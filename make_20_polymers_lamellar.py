#!/usr/bin/env python3
"""
Create a 20-chain GROMACS .gro seed structure from a single-chain .gro and .itp.

This is a morphology seed, not an X-ray-refined structure.  It places chains in an
aligned layered array so that subsequent MD annealing/NPT relaxation can be used to
approach a GIWAXS-compatible morphology.
"""
from __future__ import annotations
import argparse
import random
import re
from pathlib import Path
from typing import List, Tuple

Atom = Tuple[int, str, str, int, float, float, float]

def read_gro(path: Path) -> tuple[str, list[Atom], list[float]]:
    lines = path.read_text().splitlines()
    title = lines[0].strip()
    natoms = int(lines[1].strip())
    atoms: list[Atom] = []
    for line in lines[2:2 + natoms]:
        # Standard fixed-width .gro fields.
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        atoms.append((resid, resname, atomname, atomnr, x, y, z))
    box = [float(x) for x in lines[2 + natoms].split()]
    return title, atoms, box

def find_molecule_type(itp_path: Path) -> str:
    lines = itp_path.read_text().splitlines()
    in_mol = False
    for line in lines:
        s = line.strip()
        if not s or s.startswith(';'):
            continue
        if s.startswith('['):
            in_mol = s.lower().startswith('[ moleculetype')
            continue
        if in_mol:
            return s.split()[0]
    raise ValueError(f"Could not find [ moleculetype ] name in {itp_path}")

def write_gro(path: Path, title: str, atoms: list[Atom], box: tuple[float, float, float]) -> None:
    with path.open('w') as f:
        f.write(title[:80] + '\n')
        f.write(f"{len(atoms):5d}\n")
        for i, (resid, resname, atomname, _oldnr, x, y, z) in enumerate(atoms, start=1):
            # .gro keeps only 5 digits for atom/residue ids. This system is below 99999 atoms.
            f.write(f"{resid % 100000:5d}{resname:<5.5s}{atomname:>5.5s}{i % 100000:5d}"
                    f"{x:8.3f}{y:8.3f}{z:8.3f}\n")
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")

def write_top(path: Path, itp_name: str, molname: str, nmol: int, ff_include: str | None = None) -> None:
    with path.open('w') as f:
        f.write('; Topology for the 20-polymer seed structure\n')
        if ff_include:
            f.write(f'#include "{ff_include}"\n')
            f.write('; If your atom types are in another file, include it above the molecule .itp.\n')
        else:
            f.write('; Add your force-field/atomtypes include above this line if needed.\n')
            f.write('; Example: #include "forcefield.itp"\n')
        f.write(f'#include "{itp_name}"\n\n')
        f.write('[ system ]\n')
        f.write(f'{nmol} {molname} polymers, aligned lamellar seed\n\n')
        f.write('[ molecules ]\n')
        f.write(f'{molname:<16s} {nmol}\n')

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--gro', default='1(6).gro', help='single polymer .gro')
    p.add_argument('--itp', default='1(4).itp', help='single polymer .itp')
    p.add_argument('--nmol', type=int, default=20)
    p.add_argument('--ny', type=int, default=5, help='chains per layer along y')
    p.add_argument('--lamellar-d-nm', type=float, default=2.80,
                   help='center-to-center layer spacing in z; 2.5-3.0 nm is typical for a low-q peak near 0.21-0.25 A^-1')
    p.add_argument('--gap-y-nm', type=float, default=0.45,
                   help='extra clearance beyond the single-chain y extent')
    p.add_argument('--margin-nm', type=float, default=0.55)
    p.add_argument('--jitter-nm', type=float, default=0.03,
                   help='small random displacement to avoid a perfectly artificial lattice')
    p.add_argument('--seed', type=int, default=20260618)
    p.add_argument('--out-gro', default='20_JR4Q_lamellar_seed.gro')
    p.add_argument('--out-top', default='topol_20_JR4Q.top')
    p.add_argument('--ff-include', default=None,
                   help='optional force-field/atomtypes include file, written above the molecule itp in topol')
    args = p.parse_args()

    gro_path = Path(args.gro)
    itp_path = Path(args.itp)
    title, atoms, _box = read_gro(gro_path)
    molname = find_molecule_type(itp_path)

    xs = [a[4] for a in atoms]
    ys = [a[5] for a in atoms]
    zs = [a[6] for a in atoms]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    minz, maxz = min(zs), max(zs)
    cx = 0.5 * (minx + maxx)
    cy = 0.5 * (miny + maxy)
    cz = 0.5 * (minz + maxz)
    x_extent = maxx - minx
    y_extent = maxy - miny
    z_extent = maxz - minz

    ny = args.ny
    nz = (args.nmol + ny - 1) // ny
    if ny * nz < args.nmol:
        raise ValueError('ny*nz must fit nmol')

    y_step = y_extent + args.gap_y_nm
    z_step = max(args.lamellar_d_nm, z_extent + 0.25)
    box_x = x_extent + 2 * args.margin_nm
    box_y = (ny - 1) * y_step + y_extent + 2 * args.margin_nm
    box_z = (nz - 1) * z_step + z_extent + 2 * args.margin_nm

    rng = random.Random(args.seed)
    out_atoms: list[Atom] = []
    resid_offset = 0
    count = 0
    for iz in range(nz):
        for iy in range(ny):
            if count >= args.nmol:
                break
            # Alternate small x offsets between layers/chains, but preserve common alignment.
            x0 = args.margin_nm + x_extent / 2 + (0.15 if (iy + iz) % 2 else 0.0)
            y0 = args.margin_nm + y_extent / 2 + iy * y_step
            z0 = args.margin_nm + z_extent / 2 + iz * z_step
            dx = rng.uniform(-args.jitter_nm, args.jitter_nm)
            dy = rng.uniform(-args.jitter_nm, args.jitter_nm)
            dz = rng.uniform(-args.jitter_nm, args.jitter_nm)
            for resid, resname, atomname, atomnr, x, y, z in atoms:
                out_atoms.append((resid + resid_offset, resname, atomname, atomnr,
                                  (x - cx) + x0 + dx,
                                  (y - cy) + y0 + dy,
                                  (z - cz) + z0 + dz))
            resid_offset += max(a[0] for a in atoms)
            count += 1

    write_gro(Path(args.out_gro), f'{args.nmol} copies of {gro_path.name}: aligned lamellar seed', out_atoms, (box_x, box_y, box_z))
    write_top(Path(args.out_top), itp_path.name, molname, args.nmol, args.ff_include)

    print(f'Wrote {args.out_gro} with {len(out_atoms)} atoms, molecule type {molname}, box = {box_x:.3f} {box_y:.3f} {box_z:.3f} nm')
    print(f'Wrote {args.out_top} with [ molecules ] {molname} {args.nmol}')
    print('This is a loose ordered seed. Run energy minimization/annealing/NPT before comparing simulated scattering with GIWAXS.')

if __name__ == '__main__':
    main()
