#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse, math
import numpy as np

Atom = tuple[int,str,str,int,float,float,float]

def read_gro(path: Path):
    lines = path.read_text().splitlines()
    natoms = int(lines[1].strip())
    atoms=[]; coords=[]
    for line in lines[2:2+natoms]:
        resid=int(line[0:5]); resname=line[5:10].strip(); atomname=line[10:15].strip(); atomnr=int(line[15:20])
        x=float(line[20:28]); y=float(line[28:36]); z=float(line[36:44])
        atoms.append((resid,resname,atomname,atomnr))
        coords.append([x,y,z])
    return lines[0].strip(), atoms, np.array(coords,float)

def molecule_name_from_itp(itp: Path):
    in_mol=False
    for line in itp.read_text().splitlines():
        s=line.strip()
        if not s or s.startswith(';'):
            continue
        if s.startswith('['):
            in_mol=s.lower().startswith('[ moleculetype')
            continue
        if in_mol:
            return s.split()[0]
    raise RuntimeError('No [ moleculetype ] found')

def align_to_pca(coords):
    ctr=coords.mean(axis=0)
    X=coords-ctr
    cov=X.T@X/len(X)
    w,V=np.linalg.eigh(cov)
    idx=np.argsort(w)[::-1]
    V=V[:,idx]
    # Fix signs so PC axes resemble the original +x,+y,+z directions.
    for k in range(3):
        if V[:,k][k] < 0:
            V[:,k] *= -1
    Y=X@V
    return Y, V, ctr

def write_gro(path: Path, title: str, template_atoms, coords_all, box):
    with path.open('w') as f:
        f.write(title[:80]+'\n')
        f.write(f'{len(coords_all):5d}\n')
        nper=len(template_atoms)
        max_resid=max(a[0] for a in template_atoms)
        for i,xyz in enumerate(coords_all, start=1):
            j=(i-1)%nper
            chain=(i-1)//nper
            resid,resname,atomname,_=template_atoms[j]
            resid2=resid+chain*max_resid
            x,y,z=xyz
            f.write(f'{resid2%100000:5d}{resname:<5.5s}{atomname:>5.5s}{i%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n')
        f.write(f'{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n')

def write_top(path: Path, itp_name: str, molname: str, nmol: int):
    with path.open('w') as f:
        f.write('; Dense lamellar rectangular seed topology\n')
        f.write('; Add the force-field/atomtypes include above the molecule itp if needed.\n')
        f.write('; Example: #include "forcefield.itp"\n')
        f.write(f'#include "{itp_name}"\n\n')
        f.write('[ system ]\n')
        f.write(f'{nmol} {molname} dense lamellar rectangular seed\n\n')
        f.write('[ molecules ]\n')
        f.write(f'{molname:<16s} {nmol}\n')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--gro', default='1(6).gro')
    ap.add_argument('--itp', default='1(4).itp')
    ap.add_argument('--out-gro', default='20_JR4Q_dense_lamellar_rect.gro')
    ap.add_argument('--out-top', default='topol_20_JR4Q_dense_rect.top')
    ap.add_argument('--ny', type=int, default=4, help='number of lamellar rows')
    ap.add_argument('--nz', type=int, default=5, help='number of pi-stack rows')
    ap.add_argument('--y-spacing', type=float, default=4.50, help='nm, dense interdigitated lamellar spacing for the uploaded all-atom chain')
    ap.add_argument('--z-spacing', type=float, default=0.75, help='nm, safe stack spacing for an unrelaxed all-atom chain; raw 0.36-0.38 nm GIWAXS pi-spacing needs relaxation/side-chain interdigitation')
    ap.add_argument('--x-margin', type=float, default=0.25, help='nm total margin added to the PCA chain length')
    ap.add_argument('--yz-margin', type=float, default=0.20, help='nm total margin in y/z')
    ap.add_argument('--x-stagger', type=float, default=5.00, help='nm alternating chain shift along x for side-chain interdigitation')
    args=ap.parse_args()
    nmol=args.ny*args.nz
    title, atoms, coords = read_gro(Path(args.gro))
    molname=molecule_name_from_itp(Path(args.itp))
    Y, V, ctr = align_to_pca(coords)
    mins=Y.min(axis=0); maxs=Y.max(axis=0); ext=maxs-mins
    # Center molecule on zero in the PCA frame.
    Yc=Y-0.5*(mins+maxs)
    Lx=ext[0]+args.x_margin
    Ly=ext[1]+(args.ny-1)*args.y_spacing+args.yz_margin
    Lz=ext[2]+(args.nz-1)*args.z_spacing+args.yz_margin
    coords_all=[]
    for iy in range(args.ny):
        for iz in range(args.nz):
            m=Yc.copy()
            # Alternate orientation to promote interdigitation of side chains.
            if (iy+iz)%2:
                m[:,1]*=-1
                m[:,2]*=-1
                xshift=+0.5*args.x_stagger
            else:
                xshift=-0.5*args.x_stagger
            center=np.array([Lx/2 + xshift, args.yz_margin/2 + ext[1]/2 + iy*args.y_spacing, args.yz_margin/2 + ext[2]/2 + iz*args.z_spacing])
            coords_all.append(m+center)
    coords_all=np.vstack(coords_all)
    write_gro(Path(args.out_gro), f'{nmol} {molname}: dense rectangular lamellar seed; Lx fits one chain', atoms, coords_all, (Lx,Ly,Lz))
    write_top(Path(args.out_top), Path(args.itp).name, molname, nmol)
    vol=Lx*Ly*Lz
    print(f'Wrote {args.out_gro}: {nmol} molecules, {len(coords_all)} atoms')
    print(f'Box nm: {Lx:.3f} {Ly:.3f} {Lz:.3f}; volume {vol:.1f} nm^3')
    print(f'PCA single-chain extents nm: x={ext[0]:.3f}, y={ext[1]:.3f}, z={ext[2]:.3f}')
    print(f'Spacings nm: y={args.y_spacing:.3f}, z={args.z_spacing:.3f}, x stagger={args.x_stagger:.3f}')
    print(f'Wrote {args.out_top}; molecule type {molname}')

if __name__=='__main__':
    main()
