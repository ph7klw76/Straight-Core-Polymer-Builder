#!/usr/bin/env python3
"""
interdigitated_lamellar_builder.py

Build an interdigitated lamellar starting structure for all-atom conjugated polymer
GROMACS simulations from one complete single-chain .gro and its matching .itp.

Important assumption
--------------------
The input .gro/.itp should describe ONE chemically complete polymer chain.
For the JR4Q/OEG3 and VLBL/OC8 files used here, the .itp contains 10 residues,
so --repeat-units 10 simply checks/records that the supplied chain is a 10-repeat
polymer. This script does not guess new inter-repeat chemical bonds for an
arbitrary monomer unless those bonds are already present in the .itp.

Main output
-----------
PREFIX_wholemol.gro                  GROMACS coordinate file, molecules kept whole
PREFIX_visual_atomwrapped.gro         visualization-only coordinate file, atoms wrapped
PREFIX.top                            topology for N copies of the input molecule
PREFIX_posre_backbone_soft.itp        conditional backbone position restraints
PREFIX_minim.mdp                      restrained minimization MDP
PREFIX_npt_stage1_yz_compress.mdp     directed y/z NPT compression MDP
PREFIX_npt_stage2_1bar_relax.mdp      1 bar relaxation MDP
PREFIX_run_commands.sh                suggested commands
PREFIX_summary.txt                    build summary

Coordinate convention
---------------------
x = polymer backbone / longest PCA axis
y = lamellar / side-chain interdigitation direction
z = stacking direction

python3 interdigitated_lamellar_builder.py \
  --gro polymer20.gro \
  --itp polymer20.itp \
  --prefix polymer20_lamellar \
  --chains 20 \
  --repeat-units 20 \
  --chains-along-x 1 \
  --ny 4 \
  --nz 5 \
  --packing-mode safe \
  --x-stagger 0.0 \
  --yz-pressure 300 \
  --stage1-ps 3000 \
  --stage2-ps 3000
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

AMU_TO_G_PER_NM3_TO_G_PER_CM3 = 1.66053906660e-3  # density = amu * factor / nm^3


@dataclass
class GroData:
    title: str
    atoms: list[tuple[int, str, str, int]]  # resid, resname, atomname, atomnr
    coords: np.ndarray
    box: np.ndarray | None


@dataclass
class ItpData:
    molname: str
    atom_numbers: list[int]
    atom_types: list[str]
    atom_names: list[str]
    residue_numbers: list[int]
    residue_names: list[str]
    masses: list[float]

    @property
    def mass_amu(self) -> float:
        return float(sum(self.masses))

    @property
    def n_residues(self) -> int:
        return len(set(self.residue_numbers))


def read_gro(path: Path) -> GroData:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path} does not look like a .gro file")
    try:
        natoms = int(lines[1].strip())
    except Exception as exc:
        raise ValueError(f"Could not read atom count from {path}") from exc
    if len(lines) < 3 + natoms:
        raise ValueError(f"{path} is truncated: expected {natoms} atoms")

    atoms: list[tuple[int, str, str, int]] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + natoms]:
        # GROMACS fixed-width .gro fields. Fallback splitting is unsafe for atom names,
        # so keep fixed-width parsing.
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        atoms.append((resid, resname, atomname, atomnr))
        coords.append([x, y, z])

    box = None
    if len(lines) > 2 + natoms:
        parts = lines[2 + natoms].split()
        if len(parts) >= 3:
            box = np.array([float(parts[0]), float(parts[1]), float(parts[2])], dtype=float)
    return GroData(lines[0].strip(), atoms, np.array(coords, dtype=float), box)


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def section_name(line: str) -> str | None:
    s = line.strip()
    if not (s.startswith("[") and "]" in s):
        return None
    return s[1 : s.index("]")].strip().lower()


def read_itp(path: Path) -> ItpData:
    molname: str | None = None
    section: str | None = None
    atom_numbers: list[int] = []
    atom_types: list[str] = []
    atom_names: list[str] = []
    residue_numbers: list[int] = []
    residue_names: list[str] = []
    masses: list[float] = []

    for raw in path.read_text().splitlines():
        name = section_name(raw)
        if name is not None:
            section = name
            continue
        s = strip_comment(raw)
        if not s:
            continue
        parts = s.split()
        if section == "moleculetype" and molname is None:
            molname = parts[0]
        elif section == "atoms":
            if len(parts) < 8:
                continue
            # GROMACS .itp [ atoms ]: nr type resnr residue atom cgnr charge mass
            atom_numbers.append(int(parts[0]))
            atom_types.append(parts[1])
            residue_numbers.append(int(parts[2]))
            residue_names.append(parts[3])
            atom_names.append(parts[4])
            masses.append(float(parts[7]))

    if molname is None:
        raise ValueError(f"No [ moleculetype ] found in {path}")
    if not atom_numbers:
        raise ValueError(f"No [ atoms ] section found in {path}")
    return ItpData(molname, atom_numbers, atom_types, atom_names, residue_numbers, residue_names, masses)


def pca_align(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return centered PCA coordinates, extents, rotation matrix, and original center."""
    center = coords.mean(axis=0)
    x = coords - center
    cov = x.T @ x / len(x)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vecs = vecs[:, order]
    # Keep orientation deterministic.
    for k in range(3):
        axis = vecs[:, k]
        largest_component = np.argmax(np.abs(axis))
        if axis[largest_component] < 0:
            vecs[:, k] *= -1.0
    p = x @ vecs
    mins = p.min(axis=0)
    maxs = p.max(axis=0)
    extents = maxs - mins
    p_centered = p - 0.5 * (mins + maxs)
    return p_centered, extents, vecs, center


def choose_cross_section(nchains: int, chains_along_x: int, ny: int | None, nz: int | None) -> tuple[int, int, int]:
    nx = max(1, int(chains_along_x))
    n_cross = math.ceil(nchains / nx)
    if ny is not None and nz is not None:
        if ny * nz * nx < nchains:
            raise ValueError("ny * nz * chains_along_x is smaller than --chains")
        return nx, ny, nz
    if ny is not None:
        nz = math.ceil(n_cross / ny)
        return nx, ny, nz
    if nz is not None:
        ny = math.ceil(n_cross / nz)
        return nx, ny, nz

    # Choose a near-square factor pair. For 20 this gives 4 x 5.
    best = None
    for y in range(1, n_cross + 1):
        z = math.ceil(n_cross / y)
        empty = y * z - n_cross
        aspect_penalty = abs(y - z)
        # prefer y <= z for polymer stacks, but only weakly
        score = 10 * empty + aspect_penalty + (0.2 if y > z else 0.0)
        if best is None or score < best[0]:
            best = (score, y, z)
    assert best is not None
    return nx, best[1], best[2]


def default_spacings(ext: np.ndarray, args: argparse.Namespace) -> tuple[float, float]:
    # Safe all-atom defaults: use overlap/interdigitation but do not start from the
    # final experimental pi-spacing. The NPT compression files are meant to reduce
    # y/z afterward.
    if args.y_spacing is not None:
        y_spacing = args.y_spacing
    else:
        if args.packing_mode == "safe":
            y_spacing = max(args.min_y_spacing, args.safe_y_factor * float(ext[1]))
        else:
            y_spacing = args.compact_y_spacing
    if args.z_spacing is not None:
        z_spacing = args.z_spacing
    else:
        if args.packing_mode == "safe":
            z_spacing = max(args.min_z_spacing, args.safe_z_factor * float(ext[2]))
        else:
            z_spacing = args.compact_z_spacing
    return float(y_spacing), float(z_spacing)


def build_lamella(
    pcoords: np.ndarray,
    ext: np.ndarray,
    nchains: int,
    nx: int,
    ny: int,
    nz: int,
    y_spacing: float,
    z_spacing: float,
    x_margin: float,
    yz_margin: float,
    x_stagger: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Construct interdigitated lamella in PCA frame."""
    Lx = nx * (float(ext[0]) + x_margin)
    Ly = float(ext[1]) + (ny - 1) * y_spacing + yz_margin
    Lz = float(ext[2]) + (nz - 1) * z_spacing + yz_margin
    box = np.array([Lx, Ly, Lz], dtype=float)

    coords_all: list[np.ndarray] = []
    chain_ids: list[int] = []
    count = 0
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                if count >= nchains:
                    break
                m = pcoords.copy()
                # Alternate face-to-face/back-to-back packing for side-chain interdigitation.
                alt = (ix + iy + iz) % 2
                if alt:
                    m[:, 1] *= -1.0
                    m[:, 2] *= -1.0
                    xs = +0.5 * x_stagger
                else:
                    xs = -0.5 * x_stagger

                # Keep the box x-length close to nx chain lengths. For x-stagger we intentionally
                # permit the whole molecule to cross periodic boundaries; this is normal for a seed.
                cx = (ix + 0.5) * (float(ext[0]) + x_margin) + xs
                cy = yz_margin / 2 + float(ext[1]) / 2 + iy * y_spacing
                cz = yz_margin / 2 + float(ext[2]) / 2 + iz * z_spacing
                placed = m + np.array([cx, cy, cz])
                coords_all.append(placed)
                chain_ids.extend([count] * len(pcoords))
                count += 1
            if count >= nchains:
                break
        if count >= nchains:
            break

    coords = np.vstack(coords_all)
    metadata = {
        "Lx": Lx,
        "Ly": Ly,
        "Lz": Lz,
        "y_spacing": y_spacing,
        "z_spacing": z_spacing,
        "x_stagger": x_stagger,
        "nx": float(nx),
        "ny": float(ny),
        "nz": float(nz),
    }
    return coords, np.array(chain_ids, dtype=int), box, metadata


def wrap_atoms(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    return coords - np.floor(coords / box) * box


def write_gro(path: Path, title: str, template_atoms: list[tuple[int, str, str, int]], coords_all: np.ndarray, box: np.ndarray) -> None:
    nper = len(template_atoms)
    if len(coords_all) % nper != 0:
        raise ValueError("Total atom count is not a multiple of template atom count")
    max_resid = max(a[0] for a in template_atoms)
    with path.open("w") as f:
        f.write(title[:80] + "\n")
        f.write(f"{len(coords_all):5d}\n")
        for i, xyz in enumerate(coords_all, start=1):
            j = (i - 1) % nper
            chain = (i - 1) // nper
            resid, resname, atomname, _ = template_atoms[j]
            resid2 = resid + chain * max_resid
            x, y, z = xyz
            f.write(f"{resid2 % 100000:5d}{resname:<5.5s}{atomname:>5.5s}{i % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


def select_backbone_indices(itp: ItpData, backbone_types: Iterable[str]) -> list[int]:
    types = {t.strip() for t in backbone_types if t.strip()}
    idx = [nr for nr, typ in zip(itp.atom_numbers, itp.atom_types) if typ in types]
    return idx


def write_posre(path: Path, atom_indices: list[int], fcx: float, fcy: float, fcz: float) -> None:
    with path.open("w") as f:
        f.write("; Backbone position restraints generated by interdigitated_lamellar_builder.py\n")
        f.write("; Include this file after the molecule .itp and activate with -DPOSRES_BACKBONE_SOFT.\n")
        f.write("[ position_restraints ]\n")
        f.write(";  ai  funct      fcx      fcy      fcz\n")
        for ai in atom_indices:
            f.write(f"{ai:6d}     1  {fcx:8.1f} {fcy:8.1f} {fcz:8.1f}\n")


def write_top(path: Path, itp_include: str, molname: str, nchains: int, posre_include: str, forcefield_include: str | None) -> None:
    with path.open("w") as f:
        f.write("; Topology generated by interdigitated_lamellar_builder.py\n")
        if forcefield_include:
            f.write(f'#include "{forcefield_include}"\n\n')
        else:
            f.write("; Add your force-field/atomtypes include here, before the molecule .itp.\n")
            f.write('; Example: #include "gromos54a7_atb.ff/forcefield.itp"\n\n')
        f.write(f'#include "{itp_include}"\n\n')
        f.write("#ifdef POSRES_BACKBONE_SOFT\n")
        f.write(f'#include "{posre_include}"\n')
        f.write("#endif\n\n")
        f.write("[ system ]\n")
        f.write(f"{nchains} {molname} interdigitated lamellar polymer seed\n\n")
        f.write("[ molecules ]\n")
        f.write(f"{molname:<16s} {nchains}\n")


def mdp_common() -> str:
    return """; Common nonbonded settings
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
periodic-molecules      = no
rlist                   = 1.2
coulombtype             = PME
rcoulomb                = 1.2
pme-order               = 4
fourierspacing          = 0.12
vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
constraints             = none
constraint-algorithm    = lincs
lincs-iter              = 1
lincs-order             = 4
"""


def write_minim_mdp(path: Path) -> None:
    path.write_text(
        """; Restrained minimization for interdigitated lamellar seed
; Run with: gmx_mpi grompp -f PREFIX_minim.mdp -c PREFIX_wholemol.gro -r PREFIX_wholemol.gro -p PREFIX.top -o PREFIX_em.tpr -maxwarn 3

define                  = -DPOSRES_BACKBONE_SOFT
integrator              = steep
emtol                   = 500.0
emstep                  = 0.001
nsteps                  = 100000
nstenergy               = 500
nstlog                  = 500
"""
        + mdp_common()
    )


def write_npt_stage1(path: Path, yz_pressure: float, duration_ps: float, timestep_ps: float) -> None:
    nsteps = int(round(duration_ps / timestep_ps))
    t1 = 0
    t2 = int(0.1 * duration_ps)
    t3 = int(0.6 * duration_ps)
    t4 = int(0.85 * duration_ps)
    t5 = int(duration_ps)
    path.write_text(
        f"""; Stage 1: directed y/z compaction while keeping chain-length x almost fixed
; Purpose: remove white space without melting the lamellar/backbone template.
; Use after restrained minimization.
; Run with -r reference coordinates because position restraints are active.

define                  = -DPOSRES_BACKBONE_SOFT
integrator              = md
dt                      = {timestep_ps:.4f}
nsteps                  = {nsteps}
continuation            = no

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 5000
compressed-x-precision  = 1000

{mdp_common()}

tcoupl                  = V-rescale
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = 300

; Mild annealing. Avoid 500 K because that melted the seed in earlier tests.
annealing               = single
annealing-npoints       = 5
annealing-time          = {t1} {t2} {t3} {t4} {t5}
annealing-temp          = 300 360 360 320 300

; Berendsen anisotropic is used only for pre-compression.
; x compressibility = 0 keeps the box length close to one chain length.
; y/z pressure removes lamellar/stacking white space.
pcoupl                  = Berendsen
pcoupltype              = anisotropic
tau-p                   = 5.0
compressibility         = 0.0      4.5e-5   4.5e-5   0   0   0
ref-p                   = 1.0      {yz_pressure:.1f}    {yz_pressure:.1f}    0   0   0
refcoord-scaling        = all

gen-vel                 = yes
gen-temp                = 300
gen-seed                = -1

comm-mode               = Linear
nstcomm                 = 100
comm-grps               = System
free-energy             = no
"""
    )


def write_npt_stage2(path: Path, duration_ps: float, timestep_ps: float) -> None:
    nsteps = int(round(duration_ps / timestep_ps))
    path.write_text(
        f"""; Stage 2: 1 bar relaxation after y/z pre-compression
; This is still equilibration, not final production.
; If the lamellar template drifts too much, keep -DPOSRES_BACKBONE_SOFT for this stage.

define                  = -DPOSRES_BACKBONE_SOFT
integrator              = md
dt                      = {timestep_ps:.4f}
nsteps                  = {nsteps}
continuation            = yes

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 5000
compressed-x-precision  = 1000

{mdp_common()}

tcoupl                  = V-rescale
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = 300
annealing               = no

; C-rescale is kept isotropic here because some GROMACS builds do not support
; C-rescale + anisotropic pressure coupling. For production, validate the ensemble.
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau-p                   = 10.0
compressibility         = 4.5e-5
ref-p                   = 1.0
refcoord-scaling        = all

gen-vel                 = no
comm-mode               = Linear
nstcomm                 = 100
comm-grps               = System
free-energy             = no
"""
    )


def write_commands(path: Path, prefix: str) -> None:
    text = f"""#!/usr/bin/env bash
set -euo pipefail

# 1) Restrained minimization
# Position restraints require -r.
gmx_mpi grompp \\
  -f {prefix}_minim.mdp \\
  -c {prefix}_wholemol.gro \\
  -r {prefix}_wholemol.gro \\
  -p {prefix}.top \\
  -o {prefix}_em.tpr \\
  -maxwarn 3

gmx_mpi mdrun -deffnm {prefix}_em

# 2) Directed y/z compression. Monitor Box-X, Box-Y, Box-Z, Density.
gmx_mpi grompp \\
  -f {prefix}_npt_stage1_yz_compress.mdp \\
  -c {prefix}_em.gro \\
  -r {prefix}_em.gro \\
  -p {prefix}.top \\
  -o {prefix}_npt_yz_compress.tpr \\
  -maxwarn 3

gmx_mpi mdrun -deffnm {prefix}_npt_yz_compress

# 3) 1 bar relaxation.
gmx_mpi grompp \\
  -f {prefix}_npt_stage2_1bar_relax.mdp \\
  -c {prefix}_npt_yz_compress.gro \\
  -r {prefix}_npt_yz_compress.gro \\
  -p {prefix}.top \\
  -o {prefix}_npt_1bar_relax.tpr \\
  -maxwarn 3

gmx_mpi mdrun -deffnm {prefix}_npt_1bar_relax

# Optional checks:
# gmx_mpi energy -f {prefix}_npt_yz_compress.edr -o {prefix}_box_density.xvg
# Choose: Box-X Box-Y Box-Z Density Pressure
"""
    path.write_text(text)
    path.chmod(0o755)


def main() -> None:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Create dense/interdigitated lamellar GROMACS seeds and supporting topology/MDP files.",
    )
    ap.add_argument("--gro", required=True, help="input .gro for ONE complete polymer chain")
    ap.add_argument("--itp", required=True, help="matching .itp for ONE complete polymer chain")
    ap.add_argument("--prefix", default="lamellar_seed", help="output prefix")
    ap.add_argument("--chains", "-x", type=int, default=20, help="number of polymer chains")
    ap.add_argument("--repeat-units", "-y", type=int, default=10, help="expected repeat/residue count in the input chain")
    ap.add_argument("--chains-along-x", "-z", type=int, default=1, help="number of chain lengths along the longest box direction")
    ap.add_argument("--ny", type=int, default=None, help="number of lamellar rows; auto if omitted")
    ap.add_argument("--nz", type=int, default=None, help="number of stacking rows; auto if omitted")
    ap.add_argument("--packing-mode", choices=["safe", "compact"], default="safe", help="safe avoids extreme initial overlaps; compact starts closer but may require careful minimization")
    ap.add_argument("--y-spacing", type=float, default=None, help="manual lamellar row spacing in nm")
    ap.add_argument("--z-spacing", type=float, default=None, help="manual stacking row spacing in nm")
    ap.add_argument("--safe-y-factor", type=float, default=0.78, help="safe y spacing = this factor times molecule PCA y-width")
    ap.add_argument("--safe-z-factor", type=float, default=1.00, help="safe z spacing = this factor times molecule PCA z-width")
    ap.add_argument("--min-y-spacing", type=float, default=4.50, help="minimum safe y spacing in nm when auto-selected")
    ap.add_argument("--min-z-spacing", type=float, default=0.75, help="minimum safe z spacing in nm when auto-selected")
    ap.add_argument("--compact-y-spacing", type=float, default=3.00, help="compact-mode y spacing in nm")
    ap.add_argument("--compact-z-spacing", type=float, default=0.70, help="compact-mode z spacing in nm")
    ap.add_argument("--x-margin", type=float, default=0.25, help="extra box length per chain along x in nm")
    ap.add_argument("--yz-margin", type=float, default=0.25, help="extra box margin in y and z in nm")
    ap.add_argument("--x-stagger", type=float, default=0.0, help="alternating x stagger in nm; use 0 for clean visualization, 2-5 for stronger interdigitation")
    ap.add_argument("--backbone-types", default="CAro,CPos,S,NOpt,NTer,OEOpt", help="comma-separated atom types to restrain")
    ap.add_argument("--posre-fc", default="20,20,50", help="fcx,fcy,fcz for soft backbone restraints in kJ mol^-1 nm^-2")
    ap.add_argument("--forcefield-include", default="gromos54a7_atb.ff/forcefield.itp", help="force-field include line for .top; use empty string to omit")
    ap.add_argument("--yz-pressure", type=float, default=300.0, help="stage-1 y/z compression pressure in bar")
    ap.add_argument("--stage1-ps", type=float, default=2000.0, help="stage-1 NPT duration in ps")
    ap.add_argument("--stage2-ps", type=float, default=2000.0, help="stage-2 NPT duration in ps")
    ap.add_argument("--dt", type=float, default=0.0005, help="MD timestep in ps for NPT stages")
    args = ap.parse_args()

    gro_path = Path(args.gro)
    itp_path = Path(args.itp)
    gro = read_gro(gro_path)
    itp = read_itp(itp_path)
    if len(gro.atoms) != len(itp.atom_numbers):
        raise ValueError(
            f"Atom count mismatch: {gro_path} has {len(gro.atoms)} atoms but {itp_path} has {len(itp.atom_numbers)} atoms"
        )

    if itp.n_residues != args.repeat_units:
        print(
            f"WARNING: --repeat-units is {args.repeat_units}, but {itp_path.name} contains {itp.n_residues} unique residues. "
            "The script will still pack the supplied chain as one molecule; it will not rewrite the chemistry.",
            file=sys.stderr,
        )

    pcoords, ext, _, _ = pca_align(gro.coords)
    nx, ny, nz = choose_cross_section(args.chains, args.chains_along_x, args.ny, args.nz)
    y_spacing, z_spacing = default_spacings(ext, args)
    coords_all, chain_ids, box, meta = build_lamella(
        pcoords,
        ext,
        args.chains,
        nx,
        ny,
        nz,
        y_spacing,
        z_spacing,
        args.x_margin,
        args.yz_margin,
        args.x_stagger,
    )

    prefix = args.prefix
    whole_gro = Path(f"{prefix}_wholemol.gro")
    visual_gro = Path(f"{prefix}_visual_atomwrapped.gro")
    posre = Path(f"{prefix}_posre_backbone_soft.itp")
    top = Path(f"{prefix}.top")
    minim = Path(f"{prefix}_minim.mdp")
    npt1 = Path(f"{prefix}_npt_stage1_yz_compress.mdp")
    npt2 = Path(f"{prefix}_npt_stage2_1bar_relax.mdp")
    commands = Path(f"{prefix}_run_commands.sh")
    summary = Path(f"{prefix}_summary.txt")

    title = f"{args.chains} {itp.molname} chains; interdigitated lamellar seed; repeat-units={args.repeat_units}"
    write_gro(whole_gro, title, gro.atoms, coords_all, box)
    write_gro(visual_gro, title + " atom-wrapped", gro.atoms, wrap_atoms(coords_all, box), box)

    backbone_types = [x.strip() for x in args.backbone_types.split(",")]
    bidx = select_backbone_indices(itp, backbone_types)
    if not bidx:
        print("WARNING: no backbone atoms selected. Check --backbone-types.", file=sys.stderr)
    fcx, fcy, fcz = [float(x.strip()) for x in args.posre_fc.split(",")]
    write_posre(posre, bidx, fcx, fcy, fcz)

    ff = args.forcefield_include.strip() or None
    write_top(top, itp_path.name, itp.molname, args.chains, posre.name, ff)
    write_minim_mdp(minim)
    write_npt_stage1(npt1, args.yz_pressure, args.stage1_ps, args.dt)
    write_npt_stage2(npt2, args.stage2_ps, args.dt)
    write_commands(commands, prefix)

    volume = float(np.prod(box))
    density = itp.mass_amu * args.chains * AMU_TO_G_PER_NM3_TO_G_PER_CM3 / volume
    summary_text = f"""Build summary for {prefix}

Input GRO: {gro_path.name}
Input ITP: {itp_path.name}
Molecule type: {itp.molname}
Input atoms per chain: {len(gro.atoms)}
Residues/repeat units detected in ITP: {itp.n_residues}
Requested repeat-units label/check: {args.repeat_units}
Number of chains: {args.chains}
Chains along x: {nx}
Lamellar rows y: {ny}
Stacking rows z: {nz}

PCA single-chain extents / nm:
  x backbone length = {ext[0]:.4f}
  y side-chain width = {ext[1]:.4f}
  z thickness        = {ext[2]:.4f}

Box / nm:
  Lx = {box[0]:.4f}
  Ly = {box[1]:.4f}
  Lz = {box[2]:.4f}
  Volume = {volume:.2f} nm^3
Estimated initial density = {density:.3f} g/cm^3

Packing:
  y_spacing = {y_spacing:.4f} nm
  z_spacing = {z_spacing:.4f} nm
  x_stagger = {args.x_stagger:.4f} nm
  packing_mode = {args.packing_mode}

Backbone restraint atom types: {', '.join(backbone_types)}
Backbone restraint atoms selected: {len(bidx)}
Position restraint force constants: fcx={fcx}, fcy={fcy}, fcz={fcz} kJ mol^-1 nm^-2

Recommended sequence:
  1. {prefix}_minim.mdp
  2. {prefix}_npt_stage1_yz_compress.mdp
  3. {prefix}_npt_stage2_1bar_relax.mdp

Use {prefix}_wholemol.gro for grompp and {prefix}_visual_atomwrapped.gro for quick visualization.
"""
    summary.write_text(summary_text)

    print(summary_text)
    print("Created files:")
    for p in [whole_gro, visual_gro, top, posre, minim, npt1, npt2, commands, summary]:
        print(f"  {p}")


if __name__ == "__main__":
    main()
