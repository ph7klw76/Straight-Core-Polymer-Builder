#!/usr/bin/env python3
"""
Rotate and translate a GROMACS .gro molecule so that its long molecular/core axis
is aligned with the x-axis, then fit it into a rectangular box with a chosen
padding from all box edges.

Default behavior:
    - Uses heavy atoms as a practical proxy for the aromatic/core axis.
    - Performs PCA/SVD alignment so the longest heavy-atom axis becomes x.
    - Translates the molecule so min(x,y,z) = padding.
    - Sets the box size to molecular span + 2 * padding.
    - Preserves atom/residue columns and velocity columns if present.

Example:
    python align_gro_xaxis_fitbox.py input.gro -o output.gro

Example with 0.2 nm padding:
    python align_gro_xaxis_fitbox.py input.gro -o output.gro --padding 0.2

Example using explicit aromatic-core atom names:
    python align_gro_xaxis_fitbox.py input.gro -o output.gro --core-names C1,C2,C3,C4,C5,C6,N1,N2
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class GroAtom:
    line: str
    prefix: str
    suffix: str
    atomname: str
    atomnr: int
    coord: np.ndarray


def is_hydrogen_atom(atomname: str) -> bool:
    """
    Detect hydrogen-like atom names.

    Handles names such as:
        H, H1, H12, HA, HB2, 1H, 2HA

    For normal organic .gro files this is usually sufficient.
    """
    name = atomname.strip().upper()
    return bool(re.match(r"^\d*H", name))


def parse_gro_atom_line(line: str, fallback_atomnr: int) -> GroAtom:
    """
    Parse a standard .gro atom line.

    Standard GRO fixed-width columns:
        0:5    residue number
        5:10   residue name
        10:15  atom name
        15:20  atom number
        20:28  x
        28:36  y
        36:44  z
        44:     optional velocities
    """
    if len(line) >= 44:
        prefix = line[:20]
        suffix = line[44:]

        atomname = line[10:15].strip()

        try:
            atomnr = int(line[15:20])
        except ValueError:
            atomnr = fallback_atomnr

        try:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError as exc:
            raise ValueError(f"Could not parse coordinates from line:\n{line}") from exc

        return GroAtom(
            line=line,
            prefix=prefix,
            suffix=suffix,
            atomname=atomname,
            atomnr=atomnr,
            coord=np.array([x, y, z], dtype=float),
        )

    raise ValueError(
        "Atom line is too short for standard .gro fixed-width format:\n"
        f"{line}"
    )


def read_gro(path: Path) -> tuple[str, list[GroAtom], str]:
    lines = path.read_text().splitlines()

    if len(lines) < 3:
        raise ValueError("This file is too short to be a valid .gro file.")

    title = lines[0]

    try:
        natoms = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError("Could not read atom count from the second line.") from exc

    expected_line_count = natoms + 3
    if len(lines) < expected_line_count:
        raise ValueError(
            f"Atom count mismatch: header says {natoms} atoms, "
            f"but file has only {len(lines) - 3 + 1} possible atom lines."
        )

    atom_lines = lines[2 : 2 + natoms]
    box_line = lines[2 + natoms]

    atoms = [
        parse_gro_atom_line(line, fallback_atomnr=i + 1)
        for i, line in enumerate(atom_lines)
    ]

    return title, atoms, box_line


def choose_core_mask(
    atoms: list[GroAtom],
    selection: str = "heavy",
    core_names: str | None = None,
) -> np.ndarray:
    """
    Select atoms used to define the long axis.

    Recommended:
        - default: heavy atoms
        - if you know the exact aromatic-core atom names, use --core-names
    """
    if core_names:
        requested = {
            name.strip()
            for name in core_names.split(",")
            if name.strip()
        }

        mask = np.array(
            [atom.atomname.strip() in requested for atom in atoms],
            dtype=bool,
        )

        if mask.sum() < 3:
            raise ValueError(
                "Fewer than 3 atoms matched --core-names. "
                "Please check the atom names in your .gro file."
            )

        return mask

    if selection == "all":
        return np.ones(len(atoms), dtype=bool)

    if selection == "heavy":
        return np.array(
            [not is_hydrogen_atom(atom.atomname) for atom in atoms],
            dtype=bool,
        )

    raise ValueError(f"Unknown selection mode: {selection}")


def align_long_axis_to_x(
    coords: np.ndarray,
    core_mask: np.ndarray,
    atom_numbers: np.ndarray,
    orient_by_atom_number: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Align the longest PCA axis of the selected core atoms to the x-axis.

    Returns:
        rotated_coords
        rotation_axes
        singular_values
        core_center
    """
    core_coords = coords[core_mask]

    if core_coords.shape[0] < 3:
        raise ValueError("Need at least 3 selected atoms for PCA alignment.")

    core_center = core_coords.mean(axis=0)
    centered_core = core_coords - core_center

    # SVD gives principal directions.
    # Vt[0] is the longest axis in the original coordinate basis.
    _, singular_values, Vt = np.linalg.svd(centered_core, full_matrices=False)

    # Transform coordinates into principal-component coordinates.
    # After this, PC1 -> x, PC2 -> y, PC3 -> z.
    rotated = (coords - core_center) @ Vt.T

    # Make the x-axis direction deterministic.
    # This does not change the fit, but makes output orientation reproducible.
    if orient_by_atom_number:
        core_indices = np.where(core_mask)[0]
        sorted_local = np.argsort(atom_numbers[core_mask])

        ncore = len(sorted_local)
        if ncore >= 2:
            k = max(1, min(max(10, int(0.05 * ncore)), ncore // 2))

            early_mean_x = rotated[core_indices[sorted_local[:k]], 0].mean()
            late_mean_x = rotated[core_indices[sorted_local[-k:]], 0].mean()

            if late_mean_x < early_mean_x:
                rotated[:, 0] *= -1.0
                Vt[0, :] *= -1.0

    return rotated, Vt, singular_values, core_center


def fit_box_with_padding(
    rotated_coords: np.ndarray,
    padding: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Translate coordinates so all coordinates are positive and the molecule has
    exactly `padding` nm clearance from each box edge.

    Returns:
        final_coords
        box_lengths
        final_min
        final_max
    """
    if padding <= 0:
        raise ValueError("Padding must be greater than zero.")

    mins = rotated_coords.min(axis=0)
    maxs = rotated_coords.max(axis=0)

    span = maxs - mins
    box_lengths = span + 2.0 * padding

    shift = padding - mins
    final_coords = rotated_coords + shift

    final_min = final_coords.min(axis=0)
    final_max = final_coords.max(axis=0)

    return final_coords, box_lengths, final_min, final_max


def write_gro(
    output_path: Path,
    title: str,
    atoms: list[GroAtom],
    new_coords: np.ndarray,
    box_lengths: np.ndarray,
    padding: float,
) -> None:
    output_lines: list[str] = []

    output_lines.append(
        f"{title} | long-axis aligned to x; fitted box; padding {padding:g} nm"
    )
    output_lines.append(f"{len(atoms):5d}")

    for atom, coord in zip(atoms, new_coords):
        x, y, z = coord

        # GRO coordinates use 8.3f formatting.
        # Preserve original prefix and suffix, including velocities if present.
        output_lines.append(
            f"{atom.prefix}{x:8.3f}{y:8.3f}{z:8.3f}{atom.suffix}"
        )

    # Orthorhombic box.
    output_lines.append(
        f"{box_lengths[0]:10.5f}{box_lengths[1]:10.5f}{box_lengths[2]:10.5f}"
    )

    output_path.write_text("\n".join(output_lines) + "\n")


def process_gro(
    input_path: Path,
    output_path: Path,
    padding: float = 0.2,
    selection: str = "heavy",
    core_names: str | None = None,
    orient_by_atom_number: bool = True,
) -> None:
    title, atoms, _old_box = read_gro(input_path)

    coords = np.array([atom.coord for atom in atoms], dtype=float)
    atom_numbers = np.array([atom.atomnr for atom in atoms], dtype=int)

    core_mask = choose_core_mask(
        atoms=atoms,
        selection=selection,
        core_names=core_names,
    )

    rotated_coords, rotation_axes, singular_values, core_center = align_long_axis_to_x(
        coords=coords,
        core_mask=core_mask,
        atom_numbers=atom_numbers,
        orient_by_atom_number=orient_by_atom_number,
    )

    final_coords, box_lengths, final_min, final_max = fit_box_with_padding(
        rotated_coords=rotated_coords,
        padding=padding,
    )

    write_gro(
        output_path=output_path,
        title=title,
        atoms=atoms,
        new_coords=final_coords,
        box_lengths=box_lengths,
        padding=padding,
    )

    # Re-check rounded values as actually written to .gro.
    rounded_coords = np.round(final_coords, 3)
    rounded_min = rounded_coords.min(axis=0)
    rounded_max = rounded_coords.max(axis=0)

    variance = singular_values**2
    variance_percent = variance / variance.sum() * 100.0

    print("Finished.")
    print(f"Input file:  {input_path}")
    print(f"Output file: {output_path}")
    print()
    print(f"Number of atoms:       {len(atoms)}")
    print(f"Atoms used for axis:   {int(core_mask.sum())}")
    print(f"Selection mode:        {selection if core_names is None else 'core-names'}")
    print()
    print("Principal-axis variance percentage:")
    print(f"  PC1 / x-axis: {variance_percent[0]:.3f}%")
    print(f"  PC2 / y-axis: {variance_percent[1]:.3f}%")
    print(f"  PC3 / z-axis: {variance_percent[2]:.3f}%")
    print()
    print("Original-basis direction of the new x-axis:")
    print(f"  [{rotation_axes[0,0]: .6f}, {rotation_axes[0,1]: .6f}, {rotation_axes[0,2]: .6f}]")
    print()
    print("Final fitted box / nm:")
    print(f"  x: {box_lengths[0]:.5f}")
    print(f"  y: {box_lengths[1]:.5f}")
    print(f"  z: {box_lengths[2]:.5f}")
    print()
    print("Final coordinate range after rounding to .gro precision / nm:")
    print(f"  min: [{rounded_min[0]:.3f}, {rounded_min[1]:.3f}, {rounded_min[2]:.3f}]")
    print(f"  max: [{rounded_max[0]:.3f}, {rounded_max[1]:.3f}, {rounded_max[2]:.3f}]")
    print()
    print("Final box line:")
    print(f"  {box_lengths[0]:10.5f}{box_lengths[1]:10.5f}{box_lengths[2]:10.5f}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Rotate a .gro molecule so its long core axis lies along x, "
            "then translate it into a fitted positive-coordinate box."
        )
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Input .gro file.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output .gro file. Default: <input_stem>_xaxis_centered_0p2nm.gro",
    )

    parser.add_argument(
        "--padding",
        type=float,
        default=0.2,
        help="Padding from molecule to box edge in nm. Default: 0.2",
    )

    parser.add_argument(
        "--selection",
        choices=["heavy", "all"],
        default="heavy",
        help=(
            "Atoms used to define the long axis. "
            "'heavy' excludes hydrogen-like atom names. Default: heavy"
        ),
    )

    parser.add_argument(
        "--core-names",
        type=str,
        default=None,
        help=(
            "Comma-separated exact atom names to define the aromatic core, "
            "for example: C1,C2,C3,C4,C5,C6,N1,N2. "
            "If supplied, this overrides --selection."
        ),
    )

    parser.add_argument(
        "--no-index-orient",
        action="store_true",
        help=(
            "Do not force later atom numbers toward +x. "
            "The long axis will still be aligned with x, but the sign may flip."
        ),
    )

    args = parser.parse_args()

    input_path: Path = args.input

    if args.output is None:
        output_path = input_path.with_name(
            f"{input_path.stem}_xaxis_centered_{str(args.padding).replace('.', 'p')}nm.gro"
        )
    else:
        output_path = args.output

    process_gro(
        input_path=input_path,
        output_path=output_path,
        padding=args.padding,
        selection=args.selection,
        core_names=args.core_names,
        orient_by_atom_number=not args.no_index_orient,
    )


if __name__ == "__main__":
    main()
