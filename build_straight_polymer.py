#!/usr/bin/env python3
"""
build_straight_polymer.py

Build an N-repeat straight-core polymer from a monomer PDB and a matching
GROMACS/ATB/GROMOS-style .itp topology.

Typical use:

  python build_straight_polymer.py \
      --pdb monomer.pdb \
      --itp monomer.itp \
      --n 10 \
      --start C25 \
      --end C21 \
      --out-prefix polymer_10mer

Definitions:
  --start is the atom in each repeat that bonds to the PREVIOUS repeat.
  --end   is the atom in each repeat that bonds to the NEXT repeat.

For repeat i -> i+1, the program creates:
  end(i) -- start(i+1)

The program assumes the PDB atom order matches the [ atoms ] order in the .itp.
It aborts if atom counts or element identities disagree.

Outputs:
  <out-prefix>.pdb
  <out-prefix>.itp
  <out-prefix>_validation.txt

Important scientific caveat:
  The inter-repeat bonded terms are constructed from generic, GROMOS/ATB-like
  alkyl parameters. This is suitable for building a working starting topology,
  but should be validated with gmx grompp, minimization, and test MD.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

Vec3 = Tuple[float, float, float]

# -----------------------------------------------------------------------------
# Vector helpers
# -----------------------------------------------------------------------------

def vsub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vadd(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vscale(a: Vec3, s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a: Vec3) -> float:
    return math.sqrt(dot(a, a))


def normalize(a: Vec3) -> Vec3:
    n = norm(a)
    if n < 1.0e-12:
        raise ValueError("Cannot normalize a near-zero vector")
    return (a[0] / n, a[1] / n, a[2] / n)


def dist(a: Vec3, b: Vec3) -> float:
    return norm(vsub(a, b))

# -----------------------------------------------------------------------------
# Data models
# -----------------------------------------------------------------------------

@dataclass
class PDBAtom:
    serial: int
    name: str
    resname: str
    resseq: int
    coord: Vec3
    elem: str


@dataclass
class ITPAtom:
    nr: int
    atype: str
    resnr: int
    resid: str
    atom: str
    cgnr: int
    charge: float
    mass: float
    line: int


@dataclass
class TopItem:
    line: int
    parts: List[str]
    raw: str


@dataclass
class PolyAtom:
    repeat: int
    old_index: int  # 0-based monomer atom index
    atype: str
    resid: str
    atom: str
    elem: str
    charge: float
    mass: float
    coord: Vec3

# -----------------------------------------------------------------------------
# PDB and ITP parsers
# -----------------------------------------------------------------------------

def atom_element_from_name(name: str) -> str:
    m = re.match(r"[A-Za-z]+", name.strip())
    if not m:
        return ""
    return m.group(0)[0].upper()


def parse_pdb(path: str) -> Tuple[List[PDBAtom], Set[Tuple[int, int]]]:
    atoms: List[PDBAtom] = []
    conect: Set[Tuple[int, int]] = set()
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                name = line[12:16].strip()
                resname = line[17:21].strip() or "MOL"
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    resseq = 1
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = line[76:78].strip()
                if not elem:
                    elem = atom_element_from_name(name)
                atoms.append(PDBAtom(serial, name, resname, resseq, (x, y, z), elem.upper()))
            elif line.startswith("CONECT"):
                fields = []
                for i in range(6, len(line.rstrip()), 5):
                    s = line[i:i+5].strip()
                    if s:
                        fields.append(int(s))
                if fields:
                    a = fields[0] - 1
                    for b1 in fields[1:]:
                        b = b1 - 1
                        if a != b:
                            conect.add(tuple(sorted((a, b))))
    return atoms, conect


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def read_itp_sections(path: str):
    preamble: List[Tuple[int, str]] = []
    sections: List[dict] = []
    cur = None
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            line = line.rstrip("\n")
            m = re.match(r"^\s*\[\s*([^\]]+)\s*\]\s*$", line)
            if m:
                cur = {"name": m.group(1).strip().lower(), "start": ln, "lines": []}
                sections.append(cur)
            elif cur is not None:
                cur["lines"].append((ln, line))
            else:
                preamble.append((ln, line))
    by_name: Dict[str, List[dict]] = defaultdict(list)
    for s in sections:
        by_name[s["name"]].append(s)
    return preamble, sections, by_name


def parse_itp_atoms(section: dict) -> List[ITPAtom]:
    atoms: List[ITPAtom] = []
    for ln, line in section["lines"]:
        data = strip_comment(line)
        if not data:
            continue
        parts = data.split()
        if len(parts) >= 8 and parts[0].isdigit():
            atoms.append(
                ITPAtom(
                    nr=int(parts[0]),
                    atype=parts[1],
                    resnr=int(parts[2]),
                    resid=parts[3],
                    atom=parts[4],
                    cgnr=int(parts[5]),
                    charge=float(parts[6]),
                    mass=float(parts[7]),
                    line=ln,
                )
            )
    return atoms


def parse_numeric_items(section: dict) -> List[TopItem]:
    out: List[TopItem] = []
    for ln, line in section["lines"]:
        data = strip_comment(line)
        if not data:
            continue
        parts = data.split()
        if parts and re.match(r"^-?\d+$", parts[0]):
            out.append(TopItem(line=ln, parts=parts, raw=line))
    return out

# -----------------------------------------------------------------------------
# Topology helpers
# -----------------------------------------------------------------------------

def build_adjacency(n: int, bond_items: Sequence[TopItem]) -> List[Set[int]]:
    adj = [set() for _ in range(n)]
    for item in bond_items:
        i = int(item.parts[0]) - 1
        j = int(item.parts[1]) - 1
        adj[i].add(j)
        adj[j].add(i)
    return adj


def resolve_atom(selector: str, itp_atoms: Sequence[ITPAtom]) -> int:
    """Resolve a 1-based atom number or unique atom name to a 0-based index."""
    if selector.isdigit():
        idx = int(selector) - 1
        if idx < 0 or idx >= len(itp_atoms):
            raise ValueError(f"Atom index {selector} is out of range")
        return idx
    matches = [i for i, a in enumerate(itp_atoms) if a.atom == selector]
    if not matches:
        raise ValueError(f"No atom named {selector!r} found in ITP [ atoms ]")
    if len(matches) > 1:
        raise ValueError(f"Atom name {selector!r} is not unique; use a 1-based atom number")
    return matches[0]


def attached_hydrogens(atom_idx: int, adj: Sequence[Set[int]], itp_atoms: Sequence[ITPAtom]) -> List[int]:
    return [j for j in sorted(adj[atom_idx]) if atom_element_from_name(itp_atoms[j].atom) == "H"]


def pick_hydrogen_by_direction(
    parent_idx: int,
    hlist: Sequence[int],
    direction: Vec3,
    pdb_atoms: Sequence[PDBAtom],
) -> int:
    if not hlist:
        raise ValueError(f"Atom {parent_idx+1} has no attached hydrogens to remove")
    d = normalize(direction)
    best_h = hlist[0]
    best_score = -1.0e99
    parent_xyz = pdb_atoms[parent_idx].coord
    for h in hlist:
        v = normalize(vsub(pdb_atoms[h].coord, parent_xyz))
        score = dot(v, d)
        if score > best_score:
            best_score = score
            best_h = h
    return best_h


def canonical_pair(i: int, j: int) -> Tuple[int, int]:
    return tuple(sorted((i, j)))


def angle_key(t: Tuple[int, int, int]) -> Tuple[int, int, int]:
    i, j, k = t
    return min((i, j, k), (k, j, i))


def dih_key(t: Tuple[int, int, int, int]) -> Tuple[int, int, int, int]:
    a, b, c, d = t
    return min((a, b, c, d), (d, c, b, a))


def replicate_items(
    items: Sequence[TopItem],
    nidx: int,
    nrep: int,
    old_to_new: Dict[Tuple[int, int], int],
) -> List[Tuple[List[int], List[str]]]:
    out: List[Tuple[List[int], List[str]]] = []
    for r in range(nrep):
        for item in items:
            old_idxs = [int(x) - 1 for x in item.parts[:nidx]]
            if all((r, idx) in old_to_new for idx in old_idxs):
                new_idxs = [old_to_new[(r, idx)] + 1 for idx in old_idxs]
                rest = item.parts[nidx:]
                out.append((new_idxs, rest))
    return out

# -----------------------------------------------------------------------------
# Polymer construction
# -----------------------------------------------------------------------------

def build_polymer(args) -> None:
    pdb_atoms, _pdb_conect = parse_pdb(args.pdb)
    _pre, sections, by = read_itp_sections(args.itp)

    required = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals"]
    for name in required:
        if name not in by:
            raise ValueError(f"Missing required [{name}] section in {args.itp}")
    if len(by["dihedrals"]) < 2:
        raise ValueError("Expected two [ dihedrals ] sections: improper and proper")

    itp_atoms = parse_itp_atoms(by["atoms"][0])
    bonds_items = parse_numeric_items(by["bonds"][0])
    pairs_items = parse_numeric_items(by["pairs"][0])
    angles_items = parse_numeric_items(by["angles"][0])
    dih_imp_items = parse_numeric_items(by["dihedrals"][0])
    dih_prop_items = parse_numeric_items(by["dihedrals"][1])
    excl_items = parse_numeric_items(by["exclusions"][0]) if "exclusions" in by else []

    if len(pdb_atoms) != len(itp_atoms):
        raise ValueError(
            f"PDB atom count ({len(pdb_atoms)}) does not match ITP atom count ({len(itp_atoms)}). "
            "This program requires monomer PDB and ITP to have the same atom order."
        )

    # Check element/order consistency.
    element_mismatches = []
    for i, (pa, ia) in enumerate(zip(pdb_atoms, itp_atoms), start=1):
        e_itp = atom_element_from_name(ia.atom)
        e_pdb = pa.elem.strip().upper()[:1]
        if e_itp and e_pdb and e_itp != e_pdb:
            element_mismatches.append((i, ia.atom, ia.atype, pa.name, pa.elem))
    if element_mismatches and not args.allow_element_mismatch:
        sample = element_mismatches[:10]
        raise ValueError(
            "PDB/ITP atom order appears inconsistent. First mismatches: "
            + repr(sample)
            + ". Use a PDB in the same atom order as the ITP, or set --allow-element-mismatch only if you know what you are doing."
        )

    nmono = len(itp_atoms)
    nrep = args.n
    if nrep < 1:
        raise ValueError("--n must be >= 1")

    adj = build_adjacency(nmono, bonds_items)

    start_idx = resolve_atom(args.start, itp_atoms)
    end_idx = resolve_atom(args.end, itp_atoms)

    # Direction: start -> end defines the straight repeat axis.
    core_vec = vsub(pdb_atoms[end_idx].coord, pdb_atoms[start_idx].coord)
    axis = normalize(core_vec)
    translation = vadd(core_vec, vscale(axis, args.link_length))

    start_h = resolve_atom(args.start_h, itp_atoms) if args.start_h else None
    end_h = resolve_atom(args.end_h, itp_atoms) if args.end_h else None

    if start_h is None:
        # Removed from start atom when it links backward to previous repeat.
        start_h = pick_hydrogen_by_direction(start_idx, attached_hydrogens(start_idx, adj, itp_atoms), vscale(axis, -1.0), pdb_atoms)
    if end_h is None:
        # Removed from end atom when it links forward to next repeat.
        end_h = pick_hydrogen_by_direction(end_idx, attached_hydrogens(end_idx, adj, itp_atoms), axis, pdb_atoms)

    if start_h not in adj[start_idx]:
        raise ValueError(f"--start-h atom {start_h+1} is not bonded to --start atom {start_idx+1}")
    if end_h not in adj[end_idx]:
        raise ValueError(f"--end-h atom {end_h+1} is not bonded to --end atom {end_idx+1}")

    old_to_new: Dict[Tuple[int, int], int] = {}
    poly_atoms: List[PolyAtom] = []

    for r in range(nrep):
        shift = vscale(translation, r)
        for i, ia in enumerate(itp_atoms):
            if r < nrep - 1 and i == end_h:
                continue
            if r > 0 and i == start_h:
                continue

            charge = ia.charge
            # Redistribute removed H charge onto parent atom.
            if i == end_idx and r < nrep - 1:
                charge += itp_atoms[end_h].charge
            if i == start_idx and r > 0:
                charge += itp_atoms[start_h].charge

            new_idx = len(poly_atoms)
            old_to_new[(r, i)] = new_idx
            poly_atoms.append(
                PolyAtom(
                    repeat=r + 1,
                    old_index=i,
                    atype=ia.atype,
                    resid=args.resid or ia.resid,
                    atom=ia.atom,
                    elem=atom_element_from_name(ia.atom),
                    charge=charge,
                    mass=ia.mass,
                    coord=vadd(pdb_atoms[i].coord, shift),
                )
            )

    # Replicate topology sections.
    all_bonds = replicate_items(bonds_items, 2, nrep, old_to_new)
    all_pairs = replicate_items(pairs_items, 2, nrep, old_to_new)
    all_angles = replicate_items(angles_items, 3, nrep, old_to_new)
    all_dih_imp = replicate_items(dih_imp_items, 4, nrep, old_to_new)
    all_dih_prop = replicate_items(dih_prop_items, 4, nrep, old_to_new)
    all_excl = replicate_items(excl_items, 2, nrep, old_to_new)

    # Add inter-repeat bonds.
    link_bond_rest = args.link_bond.split()
    for r in range(nrep - 1):
        a = old_to_new[(r, end_idx)] + 1
        b = old_to_new[(r + 1, start_idx)] + 1
        all_bonds.append(([a, b], link_bond_rest))

    # Build polymer graph.
    n_poly = len(poly_atoms)
    bond_set: Set[Tuple[int, int]] = set()
    graph = [set() for _ in range(n_poly)]
    for idxs, _rest in all_bonds:
        i, j = idxs[0] - 1, idxs[1] - 1
        bond_set.add(canonical_pair(i, j))
        graph[i].add(j)
        graph[j].add(i)

    # Add missing cross-repeat angles generated by new bonds.
    angle_set = set(angle_key(tuple(x - 1 for x in idxs)) for idxs, _ in all_angles)
    angle_hcc = args.angle_hcc.split()
    angle_heavy = args.angle_heavy.split()

    def choose_angle_params(i: int, j: int, k: int) -> List[str]:
        # Default: H-C-C vs heavy-C-C-like. This is approximate but useful.
        if poly_atoms[i].elem == "H" or poly_atoms[k].elem == "H":
            return angle_hcc
        return angle_heavy

    added_angles = 0
    for j in range(n_poly):
        neigh = sorted(graph[j])
        for a_i in range(len(neigh)):
            for c_i in range(a_i + 1, len(neigh)):
                i = neigh[a_i]
                k = neigh[c_i]
                if len({poly_atoms[i].repeat, poly_atoms[j].repeat, poly_atoms[k].repeat}) <= 1:
                    continue
                key = angle_key((i, j, k))
                if key not in angle_set:
                    all_angles.append(([i + 1, j + 1, k + 1], choose_angle_params(i, j, k)))
                    angle_set.add(key)
                    added_angles += 1

    # Add missing cross-repeat 1-4 pairs and proper dihedrals.
    pair_set = set(canonical_pair(idxs[0] - 1, idxs[1] - 1) for idxs, _ in all_pairs)
    dih_prop_set = set(dih_key(tuple(x - 1 for x in idxs)) for idxs, _ in all_dih_prop)
    pair_rest = args.pair.split()
    dihedral_rest = args.dihedral.split()

    paths_seen: Set[Tuple[int, int, int, int]] = set()
    added_pairs = 0
    added_dihedrals = 0
    for a in range(n_poly):
        for b in graph[a]:
            for c in graph[b]:
                if c in (a, b):
                    continue
                for d in graph[c]:
                    if d in (a, b, c):
                        continue
                    key = dih_key((a, b, c, d))
                    if key in paths_seen:
                        continue
                    paths_seen.add(key)
                    if len({poly_atoms[x].repeat for x in (a, b, c, d)}) <= 1:
                        continue
                    pkey = canonical_pair(a, d)
                    if pkey not in pair_set:
                        all_pairs.append(([a + 1, d + 1], pair_rest))
                        pair_set.add(pkey)
                        added_pairs += 1
                    if key not in dih_prop_set:
                        all_dih_prop.append(([a + 1, b + 1, c + 1, d + 1], dihedral_rest))
                        dih_prop_set.add(key)
                        added_dihedrals += 1

    # Write output files.
    out_pdb = args.out_prefix + ".pdb"
    out_itp = args.out_prefix + ".itp"
    out_report = args.out_prefix + "_validation.txt"
    molname = args.molname or ((itp_atoms[0].resid if itp_atoms else "POL") + str(nrep))

    write_pdb(out_pdb, poly_atoms, bond_set, args, start_idx, end_idx, start_h, end_h)
    write_itp(out_itp, poly_atoms, molname, all_bonds, all_pairs, all_angles, all_dih_imp, all_dih_prop, all_excl)

    report = validate(poly_atoms, bond_set, all_bonds, all_pairs, all_angles, all_dih_prop, out_pdb)
    report.update({
        "input_pdb": args.pdb,
        "input_itp": args.itp,
        "n_repeats": nrep,
        "start_atom": f"{start_idx+1}:{itp_atoms[start_idx].atom}",
        "end_atom": f"{end_idx+1}:{itp_atoms[end_idx].atom}",
        "start_h_removed": f"{start_h+1}:{itp_atoms[start_h].atom}",
        "end_h_removed": f"{end_h+1}:{itp_atoms[end_h].atom}",
        "link_length_A": args.link_length,
        "added_cross_angles": added_angles,
        "added_cross_pairs": added_pairs,
        "added_cross_dihedrals": added_dihedrals,
        "out_pdb": out_pdb,
        "out_itp": out_itp,
    })
    write_report(out_report, report)

    print("Built polymer successfully.")
    print(f"  PDB: {out_pdb}")
    print(f"  ITP: {out_itp}")
    print(f"  Report: {out_report}")
    print(f"  Atoms: {len(poly_atoms)}")
    print(f"  Total charge: {report['total_charge']:.6f}")
    print(f"  Cross terms added: angles={added_angles}, pairs={added_pairs}, proper_dihedrals={added_dihedrals}")
    print(f"  Static validation: bad_refs={report['bad_refs']}, duplicates={report['duplicates']}, pdb_bond_mismatch={report['pdb_bond_mismatch']}")

# -----------------------------------------------------------------------------
# Writers and validation
# -----------------------------------------------------------------------------

def write_pdb(
    path: str,
    atoms: Sequence[PolyAtom],
    bond_set: Set[Tuple[int, int]],
    args,
    start_idx: int,
    end_idx: int,
    start_h: int,
    end_h: int,
) -> None:
    adj = [[] for _ in atoms]
    for i, j in sorted(bond_set):
        adj[i].append(j + 1)
        adj[j].append(i + 1)

    with open(path, "w") as f:
        f.write("HEADER    STRAIGHT-CORE POLYMER GENERATED FROM MONOMER PDB/ITP\n")
        f.write(f"TITLE     {args.n}-REPEAT POLYMER\n")
        f.write(f"REMARK    start atom: {args.start}; end atom: {args.end}; start-H removed: {args.start_h or 'auto'}; end-H removed: {args.end_h or 'auto'}\n")
        f.write(f"REMARK    Inter-repeat bond length: {args.link_length:.4f} Angstrom.\n")
        for idx, a in enumerate(atoms, start=1):
            x, y, z = a.coord
            name = a.atom[:4]
            resname = a.resid[:4]
            elem = a.elem.rjust(2)
            f.write(f"HETATM{idx:5d} {name:>4s} {resname:<4s}{a.repeat:5d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem}\n")
        for idx, neigh in enumerate(adj, start=1):
            neigh = sorted(set(neigh))
            for s in range(0, len(neigh), 4):
                chunk = neigh[s:s+4]
                f.write("CONECT" + f"{idx:5d}" + "".join(f"{x:5d}" for x in chunk) + "\n")
        f.write("END\n")


def write_section(f, title: str, header: str, entries: Sequence[Tuple[List[int], List[str]]]) -> None:
    f.write(f"\n[ {title} ]\n")
    if header:
        f.write(header + "\n")
    for idxs, rest in entries:
        f.write("".join(f"{x:6d}" for x in idxs))
        if rest:
            f.write("    " + "    ".join(rest))
        f.write("\n")


def write_itp(
    path: str,
    atoms: Sequence[PolyAtom],
    molname: str,
    bonds,
    pairs,
    angles,
    dih_imp,
    dih_prop,
    excl,
) -> None:
    total_charge = sum(a.charge for a in atoms)
    with open(path, "w") as f:
        f.write("; Polymer topology generated by build_straight_polymer.py\n")
        f.write("; Inter-repeat bonded parameters are generic GROMOS/ATB-like terms and must be validated.\n")
        f.write(f"; Total charge: {total_charge:.6f}\n")
        f.write("\n[ moleculetype ]\n")
        f.write("; Name   nrexcl\n")
        f.write(f"{molname}     3\n")
        f.write("\n[ atoms ]\n")
        f.write(";  nr  type  resnr  resid  atom  cgnr  charge    mass\n")
        for i, a in enumerate(atoms, start=1):
            f.write(f"{i:6d} {a.atype:>6s} {a.repeat:6d} {a.resid:>6s} {a.atom:>6s} {i:6d} {a.charge:10.6f} {a.mass:10.4f}\n")

        write_section(f, "bonds", ";  ai    aj  funct   c0         c1", bonds)
        write_section(f, "pairs", ";  ai    aj  funct", pairs)
        write_section(f, "angles", ";  ai    aj    ak  funct   angle     fc", angles)
        write_section(f, "dihedrals", "; GROMOS improper dihedrals\n;  ai    aj    ak    al  funct   angle     fc", dih_imp)
        write_section(f, "dihedrals", "; Proper dihedrals, including generated cross-repeat torsions\n;  ai    aj    ak    al  funct    ph0      cp     mult", dih_prop)
        if excl:
            write_section(f, "exclusions", ";  ai    aj", excl)


def parse_output_pdb_bonds(path: str) -> Set[Tuple[int, int]]:
    _atoms, bonds = parse_pdb(path)
    return bonds


def validate(
    atoms: Sequence[PolyAtom],
    bond_set: Set[Tuple[int, int]],
    bonds,
    pairs,
    angles,
    dih_prop,
    out_pdb: str,
) -> dict:
    n = len(atoms)
    graph = [set() for _ in range(n)]
    for i, j in bond_set:
        graph[i].add(j)
        graph[j].add(i)

    def refs_bad(entries, nidx):
        bad = 0
        for idxs, _rest in entries:
            for x in idxs[:nidx]:
                if x < 1 or x > n:
                    bad += 1
        return bad

    def dup_count(entries, canonical_func):
        keys = []
        for idxs, _rest in entries:
            key = canonical_func(tuple(x - 1 for x in idxs))
            keys.append(key)
        c = Counter(keys)
        return sum(1 for _k, v in c.items() if v > 1)

    bad_refs = {
        "bonds": refs_bad(bonds, 2),
        "pairs": refs_bad(pairs, 2),
        "angles": refs_bad(angles, 3),
        "dih_prop": refs_bad(dih_prop, 4),
    }
    duplicates = {
        "bonds": dup_count(bonds, lambda t: canonical_pair(t[0], t[1])),
        "pairs": dup_count(pairs, lambda t: canonical_pair(t[0], t[1])),
        "angles": dup_count(angles, angle_key),
        "dih_prop": dup_count(dih_prop, dih_key),
    }

    pdb_bonds = parse_output_pdb_bonds(out_pdb)
    pdb_bonds0 = set((i - 1, j - 1) for i, j in pdb_bonds)
    pdb_bond_mismatch = {
        "missing_in_pdb": len(bond_set - pdb_bonds0),
        "extra_in_pdb": len(pdb_bonds0 - bond_set),
    }

    # Cross-repeat missing terms.
    angle_set = set(angle_key(tuple(x - 1 for x in idxs)) for idxs, _ in angles)
    pair_set = set(canonical_pair(idxs[0] - 1, idxs[1] - 1) for idxs, _ in pairs)
    dih_set = set(dih_key(tuple(x - 1 for x in idxs)) for idxs, _ in dih_prop)

    missing_cross_angles = 0
    for j in range(n):
        neigh = sorted(graph[j])
        for a_i in range(len(neigh)):
            for c_i in range(a_i + 1, len(neigh)):
                i = neigh[a_i]
                k = neigh[c_i]
                if len({atoms[i].repeat, atoms[j].repeat, atoms[k].repeat}) > 1 and angle_key((i, j, k)) not in angle_set:
                    missing_cross_angles += 1

    missing_cross_pairs = 0
    missing_cross_dihedrals = 0
    seen_paths = set()
    for a in range(n):
        for b in graph[a]:
            for c in graph[b]:
                if c in (a, b):
                    continue
                for d in graph[c]:
                    if d in (a, b, c):
                        continue
                    key = dih_key((a, b, c, d))
                    if key in seen_paths:
                        continue
                    seen_paths.add(key)
                    if len({atoms[x].repeat for x in (a, b, c, d)}) <= 1:
                        continue
                    if canonical_pair(a, d) not in pair_set:
                        missing_cross_pairs += 1
                    if key not in dih_set:
                        missing_cross_dihedrals += 1

    # Simple severe clash screen, excluding bonded and 1-3 pairs. O(N^2), acceptable for small/medium polymers.
    excl = set(bond_set)
    for i in range(n):
        for j in graph[i]:
            for k in graph[j]:
                if k != i:
                    excl.add(canonical_pair(i, k))

    severe_clashes = 0
    if n <= 8000:
        for i in range(n):
            ci = atoms[i].coord
            ei = atoms[i].elem
            for j in range(i + 1, n):
                if (i, j) in excl:
                    continue
                ej = atoms[j].elem
                threshold = 0.75 if ei == "H" and ej == "H" else (0.90 if (ei == "H" or ej == "H") else 1.20)
                if dist(ci, atoms[j].coord) < threshold:
                    severe_clashes += 1

    return {
        "atoms": n,
        "total_charge": sum(a.charge for a in atoms),
        "bad_refs": bad_refs,
        "duplicates": duplicates,
        "pdb_bond_mismatch": pdb_bond_mismatch,
        "missing_cross_angles": missing_cross_angles,
        "missing_cross_pairs": missing_cross_pairs,
        "missing_cross_dihedrals": missing_cross_dihedrals,
        "severe_clashes": severe_clashes,
    }


def write_report(path: str, report: dict) -> None:
    with open(path, "w") as f:
        f.write("Straight-core polymer build validation report\n")
        f.write("============================================\n\n")
        for key, value in report.items():
            f.write(f"{key}: {value}\n")
        f.write("\nNotes:\n")
        f.write("- This is a static validation, not a replacement for gmx grompp.\n")
        f.write("- Run gmx grompp, energy minimization, and a short test MD before production.\n")

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="Build an N-repeat straight-core polymer PDB and matching GROMACS ITP from a monomer PDB/ITP."
    )
    p.add_argument("--pdb", required=True, help="Input monomer PDB. Must match ITP atom order.")
    p.add_argument("--itp", required=True, help="Input monomer ITP.")
    p.add_argument("--n", type=int, required=True, help="Total repeat units in polymer.")
    p.add_argument("--start", required=True, help="Start/link-back atom: atom name or 1-based index. This atom bonds to previous repeat.")
    p.add_argument("--end", required=True, help="End/link-forward atom: atom name or 1-based index. This atom bonds to next repeat.")
    p.add_argument("--start-h", default=None, help="Hydrogen on --start removed for backward links. Atom name or 1-based index. Default: auto.")
    p.add_argument("--end-h", default=None, help="Hydrogen on --end removed for forward links. Atom name or 1-based index. Default: auto.")
    p.add_argument("--link-length", type=float, default=1.54, help="Inter-repeat C-C link length in Angstrom for PDB coordinates. Default: 1.54.")
    p.add_argument("--out-prefix", required=True, help="Output prefix; writes .pdb, .itp, and _validation.txt.")
    p.add_argument("--molname", default=None, help="Moleculetype name in output ITP. Default: resid + N.")
    p.add_argument("--resid", default=None, help="Residue name in output PDB/ITP. Default: original ITP residue name.")

    # Generic GROMOS/ATB-like default parameters. Users can override if they know better parameters.
    p.add_argument("--link-bond", default="2 0.1530 7.1500e+06", help="Inter-repeat bond parameter fields after ai aj. Default: GROMOS-like C-C single bond.")
    p.add_argument("--pair", default="1", help="1-4 pair parameter fields after ai aj. Default: 1.")
    p.add_argument("--angle-hcc", default="2 111.00 530.00", help="Angle parameter for H-C-C-like cross angles.")
    p.add_argument("--angle-heavy", default="2 109.50 520.00", help="Angle parameter for heavy-C-C-like cross angles.")
    p.add_argument("--dihedral", default="1 0.00 5.92 3", help="Proper dihedral parameter for generated cross-repeat torsions.")

    p.add_argument("--allow-element-mismatch", action="store_true", help="Do not abort if PDB/ITP atom elements appear mismatched. Dangerous.")

    args = p.parse_args()
    build_polymer(args)


if __name__ == "__main__":
    main()
