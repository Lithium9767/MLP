#!/usr/bin/env python3
"""Map PDB 7R1C residues and secondary-structure annotations to PF00741 states."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

if __package__:
    from .hmm_coordinates import parse_alignment_columns, scan_sequences, sha256_file
else:  # Support `python bioinformatics/structure_mapping.py` from the repository root.
    from hmm_coordinates import parse_alignment_columns, scan_sequences, sha256_file


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass(frozen=True)
class PDBResidue:
    sequence_position: int
    residue_number: int
    insertion_code: str
    amino_acid: str
    secondary_structure: str


def _slice_int(line: str, start: int, end: int) -> int:
    return int(line[start:end].strip())


def parse_pdb_chain(path: Path, chain: str) -> list[PDBResidue]:
    if len(chain) != 1:
        raise ValueError("PDB chain must be one character")
    lines = path.read_text(encoding="latin-1").splitlines()
    helix_ranges: list[tuple[int, int]] = []
    sheet_ranges: list[tuple[int, int]] = []
    for line in lines:
        if line.startswith("HELIX") and len(line) >= 38 and line[19] == chain and line[31] == chain:
            helix_ranges.append((_slice_int(line, 21, 25), _slice_int(line, 33, 37)))
        elif line.startswith("SHEET") and len(line) >= 38 and line[21] == chain and line[32] == chain:
            sheet_ranges.append((_slice_int(line, 22, 26), _slice_int(line, 33, 37)))

    raw: list[tuple[int, str, str]] = []
    seen: set[tuple[int, str]] = set()
    for line in lines:
        if not line.startswith("ATOM") or len(line) < 27:
            continue
        if line[12:16].strip() != "CA" or line[21] != chain:
            continue
        residue_number = _slice_int(line, 22, 26)
        insertion_code = line[26].strip()
        key = (residue_number, insertion_code)
        if key in seen:
            continue
        seen.add(key)
        amino_acid = AA3.get(line[17:20].strip().upper(), "X")
        raw.append((residue_number, insertion_code, amino_acid))
    if not raw:
        raise ValueError(f"No CA residues found for chain {chain!r}")

    residues: list[PDBResidue] = []
    for sequence_position, (number, insertion_code, amino_acid) in enumerate(raw, start=1):
        secondary = "C"
        if any(start <= number <= end for start, end in sheet_ranges):
            secondary = "E"
        if any(start <= number <= end for start, end in helix_ranges):
            secondary = "H"
        residues.append(PDBResidue(sequence_position, number, insertion_code, amino_acid, secondary))
    return residues


def structure_sequence(residues: list[PDBResidue]) -> str:
    return "".join(residue.amino_acid for residue in residues)


def build_structure_hmm_rows(
    residues: list[PDBResidue], coordinate_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    by_position = {residue.sequence_position: residue for residue in residues}
    rows = []
    for coordinate in coordinate_rows:
        residue = by_position.get(int(coordinate["raw_position"]))
        if residue is None:
            raise ValueError("HMM coordinate references a missing structure sequence position")
        rows.append({
            "pdb_sequence_position": residue.sequence_position,
            "pdb_residue_number": residue.residue_number,
            "pdb_insertion_code": residue.insertion_code,
            "amino_acid": residue.amino_acid,
            "secondary_structure": residue.secondary_structure,
            "hmm_match_state": coordinate["hmm_match_state"],
            "is_insertion": coordinate["is_insertion"],
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--chain", default="N")
    parser.add_argument("--hmm", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=0)
    args = parser.parse_args()
    residues = parse_pdb_chain(args.pdb, args.chain)
    sequence = structure_sequence(residues)
    record = {
        "internal_id": f"PDB_{args.pdb.stem}_{args.chain}",
        "sequence": sequence,
    }
    hits, hmm_info = scan_sequences(args.hmm, [record], threads=args.threads)
    hit = hits.get(record["internal_id"])
    if hit is None:
        raise RuntimeError("No PF00741 HMM hit found for the selected PDB chain")
    coordinates = parse_alignment_columns(
        hit["hmm_sequence"], hit["target_sequence"], hit["hmm_from"], hit["target_from"]
    )
    rows = build_structure_hmm_rows(residues, coordinates)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "pdb_7r1c_pf00741_map.csv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "pdb_file": args.pdb.name,
        "pdb_sha256": sha256_file(args.pdb),
        "chain": args.chain,
        "pdb_residues": len(residues),
        "mapped_residues": len(rows),
        "hmm": hmm_info,
        "hmm_sha256": sha256_file(args.hmm),
        "mapping_sha256": sha256_file(output_path),
    }
    (args.output_dir / "pdb_7r1c_mapping_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
