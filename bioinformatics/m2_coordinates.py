"""Small reviewed HMM alignment and PDB coordinate readers for the M2 run.

The alignment and PDB parsing approach originated in PR #4; these helpers
exclude its unvalidated full-dataset scan and cross-split conservation CLI.
"""

from dataclasses import dataclass
from pathlib import Path


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_alignment_columns(hmm_sequence, target_sequence, hmm_from: int,
                            target_from: int) -> list[dict]:
    """Map HMM alignment columns to full target residues; inserts have no state."""
    hmm_sequence = hmm_sequence.decode("ascii") if isinstance(hmm_sequence, bytes) else str(hmm_sequence)
    target_sequence = (target_sequence.decode("ascii") if isinstance(target_sequence, bytes)
                       else str(target_sequence))
    if len(hmm_sequence) != len(target_sequence) or hmm_from < 1 or target_from < 1:
        raise ValueError("Invalid HMM alignment length or 1-based origin")
    state, position = hmm_from, target_from
    rows = []
    for hmm_char, target_char in zip(hmm_sequence, target_sequence):
        insertion = hmm_char == "."
        if not (hmm_char.isalpha() or insertion) or not (target_char.isalpha() or target_char in ".-"):
            raise ValueError("Invalid HMM alignment symbol")
        if target_char not in ".-":
            rows.append({"raw_position": position,
                         "hmm_match_state": None if insertion else state,
                         "is_insertion": insertion,
                         "residue": target_char.upper()})
            position += 1
        if not insertion:
            state += 1
    return rows


@dataclass(frozen=True)
class PDBResidue:
    residue_number: int
    insertion_code: str
    amino_acid: str
    secondary_structure: str


def parse_pdb_chain(path: Path, chain: str) -> list[PDBResidue]:
    """Read modeled CA residues and PDB HELIX/SHEET intervals for author chain."""
    if len(chain) != 1:
        raise ValueError("PDB chain must be one character")
    lines = path.read_text(encoding="latin-1").splitlines()
    helices = []
    sheets = []
    for line in lines:
        if line.startswith("HELIX") and len(line) >= 38 and line[19] == chain and line[31] == chain:
            helices.append((int(line[21:25]), int(line[33:37])))
        elif line.startswith("SHEET") and len(line) >= 38 and line[21] == chain and line[32] == chain:
            sheets.append((int(line[22:26]), int(line[33:37])))
    seen = set()
    rows = []
    for line in lines:
        if not line.startswith("ATOM") or len(line) < 27 or line[12:16].strip() != "CA" or line[21] != chain:
            continue
        if line[16] not in " A":
            continue
        number, insertion_code = int(line[22:26]), line[26].strip()
        if (number, insertion_code) in seen:
            continue
        seen.add((number, insertion_code))
        amino_acid = AA3.get(line[17:20].strip().upper(), "X")
        secondary = "H" if any(start <= number <= end for start, end in helices) else (
            "E" if any(start <= number <= end for start, end in sheets) else "C")
        rows.append(PDBResidue(number, insertion_code, amino_acid, secondary))
    if not rows:
        raise ValueError(f"No modeled CA residues in PDB chain {chain!r}")
    return rows
