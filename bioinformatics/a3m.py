"""Convert a protein A3M into match columns and lossless residue coordinates.

HH-suite defines uppercase residues / '-' as match states, lowercase residues
as insertions, and optional '.' as insert gaps. Case must survive parsing.
Reference: https://github.com/soedinglab/hh-suite/blob/master/scripts/reformat.pl

``full_sequences`` reconstructs every residue PRESENT IN THE A3M, including
insertions. It is only a verified original sequence when --original succeeds;
an A3M alone cannot establish that external source sequences were not trimmed.

All residue/match coordinates are 1-based. Insertions have no match position.
Their anchor is the number of preceding match columns (0 before the first;
alignment length after the last). Their rank counts residues within that
sequence's insertion, ignoring '.'. Equal insertion ranks across sequences
do not establish homology. Gaps generate no residue rows.

The match FASTA is compatible with bioinformatics.msa_analysis. That tool's
residue_position is the UNGAPPED MATCH-ONLY index, not the original position.
Join it to this module's match_residue_position or use match_position instead.
For downstream compatibility this module rejects records with no match residue,
even if their insertions contain residues. This is a project QC restriction.

Example (the supplied fixture is synthetic, not a biological result)::

    python -m bioinformatics.a3m --input bioinformatics/examples/synthetic.a3m \
        --original bioinformatics/examples/synthetic_a3m_original.fasta \
        --out-dir results/synthetic_a3m --purpose synthetic
"""

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

from bioinformatics.msa_analysis import (
    AMBIGUOUS, CANONICAL, fingerprint, git_state, read_fasta, write_csv,
)


RESIDUES = CANONICAL | AMBIGUOUS
ALLOWED = RESIDUES | frozenset(residue.lower() for residue in RESIDUES) | {"-", "."}
ANNOTATION_PREFIXES = ("ss_", "sa_", "aa_")
MAP_FIELDS = [
    "sequence_id", "original_position", "residue", "state", "match_position",
    "match_residue_position", "insertion_anchor", "insertion_rank",
]
QC_FIELDS = [
    "sequence_id", "match_columns", "full_length", "match_residues",
    "inserted_residues", "match_gaps", "insert_gaps", "ambiguous_residues",
]


@dataclass
class A3MAlignment:
    """Validated protein records; annotation records are excluded and audited."""

    sequences: dict
    headers: dict
    excluded_records: list
    comments: list

    @property
    def match_alignment(self):
        return {key: "".join(r for r in value if r in RESIDUES or r == "-")
                for key, value in self.sequences.items()}

    @property
    def full_sequences(self):
        """Reconstruct match + inserted residues, never match residues alone."""
        return {key: value.replace("-", "").replace(".", "").upper()
                for key, value in self.sequences.items()}

    @property
    def match_columns(self):
        return len(next(iter(self.match_alignment.values())))


def read_a3m(path):
    """Parse wrapped A3M, preserving case and rejecting ambiguous input.

    First header token is the stable ID. Duplicate IDs, including annotation
    IDs, are rejected. Reserved HH-suite ss_/sa_/aa_ records are excluded before
    residue validation; their lengths and content hashes remain in the audit.
    Full-line # comments are retained as metadata, never parsed as sequence.
    """
    records = {}
    headers = {}
    comments = []
    identifier = None
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            comments.append({"line": number, "text": line})
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            tokens = header.split()
            if not tokens:
                raise ValueError(f"line {number}: empty A3M identifier")
            identifier = tokens[0]
            if ">" in identifier or any(ord(c) < 33 or ord(c) == 127 for c in identifier):
                raise ValueError(f"line {number}: invalid A3M identifier")
            if identifier in records:
                raise ValueError(f"line {number}: duplicate identifier {identifier}")
            records[identifier] = []
            headers[identifier] = header
            continue
        if identifier is None:
            raise ValueError(f"line {number}: sequence before A3M header")
        records[identifier].append(line)

    sequences = {}
    excluded = []
    expected_columns = None
    for key, parts in records.items():
        sequence = "".join(parts)
        if key.startswith(ANNOTATION_PREFIXES):
            excluded.append({
                "sequence_id": key, "header": headers[key],
                "reason": "Reserved HH-suite annotation prefix (ss_, sa_, or aa_)",
                "characters": len(sequence),
                "content_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
            })
            continue
        if not sequence:
            raise ValueError(f"{key}: empty A3M sequence")
        invalid = set(sequence) - ALLOWED
        if invalid:
            raise ValueError(f"{key}: invalid A3M characters {sorted(invalid)}")
        matches = sum(residue in RESIDUES for residue in sequence)
        columns = matches + sequence.count("-")
        if not columns:
            raise ValueError(f"{key}: A3M has no match columns")
        if not matches:
            raise ValueError(f"{key}: no match residues; match-only alignment would be entirely gapped")
        if expected_columns is None:
            expected_columns = columns
        elif columns != expected_columns:
            raise ValueError(f"{key}: match column count {columns} differs from {expected_columns}")
        sequences[key] = sequence
    if not sequences:
        raise ValueError("A3M contains no protein sequences after annotation exclusion")
    return A3MAlignment(sequences, {key: headers[key] for key in sequences}, excluded, comments)


def validate_originals(alignment, originals):
    """Require exact IDs and ALL A3M residues, including inserted residues."""
    reconstructed = alignment.full_sequences
    if reconstructed.keys() != originals.keys():
        missing = sorted(reconstructed.keys() - originals.keys())
        extra = sorted(originals.keys() - reconstructed.keys())
        raise ValueError(f"original FASTA and A3M identifiers differ; missing={missing}, extra={extra}")
    for key, sequence in reconstructed.items():
        if sequence != originals[key]:
            raise ValueError(f"{key}: full A3M sequence (including inserts) does not reconstruct original")


def coordinate_rows(alignment):
    """One row per actual residue, with explicit match or insertion coordinates."""
    for key, sequence in alignment.sequences.items():
        original_position = match_position = match_residue_position = insertion_rank = 0
        for residue in sequence:
            if residue == ".":
                continue
            if residue == "-":
                match_position += 1
                insertion_rank = 0
                continue
            original_position += 1
            is_match = residue in RESIDUES
            if is_match:
                match_position += 1
                match_residue_position += 1
                insertion_rank = 0
            else:
                insertion_rank += 1
            yield {
                "sequence_id": key, "original_position": original_position,
                "residue": residue.upper(), "state": "match" if is_match else "insert",
                "match_position": match_position if is_match else None,
                "match_residue_position": match_residue_position if is_match else None,
                "insertion_anchor": None if is_match else match_position,
                "insertion_rank": None if is_match else insertion_rank,
            }


def sequence_qc(alignment):
    for key, sequence in alignment.sequences.items():
        matches = sum(residue in RESIDUES for residue in sequence)
        inserts = sum(residue.islower() for residue in sequence)
        yield {
            "sequence_id": key, "match_columns": matches + sequence.count("-"),
            "full_length": matches + inserts, "match_residues": matches,
            "inserted_residues": inserts, "match_gaps": sequence.count("-"),
            "insert_gaps": sequence.count("."),
            "ambiguous_residues": sum(residue.upper() in AMBIGUOUS for residue in sequence),
        }


def write_fasta(path, sequences):
    with Path(path).open("w", encoding="utf-8", newline="\n") as stream:
        for identifier, sequence in sequences.items():
            stream.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                stream.write(sequence[start:start + 80] + "\n")


def run(args):
    """Validate all inputs before creating the new, never-overwritten output dir."""
    if args.purpose not in ("synthetic", "exploratory"):
        raise ValueError("purpose must be synthetic or exploratory")
    alignment = read_a3m(args.input)
    if args.original is not None:
        validate_originals(alignment, read_fasta(args.original, aligned=False))
    qc = list(sequence_qc(alignment))
    destination = Path(args.out_dir)
    # Hash provenance before writing. Inputs are never modified by this tool.
    input_info = fingerprint(args.input)
    original_info = fingerprint(args.original) if args.original is not None else None
    implementation = fingerprint(__file__)
    destination.mkdir(parents=True, exist_ok=False)
    try:
        write_fasta(destination / "match_alignment.fasta", alignment.match_alignment)
        write_fasta(destination / "full_sequences.fasta", alignment.full_sequences)
        write_csv(destination / "coordinate_map.csv", coordinate_rows(alignment), MAP_FIELDS)
        write_csv(destination / "sequence_qc.csv", qc, QC_FIELDS)
        manifest = {
            "schema_version": "1.0", "created_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": args.purpose, "synthetic": args.purpose == "synthetic",
            "synthetic_flag_source": "Caller-declared --purpose; not inferred biological evidence",
            "input": input_info, "original": original_info,
            "original_verified": args.original is not None,
            "parameters": {"input": str(args.input), "original": str(args.original) if args.original else None,
                           "out_dir": str(args.out_dir), "purpose": args.purpose},
            "n_sequences": len(alignment.sequences), "match_columns": alignment.match_columns,
            "n_full_residues": sum(row["full_length"] for row in qc),
            "n_inserted_residues": sum(row["inserted_residues"] for row in qc),
            "n_excluded_records": len(alignment.excluded_records),
            "excluded_records": alignment.excluded_records, "input_comments": alignment.comments,
            "source_headers": alignment.headers,
            "coordinate_convention": {
                "original_position": "1-based in full A3M reconstruction, including lowercase inserts",
                "match_position": "1-based match column; absent for inserted residues",
                "match_residue_position": "1-based ungapped match-only residue; matches msa_analysis residue_position",
                "insertion_anchor": "Number of preceding match columns, including '-' (0 to match_columns)",
                "insertion_rank": "1-based actual residue within this sequence's insertion; ignores '.'",
                "gaps": "No residue rows for '-' or '.'",
            },
            "python": platform.python_version(), "git": git_state(),
            "implementation": implementation, "command_argv": sys.argv,
            "outputs": {p.name: fingerprint(p)["sha256"] for p in sorted(destination.iterdir())},
            "limitations": [
                "Full sequences reconstruct A3M residues; untrimmed source sequences are verified only with --original",
                "Match-only FASTA omits inserts; its ungapped residue numbering is not original_position",
                "Matching insertion anchor/rank across sequences does not establish inserted-residue homology",
                "ss_, sa_, aa_ identifiers are reserved annotation records and are excluded",
                "Records without any uppercase match residue are rejected for downstream analysis compatibility",
                "Format validation and descriptive QC do not provide functional labels or biological conclusions",
            ],
        }
        (destination / "manifest.json").write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    except Exception:
        (destination / "FAILED.txt").write_text(
            "Run incomplete. Do not use partial outputs. Rerun into a new directory.\n", encoding="utf-8")
        raise
    return manifest


def main():
    parser = argparse.ArgumentParser(description="Convert protein A3M without losing insertion coordinates")
    parser.add_argument("--input", required=True, type=Path, help="A3M; letter case defines match/insert states")
    parser.add_argument("--original", type=Path, help="Optional ungapped FASTA with exactly matching sequence IDs")
    parser.add_argument("--out-dir", required=True, type=Path, help="New directory; never overwritten")
    parser.add_argument("--purpose", required=True, choices=["synthetic", "exploratory"])
    args = parser.parse_args()
    try:
        manifest = run(args)
    except (ValueError, OSError) as error:
        parser.exit(2, f"error: {error}\n")
    print(f"Completed: {manifest['n_sequences']} sequences, {manifest['match_columns']} match columns, "
          f"{manifest['n_inserted_residues']} inserted residues, "
          f"{manifest['n_excluded_records']} excluded annotations -> {args.out_dir}")


if __name__ == "__main__":
    main()
