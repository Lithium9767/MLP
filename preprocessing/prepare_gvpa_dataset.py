#!/usr/bin/env python3
"""Prepare the teacher-provided GvpA FASTA as an auditable candidate dataset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


FIELDS = [
    "sequence_id",
    "accession",
    "database",
    "source_url",
    "retrieval_date",
    "sequence",
    "sequence_length",
    "sequence_hash",
    "valid_residues",
    "organism",
    "taxonomy_id",
    "lineage",
    "pfam_start",
    "pfam_end",
    "label",
    "label_definition",
    "label_source",
    "evidence_level",
    "duplicate_cluster",
    "similarity_cluster",
    "split",
    "exclusion_reason",
]

ALLOWED_AA = frozenset("ACDEFGHIKLMNPQRSTVWYXBZJUO")
ORGANISM_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts).upper()))
            header = line[1:].strip()
            parts = []
        elif header is None:
            raise ValueError(f"Sequence before first header at {path}:{line_number}")
        else:
            parts.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(parts).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def build_rows(path: Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    raw_records = parse_fasta(path)
    rows: list[dict[str, str]] = []
    seen_hashes: dict[str, str] = {}
    accessions: Counter[str] = Counter()

    for index, (header, sequence) in enumerate(raw_records, start=1):
        accession = header.split(maxsplit=1)[0] if header else ""
        accessions[accession] += 1
        digest = hashlib.sha256(sequence.encode("ascii", errors="ignore")).hexdigest()
        invalid = sorted(set(sequence) - ALLOWED_AA)
        reasons = []
        if not accession:
            reasons.append("missing_accession")
        if not sequence:
            reasons.append("empty_sequence")
        if invalid:
            reasons.append("invalid_residues:" + "".join(invalid))
        if digest in seen_hashes:
            reasons.append("exact_duplicate_of:" + seen_hashes[digest])
        else:
            seen_hashes[digest] = accession or f"record_{index}"

        organism_match = ORGANISM_RE.search(header)
        organism = organism_match.group(1).strip() if organism_match else ""
        rows.append(
            {
                "sequence_id": accession or f"record_{index}",
                "accession": accession,
                "database": "NCBI Protein (inferred from supplied accession)",
                "source_url": (
                    f"https://www.ncbi.nlm.nih.gov/protein/{accession}" if accession else ""
                ),
                "retrieval_date": "",
                "sequence": sequence,
                "sequence_length": str(len(sequence)),
                "sequence_hash": digest,
                "valid_residues": str(bool(sequence) and not invalid).lower(),
                "organism": organism,
                "taxonomy_id": "",
                "lineage": "",
                "pfam_start": "",
                "pfam_end": "",
                "label": "",
                "label_definition": "pending_prediction_target_confirmation",
                "label_source": "",
                "evidence_level": "unlabelled_sequence",
                "duplicate_cluster": f"exact_{digest[:16]}",
                "similarity_cluster": "",
                "split": "",
                "exclusion_reason": ";".join(reasons),
            }
        )

    included = [row for row in rows if not row["exclusion_reason"]]
    lengths = [int(row["sequence_length"]) for row in included]
    summary = {
        "dataset_status": "candidate_not_frozen",
        "source_file": str(path),
        "raw_records": len(rows),
        "included_after_validation_and_exact_dedup": len(included),
        "excluded_records": len(rows) - len(included),
        "unique_accessions": len(accessions),
        "duplicate_accession_values": sum(count > 1 for count in accessions.values()),
        "missing_organism_records": sum(not row["organism"] for row in rows),
        "min_length": min(lengths) if lengths else None,
        "max_length": max(lengths) if lengths else None,
        "prediction_label_status": "missing_pending_M1_decision",
        "similarity_split_status": "pending_mmseqs2",
    }
    return rows, summary


def write_outputs(rows: list[dict[str, str]], summary: dict[str, object], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    with (output / "metadata.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    with (output / "sequences.fasta").open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            if row["exclusion_reason"]:
                continue
            handle.write(f'>{row["sequence_id"]}\n')
            sequence = row["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    (output / "audit_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("gv/gv/design/data/GvpA.fasta"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/gvpa_v0.1_candidate"),
    )
    args = parser.parse_args()
    rows, summary = build_rows(args.input)
    write_outputs(rows, summary, args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
