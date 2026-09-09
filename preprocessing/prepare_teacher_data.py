#!/usr/bin/env python3
"""Convert the teacher-provided recognition JSON files into an auditable dataset.

This script deliberately creates a *candidate* dataset.  The GVP subtype is copied
from the provided directory/annotation and must not be treated as an independently
validated functional label until A/B/C finish the M1 label audit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


METADATA_FIELDS = [
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

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = frozenset("XBZJUO")
ACCESSION_RE = re.compile(r"^[A-Za-z]{1,4}_?[A-Za-z0-9]+(?:\.\d+)?$")


def normalise_sequence(value: Any) -> str:
    """Return an upper-case sequence with all whitespace removed."""
    return "".join(str(value or "").split()).upper()


def canonical_label(directory_name: str) -> str:
    """Convert a directory such as gvpa to the report label GvpA."""
    value = directory_name.strip().lower()
    if not re.fullmatch(r"gvp[a-z]", value):
        raise ValueError(f"Cannot infer GVP label from directory: {directory_name}")
    return "Gvp" + value[-1].upper()


def source_url(accession: str) -> str:
    if accession and ACCESSION_RE.fullmatch(accession):
        return f"https://www.ncbi.nlm.nih.gov/protein/{accession}"
    return ""


def load_records(source_dir: Path) -> list[dict[str, str]]:
    json_paths = sorted(source_dir.glob("*/*.json"))
    if not json_paths:
        raise FileNotFoundError(f"No */*.json files found under {source_dir}")

    rows: list[dict[str, str]] = []
    for json_path in json_paths:
        label = canonical_label(json_path.parent.name)
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        if not isinstance(payload, list):
            raise ValueError(f"Expected a JSON list in {json_path}")

        for index, item in enumerate(payload):
            if not isinstance(item, dict):
                raise ValueError(f"Record {index} in {json_path} is not an object")

            accession = str(item.get("unique_sequence_id") or "").strip()
            sequence = normalise_sequence(item.get("sequence"))
            annotation = item.get("representative_annotation") or {}
            if not isinstance(annotation, dict):
                annotation = {}

            invalid = sorted(set(sequence) - STANDARD_AA - AMBIGUOUS_AA)
            reasons = []
            if not accession:
                reasons.append("missing_accession")
            if not sequence:
                reasons.append("empty_sequence")
            if invalid:
                reasons.append("invalid_residues:" + "".join(invalid))

            digest = hashlib.sha256(sequence.encode("ascii", errors="ignore")).hexdigest()
            rows.append(
                {
                    "sequence_id": accession or f"{label}_{index + 1}",
                    "accession": accession,
                    "database": "NCBI Protein (inferred from supplied accession)",
                    "source_url": source_url(accession),
                    "retrieval_date": "",
                    "sequence": sequence,
                    "sequence_length": str(len(sequence)),
                    "sequence_hash": digest,
                    "valid_residues": str(not invalid and bool(sequence)).lower(),
                    "organism": str(annotation.get("organism") or "").strip(),
                    "taxonomy_id": "",
                    "lineage": "",
                    "pfam_start": "",
                    "pfam_end": "",
                    "label": label,
                    "label_definition": "GVP subtype supplied by directory and representative annotation",
                    "label_source": str(json_path),
                    "evidence_level": "database_annotation_unverified",
                    "duplicate_cluster": f"exact_{digest[:16]}",
                    "similarity_cluster": "",
                    "split": "",
                    "exclusion_reason": ";".join(reasons),
                }
            )
    return rows


def audit(rows: list[dict[str, str]]) -> dict[str, Any]:
    label_counts = Counter(row["label"] for row in rows)
    lengths: dict[str, list[int]] = defaultdict(list)
    sequence_labels: dict[str, set[str]] = defaultdict(set)
    accession_labels: dict[str, set[str]] = defaultdict(set)
    invalid_count = 0
    missing_organism_count = 0

    for row in rows:
        lengths[row["label"]].append(int(row["sequence_length"]))
        sequence_labels[row["sequence_hash"]].add(row["label"])
        if row["accession"]:
            accession_labels[row["accession"]].add(row["label"])
        invalid_count += bool(row["exclusion_reason"])
        missing_organism_count += not bool(row["organism"])

    per_label = {}
    for label in sorted(label_counts):
        values = sorted(lengths[label])
        middle = len(values) // 2
        median = (
            values[middle]
            if len(values) % 2
            else (values[middle - 1] + values[middle]) / 2
        )
        per_label[label] = {
            "count": label_counts[label],
            "min_length": values[0],
            "median_length": median,
            "max_length": values[-1],
        }

    counts = list(label_counts.values())
    return {
        "dataset_status": "candidate_not_frozen",
        "label_evidence": "database_annotation_unverified",
        "total_records": len(rows),
        "unique_sequence_hashes": len(sequence_labels),
        "unique_accessions": len({row["accession"] for row in rows if row["accession"]}),
        "invalid_or_incomplete_records": invalid_count,
        "missing_organism_records": missing_organism_count,
        "cross_label_exact_sequence_conflicts": sum(
            len(labels) > 1 for labels in sequence_labels.values()
        ),
        "cross_label_accession_conflicts": sum(
            len(labels) > 1 for labels in accession_labels.values()
        ),
        "class_imbalance_ratio_max_to_min": round(max(counts) / min(counts), 3),
        "per_label": per_label,
        "blocking_items": [
            "Confirm that GVP-subtype recognition is the intended prediction task.",
            "Trace the original database query, retrieval date, version, and licence.",
            "Independently review label evidence; unknown annotations are not negatives.",
            "Create sequence-similarity clusters before train/validation/test splitting.",
            "Keep Pfam validation independent from label construction to avoid circular reasoning.",
        ],
    }


def write_outputs(rows: list[dict[str, str]], output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = output_dir / "metadata.csv"
    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    fasta_path = output_dir / "sequences.fasta"
    with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(f'>{row["sequence_id"]}|label={row["label"]}\n')
            sequence = row["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    summary = audit(rows)
    (output_dir / "audit_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("gv/gv/recognition/data/gvp"),
        help="Directory containing one JSON file per GVP subtype directory.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/teacher_recognition_v0.1_candidate"),
        help="Output directory (normally ignored by Git).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = load_records(args.source)
    summary = write_outputs(rows, args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
