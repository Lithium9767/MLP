#!/usr/bin/env python3
"""Check M2 data-audit, MMseqs2, and homology-split receipts.

Exit status is 0 only when the committed receipts are complete. A missing
receipt, a missing SHA-256, or homology-cluster leakage fails the run.
This script does not replace the frozen split and does not read validation
sequences for candidate discovery.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_AUDIT = ROOT / "results" / "data_audit" / "gvpa_v1_audit_summary.json"
DEFAULT_MMSEQS = ROOT / "results" / "data_audit" / "gvpa_v1_mmseqs_run.json"
DEFAULT_SPLIT = ROOT / "results" / "data_audit" / "gvpa_v1_split_summary.json"

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

AUDIT_FIELDS = (
    "dataset_version",
    "candidate_records",
    "sequence_qc_eligible",
    "primary_cohort",
    "unique_internal_ids",
    "output_sha256",
    "source",
)
AUDIT_HASHES = (
    ("output_sha256", "metadata.csv"),
    ("output_sha256", "members.csv"),
    ("output_sha256", "sequences_for_clustering.fasta"),
    ("source", "source_member_sha256"),
    ("source", "source_zip_sha256"),
)
MMSEQS_FIELDS = (
    "status",
    "mmseqs_version",
    "input_sequences",
    "output_clusters",
    "min_sequence_identity",
    "coverage",
    "coverage_mode",
    "input_fasta_sha256",
    "mmseqs_cluster_tsv_sha256",
)
MMSEQS_HASHES = (
    ("input_fasta_sha256",),
    ("mmseqs_cluster_tsv_sha256",),
)
SPLIT_FIELDS = (
    "split_version",
    "schema_version",
    "seed",
    "cluster_leakage",
    "total_sequences",
    "total_clusters",
    "sequence_counts",
    "cluster_counts",
    "metadata_sha256",
    "mmseqs_cluster_tsv_sha256",
    "split_manifest_sha256",
)
SPLIT_HASHES = (
    ("metadata_sha256",),
    ("mmseqs_cluster_tsv_sha256",),
    ("split_manifest_sha256",),
)
FROZEN = {
    "dataset_version": "gvpa-recognition-c7f6f005d717",
    "split_version": "homology-8b9005e2d9-s42",
    "candidate_records": 2078,
    "unique_internal_ids": 2078,
    "sequence_qc_eligible": 2076,
    "primary_cohort": 1721,
    "input_sequences": 2076,
    "output_clusters": 478,
    "total_sequences": 2076,
    "total_clusters": 478,
    "discovery_sequences": 1453,
    "validation_sequences": 623,
    "discovery_clusters": 335,
    "validation_clusters": 143,
}


def load_receipt(path: Path) -> tuple[dict[str, Any] | None, list[str]]:
    if not path.is_file():
        return None, [f"missing receipt: {path}"]
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return None, [f"unreadable receipt: {path}: {exc}"]
    if not isinstance(payload, dict):
        return None, [f"missing receipt: {path} is not a JSON object"]
    return payload, []


def missing_fields(payload: dict[str, Any], fields: tuple[str, ...], receipt_name: str) -> list[str]:
    errors = []
    for field in fields:
        if field not in payload or payload[field] is None or payload[field] == "":
            errors.append(f"missing receipt field: {receipt_name}.{field}")
    return errors


def lookup(payload: dict[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = payload
    for part in path:
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def hash_errors(payload: dict[str, Any], paths: tuple[tuple[str, ...], ...], receipt_name: str) -> list[str]:
    errors = []
    for path in paths:
        label = f"{receipt_name}.{'.'.join(path)}"
        value = lookup(payload, path)
        if value is None or value == "":
            errors.append(f"missing hash: {label}")
        elif not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
            errors.append(f"invalid hash: {label}")
    return errors


def cluster_leakage_errors(split_receipt: dict[str, Any], manifest: Path | None) -> list[str]:
    errors = []
    if split_receipt.get("cluster_leakage") is not False:
        errors.append("cluster leakage: split receipt cluster_leakage is not false")
    if manifest is None:
        return errors
    if not manifest.is_file():
        return errors + [f"missing receipt: {manifest}"]
    try:
        leaked = leaking_clusters(manifest)
    except ValueError as exc:
        return errors + [str(exc)]
    errors.extend(
        f"cluster leakage: homology cluster {cluster_id} is assigned to more than one split"
        for cluster_id in leaked
    )
    return errors


def leaking_clusters(manifest: Path) -> list[str]:
    splits_by_cluster: dict[str, set[str]] = defaultdict(set)
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        required = {"homology_cluster", "split"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise ValueError(f"missing receipt field: split_manifest.{missing}")
        for line_number, row in enumerate(reader, start=2):
            cluster_id = (row.get("homology_cluster") or "").strip()
            split_name = (row.get("split") or "").strip()
            if not cluster_id or not split_name:
                raise ValueError(f"missing receipt field: split_manifest row {line_number}")
            splits_by_cluster[cluster_id].add(split_name)
    return sorted(cluster_id for cluster_id, names in splits_by_cluster.items() if len(names) > 1)


def agreement_errors(audit: dict[str, Any], mmseqs: dict[str, Any], split: dict[str, Any]) -> list[str]:
    errors = []
    fasta_hash = lookup(audit, ("output_sha256", "sequences_for_clustering.fasta"))
    metadata_hash = lookup(audit, ("output_sha256", "metadata.csv"))
    if isinstance(fasta_hash, str) and isinstance(mmseqs.get("input_fasta_sha256"), str):
        if fasta_hash != mmseqs["input_fasta_sha256"]:
            errors.append("hash mismatch: audit sequences_for_clustering.fasta != mmseqs input_fasta_sha256")
    if isinstance(metadata_hash, str) and isinstance(split.get("metadata_sha256"), str):
        if metadata_hash != split["metadata_sha256"]:
            errors.append("hash mismatch: audit metadata.csv != split metadata_sha256")
    mmseqs_clusters = mmseqs.get("mmseqs_cluster_tsv_sha256")
    split_clusters = split.get("mmseqs_cluster_tsv_sha256")
    if isinstance(mmseqs_clusters, str) and isinstance(split_clusters, str) and mmseqs_clusters != split_clusters:
        errors.append("hash mismatch: mmseqs cluster TSV != split mmseqs_cluster_tsv_sha256")
    return errors


def baseline_errors(audit: dict[str, Any], mmseqs: dict[str, Any], split: dict[str, Any]) -> list[str]:
    sequence_counts = split.get("sequence_counts") if isinstance(split.get("sequence_counts"), dict) else {}
    cluster_counts = split.get("cluster_counts") if isinstance(split.get("cluster_counts"), dict) else {}
    observed = {
        "dataset_version": audit.get("dataset_version"),
        "split_version": split.get("split_version"),
        "candidate_records": audit.get("candidate_records"),
        "unique_internal_ids": audit.get("unique_internal_ids"),
        "sequence_qc_eligible": audit.get("sequence_qc_eligible"),
        "primary_cohort": audit.get("primary_cohort"),
        "input_sequences": mmseqs.get("input_sequences"),
        "output_clusters": mmseqs.get("output_clusters"),
        "total_sequences": split.get("total_sequences"),
        "total_clusters": split.get("total_clusters"),
        "discovery_sequences": sequence_counts.get("discovery"),
        "validation_sequences": sequence_counts.get("validation"),
        "discovery_clusters": cluster_counts.get("discovery"),
        "validation_clusters": cluster_counts.get("validation"),
    }
    return [
        f"frozen baseline mismatch: {name} is {observed[name]!r}, expected {expected!r}"
        for name, expected in FROZEN.items()
        if observed[name] != expected
    ]


def validate_release(
    audit_path: Path,
    mmseqs_path: Path,
    split_path: Path,
    manifest_path: Path | None = None,
) -> list[str]:
    audit, audit_errors = load_receipt(audit_path)
    mmseqs, mmseqs_errors = load_receipt(mmseqs_path)
    split, split_errors = load_receipt(split_path)
    errors = audit_errors + mmseqs_errors + split_errors
    if audit is not None:
        errors.extend(missing_fields(audit, AUDIT_FIELDS, "audit"))
        errors.extend(hash_errors(audit, AUDIT_HASHES, "audit"))
    if mmseqs is not None:
        errors.extend(missing_fields(mmseqs, MMSEQS_FIELDS, "mmseqs"))
        errors.extend(hash_errors(mmseqs, MMSEQS_HASHES, "mmseqs"))
    if split is not None:
        errors.extend(missing_fields(split, SPLIT_FIELDS, "split"))
        errors.extend(hash_errors(split, SPLIT_HASHES, "split"))
        errors.extend(cluster_leakage_errors(split, manifest_path))
    elif manifest_path is not None:
        errors.extend(cluster_leakage_errors({"cluster_leakage": False}, manifest_path))
    if audit is not None and mmseqs is not None and split is not None:
        errors.extend(agreement_errors(audit, mmseqs, split))
        errors.extend(baseline_errors(audit, mmseqs, split))
    return errors


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--mmseqs", type=Path, default=DEFAULT_MMSEQS)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Optional split_manifest.csv. When given, a cluster in more than one split fails validation.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    errors = validate_release(args.audit, args.mmseqs, args.split, args.manifest)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print("M2 release validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
