#!/usr/bin/env python3
"""Build and check the frozen M3 input manifests without publishing sequences."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

DATASET_VERSION = "gvpa-recognition-c7f6f005d717"
SPLIT_VERSION = "homology-8b9005e2d9-s42"
METADATA_SHA256 = "62893e31a2dbc245f29c0e9af0a282dcb5099d0e95b96db1a9cead2ac933004e"
SPLIT_SHA256 = "0a6bd3b44d284a407e1cc9d5be758f0f9439c28395deb6d62d2c787f95341e9c"
FIELDS = ["internal_id", "sequence_id", "sequence_sha256", "sequence_length",
          "analysis_cohort", "split", "homology_cluster"]
EXPECTED = {"discovery_primary": 1202, "validation_primary": 519, "sensitivity": 355}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing columns {sorted(required - set(reader.fieldnames or []))}")
        return list(reader)


def unique_by_id(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    ids = [row["internal_id"] for row in rows]
    if not all(ids) or len(ids) != len(set(ids)):
        raise ValueError(f"{label}: missing or duplicate internal_id")
    return {row["internal_id"]: row for row in rows}


def source_rows(metadata: Path, split: Path) -> dict[str, list[dict[str, str]]]:
    if sha256(metadata) != METADATA_SHA256 or sha256(split) != SPLIT_SHA256:
        raise ValueError("M2 input hash differs from the frozen release")
    meta = unique_by_id(read_csv(metadata, {
        "internal_id", "sequence_id", "sequence", "sequence_sha256", "sequence_length",
        "sequence_qc_eligible", "primary_analysis_eligible", "analysis_cohort"}), "metadata")
    assignments = unique_by_id(read_csv(split, {
        "internal_id", "sequence_id", "analysis_cohort", "primary_analysis_eligible",
        "homology_cluster", "split"}), "split")
    eligible = {key for key, row in meta.items() if row["sequence_qc_eligible"] == "true"}
    if eligible != set(assignments) or len(eligible) != 2076:
        raise ValueError("split IDs do not cover exactly the 2076 QC eligible sequences")
    if len({row["sequence_id"] for row in meta.values()}) != len(meta):
        raise ValueError("metadata contains duplicate sequence_id")
    owners: dict[str, str] = {}
    groups: dict[str, list[dict[str, str]]] = {name: [] for name in EXPECTED}
    for key, assignment in assignments.items():
        row = meta[key]
        sequence = row["sequence"].upper()
        if hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
            raise ValueError(f"{key}: sequence hash mismatch")
        if len(sequence) != int(row["sequence_length"]):
            raise ValueError(f"{key}: sequence length mismatch")
        if assignment["sequence_id"] != row["sequence_id"] or assignment["analysis_cohort"] != row["analysis_cohort"]:
            raise ValueError(f"{key}: metadata and split disagree")
        if assignment["primary_analysis_eligible"] != row["primary_analysis_eligible"]:
            raise ValueError(f"{key}: primary eligibility differs")
        side = assignment["split"]
        cluster = assignment["homology_cluster"]
        if side not in {"discovery", "validation"} or not cluster:
            raise ValueError(f"{key}: invalid split or cluster")
        if cluster in owners and owners[cluster] != side:
            raise ValueError(f"{cluster}: homology cluster crosses splits")
        owners[cluster] = side
        primary = row["primary_analysis_eligible"] == "true"
        if primary != (row["analysis_cohort"] == "primary"):
            raise ValueError(f"{key}: inconsistent primary cohort")
        group = f"{side}_primary" if primary else "sensitivity"
        groups[group].append({
            "internal_id": key, "sequence_id": row["sequence_id"],
            "sequence_sha256": row["sequence_sha256"],
            "sequence_length": row["sequence_length"],
            "analysis_cohort": row["analysis_cohort"], "split": side,
            "homology_cluster": cluster,
        })
    for name, rows in groups.items():
        if len(rows) != EXPECTED[name]:
            raise ValueError(f"{name}: expected {EXPECTED[name]}, found {len(rows)}")
        rows.sort(key=lambda row: row["internal_id"])
    return groups


def write_manifests(groups: dict[str, list[dict[str, str]]], output: Path) -> dict[str, str]:
    output.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name, rows in groups.items():
        path = output / f"{name}_manifest.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        hashes[path.name] = sha256(path)
    return hashes


def verify_manifests(groups: dict[str, list[dict[str, str]]], output: Path) -> dict[str, str]:
    hashes = {}
    for name, expected in groups.items():
        path = output / f"{name}_manifest.csv"
        actual = read_csv(path, set(FIELDS))
        if actual != expected:
            raise ValueError(f"{path}: content differs from frozen M2 inputs")
        hashes[path.name] = sha256(path)
    return hashes


def verify_embedding_ids(path: Path, groups: dict[str, list[dict[str, str]]], cohort: str) -> None:
    """Check D's one-row-per-sequence embedding manifest before accepting handoff."""
    rows = read_csv(path, {"internal_id", "sequence_sha256", "analysis_cohort", "split"})
    actual = unique_by_id(rows, "embedding manifest")
    expected = {row["internal_id"]: row for row in groups[cohort]}
    if set(actual) != set(expected):
        raise ValueError(f"embedding IDs: missing {len(set(expected)-set(actual))}, extra {len(set(actual)-set(expected))}")
    for key, row in actual.items():
        for field in ("sequence_sha256", "analysis_cohort", "split"):
            if row[field] != expected[key][field]:
                raise ValueError(f"{key}: embedding {field} differs from input")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["build", "verify"])
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--embeddings-manifest", type=Path)
    parser.add_argument("--cohort", choices=["discovery_primary", "validation_primary", "sensitivity"])
    args = parser.parse_args()
    groups = source_rows(args.metadata, args.split_manifest)
    hashes = write_manifests(groups, args.output_dir) if args.command == "build" else verify_manifests(groups, args.output_dir)
    if args.embeddings_manifest:
        if args.command != "verify" or not args.cohort:
            parser.error("--embeddings-manifest requires verify and --cohort")
        verify_embedding_ids(args.embeddings_manifest, groups, args.cohort)
    print(json.dumps({"dataset_version": DATASET_VERSION, "split_version": SPLIT_VERSION,
                      "counts": {key: len(value) for key, value in groups.items()},
                      "sensitivity_by_cohort": dict(Counter(row["analysis_cohort"] for row in groups["sensitivity"])),
                      "output_sha256": hashes}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
