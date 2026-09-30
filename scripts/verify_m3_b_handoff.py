#!/usr/bin/env python3
"""Check the frozen M3 B handoff without extraction or model validation.

This is an engineering check. It cannot establish who prepared the package or
substitute for B's own handoff and a non-author review.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import subprocess
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INNER = "m3_b_handoff/m3_b_handoff_bundle.zip"
TABLES = (
    "discovery_primary_manifest.csv",
    "validation_primary_manifest.csv",
    "sensitivity_manifest.csv",
)
INNER_FILES = set(TABLES) | {"metadata.csv", "split_manifest.csv", "sequences_for_clustering.fasta"}
EXPECTED_SHA256 = {
    # Previously frozen M2 inputs, also referenced in the M3 discovery config.
    "metadata.csv": "62893e31a2dbc245f29c0e9af0a282dcb5099d0e95b96db1a9cead2ac933004e",
    "split_manifest.csv": "0a6bd3b44d284a407e1cc9d5be758f0f9439c28395deb6d62d2c787f95341e9c",
    "sequences_for_clustering.fasta": "4ff33c04be1d28f42ef5afa7d67261229d19c8e301784ea59750030196e8b122",
}
EXPECTED_COHORTS = {
    "primary": 1721,
    "sensitivity_type_conflict": 305,
    "sensitivity_partial_and_type_conflict": 30,
    "sensitivity_partial": 20,
    "excluded_sequence_qc": 2,
}
EXPECTED_SPLITS = {"discovery": 1453, "validation": 623}
EXPECTED_CLUSTERS = {"discovery": 335, "validation": 143}
EXPECTED_RECORDS = {"metadata": 2078, "split": 2076}
EXPECTED_MANIFESTS = {
    "discovery_primary_manifest.csv": 1202,
    "validation_primary_manifest.csv": 519,
    "sensitivity_manifest.csv": 355,
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def table(data: bytes, name: str) -> list[dict[str, str]]:
    try:
        return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
    except UnicodeError as exc:
        raise ValueError(f"Invalid UTF-8 in {name}") from exc


def unique(rows: list[dict[str, str]], name: str) -> dict[str, dict[str, str]]:
    require(all(row.get("internal_id") for row in rows), f"Missing internal_id: {name}")
    indexed = {row["internal_id"]: row for row in rows}
    require(len(indexed) == len(rows), f"Duplicate internal_id: {name}")
    return indexed


def checked_zip_members(archive: zipfile.ZipFile, expected: set[str], label: str) -> None:
    names = [item.filename for item in archive.infolist() if not item.is_dir()]
    require(len(names) == len(set(names)), f"Duplicate ZIP member: {label}")
    require(set(names) == expected, f"Unexpected {label} ZIP members: {set(names) ^ expected}")


def parse_fasta(data: bytes) -> dict[str, str]:
    records: dict[str, str] = {}
    key = None
    residues: list[str] = []
    for line in data.decode("ascii").splitlines():
        if line.startswith(">"):
            if key is not None:
                require(key not in records, f"Duplicate FASTA ID: {key}")
                records[key] = "".join(residues)
            key, residues = line[1:].split()[0], []
            require(bool(key), "Blank FASTA ID")
        else:
            require(key is not None, "FASTA sequence before header")
            residues.append(line.strip())
    if key is not None:
        require(key not in records, f"Duplicate FASTA ID: {key}")
        records[key] = "".join(residues)
    return records


def check_rows(metadata: list[dict[str, str]], split: list[dict[str, str]],
               manifests: dict[str, list[dict[str, str]]], fasta: dict[str, str]) -> dict:
    meta_by = unique(metadata, "metadata.csv")
    split_by = unique(split, "split_manifest.csv")
    require(len(metadata) == EXPECTED_RECORDS["metadata"] and
            len(split) == EXPECTED_RECORDS["split"], "Frozen record counts changed")
    cohorts = Counter(row["analysis_cohort"] for row in metadata)
    require(cohorts == EXPECTED_COHORTS, f"Unexpected metadata cohorts: {cohorts}")
    require(set(split_by) == set(fasta), "FASTA/split ID mismatch")
    require(set(split_by) == {key for key, row in meta_by.items()
                              if row["sequence_qc_eligible"] == "true"}, "Eligible/split ID mismatch")
    for key, row in meta_by.items():
        sequence = row["sequence"]
        require(digest(sequence.encode("ascii")) == row["sequence_sha256"], f"Sequence hash mismatch: {key}")
        require(len(sequence) == int(row["sequence_length"]), f"Sequence length mismatch: {key}")
        if key in split_by:
            require(fasta[key] == sequence, f"FASTA/metadata sequence mismatch: {key}")
            s = split_by[key]
            require(s["sequence_id"] == row["sequence_id"] and
                    s["analysis_cohort"] == row["analysis_cohort"] and
                    s["primary_analysis_eligible"] == row["primary_analysis_eligible"],
                    f"Metadata/split mismatch: {key}")
    split_counts = Counter(row["split"] for row in split)
    require(split_counts == EXPECTED_SPLITS, f"Unexpected split counts: {split_counts}")
    cluster_sets = {name: {row["homology_cluster"] for row in split if row["split"] == name}
                    for name in EXPECTED_SPLITS}
    require({k: len(v) for k, v in cluster_sets.items()} == EXPECTED_CLUSTERS,
            "Unexpected homology cluster counts")
    require(not cluster_sets["discovery"] & cluster_sets["validation"], "Homology cluster leakage")
    seen: set[str] = set()
    for name, rows in manifests.items():
        indexed = unique(rows, name)
        require(len(rows) == EXPECTED_MANIFESTS[name], f"Unexpected manifest count: {name}")
        require(not seen & set(indexed), f"Manifest overlap: {name}")
        seen.update(indexed)
        for key, row in indexed.items():
            require(key in split_by, f"Manifest ID absent from split: {key}")
            source = split_by[key]
            meta = meta_by[key]
            require(all(row[field] == source[field] for field in
                        ("sequence_id", "analysis_cohort", "split", "homology_cluster")),
                    f"Manifest/split mismatch: {key}")
            require(row["sequence_sha256"] == meta["sequence_sha256"] and
                    int(row["sequence_length"]) == len(meta["sequence"]),
                    f"Manifest/metadata mismatch: {key}")
            if name == "discovery_primary_manifest.csv":
                require(row["analysis_cohort"] == "primary" and row["split"] == "discovery",
                        f"Wrong discovery cohort: {key}")
            elif name == "validation_primary_manifest.csv":
                require(row["analysis_cohort"] == "primary" and row["split"] == "validation",
                        f"Wrong validation cohort: {key}")
            else:
                require(row["analysis_cohort"].startswith("sensitivity_"), f"Wrong sensitivity cohort: {key}")
    require(seen == set(split_by), "Manifests do not partition split IDs")
    return {
        "metadata_rows": len(metadata), "split_rows": len(split),
        "metadata_cohorts": dict(sorted(cohorts.items())),
        "split_counts": dict(sorted(split_counts.items())),
        "cluster_counts": {k: len(v) for k, v in cluster_sets.items()},
        "cross_split_clusters": 0,
        "manifest_rows": {k: len(v) for k, v in manifests.items()},
        "fasta_records": len(fasta),
    }


def verify_bundle(path: Path) -> dict:
    with zipfile.ZipFile(path) as outer:
        expected_outer = {INNER} | {f"m3_b_handoff/{name}" for name in TABLES}
        checked_zip_members(outer, expected_outer, "outer")
        require(outer.testzip() is None, "Outer ZIP CRC failure")
        inner_bytes = outer.read(INNER)
        with zipfile.ZipFile(io.BytesIO(inner_bytes)) as inner:
            checked_zip_members(inner, INNER_FILES, "inner")
            require(sum(item.file_size for item in inner.infolist()) < 10_000_000, "Excessive uncompressed size")
            require(inner.testzip() is None, "Inner ZIP CRC failure")
            data = {name: inner.read(name) for name in INNER_FILES}
            for name in TABLES:
                require(outer.read(f"m3_b_handoff/{name}") == data[name], f"Outer/inner mismatch: {name}")
    hashes = {name: digest(data[name]) for name in sorted(data)}
    for name, expected in EXPECTED_SHA256.items():
        require(hashes[name] == expected, f"Frozen M2 input hash mismatch: {name}")
    counts = check_rows(table(data["metadata.csv"], "metadata.csv"),
                        table(data["split_manifest.csv"], "split_manifest.csv"),
                        {name: table(data[name], name) for name in TABLES},
                        parse_fasta(data["sequences_for_clustering.fasta"]))
    return {"outer_zip_sha256": digest(path.read_bytes()), "inner_zip_sha256": digest(inner_bytes),
            "inner_file_sha256": hashes, "checks": counts,
            "dataset_version": "gvpa-recognition-c7f6f005d717",
            "split_version": "homology-8b9005e2d9-s42"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()
    require(not status, "Formal receipt requires a clean committed worktree")
    run_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    result = verify_bundle(args.bundle)
    result.update({"status": "assistant_technical_precheck_passed_pending_B_attestation",
                   "run_commit": run_sha, "dirty": False,
                   "run_utc": datetime.now(timezone.utc).isoformat(),
                   "python": sys.version.split()[0], "platform": platform.platform(),
                   "command": "python scripts/verify_m3_b_handoff.py --bundle <local-bundle-path> --receipt <receipt-path>",
                   "bundle_filename": args.bundle.name,
                   "limitations": ["Package creator and B's independent check are unverified",
                                   "Original large-file shared-storage URI and access rights are unverified",
                                   "Validation sequences were checked for integrity only; no M3 model validation ran"]})
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "run_commit": run_sha,
                      "counts": result["checks"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
