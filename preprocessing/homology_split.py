#!/usr/bin/env python3
"""Create a deterministic discovery/validation split from MMseqs2 clusters."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_ratios(value: str) -> dict[str, float]:
    result: dict[str, float] = {}
    try:
        for item in value.split(","):
            name, raw = item.split("=", 1)
            name = name.strip()
            ratio = float(raw)
            if not name or name in result or ratio <= 0:
                raise ValueError
            result[name] = ratio
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Ratios must look like discovery=0.7,validation=0.3") from exc
    if len(result) < 2:
        raise argparse.ArgumentTypeError("At least two positive split ratios are required")
    total = sum(result.values())
    return {name: ratio / total for name, ratio in result.items()}


def read_eligible_ids(metadata_path: Path) -> set[str]:
    with metadata_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sequence_id", "sequence_qc_eligible"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Metadata must contain {sorted(required)}")
        ids = {
            row["sequence_id"].strip()
            for row in reader
            if row["sequence_qc_eligible"].strip().casefold() == "true"
        }
    if not ids or "" in ids:
        raise ValueError("No valid eligible sequence IDs found")
    return ids


def read_clusters(path: Path, expected_ids: set[str]) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = defaultdict(list)
    owner: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
                raise ValueError(f"Malformed cluster row at {path}:{line_number}")
            representative, member = (part.strip() for part in parts)
            if member not in expected_ids:
                raise ValueError(f"Cluster member {member!r} is not an eligible sequence")
            previous = owner.get(member)
            if previous is not None:
                raise ValueError(
                    f"Cluster member {member!r} appears more than once ({previous!r}, {representative!r})"
                )
            owner[member] = representative
            clusters[representative].append(member)
    missing = sorted(expected_ids - owner.keys())
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"Cluster file is missing {len(missing)} eligible IDs; first: {preview}")
    if not clusters:
        raise ValueError("No clusters found")
    return {representative: sorted(members) for representative, members in sorted(clusters.items())}


def _tie_key(seed: int, representative: str) -> str:
    return hashlib.sha256(f"{seed}:{representative}".encode("utf-8")).hexdigest()


def assign_clusters(clusters: dict[str, list[str]], ratios: dict[str, float], seed: int) -> dict[str, str]:
    total = sum(map(len, clusters.values()))
    targets = {name: ratio * total for name, ratio in ratios.items()}
    counts = {name: 0 for name in ratios}
    ordered = sorted(clusters, key=lambda representative: (-len(clusters[representative]), _tie_key(seed, representative)))
    assignments: dict[str, str] = {}
    for representative in ordered:
        split = min(
            ratios,
            key=lambda name: (counts[name] / targets[name], -ratios[name], name),
        )
        assignments[representative] = split
        counts[split] += len(clusters[representative])
    return assignments


def write_outputs(
    clusters: dict[str, list[str]],
    assignments: dict[str, str],
    metadata_path: Path,
    cluster_path: Path,
    ratios: dict[str, float],
    seed: int,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "split_manifest.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence_id", "homology_cluster", "split"])
        writer.writeheader()
        for representative in sorted(clusters):
            for member in clusters[representative]:
                writer.writerow({
                    "sequence_id": member,
                    "homology_cluster": representative,
                    "split": assignments[representative],
                })

    sequence_counts = Counter()
    cluster_counts = Counter()
    for representative, members in clusters.items():
        split = assignments[representative]
        sequence_counts[split] += len(members)
        cluster_counts[split] += 1
    cluster_sha = sha256_file(cluster_path)
    summary = {
        "schema_version": "1.0",
        "split_version": f"homology-{cluster_sha[:10]}-s{seed}",
        "method": "whole-homology-cluster greedy allocation",
        "seed": seed,
        "target_ratios": ratios,
        "sequence_counts": dict(sequence_counts),
        "cluster_counts": dict(cluster_counts),
        "total_sequences": sum(sequence_counts.values()),
        "total_clusters": len(clusters),
        "metadata_sha256": sha256_file(metadata_path),
        "mmseqs_cluster_tsv_sha256": cluster_sha,
        "split_manifest_sha256": sha256_file(manifest_path),
        "cluster_leakage": False,
    }
    (output_dir / "split_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--clusters-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ratios", type=parse_ratios, default=parse_ratios("discovery=0.7,validation=0.3"))
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    expected_ids = read_eligible_ids(args.metadata)
    clusters = read_clusters(args.clusters_tsv, expected_ids)
    assignments = assign_clusters(clusters, args.ratios, args.seed)
    summary = write_outputs(
        clusters, assignments, args.metadata, args.clusters_tsv, args.ratios, args.seed, args.output_dir
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
