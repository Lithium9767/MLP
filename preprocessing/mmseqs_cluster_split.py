#!/usr/bin/env python3
"""Cluster sequences with MMseqs2 and create a deterministic group-wise split."""

from __future__ import annotations

import argparse
import csv
import random
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


def run_mmseqs(
    executable: str,
    fasta: Path,
    work_dir: Path,
    min_seq_id: float,
    coverage: float,
    threads: int,
) -> Path:
    resolved = shutil.which(executable)
    if not resolved:
        raise SystemExit(
            f"MMseqs2 executable '{executable}' was not found. Install it first or pass --mmseqs."
        )
    work_dir.mkdir(parents=True, exist_ok=True)
    prefix = work_dir / "gvpa"
    tmp_dir = work_dir / "tmp"
    command = [
        resolved,
        "easy-cluster",
        str(fasta),
        str(prefix),
        str(tmp_dir),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        "1",
        "--threads",
        str(threads),
    ]
    subprocess.run(command, check=True)
    cluster_tsv = Path(f"{prefix}_cluster.tsv")
    if not cluster_tsv.is_file():
        raise RuntimeError(f"MMseqs2 did not create {cluster_tsv}")
    return cluster_tsv


def read_clusters(path: Path) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise ValueError(f"Malformed cluster row at {path}:{line_number}")
            representative, member = parts[0], parts[1]
            clusters[representative].append(member)
    if not clusters:
        raise ValueError(f"No clusters found in {path}")
    return dict(clusters)


def assign_clusters(
    clusters: dict[str, list[str]], ratios: tuple[float, float, float], seed: int
) -> dict[str, str]:
    names = ("train", "val", "test")
    total = sum(len(members) for members in clusters.values())
    targets = {name: ratio * total for name, ratio in zip(names, ratios)}
    counts = {name: 0 for name in names}
    rng = random.Random(seed)
    items = list(clusters.items())
    rng.shuffle(items)
    items.sort(key=lambda item: len(item[1]), reverse=True)

    cluster_split: dict[str, str] = {}
    for representative, members in items:
        split = max(names, key=lambda name: (targets[name] - counts[name], -counts[name]))
        cluster_split[representative] = split
        counts[split] += len(members)
    return cluster_split


def write_split(
    clusters: dict[str, list[str]], cluster_split: dict[str, str], output: Path
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sequence_id", "similarity_cluster", "split"]
        )
        writer.writeheader()
        for representative in sorted(clusters):
            for member in sorted(clusters[representative]):
                writer.writerow(
                    {
                        "sequence_id": member,
                        "similarity_cluster": representative,
                        "split": cluster_split[representative],
                    }
                )


def parse_ratios(value: str) -> tuple[float, float, float]:
    values = tuple(float(part) for part in value.split(","))
    if len(values) != 3 or any(part <= 0 for part in values):
        raise argparse.ArgumentTypeError("Ratios must be three positive comma-separated numbers")
    total = sum(values)
    return tuple(part / total for part in values)  # type: ignore[return-value]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--mmseqs", default="mmseqs")
    parser.add_argument("--min-seq-id", type=float, default=0.8)
    parser.add_argument("--coverage", type=float, default=0.8)
    parser.add_argument("--ratios", type=parse_ratios, default=parse_ratios("0.8,0.1,0.1"))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    cluster_tsv = run_mmseqs(
        args.mmseqs,
        args.fasta,
        args.work_dir,
        args.min_seq_id,
        args.coverage,
        args.threads,
    )
    clusters = read_clusters(cluster_tsv)
    assignments = assign_clusters(clusters, args.ratios, args.seed)
    write_split(clusters, assignments, args.output)
    print(f"Wrote {sum(map(len, clusters.values()))} sequences in {len(clusters)} clusters")


if __name__ == "__main__":
    main()
