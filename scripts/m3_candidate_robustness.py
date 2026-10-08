#!/usr/bin/env python3
"""Discovery-only sequence-level diagnostics for an exploratory M3 cluster.

These post-selection diagnostics do not establish function or A approval.
No validation record is read or used for parameter selection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import random
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def validate_discovery_rows(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("No discovery windows")
    keys = set()
    for row in rows:
        if row["split"] != "discovery":
            raise ValueError("Only discovery windows may be analyzed")
        start = int(row["raw_start"])
        end = int(row["raw_end"])
        length = int(row["sequence_length"])
        if not 1 <= start <= end <= length or end - start + 1 != int(row["window_length"]):
            raise ValueError("Invalid discovery window bounds")
        key = (row["internal_id"], start, end)
        if key in keys:
            raise ValueError("Duplicate discovery window")
        keys.add(key)


def match_same_sequence_position(
    rows: list[dict[str, str]], cluster: int, tolerance: float
) -> tuple[list[tuple[dict[str, str], dict[str, str]]], int]:
    """Choose the nearest noncandidate window in the same sequence and position band."""
    controls = defaultdict(list)
    targets = []
    for row in rows:
        if int(row["cluster"]) == cluster:
            targets.append(row)
        else:
            controls[row["internal_id"]].append(row)
    pairs = []
    unmatched = 0
    for target in targets:
        position = float(target["normalized_start"])
        nearby = [row for row in controls[target["internal_id"]]
                  if abs(float(row["normalized_start"]) - position) <= tolerance]
        if not nearby:
            unmatched += 1
            continue
        control = min(nearby, key=lambda row: (
            abs(float(row["normalized_start"]) - position), int(row["raw_start"])))
        pairs.append((target, control))
    return pairs, unmatched


def state_present(row: dict[str, str], state: int) -> bool:
    return str(state) in row["hmm_states"].split(";")


def sequence_level_delta(
    pairs: list[tuple[dict[str, str], dict[str, str]]], state: int
) -> dict[str, float]:
    by_id = defaultdict(list)
    for target, control in pairs:
        if target["internal_id"] != control["internal_id"]:
            raise ValueError("Cross-sequence control")
        by_id[target["internal_id"]].append(int(state_present(target, state)) -
                                                int(state_present(control, state)))
    return {key: statistics.mean(values) for key, values in by_id.items()}


def bootstrap_mean_ci(values: list[float], seed: int, repeats: int = 2000) -> tuple[float, float]:
    if not values:
        raise ValueError("No matched sequences for bootstrap")
    rng = random.Random(seed)
    estimates = sorted(statistics.mean(rng.choices(values, k=len(values))) for _ in range(repeats))
    return estimates[int(.025 * repeats)], estimates[int(.975 * repeats)]


def length_bin(length: int) -> str:
    if length <= 100:
        return "<=100"
    if length <= 130:
        return "101-130"
    if length <= 160:
        return "131-160"
    return ">160"


def shorter_length_auc(candidate_lengths: list[int], other_lengths: list[int]) -> float:
    """Fraction of candidate/other sequence pairs ordered by shorter length."""
    if not candidate_lengths or not other_lengths:
        raise ValueError("Both sequence groups are required for length AUC")
    wins = sum(a < b for a in candidate_lengths for b in other_lengths)
    ties = sum(a == b for a in candidate_lengths for b in other_lengths)
    return (wins + .5 * ties) / (len(candidate_lengths) * len(other_lengths))


def analyze(rows: list[dict[str, str]], cluster: int, state: int, seed: int) -> tuple[dict, list[dict], list[dict]]:
    validate_discovery_rows(rows)
    targets = [row for row in rows if int(row["cluster"]) == cluster]
    if not targets:
        raise ValueError("Cluster has no discovery windows")
    lengths = {}
    for row in rows:
        key = row["internal_id"]
        length = int(row["sequence_length"])
        if key in lengths and lengths[key] != length:
            raise ValueError("Inconsistent sequence length")
        lengths[key] = length
    participant_ids = {row["internal_id"] for row in targets}
    other_ids = set(lengths) - participant_ids
    if not other_ids:
        raise ValueError("No noncandidate sequences for comparison")
    bins = []
    for label in ("<=100", "101-130", "131-160", ">160"):
        members = {key for key, length in lengths.items() if length_bin(length) == label}
        participants = len(members & participant_ids)
        bins.append({"length_bin": label, "n_sequences": len(members),
                     "n_with_candidate": participants,
                     "candidate_sequence_fraction": participants / len(members) if members else ""})
    diagnostics = {}
    sequence_rows = []
    for tolerance in (.05, .10):
        pairs, unmatched = match_same_sequence_position(rows, cluster, tolerance)
        deltas = sequence_level_delta(pairs, state)
        key = f"same_sequence_position_{tolerance:.2f}"
        diagnostics[key] = {"matched_windows": len(pairs), "unmatched_windows": unmatched,
                            "matched_sequences": len(deltas),
                            "target_state_fraction": statistics.mean(state_present(a, state) for a, _ in pairs) if pairs else None,
                            "control_state_fraction": statistics.mean(state_present(b, state) for _, b in pairs) if pairs else None,
                            "mean_sequence_delta": statistics.mean(deltas.values()) if deltas else None,
                            "bootstrap_sequence_ci95": bootstrap_mean_ci(list(deltas.values()), seed) if deltas else None}
        if tolerance == .10:
            sequence_rows = [{"internal_id": seq_id, "matched_window_count": sum(a["internal_id"] == seq_id for a, _ in pairs),
                              "state_fraction_delta": value} for seq_id, value in sorted(deltas.items())]
    summary = {"cluster": cluster, "post_selection_state": state,
               "n_discovery_sequences": len(lengths), "n_discovery_windows": len(rows),
               "n_candidate_sequences": len(participant_ids), "n_candidate_windows": len(targets),
               "candidate_state_fraction": statistics.mean(state_present(row, state) for row in targets),
               "candidate_median_length": statistics.median(lengths[key] for key in participant_ids),
               "noncandidate_median_length": statistics.median(lengths[key] for key in other_ids),
               "shorter_length_only_auc_for_cluster_participation": shorter_length_auc(
                   [lengths[key] for key in participant_ids], [lengths[key] for key in other_ids]),
               "length_auc_interpretation": "Descriptive ability of sequence length to distinguish cluster participation; not biological prediction performance.",
               "candidate_median_normalized_start": statistics.median(float(row["normalized_start"]) for row in targets),
               "control_diagnostics": diagnostics,
               "interpretation": "Post-selected discovery diagnostic; windows overlap and are not independent experimental units. No functional p-value or validation claim."}
    return summary, bins, sequence_rows


def clean_git() -> str:
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip():
        raise ValueError("Formal run requires clean committed worktree")
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--source-receipt", type=Path, required=True)
    parser.add_argument("--cluster", type=int, required=True)
    parser.add_argument("--state", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run_commit = clean_git()
    source = json.loads(args.source_receipt.read_text(encoding="utf-8"))
    if source["result"]["validation_used"] or source["config"]["fit_split"] != "discovery":
        raise ValueError("Source run used validation")
    if sha256(args.windows) != source["outputs_sha256"]["candidate_windows.csv"]:
        raise ValueError("Candidate window input hash mismatch")
    if args.output_dir.exists():
        raise ValueError("Choose a fresh output directory")
    summary, bins, sequence_rows = analyze(read_csv(args.windows), args.cluster, args.state, args.seed)
    if not sequence_rows:
        raise ValueError("No matched sequence controls at tolerance 0.10")
    args.output_dir.mkdir(parents=True)
    write_json(args.output_dir / "summary.json", summary)
    write_csv(args.output_dir / "length_bins.csv", bins)
    write_csv(args.output_dir / "sequence_deltas.csv", sequence_rows)
    receipt = {"status": "completed_discovery_diagnostic", "run_commit": run_commit,
               "dirty": False, "utc": datetime.now(timezone.utc).isoformat(),
               "python": sys.version, "platform": platform.platform(),
               "command": sys.argv, "seed": args.seed,
               "source_run_commit": source["provenance"]["commit"],
               "source_receipt_sha256": sha256(args.source_receipt),
               "candidate_windows_sha256": sha256(args.windows),
               "validation_used": False,
               "outputs_sha256": {p.name: sha256(p) for p in args.output_dir.iterdir() if p.is_file()}}
    write_json(args.output_dir / "run_receipt.json", receipt)
    print(json.dumps({"summary": summary, "run_commit": run_commit}, ensure_ascii=False))


if __name__ == "__main__":
    main()
