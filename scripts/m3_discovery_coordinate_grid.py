"""Publish sequence-free coordinates for discovery primary window grids.

Uses the existing M2 scan only. No ESM inference, candidate selection,
conservation calculation or validation-window mapping is performed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import subprocess
import sys

from bioinformatics.m2_pf00741 import digest
from scripts.m3_coordinate_handoff import (
    DEFAULTS, ROOT, WINDOW_COLUMNS, RESIDUE_COLUMNS, load_verified_index,
    map_requests,
)

PUBLIC_RESIDUE_COLUMNS = tuple(
    column for column in RESIDUE_COLUMNS
    if column not in {"residue", "7r1c_reference_residue"}
)


def discovery_windows(records: dict, width: int, step: int):
    """Match D's 1-based inclusive iter_windows; omit incomplete tail windows."""
    if width <= 0 or step <= 0:
        raise ValueError("Window length and step must be positive")
    for identifier in sorted(records):
        record = records[identifier]
        if record["split"] != "discovery" or not record["primary_analysis_eligible"]:
            continue
        for start0 in range(0, len(record["sequence"]) - width + 1, step):
            start, end = start0 + 1, start0 + width
            yield {"window_id": f"{identifier}:{start}-{end}",
                   "internal_id": identifier, "start_1based": start,
                   "end_1based_inclusive": end}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--window-length", type=int, default=30)
    parser.add_argument("--step", type=int, default=5)
    args = parser.parse_args()
    if args.window_length <= 0 or args.step <= 0:
        parser.error("Window length and step must be positive")
    if git("status", "--porcelain"):
        raise ValueError("Coordinate run requires a clean committed working tree")
    if args.out_dir.exists():
        raise FileExistsError("Use a new output directory; existing runs are not overwritten")
    index, summary = load_verified_index({k: ROOT / v for k, v in DEFAULTS.items()})
    args.out_dir.mkdir(parents=True)
    paths = {name: args.out_dir / name for name in (
        "windows.csv", "window_coordinate_map.csv", "residue_coordinate_map.csv")}
    counts, residue_counts, source_ids = Counter(), Counter(), set()
    n_windows = n_residues = 0
    with paths["windows.csv"].open("w", newline="", encoding="utf-8") as requests_file, \
         paths["window_coordinate_map.csv"].open("w", newline="", encoding="utf-8") as window_file, \
         paths["residue_coordinate_map.csv"].open("w", newline="", encoding="utf-8") as residue_file:
        requests_writer = csv.DictWriter(requests_file, fieldnames=(
            "window_id", "internal_id", "start_1based", "end_1based_inclusive"))
        window_writer = csv.DictWriter(window_file, fieldnames=WINDOW_COLUMNS)
        residue_writer = csv.DictWriter(residue_file, fieldnames=PUBLIC_RESIDUE_COLUMNS,
                                       extrasaction="ignore")
        for writer in (requests_writer, window_writer, residue_writer):
            writer.writeheader()
        for request in discovery_windows(index.records, args.window_length, args.step):
            windows, residues = map_requests(index, [request], 0.8)
            assert len(windows) == 1 and len(residues) == args.window_length
            requests_writer.writerow(request)
            window_writer.writerows(windows)
            residue_writer.writerows(residues)
            counts.update([windows[0]["coverage_status"]])
            residue_counts.update(row["mapping_status"] for row in residues)
            source_ids.add(request["internal_id"])
            n_windows += 1
            n_residues += len(residues)
    eligible = {identifier for identifier, record in index.records.items()
                if record["split"] == "discovery" and record["primary_analysis_eligible"]}
    summary.update({
        "schema_version": "1.0", "experiment_id": f"M3-C-GRID-W{args.window_length}-001",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "git": {"commit": git("rev-parse", "HEAD"), "dirty": False},
        "command": sys.argv, "python": platform.python_version(),
        "purpose": "Complete discovery-primary coordinate grid; not selected candidates",
        "window_length": args.window_length, "window_step": args.step,
        "coordinate_convention": "natural/window 1-based inclusive; HMM 1-39; 7R1C author N",
        "join_key": ["internal_id", "start_1based", "end_1based_inclusive"],
        "d_equivalent_join_key": ["internal_id", "raw_start", "raw_end"],
        "d_reference_commit": "79009c83c191a77144181346693b5c2f82e514fc",
        "window_policy": "start0 in range(0, sequence_length-width+1, step); no padded tail",
        "n_primary_discovery_with_windows": len(source_ids),
        "n_primary_discovery_shorter_than_window": len(eligible - source_ids),
        "n_windows": n_windows, "n_residue_rows": n_residues,
        "coverage_status_counts": dict(counts), "residue_status_counts": dict(residue_counts),
        "technical_min_match_fraction": 0.8,
        "candidate_selection_performed": False, "validation_windows_mapped": 0,
        "sequence_residues_in_published_tables": False,
        "script_sha256": digest(Path(__file__)),
        "mapping_module_sha256": digest(ROOT / "bioinformatics/m3_window_coordinates.py"),
        "output_sha256": {name: digest(path) for name, path in paths.items()},
        "limits": ["Coverage threshold only marks technical quality, not candidate acceptance",
                   "HMM/PDB correspondence is homologous position, not a functional label",
                   "Cluster labels/candidate selection require D's actual per-window results",
                   "Validation and post-freeze biological checks await A's real freeze record"],
    })
    (args.out_dir / "mapping_receipt.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"windows": n_windows, "residue_rows": n_residues,
                      "coverage": dict(counts), "out_dir": str(args.out_dir)}))


if __name__ == "__main__":
    main()
