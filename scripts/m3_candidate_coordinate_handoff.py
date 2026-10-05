"""Verify D's original discovery windows and attach C's frozen coordinates.

No candidate selection, model inference or validation analysis is performed.
Cluster -1 (noise) and all unmapped windows remain in the output.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

from bioinformatics.m2_pf00741 import digest
from scripts.m3_coordinate_handoff import DEFAULTS, ROOT, WINDOW_COLUMNS, load_verified_index, map_requests
from scripts.m3_discovery_coordinate_grid import PUBLIC_RESIDUE_COLUMNS

EXPECTED_WINDOWS_SHA = "d8b82d04cdee709507567a7e626e69f985592537a2580bc60485224d0fbec328"
EXPECTED_D_COMMIT = "dc066b76c637cdb1b90b7c8490e28ad81c4cff52"
D_FIELDS = ("d_cluster", "d_probability", "d_review_status", "d_interpretation")


def validate_window(row, record, mapped):
    """Reject wrong split/hash/bounds and discrepancies with C's raw mapping."""
    if row["split"] != "discovery" or record["split"] != "discovery" or not record["primary_analysis_eligible"]:
        raise ValueError("Only discovery primary windows are permitted")
    if row["sequence_sha256"] != record["sequence_sha256"]:
        raise ValueError("Sequence SHA mismatch")
    start, end = int(row["raw_start"]), int(row["raw_end"])
    if not (1 <= start <= end <= len(record["sequence"])) or end-start+1 != 30 or (start-1) % 5:
        raise ValueError("Invalid window bounds/grid")
    if int(row["window_length"]) != 30 or int(row["sequence_length"]) != len(record["sequence"]):
        raise ValueError("Window/sequence length mismatch")
    states = [int(x) for x in row["hmm_states"].split(";") if x]
    if states != json.loads(mapped["hmm_match_states_json"]):
        raise ValueError("D/C HMM state mismatch")
    if abs(float(row["mapped_fraction"]) - float(mapped["match_fraction"])) > 1e-12:
        raise ValueError("D/C mapped fraction mismatch")
    probability = float(row["probability"])
    if not 0 <= probability <= 1 or int(row["cluster"]) not in {-1, 0, 1, 2, 3}:
        raise ValueError("Invalid D cluster/probability")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zip", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise ValueError("Run from a clean committed working tree")
    if args.out_dir.exists():
        raise FileExistsError("Use a new output directory")
    with zipfile.ZipFile(args.zip) as bundle:
        if len(bundle.namelist()) != len(set(bundle.namelist())):
            raise ValueError("Duplicate ZIP members")
        payload = {name: bundle.read(name) for name in (
            "candidate_windows.csv", "run_receipt.json", "candidate_review.csv", "HANDOFF.txt")}
    sha = lambda b: hashlib.sha256(b).hexdigest()
    receipt = json.loads(payload["run_receipt.json"])
    if sha(payload["candidate_windows.csv"]) != EXPECTED_WINDOWS_SHA or receipt["outputs_sha256"]["candidate_windows.csv"] != EXPECTED_WINDOWS_SHA:
        raise ValueError("D window payload differs from published run hash")
    if receipt["provenance"]["commit"] != EXPECTED_D_COMMIT or receipt["provenance"]["dirty"] or receipt["result"]["validation_used"]:
        raise ValueError("Unexpected D run provenance or validation use")
    index, summary = load_verified_index({k: ROOT/v for k, v in DEFAULTS.items()})
    if receipt["config"]["coordinates_sha256"] != digest(ROOT/DEFAULTS["coordinates"]):
        raise ValueError("D used different source coordinates")
    reviews = {r["cluster"]: r for r in csv.DictReader(io.StringIO(payload["candidate_review.csv"].decode("utf-8-sig")))}
    if set(reviews) != {"0", "1", "2", "3"}:
        raise ValueError("Missing cluster review")
    args.out_dir.mkdir(parents=True)
    for name, data in payload.items():
        (args.out_dir/name).write_bytes(data)
    window_path = args.out_dir/"candidate_window_coordinate_map.csv"
    residue_path = args.out_dir/"candidate_residue_coordinate_map.csv"
    seen, identifiers, counts, cluster_counts = set(), set(), Counter(), Counter()
    cluster_summary = {}
    n_residues = 0
    with window_path.open("w", encoding="utf-8", newline="") as wf, residue_path.open("w", encoding="utf-8", newline="") as rf:
        ww = csv.DictWriter(wf, fieldnames=(*WINDOW_COLUMNS, *D_FIELDS))
        rw = csv.DictWriter(rf, fieldnames=(*PUBLIC_RESIDUE_COLUMNS, "d_cluster"), extrasaction="ignore")
        ww.writeheader(); rw.writeheader()
        for row in csv.DictReader(io.StringIO(payload["candidate_windows.csv"].decode("utf-8-sig"))):
            identifier = row["internal_id"]
            if identifier not in index.records:
                raise ValueError("Unknown internal_id")
            start, end = int(row["raw_start"]), int(row["raw_end"])
            key = (identifier, start, end)
            if key in seen:
                raise ValueError("Duplicate D window")
            seen.add(key); identifiers.add(identifier)
            request = {"window_id": f"{identifier}:{start}-{end}", "internal_id": identifier,
                       "start_1based": start, "end_1based_inclusive": end}
            windows, residues = map_requests(index, [request], 0.8)
            mapped = windows[0]
            validate_window(row, index.records[identifier], mapped)
            cluster = row["cluster"]
            review = reviews.get(cluster, {"review_status": "noise_not_candidate", "interpretation": "HDBSCAN noise; retained without assigning a candidate cluster"})
            mapped.update(d_cluster=cluster, d_probability=row["probability"],
                          d_review_status=review["review_status"], d_interpretation=review["interpretation"])
            ww.writerow(mapped)
            for residue in residues:
                rw.writerow({**residue, "d_cluster": cluster})
            n_residues += len(residues)
            counts[mapped["coverage_status"]] += 1
            cluster_counts[cluster] += 1
            stats = cluster_summary.setdefault(cluster, {"sequences": set(), "coverage": Counter(), "states": set(), "author_positions": set()})
            stats["sequences"].add(identifier); stats["coverage"][mapped["coverage_status"]] += 1
            stats["states"].update(json.loads(mapped["hmm_match_states_json"]))
            stats["author_positions"].update(int(r["7r1c_author_residue_number"]) for r in residues if r["7r1c_author_residue_number"] not in (None, ""))
    expected = {(identifier, start+1, start+30) for identifier, record in index.records.items()
                if record["split"] == "discovery" and record["primary_analysis_eligible"]
                for start in range(0, len(record["sequence"])-30+1, 5)}
    if seen != expected or len(seen) != 24581 or len(identifiers) != 1202:
        raise ValueError("D payload is not the complete frozen discovery grid")
    summary.update(experiment_id="M3-C-CANDIDATE-001", status="completed_discovery_coordinate_integration",
        run_utc=datetime.now(timezone.utc).isoformat(), command=sys.argv,
        git={"commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), "dirty": False},
        d_run_commit=EXPECTED_D_COMMIT, d_zip_sha256=digest(args.zip),
        d_files_sha256={k: sha(v) for k,v in payload.items()},
        n_windows=len(seen), n_sequences=len(identifiers), n_residue_rows=n_residues,
        n_clustered_windows=sum(v for k,v in cluster_counts.items() if k != "-1"),
        cluster_window_counts=dict(cluster_counts), coverage_status_counts=dict(counts),
        all_window_hashes_bounds_hmm_states_and_match_fractions_checked=True,
        validation_used=False, A_candidate_freeze_confirmed=False,
        interpretation="Exploratory D cluster membership only; PF00741 hit is not a functional label",
        cluster_summary={k: {"n_sequences": len(s["sequences"]), "coverage": dict(s["coverage"]),
            "hmm_states_union": sorted(s["states"]), "7r1c_author_positions_union": sorted(s["author_positions"]),
            "interpretation": reviews.get(k, {}).get("interpretation", "noise_not_candidate")}
            for k,s in cluster_summary.items()},
        outputs_sha256={p.name: digest(p) for p in (window_path, residue_path)})
    (args.out_dir/"mapping_receipt.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"windows": len(seen), "residue_rows": n_residues, "clusters": dict(cluster_counts)}))


if __name__ == "__main__":
    main()
