"""Verify frozen M2 coordinates and hand M3 windows to PF00741/7R1C.

Without ``--windows`` this is a prelaunch integrity check only. With a CSV of
1-based inclusive windows it emits per-window and per-residue maps in an ignored
shared directory. It never ranks candidates or estimates conservation. The
validation mode requires an A-owned freeze record before opening any windows.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from bioinformatics.m2_pf00741 import HMM_SHA256, PDB_SHA256, digest, validate_b_handoff
from bioinformatics.m3_window_coordinates import load_coordinate_index


ROOT = Path(__file__).resolve().parents[1]
DEFAULTS = {
    "metadata": "data/processed/gvpa_v1_reproduction/metadata.csv",
    "fasta": "data/processed/gvpa_v1_reproduction/sequences_for_clustering.fasta",
    "split_manifest": "data/processed/gvpa_v1_reproduction/split/split_manifest.csv",
    "audit_summary": "results/data_audit/gvpa_v1_audit_summary.json",
    "split_summary": "results/data_audit/gvpa_v1_split_summary.json",
    "status": "data/processed/m2_c/pf00741_scan_002/per_sequence_status.csv",
    "coordinates": "data/processed/m2_c/pf00741_scan_002/hmm_coordinate_map.csv",
    "structure": "data/processed/m2_c/7r1c_002/7r1c_residue_hmm_map.csv",
    "scan_receipt": "results/bioinformatics/pf00741_scan_002_summary.json",
    "coordinate_receipt": "results/bioinformatics/hmm_coordinate_002_summary.json",
    "structure_receipt": "results/bioinformatics/structure_7r1c_summary.json",
}
PINNED_M2_RECEIPT_SHA256 = {
    "scan_receipt": "befbebcc8358a35292839b2d32375684effee14faba868f9bfc6ce94ffacb781",
    "coordinate_receipt": "71e04b8e41d7d1a1e4cae56dfb9213849c51e465ae335bb22e6e6353987d2df6",
    "structure_receipt": "81d3e3cc793d30d7fbad4ba71dd91a58eff58ded086f315f4a22b8565a478ffc",
}
WINDOW_FIELDS = ("internal_id", "start_1based", "end_1based_inclusive")
WINDOW_COLUMNS = ("window_id", "internal_id", "sequence_sha256", "split",
                  "analysis_cohort", "primary_analysis_eligible", "start_1based",
                  "end_1based_inclusive", "window_length", "m2_hmm_status",
                  "hmm_accession", "7r1c_author_chain",
                  "coverage_status", "match_fraction", "aligned_fraction",
                  "n_match_residues", "n_insertion_residues",
                  "n_unaligned_residues", "n_ambiguous_residues", "n_domains",
                  "n_reference_modeled", "n_reference_unmodeled",
                  "domain_numbers_json", "hmm_match_states_json",
                  "hmm_deleted_states_internal_json")
RESIDUE_COLUMNS = ("window_id", "internal_id", "raw_position", "residue",
                   "mapping_status", "domain_number", "hmm_match_state",
                   "coordinate_options_json", "reference_status",
                   "hmm_accession", "7r1c_author_chain",
                   "7r1c_deposited_position", "7r1c_author_residue_number",
                   "7r1c_author_insertion_code", "7r1c_reference_residue",
                   "7r1c_modeled", "7r1c_secondary_structure")


def _json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _verify_sha(path: Path, expected: str) -> str:
    actual = digest(path)
    if actual != expected:
        raise ValueError(f"Frozen SHA-256 mismatch: {path}")
    return actual


def load_verified_index(paths: dict[str, Path]):
    """Reuse M2's B-input validator, then verify M2 raw/structure receipts."""
    records, handoff = validate_b_handoff(
        paths["metadata"], paths["fasta"], paths["split_manifest"],
        paths["audit_summary"], paths["split_summary"])
    for name, expected in PINNED_M2_RECEIPT_SHA256.items():
        _verify_sha(paths[name], expected)
    scan = _json(paths["scan_receipt"])
    coordinate = _json(paths["coordinate_receipt"])
    structure = _json(paths["structure_receipt"])
    if any(item.get("git", {}).get("dirty") is not False
           for item in (scan, coordinate, structure)):
        raise ValueError("M2 scan/structure must have clean-run receipts")
    if scan.get("input") != coordinate.get("input") or scan.get("input") != handoff:
        raise ValueError("M2 scan and coordinate receipts do not match frozen B handoff")
    if (scan.get("hmm", {}).get("accession") != "PF00741.24"
            or scan.get("hmm", {}).get("length") != 39
            or scan.get("hmm", {}).get("sha256") != HMM_SHA256
            or coordinate.get("hmm") != scan.get("hmm")
            or structure.get("hmm") != scan.get("hmm")
            or structure.get("pdb_id") != "7R1C" or structure.get("chain") != "N"
            or structure.get("pdb_sha256") != PDB_SHA256):
        raise ValueError("M2 HMM/7R1C reference identity differs from frozen receipt")
    raw_hashes = scan.get("raw_output_sha256", {})
    if coordinate.get("raw_output_sha256") != raw_hashes:
        raise ValueError("M2 scan/coordinate raw-output hashes disagree")
    hashes = {
        "metadata.csv": _verify_sha(paths["metadata"], handoff["metadata_sha256"]),
        "sequences_for_clustering.fasta": _verify_sha(paths["fasta"], handoff["fasta_sha256"]),
        "split_manifest.csv": _verify_sha(paths["split_manifest"], handoff["split_manifest_sha256"]),
        "per_sequence_status.csv": _verify_sha(paths["status"], raw_hashes["per_sequence_status.csv"]),
        "hmm_coordinate_map.csv": _verify_sha(paths["coordinates"], raw_hashes["hmm_coordinate_map.csv"]),
        "7r1c_residue_hmm_map.csv": _verify_sha(paths["structure"], structure["coordinate_map_sha256"]),
    }
    index = load_coordinate_index(records, paths["status"], paths["coordinates"],
                                  paths["structure"], hmm_length=39)
    counts = Counter((row["split"], row["primary_analysis_eligible"]) for row in records)
    if (len(records) != 2076 or counts[("discovery", "true")] != 1202
            or counts[("validation", "true")] != 519
            or len(index.reference_by_state) != 39
            or coordinate.get("n_coordinate_rows") != 80219):
        raise ValueError("Frozen M2 coordinate population/structure count differs")
    return index, {
        "dataset_version": handoff["dataset_version"],
        "split_version": handoff["split_version"],
        "m2_scan_run_commit": scan["git"]["commit"],
        "m2_scan_receipt_sha256": digest(paths["scan_receipt"]),
        "m2_coordinate_receipt_sha256": digest(paths["coordinate_receipt"]),
        "m2_structure_receipt_sha256": digest(paths["structure_receipt"]),
        "input_and_raw_table_sha256": hashes,
        "n_qc_eligible_sequences": len(records),
        "n_primary_discovery": counts[("discovery", "true")],
        "n_primary_validation": counts[("validation", "true")],
        "n_m2_coordinate_rows": coordinate["n_coordinate_rows"],
        "n_7r1c_deposited_residues": structure["n_deposited_residues"],
        "n_7r1c_modeled_residues": structure["n_modeled_residues"],
        "n_7r1c_unmodeled_residues": structure["n_unmodeled_residues"],
        "n_7r1c_hmm_match_states": len(index.reference_by_state),
    }


def require_validation_freeze(split: str, freeze_record: Path | None) -> str | None:
    """Require A's machine-readable Markdown front matter before validation.

    This checks that the record *claims* a frozen candidate/rule artifact. A
    still owns the temporal audit proving it preceded validation inspection.
    """
    if split == "validation":
        if freeze_record is None or not freeze_record.is_file() or freeze_record.stat().st_size == 0:
            raise ValueError("Validation mapping requires A's nonempty candidate freeze record")
        content = freeze_record.read_text(encoding="utf-8-sig")
        front = re.match(r"\A---\s*\n(.*?)\n---(?:\s*\n|\s*\Z)", content, flags=re.DOTALL)
        if front is None:
            raise ValueError("A's freeze record needs Markdown front matter")
        fields = {}
        for line in front.group(1).splitlines():
            if not line.strip():
                continue
            if ":" not in line:
                raise ValueError("Malformed freeze record front matter")
            key, value = (part.strip() for part in line.split(":", 1))
            if key in fields or not key:
                raise ValueError("Duplicate or empty freeze record field")
            fields[key] = value
        if fields.get("status") != "frozen":
            raise ValueError("Candidate freeze status is not frozen")
        for key in ("candidate_manifest_sha256", "evaluation_rules_sha256"):
            if not re.fullmatch(r"[0-9a-f]{64}", fields.get(key, "")):
                raise ValueError(f"Candidate freeze lacks valid {key}")
        try:
            frozen_at = datetime.fromisoformat(fields.get("frozen_at_utc", "").replace("Z", "+00:00"))
        except ValueError as error:
            raise ValueError("Candidate freeze lacks valid frozen_at_utc") from error
        if frozen_at.tzinfo is None or frozen_at.utcoffset() != timezone.utc.utcoffset(frozen_at):
            raise ValueError("Candidate freeze time must be UTC with timezone")
        return digest(freeze_record)
    if freeze_record is not None:
        raise ValueError("A validation freeze record is only valid in validation mode")
    return None


def read_windows(path: Path, index, *, split: str, cohort: str) -> list[dict]:
    """Require explicit 1-based coordinates and the requested frozen queue."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)) or not set(WINDOW_FIELDS) <= set(reader.fieldnames):
            raise ValueError("Windows CSV needs unique internal_id,start_1based,end_1based_inclusive headers")
        rows = list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError("Windows CSV is empty or malformed")
    seen_ids, seen_coordinates = set(), set()
    for row in rows:
        identifier = row["internal_id"]
        if identifier not in index.records:
            raise ValueError(f"Unknown window internal_id: {identifier}")
        record = index.records[identifier]
        if record["split"] != split:
            raise ValueError(f"Window {identifier} belongs to {record['split']}, not {split}")
        if cohort == "primary" and not record["primary_analysis_eligible"]:
            raise ValueError(f"Window {identifier} is outside the primary queue")
        if cohort == "sensitivity" and record["primary_analysis_eligible"]:
            raise ValueError(f"Window {identifier} is not a sensitivity record")
        try:
            start = int(row["start_1based"])
            end = int(row["end_1based_inclusive"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"Noninteger window boundary for {identifier}") from error
        if str(start) != row["start_1based"] or str(end) != row["end_1based_inclusive"]:
            raise ValueError(f"Noncanonical integer window boundary for {identifier}")
        if start < 1 or end < start or end > len(record["sequence"]):
            raise ValueError(f"Out-of-range window for {identifier}")
        key = (identifier, start, end)
        if key in seen_coordinates:
            raise ValueError(f"Duplicate window coordinates: {key}")
        seen_coordinates.add(key)
        row["start_1based"], row["end_1based_inclusive"] = start, end
        row["window_id"] = row.get("window_id") or f"{identifier}:{start}-{end}"
        if row["window_id"] in seen_ids:
            raise ValueError(f"Duplicate window_id: {row['window_id']}")
        seen_ids.add(row["window_id"])
    return rows


def _write_csv(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def map_requests(index, requests: list[dict], min_match_fraction: float) -> tuple[list[dict], list[dict]]:
    window_rows, residue_rows = [], []
    for request in requests:
        mapped = index.map_window(request["internal_id"], request["start_1based"],
                                  request["end_1based_inclusive"],
                                  min_match_fraction=min_match_fraction)
        window_id = request["window_id"]
        window_rows.append({
            **{key: mapped[key] for key in WINDOW_COLUMNS if key in mapped},
            "window_id": window_id,
            "hmm_accession": "PF00741.24", "7r1c_author_chain": "N",
            "n_reference_modeled": sum(item["reference_status"] == "reference_modeled"
                                       for item in mapped["residues"]),
            "n_reference_unmodeled": sum(item["reference_status"] == "reference_unmodeled"
                                         for item in mapped["residues"]),
            "domain_numbers_json": json.dumps(mapped["domain_numbers"]),
            "hmm_match_states_json": json.dumps(mapped["hmm_match_states"]),
            "hmm_deleted_states_internal_json": json.dumps(mapped["hmm_deleted_states_internal"]),
        })
        for residue in mapped["residues"]:
            reference = residue["reference_7r1c"] or {}
            residue_rows.append({
                "window_id": window_id, "internal_id": mapped["internal_id"],
                "raw_position": residue["raw_position"], "residue": residue["residue"],
                "mapping_status": residue["mapping_status"],
                "domain_number": residue["domain_number"],
                "hmm_match_state": residue["hmm_match_state"],
                "coordinate_options_json": json.dumps(residue["coordinate_options"]),
                "reference_status": residue["reference_status"],
                "hmm_accession": "PF00741.24", "7r1c_author_chain": "N",
                "7r1c_deposited_position": reference.get("deposited_position"),
                "7r1c_author_residue_number": reference.get("pdb_residue_number"),
                "7r1c_author_insertion_code": reference.get("pdb_insertion_code"),
                "7r1c_reference_residue": reference.get("reference_residue"),
                "7r1c_modeled": reference.get("modeled"),
                "7r1c_secondary_structure": reference.get("secondary_structure"),
            })
    return window_rows, residue_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in DEFAULTS.items():
        parser.add_argument("--" + name.replace("_", "-"), type=Path, default=ROOT / default)
    parser.add_argument("--windows", type=Path, help="CSV with internal_id and 1-based inclusive start/end")
    parser.add_argument("--split", choices=("discovery", "validation"), default="discovery")
    parser.add_argument("--cohort", choices=("primary", "sensitivity", "all"), default="primary")
    parser.add_argument("--validation-freeze", type=Path,
                        help="A-owned candidate freeze record, required for validation windows")
    parser.add_argument("--min-match-fraction", type=float, default=0.8,
                        help="Technical low-coverage flag only; not a candidate filter")
    parser.add_argument("--out-dir", type=Path,
                        help="New shared/ignored output directory; required with --windows")
    parser.add_argument("--summary-out", type=Path,
                        help="Default: results/M3/coordinate_summary.json without windows")
    args = parser.parse_args()
    if args.windows is None and (args.out_dir is not None or args.validation_freeze is not None):
        parser.error("--out-dir/--validation-freeze require --windows")
    if args.windows is not None and args.out_dir is None:
        parser.error("--windows requires --out-dir")
    freeze_sha = require_validation_freeze(args.split, args.validation_freeze)
    paths = {name: getattr(args, name) for name in DEFAULTS}
    index, summary = load_verified_index(paths)
    summary.update({
        "schema_version": "1.0", "run_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "M3 C prelaunch coordinate handoff, not ESM discovery or validation",
        "coordinate_convention": "natural/window 1-based inclusive; HMM match states 1-39; 7R1C author chain N",
        "coverage_policy": "Unaligned/inserted/ambiguous and unmodeled reference positions remain explicit; min_match_fraction only flags technical coverage",
        "implementation_sha256": digest(Path(__file__)),
        "mapping_module_sha256": digest(ROOT / "bioinformatics/m3_window_coordinates.py"),
        "validation_freeze_sha256": freeze_sha,
        "candidate_selection_performed": False,
        "validation_used_for_parameter_selection": False,
    })
    if args.windows is None:
        summary["window_mapping_status"] = "not_run_prelaunch"
        destination = args.summary_out or ROOT / "results/M3/coordinate_summary.json"
        if destination.exists():
            raise FileExistsError(f"Refusing to overwrite coordinate summary: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(destination)
        return
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty mapping directory: {args.out_dir}")
    summary_path = args.summary_out or args.out_dir / "mapping_receipt.json"
    if summary_path.exists():
        raise FileExistsError(f"Refusing to overwrite mapping receipt: {summary_path}")
    requests = read_windows(args.windows, index, split=args.split, cohort=args.cohort)
    window_rows, residue_rows = map_requests(index, requests, args.min_match_fraction)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    window_path = args.out_dir / "window_coordinate_map.csv"
    residue_path = args.out_dir / "residue_coordinate_map.csv"
    _write_csv(window_path, WINDOW_COLUMNS, window_rows)
    _write_csv(residue_path, RESIDUE_COLUMNS, residue_rows)
    summary.update({
        "window_mapping_status": "completed_coordinate_only",
        "requested_split": args.split, "requested_cohort": args.cohort,
        "technical_min_match_fraction": args.min_match_fraction,
        "windows_input_sha256": digest(args.windows),
        "n_windows": len(window_rows), "n_residue_rows": len(residue_rows),
        "coverage_status_counts": dict(Counter(row["coverage_status"] for row in window_rows)),
        "output_sha256": {window_path.name: digest(window_path),
                          residue_path.name: digest(residue_path)},
    })
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(summary_path)


if __name__ == "__main__":
    main()
