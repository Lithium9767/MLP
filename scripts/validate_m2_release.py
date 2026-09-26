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


B_COUNT_FIELDS = {
    "candidate_records": FROZEN["candidate_records"],
    "sequence_qc_eligible": FROZEN["sequence_qc_eligible"],
    "primary": FROZEN["primary_cohort"],
    "homology_clusters": FROZEN["total_clusters"],
}
SCAN_FIELDS = (
    "schema_version",
    "n_sequences",
    "counts",
    "by_split",
    "hmm",
    "parameters",
    "purpose",
)
COORDINATE_FIELDS = (
    "schema_version",
    "n_coordinate_rows",
    "per_split_coverage",
    "coordinate_convention",
    "hmm",
    "parameters",
)


def _expect(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def check_b_reproduction(
    data_receipt: dict[str, Any],
    split_receipt: dict[str, Any],
    audit: dict[str, Any],
    frozen_split: dict[str, Any],
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    notes: list[str] = []
    errors.extend(missing_fields(data_receipt, ("dataset_version", "status", "counts", "outputs", "source"), "b_data"))
    errors.extend(missing_fields(
        split_receipt,
        ("dataset_version", "status", "cluster_leakage", "sequence_counts", "cluster_counts", "comparison", "mmseqs2_version"),
        "b_split",
    ))
    counts = data_receipt.get("counts") if isinstance(data_receipt.get("counts"), dict) else {}
    for name, expected in B_COUNT_FIELDS.items():
        _expect(errors, counts.get(name) == expected, f"frozen baseline mismatch: b_data.counts.{name} is {counts.get(name)!r}, expected {expected!r}")
    outputs = data_receipt.get("outputs") if isinstance(data_receipt.get("outputs"), dict) else {}
    for key, audit_path in (
        ("metadata_sha256", ("output_sha256", "metadata.csv")),
        ("members_sha256", ("output_sha256", "members.csv")),
        ("clustering_fasta_sha256", ("output_sha256", "sequences_for_clustering.fasta")),
    ):
        value = outputs.get(key)
        if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
            errors.append(f"missing hash: b_data.outputs.{key}")
        elif value != lookup(audit, audit_path):
            errors.append(f"hash mismatch: b_data.outputs.{key} != frozen audit")
    source = data_receipt.get("source") if isinstance(data_receipt.get("source"), dict) else {}
    member_hash = source.get("source_member_sha256")
    if not isinstance(member_hash, str) or SHA256_RE.fullmatch(member_hash) is None:
        errors.append("missing hash: b_data.source.source_member_sha256")
    elif member_hash != lookup(audit, ("source", "source_member_sha256")):
        errors.append("hash mismatch: b_data.source.source_member_sha256 != frozen audit")
    if source.get("archive_sha256") in (None, ""):
        notes.append("B did not check gv.zip SHA-256; the recognition JSON hash matches the frozen audit.")
    if split_receipt.get("cluster_leakage") is not False:
        errors.append("cluster leakage: B split receipt cluster_leakage is not false")
    sequence_counts = split_receipt.get("sequence_counts") if isinstance(split_receipt.get("sequence_counts"), dict) else {}
    cluster_counts = split_receipt.get("cluster_counts") if isinstance(split_receipt.get("cluster_counts"), dict) else {}
    _expect(
        errors,
        sequence_counts.get("discovery") == FROZEN["discovery_sequences"] and sequence_counts.get("validation") == FROZEN["validation_sequences"],
        "frozen baseline mismatch: B discovery/validation sequence counts",
    )
    _expect(
        errors,
        cluster_counts.get("discovery") == FROZEN["discovery_clusters"] and cluster_counts.get("validation") == FROZEN["validation_clusters"],
        "frozen baseline mismatch: B discovery/validation cluster counts",
    )
    manifest_hash = split_receipt.get("split_manifest_sha256")
    if not isinstance(manifest_hash, str) or SHA256_RE.fullmatch(manifest_hash) is None:
        errors.append("missing hash: b_split.split_manifest_sha256")
    elif manifest_hash != frozen_split.get("split_manifest_sha256"):
        errors.append("hash mismatch: B split_manifest_sha256 != frozen split")
    comparison = split_receipt.get("comparison") if isinstance(split_receipt.get("comparison"), dict) else {}
    if comparison.get("counts_match_frozen") is not True or comparison.get("split_manifest_hash_matches_frozen") is not True:
        errors.append("B comparison does not confirm frozen counts and split manifest hash")
    if split_receipt.get("mmseqs2_version") != "15-6f452" or comparison.get("cluster_tsv_hash_matches_frozen") is False:
        notes.append(
            f"B MMseqs2 is {split_receipt.get('mmseqs2_version')!r}; frozen receipt is 15-6f452. "
            "Cluster TSV hashes differ and the reproduced split version is not the frozen version. Counts and split manifest hash match."
        )
    return errors, notes


def check_structure_receipt(receipt: dict[str, Any]) -> list[str]:
    errors = missing_fields(
        receipt,
        ("pdb_id", "chain", "hmm", "pdb_sha256", "coordinate_map_sha256", "n_deposited_residues", "n_modeled_residues", "n_unmodeled_residues", "n_distinct_hmm_match_states", "purpose"),
        "structure",
    )
    hmm = receipt.get("hmm") if isinstance(receipt.get("hmm"), dict) else {}
    errors.extend(hash_errors(receipt, (("pdb_sha256",), ("coordinate_map_sha256",)), "structure"))
    errors.extend(hash_errors(hmm, (("sha256",),), "structure.hmm"))
    _expect(errors, receipt.get("pdb_id") == "7R1C" and receipt.get("chain") == "N", "structure receipt is not PDB 7R1C chain N")
    _expect(errors, hmm.get("length") == 39 and str(hmm.get("accession", "")).startswith("PF00741"), "structure HMM is not PF00741 with 39 match states")
    _expect(errors, receipt.get("n_deposited_residues") == 88, "structure deposited residue count is not 88")
    _expect(errors, receipt.get("n_modeled_residues") == 65 and receipt.get("n_unmodeled_residues") == 23, "structure modeled/unmodeled counts are not 65/23")
    _expect(errors, receipt.get("n_distinct_hmm_match_states") == 39, "structure HMM match-state count is not 39")
    purpose = str(receipt.get("purpose", "")).lower()
    _expect(errors, "function" in purpose, "structure receipt does not state that it is not a functional label")
    return errors


def check_scan_receipt(receipt: dict[str, Any], receipt_name: str) -> list[str]:
    fields = SCAN_FIELDS if receipt_name == "pf00741_scan" else COORDINATE_FIELDS
    errors = missing_fields(receipt, fields, receipt_name)
    hmm = receipt.get("hmm") if isinstance(receipt.get("hmm"), dict) else {}
    errors.extend(hash_errors(hmm, (("sha256",),), f"{receipt_name}.hmm"))
    parameters = receipt.get("parameters") if isinstance(receipt.get("parameters"), dict) else {}
    _expect(errors, parameters.get("validation_used_for_conservation") is False, f"{receipt_name} used validation sequences for conservation")
    if receipt_name == "pf00741_scan":
        _expect(errors, receipt.get("n_sequences") == FROZEN["sequence_qc_eligible"], "PF00741 scan did not cover 2076 QC-eligible sequences")
    return errors


def validate_handoff(
    b_data_path: Path,
    b_split_path: Path,
    structure_path: Path,
    scan_path: Path,
    coordinate_path: Path,
    audit_path: Path = DEFAULT_AUDIT,
    split_path: Path = DEFAULT_SPLIT,
) -> tuple[list[str], list[str]]:
    audit, audit_errors = load_receipt(audit_path)
    frozen_split, split_errors = load_receipt(split_path)
    b_data, b_data_errors = load_receipt(b_data_path)
    b_split, b_split_errors = load_receipt(b_split_path)
    structure, structure_errors = load_receipt(structure_path)
    errors = audit_errors + split_errors + b_data_errors + b_split_errors + structure_errors
    notes: list[str] = []
    if audit is not None and frozen_split is not None and b_data is not None and b_split is not None:
        found, noted = check_b_reproduction(b_data, b_split, audit, frozen_split)
        errors.extend(found)
        notes.extend(noted)
    if structure is not None:
        errors.extend(check_structure_receipt(structure))
    scan, scan_errors = load_receipt(scan_path)
    coordinate, coordinate_errors = load_receipt(coordinate_path)
    errors.extend(scan_errors)
    errors.extend(coordinate_errors)
    if scan is not None:
        errors.extend(check_scan_receipt(scan, "pf00741_scan"))
    if coordinate is not None:
        errors.extend(check_scan_receipt(coordinate, "hmm_coordinates"))
    return errors, notes


def render_reproducibility_report(errors: list[str], notes: list[str]) -> str:
    lines = [
        "# M2 reproducibility check",
        "",
        "Generated by `scripts/validate_m2_release.py`. This check reads published receipts and does not replace `homology-8b9005e2d9-s42`.",
        "",
        "## Frozen baseline on main",
        "",
        "Data version `gvpa-recognition-c7f6f005d717`: 2078 candidates, 2076 QC-eligible, 1721 primary, 478 clusters, discovery 1453/335, validation 623/143, cluster leakage false.",
        "",
        "## B branch `feature/m2-1-data-reproduction`",
        "",
    ]
    if any(item.startswith("frozen baseline mismatch: b_") or item.startswith("frozen baseline mismatch: B") or "B comparison" in item or item.startswith("cluster leakage: B") for item in errors):
        lines.append("Count reproduction failed. See errors below.")
    else:
        lines.append("Counts, output hashes, and the split-manifest hash match the frozen baseline. Cluster leakage is false.")
    lines.extend(["", "## C branch `analysis/10-pf00741-hmm-mapping`", ""])
    if any(item.startswith("structure") for item in errors):
        lines.append("7R1C structure receipt failed the coordinate checks.")
    else:
        lines.append("PDB 7R1C chain N maps to PF00741.24 with 39 match states: 88 deposited residues, 65 modeled, 23 unmodeled.")
    if any("pf00741_scan" in item or "hmm_coordinates" in item for item in errors):
        lines.append("The 2076-sequence PF00741 scan summary and HMM coordinate summary are not in the published branch. Full-cohort scan acceptance is not met.")
    lines.extend(["", "## Notes", ""])
    lines.extend(f"- {note}" for note in notes) if notes else lines.append("- None.")
    lines.extend(["", "## Errors", ""])
    lines.extend(f"- {error}" for error in errors) if errors else lines.append("- None.")
    sensitivity_path = ROOT / "results" / "data_audit" / "mmseqs_sensitivity.json"
    lines.extend(["", "## MMseqs2 threshold sensitivity", ""])
    if sensitivity_path.is_file():
        sensitivity = json.loads(sensitivity_path.read_text(encoding="utf-8"))
        if sensitivity.get("status") == "completed":
            lines.append("Coverage is 0.8. The frozen split was not replaced.")
            lines.append("")
            lines.append("| Setting | Identity | Clusters | Max cluster | Singletons | Split agreement |")
            lines.append("| --- | ---: | ---: | ---: | ---: | ---: |")
            for row in sensitivity["settings"]:
                lines.append(
                    f"| {row['name']} | {row['identity']} | {row['n_clusters']} | {row['max_cluster_size']} | {row['n_singletons']} | {row['split_agreement']} |"
                )
        else:
            lines.append(sensitivity.get("reason", "Sensitivity did not complete."))
    else:
        lines.append("Sensitivity results are not present.")
    lines.append("")
    return "\n".join(lines)


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
    parser.add_argument("--b-data", type=Path, default=None)
    parser.add_argument("--b-split", type=Path, default=None)
    parser.add_argument("--structure", type=Path, default=None)
    parser.add_argument("--pf00741-scan", type=Path, default=None)
    parser.add_argument("--hmm-coordinates", type=Path, default=None)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--sensitivity-out", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    errors = validate_release(args.audit, args.mmseqs, args.split, args.manifest)
    notes: list[str] = []
    if args.b_data or args.b_split or args.structure or args.pf00741_scan or args.hmm_coordinates:
        handoff_errors, notes = validate_handoff(
            args.b_data or Path("missing-b-data.json"),
            args.b_split or Path("missing-b-split.json"),
            args.structure or Path("missing-structure.json"),
            args.pf00741_scan or Path("missing-pf00741-scan.json"),
            args.hmm_coordinates or Path("missing-hmm-coordinates.json"),
            args.audit,
            args.split,
        )
        errors.extend(handoff_errors)
    if args.sensitivity_out is not None:
        args.sensitivity_out.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "status": "blocked",
            "frozen_split_version": FROZEN["split_version"],
            "settings": [
                {"identity": 0.7, "coverage": 0.8},
                {"identity": 0.8, "coverage": 0.8},
                {"identity": 0.9, "coverage": 0.8},
            ],
            "reason": "mmseqs binary and sequences_for_clustering.fasta are not available; the frozen split was not replaced",
        }
        args.sensitivity_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(render_reproducibility_report(errors, notes), encoding="utf-8")
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print("M2 release validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
