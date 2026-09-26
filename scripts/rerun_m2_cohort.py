"""Regenerate only the M2 cohort figure from a clean commit, with provenance."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import m2_visual_summary as visual


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(metadata_path: Path, audit_path: Path, run_id: str) -> Path:
    if not re.fullmatch(r"M2-E-FIGURES-[0-9]{3}", run_id):
        raise ValueError("Use a new M2-E-FIGURES-NNN run ID")
    directory = ROOT / "reports/M2/runs" / run_id
    if directory.exists():
        raise ValueError("Run ID already exists; choose a new one")
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True)
    if status.strip():
        raise ValueError("Commit code first: working tree must be clean before the run")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if sha256(metadata_path) != audit["output_sha256"]["metadata.csv"]:
        raise ValueError("Metadata SHA256 differs from frozen audit")
    metadata = visual.read_csv(metadata_path)
    counts = dict(visual.cohort_counts(metadata))
    if counts != audit["cohort_counts"] or len(metadata) != audit["candidate_records"]:
        raise ValueError("All cohort counts must match frozen audit")
    manifest_path = ROOT / "reports/M2/figure_manifest.csv"
    manifest = visual.read_csv(manifest_path)
    entries = [row for row in manifest if row["figure"] == "cohort_counts.png"]
    if len(entries) != 1:
        raise ValueError("Expected exactly one cohort figure manifest row")
    registry_path = ROOT / "experiments/registry.csv"
    registry = visual.read_csv(registry_path)
    if any(row["experiment_id"] == run_id for row in registry):
        raise ValueError("Run ID already registered")
    figure = ROOT / "reports/M2/figures/cohort_counts.png"
    previous_sha256 = sha256(figure)
    visual.configure_style()
    rows = visual.render_cohort_figure(metadata, figure)
    directory.mkdir(parents=True)
    visual.write_csv(directory / "cohort_counts.csv", rows)
    receipt_path = directory / "run_receipt.json"
    for row in manifest:
        for field in ("git_commit", "output_sha256", "input_metadata_sha256", "receipt"):
            row.setdefault(field, "")
    entries[0].update(experiment_id=run_id, script="scripts/rerun_m2_cohort.py",
                      git_commit=commit, output_sha256=sha256(figure),
                      input_metadata_sha256=sha256(metadata_path),
                      receipt=receipt_path.relative_to(ROOT).as_posix())
    visual.write_csv(manifest_path, manifest)
    # Read back persisted artifacts: verify source, exported counts and manifest agree.
    exported = visual.read_csv(directory / "cohort_counts.csv")
    if {row["analysis_cohort"]: int(row["count"]) for row in exported} != counts:
        raise ValueError("Exported counts disagree with metadata")
    manifest_row = next(row for row in visual.read_csv(manifest_path) if row["figure"] == figure.name)
    if manifest_row["output_sha256"] != sha256(figure) or manifest_row["git_commit"] != commit:
        raise ValueError("Figure manifest provenance mismatch")
    receipt = {
        "experiment_id": run_id, "executor": "Codex engineering assistant; not independent B/C review",
        "run_utc": datetime.now(timezone.utc).isoformat(), "git_commit": commit, "dirty_at_start": False,
        "command": ["python", *sys.argv], "dataset_version": audit["dataset_version"],
        "split_version": visual.SPLIT_VERSION, "seed": None,
        "software": {"python": platform.python_version(), "platform": platform.platform(),
                     "matplotlib": visual.matplotlib.__version__},
        "input_sha256": {"metadata.csv": sha256(metadata_path), "audit_summary.json": sha256(audit_path)},
        "implementation_sha256": {"rerun_m2_cohort.py": sha256(Path(__file__)),
                                  "m2_visual_summary.py": sha256(ROOT / "scripts/m2_visual_summary.py")},
        "cohort_counts": counts, "previous_figure_sha256": previous_sha256,
        "output_sha256": {p.relative_to(ROOT).as_posix(): sha256(p)
                          for p in (figure, directory / "cohort_counts.csv", manifest_path)},
        "checks": {"all_frozen_cohorts_match": True, "rendered_bars_match_source": True,
                   "exported_counts_match_source": True, "figure_manifest_hash_matches": True},
        "limitations": ["Only cohort_counts.png regenerated; other six figures not revalidated",
                        "Old M2-E-FIGURES-001 execution commit remains unknown",
                        "Assistant precheck is not independent review or A acceptance"],
    }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    registry.append({"experiment_id": run_id, "owner": "Codex-assistant", "status": "completed_pending_independent_review",
                     "rq": "M2-figures", "dataset_version": audit["dataset_version"], "split_version": visual.SPLIT_VERSION,
                     "config_path": "results/data_audit/gvpa_v1_audit_summary.json", "seed": "",
                     "git_commit": commit, "metrics_path": (directory / "cohort_counts.csv").relative_to(ROOT).as_posix(),
                     "artifact_uri": directory.relative_to(ROOT).as_posix(),
                     "conclusion": "Cohort figure regenerated from verified metadata: 1721;20;305;30;2",
                     "limitations": "Assistant run; only cohort figure rerun; independent review pending; receipt=" + receipt_path.relative_to(ROOT).as_posix()})
    visual.write_csv(registry_path, registry)
    return receipt_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--audit-summary", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    print(run(args.metadata, args.audit_summary, args.run_id).relative_to(ROOT))


if __name__ == "__main__":
    main()
