"""Assistant/read-only byte-level check of M2 receipts and available local files.

No download, shared-access verification, or independent-review signature is implied.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=False)


def audit(pdb_candidate: Path | None = None) -> dict:
    scan = json.loads((ROOT / "results/bioinformatics/pf00741_scan_summary.json").read_text(encoding="utf-8"))
    structure = json.loads((ROOT / "results/bioinformatics/structure_7r1c_summary.json").read_text(encoding="utf-8"))
    rows = []
    def check(label: str, path: Path, expected: str, source: str) -> None:
        found = path.is_file()
        actual = digest(path.read_bytes()) if found else None
        row = {"artifact": label, "location_checked": path.as_posix(), "expected_sha256": expected,
               "actual_sha256": actual, "expected_source": source,
               "status": "match" if actual == expected else "mismatch" if found else "unverified_missing",
               "shared_uri": None, "shared_access": "unverified; no location/access procedure supplied"}
        if found and actual != expected:
            normalized = digest(path.read_bytes().replace(b"\r\n", b"\n"))
            row["lf_normalized_sha256_diagnostic_only"] = normalized
            row["lf_matches_expected"] = normalized == expected
            crlf = digest(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
            row["crlf_normalized_sha256_diagnostic_only"] = crlf
            row["crlf_matches_expected"] = crlf == expected
        rows.append(row)

    command = scan["command"]
    for flag, key in [("--metadata", "metadata_sha256"), ("--fasta", "fasta_sha256"),
                      ("--split-manifest", "split_manifest_sha256"), ("--audit-summary", "audit_summary_sha256"),
                      ("--split-summary", "split_summary_sha256")]:
        path = Path(command[command.index(flag) + 1])
        check(path.as_posix(), path, scan["input"][key], "pf00741_scan_summary.input." + key)
    out_dir = Path(command[command.index("--out-dir") + 1])
    for name, value in scan["raw_output_sha256"].items():
        check(name, out_dir / name, value, "pf00741_scan_summary.raw_output_sha256")
    for receipt, source in [(scan, "pf00741_scan_summary"), (structure, "structure_7r1c_summary")]:
        cmd = receipt["command"]
        for flag, expected in [("--hmm", receipt["hmm"]["sha256"]),
                               ("--hmm-receipt", receipt["hmm_download_receipt_sha256"])]:
            check(source + flag, Path(cmd[cmd.index(flag) + 1]), expected, source)
    cmd = structure["command"]
    check("PDB recorded location", Path(cmd[cmd.index("--pdb") + 1]), structure["pdb_sha256"], "structure_7r1c_summary")
    check("PDB download receipt", Path(cmd[cmd.index("--pdb-receipt") + 1]), structure["pdb_download_receipt_sha256"], "structure_7r1c_summary")
    check("7R1C residue map", Path(cmd[cmd.index("--out-dir") + 1]) / structure["coordinate_map_filename"],
          structure["coordinate_map_sha256"], "structure_7r1c_summary")
    if pdb_candidate:
        check("PDB alternative local snapshot; not the recorded/shared location", pdb_candidate,
              structure["pdb_sha256"], "structure_7r1c_summary")
    implementation = []
    base = scan["git"]["commit"]
    for key, path in [("workflow", "bioinformatics/m2_pf00741.py"), ("coordinates", "bioinformatics/m2_coordinates.py")]:
        at_base = git("show", base + ":" + path)
        current_blob = git("show", "HEAD:" + path)
        current_bytes = (ROOT / path).read_bytes()
        implementation.append({"path": path, "expected_sha256": scan["implementation_sha256"][key],
                               "recorded_base_contains_file": at_base.returncode == 0,
                               "base_blob_sha256": digest(at_base.stdout) if at_base.returncode == 0 else None,
                               "current_blob_sha256": digest(current_blob.stdout),
                               "current_worktree_sha256": digest(current_bytes),
                               "worktree_matches_receipt": digest(current_bytes) == scan["implementation_sha256"][key],
                               "blob_matches_receipt": digest(current_blob.stdout) == scan["implementation_sha256"][key]})
    return {"kind": "assistant_technical_precheck_not_independent_review", "checked_utc": datetime.now(timezone.utc).isoformat(),
            "checked_sha": git("rev-parse", "HEAD").stdout.decode().strip(),
            "working_tree_dirty_at_check": bool(git("status", "--porcelain").stdout.strip()),
            "command": ["python", *sys.argv], "python": platform.python_version(),
            "artifacts": rows, "scan_git_receipt": scan["git"],
            "recorded_scan_commit_exists": git("cat-file", "-e", base + "^{commit}").returncode == 0,
            "implementation_checks": implementation, "structure_git_receipt": structure["git"],
            "dirty_run_provenance": "unverified: no original dirty patch/full working tree or environment lock supplied",
            "shared_storage_access": "unverified: only local/shared placeholder in receipt; Issue 15 has no access instructions",
            "limits": ["Missing files are not verified by copying receipt hashes", "LF-normalized equality is diagnostic, not byte equality",
                       "A matching alternative PDB proves file bytes only, not shared access or scan reproducibility"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-candidate", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.pdb_candidate)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    from collections import Counter
    print(dict(Counter(row["status"] for row in result["artifacts"])))


if __name__ == "__main__":
    main()
