import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import validate_m2_release as release


AUDIT = release.DEFAULT_AUDIT
MMSEQS = release.DEFAULT_MMSEQS
SPLIT = release.DEFAULT_SPLIT


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


class ReleaseValidationTest(unittest.TestCase):
    def test_frozen_receipts_pass(self):
        self.assertEqual(release.validate_release(AUDIT, MMSEQS, SPLIT), [])

    def test_missing_receipt_fails(self):
        missing = ROOT / "results" / "data_audit" / "does_not_exist.json"
        errors = release.validate_release(missing, MMSEQS, SPLIT)
        self.assertTrue(any(item.startswith("missing receipt:") for item in errors))

    def test_missing_receipt_field_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "split.json"
            payload = json.loads(SPLIT.read_text(encoding="utf-8"))
            del payload["split_version"]
            write_json(path, payload)
            errors = release.validate_release(AUDIT, MMSEQS, path)
        self.assertIn("missing receipt field: split.split_version", errors)

    def test_missing_hash_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "audit.json"
            payload = json.loads(AUDIT.read_text(encoding="utf-8"))
            del payload["output_sha256"]["metadata.csv"]
            write_json(path, payload)
            errors = release.validate_release(path, MMSEQS, SPLIT)
        self.assertIn("missing hash: audit.output_sha256.metadata.csv", errors)

    def test_empty_hash_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mmseqs.json"
            payload = json.loads(MMSEQS.read_text(encoding="utf-8"))
            payload["mmseqs_cluster_tsv_sha256"] = ""
            write_json(path, payload)
            errors = release.validate_release(AUDIT, path, SPLIT)
        self.assertTrue(any(item.startswith("missing hash: mmseqs.mmseqs_cluster_tsv_sha256") for item in errors))

    def test_cluster_leakage_flag_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "split.json"
            payload = json.loads(SPLIT.read_text(encoding="utf-8"))
            payload["cluster_leakage"] = True
            write_json(path, payload)
            errors = release.validate_release(AUDIT, MMSEQS, path)
        self.assertIn("cluster leakage: split receipt cluster_leakage is not false", errors)

    def test_manifest_cluster_leakage_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "split_manifest.csv"
            manifest.write_text(
                "internal_id,homology_cluster,split\n"
                "GVPA_000001,cluster-a,discovery\n"
                "GVPA_000002,cluster-a,validation\n",
                encoding="utf-8",
            )
            errors = release.validate_release(AUDIT, MMSEQS, SPLIT, manifest)
        self.assertIn(
            "cluster leakage: homology cluster cluster-a is assigned to more than one split",
            errors,
        )

    def test_clean_manifest_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "split_manifest.csv"
            manifest.write_text(
                "internal_id,homology_cluster,split\n"
                "GVPA_000001,cluster-a,discovery\n"
                "GVPA_000002,cluster-b,validation\n",
                encoding="utf-8",
            )
            self.assertEqual(release.validate_release(AUDIT, MMSEQS, SPLIT, manifest), [])

    def test_command_exits_nonzero_when_receipt_missing(self):
        completed = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "validate_m2_release.py"), "--audit", "missing-audit.json"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 1)
        self.assertIn("missing receipt:", completed.stderr)


if __name__ == "__main__":
    unittest.main()
