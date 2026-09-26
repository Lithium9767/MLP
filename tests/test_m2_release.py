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

    def test_b_reproduction_counts_match_and_record_mmseqs_version(self):
        audit = json.loads(AUDIT.read_text(encoding="utf-8"))
        frozen_split = json.loads(SPLIT.read_text(encoding="utf-8"))
        data = {
            "dataset_version": release.FROZEN["dataset_version"],
            "status": "completed_with_version_note",
            "counts": {
                "candidate_records": 2078,
                "sequence_qc_eligible": 2076,
                "primary": 1721,
                "homology_clusters": 478,
            },
            "outputs": {
                "metadata_sha256": audit["output_sha256"]["metadata.csv"],
                "members_sha256": audit["output_sha256"]["members.csv"],
                "clustering_fasta_sha256": audit["output_sha256"]["sequences_for_clustering.fasta"],
            },
            "source": {"source_member_sha256": audit["source"]["source_member_sha256"], "archive_sha256": None},
        }
        split = {
            "dataset_version": release.FROZEN["dataset_version"],
            "status": "completed_with_version_note",
            "cluster_leakage": False,
            "sequence_counts": {"discovery": 1453, "validation": 623},
            "cluster_counts": {"discovery": 335, "validation": 143},
            "split_manifest_sha256": frozen_split["split_manifest_sha256"],
            "mmseqs2_version": "18.8cc5c",
            "comparison": {
                "counts_match_frozen": True,
                "split_manifest_hash_matches_frozen": True,
                "cluster_tsv_hash_matches_frozen": False,
            },
        }
        errors, notes = release.check_b_reproduction(data, split, audit, frozen_split)
        self.assertEqual(errors, [])
        self.assertTrue(any("18.8cc5c" in note for note in notes))

    def test_b_cluster_leakage_fails(self):
        audit = json.loads(AUDIT.read_text(encoding="utf-8"))
        frozen_split = json.loads(SPLIT.read_text(encoding="utf-8"))
        data = {
            "dataset_version": release.FROZEN["dataset_version"],
            "status": "completed",
            "counts": {"candidate_records": 2078, "sequence_qc_eligible": 2076, "primary": 1721, "homology_clusters": 478},
            "outputs": {
                "metadata_sha256": audit["output_sha256"]["metadata.csv"],
                "members_sha256": audit["output_sha256"]["members.csv"],
                "clustering_fasta_sha256": audit["output_sha256"]["sequences_for_clustering.fasta"],
            },
            "source": {"source_member_sha256": audit["source"]["source_member_sha256"], "archive_sha256": "a" * 64},
        }
        split = {
            "dataset_version": release.FROZEN["dataset_version"],
            "status": "completed",
            "cluster_leakage": True,
            "sequence_counts": {"discovery": 1453, "validation": 623},
            "cluster_counts": {"discovery": 335, "validation": 143},
            "split_manifest_sha256": frozen_split["split_manifest_sha256"],
            "mmseqs2_version": "15-6f452",
            "comparison": {"counts_match_frozen": True, "split_manifest_hash_matches_frozen": True, "cluster_tsv_hash_matches_frozen": True},
        }
        errors, _notes = release.check_b_reproduction(data, split, audit, frozen_split)
        self.assertTrue(any(item.startswith("cluster leakage:") for item in errors))

    def test_missing_pf00741_scan_receipt_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            errors, _notes = release.validate_handoff(
                AUDIT,
                SPLIT,
                AUDIT,
                root / "missing-scan.json",
                root / "missing-coordinates.json",
            )
        self.assertTrue(any(item.startswith("missing receipt:") for item in errors))

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
