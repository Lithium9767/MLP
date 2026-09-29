import csv
import tempfile
import unittest
from pathlib import Path

from scripts.m3_data_handoff import FIELDS, verify_embedding_ids, verify_manifests, write_manifests


class HandoffChecksTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.row = dict(zip(FIELDS, ["GVPA_000001", "seq1", "a" * 64, "100", "primary", "discovery", "GVPA_000001"]))
        self.groups = {"discovery_primary": [self.row], "validation_primary": [], "sensitivity": []}
        write_manifests(self.groups, self.root)

    def test_manifest_tampering_is_rejected(self):
        verify_manifests(self.groups, self.root)
        path = self.root / "discovery_primary_manifest.csv"
        text = path.read_text().replace("a" * 64, "b" * 64)
        path.write_text(text)
        with self.assertRaisesRegex(ValueError, "content differs"):
            verify_manifests(self.groups, self.root)

    def test_embedding_missing_duplicate_and_wrong_hash_rejected(self):
        path = self.root / "embeddings.csv"
        fields = ["internal_id", "sequence_sha256", "analysis_cohort", "split"]

        def write(rows):
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fields)
                writer.writeheader()
                writer.writerows(rows)

        valid = {field: self.row[field] for field in fields}
        write([valid])
        verify_embedding_ids(path, self.groups, "discovery_primary")
        write([valid, valid])
        with self.assertRaisesRegex(ValueError, "duplicate internal_id"):
            verify_embedding_ids(path, self.groups, "discovery_primary")
        write([])
        with self.assertRaisesRegex(ValueError, "missing 1"):
            verify_embedding_ids(path, self.groups, "discovery_primary")
        write([{**valid, "sequence_sha256": "b" * 64}])
        with self.assertRaisesRegex(ValueError, "sequence_sha256 differs"):
            verify_embedding_ids(path, self.groups, "discovery_primary")
        write([{**valid, "internal_id": "unknown"}])
        with self.assertRaisesRegex(ValueError, "extra 1"):
            verify_embedding_ids(path, self.groups, "discovery_primary")


if __name__ == "__main__":
    unittest.main()
