import hashlib
import unittest
from unittest.mock import patch

from scripts import verify_m3_b_handoff as handoff


def row(identifier, sequence, cohort, split, cluster):
    return {
        "internal_id": identifier,
        "sequence_id": identifier,
        "sequence": sequence,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "sequence_length": str(len(sequence)),
        "sequence_qc_eligible": "true",
        "primary_analysis_eligible": str(cohort == "primary").lower(),
        "analysis_cohort": cohort,
        "split": split,
        "homology_cluster": cluster,
    }


class HandoffCheckTest(unittest.TestCase):
    def setUp(self):
        self.metadata = [row("a", "ACD", "primary", "discovery", "g1"),
                         row("b", "ACE", "primary", "validation", "g2"),
                         row("c", "ACF", "sensitivity_partial", "discovery", "g1")]
        self.split = [{key: value for key, value in r.items() if key not in
                       ("sequence", "sequence_sha256", "sequence_length", "sequence_qc_eligible")}
                      for r in self.metadata]
        self.manifests = {
            "discovery_primary_manifest.csv": [self.metadata[0].copy()],
            "validation_primary_manifest.csv": [self.metadata[1].copy()],
            "sensitivity_manifest.csv": [self.metadata[2].copy()],
        }
        self.fasta = {r["internal_id"]: r["sequence"] for r in self.metadata}
        self.patches = [
            patch.object(handoff, "EXPECTED_RECORDS", {"metadata": 3, "split": 3}),
            patch.object(handoff, "EXPECTED_COHORTS", {"primary": 2, "sensitivity_partial": 1}),
            patch.object(handoff, "EXPECTED_SPLITS", {"discovery": 2, "validation": 1}),
            patch.object(handoff, "EXPECTED_CLUSTERS", {"discovery": 1, "validation": 1}),
            patch.object(handoff, "EXPECTED_MANIFESTS", {name: 1 for name in self.manifests}),
        ]
        for item in self.patches:
            item.start()
            self.addCleanup(item.stop)

    def test_valid_partition(self):
        result = self._check()
        self.assertEqual(result["cross_split_clusters"], 0)

    def test_cross_split_homology_leakage_rejected(self):
        self.split[1]["homology_cluster"] = "g1"
        with self.assertRaisesRegex(ValueError, "Homology cluster leakage"):
            self._check()

    def test_nonmatching_sequence_hash_rejected(self):
        self.metadata[0]["sequence_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "Sequence hash mismatch"):
            self._check()

    def test_manifest_split_mismatch_rejected(self):
        self.manifests["discovery_primary_manifest.csv"][0]["split"] = "validation"
        with self.assertRaisesRegex(ValueError, "Manifest/split mismatch"):
            self._check()

    def _check(self):
        return handoff.check_rows(self.metadata, self.split, self.manifests, self.fasta)


if __name__ == "__main__":
    unittest.main()
