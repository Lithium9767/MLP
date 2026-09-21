import csv
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

from preprocessing.homology_split import assign_clusters, parse_ratios, read_clusters, read_eligible_ids
from preprocessing.prepare_gvpa_dataset import DEFAULT_MEMBER, build_rows, load_source, write_outputs


def annotation(identifier, *, types=("GvpA",), description="protein", source="source.fasta"):
    return {
        "id": identifier,
        "description": description,
        "full_header": description,
        "organism": "Test organism",
        "source_file": source,
        "gene": "gvpa",
        "gvp_types": list(types),
    }


def record(identifier, sequence, representative=None, members=None):
    representative = representative or annotation(identifier)
    members = [representative] if members is None else members
    return {
        "unique_sequence_id": identifier,
        "sequence": sequence,
        "representative_annotation": representative,
        "cluster_statistics": {
            "total_original_sequences": len(members) + 1,
            "redundant_count": len(members),
        },
        "redundant_members_details": members,
    }


class PrepareDatasetTest(unittest.TestCase):
    def test_quality_flags_and_outputs(self):
        partial = annotation("a", description="gas vesicle protein, partial")
        conflict = annotation("b_member", types=("GvpJ",))
        records = [
            record("a", "ACDEFG", partial, [partial]),
            record("b", "ACDEFA", members=[annotation("b"), conflict]),
            record("c", "ACDEXG"),
            record("d", "ACDEFG"),
        ]
        rows, members, summary = build_rows(records)
        self.assertEqual(summary["candidate_records"], 4)
        self.assertEqual(summary["total_original_sequence_records"], 9)
        self.assertEqual(summary["redundant_member_records"], 5)
        self.assertEqual(summary["partial_any_member"], 1)
        self.assertEqual(summary["type_conflict"], 1)
        self.assertEqual(summary["nonstandard_sequence"], 1)
        self.assertEqual(summary["exact_duplicate_sequence"], 1)
        self.assertEqual(rows[0]["analysis_cohort"], "sensitivity_partial")
        self.assertEqual(rows[1]["analysis_cohort"], "sensitivity_type_conflict")
        self.assertEqual(rows[2]["sequence_qc_eligible"], "false")
        self.assertEqual(rows[3]["exact_duplicate_of"], "a")

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            completed = write_outputs(
                rows,
                members,
                summary,
                {
                    "source_zip": "test.zip",
                    "source_zip_sha256": "0" * 64,
                    "source_member": DEFAULT_MEMBER,
                    "source_member_size": 1,
                    "source_member_sha256": "1" * 64,
                },
                output,
            )
            fasta = (output / "sequences_for_clustering.fasta").read_text()
            self.assertIn(">a", fasta)
            self.assertIn(">b", fasta)
            self.assertNotIn(">c", fasta)
            self.assertNotIn(">d", fasta)
            self.assertEqual(completed["dataset_version"], "gvpa-recognition-111111111111")

    def test_loads_expected_zip_member(self):
        with tempfile.TemporaryDirectory() as directory:
            archive_path = Path(directory) / "gv.zip"
            payload = [record("a", "ACDEFG")]
            with zipfile.ZipFile(archive_path, "w") as archive:
                archive.writestr(DEFAULT_MEMBER, json.dumps(payload))
            records, provenance = load_source(archive_path)
            self.assertEqual(records[0]["unique_sequence_id"], "a")
            self.assertEqual(provenance["source_member"], DEFAULT_MEMBER)
            self.assertEqual(len(provenance["source_member_sha256"]), 64)

    def test_missing_required_field_rejected(self):
        with self.assertRaisesRegex(ValueError, "missing required field sequence"):
            build_rows([{"unique_sequence_id": "a"}])


class HomologySplitTest(unittest.TestCase):
    def _metadata(self, root: Path, ids=("a", "b", "c", "d")) -> Path:
        path = root / "metadata.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sequence_id", "sequence_qc_eligible"])
            writer.writeheader()
            writer.writerows({"sequence_id": item, "sequence_qc_eligible": "true"} for item in ids)
        return path

    def test_whole_clusters_stay_together_and_are_deterministic(self):
        clusters = {"r1": ["a", "b"], "r2": ["c"], "r3": ["d"]}
        ratios = parse_ratios("discovery=0.7,validation=0.3")
        first = assign_clusters(clusters, ratios, seed=42)
        second = assign_clusters(clusters, ratios, seed=42)
        self.assertEqual(first, second)
        self.assertEqual(set(first), set(clusters))
        self.assertGreaterEqual(len(set(first.values())), 2)

    def test_cluster_coverage_and_duplicate_assignment_are_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = read_eligible_ids(self._metadata(root, ("a", "b")))
            missing = root / "missing.tsv"
            missing.write_text("r1\ta\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing 1 eligible"):
                read_clusters(missing, expected)
            duplicate = root / "duplicate.tsv"
            duplicate.write_text("r1\ta\nr2\ta\nr2\tb\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "appears more than once"):
                read_clusters(duplicate, expected)


if __name__ == "__main__":
    unittest.main()
