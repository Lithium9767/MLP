import tempfile
import unittest
from pathlib import Path

from preprocessing.mmseqs_cluster_split import assign_clusters, read_clusters, write_split
from preprocessing.prepare_gvpa_dataset import build_rows, write_outputs


class PrepareGvpADatasetTest(unittest.TestCase):
    def test_build_rows_and_exact_dedup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fasta"
            fasta.write_text(
                ">WP_1.1 protein [Species one]\nACDEFG\n"
                ">WP_2.1 duplicate [Species two]\nACDEFG\n"
                ">WP_3.1 protein [Species three]\nACDX\n"
            )
            rows, summary = build_rows(fasta)
            write_outputs(rows, summary, root / "output")

            self.assertEqual(summary["raw_records"], 3)
            self.assertEqual(summary["included_after_validation_and_exact_dedup"], 2)
            self.assertIn("exact_duplicate_of:WP_1.1", rows[1]["exclusion_reason"])
            self.assertEqual(rows[0]["organism"], "Species one")
            self.assertTrue((root / "output" / "metadata.csv").is_file())

    def test_group_split_keeps_cluster_together(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cluster_file = root / "clusters.tsv"
            cluster_file.write_text("r1\tr1\nr1\ta\nr2\tr2\nr3\tr3\n")
            clusters = read_clusters(cluster_file)
            assignment = assign_clusters(clusters, (0.5, 0.25, 0.25), seed=7)
            output = root / "split.csv"
            write_split(clusters, assignment, output)

            self.assertEqual(len(clusters), 3)
            self.assertIn(assignment["r1"], {"train", "val", "test"})
            lines = output.read_text().splitlines()
            a_rows = [line for line in lines if line.startswith(("r1,", "a,"))]
            self.assertEqual(len({line.rsplit(",", 1)[-1] for line in a_rows}), 1)


if __name__ == "__main__":
    unittest.main()
