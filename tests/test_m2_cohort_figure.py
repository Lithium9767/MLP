import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

from scripts import m2_visual_summary as visual
from scripts import rerun_m2_cohort as rerun


class CohortFigureTest(unittest.TestCase):
    def setUp(self):
        # Real schema, deliberately unequal nonzero counts: aliases would all be zero.
        self.counts = {"primary": 7, "sensitivity_partial": 2,
                       "sensitivity_type_conflict": 3,
                       "sensitivity_partial_and_type_conflict": 4,
                       "excluded_sequence_qc": 1}
        self.metadata = [{"analysis_cohort": key} for key, n in self.counts.items() for _ in range(n)]

    def test_actual_matplotlib_bars_match_source_categories(self):
        captured = []
        def inspect_figure(fig, path):
            ax = fig.axes[0]
            horizontal = ax.get_xlabel() == "Sequences"
            captured.extend(p.get_width() if horizontal else p.get_height() for p in ax.patches)
            visual.plt.close(fig)
        visual.configure_style()
        with patch.object(visual, "save_figure", side_effect=inspect_figure):
            rows = visual.render_cohort_figure(self.metadata, Path("unused.png"))
        source = Counter(row["analysis_cohort"] for row in self.metadata)
        self.assertEqual(captured, [source[row["analysis_cohort"]] for row in rows])
        self.assertEqual(captured, [7, 2, 3, 4, 1])

    def test_nonzero_source_with_zero_rendered_bar_is_rejected(self):
        with patch.object(visual, "save_bar", return_value=[7, 0, 0, 0, 0]):
            with self.assertRaisesRegex(ValueError, "Rendered/source cohort mismatch"):
                visual.render_cohort_figure(self.metadata, Path("unused.png"))

    def test_unknown_category_is_not_silently_hidden(self):
        with self.assertRaisesRegex(ValueError, "unknown analysis_cohort"):
            visual.render_cohort_figure(self.metadata + [{"analysis_cohort": "typo"}], Path("unused.png"))

    def test_rerun_rejects_dirty_tree_before_touching_inputs(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(rerun, "ROOT", Path(directory)):
            with patch.object(rerun.subprocess, "check_output", return_value=" M scripts/code.py\n"):
                with self.assertRaisesRegex(ValueError, "working tree must be clean"):
                    rerun.run(Path("missing.csv"), Path("missing.json"), "M2-E-FIGURES-002")
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_existing_run_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(rerun, "ROOT", Path(directory)):
            run_dir = Path(directory) / "reports/M2/runs/M2-E-FIGURES-002"
            run_dir.mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "Run ID already exists"):
                rerun.run(Path("missing.csv"), Path("missing.json"), "M2-E-FIGURES-002")


if __name__ == "__main__":
    unittest.main()
