import csv
import math

import pytest

from bioinformatics.regions import extract_regions, read_scores, run


def rows(values):
    return [{"position": i, "conservation": v, "sufficient_support": v is not None}
            for i, v in enumerate(values, 1)]


def test_runs_do_not_bridge_missing_or_low_support():
    data = rows([1, 0.9, None, 1, 1, 1, 1, 0.8, 1])
    data[4]["sufficient_support"] = False
    result = extract_regions(data, 0.9, 2)
    assert [(r["start"], r["end"]) for r in result] == [(1, 2), (6, 7)]
    assert result[0]["mean_conservation"] == pytest.approx(0.95)


def test_terminal_run_and_minimum_length():
    assert extract_regions(rows([0.8, 1, 1]), 0.9, 3) == []
    result = extract_regions(rows([0.8, 1, 1]), 0.9, 2)
    assert (result[0]["start"], result[0]["end"], result[0]["length"]) == (2, 3, 2)


@pytest.mark.parametrize("threshold,length", [(math.nan, 1), (math.inf, 1), (-0.1, 2), (1.1, 1), (0.9, 0)])
def test_invalid_parameters(threshold, length):
    with pytest.raises(ValueError):
        extract_regions(rows([1]), threshold, length)


@pytest.mark.parametrize("position,score,support", [(2, "1", "True"), (1, "nan", "True"),
    (1, "", "True"), (1, "1", "maybe"), (1, "2", "True")])
def test_invalid_scores(tmp_path, position, score, support):
    path = tmp_path / "scores.csv"
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["msa_position", "conservation", "sufficient_support"])
        writer.writerow([position, score, support])
    with pytest.raises(ValueError):
        read_scores(path)


@pytest.mark.parametrize("content", [
    "msa_position,conservation,sufficient_support,conservation\n1,0.0,True,1.0\n",
    "msa_position,conservation,sufficient_support\n1,0.0,True,1.0\n",
    "msa_position,conservation,sufficient_support\n1,0.0\n",
])
def test_rejects_ambiguous_csv(tmp_path, content):
    path = tmp_path / "scores.csv"
    path.write_text(content)
    with pytest.raises(ValueError):
        read_scores(path)


def test_direct_api_rejects_missing_column():
    values = rows([1, 1])
    values[1]["position"] = 3
    with pytest.raises(ValueError, match="consecutive"):
        extract_regions(values, 0.9, 1)


def test_write_failure_has_marker(tmp_path, monkeypatch):
    from types import SimpleNamespace
    path = tmp_path / "scores.csv"
    path.write_text("msa_position,conservation,sufficient_support\n1,1,True\n")
    def fail(*args, **kwargs):
        raise OSError("simulated disk error")
    monkeypatch.setattr("bioinformatics.regions.write_csv", fail)
    out = tmp_path / "failed"
    with pytest.raises(OSError, match="disk"):
        run(SimpleNamespace(scores=path, min_conservation=0.9, min_length=1, out_dir=out))
    assert (out / "FAILED.txt").exists()
    assert not (out / "manifest.json").exists()
