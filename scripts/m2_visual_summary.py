#!/usr/bin/env python3
"""Generate M2 visual summary and a scripted spot check for role E.

Figures and the spot-check table are written from the local metadata, split,
PF00741 scan, and 7R1C map. HMM occupancy is family-coordinate coverage, not
evidence of biochemical function.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATASET_VERSION = "gvpa-recognition-c7f6f005d717"
SPLIT_VERSION = "homology-8b9005e2d9-s42"
EXPERIMENT_ID = "M2-E-FIGURES-001"
LIMITATION = "HMM occupancy and PF00741 coverage are homology coordinates, not functional proof."
FONT_SIZES = {"label": 9.5, "tick": 8, "legend": 8.5, "title": 9.5}
FIGSIZE = (3.4, 2.5)


def configure_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.labelweight": "normal",
            "legend.frameon": False,
            "figure.dpi": 300,
            "axes.unicode_minus": False,
        }
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style_axes(ax) -> None:
    ax.tick_params(axis="both", labelsize=FONT_SIZES["tick"], pad=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.3)
    ax.spines["bottom"].set_color("#333333")
    ax.spines["left"].set_linewidth(1.3)
    ax.spines["left"].set_color("#333333")


def fit_text(fig, artist, max_width: float, start: float, minimum: float = 6.0) -> None:
    """Shrink a text artist until it fits max_width in display pixels."""
    renderer = fig.canvas.get_renderer()
    size = start
    artist.set_fontsize(size)
    fig.canvas.draw()
    while artist.get_window_extent(renderer).width > max_width and size > minimum:
        size -= 0.5
        artist.set_fontsize(size)
        fig.canvas.draw()


def fit_title(fig, ax, title: str) -> None:
    artist = ax.set_title(title, fontsize=FONT_SIZES["title"])
    fig.canvas.draw()
    limit = ax.get_window_extent().width * 0.98
    fit_text(fig, artist, limit, FONT_SIZES["title"])


def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def save_bar(path: Path, labels: list[str], values: list[int], title: str, ylabel: str) -> None:
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.bar(range(len(labels)), values, color="#1f4e79")
    ax.set_xticks(range(len(labels)))
    tick_labels = ax.set_xticklabels(labels, fontsize=FONT_SIZES["tick"])
    ax.set_ylabel(ylabel, fontsize=FONT_SIZES["label"], labelpad=3)
    style_axes(ax)
    fig.tight_layout()
    fig.canvas.draw()
    slot = ax.get_window_extent().width / max(len(labels), 1) * 0.92
    widest = max(artist.get_window_extent(fig.canvas.get_renderer()).width for artist in tick_labels)
    if widest > slot:
        ax.cla()
        ax.barh(range(len(labels)), values, color="#1f4e79")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=FONT_SIZES["tick"])
        ax.invert_yaxis()
        ax.set_xlabel(ylabel, fontsize=FONT_SIZES["label"], labelpad=3)
        ceiling = max(values) * 1.18
        ax.set_xlim(0, ceiling)
        for index, value in enumerate(values):
            ax.text(value + ceiling * 0.02, index, str(value), va="center", ha="left", fontsize=FONT_SIZES["tick"])
        style_axes(ax)
        fig.tight_layout()
    fit_title(fig, ax, title)
    save_figure(fig, path)


def cohort_counts(metadata: list[dict[str, str]]) -> Counter[str]:
    return Counter(row["analysis_cohort"] for row in metadata)


def cluster_sizes(cluster_tsv: Path) -> list[int]:
    clusters: dict[str, int] = defaultdict(int)
    for line in cluster_tsv.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        representative, _member = line.split("\t")
        clusters[representative] += 1
    return list(clusters.values())


def occupancy(coordinates: list[dict[str, str]]) -> list[int]:
    occupied: dict[int, set[str]] = defaultdict(set)
    for row in coordinates:
        if row["is_insertion"] == "True" or not row["hmm_match_state"]:
            continue
        occupied[int(row["hmm_match_state"])].add(row["internal_id"])
    return [len(occupied.get(state, ())) for state in range(1, 40)]


def first_where(rows: list[dict[str, str]], predicate) -> dict[str, str] | None:
    return next((row for row in rows if predicate(row)), None)


def spot_rows(metadata: list[dict[str, str]], coordinates: list[dict[str, str]], structure: list[dict[str, str]]) -> list[dict[str, str]]:
    by_id = {row["internal_id"]: row for row in metadata}
    checks: list[tuple[str, dict[str, str] | None, str]] = []
    checks.append(("primary", first_where(metadata, lambda row: row["analysis_cohort"] == "primary"), "Primary cohort row is QC-eligible and marked primary."))
    checks.append(("partial", first_where(metadata, lambda row: row["is_partial_any_member"] == "true" and row["has_type_conflict"] == "false"), "Partial flag comes from the annotation text, not from a filename."))
    checks.append(("type_conflict", first_where(metadata, lambda row: row["has_type_conflict"] == "true" and row["is_partial_any_member"] == "false"), "A non-GvpA member type is flagged and kept out of the primary cohort."))
    checks.append(("partial_and_type_conflict", first_where(metadata, lambda row: row["analysis_cohort"] == "sensitivity_partial_and_type_conflict"), "Combined partial and type-conflict row stays in a sensitivity cohort."))
    excluded = first_where(metadata, lambda row: row["sequence_qc_eligible"] == "false")
    checks.append(("nonstandard_or_excluded", excluded, "Sequence-QC exclusion is explicit and the row is absent from clustering."))
    outlier = first_where(metadata, lambda row: row["is_length_outlier"] == "true")
    checks.append(("length_outlier", outlier, "No length outlier was flagged under the mean +/- 3 SD rule." if outlier is None else "Length outlier is flagged rather than silently dropped."))

    modeled = first_where(structure, lambda row: row["modeled"] == "True" and row["hmm_match_state"])
    unmodeled = first_where(structure, lambda row: row["modeled"] == "False")
    if modeled is not None:
        checks.append(("structure_modeled", {"internal_id": "7R1C_N", "sequence_id": "7R1C"}, f"Deposited position {modeled['deposited_position']} residue {modeled['residue']} is modeled as PDB residue {modeled['pdb_residue_number']} and HMM state {modeled['hmm_match_state']}."))
    if unmodeled is not None:
        checks.append(("structure_unmodeled", {"internal_id": "7R1C_N", "sequence_id": "7R1C"}, f"Deposited position {unmodeled['deposited_position']} residue {unmodeled['residue']} has no experimental CA atom."))

    output = []
    for name, row, conclusion in checks:
        output.append({
            "check": name,
            "internal_id": "" if row is None else row["internal_id"],
            "sequence_id": "" if row is None else row.get("sequence_id", ""),
            "conclusion": conclusion,
        })
    return output


def hmm_spot_checks(metadata: list[dict[str, str]], coordinates: list[dict[str, str]], status: list[dict[str, str]]) -> list[dict[str, str]]:
    by_id = {row["internal_id"]: row for row in metadata}
    rows = []
    sample = next(row for row in coordinates if row["is_insertion"] == "False")
    sequence = by_id[sample["internal_id"]]["sequence"]
    position = int(sample["raw_position"])
    matches = sequence[position - 1] == sample["residue"]
    rows.append({
        "check": "hmm_residue_identity",
        "internal_id": sample["internal_id"],
        "sequence_id": by_id[sample["internal_id"]]["sequence_id"],
        "conclusion": f"Raw position {position} is {sample['residue']} in both the HMM map and the metadata sequence." if matches else "Raw position residue disagrees with the metadata sequence.",
    })
    deleted = next((row for row in status if int(row["n_deleted_match_states"]) > 0), None)
    if deleted is None:
        rows.append({"check": "hmm_deletion", "internal_id": "", "sequence_id": "", "conclusion": "No scanned sequence has a deleted PF00741 match state."})
    else:
        meta = by_id[deleted["internal_id"]]
        rows.append({"check": "hmm_deletion", "internal_id": deleted["internal_id"], "sequence_id": meta["sequence_id"], "conclusion": f"{deleted['n_deleted_match_states']} HMM match states are deleted and were not assigned a raw residue."})
    insertion = next((row for row in coordinates if row["is_insertion"] == "True"), None)
    if insertion is None:
        rows.append({"check": "hmm_insertion", "internal_id": "", "sequence_id": "", "conclusion": "The GA alignment contains no insertion, so no raw residue lacks a match state for that reason."})
    else:
        meta = by_id[insertion["internal_id"]]
        rows.append({"check": "hmm_insertion", "internal_id": insertion["internal_id"], "sequence_id": meta["sequence_id"], "conclusion": f"Raw position {insertion['raw_position']} is an insertion and has no HMM match state."})
    return rows


def render_spot_check(rows: list[dict[str, str]], questions: list[tuple[str, str]]) -> str:
    lines = [
        "# M2 manual spot check",
        "",
        f"Generated by `scripts/m2_visual_summary.py`. Dataset `{DATASET_VERSION}`, split `{SPLIT_VERSION}`.",
        "",
        "HMM coverage and occupancy describe PF00741 match states. They are not functional labels.",
        "",
        "| Check | Internal ID | Sequence ID | Conclusion |",
        "| --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(f"| {row['check']} | {row['internal_id']} | {row['sequence_id']} | {row['conclusion']} |")
    lines.extend(["", "## Questions for the M2 review", ""])
    for question, evidence in questions:
        lines.append(f"- {question} Evidence: {evidence}")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=ROOT / "data/processed/gvpa_v1/metadata.csv")
    parser.add_argument("--clusters", type=Path, default=ROOT / "work/gvpa_cluster.tsv")
    parser.add_argument("--coordinates", type=Path, default=ROOT / "data/processed/m2_c/pf00741_scan_001/hmm_coordinate_map.csv")
    parser.add_argument("--status", type=Path, default=ROOT / "data/processed/m2_c/pf00741_scan_001/per_sequence_status.csv")
    parser.add_argument("--structure", type=Path, default=ROOT / "data/processed/m2_c/7r1c_e/7r1c_residue_hmm_map.csv")
    parser.add_argument("--scan-summary", type=Path, default=ROOT / "results/bioinformatics/pf00741_scan_summary.json")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "reports/M2/figures")
    args = parser.parse_args()

    configure_style()
    metadata = read_csv(args.metadata)
    coordinates = read_csv(args.coordinates)
    status = read_csv(args.status)
    structure = read_csv(args.structure)
    scan = json.loads(args.scan_summary.read_text(encoding="utf-8"))
    counts = cohort_counts(metadata)
    if sum(counts.values()) != 2078 or counts["primary"] != 1721:
        raise SystemExit(f"Cohort counts {dict(counts)} disagree with the frozen audit")
    if scan["counts"]["accepted"] != 2076:
        raise SystemExit("PF00741 accepted count disagrees with the scan summary")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    labels = ["primary", "partial", "type conflict", "partial+type", "excluded"]
    save_bar(args.out_dir / "cohort_counts.png", labels, [counts[label] for label in labels], "GvpA analysis cohorts (n=2078)", "Sequences")

    lengths = [int(row["sequence_length"]) for row in metadata]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.hist(lengths, bins=range(50, 185, 5), color="#1f4e79")
    ax.set_xlabel("Length (aa)", fontsize=FONT_SIZES["label"], labelpad=3)
    ax.set_ylabel("Sequences", fontsize=FONT_SIZES["label"], labelpad=3)
    style_axes(ax)
    fig.tight_layout()
    fit_title(fig, ax, "GvpA sequence length")
    save_figure(fig, args.out_dir / "length_distribution.png")

    sizes = cluster_sizes(args.clusters)
    if len(sizes) != 478 or sum(sizes) != 2076:
        raise SystemExit(f"Cluster table has {len(sizes)} clusters and {sum(sizes)} sequences")
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.hist(sizes, bins=30, color="#1f4e79")
    ax.set_xlabel("Sequences in cluster", fontsize=FONT_SIZES["label"], labelpad=3)
    ax.set_ylabel("Clusters", fontsize=FONT_SIZES["label"], labelpad=3)
    style_axes(ax)
    fig.tight_layout()
    fit_title(fig, ax, "Homology cluster sizes at identity 0.8")
    save_figure(fig, args.out_dir / "cluster_sizes.png")

    split_counts = Counter(row["split"] for row in status)
    save_bar(args.out_dir / "discovery_validation.png", ["discovery", "validation"], [split_counts["discovery"], split_counts["validation"]], "Frozen discovery/validation split", "QC-eligible sequences")
    status_counts = Counter(row["status"] for row in status)
    save_bar(args.out_dir / "pf00741_status.png", list(status_counts), list(status_counts.values()), "PF00741 GA status for 2076 sequences", "Sequences")

    covered = occupancy(coordinates)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.bar(range(1, 40), covered, color="#1f4e79")
    ax.set_xlabel("HMM match state", fontsize=FONT_SIZES["label"], labelpad=3)
    ax.set_ylabel("Sequences", fontsize=FONT_SIZES["label"], labelpad=3)
    ax.set_ylim(0, 2076)
    style_axes(ax)
    fig.tight_layout()
    fit_title(fig, ax, "Sequences occupying each PF00741 match state")
    save_figure(fig, args.out_dir / "hmm_occupancy.png")

    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis("off")
    steps = ["2078 candidates", "2076 QC-eligible", "478 clusters", "1453/623 split"]
    texts = []
    for index, step in enumerate(steps):
        y = 6.6 - index * 1.7
        ax.add_patch(plt.Rectangle((1.2, y), 7.6, 1.2, fill=False, linewidth=1.2))
        texts.append(ax.text(5.0, y + 0.6, step, ha="center", va="center", fontsize=FONT_SIZES["label"]))
    fig.tight_layout()
    fig.canvas.draw()
    for artist in texts:
        fit_text(fig, artist, ax.get_window_extent().width * 0.7, FONT_SIZES["label"])
    fit_title(fig, ax, "Data filter counted from the frozen audit")
    save_figure(fig, args.out_dir / "filter_flow.png")

    figures = [
        "filter_flow.png", "cohort_counts.png", "length_distribution.png", "cluster_sizes.png",
        "discovery_validation.png", "pf00741_status.png", "hmm_occupancy.png",
    ]
    write_csv(ROOT / "reports/M2/figure_manifest.csv", [
        {
            "figure": name,
            "path": f"reports/M2/figures/{name}",
            "experiment_id": EXPERIMENT_ID,
            "dataset_version": DATASET_VERSION,
            "split_version": SPLIT_VERSION,
            "script": "scripts/m2_visual_summary.py",
            "limitation": LIMITATION,
        }
        for name in figures
    ])
    rows = spot_rows(metadata, coordinates, structure) + hmm_spot_checks(metadata, coordinates, status)
    questions = [
        ("Why is there no supervised label?", "reports/M2/m2_reproducibility.md and the audit label_status field."),
        ("Did homology clusters cross discovery and validation?", "Frozen split cluster_leakage is false; D's validator checks the manifest."),
        ("Does PF00741 acceptance prove function?", "No. The scan summary states that a hit is not a functional label."),
        ("Why do 0.7 and 0.9 change the split?", "results/data_audit/mmseqs_sensitivity.json; the frozen 0.8 split was not replaced."),
    ]
    (ROOT / "reports/M2/manual_spot_check.md").write_text(render_spot_check(rows, questions), encoding="utf-8")
    print(f"wrote {len(figures)} figures and {len(rows)} spot checks")


if __name__ == "__main__":
    main()
