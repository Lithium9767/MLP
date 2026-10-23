"""Frozen sequence, coordinate and window contracts for spectral label transfer.

The 7R1C extension rule selects sufficiently mapped windows; it never supplies
region or function labels to the embedding features. Coordinates are one-based
and inclusive. Alignment thresholds reproduce the discovery audit unchanged.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np

from bioinformatics.m2_pf00741 import CANONICAL, csv_rows, sequence_digest


EMBEDDING_DIMENSION = 480
ALIGNMENT_PARAMETERS = ((-8, -0.5), (-10, -0.5), (-12, -1), (-10, -2))
MAX_ALIGNMENT_TIES = 32
VALID_SPLITS = frozenset({"discovery", "validation"})


def _sequence(row: Mapping) -> str:
    sequence = row["sequence"]
    if (not sequence or set(sequence) - CANONICAL
            or len(sequence) != int(row["sequence_length"])
            or sequence_digest(sequence) != row["sequence_sha256"]):
        raise ValueError(f"Sequence, length or SHA-256 mismatch: {row['internal_id']}")
    return sequence


def load_primary_records(metadata_path: Path, split_path: Path, split: str) -> list[dict]:
    """Join the frozen split to QC-eligible primary records, in internal-ID order.

    Duplicate IDs, inconsistent eligibility, unknown assignments and homology
    groups crossing the split are errors rather than silently dropped rows.
    """
    if split not in VALID_SPLITS:
        raise ValueError(f"Unknown split: {split}")
    metadata = csv_rows(Path(metadata_path), {
        "internal_id", "sequence_id", "sequence", "sequence_length", "sequence_sha256",
        "sequence_qc_eligible", "primary_analysis_eligible", "analysis_cohort"})
    by_id, accessions, sequence_hashes = {}, set(), set()
    for row in metadata:
        sid = row["internal_id"]
        if (not sid or sid in by_id or not row["sequence_id"]
                or row["sequence_id"] in accessions or row["sequence_sha256"] in sequence_hashes):
            raise ValueError("Missing or duplicate sequence ID in metadata")
        if (row["sequence_qc_eligible"] not in {"true", "false"}
                or row["primary_analysis_eligible"] not in {"true", "false"}):
            raise ValueError(f"Invalid eligibility flags: {sid}")
        if row["primary_analysis_eligible"] == "true" and (
                row["sequence_qc_eligible"] != "true" or row["analysis_cohort"] != "primary"):
            raise ValueError(f"Inconsistent primary eligibility: {sid}")
        by_id[sid] = row
        accessions.add(row["sequence_id"])
        sequence_hashes.add(row["sequence_sha256"])
    assignments, group_split = {}, {}
    for row in csv_rows(Path(split_path), {
            "internal_id", "sequence_id", "split", "homology_cluster",
            "analysis_cohort", "primary_analysis_eligible"}):
        sid, subset, group = row["internal_id"], row["split"], row["homology_cluster"]
        if sid not in by_id or sid in assignments or subset not in VALID_SPLITS or not group:
            raise ValueError("Unknown/duplicate split ID, invalid split or missing homology group")
        source = by_id[sid]
        if source["sequence_qc_eligible"] != "true" or any(
                source[key] != row[key] for key in (
                    "sequence_id", "analysis_cohort", "primary_analysis_eligible")):
            raise ValueError(f"Split and QC/metadata disagree: {sid}")
        if group in group_split and group_split[group] != subset:
            raise ValueError(f"Homology group crosses discovery/validation: {group}")
        assignments[sid], group_split[group] = row, subset
    eligible = {sid for sid, row in by_id.items() if row["sequence_qc_eligible"] == "true"}
    if set(assignments) != eligible:
        raise ValueError("Split manifest does not cover every QC-eligible sequence")
    records = []
    for sid in sorted(assignments):
        row, assignment = by_id[sid], assignments[sid]
        if assignment["split"] == split and row["primary_analysis_eligible"] == "true":
            _sequence(row)
            records.append({**row, "split": split,
                            "homology_cluster": assignment["homology_cluster"]})
    if not records:
        raise ValueError(f"No eligible primary records in {split}")
    return records


def _record_index(records: list[dict]) -> dict[str, dict]:
    result = {}
    for row in records:
        sid = row["internal_id"]
        if (not sid or sid in result or row["split"] not in VALID_SPLITS
                or not row["homology_cluster"]):
            raise ValueError("Missing/duplicate ID or invalid split/homology group")
        _sequence(row)
        result[sid] = row
    if not result:
        raise ValueError("No sequence records")
    return result


def _validate_anchors(sequence: str, anchors: Mapping[int, int]) -> None:
    if (len(anchors) < 20 or any(type(h) is not int or not 1 <= h <= 39 for h in anchors)
            or any(type(q) is not int or not 1 <= q <= len(sequence) for q in anchors.values())):
        raise ValueError("Require at least 20 actual HMM states in 1..39 at valid raw positions")
    positions = [anchors[h] for h in sorted(anchors)]
    if any(a >= b for a, b in zip(positions, positions[1:])):
        raise ValueError("HMM anchors must have unique, monotonically increasing raw positions")


def load_hmm_anchors(coordinates: Path, records: list[dict]) -> dict[str, dict[int, int]]:
    """Read actual HMM matches, checking row identity; never impute deleted states."""
    by_id = _record_index(records)
    anchors = {sid: {} for sid in by_id}
    seen = {sid: set() for sid in by_id}
    for row in csv_rows(Path(coordinates), {
            "internal_id", "sequence_sha256", "split", "analysis_cohort",
            "raw_position", "residue", "hmm_match_state"}):
        sid = row["internal_id"]
        if sid not in by_id:
            continue
        source = by_id[sid]
        position = int(row["raw_position"])
        if (row["sequence_sha256"] != source["sequence_sha256"]
                or row["split"] != source["split"] or row["analysis_cohort"] != "primary"
                or not 1 <= position <= len(source["sequence"])
                or source["sequence"][position - 1] != row["residue"]):
            raise ValueError(f"Coordinate sequence SHA, residue, split or position mismatch: {sid}")
        if position in seen[sid]:
            raise ValueError(f"Duplicate coordinate raw position: {sid}:{position}")
        seen[sid].add(position)
        if row["hmm_match_state"]:
            state = int(row["hmm_match_state"])
            if state in anchors[sid]:
                raise ValueError(f"Duplicate HMM state: {sid}:{state}")
            anchors[sid][state] = position
    for sid, present in anchors.items():
        _validate_anchors(by_id[sid]["sequence"], present)
    return anchors


def _extension_maps(ref_positions, query_positions, reference, sequence, matrix):
    """Return only pairs shared by all optima under all four gap settings."""
    from Bio.Align import PairwiseAligner

    target = "".join(reference[p - 1] for p in ref_positions)
    query = "".join(sequence[p - 1] for p in query_positions)
    if not target or not query:
        return {}, {"max_optimal_ties": 0, "truncated_ties": False, "scores": []}
    maps, ties, scores = [], [], []
    for opening, extending in ALIGNMENT_PARAMETERS:
        aligner = PairwiseAligner()
        aligner.mode, aligner.substitution_matrix = "global", matrix
        aligner.open_gap_score, aligner.extend_gap_score = opening, extending
        # Both strands pay gaps adjacent to the anchor; only distal gaps are free.
        aligner.open_right_insertion_score = aligner.extend_right_insertion_score = 0
        aligner.open_right_deletion_score = aligner.extend_right_deletion_score = 0
        alignments = aligner.align(target, query)
        scores.append(float(alignments.score))
        try:
            count = len(alignments)
        except OverflowError:
            count = MAX_ALIGNMENT_TIES + 1
        ties.append(count)
        for index, alignment in enumerate(alignments):
            if index >= MAX_ALIGNMENT_TIES:
                break
            pairs = {}
            for ref_block, query_block in zip(*alignment.aligned):
                for ref_index, query_index in zip(range(*ref_block), range(*query_block)):
                    pairs[ref_positions[ref_index]] = query_positions[query_index]
            maps.append(pairs)
    truncated = max(ties) > MAX_ALIGNMENT_TIES
    stable = {} if truncated else {
        p: q for p, q in maps[0].items() if all(mapping.get(p) == q for mapping in maps)}
    return stable, {"max_optimal_ties": max(ties), "truncated_ties": truncated, "scores": scores}


def _quality(positions, mapping, reference, sequence, matrix, minimum):
    mapped = [p for p in positions if p in mapping]
    count = len(mapped)
    identity = sum(reference[p - 1] == sequence[mapping[p] - 1] for p in mapped) / count if count else 0.
    positive = sum(matrix[reference[p - 1], sequence[mapping[p] - 1]] > 0 for p in mapped) / count if count else 0.
    return {"mapped_residues": count, "region_residues": len(positions), "identity": identity,
            "positive_substitution_fraction": positive,
            "passes_descriptive_sequence_support": bool(count >= minimum and identity >= .25 and positive >= .6)}


def accepted_reference_mapping(sequence: str, anchors: dict[int, int], reference: str) -> tuple[dict, dict]:
    """Transfer the frozen 88-residue 7R1C correspondence using discovery rules.

    All actual HMM states retain reference position ``state + 10``. Only modeled
    extensions 2..10 and 50..66 can be added. Missing core states stay missing.
    Quality is descriptive alignment support, not a correctness probability.
    """
    from Bio.Align import substitution_matrices

    if not sequence or set(sequence) - CANONICAL or len(reference) != 88 or set(reference) - CANONICAL:
        raise ValueError("Require canonical sequence and frozen 88-residue reference")
    _validate_anchors(sequence, anchors)
    matrix = substitution_matrices.load("BLOSUM62")
    core = {h + 10: q for h, q in anchors.items()}
    first, last = min(core), max(core)
    nmap, nd = _extension_maps(list(range(first - 1, 0, -1)),
                              list(range(core[first] - 1, 0, -1)), reference, sequence, matrix)
    cmap, cd = _extension_maps(list(range(last + 1, 89)),
                              list(range(core[last] + 1, len(sequence) + 1)), reference, sequence, matrix)
    nmap, cmap = {p: q for p, q in nmap.items() if p < 11}, {p: q for p, q in cmap.items() if p > 49}
    nq = _quality(range(2, 11), nmap, reference, sequence, matrix, 7)
    aq = _quality(range(50, 63), cmap, reference, sequence, matrix, 10)
    cq = _quality(range(50, 67), cmap, reference, sequence, matrix, 13)
    accepted = dict(core)
    if nq["passes_descriptive_sequence_support"]:
        accepted.update({p: q for p, q in nmap.items() if 2 <= p <= 66})
    if aq["passes_descriptive_sequence_support"]:
        accepted.update({p: q for p, q in cmap.items() if 50 <= p <= 62})
        if cq["passes_descriptive_sequence_support"]:
            accepted.update({p: q for p, q in cmap.items() if 63 <= p <= 66})
    reverse = {q: p for p, q in accepted.items()}
    if len(reverse) != len(accepted):
        raise ValueError("Reference transfer contains duplicate raw positions")
    return reverse, {"N": nq, "alpha2": aq, "C_context": cq,
                     "N_alignment": nd, "C_alignment": cd}


def build_window_features(records: list[dict], embeddings: dict[str, np.ndarray],
                          anchors: dict, mappings: dict, *, width: int = 10,
                          step: int = 5, min_accepted: int = 8) -> tuple[list[dict], np.ndarray, list[dict]]:
    """Build full-protein ESM windows minus the mean of actual HMM-match residues.

    Embeddings must be residue arrays ``(sequence_length, 480)`` keyed by internal
    ID, with provenance checked by the caller's embedding loader. Returned rows
    include all windows, but only accepted windows enter the feature matrix.
    """
    if (any(type(value) is not int for value in (width, step, min_accepted))
            or width < 1 or step < 1 or not 1 <= min_accepted <= width):
        raise ValueError("Invalid width, step or minimum accepted residue count")
    by_id = _record_index(records)
    required = set(by_id)
    if any(not required <= set(source) for source in (embeddings, anchors, mappings)):
        raise ValueError("Missing sequence embedding, HMM anchors or accepted mapping")
    selected_rows, features, all_rows = [], [], []
    for sid in sorted(by_id):
        row, embedding = by_id[sid], np.asarray(embeddings[sid])
        sequence = row["sequence"]
        if (embedding.shape != (len(sequence), EMBEDDING_DIMENSION)
                or embedding.dtype.kind != "f" or not np.isfinite(embedding).all()):
            raise ValueError(f"Invalid residue embedding shape, dtype or finite values: {sid}")
        _validate_anchors(sequence, anchors[sid])
        mapping = mappings[sid]
        if (not isinstance(mapping, Mapping)
                or any(type(q) is not int or not 1 <= q <= len(sequence) for q in mapping)
                or any(type(p) is not int or not 2 <= p <= 66 for p in mapping.values())):
            raise ValueError(f"Invalid accepted reference mapping: {sid}")
        ordered_reference = [mapping[q] for q in sorted(mapping)]
        if (len(set(mapping.values())) != len(mapping)
                or any(a >= b for a, b in zip(ordered_reference, ordered_reference[1:]))
                or any(mapping.get(q) != h + 10 for h, q in anchors[sid].items())
                or any(11 <= p <= 49 and anchors[sid].get(p - 10) != q for q, p in mapping.items())):
            raise ValueError(f"Accepted mapping is invalid or changes a frozen HMM anchor: {sid}")
        background = embedding[np.array([anchors[sid][h] for h in sorted(anchors[sid])]) - 1].mean(axis=0)
        for start in range(1, len(sequence) - width + 2, step):
            end = start + width - 1
            reference_positions = [mapping[q] for q in range(start, end + 1) if q in mapping]
            selected = len(reference_positions) >= min_accepted
            result = {key: row[key] for key in (
                "internal_id", "sequence_sha256", "homology_cluster", "split", "sequence_length")}
            result.update(window_length=width, raw_start=start, raw_end=end,
                          accepted_mapping_count=len(reference_positions),
                          reference_positions=";".join(map(str, reference_positions)),
                          fragment=sequence[start - 1:end], selected=selected)
            all_rows.append(result)
            if selected:
                selected_rows.append(dict(result))
                features.append(embedding[start - 1:end].mean(axis=0) - background)
    matrix = np.asarray(features, dtype=np.float32).reshape(-1, EMBEDDING_DIMENSION)
    if not np.isfinite(matrix).all():
        raise ValueError("Non-finite window features")
    return selected_rows, matrix, all_rows
