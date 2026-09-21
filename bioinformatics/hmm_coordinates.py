#!/usr/bin/env python3
"""Map GvpA raw residue positions to PF00741 HMM match states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_text(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def parse_alignment_columns(
    hmm_sequence: str,
    target_sequence: str,
    hmm_from: int,
    target_from: int,
) -> list[dict[str, Any]]:
    """Return one row per aligned target residue, including insertion residues."""
    hmm_sequence = as_text(hmm_sequence)
    target_sequence = as_text(target_sequence)
    if len(hmm_sequence) != len(target_sequence):
        raise ValueError("HMM and target alignment strings have different lengths")
    if hmm_from < 1 or target_from < 1:
        raise ValueError("Alignment coordinates are 1-based positive integers")
    hmm_state = int(hmm_from)
    sequence_position = int(target_from)
    rows: list[dict[str, Any]] = []
    for hmm_char, target_char in zip(hmm_sequence, target_sequence):
        insertion = hmm_char == "."
        target_present = target_char not in ".-"
        if target_present:
            rows.append({
                "raw_position": sequence_position,
                "hmm_match_state": None if insertion else hmm_state,
                "is_insertion": insertion,
                "residue": target_char.upper(),
            })
            sequence_position += 1
        if not insertion:
            hmm_state += 1
    return rows


def map_window_to_hmm_states(
    coordinate_rows: list[dict[str, Any]], start: int, end: int
) -> list[int]:
    if start < 1 or end < start:
        raise ValueError("Window coordinates must be valid 1-based inclusive positions")
    return sorted({
        int(row["hmm_match_state"])
        for row in coordinate_rows
        if start <= int(row["raw_position"]) <= end and row["hmm_match_state"] is not None
    })


def shannon_conservation(residues: list[str]) -> tuple[float, float, str, float]:
    occupied = [residue for residue in residues if residue in CANONICAL_AA]
    occupancy = len(occupied) / len(residues) if residues else 0.0
    if not occupied:
        return 0.0, 0.0, "", occupancy
    counts = Counter(occupied)
    total = len(occupied)
    entropy = -sum((count / total) * math.log2(count / total) for count in counts.values())
    denominator = math.log2(min(20, total))
    conservation = 1.0 if denominator == 0 else max(0.0, 1.0 - entropy / denominator)
    consensus, consensus_count = counts.most_common(1)[0]
    return conservation, consensus_count / total, consensus, occupancy


def conservation_by_state(
    coordinate_rows: list[dict[str, Any]], total_sequences: int
) -> list[dict[str, Any]]:
    if total_sequences <= 0:
        raise ValueError("total_sequences must be positive")
    residues: dict[int, list[str]] = defaultdict(list)
    for row in coordinate_rows:
        state = row.get("hmm_match_state")
        if state is not None:
            residues[int(state)].append(str(row.get("residue", "")).upper())
    result = []
    for state in sorted(residues):
        observed = residues[state]
        padded = observed + [""] * (total_sequences - len(observed))
        conservation, identity, consensus, occupancy = shannon_conservation(padded)
        result.append({
            "hmm_match_state": state,
            "conservation": conservation,
            "identity": identity,
            "consensus": consensus,
            "occupancy": occupancy,
            "n_occupied": len(observed),
        })
    return result


def read_sequences(metadata_path: Path) -> list[dict[str, str]]:
    with metadata_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"internal_id", "sequence", "sequence_sha256", "sequence_qc_eligible"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Metadata must contain {sorted(required)}")
        rows = [
            row for row in reader
            if row["sequence_qc_eligible"].strip().casefold() == "true"
        ]
    if not rows:
        raise ValueError("No sequence-QC-eligible records")
    seen: set[str] = set()
    for row in rows:
        internal_id = row["internal_id"]
        sequence = row["sequence"].replace(" ", "").upper()
        if not internal_id or internal_id in seen:
            raise ValueError("Missing or duplicate internal ID")
        seen.add(internal_id)
        if not sequence or set(sequence) - CANONICAL_AA:
            raise ValueError(f"Invalid sequence for {internal_id}")
        if hashlib.sha256(sequence.encode("utf-8")).hexdigest() != row["sequence_sha256"]:
            raise ValueError(f"Sequence hash mismatch for {internal_id}")
        row["sequence"] = sequence
    return rows


def scan_sequences(
    hmm_path: Path,
    records: list[dict[str, str]],
    *,
    threads: int = 0,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    try:
        import pyhmmer
        from pyhmmer.easel import Alphabet, TextSequence
        from pyhmmer.plan7 import HMMFile
    except ImportError as exc:  # pragma: no cover - optional M3 environment
        raise RuntimeError("Install requirements-m3.txt before PF00741 HMM mapping") from exc
    with HMMFile(hmm_path) as handle:
        hmm = handle.read()
    if hmm is None:
        raise ValueError(f"No HMM found in {hmm_path}")
    alphabet = Alphabet.amino()
    sequences = [
        TextSequence(name=row["internal_id"].encode("utf-8"), sequence=row["sequence"]).digitize(alphabet)
        for row in records
    ]
    best_hits: dict[str, dict[str, Any]] = {}
    for top_hits in pyhmmer.hmmsearch([hmm], sequences, cpus=threads):
        for hit in top_hits:
            candidates = []
            for domain in hit.domains:
                alignment = domain.alignment
                candidates.append({
                    "score": float(domain.score),
                    "i_evalue": float(domain.i_evalue),
                    "hmm_from": int(alignment.hmm_from),
                    "hmm_to": int(alignment.hmm_to),
                    "target_from": int(alignment.target_from),
                    "target_to": int(alignment.target_to),
                    "envelope_from": int(domain.env_from),
                    "envelope_to": int(domain.env_to),
                    "hmm_sequence": as_text(alignment.hmm_sequence),
                    "target_sequence": as_text(alignment.target_sequence),
                })
            if candidates:
                best_hits[as_text(hit.name)] = max(candidates, key=lambda item: item["score"])
    hmm_info = {
        "name": as_text(hmm.name),
        "accession": as_text(hmm.accession) if hmm.accession else "",
        "length": int(hmm.M),
        "gathering_score": float(hmm.cutoffs.gathering1) if hmm.cutoffs.gathering_available else None,
        "pyhmmer_version": pyhmmer.__version__,
    }
    return best_hits, hmm_info


def write_mapping_outputs(
    records: list[dict[str, str]],
    hits: dict[str, dict[str, Any]],
    hmm_info: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    domain_rows: list[dict[str, Any]] = []
    coordinate_rows: list[dict[str, Any]] = []
    for record in records:
        internal_id = record["internal_id"]
        hit = hits.get(internal_id)
        if hit is None:
            domain_rows.append({
                "internal_id": internal_id, "hit": False, "score": "", "i_evalue": "",
                "hmm_from": "", "hmm_to": "", "target_from": "", "target_to": "",
                "n_mapped_residues": 0,
            })
            continue
        parsed = parse_alignment_columns(
            hit["hmm_sequence"], hit["target_sequence"], hit["hmm_from"], hit["target_from"]
        )
        for row in parsed:
            coordinate_rows.append({"internal_id": internal_id, **row})
        domain_rows.append({
            "internal_id": internal_id,
            "hit": True,
            "score": hit["score"],
            "i_evalue": hit["i_evalue"],
            "hmm_from": hit["hmm_from"],
            "hmm_to": hit["hmm_to"],
            "target_from": hit["target_from"],
            "target_to": hit["target_to"],
            "n_mapped_residues": len(parsed),
        })
    domain_path = output_dir / "pf00741_domains.csv"
    coordinate_path = output_dir / "pf00741_coordinate_map.csv"
    for path, rows in ((domain_path, domain_rows), (coordinate_path, coordinate_rows)):
        with path.open("w", encoding="utf-8", newline="") as handle:
            if not rows:
                raise ValueError(f"No rows generated for {path.name}")
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    conservation = conservation_by_state(coordinate_rows, len(records))
    conservation_path = output_dir / "pf00741_conservation.csv"
    with conservation_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(conservation[0]))
        writer.writeheader()
        writer.writerows(conservation)
    return {
        "hmm": hmm_info,
        "n_sequences": len(records),
        "n_hits": len(hits),
        "n_coordinate_rows": len(coordinate_rows),
        "domain_table_sha256": sha256_file(domain_path),
        "coordinate_map_sha256": sha256_file(coordinate_path),
        "conservation_sha256": sha256_file(conservation_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--hmm", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=0)
    args = parser.parse_args()
    records = read_sequences(args.metadata)
    hits, hmm_info = scan_sequences(args.hmm, records, threads=args.threads)
    receipt = write_mapping_outputs(records, hits, hmm_info, args.output_dir)
    receipt.update({
        "metadata_sha256": sha256_file(args.metadata),
        "hmm_sha256": sha256_file(args.hmm),
    })
    (args.output_dir / "pf00741_mapping_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
