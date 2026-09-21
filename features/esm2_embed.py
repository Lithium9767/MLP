#!/usr/bin/env python3
"""Extract split-aware ESM-2 residue embeddings for the GvpA pipeline.

The validation split is locked by default. Embeddings are cached outside Git by
sequence hash and model ID; raw attention is intentionally not treated as a
functional attribution score in this module.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


MODEL_ID = "facebook/esm2_t12_35M_UR50D"
CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class SequenceRecord:
    internal_id: str
    sequence_id: str
    sequence: str
    sequence_sha256: str
    split: str
    analysis_cohort: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalise_sequence(value: str) -> str:
    sequence = re.sub(r"\s+", "", str(value)).upper()
    invalid = sorted(set(sequence) - CANONICAL_AA)
    if not sequence:
        raise ValueError("Empty protein sequence")
    if invalid:
        raise ValueError(f"Non-canonical residues: {''.join(invalid)}")
    return sequence


def read_embedding_records(
    metadata_path: Path,
    split_path: Path,
    *,
    selected_split: str = "discovery",
    primary_only: bool = True,
) -> list[SequenceRecord]:
    with metadata_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "internal_id", "sequence_id", "sequence", "sequence_sha256",
            "sequence_qc_eligible", "primary_analysis_eligible", "analysis_cohort",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Metadata must contain {sorted(required)}")
        metadata_rows = list(reader)
    metadata = {row["internal_id"]: row for row in metadata_rows}
    if "" in metadata or len(metadata) != len(metadata_rows):
        raise ValueError("Metadata contains empty or duplicate internal IDs")

    with split_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"internal_id", "sequence_id", "split"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Split manifest must contain {sorted(required)}")
        split_rows = list(reader)
    split_by_id = {row["internal_id"]: row for row in split_rows}
    if "" in split_by_id or len(split_by_id) != len(split_rows):
        raise ValueError("Split manifest contains empty or duplicate internal IDs")

    records: list[SequenceRecord] = []
    for internal_id, split_row in split_by_id.items():
        if split_row["split"] != selected_split:
            continue
        row = metadata.get(internal_id)
        if row is None:
            raise ValueError(f"Split manifest references unknown internal ID {internal_id}")
        if row["sequence_id"] != split_row["sequence_id"]:
            raise ValueError(f"Sequence ID mismatch for {internal_id}")
        if row["sequence_qc_eligible"].casefold() != "true":
            raise ValueError(f"Split manifest includes sequence-QC-ineligible ID {internal_id}")
        if primary_only and row["primary_analysis_eligible"].casefold() != "true":
            continue
        sequence = normalise_sequence(row["sequence"])
        digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        if digest != row["sequence_sha256"]:
            raise ValueError(f"Sequence hash mismatch for {internal_id}")
        records.append(
            SequenceRecord(
                internal_id=internal_id,
                sequence_id=row["sequence_id"],
                sequence=sequence,
                sequence_sha256=digest,
                split=selected_split,
                analysis_cohort=row["analysis_cohort"],
            )
        )
    if not records:
        raise ValueError(f"No records selected for split={selected_split!r}")
    return sorted(records, key=lambda record: record.internal_id)


def iter_windows(sequence_length: int, window_length: int = 30, step: int = 5) -> list[tuple[int, int]]:
    if window_length <= 0 or step <= 0:
        raise ValueError("Window length and step must be positive")
    if sequence_length < window_length:
        return []
    return [
        (start + 1, start + window_length)
        for start in range(0, sequence_length - window_length + 1, step)
    ]


def model_cache_key(model_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", model_id).strip("_")


def load_backend(model_id: str, device: str | None = None):
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:  # pragma: no cover - depends on optional M3 environment
        raise RuntimeError("Install requirements-m3.txt before extracting ESM-2 embeddings") from exc
    resolved_device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = AutoModel.from_pretrained(model_id)
    model.to(resolved_device)
    model.eval()
    return tokenizer, model, torch, resolved_device


def _batches(values: list[SequenceRecord], size: int) -> Iterable[list[SequenceRecord]]:
    if size <= 0:
        raise ValueError("Batch size must be positive")
    for start in range(0, len(values), size):
        yield values[start:start + size]


def embed_records(
    records: list[SequenceRecord],
    *,
    output_dir: Path,
    model_id: str = MODEL_ID,
    batch_size: int = 4,
    device: str | None = None,
) -> dict[str, Any]:
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - depends on optional M3 environment
        raise RuntimeError("Install requirements-m3.txt before extracting ESM-2 embeddings") from exc
    tokenizer, model, torch, resolved_device = load_backend(model_id, device)
    cache_dir = output_dir / model_cache_key(model_id)
    cache_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str]] = []

    for batch in _batches(records, batch_size):
        sequences = [record.sequence for record in batch]
        encoded = tokenizer(
            sequences,
            return_tensors="pt",
            padding=True,
            add_special_tokens=True,
            truncation=False,
        )
        encoded = {name: tensor.to(resolved_device) for name, tensor in encoded.items()}
        with torch.no_grad():
            output = model(**encoded, output_attentions=False)
        for index, record in enumerate(batch):
            residue_embedding = (
                output.last_hidden_state[index, 1:1 + len(record.sequence)]
                .detach().cpu().numpy().astype(np.float32)
            )
            if residue_embedding.shape[0] != len(record.sequence):
                raise RuntimeError(f"Unexpected token/residue length for {record.internal_id}")
            path = cache_dir / f"{record.sequence_sha256}.npz"
            np.savez_compressed(
                path,
                residue_embedding=residue_embedding,
                full_length_mean=residue_embedding.mean(axis=0),
                sequence_sha256=np.array(record.sequence_sha256),
                model_id=np.array(model_id),
            )
            manifest_rows.append({
                "internal_id": record.internal_id,
                "sequence_id": record.sequence_id,
                "split": record.split,
                "analysis_cohort": record.analysis_cohort,
                "sequence_length": str(len(record.sequence)),
                "sequence_sha256": record.sequence_sha256,
                "model_id": model_id,
                "embedding_path": str(path),
                "embedding_sha256": sha256_file(path),
                "residue_dimension": str(residue_embedding.shape[1]),
            })

    manifest_path = output_dir / "embedding_manifest.csv"
    fields = list(manifest_rows[0])
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(manifest_rows)
    return {
        "model_id": model_id,
        "device": resolved_device,
        "n_sequences": len(records),
        "split": records[0].split,
        "primary_only": all(record.analysis_cohort == "primary" for record in records),
        "embedding_manifest_sha256": sha256_file(manifest_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--split", choices=("discovery", "validation"), default="discovery")
    parser.add_argument("--allow-validation", action="store_true")
    parser.add_argument("--include-sensitivity", action="store_true")
    parser.add_argument("--model-id", default=MODEL_ID)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--device")
    args = parser.parse_args()
    if args.split == "validation" and not args.allow_validation:
        parser.error("Validation is locked; pass --allow-validation only after discovery rules are frozen")
    records = read_embedding_records(
        args.metadata,
        args.split_manifest,
        selected_split=args.split,
        primary_only=not args.include_sensitivity,
    )
    receipt = embed_records(
        records,
        output_dir=args.output_dir,
        model_id=args.model_id,
        batch_size=args.batch_size,
        device=args.device,
    )
    receipt.update({
        "metadata_sha256": sha256_file(args.metadata),
        "split_manifest_sha256": sha256_file(args.split_manifest),
    })
    (args.output_dir / "embedding_run.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
