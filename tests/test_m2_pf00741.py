"""Scientific coordinate and frozen-input failure cases for the M2 C workflow."""

import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from bioinformatics import m2_pf00741 as m2


class Domains(list):
    @property
    def included(self):
        return (domain for domain in self if domain.included)


def fake_hit(score, domain_scores, included=True):
    return SimpleNamespace(score=score, domains=Domains(
        [SimpleNamespace(score=value, included=included) for value in domain_scores]))


def test_ga_and_boundary_are_distinct_from_no_hit():
    assert m2.classify_hit(fake_hit(30, [28]), None, 25, 25, 2) == "accepted"
    assert m2.classify_hit(fake_hit(55, [28, 27]), None, 25, 25, 2) == "multi_hit"
    assert m2.classify_hit(None, fake_hit(24, [24]), 25, 25, 2) == "boundary"
    assert m2.classify_hit(None, fake_hit(20, [20]), 25, 25, 2) == "weak_hit"
    assert m2.classify_hit(None, None, 25, 25, 2) == "no_hit"
    assert m2.classify_hit(fake_hit(28, [28], included=False), None, 25, 25, 2) == "failed"


def test_alignment_checks_full_residue_and_insertion_state():
    alignment = SimpleNamespace(hmm_sequence="AB.CD", target_sequence="AXY-D",
                                hmm_from=5, target_from=10, target_to=13)
    domain = SimpleNamespace(alignment=alignment)
    rows = m2.alignment_rows(domain, "Q" * 9 + "AXYD", 10)
    assert [(r["raw_position"], r["hmm_match_state"]) for r in rows] == [
        (10, 5), (11, 6), (12, None), (13, 8)]
    with pytest.raises(ValueError, match="full input sequence"):
        m2.alignment_rows(domain, "Q" * 9 + "AXWD", 10)
    with pytest.raises(ValueError, match="out of range"):
        m2.alignment_rows(domain, "Q" * 9 + "AXYD", 7)


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_frozen_handoff_checks_sequence_and_split_identity(tmp_path, monkeypatch):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    def sequence(index):
        digits = []
        for _ in range(3):
            index, remainder = divmod(index, 20)
            digits.append(alphabet[remainder])
        return "A" + "".join(digits) + "C"

    metadata = tmp_path / "metadata.csv"
    fasta = tmp_path / "sequences_for_clustering.fasta"
    split = tmp_path / "split_manifest.csv"
    audit_path = tmp_path / "audit.json"
    split_summary_path = tmp_path / "split_summary.json"
    metadata_rows, split_rows, fasta_lines = [], [], []
    for index in range(2078):
        identifier = f"GVPA_{index+1:06d}"
        seq = sequence(index)
        eligible = index < 2076
        subset = "discovery" if index < 1453 else "validation"
        metadata_rows.append({"internal_id": identifier, "sequence": seq,
                              "sequence_length": len(seq),
                              "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
                              "sequence_qc_eligible": str(eligible).lower(),
                              "primary_analysis_eligible": str(eligible).lower(),
                              "analysis_cohort": "primary" if eligible else "excluded_sequence_qc"})
        if eligible:
            fasta_lines.extend([f">{identifier}", seq])
            split_rows.append({"internal_id": identifier, "split": subset,
                               "homology_cluster": f"CLUSTER_{index}",
                               "analysis_cohort": "primary",
                               "primary_analysis_eligible": "true"})
    write_csv(metadata, list(metadata_rows[0]), metadata_rows)
    fasta.write_text("\n".join(fasta_lines) + "\n", encoding="ascii")
    write_csv(split, list(split_rows[0]), split_rows)

    def update_references():
        audit = {"dataset_version": m2.DATA_VERSION, "unique_internal_ids": 2078,
                 "sequence_qc_eligible": 2076,
                 "output_sha256": {"metadata.csv": m2.digest(metadata),
                                   "sequences_for_clustering.fasta": m2.digest(fasta)}}
        summary = {"split_version": m2.SPLIT_VERSION,
                   "split_manifest_sha256": m2.digest(split),
                   "sequence_counts": {"discovery": 1453, "validation": 623},
                   "cluster_counts": {"discovery": 1453, "validation": 623},
                   "primary_sequence_counts": {"discovery": 1453, "validation": 623}}
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        split_summary_path.write_text(json.dumps(summary), encoding="utf-8")
        monkeypatch.setattr(m2, "FROZEN_AUDIT_SHA256", m2.digest(audit_path))
        monkeypatch.setattr(m2, "FROZEN_SPLIT_SHA256", m2.digest(split_summary_path))

    update_references()
    records, provenance = m2.validate_b_handoff(metadata, fasta, split, audit_path, split_summary_path)
    assert len(records) == 2076
    assert provenance["n_sequences"] == 2076

    # Preserve all file/reference hashes while introducing a semantic sequence mismatch.
    fasta_lines[1] = "CCCCA"
    fasta.write_text("\n".join(fasta_lines) + "\n", encoding="ascii")
    update_references()
    with pytest.raises(ValueError, match="Sequence identity"):
        m2.validate_b_handoff(metadata, fasta, split, audit_path, split_summary_path)

    fasta_lines[1] = metadata_rows[0]["sequence"]
    fasta.write_text("\n".join(fasta_lines) + "\n", encoding="ascii")
    split_rows[1]["homology_cluster"] = split_rows[0]["homology_cluster"]
    split_rows[1]["split"] = "validation"
    write_csv(split, list(split_rows[0]), split_rows)
    update_references()
    with pytest.raises(ValueError, match="crosses discovery/validation"):
        m2.validate_b_handoff(metadata, fasta, split, audit_path, split_summary_path)


def test_structure_sequence_checks_deposited_count(tmp_path):
    pdb = tmp_path / "7R1C.pdb"
    pdb.write_text("SEQRES   1 N    3  ALA CYS ASP\n", encoding="ascii")
    assert m2.deposited_sequence(pdb, "N") == "ACD"
    pdb.write_text("SEQRES   1 N    4  ALA CYS ASP\n", encoding="ascii")
    with pytest.raises(ValueError, match="SEQRES"):
        m2.deposited_sequence(pdb, "N")


def test_source_receipt_needs_matching_hash_and_url(tmp_path):
    hmm = tmp_path / "PF00741.hmm"
    hmm.write_bytes(b"test")
    receipt = tmp_path / "receipt.json"
    payload = {"url": m2.HMM_URL, "sha256": m2.digest(hmm),
               "retrieved_utc": "2026-09-23T00:00:00Z"}
    receipt.write_text(json.dumps(payload), encoding="utf-8")
    assert m2.check_source_receipt(receipt, hmm, m2.HMM_URL) == payload
    hmm.write_bytes(b"changed")
    with pytest.raises(ValueError, match="receipt"):
        m2.check_source_receipt(receipt, hmm, m2.HMM_URL)


def test_full_scan_assigns_all_2076_synthetic_statuses_without_validation_conservation(tmp_path, monkeypatch):
    pyhmmer = pytest.importorskip("pyhmmer")
    alphabet = pyhmmer.easel.Alphabet.amino()
    motif = "ACDEFGHIKLMNPQRSTVWY"
    digital = pyhmmer.easel.TextSequence(name=b"toy", sequence=motif).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(digital, pyhmmer.plan7.Background(alphabet))
    hmm.accession = b"PF00741.1"
    hmm.cutoffs.gathering = (5.0, 5.0)
    hmm_path = tmp_path / "toy.hmm"
    with hmm_path.open("wb") as handle:
        hmm.write(handle)
    monkeypatch.setattr(m2, "HMM_ACCESSION", "PF00741.1")
    monkeypatch.setattr(m2, "HMM_SHA256", m2.digest(hmm_path))
    receipt_path = tmp_path / "download_receipt.json"
    receipt_path.write_text(json.dumps({"url": m2.HMM_URL, "sha256": m2.digest(hmm_path),
                                        "retrieved_utc": "2026-09-23T00:00:00Z",
                                        "model_accession": "PF00741.1"}), encoding="utf-8")
    records = []
    for index in range(2076):
        sequence = motif if index == 0 else "A" * len(motif)
        records.append({"internal_id": f"GVPA_{index+1:06d}", "sequence": sequence,
                        "sequence_sha256": m2.sequence_digest(sequence),
                        "split": "discovery" if index < 1453 else "validation",
                        "analysis_cohort": "primary", "primary_analysis_eligible": "true"})
    monkeypatch.setattr(m2, "validate_b_handoff",
                        lambda *_: (records, {"n_sequences": 2076, "dataset_version": "synthetic-test"}))
    args = SimpleNamespace(metadata=tmp_path / "metadata.csv", fasta=tmp_path / "input.fasta",
                           split_manifest=tmp_path / "split.csv", audit_summary=tmp_path / "audit.json",
                           split_summary=tmp_path / "summary.json", hmm=hmm_path,
                           hmm_receipt=receipt_path, out_dir=tmp_path / "raw",
                           scan_summary=tmp_path / "scan_summary.json",
                           coordinate_summary=tmp_path / "coordinate_summary.json",
                           threads=1, boundary_bits=2.0)
    result = m2.run_scan(args)
    assert result["n_sequences"] == 2076
    assert sum(result["counts"].values()) == 2076
    assert result["by_split"]["discovery"] and result["by_split"]["validation"]
    assert result["parameters"]["validation_used_for_conservation"] is False
    coordinate = json.loads(args.coordinate_summary.read_text(encoding="utf-8"))
    assert coordinate["n_coordinate_rows"] > 0
    assert coordinate["n_distinct_mapped_sequences"] >= 1
    with (args.out_dir / "per_sequence_status.csv").open(encoding="utf-8", newline="") as handle:
        statuses = list(csv.DictReader(handle))
    assert len(statuses) == 2076
    assert set(row["status"] for row in statuses) <= {
        "accepted", "multi_hit", "boundary", "weak_hit", "no_hit", "failed"}
