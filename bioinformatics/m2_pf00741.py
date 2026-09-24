"""Audited M2 PF00741 scan and 7R1C mapping; large outputs stay outside Git.

The `scan` command requires the complete B handoff. `structure` can run
independently against the public PF00741 model and PDB 7R1C.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from .m2_coordinates import AA3, parse_alignment_columns, parse_pdb_chain


CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")
DATA_VERSION = "gvpa-recognition-c7f6f005d717"
SPLIT_VERSION = "homology-8b9005e2d9-s42"
HMM_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00741?annotation=hmm"
PDB_URL = "https://files.rcsb.org/download/7R1C.pdb"
HMM_ACCESSION = "PF00741.24"
HMM_SHA256 = "e44f933c618426c2d0465657a64f1890f50feaeaea626ecb0b8dedb734f38c2f"
PDB_SHA256 = "c83776e7091d1342a4c6dbc40fbcc7be2e2beff87b12f1c52c9d6de5be0be870"
SPLITS = frozenset({"discovery", "validation"})
FROZEN_AUDIT_SHA256 = "908ce3dbf0e38731980c166f83b3921956097bc1498aff091e625084aacc449c"
FROZEN_SPLIT_SHA256 = "e192f67a4e5534552cdcdfcb5346a23c4bcc365f13569f534be86fe850e6efb3"


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def sequence_digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def load_json(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return data


def csv_rows(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or repeated CSV headers: {path}")
        if not required <= set(reader.fieldnames):
            raise ValueError(f"Missing {sorted(required - set(reader.fieldnames))} in {path}")
        rows = list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError(f"Empty or malformed CSV: {path}")
    return rows


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = None
    for number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        if line.startswith(">"):
            identifier = line[1:].strip()
            if not identifier or any(char.isspace() for char in identifier) or identifier in records:
                raise ValueError(f"Invalid or repeated FASTA ID on line {number}")
            records[identifier] = ""
            current = identifier
        elif line:
            if current is None or set(line) - CANONICAL:
                raise ValueError(f"Invalid FASTA sequence on line {number}")
            records[current] += line
    if not records or any(not sequence for sequence in records.values()):
        raise ValueError("FASTA has no records or contains an empty sequence")
    return records


def validate_b_handoff(metadata_path: Path, fasta_path: Path, split_path: Path,
                       audit_path: Path, split_summary_path: Path) -> tuple[list[dict], dict]:
    """Require exact frozen file hashes and one-to-one QC/split/FASTA coverage."""
    if digest(audit_path) != FROZEN_AUDIT_SHA256 or digest(split_summary_path) != FROZEN_SPLIT_SHA256:
        raise ValueError("Audit/split references differ from the committed frozen M2 summaries")
    audit = load_json(audit_path)
    split_summary = load_json(split_summary_path)
    if audit.get("dataset_version") != DATA_VERSION or split_summary.get("split_version") != SPLIT_VERSION:
        raise ValueError("Dataset or split version is not the frozen M2 version")
    expected = {
        metadata_path: audit["output_sha256"]["metadata.csv"],
        fasta_path: audit["output_sha256"]["sequences_for_clustering.fasta"],
        split_path: split_summary["split_manifest_sha256"],
    }
    for path, sha in expected.items():
        if digest(path) != sha:
            raise ValueError(f"B handoff hash mismatch: {path}")
    metadata = csv_rows(metadata_path, {"internal_id", "sequence", "sequence_length",
                                        "sequence_sha256", "sequence_qc_eligible",
                                        "primary_analysis_eligible", "analysis_cohort"})
    if len(metadata) != audit["unique_internal_ids"] or len(metadata) != 2078:
        raise ValueError("Metadata record count differs from frozen audit")
    by_id: dict[str, dict] = {}
    for row in metadata:
        identifier = row["internal_id"]
        if not identifier or identifier in by_id:
            raise ValueError("Missing or duplicate internal ID in metadata")
        by_id[identifier] = row
    eligible = {identifier: row for identifier, row in by_id.items()
                if row["sequence_qc_eligible"] == "true"}
    if len(eligible) != audit["sequence_qc_eligible"] or len(eligible) != 2076:
        raise ValueError("QC-eligible count differs from frozen audit")
    fasta = read_fasta(fasta_path)
    if set(fasta) != set(eligible):
        raise ValueError("FASTA and QC-eligible metadata IDs differ")
    for identifier, row in eligible.items():
        sequence = row["sequence"]
        if (not sequence or set(sequence) - CANONICAL or sequence != fasta[identifier]
                or len(sequence) != int(row["sequence_length"])
                or sequence_digest(sequence) != row["sequence_sha256"]):
            raise ValueError(f"Sequence identity, length or SHA-256 mismatch for {identifier}")
    split_rows = csv_rows(split_path, {"internal_id", "split", "homology_cluster",
                                        "analysis_cohort", "primary_analysis_eligible"})
    if len(split_rows) != len(eligible):
        raise ValueError("Split manifest has the wrong row count")
    assignments: dict[str, dict] = {}
    cluster_split: dict[str, str] = {}
    counts = Counter()
    clusters: dict[str, set[str]] = defaultdict(set)
    primary = Counter()
    for row in split_rows:
        identifier, subset, cluster = row["internal_id"], row["split"], row["homology_cluster"]
        if identifier not in eligible or identifier in assignments or subset not in SPLITS or not cluster:
            raise ValueError("Unknown/duplicate ID, split name or missing cluster")
        source = eligible[identifier]
        if (row["analysis_cohort"] != source["analysis_cohort"]
                or row["primary_analysis_eligible"] != source["primary_analysis_eligible"]):
            raise ValueError(f"Split and metadata fields disagree for {identifier}")
        if cluster in cluster_split and cluster_split[cluster] != subset:
            raise ValueError(f"Homology cluster crosses discovery/validation: {cluster}")
        cluster_split[cluster] = subset
        clusters[subset].add(cluster)
        counts[subset] += 1
        primary[subset] += row["primary_analysis_eligible"] == "true"
        assignments[identifier] = row
    if set(assignments) != set(eligible):
        raise ValueError("Split does not cover every QC-eligible ID")
    if (dict(counts) != split_summary["sequence_counts"]
            or {key: len(value) for key, value in clusters.items()} != split_summary["cluster_counts"]
            or dict(primary) != split_summary["primary_sequence_counts"]):
        raise ValueError("Split counts do not match frozen M2 summary")
    ordered = [{**eligible[identifier], "split": assignments[identifier]["split"],
                "homology_cluster": assignments[identifier]["homology_cluster"]}
               for identifier in sorted(eligible)]
    return ordered, {"dataset_version": DATA_VERSION, "split_version": SPLIT_VERSION,
                     "metadata_sha256": expected[metadata_path], "fasta_sha256": expected[fasta_path],
                     "split_manifest_sha256": expected[split_path],
                     "audit_summary_sha256": digest(audit_path),
                     "split_summary_sha256": digest(split_summary_path),
                     "n_sequences": len(ordered), "n_clusters": len(cluster_split)}


def load_hmm(path: Path):
    try:
        import pyhmmer
        from pyhmmer.plan7 import HMMFile
    except ImportError as error:
        raise RuntimeError("Install PyHMMER before running this command") from error
    with HMMFile(path) as handle:
        hmm = handle.read()
        if hmm is None or handle.read() is not None:
            raise ValueError("Expected exactly one profile HMM")
    accession = (hmm.accession.decode("ascii") if isinstance(hmm.accession, bytes)
                 else hmm.accession or "")
    if accession != HMM_ACCESSION or digest(path) != HMM_SHA256:
        raise ValueError(f"PF00741 model accession/hash differs from frozen {HMM_ACCESSION}")
    if not hmm.cutoffs.gathering_available:
        raise ValueError("PF00741 HMM has no sequence/domain gathering cutoffs")
    return hmm, {"accession": accession, "length": int(hmm.M),
                 "ga_sequence_bits": float(hmm.cutoffs.gathering1),
                 "ga_domain_bits": float(hmm.cutoffs.gathering2),
                 "sha256": digest(path), "source_url": HMM_URL,
                 "pyhmmer_version": pyhmmer.__version__}


def scan_hits(hmm, records: list[dict], threads: int):
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextSequence
    alphabet = Alphabet.amino()
    sequences = [TextSequence(name=row["internal_id"].encode("ascii"),
                              sequence=row["sequence"]).digitize(alphabet)
                 for row in records]
    # GA gives the official family cutoff; a second permissive search records
    # below-GA candidates so a near-threshold result is not silently a no-hit.
    accepted = next(pyhmmer.hmmsearch([hmm], sequences, cpus=threads, bit_cutoffs="gathering"))
    permissive = next(pyhmmer.hmmsearch([hmm], sequences, cpus=threads,
                                        T=0.0, domT=0.0, incT=0.0, incdomT=0.0))
    def keyed(top_hits):
        result = {}
        for hit in top_hits:
            identifier = hit.name.decode("ascii") if isinstance(hit.name, bytes) else hit.name
            if identifier in result:
                raise ValueError(f"Duplicate HMM hit for {identifier}")
            result[identifier] = hit
        return result
    ga, weak = keyed(accepted), keyed(permissive)
    expected = {row["internal_id"] for row in records}
    if not set(ga) <= expected or not set(weak) <= expected:
        raise ValueError("HMM search returned an unknown internal ID")
    return ga, weak


def alignment_rows(domain, sequence: str, hmm_length: int) -> list[dict]:
    alignment = domain.alignment
    rows = parse_alignment_columns(alignment.hmm_sequence, alignment.target_sequence,
                                   int(alignment.hmm_from), int(alignment.target_from))
    if not rows:
        raise ValueError("Reported HMM domain has no aligned residues")
    for row in rows:
        position, state = row["raw_position"], row["hmm_match_state"]
        if position > len(sequence) or sequence[position - 1] != row["residue"]:
            raise ValueError("HMM alignment does not match the full input sequence")
        if state is not None and not 1 <= state <= hmm_length:
            raise ValueError("HMM match state is out of range")
    if rows[-1]["raw_position"] != alignment.target_to:
        raise ValueError("HMM alignment terminal coordinate disagrees with target_to")
    return rows


def git_state() -> dict:
    def run(*args):
        result = subprocess.run(["git", *args], capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else ""
    return {"commit": run("rev-parse", "HEAD"), "dirty": bool(run("status", "--porcelain"))}


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def check_source_receipt(path: Path, source_path: Path, url: str) -> dict:
    content = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(content, list):
        if len(content) != 1:
            raise ValueError("Expected a single downloaded source receipt")
        content = content[0]
    if (not isinstance(content, dict) or content.get("url", content.get("source_url")) != url
            or content.get("sha256") != digest(source_path)
            or not content.get("retrieved_utc")):
        raise ValueError(f"Source receipt does not match {source_path}")
    return content


def hit_detail(hit) -> dict:
    domains = []
    for domain in hit.domains:
        alignment = domain.alignment
        domains.append({
            "score_bits": float(domain.score), "i_evalue": float(domain.i_evalue),
            "included": bool(domain.included), "reported": bool(domain.reported),
            "hmm_from": int(alignment.hmm_from), "hmm_to": int(alignment.hmm_to),
            "target_from": int(alignment.target_from), "target_to": int(alignment.target_to),
            "hmm_sequence": str(alignment.hmm_sequence),
            "target_sequence": str(alignment.target_sequence),
        })
    return {"score_bits": float(hit.score), "evalue": float(hit.evalue),
            "included": bool(hit.included), "reported": bool(hit.reported),
            "domains": domains}


def classify_hit(ga_hit, permissive_hit, ga_sequence_bits: float, ga_domain_bits: float,
                 boundary_bits: float) -> str:
    if ga_hit is not None:
        included_domains = list(ga_hit.domains.included)
        if not included_domains:
            return "failed"
        return "multi_hit" if len(included_domains) > 1 else "accepted"
    if permissive_hit is None:
        return "no_hit"
    domain_scores = [float(domain.score) for domain in permissive_hit.domains]
    distance = min(float(permissive_hit.score) - ga_sequence_bits,
                   max(domain_scores, default=float("-inf")) - ga_domain_bits)
    return "boundary" if distance >= -boundary_bits else "weak_hit"


def run_scan(args) -> dict:
    if args.threads < 1 or not math.isfinite(args.boundary_bits) or args.boundary_bits < 0:
        raise ValueError("Threads must be positive and boundary distance nonnegative")
    records, handoff = validate_b_handoff(args.metadata, args.fasta, args.split_manifest,
                                           args.audit_summary, args.split_summary)
    hmm, hmm_info = load_hmm(args.hmm)
    hmm_receipt = check_source_receipt(args.hmm_receipt, args.hmm, HMM_URL)
    if hmm_receipt.get("model_accession") != hmm_info["accession"]:
        raise ValueError("HMM download receipt and model accession disagree")
    if args.out_dir.exists() or args.scan_summary.exists() or args.coordinate_summary.exists():
        raise FileExistsError("Refusing to overwrite an existing M2 scan output")
    args.out_dir.mkdir(parents=True)
    try:
        ga, weak = scan_hits(hmm, records, args.threads)
        status_rows = []
        coordinate_rows = []
        raw_hits = []
        counts_by_split = {subset: Counter() for subset in sorted(SPLITS)}
        coverage_by_split = {subset: {"accepted_residues": 0, "accepted_match_states": 0,
                                      "sequence_residues": 0} for subset in sorted(SPLITS)}
        failures = []
        for record in records:
            identifier, sequence, subset = record["internal_id"], record["sequence"], record["split"]
            ga_hit, weak_hit = ga.get(identifier), weak.get(identifier)
            status = classify_hit(ga_hit, weak_hit, hmm_info["ga_sequence_bits"],
                                  hmm_info["ga_domain_bits"], args.boundary_bits)
            if ga_hit is not None:
                raw_hits.append({"internal_id": identifier, "search": "gathering", "hit": hit_detail(ga_hit)})
            if weak_hit is not None:
                raw_hits.append({"internal_id": identifier, "search": "permissive_T0_domT0",
                                 "hit": hit_detail(weak_hit)})
            mapped: set[int] = set()
            states: set[int] = set()
            deletion_states: set[int] = set()
            included_domains = list(ga_hit.domains.included) if ga_hit is not None else []
            if status in {"accepted", "multi_hit"}:
                try:
                    for number, domain in enumerate(included_domains, 1):
                        rows = alignment_rows(domain, sequence, hmm_info["length"])
                        match_states = {row["hmm_match_state"] for row in rows
                                        if row["hmm_match_state"] is not None}
                        deletion_states.update(set(range(int(domain.alignment.hmm_from),
                                                         int(domain.alignment.hmm_to) + 1)) - match_states)
                        for row in rows:
                            mapped.add(row["raw_position"])
                            if row["hmm_match_state"] is not None:
                                states.add(row["hmm_match_state"])
                            coordinate_rows.append({
                                "internal_id": identifier, "sequence_sha256": record["sequence_sha256"],
                                "split": subset, "analysis_cohort": record["analysis_cohort"],
                                "domain_number": number, "raw_position": row["raw_position"],
                                "residue": row["residue"], "hmm_match_state": row["hmm_match_state"],
                                "is_insertion": row["is_insertion"],
                            })
                except (ValueError, KeyError, TypeError) as error:
                    status = "failed"
                    failures.append({"internal_id": identifier, "reason": str(error)})
                    coordinate_rows = [row for row in coordinate_rows if row["internal_id"] != identifier]
                    mapped.clear()
                    states.clear()
                    deletion_states.clear()
            counts_by_split[subset][status] += 1
            if status in {"accepted", "multi_hit"}:
                coverage_by_split[subset]["accepted_residues"] += len(mapped)
                coverage_by_split[subset]["accepted_match_states"] += len(states)
            coverage_by_split[subset]["sequence_residues"] += len(sequence)
            status_rows.append({
                "internal_id": identifier, "sequence_sha256": record["sequence_sha256"],
                "split": subset, "analysis_cohort": record["analysis_cohort"],
                "is_partial_or_conflict": record["primary_analysis_eligible"] != "true",
                "status": status, "ga_sequence_bits": hmm_info["ga_sequence_bits"],
                "ga_domain_bits": hmm_info["ga_domain_bits"],
                "best_score_bits": (float(ga_hit.score) if ga_hit is not None else
                                    float(weak_hit.score) if weak_hit is not None else ""),
                "n_ga_domains": len(included_domains), "n_mapped_residues": len(mapped),
                "n_match_states": len(states), "n_deleted_match_states": len(deletion_states),
                "residue_coverage": len(mapped) / len(sequence) if mapped else "",
                "hmm_coverage": len(states) / hmm_info["length"] if states else "",
            })
        if len(status_rows) != 2076:
            raise ValueError("Not every QC-eligible input received a scan status")
        raw_dir = args.out_dir
        write_csv(raw_dir / "per_sequence_status.csv", list(status_rows[0]), status_rows)
        write_csv(raw_dir / "hmm_coordinate_map.csv",
                  ["internal_id", "sequence_sha256", "split", "analysis_cohort", "domain_number",
                   "raw_position", "residue", "hmm_match_state", "is_insertion"], coordinate_rows)
        with (raw_dir / "raw_hits.jsonl").open("w", encoding="utf-8") as handle:
            for row in raw_hits:
                handle.write(json.dumps(row, ensure_ascii=False) + "\n")
        write_json(raw_dir / "failures.json", {"records": failures})
        totals = Counter(row["status"] for row in status_rows)
        base = {"schema_version": "1.0", "run_utc": datetime.now(timezone.utc).isoformat(),
                "purpose": "M2 real PF00741 family/coordinate audit, not functional labels",
                "git": git_state(), "input": handoff, "hmm": hmm_info,
                "command": ["python", "-m", "bioinformatics.m2_pf00741", *sys.argv[1:]],
                "implementation_sha256": {"workflow": digest(Path(__file__)),
                                          "coordinates": digest(Path(__file__).with_name("m2_coordinates.py"))},
                "hmm_download_receipt_sha256": digest(args.hmm_receipt),
                "software": {"python": sys.version.split()[0],
                             "pyhmmer": hmm_info["pyhmmer_version"]},
                "parameters": {"threads": args.threads, "search": "GA and permissive T=0/domT=0",
                               "boundary_distance_bits": args.boundary_bits,
                               "validation_used_for_conservation": False},
                "raw_output_sha256": {p.name: digest(p) for p in sorted(raw_dir.iterdir()) if p.is_file()},
                "raw_output_storage": "ignored local/shared directory; request files by hashes"}
        scan_summary = {**base, "n_sequences": len(status_rows), "counts": dict(totals),
                        "by_split": {key: dict(value) for key, value in counts_by_split.items()},
                        "n_multihit": totals["multi_hit"], "n_failed": totals["failed"],
                        "limitations": ["No PF00741 hit is a functional negative label",
                                        "Permissive near-threshold search reports candidates, not accepted family hits",
                                        "Only GA-included domains receive HMM coordinates",
                                        "No conservation estimate uses validation sequences"]}
        coordinate_summary = {**base, "n_coordinate_rows": len(coordinate_rows),
                              "per_split_coverage": coverage_by_split,
                              "n_insertions": sum(row["is_insertion"] for row in coordinate_rows),
                              "n_distinct_mapped_sequences": len({row["internal_id"] for row in coordinate_rows}),
                              "coordinate_convention": "1-based full natural sequence position; insertion has null HMM match state"}
        write_json(args.scan_summary, scan_summary)
        write_json(args.coordinate_summary, coordinate_summary)
        return scan_summary
    except Exception:
        (args.out_dir / "FAILED.txt").write_text("Incomplete scan; do not use partial outputs.\n", encoding="utf-8")
        raise


def deposited_sequence(path: Path, chain: str) -> str:
    residues = []
    declared_lengths = set()
    for line in path.read_text(encoding="latin-1").splitlines():
        if line.startswith("SEQRES") and line.split()[2] == chain:
            fields = line.split()
            declared_lengths.add(int(fields[3]))
            residues.extend(AA3.get(code.upper(), "X") for code in fields[4:])
    if len(declared_lengths) != 1 or len(residues) != declared_lengths.pop() or set(residues) - CANONICAL:
        raise ValueError("PDB SEQRES is missing, inconsistent or noncanonical")
    return "".join(residues)


def run_structure(args) -> dict:
    if args.threads < 1 or args.chain != "N":
        raise ValueError("7R1C M2 mapping requires chain N and positive thread count")
    if args.out_dir.exists() or args.summary.exists():
        raise FileExistsError("Refusing to overwrite an existing structure output")
    check_source_receipt(args.pdb_receipt, args.pdb, PDB_URL)
    if digest(args.pdb) != PDB_SHA256:
        raise ValueError("7R1C PDB bytes differ from the M2 frozen source")
    hmm, hmm_info = load_hmm(args.hmm)
    hmm_receipt = check_source_receipt(args.hmm_receipt, args.hmm, HMM_URL)
    if hmm_receipt.get("model_accession") != hmm_info["accession"]:
        raise ValueError("HMM accession differs from its source receipt")
    header = args.pdb.read_text(encoding="latin-1").splitlines()[0]
    if not header.startswith("HEADER") or "7R1C" not in header:
        raise ValueError("PDB source is not 7R1C")
    deposited = deposited_sequence(args.pdb, args.chain)
    modeled = parse_pdb_chain(args.pdb, args.chain)
    modeled_by_number = {}
    for row in modeled:
        if (row.insertion_code or row.residue_number not in range(1, len(deposited) + 1)
                or deposited[row.residue_number - 1] != row.amino_acid
                or row.residue_number in modeled_by_number):
            raise ValueError("PDB modeled residue does not match deposited sequence numbering")
        modeled_by_number[row.residue_number] = row
    ga, _ = scan_hits(hmm, [{"internal_id": "7R1C_N", "sequence": deposited}], args.threads)
    hit = ga.get("7R1C_N")
    if hit is None:
        raise ValueError("7R1C deposited sequence has no GA-qualified PF00741 hit")
    domains = list(hit.domains.included)
    if len(domains) != 1:
        raise ValueError("7R1C must have exactly one included PF00741 domain for unambiguous mapping")
    aligned = alignment_rows(domains[0], deposited, hmm_info["length"])
    by_position = {row["raw_position"]: row for row in aligned}
    if len(by_position) != len(aligned):
        raise ValueError("7R1C HMM alignment has duplicate residue positions")
    table = []
    for position, amino_acid in enumerate(deposited, 1):
        model = modeled_by_number.get(position)
        alignment = by_position.get(position)
        table.append({"deposited_position": position,
                      "pdb_residue_number": model.residue_number if model else "",
                      "pdb_insertion_code": model.insertion_code if model else "",
                      "residue": amino_acid, "modeled": model is not None,
                      "secondary_structure": model.secondary_structure if model else "",
                      "hmm_match_state": alignment["hmm_match_state"] if alignment else "",
                      "is_hmm_insertion": alignment["is_insertion"] if alignment else "",
                      "hmm_aligned": alignment is not None})
    args.out_dir.mkdir(parents=True)
    try:
        output = args.out_dir / "7r1c_residue_hmm_map.csv"
        write_csv(output, list(table[0]), table)
        summary = {"schema_version": "1.0", "run_utc": datetime.now(timezone.utc).isoformat(),
                   "purpose": "real 7R1C structural reference mapping, not functional labels",
                   "git": git_state(), "pdb_id": "7R1C", "chain": args.chain,
                   "command": ["python", "-m", "bioinformatics.m2_pf00741", *sys.argv[1:]],
                   "implementation_sha256": {"workflow": digest(Path(__file__)),
                                             "coordinates": digest(Path(__file__).with_name("m2_coordinates.py"))},
                   "pdb_source_url": PDB_URL, "pdb_sha256": digest(args.pdb),
                   "pdb_download_receipt_sha256": digest(args.pdb_receipt),
                   "hmm": hmm_info, "hmm_download_receipt_sha256": digest(args.hmm_receipt),
                   "software": {"python": sys.version.split()[0],
                                "pyhmmer": hmm_info["pyhmmer_version"]},
                   "n_deposited_residues": len(deposited), "n_modeled_residues": len(modeled),
                   "n_unmodeled_residues": len(deposited) - len(modeled),
                   "n_hmm_aligned_residues": len(aligned),
                   "n_modeled_hmm_aligned_residues": sum(row["modeled"] and row["hmm_aligned"] for row in table),
                   "n_hmm_match_residues": sum(row["hmm_match_state"] not in ("", None) for row in table),
                   "n_hmm_insertions": sum(row["is_hmm_insertion"] is True for row in table),
                   "n_distinct_hmm_match_states": len({row["hmm_match_state"] for row in table
                                                       if row["hmm_match_state"] not in ("", None)}),
                   "ga_hit_score_bits": float(hit.score),
                   "ga_domain_score_bits": float(domains[0].score),
                   "coordinate_convention": "1-based deposited position; PDB author number only for modeled residues; HMM insertion has no match state",
                   "coordinate_map_sha256": digest(output), "coordinate_map_filename": output.name,
                   "limitations": ["Only 65 of 88 deposited residues have experimental coordinates",
                                   "HMM match states are homology coordinates, not experimental evidence of function",
                                   "Secondary structure uses PDB HELIX/SHEET records for modeled residues"]}
        write_json(args.summary, summary)
        return summary
    except Exception:
        (args.out_dir / "FAILED.txt").write_text("Incomplete 7R1C mapping.\n", encoding="utf-8")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    structure = modes.add_parser("structure")
    structure.add_argument("--hmm", type=Path, required=True)
    structure.add_argument("--hmm-receipt", type=Path, required=True)
    structure.add_argument("--pdb", type=Path, required=True)
    structure.add_argument("--pdb-receipt", type=Path, required=True)
    structure.add_argument("--chain", default="N")
    structure.add_argument("--threads", type=int, default=1)
    structure.add_argument("--out-dir", type=Path, required=True)
    structure.add_argument("--summary", type=Path, required=True)
    scan = modes.add_parser("scan")
    scan.add_argument("--metadata", type=Path, required=True)
    scan.add_argument("--fasta", type=Path, required=True)
    scan.add_argument("--split-manifest", type=Path, required=True)
    scan.add_argument("--audit-summary", type=Path, required=True)
    scan.add_argument("--split-summary", type=Path, required=True)
    scan.add_argument("--hmm", type=Path, required=True)
    scan.add_argument("--hmm-receipt", type=Path, required=True)
    scan.add_argument("--out-dir", type=Path, required=True)
    scan.add_argument("--scan-summary", type=Path, required=True)
    scan.add_argument("--coordinate-summary", type=Path, required=True)
    scan.add_argument("--threads", type=int, default=1)
    scan.add_argument("--boundary-bits", type=float, default=2.0)
    args = parser.parse_args()
    try:
        result = run_structure(args) if args.mode == "structure" else run_scan(args)
    except (ValueError, OSError, RuntimeError, KeyError, UnicodeError) as error:
        parser.exit(2, f"error: {error}\n")
    print(json.dumps({"status": "completed", "mode": args.mode,
                      "result": str(args.summary if args.mode == "structure" else args.scan_summary)},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
