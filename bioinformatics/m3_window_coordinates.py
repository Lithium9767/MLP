"""M3 window-to-PF00741/7R1C coordinates from *existing* M2 tables.

This is a coordinate reader, not a new HMM scan or a functional classifier.
Natural-sequence and window coordinates are 1-based and inclusive. A PF00741
match state is a homology coordinate; a 7R1C author residue number exists only
where the deposited reference residue was experimentally modeled.

Call ``load_coordinate_index`` after verifying the M2 input/output receipts,
or construct an index from already checked rows with ``from_rows``. The
``min_match_fraction`` argument only flags technical coverage; candidate
selection and its threshold must be frozen separately by M3 A.
"""

from __future__ import annotations

from collections import defaultdict
import csv
from pathlib import Path

from .m2_pf00741 import CANONICAL, sequence_digest


RECORD_FIELDS = frozenset({"internal_id", "sequence", "sequence_sha256", "split",
                           "analysis_cohort", "primary_analysis_eligible"})
STATUS_FIELDS = frozenset({"internal_id", "sequence_sha256", "split", "analysis_cohort",
                           "is_partial_or_conflict", "status", "n_ga_domains",
                           "n_mapped_residues", "n_match_states"})
COORDINATE_FIELDS = frozenset({"internal_id", "sequence_sha256", "split", "analysis_cohort",
                               "domain_number", "raw_position", "residue",
                               "hmm_match_state", "is_insertion"})
STRUCTURE_FIELDS = frozenset({"deposited_position", "pdb_residue_number",
                              "pdb_insertion_code", "residue", "modeled",
                              "secondary_structure", "hmm_match_state",
                              "is_hmm_insertion", "hmm_aligned"})
HIT_STATUSES = frozenset({"accepted", "multi_hit"})
ALL_STATUSES = HIT_STATUSES | {"boundary", "weak_hit", "no_hit", "failed"}


def _bool(value: object, name: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in {"true", "false"}:
        return value.lower() == "true"
    raise ValueError(f"{name} must be a true/false value")


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")
    try:
        result = int(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"{name} must be a positive integer") from error
    if result < 1 or str(value).strip() != str(result):
        raise ValueError(f"{name} must be a positive integer")
    return result


def _nonnegative_int(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        result = int(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"{name} must be a nonnegative integer") from error
    if result < 0 or str(value).strip() != str(result):
        raise ValueError(f"{name} must be a nonnegative integer")
    return result


def _optional_state(value: object, hmm_length: int) -> int | None:
    if value is None or value == "":
        return None
    state = _positive_int(value, "hmm_match_state")
    if state > hmm_length:
        raise ValueError("HMM match state exceeds model length")
    return state


def _csv_rows(path: Path, required: frozenset[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)) or not required <= set(fields):
            raise ValueError(f"Missing/duplicate CSV headers in {path}")
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"Malformed CSV row in {path}")
    return rows


class CoordinateIndex:
    """Validated M2 natural-residue and 7R1C reference coordinate index."""

    @classmethod
    def from_rows(cls, records: list[dict], status_rows: list[dict],
                  coordinate_rows: list[dict], structure_rows: list[dict],
                  *, hmm_length: int = 39) -> CoordinateIndex:
        """Validate one-to-one record/status and all aligned rows before mapping.

        ``records`` are the M2/B frozen records, including full sequences. The
        existing M2 scan's GA-qualified domains are the only coordinate source.
        ``structure_rows`` are the M2 7R1C deposited-residue table, including
        unmodeled residues. This method does not itself verify file SHA-256.
        """
        self = cls()
        self.hmm_length = _positive_int(hmm_length, "hmm_length")
        self.records = {}
        for row in records:
            if not RECORD_FIELDS <= row.keys():
                raise ValueError("Frozen input record lacks required fields")
            identifier, sequence = row["internal_id"], row["sequence"]
            if not identifier or identifier in self.records or not sequence or set(sequence) - CANONICAL:
                raise ValueError("Duplicate ID or invalid natural sequence")
            if sequence_digest(sequence) != row["sequence_sha256"]:
                raise ValueError(f"Sequence SHA-256 mismatch for {identifier}")
            if row["split"] not in {"discovery", "validation"} or not row["analysis_cohort"]:
                raise ValueError("Invalid split/cohort in frozen input")
            self.records[identifier] = {**row, "primary_analysis_eligible":
                                        _bool(row["primary_analysis_eligible"],
                                              "primary_analysis_eligible")}
        if not self.records:
            raise ValueError("No frozen input records")

        self.statuses = {}
        for row in status_rows:
            if not STATUS_FIELDS <= row.keys():
                raise ValueError("M2 status row lacks required fields")
            identifier = row["internal_id"]
            if identifier not in self.records or identifier in self.statuses:
                raise ValueError("Unknown or duplicate status ID")
            record = self.records[identifier]
            if any(row[key] != record[key] for key in ("sequence_sha256", "split", "analysis_cohort")):
                raise ValueError(f"M2 status identity differs from input for {identifier}")
            if row["status"] not in ALL_STATUSES:
                raise ValueError(f"Unknown M2 HMM status for {identifier}")
            if _bool(row["is_partial_or_conflict"], "is_partial_or_conflict") == record["primary_analysis_eligible"]:
                raise ValueError(f"M2 sensitivity/primary flag disagrees for {identifier}")
            self.statuses[identifier] = {**row,
                                         "n_ga_domains": _nonnegative_int(row["n_ga_domains"], "n_ga_domains"),
                                         "n_mapped_residues": _nonnegative_int(row["n_mapped_residues"], "n_mapped_residues"),
                                         "n_match_states": _nonnegative_int(row["n_match_states"], "n_match_states")}
        if set(self.statuses) != set(self.records):
            raise ValueError("M2 status table does not cover every frozen input ID")

        self.coordinates: dict[str, dict[int, list[dict]]] = defaultdict(lambda: defaultdict(list))
        domains: dict[str, dict[int, list[dict]]] = defaultdict(lambda: defaultdict(list))
        used_keys = set()
        for row in coordinate_rows:
            if not COORDINATE_FIELDS <= row.keys():
                raise ValueError("M2 coordinate row lacks required fields")
            identifier = row["internal_id"]
            if identifier not in self.records or self.statuses[identifier]["status"] not in HIT_STATUSES:
                raise ValueError("Coordinate row refers to unknown or non-GA-qualified input")
            record = self.records[identifier]
            if any(row[key] != record[key] for key in ("sequence_sha256", "split", "analysis_cohort")):
                raise ValueError(f"M2 coordinate identity differs from input for {identifier}")
            domain = _positive_int(row["domain_number"], "domain_number")
            position = _positive_int(row["raw_position"], "raw_position")
            key = (identifier, domain, position)
            if key in used_keys or position > len(record["sequence"]):
                raise ValueError("Repeated domain/position or out-of-range natural coordinate")
            used_keys.add(key)
            if row["residue"] != record["sequence"][position - 1]:
                raise ValueError(f"M2 coordinate residue differs from natural sequence for {identifier}")
            state = _optional_state(row["hmm_match_state"], self.hmm_length)
            insertion = _bool(row["is_insertion"], "is_insertion")
            if insertion != (state is None):
                raise ValueError("Insertion must have null HMM match state, and vice versa")
            item = {"domain_number": domain, "raw_position": position,
                    "hmm_match_state": state, "is_insertion": insertion}
            self.coordinates[identifier][position].append(item)
            domains[identifier][domain].append(item)
        for identifier, record in self.records.items():
            status = self.statuses[identifier]
            positions = self.coordinates[identifier]
            observed_states = {item["hmm_match_state"] for items in positions.values()
                               for item in items if item["hmm_match_state"] is not None}
            if (status["n_mapped_residues"] != len(positions)
                    or status["n_match_states"] != len(observed_states)
                    or (status["status"] in HIT_STATUSES
                        and status["n_ga_domains"] != len(domains[identifier]))):
                raise ValueError(f"M2 coordinate/status counts disagree for {identifier}")
            if status["status"] in HIT_STATUSES and not positions:
                raise ValueError(f"GA-qualified hit has no coordinate rows for {identifier}")
            if status["status"] == "accepted" and len(domains[identifier]) != 1:
                raise ValueError(f"Single accepted hit must have one domain for {identifier}")
            if status["status"] == "multi_hit" and len(domains[identifier]) < 2:
                raise ValueError(f"Multi-hit status must have multiple domains for {identifier}")
            if status["status"] in HIT_STATUSES and set(domains[identifier]) != set(
                    range(1, status["n_ga_domains"] + 1)):
                raise ValueError(f"M2 domain numbering has gaps for {identifier}")
            for domain, items in domains[identifier].items():
                items.sort(key=lambda item: item["raw_position"])
                last_state = 0
                for previous, current in zip(items, items[1:]):
                    if current["raw_position"] != previous["raw_position"] + 1:
                        raise ValueError(f"Gap in M2 alignment rows for {identifier} domain {domain}")
                for item in items:
                    state = item["hmm_match_state"]
                    if state is not None:
                        if state <= last_state:
                            raise ValueError(f"Non-increasing HMM states for {identifier} domain {domain}")
                        last_state = state

        self.reference_by_state = {}
        seen_positions = set()
        for row in structure_rows:
            if not STRUCTURE_FIELDS <= row.keys():
                raise ValueError("7R1C coordinate row lacks required fields")
            position = _positive_int(row["deposited_position"], "deposited_position")
            if position in seen_positions:
                raise ValueError("Repeated 7R1C deposited position")
            seen_positions.add(position)
            modeled = _bool(row["modeled"], "modeled")
            aligned = _bool(row["hmm_aligned"], "hmm_aligned")
            state = _optional_state(row["hmm_match_state"], self.hmm_length)
            insertion = (None if row["is_hmm_insertion"] in (None, "") else
                         _bool(row["is_hmm_insertion"], "is_hmm_insertion"))
            if (not aligned and (state is not None or insertion is not None)) or (aligned and insertion != (state is None)):
                raise ValueError("7R1C HMM alignment/insertion flags disagree")
            if modeled:
                pdb_number = _positive_int(row["pdb_residue_number"], "pdb_residue_number")
                if row["secondary_structure"] not in {"H", "E", "C"}:
                    raise ValueError("Modeled 7R1C residue lacks secondary structure")
            else:
                pdb_number = None
                if row["pdb_residue_number"] not in (None, "") or row["secondary_structure"] not in (None, ""):
                    raise ValueError("Unmodeled 7R1C residue must not have PDB coordinates")
            if not isinstance(row["residue"], str) or len(row["residue"]) != 1 or row["residue"] not in CANONICAL:
                raise ValueError("Invalid 7R1C deposited residue")
            if state is not None:
                if state in self.reference_by_state:
                    raise ValueError("One HMM state maps to multiple 7R1C residues")
                self.reference_by_state[state] = {
                    "deposited_position": position,
                    "pdb_residue_number": pdb_number,
                    "pdb_insertion_code": row["pdb_insertion_code"] or None,
                    "modeled": modeled,
                    "secondary_structure": row["secondary_structure"] or None,
                    "reference_residue": row["residue"],
                }
        if not seen_positions or seen_positions != set(range(1, max(seen_positions) + 1)):
            raise ValueError("7R1C deposited positions must be complete and contiguous")
        return self

    def map_window(self, internal_id: str, start: int, end: int, *,
                   min_match_fraction: float = 0.8) -> dict:
        """Map an inclusive natural-sequence window without filling missing states.

        A low-coverage flag is a technical annotation, not a discovery decision.
        M3 A must freeze any acceptance threshold before validation is opened.
        """
        if internal_id not in self.records:
            raise KeyError(f"Unknown internal_id: {internal_id}")
        if isinstance(start, bool) or isinstance(end, bool) or not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("Window start/end must be integer 1-based positions")
        record = self.records[internal_id]
        if start < 1 or end < start or end > len(record["sequence"]):
            raise ValueError("Window must be within the full natural sequence, 1-based inclusive")
        if not isinstance(min_match_fraction, (int, float)) or isinstance(min_match_fraction, bool) or not 0 <= min_match_fraction <= 1:
            raise ValueError("min_match_fraction must lie in [0, 1]")
        residues = []
        counts = {"match": 0, "insertion": 0, "unaligned": 0, "ambiguous_multi_domain": 0}
        states = set()
        for position in range(start, end + 1):
            options = sorted(self.coordinates[internal_id].get(position, []),
                             key=lambda item: item["domain_number"])
            if len(options) > 1:
                kind, state, domain = "ambiguous_multi_domain", None, None
            elif options:
                state, domain = options[0]["hmm_match_state"], options[0]["domain_number"]
                kind = "insertion" if state is None else "match"
            else:
                kind, state, domain = "unaligned", None, None
            counts[kind] += 1
            reference = self.reference_by_state.get(state) if state is not None else None
            reference_status = ("ambiguous_multi_domain" if kind == "ambiguous_multi_domain"
                                else "no_match_state" if state is None
                                else "not_mapped_in_7r1c" if reference is None
                                else "reference_modeled" if reference["modeled"]
                                else "reference_unmodeled")
            if state is not None:
                states.add(state)
            residues.append({"raw_position": position, "residue": record["sequence"][position - 1],
                             "mapping_status": kind, "domain_number": domain,
                             "hmm_match_state": state, "coordinate_options": options,
                             "reference_status": reference_status,
                             "reference_7r1c": reference})
        length = end - start + 1
        match_fraction = counts["match"] / length
        domain_numbers = sorted({item["domain_number"] for residue in residues
                                 for item in residue["coordinate_options"]})
        if counts["ambiguous_multi_domain"]:
            coverage_status = "ambiguous_multi_domain"
        elif len(domain_numbers) > 1:
            coverage_status = "multi_domain_window"
        elif counts["match"] + counts["insertion"] == 0:
            coverage_status = "unmapped"
        elif match_fraction < min_match_fraction:
            coverage_status = "low_match_coverage"
        else:
            coverage_status = "mapped"

        # Infer only deletions *internal* to a continuously represented M2
        # alignment segment. Do not invent HMM states at unaligned window edges.
        deleted_states = set()
        by_domain: dict[int, list[dict]] = defaultdict(list)
        for residue in residues:
            if residue["domain_number"] is not None:
                by_domain[residue["domain_number"]].append(residue)
        for domain_residues in by_domain.values():
            previous_state = None
            previous_position = None
            for residue in domain_residues:
                state = residue["hmm_match_state"]
                if state is not None:
                    if previous_state is not None and all(
                            item["domain_number"] == residue["domain_number"]
                            for item in residues[previous_position - start + 1:residue["raw_position"] - start]):
                        deleted_states.update(range(previous_state + 1, state))
                    previous_state, previous_position = state, residue["raw_position"]
        return {
            "internal_id": internal_id, "sequence_sha256": record["sequence_sha256"],
            "split": record["split"], "analysis_cohort": record["analysis_cohort"],
            "primary_analysis_eligible": record["primary_analysis_eligible"],
            "m2_hmm_status": self.statuses[internal_id]["status"],
            "start_1based": start, "end_1based_inclusive": end, "window_length": length,
            "min_match_fraction": float(min_match_fraction),
            "match_fraction": match_fraction,
            "aligned_fraction": (counts["match"] + counts["insertion"]) / length,
            "coverage_status": coverage_status,
            "n_match_residues": counts["match"], "n_insertion_residues": counts["insertion"],
            "n_unaligned_residues": counts["unaligned"],
            "n_ambiguous_residues": counts["ambiguous_multi_domain"],
            "domain_numbers": domain_numbers, "n_domains": len(domain_numbers),
            "hmm_match_states": sorted(states),
            "hmm_deleted_states_internal": sorted(deleted_states),
            "residues": residues,
        }


def load_coordinate_index(records: list[dict], status_path: Path, coordinate_path: Path,
                          structure_path: Path, *, hmm_length: int = 39) -> CoordinateIndex:
    """Read M2 CSV tables after caller verifies frozen receipts and SHA-256.

    No HMM rescanning, conservation estimation, or validation selection occurs.
    """
    return CoordinateIndex.from_rows(
        records,
        _csv_rows(status_path, STATUS_FIELDS),
        _csv_rows(coordinate_path, COORDINATE_FIELDS),
        _csv_rows(structure_path, STRUCTURE_FIELDS),
        hmm_length=hmm_length,
    )
