# M3 Completion Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Complete the remaining M3 engineering work in dependency order while preserving real A/member/teacher approvals and the discovery/validation boundary.

**Architecture:** Audit discovery candidates with sequence-level and length/position-matched diagnostics before any freeze. Bind any A-approved candidate to immutable proposal/receipt hashes, then run the existing validation gate once without refitting. Integrate C/E evidence and B provenance through separate reviewable PRs; do not treat assistant checks as member reviews.

**Tech Stack:** Python 3.13, CSV/JSON, NumPy/scikit-learn/hdbscan, ESM-2, GitHub PRs.

---

### Task 1: Audit discovery robustness (A/D)

**Files:** Add `scripts/m3_candidate_robustness.py`, `tests/test_m3_candidate_robustness.py`; archive a run receipt and small statistics under `reports/M3/runs/`.

1. Add tests for per-sequence aggregation, same-sequence null sampling, and length-bin confounding checks. Run tests and observe failure before implementation.
2. Implement discovery-only analysis of `candidate_windows.csv` and frozen metadata; reject validation records and hash mismatches.
3. Compare candidate coverage at the sequence level, matched noncandidate windows, normalized position, and length distribution. Preserve negative findings.
4. Commit code on PR #22; run from a clean commit with a new run ID; archive receipt and output hashes in a later commit.
5. Review whether any specific cluster is sufficiently stable to offer A a concrete freeze choice. Do not infer A approval from the analysis.

### Task 2: Freeze gate and validation (A then D)

**Files:** `reports/M3/runs/M3-D-DISCOVERY-001/candidate_freeze_proposal.json`, A's approval record, `scripts/m3_validate_frozen.py`, a new validation run receipt.

1. Present A with the exact candidate cluster(s), parameters, proposal hash, receipt hash, and limitations. If none is defensible, keep validation locked and report that outcome.
2. Only after a genuine A decision, bind `approved_by`, `approval_evidence`, selected clusters, and current hashes to a reviewable record.
3. Run validation exactly once from clean code and fixed inputs. Do not refit PCA/HDBSCAN or select thresholds from validation.
4. Record actual SHA, commands, output hashes, success/failure, and comparison with discovery.

### Task 3: C/E evidence (after freeze; independent PRs)

1. Review PR #25 against D's original window table and M2 HMM/7R1C coordinates; keep C's own authorship and review history.
2. Align C's coordinate interface with the frozen candidate. Add conservation/structure evidence only where actually computed.
3. Implement E's residue perturbation and matched random baseline with token/residue/HMM coordinate tests; keep claims exploratory.

### Task 4: B provenance and repository integration

1. Resolve Issue #24 only with B's actual generation method, run identity, shared-storage URI/access and measured hash. Assistant precheck #23 remains labeled as such.
2. Reconcile overlapping #17/#22 ESM-2 changes, review #7/#17/#22/#23/#25 with genuine non-author members, then merge per `CONTRIBUTING.md`.
3. Run the complete test suite and final M3 validator/checks on the proposed final commit; update M3 report, registry and midterm materials.

**Stop conditions:** Missing human A approval, member attestation, teacher confirmation or non-author Review remains explicitly pending. Do not invent them, publish a member Review, or mark M3 complete while they are absent.
