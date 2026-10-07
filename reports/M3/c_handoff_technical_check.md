# C M3 coordinate handoff: assistant technical precheck

Status: **technical package checks passed; PR #25 remains unmerged and unreviewed by a non-author member**. This check does not sign for C or approve a biological candidate.

- PR: [#25](https://github.com/Lithium9767/MLP/pull/25), head `9d0f705c8e7eeb8374f221e990e38a0f772424aa`, author `fengbujue777` at check time.
- Actual local copy of `MLP_C_M3_candidate_coordinates_001.zip` was opened, CRC checked, and independently SHA-256 measured as `6b06159fe0a3c0ea01e2e1d34e41ff3d4fed6cbefb0cf1719f63596e34662473` (3,875,368 bytes).
- Every ZIP member listed in its `manifest.json` was streamed and matched against its declared byte length and SHA-256. There were zero mismatches.
- `candidate_windows.csv` SHA-256 is `d8b82d04cdee709507567a7e626e69f985592537a2580bc60485224d0fbec328`, identical to D's original discovery receipt. The coordinate table has 24,581 windows; the residue table has 737,430 rows, exactly 30 per window, with all 24,581 window IDs represented. Its manifest states `validation_used=false` and `candidate_freeze_status=exploratory_not_A_frozen`.
- From a detached checkout at PR head, `python -m pytest -q tests` produced `100 passed, 1 skipped, 4 subtests passed` with Python 3.13 and the local dependency set. The skipped test is not counted as verification.
- The five M2 JSON summaries appearing changed in PR #25 differ only by line endings (`git diff --ignore-space-at-eol` is empty). They should be removed from the PR diff before merge to reduce unrelated churn. The PR also edits the registry, which overlaps A/D's PR #22 and will need reconciliation.

**Public-record hygiene:** the tracked candidate mapping receipt in PR #25 contains a personal local filesystem/messaging path. Replace or redact that field in a new C-authored commit while preserving original execution provenance in a restricted handoff; do not invent a clean historical command. The current ZIP may carry the same path, so do not circulate it further without checking the intended audience and preparing a sanitized, hash-versioned replacement if public distribution is required. This issue is separate from coordinate correctness.

This precheck does not independently establish every HMM/7R1C coordinate value, C's authorship of the raw mapping run, a shared-storage access path, conservation scores, or a frozen candidate. C's own mapping claims and prior automated verification remain traceable in PR #25, while biological analysis and validation stay pending A's candidate decision.
