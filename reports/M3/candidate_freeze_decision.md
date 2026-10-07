# M3 discovery candidate audit (pending A decision)

Status: **no candidate frozen; validation remains locked**. This is an assistant-executed technical audit, not A's approval or a member review.

The frozen discovery run `M3-D-DISCOVERY-001` used 1,202 primary discovery sequences and produced 24,581 overlapping 30-aa windows. Its exploratory cluster 3 covers 2,216 windows from 265 sequences, with 63.94% of its windows touching HMM state 25. That descriptive enrichment is insufficient to call the cluster a functional region.

The clean diagnostic run `M3-D-ROBUST-001` used code commit `8a0bf23715783de9fc1c5fe9fd3e6a5bb83644ee`, source window SHA-256 `d8b82d04cdee709507567a7e626e69f985592537a2580bc60485224d0fbec328`, and discovery data only. Its [receipt](runs/M3-D-ROBUST-001/run_receipt.json) records the exact command, environment and output hashes; [summary](runs/M3-D-ROBUST-001/summary.json), [length bins](runs/M3-D-ROBUST-001/length_bins.csv), and [sequence-level matched differences](runs/M3-D-ROBUST-001/sequence_deltas.csv) preserve the actual results.

| Discovery diagnostic | Observed result |
| --- | ---: |
| Median sequence length, cluster 3 participants / other sequences | 87 / 144 aa |
| Length alone distinguishing cluster 3 participation (shorter-length AUC) | 0.9688 |
| Sequences <=100 aa with cluster 3 windows | 245 / 298 (82.2%) |
| Sequences 131–160 aa with cluster 3 windows | 1 / 845 (0.12%) |
| Same-sequence position-matched windows, normalized-start tolerance 0.05 | 32 / 2,216 (1.4%) |
| Same-sequence position-matched windows, tolerance 0.10 | 294 / 2,216 (13.3%) |

At tolerance 0.10, the mean per-sequence difference in state-25 presence between matched candidate and noncandidate windows was +0.199 (descriptive bootstrap 95% interval 0.154–0.245 across 197 matched sequences). Most candidate windows had no suitable same-sequence control; the matched subset is selected and does not establish a cluster-wide effect. Windows overlap, state 25 was inspected after selecting the cluster, and these are **not** functional p-values or validation results. The length AUC measures a technical confound, not predictive biological accuracy.

Other clusters do not currently rescue the freeze: cluster 0 lacks usable HMM mapping; cluster 1 did not exceed the prior matched-random concentration diagnostic; cluster 2 is diffuse and has no majority common HMM region. Parameter sensitivity also varies substantially across window sizes. See [the existing discovery report](m3_summary.md) and [candidate review](candidate_review.csv).

**Decision needed from A:** review a revised discovery-only candidate definition or explicitly approve a narrowly labeled exploratory cluster despite these limitations. No approval is recorded here. D must not run `scripts/m3_validate_frozen.py` until a concrete proposal and genuine A approval satisfy its hash-checked gate. Do not use validation to select a better cluster or change its thresholds.

Independent work may proceed on C's coordinate checks, B's provenance, and E's unselected-method tests. Biological conservation, mutation and structure conclusions remain unverified until their respective real evidence is produced.
