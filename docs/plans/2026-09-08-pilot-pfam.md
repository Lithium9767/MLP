# Pfam pilot implementation plan

Goal: reproducible 30/15/15 recognition sampling and local-only InterProScan preparation, followed by conservative parsing of real results.

Architecture: a standard-library `scripts/pilot_pfam.py` command with select/run/parse/summarize subcommands shares hashes, manifests and rules. This replaces four thin wrapper files to avoid duplicated validation. Changes stay on a local feature branch; no remote publishing or sequence upload.

1. Inspect actual JSON schema and native/WSL tools. Explicitly select GvpA input and separate other-family background inputs.
2. Implement deterministic, provenance-preserving selection with QC and approximate near-duplicate screening, and record limits of that screen.
3. Provide local runner that checks installed help/version, disables remote precalculated lookup and records receipts/logs. When absent, emit scan_not_run and zero submitted counts.
4. Parse XML with complete sequence identity checks and successful run receipt. TSV omission alone can never establish absence. Freeze PF00741 rule before execution.
5. Test schema failures, duplicates/hashes/determinism, false/missing signatures, incomplete runs, boundary hits and accounting. Test fixtures are explicitly synthetic and never used as scan outputs.
6. Generate selected inputs and pending-status report, document Linux execution commands, stop at preparation if no local scanner exists.

Environment inspection: Windows PATH has no InterProScan/HMMER; Ubuntu WSL is present, but PATH and bounded /opt,/usr/local,/home searches found no interproscan.sh or Pfam-A.hmm. No remote submission authorized. Input remains gv.zip::gv/recognition/data/gvp/gvpa/GvpA_sequences.json.
