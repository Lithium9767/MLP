#!/usr/bin/env python3
"""Download Pfam HMMs with a provenance receipt."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import urllib.request
from pathlib import Path


DEFAULT_ACCESSIONS = ("PF00741",)
API = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}?annotation=hmm"


def download_hmm(accession: str, output_path: Path, *, timeout: int = 120) -> dict[str, object]:
    if not accession.startswith("PF") or not accession[2:].isdigit():
        raise ValueError(f"Invalid Pfam accession {accession!r}")
    url = API.format(accession=accession)
    request = urllib.request.Request(url, headers={"User-Agent": "MLP-GvpA/2.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        data = response.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    if not data.startswith(b"HMMER"):
        raise RuntimeError(f"{accession} response is not an HMM")
    headers = {}
    for line in data.splitlines():
        if line.startswith(b"HMM "):
            break
        parts = line.decode("ascii").split(maxsplit=1)
        if len(parts) == 2:
            headers[parts[0]] = parts[1].strip()
    model_accession = headers.get("ACC", "")
    if model_accession.split(".")[0] != accession or not model_accession[len(accession):].startswith("."):
        raise ValueError(f"Downloaded HMM accession {model_accession!r} does not match {accession}")
    if not headers.get("GA") or not headers.get("LENG"):
        raise ValueError("Downloaded HMM lacks length or gathering thresholds")
    if output_path.exists():
        raise FileExistsError(f"Will not overwrite existing HMM: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(data)
    return {
        "accession": accession,
        "model_accession": model_accession,
        "model_date": headers.get("DATE", ""),
        "gathering_thresholds": headers["GA"],
        "url": url,
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "size": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "output": str(output_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("data/raw/pfam"))
    parser.add_argument("--accessions", nargs="+", default=list(DEFAULT_ACCESSIONS))
    args = parser.parse_args()
    receipts = [
        download_hmm(accession, args.output_dir / f"{accession}.hmm")
        for accession in args.accessions
    ]
    receipt_path = args.output_dir / "download_receipt.json"
    receipt_path.write_text(
        json.dumps(receipts, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipts, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
