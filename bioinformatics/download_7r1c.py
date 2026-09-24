"""Download PDB 7R1C with an immutable-source SHA-256 receipt."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request


URL = "https://files.rcsb.org/download/7R1C.pdb"


def download(output_dir: Path) -> dict:
    pdb_path = output_dir / "7R1C.pdb"
    receipt_path = output_dir / "7R1C_download_receipt.json"
    if pdb_path.exists() or receipt_path.exists():
        raise FileExistsError("7R1C source or receipt already exists; use a new directory")
    request = urllib.request.Request(URL, headers={"User-Agent": "MLP-C-M2/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read()
    if not content.startswith(b"HEADER") or b"7R1C" not in content[:100] or b"\nATOM " not in content:
        raise ValueError("Downloaded response is not PDB 7R1C")
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_path.write_bytes(content)
    receipt = {"source_url": URL, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
               "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content),
               "path": str(pdb_path)}
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        receipt = download(args.output_dir)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
