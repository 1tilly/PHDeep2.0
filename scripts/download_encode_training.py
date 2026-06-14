#!/usr/bin/env python3
"""
Download ENCODE BED files for training.

Workflow
--------
1. Go to https://www.encodeproject.org/search/
2. Filter: Assay = DNase-seq, TF ChIP-seq, Histone ChIP-seq
           Output type = peaks
           File format = bed narrowPeak
           Organism = Homo sapiens
           Assembly = GRCh38  (or hg19)
3. Click "Download → Metadata only" → save as encode_metadata.tsv
4. Run this script:

    python scripts/download_encode_training.py \\
        --metadata data/encode_metadata.tsv \\
        --output   data/encode_beds/ \\
        --assembly GRCh38 \\
        --workers  8

The script downloads only the files referenced in the metadata TSV,
skips already-downloaded files (resume-safe), and logs failures.
"""
from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

ASSAYS = {"DNase-seq", "TF ChIP-seq", "Histone ChIP-seq", "ATAC-seq"}
FILE_FORMATS = {"bed narrowPeak"}


def _download_file(url: str, dest: Path) -> tuple[str, bool, str]:
    """Download a single file. Returns (url, success, message)."""
    if dest.exists() and dest.stat().st_size > 0:
        return url, True, "already exists"
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
        return url, True, f"downloaded {dest.stat().st_size // 1024} KB"
    except Exception as exc:
        return url, False, str(exc)


def main():
    parser = argparse.ArgumentParser(description="Download ENCODE BED training files.")
    parser.add_argument("--metadata", required=True, help="Path to ENCODE metadata TSV.")
    parser.add_argument("--output", required=True, help="Output directory for BED files.")
    parser.add_argument("--assembly", default="GRCh38", choices=["GRCh38", "hg19"])
    parser.add_argument("--workers", type=int, default=4, help="Parallel download workers.")
    parser.add_argument(
        "--assays",
        nargs="+",
        default=sorted(ASSAYS),
        help="Assay types to include (space-separated).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print files that would be downloaded without fetching them.",
    )
    args = parser.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Reading metadata from %s", args.metadata)
    meta = pd.read_csv(args.metadata, sep="\t", low_memory=False)

    mask = (
        (meta["Biosample organism"] == "Homo sapiens")
        & (meta["Assay"].isin(set(args.assays)))
        & (meta["File format"].isin(FILE_FORMATS))
        & (meta["File assembly"] == args.assembly)
        & (meta["Output type"] == "peaks")
    )
    filtered = meta[mask].dropna(subset=["File download URL", "File accession"])

    logger.info(
        "%d files selected (assembly=%s, assays=%s)",
        len(filtered), args.assembly, args.assays,
    )

    if len(filtered) == 0:
        logger.error("No files matched the filter criteria. Check your metadata TSV and filters.")
        sys.exit(1)

    tasks: list[tuple[str, Path]] = []
    for _, row in filtered.iterrows():
        url = row["File download URL"]
        accession = row["File accession"]
        file_type = str(row.get("File type", "bed.gz")).replace(" ", ".")
        dest = out_dir / f"{accession}.{file_type}.gz"
        tasks.append((url, dest))

    if args.dry_run:
        for url, dest in tasks:
            print(f"WOULD DOWNLOAD: {url} → {dest}")
        return

    n_ok = n_fail = 0
    fail_log = out_dir / "download_errors.txt"
    with fail_log.open("w") as err_fh:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(_download_file, url, dest): (url, dest) for url, dest in tasks}
            for i, future in enumerate(as_completed(futures), 1):
                url, success, msg = future.result()
                if success:
                    n_ok += 1
                else:
                    n_fail += 1
                    err_fh.write(f"FAIL\t{url}\t{msg}\n")
                    logger.warning("FAIL [%d/%d]: %s — %s", i, len(tasks), url, msg)
                if i % 50 == 0:
                    logger.info("Progress: %d/%d (ok=%d, fail=%d)", i, len(tasks), n_ok, n_fail)

    logger.info("Done. Downloaded: %d  Failed: %d  Error log: %s", n_ok, n_fail, fail_log)


if __name__ == "__main__":
    main()
