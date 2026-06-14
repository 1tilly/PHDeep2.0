#!/usr/bin/env python3
"""
Download Roadmap Epigenomics consolidated narrowPeak files.

Downloads histone modification and DNase peaks for all (or selected)
Roadmap epigenomes from the consolidated peak repository at WUSTL.

Usage
-----
    # Download H3K27ac and H3K4me3 peaks for all epigenomes
    python scripts/download_roadmap_epigenomics.py \\
        --marks H3K27ac H3K4me3 DNase \\
        --output data/roadmap_beds/ \\
        --workers 8

    # Download only specific epigenomes (E062 = Primary mononuclear cells
    #  from peripheral blood; E116 = GM12878)
    python scripts/download_roadmap_epigenomics.py \\
        --marks H3K27ac \\
        --epigenomes E062 E116 \\
        --output data/roadmap_beds/

Marks
-----
Histone: H3K4me3, H3K4me1, H3K36me3, H3K27me3, H3K9me3, H3K27ac
DNase:   DNase  (uses a slightly different filename pattern)

Data source
-----------
https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
Reference: Roadmap Epigenomics Consortium, Nature 2015.
"""
from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

BASE_URL = "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak"

# All 127 Roadmap reference epigenomes (E001–E129, with gaps)
ALL_EPIGENOMES = [
    f"E{n:03d}" for n in range(1, 130)
    if n not in {60, 64}  # E060 and E064 were withdrawn
]

HISTONE_MARKS = {"H3K4me3", "H3K4me1", "H3K36me3", "H3K27me3", "H3K9me3", "H3K27ac"}
ALL_MARKS = HISTONE_MARKS | {"DNase"}


def _build_url(epigenome: str, mark: str) -> str:
    if mark == "DNase":
        filename = f"{epigenome}-DNase.macs2.narrowPeak.gz"
    else:
        filename = f"{epigenome}-{mark}.narrowPeak.gz"
    return f"{BASE_URL}/{filename}"


def _download_file(url: str, dest: Path) -> tuple[str, bool, str]:
    if dest.exists() and dest.stat().st_size > 0:
        return url, True, "already exists"
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(url, stream=True, timeout=120) as r:
            if r.status_code == 404:
                return url, False, "404 not found (epigenome/mark combination may not exist)"
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
        return url, True, f"downloaded {dest.stat().st_size // 1024} KB"
    except Exception as exc:
        return url, False, str(exc)


def main():
    parser = argparse.ArgumentParser(description="Download Roadmap Epigenomics narrowPeak files.")
    parser.add_argument(
        "--marks",
        nargs="+",
        default=["H3K27ac", "H3K4me3", "DNase"],
        help=f"Histone marks / assays to download. Available: {sorted(ALL_MARKS)}",
    )
    parser.add_argument(
        "--epigenomes",
        nargs="*",
        default=None,
        help="Subset of epigenome IDs (e.g. E062 E116). Default: all 127.",
    )
    parser.add_argument("--output", required=True, help="Output directory.")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    invalid_marks = set(args.marks) - ALL_MARKS
    if invalid_marks:
        logger.error("Unknown marks: %s. Valid: %s", invalid_marks, sorted(ALL_MARKS))
        sys.exit(1)

    epigenomes = args.epigenomes or ALL_EPIGENOMES
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    tasks: list[tuple[str, Path]] = []
    for epigenome in epigenomes:
        for mark in args.marks:
            url = _build_url(epigenome, mark)
            dest = out_dir / f"{epigenome}-{mark}.narrowPeak.gz"
            tasks.append((url, dest))

    logger.info(
        "%d files to download (%d epigenomes × %d marks)",
        len(tasks), len(epigenomes), len(args.marks),
    )

    if args.dry_run:
        for url, dest in tasks:
            print(f"WOULD DOWNLOAD: {url} → {dest}")
        return

    n_ok = n_fail = n_miss = 0
    fail_log = out_dir / "download_errors.txt"
    with fail_log.open("w") as err_fh:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(_download_file, url, dest): url for url, dest in tasks}
            for i, future in enumerate(as_completed(futures), 1):
                url, success, msg = future.result()
                if success:
                    n_ok += 1
                elif "404" in msg:
                    n_miss += 1
                    logger.debug("MISSING: %s", url)
                else:
                    n_fail += 1
                    err_fh.write(f"FAIL\t{url}\t{msg}\n")
                    logger.warning("FAIL [%d/%d]: %s", i, len(tasks), msg)
                if i % 100 == 0:
                    logger.info(
                        "Progress: %d/%d (ok=%d, missing=%d, fail=%d)",
                        i, len(tasks), n_ok, n_miss, n_fail,
                    )

    logger.info(
        "Done. Downloaded: %d  Not found (404): %d  Failed: %d",
        n_ok, n_miss, n_fail,
    )
    if n_fail > 0:
        logger.info("Error log: %s", fail_log)


if __name__ == "__main__":
    main()
