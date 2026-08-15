#!/usr/bin/env python3
"""Freeze a one-off golden genome fixture for regression tests (PH2-013).

Fetches GRCh38 chr22:20,000,000-20,100,000 (FASTA) and its Ensembl
Regulatory Build features from the Ensembl REST API and writes them under
``tests/data/golden_pipeline/``, so tests no longer depend on a live
network call (see ``tests/conftest.py``'s former ``real_fasta``/
``real_training_dir`` fixtures, which this script's fetch logic is ported
from verbatim).

This is a one-off/rerunnable authoring script, not part of the test suite
itself — it is meant to be run manually (or via
``scripts/regenerate_golden_fixtures.py``) when the frozen fixture needs to
be regenerated from scratch, which should be rare since the whole point of
freezing it is to stop depending on the live API.
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import requests

_ENSEMBL = "https://rest.ensembl.org"
_TIMEOUT = 60  # seconds


def fetch_sequence(chrom: str, start0: int, end: int, timeout: int = _TIMEOUT) -> str:
    """Fetch a chromosome slice from Ensembl REST as an uppercase string.

    ``start0`` is 0-based (BED/pyfaidx convention); ``end`` is exclusive.
    Ensembl's REST API itself uses 1-based inclusive coordinates, so the
    conversion happens here.
    """
    chrom_key = chrom.lstrip("chr")
    url = (
        f"{_ENSEMBL}/sequence/region/human/"
        f"{chrom_key}:{start0 + 1}..{end}"
        "?content-type=text/plain"
    )
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text.strip().upper()


def fetch_regulatory_bed(
    chrom: str, start0: int, end: int, timeout: int = _TIMEOUT
) -> list[tuple[str, int, int, str]]:
    """Fetch Ensembl Regulatory Build features overlapping a region.

    Returns rows as ``(chrom, start, end, label)`` tuples with coordinates
    converted to 0-based half-open and offset relative to ``start0`` (so
    position 0 in the row corresponds to ``start0`` in the real genome).
    """
    chrom_key = chrom.lstrip("chr")
    url = (
        f"{_ENSEMBL}/overlap/region/human/"
        f"{chrom_key}:{start0 + 1}..{end}"
        "?feature=regulatory"
    )
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=timeout)
    r.raise_for_status()

    rows: list[tuple[str, int, int, str]] = []
    for feat in r.json():
        label = feat.get("description", "")
        if not label:
            continue
        # Convert Ensembl 1-based -> 0-based, then subtract region start
        row_start = int(feat["start"]) - 1 - start0
        row_end = int(feat["end"]) - start0
        if row_start < 0 or row_end > (end - start0) or row_start >= row_end:
            continue
        rows.append((chrom, row_start, row_end, label))
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom", default="chr22")
    parser.add_argument("--start", type=int, default=20_000_000, help="0-based start")
    parser.add_argument("--end", type=int, default=20_100_000, help="exclusive end")
    parser.add_argument(
        "--out", default="tests/data/golden_pipeline", help="output directory"
    )
    args = parser.parse_args(argv)

    out_dir = Path(args.out)
    training_dir = out_dir / "training"
    out_dir.mkdir(parents=True, exist_ok=True)
    training_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"Fetching {args.chrom}:{args.start}-{args.end} sequence from Ensembl REST...",
        file=sys.stderr,
    )
    seq = fetch_sequence(args.chrom, args.start, args.end)
    wrapped = "\n".join(textwrap.wrap(seq, 60))
    genome_path = out_dir / "genome.fa"
    genome_path.write_text(f">{args.chrom}\n{wrapped}\n")
    print(f"Wrote {genome_path} ({len(seq)} bp)", file=sys.stderr)

    print("Fetching Ensembl Regulatory Build features...", file=sys.stderr)
    rows = fetch_regulatory_bed(args.chrom, args.start, args.end)
    print(f"Found {len(rows)} regulatory feature rows", file=sys.stderr)
    if len(rows) < 5:
        print(
            f"ERROR: only {len(rows)} regulatory features found — "
            "check network connectivity or widen the region.",
            file=sys.stderr,
        )
        return 1

    feature_types = list(dict.fromkeys(f for _, _, _, f in rows))
    (training_dir / "features.txt").write_text("\n".join(feature_types) + "\n")
    bed_lines = [f"{c}\t{s}\t{e}\t{f}" for c, s, e, f in rows]
    (training_dir / "training_regions.bed").write_text("\n".join(bed_lines) + "\n")
    print(
        f"Wrote {training_dir / 'training_regions.bed'} ({len(rows)} regions) "
        f"and {training_dir / 'features.txt'} ({len(feature_types)} features)",
        file=sys.stderr,
    )

    import pyfaidx

    pyfaidx.Faidx(str(genome_path))
    print(f"Wrote {genome_path}.fai", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
