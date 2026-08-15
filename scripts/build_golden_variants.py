#!/usr/bin/env python3
"""Build the golden pipeline's variants.tsv fixture (PH2-013).

Generates 18 variants (14 SNV, 2 insertion, 2 deletion) at the midpoints of
training regions from the frozen ``tests/data/golden_pipeline/genome.fa`` /
``training/training_regions.bed`` fixture (see
``scripts/fetch_golden_genome_fixture.py``), each with a 500bp flank clear
of both fixture edges.

Reference alleles are read directly from the frozen FASTA via pyfaidx (never
hand-typed) so the pipeline's own ref-mismatch sanity-check warning never
fires on this fixture. Coordinates are 1-based VCF POS, LOCAL to
``genome.fa`` — position 1 is NOT a real chr22 coordinate.

Run with no arguments to (re)write ``tests/data/golden_pipeline/variants.tsv``.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pyfaidx

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "golden_pipeline"

_CHROM = "chr22"

# (0-based position, variant type, gene_symbol group). Groups: GENE_A (6),
# GENE_B (5), GENE_C (4), GENE_SINGLE (1, tests min_variants filtering),
# and 2 rows with an empty group (tests "ungrouped variants retained").
_PLAN: list[tuple[int, str, str]] = [
    (780, "SNV", "GENE_A"),
    (1434, "SNV", "GENE_A"),
    (2256, "SNV", "GENE_A"),
    (2776, "DEL", "GENE_A"),
    (10695, "SNV", "GENE_A"),
    (13491, "SNV", "GENE_A"),
    (14715, "SNV", "GENE_B"),
    (16352, "SNV", "GENE_B"),
    (16918, "INS", "GENE_B"),
    (17493, "SNV", "GENE_B"),
    (20863, "SNV", "GENE_B"),
    (21767, "SNV", "GENE_C"),
    (22912, "DEL", "GENE_C"),
    (28150, "SNV", "GENE_C"),
    (32325, "SNV", "GENE_C"),
    (33458, "SNV", "GENE_SINGLE"),
    (62799, "INS", ""),
    (84972, "SNV", ""),
]

_SNV_ROTATE = {"A": "C", "C": "G", "G": "T", "T": "A"}

_HEADER = ["chromosome", "start", "reference", "alternate", "gene_symbol"]


def build_variant_rows(genome_path: Path) -> list[dict[str, str]]:
    genome = pyfaidx.Fasta(str(genome_path), as_raw=True)
    seq = str(genome[_CHROM][:])

    rows: list[dict[str, str]] = []
    for pos0, vtype, group in _PLAN:
        pos1 = pos0 + 1  # VCF POS is 1-based
        if vtype == "SNV":
            ref = seq[pos0]
            alt = _SNV_ROTATE[ref]
        elif vtype == "INS":
            ref = seq[pos0]
            alt = ref + "TT"
        elif vtype == "DEL":
            ref = seq[pos0 : pos0 + 3]
            alt = ref[0]
        else:
            raise ValueError(f"Unknown variant type {vtype!r}")

        assert alt != ref, f"alt must differ from ref at {pos1}"
        rows.append(
            {
                "chromosome": _CHROM,
                "start": str(pos1),
                "reference": ref,
                "alternate": alt,
                "gene_symbol": group,
            }
        )

    keys = [(r["chromosome"], r["start"], r["reference"], r["alternate"]) for r in rows]
    assert len(set(keys)) == len(keys), "variant keys must be unique"
    return rows


def main() -> int:
    genome_path = FIXTURE_DIR / "genome.fa"
    rows = build_variant_rows(genome_path)

    out_path = FIXTURE_DIR / "variants.tsv"
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {out_path} ({len(rows)} variants)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
