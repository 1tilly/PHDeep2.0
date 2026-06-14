#!/usr/bin/env python3
"""
Download 1000 Genomes Phase 3 variant data for PAH-relevant genes.

Rather than downloading full per-chromosome VCF files (2–10 GB each),
this script uses tabix remote streaming to extract only the genomic
regions covering the genes of interest.

Requirements
------------
- tabix on PATH  (from htslib)
- bgzipped + tabix-indexed VCFs on the 1KG FTP (accessed remotely via tabix)

Usage
-----
    python scripts/download_1kg_variants.py \\
        --output data/1kg_pah_genes/ \\
        --assembly GRCh38

    # Custom gene list (overrides built-in PAH genes)
    python scripts/download_1kg_variants.py \\
        --output data/1kg_pah_genes/ \\
        --genes BMPR2 SOX17 \\
        --assembly GRCh38

    # GRCh37/hg19 (Phase 3 release, original assembly)
    python scripts/download_1kg_variants.py \\
        --output data/1kg_pah_genes/ \\
        --assembly GRCh37

Output
------
One VCF file per gene: <output>/<gene>_<chr>_<start>-<end>.vcf.gz
A manifest TSV: <output>/manifest.tsv (gene, chromosome, start, end, path)

References
----------
1000 Genomes Phase 3: https://doi.org/10.1038/nature15393
Data: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ (GRCh37)
      https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/ (GRCh38)
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# PAH gene definitions
# Coordinates are approximate gene body ± 50 kb flanks (GRCh38 / GRCh37).
# Extend flanks to capture regulatory regions relevant to the enhancer model.
# ---------------------------------------------------------------------------

@dataclass
class GeneRegion:
    symbol: str
    chrom_38: str
    start_38: int
    end_38: int
    chrom_37: str
    start_37: int
    end_37: int
    description: str = ""


# fmt: off
PAH_GENES: list[GeneRegion] = [
    GeneRegion("BMPR2",   "chr2",  202_378_547, 202_477_696, "2",  203_238_938, 203_390_000, "Bone morphogenetic protein receptor type 2"),
    GeneRegion("ACVRL1",  "chr12",  51_932_937,  51_949_073, "12",  51_932_937,  51_949_073, "Activin A receptor like type 1 (ALK1)"),
    GeneRegion("ENG",     "chr9",  127_815_022, 127_870_497, "9",  130_576_560, 130_631_384, "Endoglin"),
    GeneRegion("SMAD9",   "chr13",  36_834_664,  36_999_994, "13",  37_387_174,  37_549_530, "SMAD family member 9"),
    GeneRegion("KCNK3",   "chr2",   26_451_076,  26_508_804, "2",   26_602_350,  26_660_064, "Potassium two-pore domain channel subfamily K member 3"),
    GeneRegion("CAV1",    "chr7",  116_503_359, 116_527_461, "7",  116_184_893, 116_212_267, "Caveolin 1"),
    GeneRegion("SMAD1",   "chr4",  145_503_695, 145_611_723, "4",  146_418_695, 146_526_723, "SMAD family member 1"),
    GeneRegion("SMAD4",   "chr18",  51_028_394,  51_085_045, "18",  48_556_583,  48_611_411, "SMAD family member 4"),
    GeneRegion("GDF2",    "chr10",  48_397_524,  48_407_893, "10",  48_956_479,  48_966_985, "Growth differentiation factor 2 (BMP9)"),
    GeneRegion("ATP13A3", "chr3",  193_388_741, 193_469_736, "3",  194_348_741, 194_429_736, "ATPase 13A3"),
    GeneRegion("AQP1",    "chr7",  31_079_823,  31_097_978, "7",   30_950_131,  30_968_286, "Aquaporin 1"),
    GeneRegion("SOX17",   "chr8",  55_370_661,  55_376_490, "8",   56_370_661,  56_376_490, "SRY-box transcription factor 17"),
    GeneRegion("TBX4",    "chr17",  59_531_096,  59_607_125, "17",  59_531_096,  59_607_125, "T-box transcription factor 4"),
    GeneRegion("FBLN2",   "chr3",  13_506_724,  13_661_882, "3",   13_506_724,  13_661_882, "Fibulin 2"),
    GeneRegion("ABCC8",   "chr11",  17_391_709,  17_479_300, "11",  17_391_709,  17_479_300, "ATP binding cassette subfamily C member 8"),
]
# fmt: on

# ---------------------------------------------------------------------------
# 1000 Genomes FTP paths
# ---------------------------------------------------------------------------

# GRCh37 — Phase 3, 2013 release (widely used for association studies)
_1KG_GRCH37_TMPL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)

# GRCh38 — 30× high coverage resequencing (2022)
_1KG_GRCH38_TMPL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)


def _vcf_url(chrom_with_prefix: str, assembly: str) -> str:
    chrom = chrom_with_prefix.lstrip("chr")
    if assembly == "GRCh38":
        return _1KG_GRCH38_TMPL.format(chrom=chrom)
    return _1KG_GRCH37_TMPL.format(chrom=chrom)


def _tabix_extract(
    remote_vcf: str,
    chrom: str,
    start: int,
    end: int,
    dest: Path,
    flank: int = 50_000,
) -> bool:
    """
    Use `tabix` to stream-extract a region from a remote indexed VCF.
    Returns True on success.
    """
    region = f"{chrom}:{max(1, start - flank)}-{end + flank}"
    cmd = ["tabix", "-h", remote_vcf, region]
    dest.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open(dest.with_suffix(""), "wb") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, timeout=600)
        if proc.returncode != 0:
            logger.error("tabix failed for %s %s: %s", remote_vcf, region, proc.stderr.decode())
            return False

        # bgzip the output
        bgzip_cmd = ["bgzip", "-f", str(dest.with_suffix(""))]
        subprocess.run(bgzip_cmd, check=True, timeout=120)

        # Index
        tabix_idx_cmd = ["tabix", "-p", "vcf", str(dest)]
        subprocess.run(tabix_idx_cmd, check=True, timeout=60)
        return True
    except FileNotFoundError:
        logger.error(
            "tabix not found on PATH. Install htslib: https://www.htslib.org/download/"
        )
        return False
    except subprocess.TimeoutExpired:
        logger.error("tabix timed out for %s %s", remote_vcf, region)
        return False
    except Exception as exc:
        logger.error("Unexpected error for %s: %s", region, exc)
        return False


def _check_tabix() -> None:
    try:
        subprocess.run(["tabix", "--version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        logger.error(
            "tabix is required but not found on PATH.\n"
            "Install via:\n"
            "  NixOS:  nix-shell -p htslib\n"
            "  Ubuntu: sudo apt install tabix\n"
            "  conda:  conda install -c bioconda htslib"
        )
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Download 1000 Genomes Phase 3 variants for PAH genes via tabix."
    )
    parser.add_argument("--output", required=True, help="Output directory.")
    parser.add_argument(
        "--assembly",
        default="GRCh38",
        choices=["GRCh38", "GRCh37"],
        help="Reference assembly. GRCh38 uses high-coverage 30× data; GRCh37 uses Phase 3.",
    )
    parser.add_argument(
        "--genes",
        nargs="*",
        default=None,
        help="Gene symbols to download (default: all PAH genes). E.g. --genes BMPR2 SOX17",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=50_000,
        help="Flanking region around each gene body in bp (default: 50000).",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    _check_tabix()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_map = {g.symbol: g for g in PAH_GENES}
    if args.genes:
        missing = set(args.genes) - set(gene_map)
        if missing:
            logger.error("Unknown gene symbols: %s. Known: %s", missing, sorted(gene_map))
            sys.exit(1)
        targets = [gene_map[s] for s in args.genes]
    else:
        targets = PAH_GENES

    manifest_rows = []
    n_ok = n_fail = 0

    for gene in targets:
        if args.assembly == "GRCh38":
            chrom, start, end = gene.chrom_38, gene.start_38, gene.end_38
        else:
            chrom, start, end = gene.chrom_37, gene.start_37, gene.end_37

        remote_vcf = _vcf_url(chrom, args.assembly)
        chrom_plain = chrom.lstrip("chr")
        dest = out_dir / f"{gene.symbol}_chr{chrom_plain}_{start}-{end}.vcf.gz"

        logger.info("Gene %-10s  %s:%d-%d", gene.symbol, chrom, start, end)

        if args.dry_run:
            print(f"WOULD EXTRACT: {remote_vcf} [{chrom}:{start}-{end}] → {dest}")
            manifest_rows.append({
                "gene": gene.symbol,
                "chromosome": chrom,
                "start": start,
                "end": end,
                "path": str(dest),
                "description": gene.description,
            })
            continue

        if dest.exists() and dest.stat().st_size > 0:
            logger.info("  ✓ already exists, skipping.")
            n_ok += 1
        else:
            ok = _tabix_extract(remote_vcf, chrom, start, end, dest, flank=args.flank)
            if ok:
                n_ok += 1
                logger.info("  ✓ saved to %s", dest)
            else:
                n_fail += 1
                logger.warning("  ✗ failed for %s", gene.symbol)

        manifest_rows.append({
            "gene": gene.symbol,
            "chromosome": chrom,
            "start": start,
            "end": end,
            "path": str(dest),
            "description": gene.description,
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = out_dir / "manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    logger.info("Manifest written to %s", manifest_path)

    if not args.dry_run:
        logger.info("Done. Success: %d  Failed: %d", n_ok, n_fail)


if __name__ == "__main__":
    main()
