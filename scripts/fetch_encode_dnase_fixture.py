#!/usr/bin/env python3
"""Download a tiny ENCODE DNase narrowPeak fixture and write matching metadata."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import requests


ENCODE_BASE = "https://www.encodeproject.org"
# Small released DNase-seq GRCh38 narrowPeak files.
DEFAULT_ACCESSIONS = ("ENCFF887HIM", "ENCFF668BIL")


def _dataset_accession(dataset_ref: str) -> str:
    return dataset_ref.strip("/").split("/")[-1]


def fetch_fixture(output_dir: Path, accessions: tuple[str, ...]) -> None:
    input_dir = output_dir / "input_bed_files"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    session = requests.Session()
    session.headers.update({"accept": "application/json"})

    for acc in accessions:
        meta_url = f"{ENCODE_BASE}/files/{acc}/?format=json"
        meta = session.get(meta_url, timeout=60)
        meta.raise_for_status()
        payload = meta.json()

        href = payload["href"]
        download_url = f"{ENCODE_BASE}{href}"
        file_out = input_dir / f"{acc}.bed.gz"

        bed = session.get(download_url, timeout=120)
        bed.raise_for_status()
        file_out.write_bytes(bed.content)

        rows.append(
            {
                "File accession": acc,
                "File format": "bed narrowPeak",
                "File type": "bed",
                "Experiment accession": _dataset_accession(payload["dataset"]),
                "Output type": payload.get("output_type", "peaks"),
                "Assay": payload.get("assay_term_name", "DNase-seq"),
                "Biosample term name": payload.get("biosample_ontology", {}).get("term_name", "unknown"),
                "Experiment target": "",
                "File download URL": download_url,
                "Biological replicate(s)": ",".join(map(str, payload.get("biological_replicates", []))),
                "Technical replicate(s)": ",".join(map(str, payload.get("technical_replicates", []))),
                "Biosample treatments": "",
                "File assembly": payload.get("assembly", "GRCh38"),
                "File analysis title": "ENCODE4 GRCh38",
                "Biosample organism": "Homo sapiens",
            }
        )

    meta_df = pd.DataFrame(rows)
    meta_path = output_dir / "metadata.tsv"
    meta_df.to_csv(meta_path, sep="\t", index=False)

    print(f"Downloaded {len(rows)} files to: {input_dir}")
    print(f"Wrote metadata: {meta_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch ENCODE DNase narrowPeak fixture files.")
    parser.add_argument(
        "--out-dir",
        default="tests/data/encode_dnase_fixture",
        help="Output directory for metadata.tsv and input_bed_files/.",
    )
    parser.add_argument(
        "--accessions",
        nargs="+",
        default=list(DEFAULT_ACCESSIONS),
        help="ENCODE file accessions to download.",
    )
    args = parser.parse_args()

    fetch_fixture(Path(args.out_dir), tuple(args.accessions))


if __name__ == "__main__":
    main()
