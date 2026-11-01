import gzip
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("pybedtools")

from src.data_loading.bed_loader import BEDParser  # noqa: E402


def test_read_bed_file_unknown_analysis_title(tmp_path: Path) -> None:
    lines = [
        "chr1\t100\t200\tpeak1\t500\t.\t1.0\t2.0\t3.0\t50\n",
        "chr2\t300\t450\tpeak2\t600\t.\t1.0\t2.0\t3.0\t60\n",
    ]
    with gzip.open(tmp_path / "ENCFFTEST.bed.gz", "wt") as f:
        f.writelines(lines)
    row = pd.Series(
        {
            "File accession": "ENCFFTEST",
            "File type": "bed",
            "Experiment accession": "ENCSR000",
            "Assay": "DNase-seq",
            "Experiment target": None,
            "Biosample term name": "K562",
            "Technical replicate(s)": "1_1",
            "Biological replicate(s)": "1",
            "File assembly": "GRCh38",
            "File analysis title": "ENCODE4 v1.10.1 GRCh38",
        }
    )
    parser = BEDParser()
    result = parser.read_bed_file(str(tmp_path), row, "GRCh38")
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 2
    assert "proj_order" not in result.columns
