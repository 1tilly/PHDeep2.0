from pathlib import Path

from src.data_loading.vcf_loader import BCFParser


def _write_query(path: Path, gts: list[tuple[str, str]]) -> None:
    header = "# [1]CHROM\t[2]POS\t[3]REF\t[4]ALT\t[5]AC\t[6]AN\t[7]LP000001:GT\t[8]LP000002:GT\t\n"
    rows = "".join(
        f"1\t{100 + i}\tA\tC,G\t2\t4\t{a}\t{b}\t\n" for i, (a, b) in enumerate(gts)
    )
    path.write_text(header + rows)


def test_multiallelic_carriers_detected(tmp_path: Path) -> None:
    vcf_path = tmp_path / "q.vcf"
    _write_query(
        vcf_path,
        [("2/2", "0/0"), ("0|2", "0/0"), ("0/0", "0/0"), ("./.", "0/1"), ("0/10", "0/0")],
    )
    out = BCFParser().persist_to_feather(
        str(vcf_path), ["LP000001", "LP000002"], str(tmp_path / "out.feather")
    )
    assert out is not None
    assert out["pah_var"].tolist() == [True, True, False, True, True]


def test_non_carriers_with_missing_not_flagged(tmp_path: Path) -> None:
    vcf_path = tmp_path / "q.vcf"
    _write_query(vcf_path, [("./.", "0/0"), ("0|0", "./.")])
    out = BCFParser().persist_to_feather(
        str(vcf_path), ["LP000001", "LP000002"], str(tmp_path / "out.feather")
    )
    assert out is not None
    assert out["pah_var"].tolist() == [False, False]
