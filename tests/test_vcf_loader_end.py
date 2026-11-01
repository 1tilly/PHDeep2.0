import importlib.util
import warnings
from pathlib import Path
from types import ModuleType


def _load_module() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "src" / "data_loading" / "vcf_loader.py"
    spec = importlib.util.spec_from_file_location("vcf_loader_under_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_persist_to_feather_computes_end_without_positional_access(tmp_path: Path) -> None:
    module = _load_module()
    vcf_path = tmp_path / "variants.vcf"
    vcf_path.write_text(
        "# [1]CHROM\t[2]POS\t[3]REF\t[4]ALT\t[5]S000001:GT\t\n"
        "1\t100\tAC\tA\t0/1\t\n"
        "1\t200\tA\tG\t0/0\t\n"
    )
    out_path = tmp_path / "out.feather"

    parser = module.BCFParser()
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        result = parser.persist_to_feather(str(vcf_path), ["S000001"], str(out_path))

    assert result is not None
    assert list(result["start"]) == [100, 200]
    assert list(result["end"]) == [101, 200]
    assert list(result["pah_var"]) == [True, False]
    assert "length" not in result.columns
    assert out_path.exists()
