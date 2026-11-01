"""--check mode must serialize with the same writer as write mode."""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path
from types import ModuleType

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "regenerate_golden_fixtures.py"


def _load_script() -> ModuleType:
    spec = importlib.util.spec_from_file_location("regen_golden_under_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _StubDF:
    def to_csv(self, *args: object, **kwargs: object) -> str:
        return "a\tb\n1.0\t2.0\n"


def _fake_write_expected_tsv(df: object, path: Path) -> None:
    Path(path).write_text("a\tb\n1.000000\t2.000000\n")


def test_diff_df_uses_write_expected_tsv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = types.ModuleType("golden_utils")
    fake.write_expected_tsv = _fake_write_expected_tsv  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "golden_utils", fake)
    mod = _load_script()

    expected_path = tmp_path / "stats_results.tsv"
    _fake_write_expected_tsv(_StubDF(), expected_path)

    assert mod._diff_df("stats_results.tsv", _StubDF(), expected_path) == []


def test_diff_df_has_no_str_cols_param() -> None:
    import inspect

    mod = _load_script()
    assert "str_cols" not in inspect.signature(mod._diff_df).parameters
