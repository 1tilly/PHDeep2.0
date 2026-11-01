from __future__ import annotations

import importlib.util
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import requests

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "download_encode_training.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("download_encode_training", _SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _FakeResp:
    def __init__(self, chunks: list[bytes], fail: bool) -> None:
        self._chunks = chunks
        self._fail = fail

    def __enter__(self) -> _FakeResp:
        return self

    def __exit__(self, *args: Any) -> None:
        return None

    def raise_for_status(self) -> None:
        return None

    def iter_content(self, chunk_size: int) -> Iterator[bytes]:
        yield from self._chunks
        if self._fail:
            raise requests.ConnectionError("boom")


def test_interrupted_download_not_treated_as_complete(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mod = _load()
    dest = tmp_path / "out" / "ENCFF000.bed.gz"

    def bad_get(*args: Any, **kwargs: Any) -> _FakeResp:
        return _FakeResp([b"partial-data"], fail=True)

    monkeypatch.setattr(mod.requests, "get", bad_get)
    _, ok, _ = mod._download_file("http://x", dest)
    assert not ok
    assert not dest.exists()
    assert not dest.with_name(dest.name + ".part").exists()

    def good_get(*args: Any, **kwargs: Any) -> _FakeResp:
        return _FakeResp([b"full-", b"data"], fail=False)

    monkeypatch.setattr(mod.requests, "get", good_get)
    _, ok, msg = mod._download_file("http://x", dest)
    assert ok
    assert msg != "already exists"
    assert dest.read_bytes() == b"full-data"
