from typing import Any

import pytest
import requests

from src.data_loading import ensembl_api
from src.data_loading.ensembl_api import EnsemblAPI


class _FakeResponse:
    def __init__(self, text: str, ok: bool = True) -> None:
        self.text = text
        self.ok = ok

    def raise_for_status(self) -> None:
        if not self.ok:
            raise requests.HTTPError("boom")


def _fake_get(calls: list[tuple[int, int]]) -> Any:
    def fake(url: str, **kwargs: Any) -> _FakeResponse:
        region = url.rsplit("/", 1)[-1]
        coords = region.split(":", 1)[1]
        s_str, e_str = coords.split("..")
        s, e = int(s_str), int(e_str)
        calls.append((s, e))
        return _FakeResponse("A" * (e - s + 1))

    return fake


def test_chunked_region_returns_full_length(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[int, int]] = []
    monkeypatch.setattr(ensembl_api.requests, "get", _fake_get(calls))
    api = EnsemblAPI()
    (chrom, start, end), ref = api.get_ref_for_region("1", 0, 15_000_000, seq_len=1000)
    assert chrom == "1"
    assert isinstance(start, int) and isinstance(end, int)
    assert len(ref) == end - start + 1
    assert len(calls) == 2
    # chunks must not overlap and must be contiguous
    for (_, prev_e), (next_s, _) in zip(calls, calls[1:]):
        assert next_s == prev_e + 1
    assert calls[0][0] == start
    assert calls[-1][1] == end


def test_http_error_raised_before_length_assertion(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake(url: str, **kwargs: Any) -> _FakeResponse:
        return _FakeResponse("error body", ok=False)

    monkeypatch.setattr(ensembl_api.requests, "get", fake)
    api = EnsemblAPI()
    with pytest.raises(requests.HTTPError):
        api.get_reference("1", 100, 1099, seq_len=1000)


def test_request_has_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    seen: dict[str, Any] = {}

    def fake(url: str, **kwargs: Any) -> _FakeResponse:
        seen.update(kwargs)
        return _FakeResponse("A" * 10)

    monkeypatch.setattr(ensembl_api.requests, "get", fake)
    EnsemblAPI().get_reference("1", 0, 9)
    assert seen.get("timeout") is not None
