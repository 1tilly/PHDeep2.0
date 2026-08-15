#!/usr/bin/env python3
"""Regenerate the golden pipeline's `expected/*.tsv` outputs (PH2-013).

Runs the real pipeline stages against the frozen fixture under
`tests/data/golden_pipeline/` and writes their outputs to
`tests/data/golden_pipeline/expected/`:

- `bed_to_training`: via the same subprocess CLI invocation used by
  `tests/test_bed_to_training_integration.py` (against the *separate*
  `tests/data/encode_dnase_fixture/`, not the golden_pipeline genome — this
  target exists to regenerate that fixture's `output/` directory, not a
  golden_pipeline output. It is included here as one entry point for "all
  fixture regeneration in one command", per the `--only` target list).
- `aggregate` + `stats`: via `tests.golden_utils.run_aggregate_stats`
  (in-process `LocalRunner`, `DeterministicSkatBackend` swapped in for the
  R-backed `RpyBackend`), writing `aggregate_scores.tsv`,
  `aggregate_genotypes.tsv`, `stats_results.tsv`.

Usage:
    python scripts/regenerate_golden_fixtures.py [--check] [--only TARGET] [-q]

Exit codes:
    0 - wrote (or, under --check, found no diff)
    1 - diff found (--check only)
    2 - missing optional dependency (pybedtools/bedtools) for a requested target
"""
from __future__ import annotations

import argparse
import difflib
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "golden_pipeline"
ENCODE_FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "encode_dnase_fixture"
EXPECTED_DIR = FIXTURE_DIR / "expected"

sys.path.insert(0, str(REPO_ROOT / "tests"))

TARGETS = ("bed_to_training", "aggregate", "stats", "all")


def _log(quiet: bool, *args: object) -> None:
    if not quiet:
        print(*args, file=sys.stderr)


def _run_bed_to_training(quiet: bool) -> dict[str, str]:
    """Run bed_to_training's CLI against tests/data/encode_dnase_fixture/
    and return {relative_output_name: text} for its three output files."""
    try:
        import pybedtools  # noqa: F401
    except ImportError:
        _log(quiet, "pybedtools not installed — skipping bed_to_training target")
        raise _MissingDependency("pybedtools")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        out_bed = tmp / "training_regions.bed"
        out_error = tmp / "read_errors.txt"
        out_features = tmp / "features.txt"

        cmd = [
            sys.executable,
            str(REPO_ROOT / "src" / "data_processing" / "bed_to_training.py"),
            "-m", str(ENCODE_FIXTURE_DIR / "metadata.tsv"),
            "-i", str(ENCODE_FIXTURE_DIR / "input_bed_files"),
            "-o", str(out_bed),
            "-e", str(out_error),
            "-f", str(out_features),
            "-a", "GRCh38",
        ]
        completed = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr + "\n" + completed.stdout)

        return {
            "training_regions.bed": out_bed.read_text(),
            "read_errors.txt": out_error.read_text(),
            "features.txt": out_features.read_text(),
        }


class _MissingDependency(Exception):
    pass


def _run_aggregate_stats(quiet: bool) -> dict[str, "pd.DataFrame"]:  # noqa: F821
    import golden_utils as gu

    with tempfile.TemporaryDirectory() as td:
        _log(quiet, "Running aggregate -> stats over the golden fixture...")
        return gu.run_aggregate_stats(FIXTURE_DIR, Path(td))


def _diff_text(name: str, actual: str, expected_path: Path) -> list[str]:
    expected = expected_path.read_text() if expected_path.exists() else ""
    if actual == expected:
        return []
    diff = difflib.unified_diff(
        expected.splitlines(keepends=True),
        actual.splitlines(keepends=True),
        fromfile=f"expected/{name}",
        tofile=f"actual/{name}",
    )
    return list(diff)


def _diff_df(name: str, actual_df, expected_path: Path, str_cols: list[str]) -> list[str]:
    actual_text = actual_df.to_csv(sep="\t", index=False, lineterminator="\n", na_rep="")
    return _diff_text(name, actual_text, expected_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="diff instead of writing")
    parser.add_argument("--only", choices=TARGETS, default="all")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args(argv)

    targets = TARGETS[:-1] if args.only == "all" else (args.only,)

    exit_code = 0
    diffs: list[str] = []

    if "bed_to_training" in targets:
        try:
            outputs = _run_bed_to_training(args.quiet)
        except _MissingDependency:
            return 2
        for name, text in outputs.items():
            expected_path = ENCODE_FIXTURE_DIR / "output" / name
            if args.check:
                diffs.extend(_diff_text(name, text, expected_path))
            else:
                expected_path.parent.mkdir(parents=True, exist_ok=True)
                expected_path.write_text(text)
                _log(args.quiet, f"Wrote {expected_path}")

    if "aggregate" in targets or "stats" in targets:
        import golden_utils as gu

        result = _run_aggregate_stats(args.quiet)
        name_map = {
            "aggregate": [("aggregate_scores.tsv", "scores"), ("aggregate_genotypes.tsv", "genotypes")],
            "stats": [("stats_results.tsv", "stats")],
        }
        for target in ("aggregate", "stats"):
            if target not in targets:
                continue
            for filename, key in name_map[target]:
                df = result[key]
                expected_path = EXPECTED_DIR / filename
                if args.check:
                    diffs.extend(_diff_df(filename, df, expected_path, str_cols=[]))
                else:
                    gu.write_expected_tsv(df, expected_path)
                    _log(args.quiet, f"Wrote {expected_path}")

    if args.check and diffs:
        print("".join(diffs))
        exit_code = 1

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
