#!/usr/bin/env python3
"""Build the golden pipeline's predictions.tsv fixture (PH2-013).

Deliberately NOT model-inferred: an untrained/randomly-initialized model
produces near-identical values across variants and features, which would
hide exactly the kind of column-mixup bug this fixture exists to catch.
Instead every variant/feature's pred/delta values come from an explicit
deterministic formula so they are all mutually distinct.

Column layout matches ``VariantEffectPredictor.predict_variants``'s real
emission order (see ``src/prediction/predict.py``): the input variant
columns (as in ``variants.tsv``) followed by one
``ref_pred_<i>``/``alt_pred_<i>``/``delta_<i>`` triple PER FEATURE (not
grouped by kind).
"""
from __future__ import annotations

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "golden_pipeline"

N_FEATURES = 3  # matches training/features.txt (enhancer, CTCF_binding_site, promoter)


def _clamp01(x: float) -> float:
    return min(1.0, max(0.0, x))


def build_prediction_rows(variant_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out_rows: list[dict[str, str]] = []
    for i, var in enumerate(variant_rows):
        row = dict(var)
        for j in range(N_FEATURES):
            ref_pred = round(0.10 + 0.01 * i + 0.05 * j, 4)
            sign = 1 if (i + j) % 2 == 0 else -1
            magnitude = 0.02 + 0.03 * i + 0.07 * j
            alt_pred = round(ref_pred + sign * magnitude, 4)

            ref_pred = _clamp01(ref_pred)
            alt_pred = _clamp01(alt_pred)

            # The formula's magnitude term is never exactly 0, so no cell
            # naturally lands on delta==0 after clamping. Force exactly one
            # (variant 0, feature 0) so the golden fixture also exercises
            # the delta==0 case downstream (l2_delta contributions of 0, a
            # zero row in abs_delta_max ties, etc).
            if i == 0 and j == 0:
                alt_pred = ref_pred

            delta = round(alt_pred - ref_pred, 4)

            row[f"ref_pred_{j}"] = f"{ref_pred:.4f}"
            row[f"alt_pred_{j}"] = f"{alt_pred:.4f}"
            row[f"delta_{j}"] = f"{delta:.4f}"
        out_rows.append(row)
    return out_rows


def main() -> int:
    variants_path = FIXTURE_DIR / "variants.tsv"
    with variants_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        variant_fieldnames = list(reader.fieldnames or [])
        variant_rows = list(reader)

    pred_rows = build_prediction_rows(variant_rows)

    fieldnames = list(variant_fieldnames)
    for j in range(N_FEATURES):
        fieldnames += [f"ref_pred_{j}", f"alt_pred_{j}", f"delta_{j}"]

    out_path = FIXTURE_DIR / "predictions.tsv"
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pred_rows)

    print(f"Wrote {out_path} ({len(pred_rows)} variants x {N_FEATURES} features)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
