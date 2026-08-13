"""Tests for the SKAT-O input-building contract. No R/rpy2 required — all
R interaction is exercised through a fake backend."""
import numpy as np
import pandas as pd
import pytest

from src.statistical_testing.skat_o_test import (
    bh_fdr,
    build_skat_o_inputs,
    run_skat_o,
)


def _genotypes():
    return pd.DataFrame({
        "variant_id": ["1:100:A:T", "1:200:AT:A", "2:300:G:C"],
        "S1": [0, 1, 0],
        "S2": [1, 0, 2],
        "S3": [2, 2, 1],
    })


def _phenotype():
    return pd.Series([1, 0, 1], index=["S1", "S2", "S3"])


def _variant_table():
    return pd.DataFrame({
        "variant_id": ["1:100:A:T", "1:200:AT:A", "2:300:G:C"],
        "group": ["G1", "G1", "G2"],
        "eis_diff": [0.3, 0.4, 0.5],
    })


class FakeBackend:
    def __init__(self):
        self.null_model_calls = []
        self.skat_calls = []

    def run_null_model(self, y, X, out_type):
        self.null_model_calls.append((list(y), X, out_type))
        return "FAKE_NULL_MODEL"

    def run_skat(self, Z, null_model, weights, method):
        self.skat_calls.append((Z.copy(), null_model, None if weights is None else weights.copy(), method))
        return {"p_value": 0.1, "p_value_burden": 0.2, "p_value_skat": 0.15}


def test_module_imports_without_rpy2():
    import src.statistical_testing.skat_o_test as mod
    assert isinstance(mod._HAS_RPY2, bool)


def test_build_skat_o_inputs_orientation():
    inputs = build_skat_o_inputs(_variant_table(), _genotypes(), _phenotype(), weight_col="eis_diff")
    g1 = next(i for i in inputs if i.group == "G1")
    assert g1.Z.shape == (3, 2)   # samples x variants, NOT variants x features
    assert list(g1.y) == [1, 0, 1]
    assert len(g1.weights) == 2
    # Z[i,j] must equal genotypes.loc[variant_j, sample_i]
    geno = _genotypes().set_index("variant_id")
    for si, sample in enumerate(["S1", "S2", "S3"]):
        for vi, vid in enumerate(["1:100:A:T", "1:200:AT:A"]):
            assert g1.Z[si, vi] == geno.loc[vid, sample]


def test_build_skat_o_inputs_sample_order_follows_phenotype():
    genotypes = _genotypes()
    # swap S1/S2 column order
    genotypes = genotypes[["variant_id", "S2", "S1", "S3"]]
    inputs = build_skat_o_inputs(_variant_table(), genotypes, _phenotype(), weight_col="eis_diff")
    g1 = next(i for i in inputs if i.group == "G1")
    assert list(g1.y) == [1, 0, 1]
    geno = _genotypes().set_index("variant_id")
    for si, sample in enumerate(["S1", "S2", "S3"]):
        for vi, vid in enumerate(["1:100:A:T", "1:200:AT:A"]):
            assert g1.Z[si, vi] == geno.loc[vid, sample]


def test_build_skat_o_inputs_raises_on_variant_missing_from_genotypes():
    variant_table = _variant_table()
    variant_table = pd.concat([
        variant_table,
        pd.DataFrame({"variant_id": ["3:400:C:G"], "group": ["G2"], "eis_diff": [0.6]}),
    ], ignore_index=True)
    with pytest.raises(ValueError, match="3:400:C:G"):
        build_skat_o_inputs(variant_table, _genotypes(), _phenotype(), weight_col="eis_diff")


def test_build_skat_o_inputs_raises_on_sample_missing_from_genotypes():
    phenotype = pd.Series([1, 0, 1, 0], index=["S1", "S2", "S3", "S4"])
    with pytest.raises(ValueError, match="S4"):
        build_skat_o_inputs(_variant_table(), _genotypes(), phenotype, weight_col="eis_diff")


def test_build_skat_o_inputs_raises_on_empty_sample_intersection():
    phenotype = pd.Series([1, 0], index=["X1", "X2"])
    with pytest.raises(ValueError):
        build_skat_o_inputs(_variant_table(), _genotypes(), phenotype, weight_col="eis_diff")


def test_build_skat_o_inputs_raises_on_invalid_genotype_value():
    genotypes = _genotypes()
    genotypes.loc[0, "S1"] = 3
    with pytest.raises(ValueError):
        build_skat_o_inputs(_variant_table(), genotypes, _phenotype(), weight_col="eis_diff")


def test_build_skat_o_inputs_raises_on_nan_weight():
    variant_table = _variant_table()
    variant_table.loc[0, "eis_diff"] = np.nan
    with pytest.raises(ValueError):
        build_skat_o_inputs(variant_table, _genotypes(), _phenotype(), weight_col="eis_diff")


def test_build_skat_o_inputs_raises_on_negative_weight():
    variant_table = _variant_table()
    variant_table.loc[0, "eis_diff"] = -0.3
    with pytest.raises(ValueError, match="square"):
        build_skat_o_inputs(variant_table, _genotypes(), _phenotype(), weight_col="eis_diff")


def test_groups_below_min_variants_are_skipped():
    # NOTE: the base fixture's G2 has only 1 variant, which would itself be
    # excluded by an explicit min_variants=2 filter — that would contradict
    # this test's intent of isolating "singleton group gets skipped, others
    # don't". So here G2 is given a second variant, and the new G3 group is
    # the one left with a single variant to demonstrate the skip behavior.
    variant_table = _variant_table()
    variant_table = pd.concat([
        variant_table,
        pd.DataFrame({
            "variant_id": ["2:301:T:C", "3:400:C:G"],
            "group": ["G2", "G3"],
            "eis_diff": [0.55, 0.6],
        }),
    ], ignore_index=True)
    genotypes = _genotypes()
    genotypes = pd.concat([
        genotypes,
        pd.DataFrame({"variant_id": ["2:301:T:C", "3:400:C:G"], "S1": [0, 1], "S2": [1, 0], "S3": [2, 2]}),
    ], ignore_index=True)
    inputs = build_skat_o_inputs(variant_table, genotypes, _phenotype(), weight_col="eis_diff", min_variants=2)
    groups = {i.group for i in inputs}
    assert "G3" not in groups
    assert "G1" in groups
    assert "G2" in groups


def test_bh_fdr_hand_computed():
    assert list(bh_fdr([0.01, 0.02, 0.03, 0.04])) == pytest.approx([0.04, 0.04, 0.04, 0.04])
    assert list(bh_fdr([0.001, 0.5])) == pytest.approx([0.002, 0.5])
    assert list(bh_fdr([0.5])) == pytest.approx([0.5])
    # monotonicity: BH-adjusted values must be non-decreasing when sorted by original p-value
    result = bh_fdr([0.04, 0.01])
    assert result[1] <= result[0]


def test_run_skat_o_with_fake_backend_passes_correct_arguments():
    backend = FakeBackend()
    result = run_skat_o(  # noqa: F841 — result unused, this test only checks backend call args
        _variant_table(), _genotypes(), _phenotype(),
        weight_col="eis_diff", group_col="group", backend=backend,
    )
    assert len(backend.null_model_calls) == 1
    assert backend.null_model_calls[0][0] == [1, 0, 1]
    assert len(backend.skat_calls) == 2  # G1 and G2
    g1_call = next(c for c in backend.skat_calls if c[0].shape[1] == 2)
    assert list(g1_call[2]) == pytest.approx([0.3, 0.4])


def test_run_skat_o_output_schema():
    result = run_skat_o(_variant_table(), _genotypes(), _phenotype(), weight_col="eis_diff", group_col="group", backend=FakeBackend())
    assert set(result.columns) >= {"feature_id", "n_variants", "n_samples", "p_value", "p_value_burden", "p_value_skat", "q_value", "weight"}
    assert "effect_size" not in result.columns
    assert list(result["p_value"]) == sorted(result["p_value"])  # sorted by p_value ascending
    assert list(result["q_value"]) == pytest.approx(list(bh_fdr(result["p_value"])))


def test_run_skat_o_unweighted():
    backend = FakeBackend()
    result = run_skat_o(_variant_table(), _genotypes(), _phenotype(), weight_col=None, group_col="group", backend=backend)
    assert backend.skat_calls[0][2] is None  # weights=None passed to backend
    assert set(result["weight"]) == {"none"}
