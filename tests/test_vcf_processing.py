"""Tests for variant-into-reference splicing."""
import pytest

from src.data_processing.vcf_processing import VariantParser, _padded_slice

# Positionally distinguishable reference: character at index i is str(i % 10),
# so any misalignment is visible by eye.
REF = "0123456789" * 4
LENGTH, FOCUS = 20, 10


def _call(pos, ref_allele, alt_allele, reference=REF, ref_start=0):
    return VariantParser.find_variant_in_reference(
        ("chr1", pos, ref_allele, alt_allele), reference, ref_start, LENGTH, FOCUS
    )


def test_snv_in_middle_splices_exact_sequences():
    assert _call(120, "0", "A", ref_start=100) == (
        "0123456789A123456789",
        "01234567890123456789",
    )


def test_insertion_in_middle_splices_exact_sequences():
    assert _call(120, "0", "AGG", ref_start=100) == (
        "123456789AGG12345678",
        "12345678901234567890",
    )


def test_deletion_in_middle_splices_exact_sequences():
    alt_seq, ref_window = _call(120, "012", "0", ref_start=100)
    assert (alt_seq, ref_window) == (
        "01234567890345678901",
        "01234567890123456789",
    )
    assert "12" not in alt_seq[11:13]


@pytest.mark.parametrize(
    "ref_allele,alt_allele",
    [("0", "A"), ("0", "AGG"), ("012", "0"), ("0", "AAAAAAAA"), ("01234", "0")],
)
def test_output_is_always_exactly_length(ref_allele, alt_allele):
    alt_seq, ref_window = _call(120, ref_allele, alt_allele, ref_start=100)
    assert len(alt_seq) == LENGTH
    assert len(ref_window) == LENGTH


def test_focus_parameter_does_not_affect_output():
    results = [
        VariantParser.find_variant_in_reference(
            ("chr1", 120, "0", "A"), REF, 100, LENGTH, focus
        )
        for focus in (10, 4, 2, 1)
    ]
    assert all(r == results[0] for r in results)


def test_alt_longer_than_focus_raises():
    with pytest.raises(AssertionError, match="longer than the focus"):
        VariantParser.find_variant_in_reference(
            ("chr1", 120, "0", "AAAAA"), REF, 100, LENGTH, 4
        )


def test_ref_window_and_alt_seq_are_aligned_outside_variant():
    alt_seq, ref_window = _call(120, "0", "A", ref_start=100)
    assert alt_seq[:10] == ref_window[:10]
    assert alt_seq[11:] == ref_window[11:]
    assert alt_seq[10] == "A"


def test_variant_at_reference_start_pads_left_with_n():
    assert _call(100, "0", "A", ref_start=100) == (
        "NNNNNNNNNNA123456789",
        "NNNNNNNNNN0123456789",
    )


def test_variant_near_reference_start_pads_left_with_n():
    assert _call(3, "3", "A", ref_start=0) == (
        "NNNNNNN012A456789012",
        "NNNNNNN0123456789012",
    )


def test_near_start_ref_window_contains_the_variant_position():
    _, ref_window = _call(3, "3", "A", ref_start=0)
    assert ref_window[7 + 3] == "3"


def test_variant_near_reference_end_pads_right_with_n():
    assert _call(36, "6", "A", ref_start=0) == (
        "6789012345A789NNNNNN",
        "67890123456789NNNNNN",
    )


def test_reference_shorter_than_window_pads_both_sides():
    assert _call(4, "4", "A", reference="01234567", ref_start=0) == (
        "NNNNNN0123A567NNNNNN",
        "NNNNNN01234567NNNNNN",
    )


def test_padding_is_uppercase_n_only():
    alt_seq, _ = _call(3, "3", "A", ref_start=0)
    assert set(alt_seq) <= set("0123456789AN")
    assert alt_seq.count("N") == 7


# --- _padded_slice ---

_PS_SEQ = "0123456789"


def test_padded_slice_in_range_no_padding():
    result = _padded_slice(_PS_SEQ, 2, 6)
    assert result == "2345"
    assert len(result) == 4


def test_padded_slice_starts_before_zero():
    result = _padded_slice(_PS_SEQ, -3, 4)
    assert result == "NNN0123"
    assert len(result) == 7


def test_padded_slice_ends_past_len():
    result = _padded_slice(_PS_SEQ, 7, 13)
    assert result == "789NNN"
    assert len(result) == 6


def test_padded_slice_spans_both_edges():
    result = _padded_slice(_PS_SEQ, -2, 12)
    assert result == "NN0123456789NN"
    assert len(result) == 14


def test_padded_slice_entirely_before_zero():
    result = _padded_slice(_PS_SEQ, -5, -1)
    assert result == "NNNN"
    assert len(result) == 4


def test_padded_slice_entirely_past_end():
    result = _padded_slice(_PS_SEQ, 12, 16)
    assert result == "NNNN"
    assert len(result) == 4


def test_padded_slice_width_zero_or_negative():
    assert _padded_slice(_PS_SEQ, 5, 5) == ""
    assert _padded_slice(_PS_SEQ, 5, 3) == ""
    assert len(_padded_slice(_PS_SEQ, 5, 5)) == 0
    assert len(_padded_slice(_PS_SEQ, 5, 3)) == 0
