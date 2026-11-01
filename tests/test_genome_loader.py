import numpy as np

from src.data_loading.genome_loader import SequenceParser


def test_one_hot_callable_on_instance() -> None:
    result = SequenceParser().one_hot("ACGT")
    np.testing.assert_array_equal(result, np.eye(4))


def test_one_hot_unknown_base_maps_to_zero_row() -> None:
    result = SequenceParser.one_hot("ACRT")
    expected = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]], dtype=float
    )
    np.testing.assert_array_equal(result, expected)


def test_decode_round_trip_on_instance() -> None:
    parser = SequenceParser()
    encoded = parser.one_hot("ACGTN")
    assert parser.decode_one_hot(encoded) == "ACGTN"


def test_get_fasta_ref_on_instance() -> None:
    class _Rec:
        seq = "ACG"

    genome = {"chr1": [_Rec()] * 0 or type("S", (), {"__getitem__": lambda self, s: _Rec()})()}
    assert SequenceParser().get_fasta_ref(1, 1, 3, genome) == "ACG"
