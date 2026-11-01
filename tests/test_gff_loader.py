import pandas as pd

from src.data_loading.gff_loader import GeneHancerParser


def _parser(df: pd.DataFrame) -> GeneHancerParser:
    parser = GeneHancerParser("unused")
    parser.df = df
    return parser


def test_duplicate_in_gene_list_does_not_match() -> None:
    parser = _parser(pd.DataFrame({"connected_gene": [["X"]]}))
    assert len(parser.check_for_gene_connections(["A", "A"])) == 0


def test_duplicate_within_row_does_not_match() -> None:
    parser = _parser(pd.DataFrame({"connected_gene": [["X", "X"]]}))
    assert len(parser.check_for_gene_connections(["A"])) == 0


def test_tuple_gene_list_matches() -> None:
    parser = _parser(pd.DataFrame({"connected_gene": [["X"], ["Y"]]}))
    result = parser.check_for_gene_connections(("X",))
    assert len(result) == 1


def test_get_region_inclusive_bounds() -> None:
    parser = _parser(
        pd.DataFrame({"chrom": ["1", "1"], "start": [100, 150], "end": [200, 250]})
    )
    result = parser.get_region((1, 100, 200))
    assert len(result) == 1
    assert int(result.iloc[0]["start"]) == 100


def test_get_gene_snplist_filters_by_gene() -> None:
    parser = _parser(
        pd.DataFrame(
            {
                "connected_gene": [["A"], ["B"], ["A", "B"]],
                "genehancer_id": ["e1", "e2", "e3"],
                "snp_list": ["s1 ", "s2 ", "s3 "],
            }
        )
    )
    result = parser.get_geneSNPlist("A")
    assert list(result.index) == ["A"]
    assert result.loc["A", "genehancer_id"] == {"e1", "e3"}
