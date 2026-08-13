from __future__ import annotations

import pandas as pd


def _padded_slice(seq: str, start: int, end: int) -> str:
    width = end - start
    if width <= 0:
        return ""
    lo, hi = max(0, start), min(len(seq), end)
    core = seq[lo:hi] if lo < hi else ""
    left = min(width, max(0, -start))
    right = width - left - len(core)
    return "N" * left + core + "N" * right


class VariantParser():

    @staticmethod
    def load_gene(fp):
        var_df = pd.read_pickle(fp)
        return var_df[var_df.columns[:5].tolist() + var_df.columns[-6:].tolist() + var_df.columns[5:-6].tolist()]


    @staticmethod
    def find_variant_in_reference(
        variant_coord, reference: str, ref_start: int, length: int = 2000, focus: int = 200
    ) -> tuple[str, str]:
        """
		@param variant_coord: (chrom, pos, ref, alt)

		Splices the variant's alternate allele into the reference sequence and returns
		a window around it, alongside the matching unmodified reference window. Both
		returned sequences are exactly `length` long and mutually aligned (same offsets
		outside the variant). Where the window falls outside the bounds of `reference`,
		it is padded with 'N'. `focus` only constrains how long the alternate allele may
		be (see the assertion below) and does not otherwise affect the output.
		"""
        index = int(variant_coord[1]) - int(ref_start)
        ref_length = len(variant_coord[2])
        alt = variant_coord[3]

        assert len(alt) <= focus, f"The alternate length ({len(alt)}) is longer than the focus of the network: {focus} "

        pre_len = -(-(length - len(alt)) // 2)   # ceil((length - len(alt)) / 2)
        post_len = length - pre_len - len(alt)
        start_index = index - pre_len
        end_index = index + ref_length + post_len

        alt_seq = (
            _padded_slice(reference, start_index, index)
            + alt
            + _padded_slice(reference, index + ref_length, end_index)
        )
        ref_window = _padded_slice(reference, start_index, start_index + length)

        return alt_seq, ref_window


    @staticmethod
    def iterate_through_mutations_in_sequence(ref_seq, ref_start, var_df, seq_len=2000):
        for index, var in var_df.iterrows():
            yield VariantParser.find_variant_in_reference((var.chromosome, var.start, var.reference, var.alternate), ref_seq, ref_start, seq_len)
