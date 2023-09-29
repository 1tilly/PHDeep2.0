import pandas as pd


class VariantParser(BaseDataLoader):

    def load_gene(fp): 
        var_df = pd.read_pickle(fp)
        return var_df[var_df.columns[:5].tolist() + var_df.columns[-6:].tolist() + var_df.columns[5:-6].tolist()]


    def find_variant_in_reference(variant_coord, reference, ref_start, length=2000, focus=200):
        """
		@param variant_coord: (chrom, pos,ref, alt)
		"""
        index = int(variant_coord[1]) - int(ref_start)
        ref_length = len(variant_coord[2])
        ref_ens = reference[index:(index + ref_length)]
  
        assert len(variant_coord[3]) <= focus, f"The alternate length ({len(variant_coord[3])}) is longer than the focus of the network: {focus} "
        result = []
        start_index = int(index - (length - focus) / 2 - (focus - len(variant_coord[3])) / 2)
        end_index = int(index + ref_length + (length - focus) / 2 + (focus - len(variant_coord[3])) / 2)
        # The following is handling all cases: Substitution, Insertion and Deletion. It takes the whole reference up to the
        #   variant start, then appends the variant itself (in case of a deletion it is shorter than before, in case of an insertion longer)
        #   then adds the last bit of the reference, from starting index+ref_length to our result. This way the length difference of
        #   InDels is handled.
        result.extend(reference[start_index:index])
        result.extend(variant_coord[3])
        result.extend(reference[index + ref_length:end_index])
        result = "".join(result)
        diff = length - len(result)
        if diff > 0:
            print(f"Adding {diff} basepairs in the end.")
            result.append(reference[end_index:end_index + diff])
        elif diff < 0:
            print(f"Erasing {diff} basepairs in the end")
            result = result[:length]
        return result, reference[start_index:start_index + length]


    def iterate_through_mutations_in_sequence(ref_seq, ref_start, var_df, seq_len=2000):
        for index, var in var_df.iterrows():
            yield VariantParser.find_variant_in_reference((var.chromosome, var.start, var.reference, var.alternate), ref_seq, ref_start, seq_len)

