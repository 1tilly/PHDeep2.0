import itertools
from typing import Iterable
import pandas as pd


class GeneHancerParser:
    """
    In my project, GFF is only used int he context of GeneHancer. This means, the code is highly specialised for this case.
    """
    def __init__(self, path):
        self.path = path
        self.df = None

    def save_df(self, path=None):
        if path is None:
            path = self.path

        self.df.write_pickle(path, sep="\t")
        return self.df

    def load_df(self, path=None):
        if path is None:
            path = self.path

        self.df = pd.read_pickle(path)
        return self.df

    def load_gff_to_df(self, path=None, split_attr=True):
        """
        Inspired and adjusted from gffpandas
        """
        if path is None:
            path = self.path

        df = pd.read_csv(path, sep="\t", header=0)
        df["chrom"] = df["#chrom"].apply(lambda x: x[3:])

        if split_attr:
            def parse_attribute(attributes):
                kvpairs = attributes.split(";")
                res = {}
                for kv in kvpairs:
                    k, v = kv.split(sep="=")
                    if k in res.keys():
                        res[k].append(v)
                    else:
                        res[k] = [v]
                return res
            attribute_df = df.copy()

            attribute_df["at_dic"] = attribute_df.attributes.apply(
                lambda attributes: parse_attribute(attributes))
            attribute_df["at_dic_keys"] = attribute_df["at_dic"].apply(
                lambda at_dic: list(at_dic.keys()))
            merged_attribute_list = list(
                itertools.chain.from_iterable(attribute_df["at_dic_keys"]))
            nonredundant_list = sorted(list(set(merged_attribute_list)))

            for atr in nonredundant_list:
                if atr == "score":
                    atr_alt = "attr_score"
                    df[atr_alt] = attribute_df["at_dic"].apply(
                        lambda at_dic: at_dic.get(atr))

                else:
                    df[atr] = attribute_df["at_dic"].apply(
                        lambda at_dic: at_dic.get(atr))
            #pd.concat(df, df_attributes)
            df = df.explode("genehancer_id")
            self.df = df.drop(columns=["#chrom", "attributes"])
            return self.df
        else:
            # TODO this needs to be handled
            raise NotImplementedError(
                "You reached a case that is not implemented!")

    def check_for_gene_connections(self, gene_list: Iterable[str]) -> pd.DataFrame:
        if self.df is None:
            self.load_gff_to_df()

        return self.df[self.df.connected_gene.apply(lambda x: len(x+gene_list) > len(set(x+gene_list)))]

    def get_region(self, region):
        if self.df is None:
            self.load_gff_to_df()

        return self.df[(self.df.chrom == str(region[0])) & (self.df.start > region[1]) & (self.df.start < region[2]) & (self.df.end > region[1]) & (self.df.end < region[2])]

    def get_geneSNPlist(self, gene=None):
        if "snp_list" not in self.df.columns:
            raise KeyError("There is no SNPlist column in the dataframe, yet") 
        if gene is not None:
            """
            Explode by connected_gene, then merge by unique value and append the snp_lists
            """
            df = self.df.explode(column="connected_gene").groupby(by=["connected_gene"]).agg({"genehancer_id":lambda x: set(x), "snp_list": lambda x:x.sum()})          
        else:
            """
            ToDo: Here I should implement a way to run this for each and every unique gene in the geneHancer DF, efficiently
            ! Efficiency matters, a for loop would take too long
            """
            raise NotImplementedError("Implementation missing in geneHancer_parser.py, function get_geneSNPlist")
        return df

    def add_to_snplist(self, snpID_list, enhancer=None):
        if "snp_list" not in self.df.columns:
            self.df.loc["snp_list"] = ""
        if enhancer is not None:
            self.df[self.df["genehancer_id"] == enhancer, "snp_list"] = self.df.loc[self.df["genehancer_id"]
                                                                                    == enhancer, "snp_list"].astype(str) + " ".join(snpID_list)
        else:
            """
            ToDo: Here I should create functionality that uses the "get_region" function for each snp and add it to the according enhancer
            ! Efficiency matters, a for loop would take too long
            """
            raise NotImplementedError("Implementation missing in geneHancer_parser.py, function add_to_snplist")

        return self.df

    def snpList_to_list(self, df=None, sep=" "):
        if df is None:
            # This way a dataframe from outside can be passed to be finalized
            df = self.df
        df = df.snp_list.apply(lambda x: x.trim().split(sep))
        return df        
