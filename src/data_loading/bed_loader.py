import os.path

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import logging

from pybedtools import BedTool

# Initialise the logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger('bed_loader')
class BEDParser():
    """
    ! This is really just an ENCODE parser! It should be split into basic bed file functions and a ENCODE parser inheriting from that
    
    """

    def __init__():

        # We use the order of projects to determine which of the assay_values to drop
        proj_order_37 = ['Lab custom hg19', 'Lab custom hg19 V7', 'Lab custom hg19 V10', 'Lab custom hg19 V19',
                        'ENCODE3 hg19', 'ENCODE3 hg19 V19']
        proj_order_38 = ['Lab custom GRCh38', 'Lab custom GRCh38 V24', 'ENCODE3 GRCh38',
                        'ENCODE4 GRCh38', 'ENCODE3 GRCh38 V24', 'ENCODE4 GRCh38 V29',
                        'ENCODE4 v1.1.0 GRCh38 V29', 'ENCODE4 v1.2.1 GRCh38 V29', 'ENCODE4 v1.4.0 GRCh38',
                        'ENCODE4 v1.5.0 GRCh38', 'ENCODE4 v1.5.1 GRCh38', 'ENCODE4 v1.6.0 GRCh38', 'ENCODE4 v1.6.1 GRCh38',
                        'ENCODE4 v1.7.0 GRCh38', 'ENCODE4 v1.8.0 GRCh38', 'ENCODE4 v1.9.1 GRCh38',
                        'ENCODE4 v3.0.0-alpha.2 GRCh38', 'ENCODE4 v3.0.0 GRCh38']

        path_meta_file = "/home/tt419/Projects/DeepLearning/data/ENCODE_v1_orig_BED/metadata_v1.tsv"
        path_bed_out = "/home/tt419/Projects/DeepLearning/data/ENCODE_v1_orig_BED/"

        organism = "Homo sapiens"
        assays = ["DNase-seq", "TF ChIP-seq", "Histone ChIP-seq", "Control ChIP-seq", "Hi-C", "ATAC-seq"]
        columns_of_interest = ["File accession", "File format", "File type", "Experiment accession", "Output type", "Assay",
                            "Biosample term name", "Experiment target", "File download URL", "Biological replicate(s)",
                            "Technical replicate(s)", "Biosample treatments", "File assembly", "File analysis title"]

        

    def load_prepare_meta(path:str, organism:str, assays:list, c_o_i:list, file_type:list=["bed narrowPeak"], assembly:str="hg19"):
        """
        @parameter c_o_i  := Columns of interest, the columns of the meta-file, that are relevant for the selection.
        """
        with open(path, "r") as f:
            file = pd.read_csv(f, sep="\t")

        
        df = file[(file["Biosample organism"] == organism) &
            (file["Assay"].isin(assays))][c_o_i]
        df_ft = df[df["File format"].isin(file_type)]
        
        df_as = df_ft[df_ft["File assembly"] == assembly]
        
        # dropping all rows, that are NA in any of these columns, as they become unusable else.
        df_na = df_as.dropna(subset=["File analysis title", "Biosample term name", "Experiment accession", "Assay"])
        # The project order makes sure, that the latest version of each experiment is being used withing the assembly. 


        
        return df_na

    def generate_feature_exp_dict(meta_df):
        """
        This function expects the meta_df to hold either DNase-seq or specifically named targets (like in TF/Histone assays, naming the TF or histone modification as target)
        """
        feature_exp_dict = {}
        for exp in meta_df["Experiment accession"].unique():
            if exp is None:
                continue
            for _, row in meta_df[meta_df["Experiment accession"] == exp].iterrows():
                if row['Assay'] == "DNase-seq":
                    target = "DNase"
                else:
                    target =  "_".join(row['Experiment target'].split(" "))
                if not pd.notna(row['Biosample treatments']): 
                    bio_treat = "None"
                else:
                    bio_treat = "_".join(row['Biosample treatments'].split(" "))
                if len(row['Biosample term name'].split(" ")) > 1:
                    biosam = "_".join(row['Biosample term name'].split(" "))
                else:
                    biosam = row['Biosample term name']
                feature_name = f"{biosam}|{target}|{bio_treat}"
                if feature_name not in feature_exp_dict.keys() or feature_exp_dict[feature_name] is None:
                    
                    feature_exp_dict[feature_name] = [exp,]
                else:
                    
                    feature_exp_dict[feature_name].append(exp)

        return feature_exp_dict



    def generate_feature_list(feature_exp_dict):
        return sorted(feature_exp_dict.keys(), key=lambda x: (x.split("|")[1], x.split("|")[0], x.split("|")[2]))


    def read_bed_file(path:str, file, read_error_log, proj_order:str):
    
        file_path = os.path.join(path, f"{file['File accession']}.{file['File type']}.gz")
        try:
            snps = BedTool(file_path)
        except FileNotFoundError:
            logger.error(f"FNF_BED: {file_path}\n")
            return 0


        clean_chrom = lambda x: f"chr{x}" if "chr" not in x else x  # chr prefix is needed by selene

        try:
            res = snps.to_dataframe(header=None)
        except EOFError:
            logger.error(f"EOF: {file_path}\n")
            return 0
        if isinstance(res.columns[1], np.int64):
            res = res.iloc[1:, [0, 1, 2, 4]]
        else:
            res = res.iloc[:, [0, 1, 2, 4]]
        res["chrom"] = res["chrom"].apply(clean_chrom)
        chromosomes = list(map(lambda x: f"chr{x}",range(1,23))) + ["chrX", "chrY"]
        res = res[res.chrom.isin(chromosomes)]
        res["origin_file"] = file_path
        res["exp_acc"] = file["Experiment accession"]
        res["assay_name"] = file["Assay"]
        res["exp_target"] = file["Experiment target"]
        res["bio_name"] = file["Biosample term name"]
        res["tech_repl"] = file["Technical replicate(s)"]
        res["bio_repl"] = file["Biological replicate(s)"]
        res["assembly"] = file["File assembly"]
        res["analysis"] = file["File analysis title"]


        res["exp_acc"] = res["exp_acc"].astype(str)
        res["assay_name"] = res["assay_name"].astype(str)
        res["exp_target"] = res["exp_target"].astype(str)
        res["bio_name"] = res["bio_name"].astype(str)
        res["tech_repl"] = res["tech_repl"].astype(str)
        res["bio_repl"] = res["bio_repl"].astype(str)
        res["assembly"] = res["assembly"].astype(str)


        res["proj_order"] = res["analysis"].apply(lambda x: proj_order.index(x))
        res["bio_len"] = res.bio_repl.str.len()
        res["tech_len"] = res.tech_repl.str.len()
        df_sort = res.sort_values(["proj_order", "bio_len", "tech_len"], ascending=True)
        df_dup = df_sort.drop_duplicates(subset=["chrom", "start", "end", "assembly"], keep="last")
        return df_dup.drop(["proj_order", "bio_len", "tech_len"], axis=1)

    def get_peaks_per_feature(meta_df, feature_list, feature_exp_dict, path_bed_out, proj_order):
        i = 0
        read_error_log = logging.get_logger('read_error_log')
        read_error_log.setLevel(logging.ERROR)
        handler = logging.FileHandler( os.path.join(path_bed_out, "peak_size_read_errors.txt"), "w+")
        formatter = logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        read_error_log.addHandler(handler)

        for feature in feature_list:
            """
            Aggregate each feature into a df and sort/clean it, then write it to the bed file and remove the dataframe from memory
            """
            feature_df = None
            if feature_exp_dict[feature] is None:
                continue
            for exp in feature_exp_dict[feature]:  # Iterate through all experiments for 1 specific feature
                # Iterate through all files in this experiment
                for _, row in meta_df[meta_df["Experiment accession"] == exp].iterrows():
                            try:
                                tmp = self.read_bed_file(path_bed_out, row, read_error_log, proj_order)
                            except FileNotFoundError:
                                read_error_log.error(f"FNF: {row['File accession']}:{feature}")
                                continue
                            if tmp is not 0:
                                
                                tmp["feature"] = feature
                                tmp["peak_size"] = 0
                                for chr in tmp.chrom.unique():
                                    tmp.loc[tmp["chrom"]==chr, "peak_size"] = tmp[tmp["chrom"]==chr]["end"] - tmp[tmp["chrom"]==chr]["start"]
                                tmp = tmp[["chrom","start","end","peak_size"]]
                                
                                if feature_df is not None:
                                    feature_df = feature_df.append(tmp, ignore_index=True)
                                else:
                                    feature_df = tmp
                                i += 1 
                                if i%100 is 0:
                                    print(f"reached experiment iteration: {i}")

            if feature_df is not None:
                full_directory_path = os.path.join(path_bed_out, "analyses", '-'.join(feature.split('|')))
                os.makedirs(full_directory_path, exist_ok=True)
                feature_df.to_csv(os.path.join(full_directory_path,f"peak_sizes.tsv"), "\t", index=False, mode="a+", header=False)
                
            del(feature_df)



    def load_peak_file(path:str):
        df = pd.read_csv(path, sep="\t", names=["chrom", "start", "end", "peak_size"])
        
        df["peaks"] = 1

        # at each position, accumulate how many peaks there are, this is not part of loading anymore and should be its own function
        

        # A plot made of a subplot for heatmap plus line plot
        for chr in df.chrom.unique():
            logger.info(f"Working on chromosome {chr}")
            df_chr = df[df.chrom == chr]
            df_chr["positions"] = df_chr.apply(lambda row: list(range(row["start"], row["end"])), axis=1)
        
            df_exp = df_chr.explode("positions")
            df_gb_chr = df_exp[["chrom", "positions", "peaks"]].groupby(by=["chrom", "positions"]).sum()
        
            
            extent = [df_gb_chr.peaks.min(), df_gb_chr.peaks.max(), 0, 1]
            fig, (ax,ax2) = plt.subplots(nrows=2, sharex=True)
            
            ax.imshow(np.expand_dims(df_gb_chr.peaks, axis=0), cmap="plasma", aspect="auto", extent=extent)
            ax.set_yticks([])
            ax.set_xlim(extent[0], extent[1])

            ax2.plot(df_gb_chr.positions,df_gb_chr.peaks)
            plt.xticks(np.arange(df_gb_chr.positions.min(), df_gb_chr.positions.max()+1))

            plt.tight_layout()
            plt.show()
            break

