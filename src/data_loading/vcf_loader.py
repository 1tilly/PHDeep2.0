import os # For path building of paths and checking of file availability
import subprocess
import sys # to stream output of bcftools directly to file
import pandas as pd # to handle the sample ids file and saving the 
from rpy2 import robjects




class BCFParser():

    vcf_format = "%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AC_PAH\t%AN_PAH\t%AC_UPAH\t%AN_UPAH\t%AC_CTRLS\t%AN_CTRLS\t%AC_UCTRLS\t%AN_UCTRLS\t%AC_WGS13K\t%AN_WGS13K\t%AC_UWGS13K\t%AN_UWGS13K\t%NS\t%AF\t%MAF\t%HWE\t"\
    "%AF_PAH\t%MAF_PAH\t%HWE_PAH\t%AF_UPAH\t"\
    "[%GT\t]\n"
        
    def load_rpy2(self, path_to_script):
        os.environ['R_HOME'] = "/home/tt419/.conda/envs/r_env/bin/R"
        r = robjects.r
        r['source'](path_to_script)
        return r

    def load_sample_info_rpy2(self, output_path, force=False):
        persist_cohort_info_function_r = robjects.globalenv['persist_cohort_info']
        path_sample_ftr = os.path.join(output_path, "sample_info.feather")
        path_oc_ftr = os.path.join(output_path, "oc_ids.feather")
        path_samplePhenoFtr = os.path.join(output_path, "sample_pheno.feather")
        path_phenoCohort_ftr = os.path.join(output_path, "pheno_cohort.feather")
        if os.path.isfile(path_sample_ftr) and os.path.isfile(path_oc_ftr) and os.path.isfile(path_samplePhenoFtr) and os.path.isfile(path_phenoCohort_ftr):
           if force:
               persist_cohort_info_function_r(output_path)
        else:
           persist_cohort_info_function_r(output_path)
        samples = pd.read_feather(path_sample_ftr)
        oc_idmap = pd.read_feather(path_oc_ftr)
        sample_pheno = pd.read_feather(path_samplePhenoFtr)
        pheno_cohort = pd.read_feather(path_phenoCohort_ftr)
        return samples, oc_idmap, sample_pheno, pheno_cohort

    def load_sample_info_feather(self, path_sampleInfoFtr):
        sample_info = pd.read_feather(path_sampleInfoFtr)
        return sample_info

    def get_sample_IDs(self, sample_info, id_out_path, label_out_path, cohort="PAH"):
        pah_samples = sample_info.loc[sample_info["group"].isin([cohort])]
        pah_sample_ids = pah_samples["WGS.ID"].to_list()
        pah_samples["WGS.ID"].to_csv(id_out_path, sep="\t", index=False)
        
        sample_info[["WGS.ID", "labels"]].to_csv(label_out_path, sep="\t", index=False)
        return pah_sample_ids

    def load_roi_rpy2(self, output_path, genelist="PAH_disease_genes_prev_reported_and_novel", grch=37):
        """
        ROI := region of interests
        """
        persist_genelist_function_r = robjects.globalenv['persist_genelist']
        persist_genelist_function_r(output_path, genelist, grch)
        ensg = pd.read_feather(output_path)
        print(f"There are {len(ensg)} regions of interest.")
        return ensg

    def load_roi_feather(self, path_ensg37ftr):
        """
        ROI := region of interests
        """
        return pd.read_feather(path_ensg37ftr)

    def iterate_through_roi(self, df_ensg, path_release, path_vcfs, suffix_for_vcfs, specific_roi:list=[], columns=["symbol", "chr", "start","end"]):
        """
        ROI := region of interests
        """
        symbols = specific_roi or df_ensg[columns[0]].unique()
        for gene_symbol in symbols:
            print(f"Running for gene {gene_symbol}")
            gene = df_ensg.loc[df_ensg[columns[0]].isin([gene_symbol])]
            gene_chr = gene[columns[1]].to_list()[0]
            gene_start = gene[columns[2]].to_list()[0]
            gene_end = gene[columns[3]].to_list()[0]
            gene_name = gene[columns[0]].to_list()[0]
            path_to_vcf = os.path.join(path_release, path_vcfs, f"chr{gene_chr}{suffix_for_vcfs}")
            yield gene_chr, gene_start, gene_end, gene_name, path_to_vcf


    
    def run_bcftools(self, path_to_vcf, sample_label_path, out_path, release, out_format, gene_name, gene_chr, gene_start, gene_end, pah_sample_ids=None):
    #   Generate BCFtools command
        """
        This code filters and annotates a BCF file for a specific set of samples 
            and a genomic region, then queries the modified BCF file to extract data 
            in a specified format.
        The command is as follows:
        "bcftools view -Ou -S ", out.dir, "/samples_included.txt ",
                        bcf.file, " ", gsub("^chr", "", roi$chr), ":", roi$start, "-", roi$end,
                        " | bcftools plugin fill-tags ",
                        " -- -S ", paste0(analysis.dir, "/sample_labels_included.txt"),
                        " | bcftools query -f \"", format, "\""

        Following, the command is broken down into its components, combined with the right information from the python environment and subsequently run.
        """

        gene_coords =  f"{gene_chr}:{gene_start}-{gene_end}"
        command_view1 = f"~/bin/bcftools/bcftools view -IOu {path_to_vcf} {gene_coords}" # -I prevents update of AC/AN/AF
        command_plugin = f"~/bin/bcftools/bcftools plugin fill-tags -Ou -- -S {sample_label_path} -t AC,AN,AF,MAF,HWE,NS "
 
        command_query = f"~/bin/bcftools/bcftools query -H -f \"{out_format}\" "
        vcf_out_path = os.path.join(out_path, f'pah_variants_{release}_{gene_name}_{gene_chr}_{gene_start}-{gene_end}.vcf')
        full_cmd = f"{command_view1} | {command_plugin}  | {command_query} > {vcf_out_path}"
        
        if pah_sample_ids is not None:  
            command_view_filter = f"~/bin/bcftools/bcftools view -IOu -s {','.join(pah_sample_ids)} "
            full_cmd = f"{command_view1} | {command_plugin} | {command_view_filter} | {command_query} > {vcf_out_path}"



        # Execute BCFtools
        with open(vcf_out_path, 'wb') as f:
            process = subprocess.Popen(full_cmd, stdout=subprocess.PIPE, shell=True)
            for c in iter(lambda: process.stdout.read(1), b''):
                sys.stdout.buffer.write(c)
                f.buffer.write(c)

        return vcf_out_path


    def persist_to_feather(self, vcf_path, pah_IDs, var_output_path): 
        """
        To match the structure of the AVRO files, the columns should be:
        [u'chromosome', u'start', u'end', u'reference', u'alternate',u'uctrl_SD', u'upah_SD', u'uctrl_AC', u'uctrl_AN', u'upah_AC',u'upah_AN']
        """

        with open(vcf_path) as vcf:
            vcf = vcf.read().split("\n")
            if len(vcf)>1:
                header = filter(lambda l: l.startswith("#"), vcf)
                head = list(filter(lambda l: not l.startswith("##"), header))[0].split("\t")
                head = list(map(lambda x: x.split("]")[1].split(":")[0], head[:-1]))
            else:
                print(f"No entries for {vcf_path}, writing ommitted")
                return None

        vcf = pd.read_csv(vcf_path, sep="\t", comment="#", names=head, index_col=False)
        if len(vcf)>0:
            vcf.loc[:,[x for x in vcf.columns if "AF" in x]] = vcf.loc[:,[x for x in vcf.columns if "AF" in x]].replace({".":0.0})
            vcf = vcf.astype({x:"int64" for x in vcf.columns if "AC" in x})
            vcf = vcf.astype({x:"float64" for x in vcf.columns if "AF" in x})
            vcf[list(vcf.filter(regex="[A-Z]*\d{6}").columns)] = vcf[list(vcf.filter(regex="[A-Z]*\d{6}").columns)].astype(str)
            vcf["pah_var"] = vcf[pah_IDs].apply(lambda col: col.str.contains("1"), axis=0).any(axis=1)


            vcf["length"] = vcf.REF.str.len()
            vcf["end"] = vcf[["POS","length"]].apply(lambda x: x[0]+(x[1]-1) if x[1]>1 else x[0], axis=1)

            name_dict= {"POS":"start", "CHROM":"chromosome", 
                    "REF":"reference", "ALT":"alternate"
                    }
            
            def rename_acnf(col):
                if "AC_" in col or "AN_" in col or "AF_" in col:
                    al,pop = col.split("_")
                    return f"{pop.lower()}_{al}"
                else:
                    return col

            vcf = vcf.drop(["length"], axis=1).rename(name_dict, axis=1).rename(rename_acnf, axis=1)
            vcf.to_feather(var_output_path)
        else:
            print(f"No entries for {vcf_path}, writing ommitted")
            return None
        return vcf
