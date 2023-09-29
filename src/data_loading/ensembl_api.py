import requests
import pandas as pd




class EnsemblAPI():
    server = "https://grch37.rest.ensembl.org"

    def __init__(self):
        pass

    def get_reference(self, chrom, start, end=None, species="human", GRCh="GRCh37", seq_len = None):
        if end is None:
            end = start + 5000
        assert start <= end, "Start has to be before end"
        assert end - start < 10000000, "Can't fetch more then 10 Mb (10 million bases)"
        sub_path = f"/sequence/region/{species}/"
        region = f"{chrom}:{int(start)}..{int(end)}"
        r = requests.get(self.server + sub_path + region, headers={"Content-Type": "text/plain", "coord_system_version": str(GRCh)})
        if seq_len is not None:
            assert len(r.text) == seq_len, "Reference is not of length seq_len"
        if r.ok:
            return r.text
        else:
            r.raise_for_status()
            return -1
    
        
    def get_ref_for_chrom_region(self, var_df, frame_length=1000):
        chrom = var_df["chromosome"][0]
        start = var_df["start"].min() - frame_length
        end = var_df["end"].max() + frame_length
        return (chrom, start, end), self.get_reference(chrom, start, end)

        
    def get_ref_for_region(self, chrom, start, end, seq_len):
        start -= seq_len/2 # TODO: needs to be adjusted for length of ref/alt
        end += seq_len/2

        api_cap = 10000000-1
        if end - start > api_cap: # the ensembl API allows maximum 10mb 
            reference = ""
            for i in range(int(start), int(end), api_cap):
                if i+api_cap > end:
                    e = end
                else:
                    e = i+api_cap
                reference += self.get_reference(chrom,i,e, seq_len = seq_len)
        else:
            reference = self.get_reference(chrom, start, end)
        return (chrom, start, end), reference

    def get_vep_annotation(self, chrom, start, end, ref, alt, species="human", GRCh="GRCh37"):
        sub_path = f"/vep/human/region/{chrom}:{start}-{end}/{alt}?"
        r = requests.get(self.server+sub_path, headers={ "Content-Type" : "application/json"})
 
        if not r.ok:
            r.raise_for_status()
            return -1
        
        decoded = r.json()
        return decoded

    def get_vep_multiSnps(self, var_records, species="human", GRCh="GRCh37"):
        sub_path = f"/vep/homo_sapiens/region"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        variants = [f"{x[0]}  {x[1]}  .  {x[3]}  {x[4]} . . ." for x in var_records] 
        # Needed to make sure that double and single quotes are correct for json
        r = requests.post(self.server+sub_path, headers=headers, data='{' + '"variants" : ' + str(variants).replace("\'", "\"")+'}' )
        if not r.ok:
            r.raise_for_status()
            return -1
        
        decoded = r.json()
        return decoded

    def get_vep_annotation_df(self, res):
        df_v = pd.json_normalize(res, record_path=["transcript_consequences"], meta=["seq_region_name","start","end", "allele_string","input","colocated_variants","most_severe_consequence"], errors="ignore")
        return df_v.groupby(["seq_region_name",	"start",	"end",	"allele_string"]).aggregate(lambda x: list(set(x)) if len(set(x)) > 1 else list(x)[0]).reset_index()  
