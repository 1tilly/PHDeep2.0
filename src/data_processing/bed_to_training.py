import argparse
from pathlib import Path

import pandas as pd

from src.data_loading.bed_loader import BEDParser


class BedTrainingSet():

    def __init__(self, bed_parser=None):
        self.bed_parser = bed_parser or BEDParser()

    def generate_feature_exp_dict(self, assay_df):
        return self.bed_parser.generate_feature_exp_dict(assay_df)

    def generate_feature_list(self, feature_exp_dict):
        return self.bed_parser.generate_feature_list(feature_exp_dict)

    def write_bed_by_feature(self, assay_df, feature_list, feature_exp_dict, in_path, out_path, error_log, assembly="hg19", exp_count=20):
        i = 0
        written_features = []
        with open(error_log, "w+", encoding="utf-8") as read_error_log:
            for feature in feature_list:
                """
                Aggregate each feature into a df and sort/clean it, then write it to the bed file and remove the dataframe from memory
                """
                feature_frames = []
                if feature_exp_dict[feature] is None:
                    continue
                for exp in feature_exp_dict[feature]:  # Iterate through all experiments for 1 specific feature
                    # Iterate through all files in this experiment
                    for _, row in assay_df[assay_df["Experiment accession"] == exp].iterrows():
                        try:
                            tmp = self.bed_parser.read_bed_file(in_path, row, assembly)
                        except FileNotFoundError:
                            read_error_log.write(f"FNF: {row['File accession']}:{feature}\n")
                            continue
                        if isinstance(tmp, pd.DataFrame):
                            tmp_1 = self.bed_parser.choose_assay_version(tmp, assembly)
                            tmp_1["feature"] = feature
                            tmp_1 = tmp_1[["chrom", "start", "end", "feature"]]
                            feature_frames.append(tmp_1)
                            i += 1

                if feature_frames:
                    feature_df = pd.concat(feature_frames, ignore_index=True)
                    feature_df.to_csv(out_path, sep="\t", index=False, mode="a+", header=False)
                    written_features.append(feature)
                if exp_count and i >= exp_count:
                    break

        return written_features


def run_bed_to_training(
    meta_path: str | Path,
    in_bed_path: str | Path,
    out_bed_path: str | Path,
    error_log_path: str | Path,
    feature_out_path: str | Path,
    assembly: str = "hg19",
    organism: str = "Homo sapiens",
) -> dict[str, int | str]:
    meta_path = Path(meta_path)
    in_bed_path = Path(in_bed_path)
    out_bed_path = Path(out_bed_path)
    error_log_path = Path(error_log_path)
    feature_out_path = Path(feature_out_path)

    out_bed_path.parent.mkdir(parents=True, exist_ok=True)
    error_log_path.parent.mkdir(parents=True, exist_ok=True)
    feature_out_path.parent.mkdir(parents=True, exist_ok=True)
    out_bed_path.unlink(missing_ok=True)

    with open(meta_path, "r", encoding="utf-8") as f:
        meta = pd.read_csv(f, sep="\t")

    bed_parser = BEDParser()
    training_builder = BedTrainingSet(bed_parser=bed_parser)

    n_exp = len(
        meta[
            (meta["Biosample organism"] == organism)
            & (meta["Assay"].isin(bed_parser.assays))
        ]["Experiment accession"].unique()
    )

    df = meta[
        (meta["Biosample organism"] == organism)
        & (meta["Assay"].isin(bed_parser.assays))
    ][bed_parser.columns_of_interest].copy()

    assay_df = df[df["File format"].isin(["bed narrowPeak"])].copy()

    assay_df["Experiment target"] = assay_df["Experiment target"].apply(
        lambda x: x.split("-")[0] if isinstance(x, str) else x
    )
    clean_df = bed_parser.clean_meta_df(assay_df, assembly=assembly)
    feature_exp_dict = training_builder.generate_feature_exp_dict(clean_df)
    feature_list = training_builder.generate_feature_list(feature_exp_dict)
    written_features = training_builder.write_bed_by_feature(
        clean_df,
        feature_list,
        feature_exp_dict,
        str(in_bed_path),
        str(out_bed_path),
        str(error_log_path),
        assembly,
        False,
    )

    with open(feature_out_path, "w+", encoding="utf-8") as f:
        f.write("\n".join(feature_list))

    return {
        "planned_features": len(feature_list),
        "written_features": len(written_features),
        "experiments_considered": n_exp,
        "assembly": assembly,
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--meta", help="Metadata file")
    parser.add_argument("-i", "--inBed", help="Path to the input bed files")
    parser.add_argument("-o", "--outBed", help="Path to the written bed file")
    parser.add_argument("-e", "--error", help="Path to error log, usually a read_errors.txt in the output directory")
    parser.add_argument("-f", "--featureOut", help="Path for the features")
    parser.add_argument(
        "-a",
        "--assembly",
        default="hg19",
        choices=["hg19", "GRCh38"],
        help="Genome assembly used in metadata filtering and assay version selection.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    result = run_bed_to_training(
        meta_path=args.meta,
        in_bed_path=args.inBed,
        out_bed_path=args.outBed,
        error_log_path=args.error,
        feature_out_path=args.featureOut,
        assembly=args.assembly,
    )
    print(f"Planned features: {result['planned_features']}")
    print(f"Written features: {result['written_features']}")
    print(f"Experiments considered: {result['experiments_considered']}")
    print(f"Assembly used: {result['assembly']}")
