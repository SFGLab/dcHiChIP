#!/usr/bin/env python
import pandas as pd
from scipy.stats import pearsonr
import argparse


def inverse_eigenvectors(comp_path, path_out):
    """
    comp_path: path to the .bed file created from Juicer output
    path_out: output .bed file path
    """
    df = pd.read_csv(comp_path, sep="\t")
    df = df[["#1_usercol", "2_usercol", "3_usercol", "4_usercol", "6_pct_gc"]]
    df.columns = ["chr", "start", "end", "eig", "gc"]
    chroms = df["chr"].unique()

    dfs_chrom = []
    for chrom in chroms:
        df_temp = df[df["chr"] == chrom]
        corr, _ = pearsonr(df_temp["eig"], df_temp["gc"])
        if corr < 0:
            df_temp["eig"] = -1 * df_temp["eig"]
        df_temp["comp"] = ['A' if v > 0 else 'B' for v in df_temp["eig"]]
        df_temp["color"] = ['#0000FF' if v == 'B' else '#FF0000' for v in df_temp["comp"]]
        df_temp["other"] = ["."] * df_temp.shape[0]
        df_temp["start2"] = df_temp["start"]
        df_temp["end2"] = df_temp["end"]
        df_temp = df_temp[["chr", "start", "end", "comp", "eig", "other", "start2", "end2", "color"]]
        dfs_chrom.append(df_temp)

    df_final = pd.concat(dfs_chrom)
    df_final.to_csv(path_out, index=False, header=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Invert eigenvectors based on GC correlation and assign compartments.")
    parser.add_argument("--path_comp", required=True, help="Path to input compartments BED file")
    parser.add_argument("--path_out", required=True, help="Path to output BED file")

    args = parser.parse_args()

    inverse_eigenvectors(args.path_comp, args.path_out)
