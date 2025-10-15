#!/usr/bin/env python

import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process and analyze sequencing data using BWA and downstream tools."
    )
    parser.add_argument("--threads", type=int, default=2, help="Number of threads to use")
    parser.add_argument("--mapq", type=int, default=30, help="Minimum mapping quality")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--qc-file", type=str, required=True, help="QC file")
    parser.add_argument("--input-bam", type=str, default=None,required=True, help="Input BAM input")
    parser.add_argument("--command", type=str, default="filter_pair", help="Command to perform")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    print(args)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if args.command == "filter_pair":
        from feather_filter_chr import filter_pair_reads
        filter_pair_reads(args.input_bam, args.mapq, args.outdir + "/filtered.bam", args.qc_file)
