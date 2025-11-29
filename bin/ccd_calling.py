#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys
import time
import multiprocessing as mp
import numpy as np

# Maximum allowed CCD length (also used to filter very long loops)
min_huge_ccd = 25_000_000  # 25 Mb


def summarize_ccds(chrom, ccd_list, loops_chr):
    """
    Print a real-time summary for each chromosome.
    """
    if not ccd_list:
        print(f"[{chrom}] CCDs: 0 (NO DOMAINS) | loops used={len(loops_chr)}",
              file=sys.stderr)
        return

    lens = [end - start for start, end in ccd_list]
    total_span = sum(lens)
    median_len = np.median(lens)
    min_len = min(lens)
    max_len = max(lens)

    print(
        f"[{chrom}] CCDs: {len(ccd_list)} | "
        f"min={min_len:,} bp | median={int(median_len):,} bp | "
        f"max={max_len:,} bp | total span={total_span/1e6:.2f} Mb | "
        f"loops used={len(loops_chr)}",
        file=sys.stderr
    )


def get_ccds(filename, output_filename, min_loops, min_diff, min_length, summits_only):
    """
    Main CCD-calling routine for MAPS output (single-pass version).

    Expected columns in input (MAPS .bedpe):
      chr1, start1, end1, chr2, start2, end2, ...

    Steps:
      - Reads MAPS bedpe (with header).
      - (Optional) if --summits_only: keep only ClusterSummit == 1 (if present).
      - Keeps intra-chromosomal loops with span < min_huge_ccd.
      - Converts to 1D segments (chr, start1, end2).
      - For each chromosome, builds coverage via event-based model.
      - Runs CCD calling once with given min_loops (coverage-only logic).
      - Writes CCDs >= min_length and < min_huge_ccd to output BED.
      - Prints per-chromosome and global CCD summaries to stderr.
    """
    # Read with header, ignore comment lines if any
    try:
        df = pd.read_csv(filename, sep="\t", header=0, comment="#")
    except Exception as e:
        print(f"Error reading {filename}: {e}", file=sys.stderr)
        sys.exit(1)

    required_cols = {"chr1", "start1", "end1", "chr2", "start2", "end2"}
    missing = required_cols.difference(df.columns)
    if missing:
        print(
            f"Input file is missing required columns: {', '.join(sorted(missing))}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Optional summit filtering
    if summits_only and "ClusterSummit" in df.columns:
        df = df.loc[df["ClusterSummit"] == 1].copy()
        print(
            f"Using only ClusterSummit == 1 interactions (n = {len(df)})",
            file=sys.stderr,
        )
    elif summits_only:
        print(
            "Warning: --summits_only given, but 'ClusterSummit' column "
            "not found. Using all loops instead.",
            file=sys.stderr,
        )

    if df.empty:
        print("No loops remaining after filtering (empty MAPS file?).", file=sys.stderr)
        with open(output_filename, "w"):
            pass
        return

    # Only intra-chromosomal interactions, avoid very large spans
    span = df["end2"] - df["start1"]
    mask_intra = df["chr1"] == df["chr2"]
    mask_span = span < min_huge_ccd

    loops_df = df.loc[mask_intra & mask_span, ["chr1", "start1", "end2"]].copy()
    loops_df.columns = ["chr", "start", "end"]

    # Ensure integer coordinates
    loops_df["start"] = loops_df["start"].astype(int)
    loops_df["end"] = loops_df["end"].astype(int)
    loops_df = loops_df.sort_values(["chr", "start", "end"]).reset_index(drop=True)

    if loops_df.empty:
        print("No loops left after intra-chromosomal/span filtering.", file=sys.stderr)
        with open(output_filename, "w"):
            pass
        return

    # Group loops by chromosome and prepare tasks
    loops_parsed = {}
    for chromosome in loops_df["chr"].unique():
        df_chr = loops_df.loc[loops_df["chr"] == chromosome]
        loops_current_chr = list(zip(df_chr["start"].tolist(),
                                     df_chr["end"].tolist()))
        if loops_current_chr:
            loops_parsed[chromosome] = loops_current_chr

    task_list = [
        (loops_chr, chromosome, min_loops, min_diff)
        for chromosome, loops_chr in loops_parsed.items()
    ]

    threads = min(len(task_list), 24)
    with mp.Pool(threads) as pool:
        results = pool.map(parse_chrom, task_list)

    print(f"Finished CCD calling with min_loops = {min_loops}", file=sys.stderr)

    # Per-chromosome summary
    print("\n=== Per-chromosome CCD summary ===", file=sys.stderr)
    for chrom, ccd_list in results:
        loops_chr = loops_parsed.get(chrom, [])
        summarize_ccds(chrom, ccd_list, loops_chr)

    # Write final CCDs that pass length filters
    all_ccd_lengths = []
    with open(output_filename, "w") as out_f:
        for chrom, ccd_list in results:
            for start, end in ccd_list:
                length = end - start
                if length >= min_length and length < min_huge_ccd:
                    out_f.write(f"{chrom}\t{start}\t{end}\n")
                    all_ccd_lengths.append(length)

    # Global summary
    if all_ccd_lengths:
        print("\n=== GLOBAL CCD SUMMARY ===", file=sys.stderr)
        print(f"Total CCDs (after length filter): {len(all_ccd_lengths)}",
              file=sys.stderr)
        print(f"Median length: {int(np.median(all_ccd_lengths)):,} bp",
              file=sys.stderr)
        print(f"Mean length: {int(np.mean(all_ccd_lengths)):,} bp",
              file=sys.stderr)
        print(f"Max length: {max(all_ccd_lengths):,} bp",
              file=sys.stderr)
    else:
        print("\nNo CCDs passed the length filters genome-wide", file=sys.stderr)


def parse_chrom(args):
    """
    CCD calling for a single chromosome (coverage-only, strict).

    Input:
      - loops_current_chr: list of (start, end) for that chromosome
      - chromosome: chromosome name
      - min_loops: coverage threshold
      - min_diff: (unused here, kept for interface compatibility)

    Logic:
      - Build event map: +1 at start, -1 at end+1
      - Scan sorted coordinates, maintaining coverage
      - OPEN CCD when coverage >= min_loops
      - CLOSE CCD when coverage drops < min_loops
      - Domains are exactly regions with sufficient loop coverage.
        No coverage ⇒ no domain.
    """
    loops_current_chr, chromosome, min_loops, _min_diff_unused = args

    # Build event map: coverage changes at starts/ends
    events = {}
    for s, e in loops_current_chr:
        if s < 0 or e < s:
            continue  # skip invalid intervals
        events[s] = events.get(s, 0) + 1
        # end is inclusive, so coverage drops after e
        events[e + 1] = events.get(e + 1, 0) - 1

    if not events:
        return chromosome, []

    coords = sorted(events.keys())
    ccds_current_chr = []

    cov = 0
    current_ccd_left = None

    for idx, pos in enumerate(coords):
        cov += events[pos]
        next_pos = coords[idx + 1] if idx + 1 < len(coords) else None

        if cov >= min_loops:
            # we are in a high-coverage region
            if current_ccd_left is None:
                current_ccd_left = pos
        else:
            # coverage below threshold: close if something is open
            if current_ccd_left is not None:
                # end at previous position (just before coverage dropped)
                end_pos = (pos - 1) if pos > current_ccd_left else current_ccd_left
                ccds_current_chr.append((current_ccd_left, end_pos))
                current_ccd_left = None

        # if this is the last coordinate and CCD is still open,
        # close it at the last coordinate - 1
        if next_pos is None and current_ccd_left is not None:
            end_pos = (pos - 1) if pos > current_ccd_left else current_ccd_left
            ccds_current_chr.append((current_ccd_left, end_pos))
            current_ccd_left = None

    return chromosome, ccds_current_chr


def main():
    parser = argparse.ArgumentParser(
        prog="CCD caller (MAPS output)",
        description=(
            "Call chromatin contact domains (CCDs) from MAPS significant loop BEDPE.\n\n"
            "Expected MAPS columns: chr1, start1, end1, chr2, start2, end2, "
            "count, expected, fdr, ClusterLabel, ClusterSize, ClusterType, "
            "ClusterNegLog10P, ClusterSummit (header required)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("input_file", help="MAPS BEDPE file with header.")
    parser.add_argument(
        "-o",
        "--output_file",
        default="ccds.bed",
        help="Output BED file (default: ccds.bed).",
    )
    parser.add_argument(
        "--min_loops",
        type=int,
        default=2,
        help="Minimum loop coverage to define CCDs (default: 2).",
    )
    parser.add_argument(
        "--min_diff",
        type=int,
        default=2,
        help="(Unused) kept for compatibility; CCDs are coverage-based only.",
    )
    parser.add_argument(
        "--min_length",
        type=int,
        default=15_000,
        help="Minimum CCD length in bp to report (default: 15000).",
    )
    parser.add_argument(
        "--summits_only",
        action="store_true",
        help="Use only ClusterSummit == 1 interactions if the column exists.",
    )

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print("Input does not exist", file=sys.stderr)
        sys.exit(1)

    start_time_total = time.time()
    get_ccds(
        args.input_file,
        args.output_file,
        args.min_loops,
        args.min_diff,
        args.min_length,
        args.summits_only,
    )
    print(
        f"--- Executed in {time.time() - start_time_total:.2f} seconds ---",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
