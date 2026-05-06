#!/usr/bin/env python3
"""
Annotate one or more TSV variant files with a boolean column indicating
whether each variant's genomic position falls within any region in a
given BED file.

The new column is inserted immediately BEFORE the last column (is_correct),
so is_correct remains the final column in the output.

Usage examples
--------------
# Single file, output written next to input with _annotated suffix:
    python annotate_lcr.py regions.bed variants.tsv

# Multiple files using a glob:
    python annotate_lcr.py regions.bed *.tsv

# Write all outputs to a specific directory:
    python annotate_lcr.py regions.bed *.tsv -o annotated/

# Custom column name:
    python annotate_lcr.py regions.bed variants.tsv --col-name in_target
"""

import os
import sys
import glob
import bisect
import argparse

import pandas as pd


# ---------------------------------------------------------------------------
# BED loading & interval lookup
# ---------------------------------------------------------------------------

def load_bed(bed_file: str) -> dict:
    """
    Read a BED file and return a dict of merged, sorted intervals per chromosome.

    BED coordinates are 0-based half-open [start, end).
    Overlapping/adjacent intervals are merged so binary search works correctly.

    Returns
    -------
    dict : {chrom: [(start, end), ...]}  sorted and merged per chromosome
    """
    raw: dict = {}
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t") if "\t" in line else line.split()
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            raw.setdefault(chrom, []).append((start, end))

    # Merge overlapping / adjacent intervals per chromosome
    merged: dict = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        stack = [list(ivs[0])]
        for start, end in ivs[1:]:
            if start <= stack[-1][1]:          # overlapping or adjacent
                stack[-1][1] = max(stack[-1][1], end)
            else:
                stack.append([start, end])
        merged[chrom] = [(s, e) for s, e in stack]

    return merged


def pos_in_bed(chrom: str, pos: int, regions: dict) -> bool:
    """
    Return True if the 1-based position *pos* on *chrom* falls inside any
    BED region.

    Coordinate arithmetic
    ---------------------
    BED interval [start, end) covers 1-based positions start+1 … end.
    Therefore pos (1-based) is inside [start, end) iff start < pos <= end.

    The lookup uses binary search on the pre-sorted, merged interval list,
    giving O(log n) per query.
    """
    if chrom not in regions:
        return False

    intervals = regions[chrom]
    starts = [iv[0] for iv in intervals]

    # Last index whose start is strictly less than pos
    idx = bisect.bisect_left(starts, pos) - 1
    if idx < 0:
        return False

    # Interval at idx has start < pos; check whether end covers pos
    return intervals[idx][1] >= pos


# ---------------------------------------------------------------------------
# Variant-key parsing
# ---------------------------------------------------------------------------

def parse_chrom_pos(variant_key: str):
    """
    Extract (chrom, 1-based_position) from a variant_key string.

    The key format is always:  chrom:pos:ref:alt
    e.g.  chr1:9908170:A:]chr6:32366795]A
    """
    parts = variant_key.split(":")
    return parts[0], int(parts[1])


# ---------------------------------------------------------------------------
# Per-file annotation
# ---------------------------------------------------------------------------

def annotate_tsv(tsv_path: str, regions: dict, out_path: str,
                 col_name: str = "in_bed") -> None:
    """
    Read *tsv_path*, add a boolean *col_name* column immediately before the
    last column (is_correct), and write the result to *out_path*.
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)

    if "variant_key" not in df.columns:
        raise ValueError(f"'variant_key' column not found in {tsv_path}")

    results = []
    for vk in df["variant_key"]:
        try:
            chrom, pos = parse_chrom_pos(vk)
            results.append(pos_in_bed(chrom, pos, regions))
        except (IndexError, ValueError):
            results.append(False)   # malformed key → False

    # Insert new column BEFORE the last column so is_correct stays last
    insert_idx = len(df.columns) - 1          # position of current last col
    df.insert(insert_idx, col_name, results)

    df.to_csv(out_path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Add an in_bed boolean column to TSV variant files based on "
            "overlap with a BED file. The new column is placed immediately "
            "before the last column (is_correct)."
        )
    )
    parser.add_argument("bed",
                        help="BED file defining regions of interest")
    parser.add_argument("tsv", nargs="+",
                        help="TSV file(s) to annotate; glob patterns are supported")
    parser.add_argument("-o", "--output-dir", default=None,
                        help=(
                            "Directory for output files. "
                            "If omitted, each output is written next to its "
                            "input with an '_annotated' suffix added before "
                            "the extension."
                        ))
    parser.add_argument("--col-name", default="in_bed",
                        help="Name of the new boolean column (default: in_bed)")
    args = parser.parse_args()

    # Expand any glob patterns supplied on the command line
    tsv_files = []
    for pattern in args.tsv:
        matches = glob.glob(pattern)
        tsv_files.extend(matches if matches else [pattern])

    if not tsv_files:
        sys.exit("Error: no TSV files found.")

    # Load BED once
    print(f"Loading BED file: {args.bed}")
    regions = load_bed(args.bed)
    n_regions = sum(len(v) for v in regions.values())
    print(f"  → {n_regions} merged region(s) across {len(regions)} chromosome(s)")

    # Prepare output directory if requested
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)

    # Annotate each TSV
    for tsv in tsv_files:
        if args.output_dir:
            out = os.path.join(args.output_dir, os.path.basename(tsv))
        else:
            base, ext = os.path.splitext(tsv)
            out = f"{base}_annotated{ext}"

        print(f"  {tsv}  →  {out}")
        try:
            annotate_tsv(tsv, regions, out_path=out, col_name=args.col_name)
        except Exception as exc:
            print(f"    ERROR: {exc}", file=sys.stderr)

    print("Done.")


if __name__ == "__main__":
    main()