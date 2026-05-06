#!/usr/bin/env python3
"""
─────────────
Annotate Manta VCF BND records with:

    ND  = log10( distance in bp to the nearest       other BND on the same chromosome )
    ND2 = log10( distance in bp to the 2nd-nearest   other BND on the same chromosome )

BNDs that are the sole breakend on their chromosome receive ND = ND2 = 6.0
(i.e. log10 of the 1 Mbp cap).  All distances are hard-capped at 1e6 bp
(log10 value capped at 6.0).

Usage:
    python annotate_nd.py input.vcf        output.vcf
    python annotate_nd.py input.vcf.gz     output.vcf.gz
"""

import sys
import math
import gzip
import argparse
from collections import defaultdict


# ── I/O helpers ───────────────────────────────────────────────────────────────

def open_vcf(path: str):
    """Open a plain or gzip-compressed VCF for reading (text mode)."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def open_out(path: str):
    """Open a plain or gzip-compressed file for writing (text mode)."""
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


# ── Helpers ───────────────────────────────────────────────────────────────────

_ND_CAP = 6.0          # log10(1_000_000) — hard ceiling for ND / ND2


def dist_to_log10(d: int) -> float:
    """Convert a bp distance to log10, clamping d<=0 to 0.0 and capping at _ND_CAP."""
    if d <= 0:
        return 0.0
    return min(round(math.log10(d), 4), _ND_CAP)


# ── Core processing ───────────────────────────────────────────────────────────

def process_vcf(input_path: str, output_path: str) -> None:

    header_lines: list[str] = []
    records: list[dict] = []        # one dict per data line, in file order

    # ── 1. Read VCF ──────────────────────────────────────────────────────────
    with open_vcf(input_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
                continue

            fields = line.split("\t")
            records.append({
                "chrom":  fields[0],
                "pos":    int(fields[1]),
                "fields": fields,
                "is_bnd": "SVTYPE=BND" in fields[7],
                "nd":     None,          # filled in step 3
                "nd2":    None,          # filled in step 3
                "solo":   False,         # True when sole BND on its chromosome
            })

    # ── 2. Build a sorted BND position index, grouped by chromosome ──────────
    #       bnd_by_chrom[chrom] = [ (pos, record_index), ... ]  (sorted by pos)
    bnd_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for i, rec in enumerate(records):
        if rec["is_bnd"]:
            bnd_by_chrom[rec["chrom"]].append((rec["pos"], i))

    for lst in bnd_by_chrom.values():
        lst.sort()          # sort ascending by position

    # ── 3. Nearest and 2nd-nearest neighbour distances (same chromosome) ─────
    #
    # Because pos_list is sorted, distances grow monotonically as we move away
    # from j in either direction.  The two smallest distances therefore always
    # come from at most the first two elements to the left AND the first two
    # to the right — four candidates at most.  Sorting those four and picking
    # index [0] / [1] gives ND / ND2 exactly.
    #
    #   n == 1  →  sole BND on chromosome  →  ND = ND2 = _ND_CAP (6.0)
    #   n == 2  →  exactly one candidate   →  ND set,   ND2 = None
    #   n >= 3  →  at least two candidates →  ND set,   ND2 set
    #
    # All values are additionally capped at _ND_CAP by dist_to_log10().
    #
    for chrom, pos_list in bnd_by_chrom.items():
        n = len(pos_list)
        for j, (pos, idx) in enumerate(pos_list):

            # Collect up to 2 neighbours from each side
            dists: list[int] = []
            if j > 0:
                dists.append(pos - pos_list[j - 1][0])       # left  1st
            if j < n - 1:
                dists.append(pos_list[j + 1][0] - pos)       # right 1st
            if j > 1:
                dists.append(pos - pos_list[j - 2][0])       # left  2nd
            if j < n - 2:
                dists.append(pos_list[j + 2][0] - pos)       # right 2nd

            if not dists:
                # Sole BND on this chromosome — assign the cap value to both fields
                records[idx]["nd"]   = _ND_CAP
                records[idx]["nd2"]  = _ND_CAP
                records[idx]["solo"] = True
                continue

            dists.sort()

            records[idx]["nd"] = dist_to_log10(dists[0])
            if len(dists) >= 2:
                records[idx]["nd2"] = dist_to_log10(dists[1])

    # ── 4. Write annotated VCF ────────────────────────────────────────────────
    nd_header = (
        "##INFO=<ID=ND,Number=1,Type=Float,"
        'Description="Log10 of the distance (bp) to the nearest other '
        'BND breakend on the same chromosome; capped at 6.0 (1 Mbp). '
        'Sole BNDs on a chromosome receive the cap value 6.0">'
    )
    nd2_header = (
        "##INFO=<ID=ND2,Number=1,Type=Float,"
        'Description="Log10 of the distance (bp) to the 2nd-nearest other '
        'BND breakend on the same chromosome; capped at 6.0 (1 Mbp). '
        'Sole BNDs on a chromosome receive the cap value 6.0">'
    )

    with open_out(output_path) as out:
        for line in header_lines:
            if line.startswith("#CHROM"):          # insert ND/ND2 meta-lines first
                out.write(nd_header  + "\n")
                out.write(nd2_header + "\n")
            out.write(line + "\n")

        for rec in records:
            fields = rec["fields"]
            if rec["is_bnd"]:
                if rec["nd"]  is not None:
                    fields[7] += f';ND={rec["nd"]}'
                if rec["nd2"] is not None:
                    fields[7] += f';ND2={rec["nd2"]}'
            out.write("\t".join(fields) + "\n")

    # ── 5. Summary to stderr ─────────────────────────────────────────────────
    n_bnd  = sum(1 for r in records if r["is_bnd"])
    n_nd   = sum(1 for r in records if r["is_bnd"] and r["nd"]  is not None)
    n_nd2  = sum(1 for r in records if r["is_bnd"] and r["nd2"] is not None)
    n_solo = sum(1 for r in records if r["is_bnd"] and r["solo"])
    # Solos now carry nd2=cap, so BNDs without nd2 are exclusively "paired" ones
    # (exactly 2 BNDs on their chromosome → ND set, ND2 still None)
    n_pair = n_bnd - n_nd2
    print(
        f"[add_bnd_nd] {n_bnd} BND records found — "
        f"{n_nd} tagged with ND, {n_nd2} tagged with ND2"
        + (f" | {n_solo} assigned cap value ({_ND_CAP}) for ND+ND2 (sole BND on chromosome)" if n_solo else "")
        + (f" | {n_pair} tagged with ND only (paired BNDs on chromosome)" if n_pair else "")
        + ".",
        file=sys.stderr,
    )


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("input",  help="Input VCF or VCF.gz")
    ap.add_argument("output", help="Output VCF or VCF.gz")
    args = ap.parse_args()
    process_vcf(args.input, args.output)


if __name__ == "__main__":
    main()