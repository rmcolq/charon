#!/usr/bin/env python3
"""Score `charon dehost` per-read output against a truth file.

score_dehost.py --truth ont_60k.truth.tsv --calls run.tsv --run E0_ont [--host-label human]
                [--table results.tsv]

Truth labels other than --host-label are treated as microbial. Charon's output lines are
"C|U <tab> read_id <tab> call <tab> length <tab> num_hashes ...". Reads absent from the
output (filtered by length/quality) are reported as `filtered`.
Prints a readable summary and, with --table, appends one row per (run, stratum).
"""
import argparse
import os
from collections import Counter, defaultdict

STRATA = [(0, 1000, "<1k"), (1000, 5000, "1-5k"), (5000, 20000, "5-20k"), (20000, 10**12, ">20k")]


def pair_id(read_id):
    return read_id[:-2] if read_id[-2:] in ("/1", "/2") else read_id


def stratum(length):
    for lo, hi, name in STRATA:
        if lo <= length < hi:
            return name
    return "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--host-label", default="human")
    ap.add_argument("--host-category", default="human", help="category name of host in the index")
    ap.add_argument("--table")
    ap.add_argument("--no-strata", action="store_true")
    args = ap.parse_args()

    truth = {}
    with open(args.truth) as f:
        for line in f:
            rid, label = line.rstrip("\n").split("\t")[:2]
            truth[pair_id(rid)] = label

    # counts[stratum][(truth_label, outcome)]
    counts = defaultdict(Counter)
    seen = set()
    with open(args.calls) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 or parts[0] not in ("C", "U"):
                continue
            rid = pair_id(parts[1])
            if rid not in truth or rid in seen:
                continue
            seen.add(rid)
            call = parts[2]
            outcome = "unclassified" if parts[0] == "U" or call == "" else (
                "host" if call == args.host_category else "microbial")
            s = "all" if args.no_strata else stratum(int(parts[3]))
            for st in {s, "all"}:
                counts[st][(truth[rid], outcome)] += 1
    for rid, label in truth.items():
        if rid not in seen:
            counts["all"][(label, "filtered")] += 1

    labels = sorted({l for l in truth.values()}, key=lambda l: (l != args.host_label, l))
    outcomes = ["host", "microbial", "unclassified", "filtered"]
    rows = []
    order = ["all"] + [name for _, _, name in STRATA if name in counts]
    for st in order:
        c = counts[st]
        host_total = sum(c[(args.host_label, o)] for o in outcomes)
        micro = [l for l in labels if l != args.host_label]
        micro_total = sum(c[(l, o)] for l in micro for o in outcomes)
        if host_total + micro_total == 0:
            continue
        row = {
            "run": args.run, "stratum": st,
            "host_reads": host_total,
            # safety metric: host reads that would end up in the non-host output
            "host_leak_to_microbial": c[(args.host_label, "microbial")] / host_total if host_total else float("nan"),
            "host_recall": c[(args.host_label, "host")] / host_total if host_total else float("nan"),
            "microbial_reads": micro_total,
            "microbial_called_host": sum(c[(l, "host")] for l in micro) / micro_total if micro_total else float("nan"),
            "microbial_recall": sum(c[(l, "microbial")] for l in micro) / micro_total if micro_total else float("nan"),
        }
        for l in labels:
            tot = sum(c[(l, o)] for o in outcomes)
            for o in outcomes:
                row[f"{l}_{o}"] = c[(l, o)] / tot if tot else float("nan")
        rows.append(row)

    for r in rows:
        print(f"[{r['run']}] {r['stratum']:>6}  host n={r['host_reads']:>6} recall={r['host_recall']:.4f} "
              f"leak->microbial={r['host_leak_to_microbial']:.5f} | microbial n={r['microbial_reads']:>6} "
              f"recall={r['microbial_recall']:.4f} called_host={r['microbial_called_host']:.5f}")
    if rows:
        print("  per label (host / microbial / unclassified / filtered):")
        r = rows[0]
        for l in labels:
            print(f"    {l:>8}: " + " / ".join(f"{r[f'{l}_{o}']:.4f}" for o in outcomes))

    if args.table:
        cols = list(rows[0].keys())
        new = not os.path.exists(args.table)
        with open(args.table, "a") as t:
            if new:
                t.write("\t".join(cols) + "\n")
            for r in rows:
                t.write("\t".join(str(r.get(k, "")) for k in cols) + "\n")


if __name__ == "__main__":
    main()
