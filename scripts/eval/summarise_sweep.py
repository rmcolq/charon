#!/usr/bin/env python3
"""Summarise a sweep_balanced.sh run directory: index size/build cost joined with accuracy.

summarise_sweep.py <run_dir> <subset_dir> [--confidence 7,5,3]

For Illumina runs, also re-calls each run offline at other --confidence thresholds (call_host rule;
valid because thresholds don't affect model training).
"""
import argparse
import csv
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("run_dir")
ap.add_argument("subset_dir")
ap.add_argument("--confidence", default="7,5,3")
args = ap.parse_args()

indexes = {}
with open(os.path.join(args.run_dir, "indexes.tsv")) as f:
    for row in f:
        name, tb, bits, gb, wall, peak = row.rstrip("\n").split("\t")
        indexes[name] = dict(tb=tb, gb=gb, build=wall, peak=peak)

results = {}
with open(os.path.join(args.run_dir, "results.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["stratum"] == "all":
            results[row["run"]] = row

resources = {}
with open(os.path.join(args.run_dir, "resources.tsv")) as f:
    for row in f:
        run, wall, rss = row.rstrip("\n").split("\t")
        resources[run] = (wall, rss)

print("index\tdata\ttech_bins\tibf_GB\tbuild_s\tbuild_peak_GB\tclassify_s\tclassify_peak_GB\t"
      "human_host\thuman_to_microbial\tmtb\tzymo\tmicrobial_to_host")
for name, ix in indexes.items():
    for data in ("ont", "ilmn"):
        run = f"{name}_{data}"
        if run not in results:
            continue
        r = results[run]
        wall, rss = resources.get(run, ("", ""))
        print("\t".join([name, data, ix["tb"], ix["gb"], ix["build"], ix["peak"], wall, rss,
                         f"{float(r['host_recall']):.4f}", f"{float(r['host_leak_to_microbial']):.5f}",
                         f"{float(r['mtb_microbial']):.4f}", f"{float(r['zymo_microbial']):.4f}",
                         f"{float(r['microbial_called_host']):.5f}"]))

print("\noffline confidence sweep (Illumina, call_host rule, host_lo 0.05):")
for name in indexes:
    calls = os.path.join(args.run_dir, f"{name}_ilmn.calls.tsv")
    if not os.path.exists(calls):
        continue
    out = subprocess.run([sys.executable, "-I", os.path.join(here, "sweep_host_thresholds.py"), "--calls", calls,
                          "--truth", os.path.join(args.subset_dir, "ilmn_300k.truth.tsv"), "--min-len", "50",
                          "--rule", "host", "--confidence", args.confidence, "--min-prop-diff", "0.04",
                          "--host-lo", "0.05"], capture_output=True, text=True, check=True).stdout
    for line in out.splitlines()[1:]:
        print(f"{name}\t{line}")
