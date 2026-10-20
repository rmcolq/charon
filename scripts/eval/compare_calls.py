#!/usr/bin/env python3
"""Compare two `charon dehost` outputs read by read.

compare_calls.py a.calls.tsv b.calls.tsv

"Evidence" = length, num_hashes, mean quality, confidence and per-category count/proportion/unique
proportion. These don't depend on the order reads reach the model, so they must match exactly
when only performance changed. Calls and probabilities can differ when multithreading changes
which reads train the model; those are reported separately.
"""
import sys
from collections import Counter


def load(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in ("C", "U"):
                continue
            cats = []
            for field in p[8].split():
                name, count, prop, uprop, prob = field.split(":")
                cats.append((name, count, prop, uprop))
            evidence = (p[3], p[4], p[5], p[6], tuple(cats))
            out[p[1]] = (evidence, p[0] + ":" + p[2])
    return out


a, b = load(sys.argv[1]), load(sys.argv[2])
common = a.keys() & b.keys()
ev_diff = [r for r in common if a[r][0] != b[r][0]]
call_diff = Counter((a[r][1], b[r][1]) for r in common if a[r][1] != b[r][1])
print(f"reads: a={len(a)} b={len(b)} common={len(common)} only_a={len(a.keys() - b.keys())} only_b={len(b.keys() - a.keys())}")
print(f"evidence differs: {len(ev_diff)}")
for r in ev_diff[:5]:
    print("  ", r, a[r][0], "|", b[r][0])
print(f"calls differ: {sum(call_diff.values())}")
for k, v in call_diff.most_common(10):
    print(f"   {k[0]:>14} -> {k[1]:<14} {v}")
