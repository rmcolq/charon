#!/usr/bin/env python3
"""Stratified, shuffled subsample of a read set with a read->label truth file.

Single-end:  subsample_truth.py --truth read2truth.ont.tsv --reads all.ont.fq.gz \
                 --per-class 20000 --out-prefix ont_60k
Paired-end:  subsample_truth.py --truth read2truth.illumina.tsv \
                 --reads all_1.illumina.fq.gz --reads2 all_2.illumina.fq.gz \
                 --per-class 100000 --out-prefix ilmn_300k

Writes <prefix>.fq.gz (or <prefix>_1/_2.fq.gz) in shuffled order and <prefix>.truth.tsv.
Paired truth ids ending in /1 or /2 are collapsed to the pair id.
"""
import argparse
import gzip
import random
import sys
from collections import defaultdict


def pair_id(read_id):
    return read_id[:-2] if read_id[-2:] in ("/1", "/2") else read_id


def read_fastq(path):
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            yield header, seq, plus, qual


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--reads2")
    ap.add_argument("--per-class", type=int, required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    # Reservoir sample per label so memory scales with --per-class, not the truth file.
    rng = random.Random(args.seed)
    reservoirs = defaultdict(list)
    seen = defaultdict(int)
    with open(args.truth) as f:
        for line in f:
            rid, label = line.rstrip("\n").split("\t")[:2]
            if args.reads2:
                if rid.endswith("/2"):
                    continue
                rid = pair_id(rid)
            seen[label] += 1
            res = reservoirs[label]
            if len(res) < args.per_class:
                res.append(rid)
            else:
                j = rng.randrange(seen[label])
                if j < args.per_class:
                    res[j] = rid
    chosen = {}
    for label in sorted(reservoirs):
        for rid in reservoirs[label]:
            chosen[rid] = label
        print(f"{label}: {len(reservoirs[label])} of {seen[label]}", file=sys.stderr)
    del reservoirs

    def key(header):
        rid = header[1:].split()[0]
        return pair_id(rid) if args.reads2 else rid

    kept = []
    if args.reads2:
        for r1, r2 in zip(read_fastq(args.reads), read_fastq(args.reads2)):
            k = key(r1[0])
            if k in chosen:
                assert key(r2[0]) == k, f"pair mismatch {r1[0]} {r2[0]}"
                kept.append((k, "".join(r1), "".join(r2)))
    else:
        for r in read_fastq(args.reads):
            k = key(r[0])
            if k in chosen:
                kept.append((k, "".join(r), None))

    rng.shuffle(kept)
    if args.reads2:
        o1 = gzip.open(args.out_prefix + "_1.fq.gz", "wt", compresslevel=1)
        o2 = gzip.open(args.out_prefix + "_2.fq.gz", "wt", compresslevel=1)
    else:
        o1 = gzip.open(args.out_prefix + ".fq.gz", "wt", compresslevel=1)
        o2 = None
    with open(args.out_prefix + ".truth.tsv", "w") as t:
        for k, a, b in kept:
            o1.write(a)
            if o2:
                o2.write(b)
            t.write(f"{k}\t{chosen[k]}\n")
    o1.close()
    if o2:
        o2.close()
    print(f"wrote {len(kept)} of {len(chosen)} chosen", file=sys.stderr)


if __name__ == "__main__":
    main()
