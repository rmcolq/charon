#!/usr/bin/env python3
"""Re-call `charon dehost` output offline under different thresholds (mirrors ReadEntry::call_host).

sweep_host_thresholds.py --calls run.calls.tsv --truth truth.tsv [--check]
    --confidence 7,4,2 --min-prop-diff 0.04,0.02 --host-lo 0.05,0.1

Only valid when probabilities are fixed for the run (i.e. the model isn't refitted per threshold),
which holds for runs without --extract (default KDE). --check verifies that the defaults reproduce
Charon's own calls exactly before sweeping.
"""
import argparse
import itertools


def pair_id(read_id):
    return read_id[:-2] if read_id[-2:] in ("/1", "/2") else read_id


def load(path, host="human"):
    reads = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in ("C", "U"):
                continue
            cats = {}
            for field in p[8].split():
                name, count, prop, uprop, prob = field.split(":")
                cats[name] = (int(count), float(prop), float(uprop), float(prob))
            other = next(c for c in cats if c != host)
            reads.append(dict(id=pair_id(p[1]), call=p[2] if p[0] == "C" else "", length=int(p[3]),
                              n=int(p[4]), q=float(p[5]), conf=int(p[6]), comp=float(p[7]),
                              host=cats[host], other=cats[other], other_name=other,
                              cats=cats, order=list(cats)))
    return reads


def call_host(r, conf_thr, mpd, lo, min_q=15.0, min_len=0, min_comp=0.0, min_prob_diff=0.0, conf_prob_thr=0.0,
              host="human"):
    h_up, o_up = r["host"][2], r["other"][2]
    h_p, o_p = r["host"][3], r["other"][3]
    conf = r["conf"]
    if conf < conf_thr or r["q"] < min_q or r["length"] < min_len or r["comp"] < min_comp:
        return ""
    if (h_up > o_up and h_up - o_up > mpd and h_p > o_p and h_p - o_p > min_prob_diff
            and max(h_p * conf, conf) >= conf_prob_thr):
        return host
    if (h_up < lo and h_up < o_up and o_up - h_up > mpd and h_p < o_p and o_p - h_p > min_prob_diff
            and max(o_p * conf, conf) >= conf_prob_thr):
        return r["other_name"]
    return ""


def call_category(r, conf_thr, mpd, min_q=15.0, min_len=0, min_comp=0.0, min_hits=0):
    """Mirrors ReadEntry::call_category (what paired-end dehost currently uses). Category order as printed."""
    cats = [r["cats"][name] for name in r["order"]]
    uc = [round(c[2] * r["n"]) for c in cats]
    first, second = 0, 1
    if uc[second] > uc[first]:
        first, second = second, first
    conf = min(uc[first] - uc[second], 255)
    if r["q"] < min_q or r["length"] < min_len or r["comp"] < min_comp:
        return ""
    call = ""
    probs = [c[3] for c in cats]
    if probs[second] == 0 and probs[first] > 0:
        call = r["order"][first]
    elif conf > conf_thr and probs[first] > probs[second]:
        call = r["order"][first]
    counts = [c[0] for c in cats]
    if counts[second] > counts[first] or counts[first] - counts[second] < min_hits:
        call = ""
    props = [c[1] for c in cats]
    if props[second] > props[first] or props[first] - props[second] < mpd:
        call = ""
    return call


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--confidence", default="7")
    ap.add_argument("--min-prop-diff", default="0.04")
    ap.add_argument("--host-lo", default="0.05")
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--rule", choices=["host", "category"], default="host",
                    help="host = ReadEntry::call_host (single-end dehost); category = call_category (current paired dehost)")
    args = ap.parse_args()

    truth = {}
    with open(args.truth) as f:
        for line in f:
            rid, label = line.rstrip("\n").split("\t")[:2]
            truth[pair_id(rid)] = label
    reads = load(args.calls)

    def decide(r, c, m, lo):
        if args.rule == "category":
            return call_category(r, c, m, min_len=args.min_len)
        return call_host(r, c, m, lo, min_len=args.min_len)

    if args.check:
        bad = [r["id"] for r in reads if decide(r, 7, 0.04, 0.05) != r["call"]]
        print(f"check: {len(reads) - len(bad)}/{len(reads)} reproduce Charon's calls" + (f"; e.g. {bad[:3]}" if bad else ""))

    n_by_label = {}
    for label in truth.values():
        n_by_label[label] = n_by_label.get(label, 0) + 1
    labels = sorted(n_by_label, key=lambda l: (l != "human", l))
    print("conf\tmpd\thost_lo\t" + "\t".join(f"{l}_recall" for l in labels) + "\thuman_leak\tmicro_called_host")
    for c, m, lo in itertools.product([int(x) for x in args.confidence.split(",")],
                                      [float(x) for x in args.min_prop_diff.split(",")],
                                      [float(x) for x in args.host_lo.split(",")]):
        correct = {l: 0 for l in labels}
        leak = micro_host = 0
        for r in reads:
            label = truth.get(r["id"])
            if label is None:
                continue
            call = decide(r, c, m, lo)
            is_host = label == "human"
            if call == ("human" if is_host else r["other_name"]):
                correct[label] += 1
            if is_host and call == r["other_name"]:
                leak += 1
            if not is_host and call == "human":
                micro_host += 1
        n_host = n_by_label.get("human", 0)
        n_micro = sum(v for k, v in n_by_label.items() if k != "human")
        print(f"{c}\t{m}\t{lo}\t" + "\t".join(f"{correct[l] / n_by_label[l]:.4f}" for l in labels)
              + f"\t{leak / n_host:.6f}\t{micro_host / n_micro:.6f}")


if __name__ == "__main__":
    main()
