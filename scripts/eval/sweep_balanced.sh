#!/usr/bin/env bash
# Build balanced indexes over a grid of (w, k, fpr) and evaluate each on the Hall subsets.
#
# sweep_balanced.sh <charon_binary> <index.tab> <index_dir> <run_dir> <subset_dir> <datasets> <host_fpr> <w,k,fpr> ...
#   datasets: "ont", "ilmn" or "ont,ilmn"
# Appends one line per index to <run_dir>/indexes.tsv (name, technical bins, bin bits, IBF GB, build s, build peak GB).
set -euo pipefail
bin=$1; tab=$2; idx_dir=$3; run_dir=$4; subsets=$5; datasets=$6; host_fpr=$7; shift 7
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$idx_dir" "$run_dir"
for spec in "$@"; do
    IFS=, read -r w k fpr <<< "$spec"
    name=bal_w${w}k${k}_f${fpr}_h${host_fpr}
    threads=8
    (( w <= 27 )) && threads=4
    (( w <= 23 )) && threads=3
    if [[ ! -f $idx_dir/$name.idx ]]; then
        /usr/bin/time -l "$bin" index -t $threads -w "$w" -k "$k" --fpr "$fpr" --balance --category-fpr "human=$host_fpr" \
            -p "$idx_dir/$name" --temp "$idx_dir/tmp_$name" --log "$idx_dir/$name.index.log" "$tab" \
            > "$idx_dir/$name.index.time" 2>&1
    fi
    ibf=$(grep -h "Create new IBF" "$idx_dir/$name.index.log" | sed -E 's/.*with ([0-9]+) technical bins of ([0-9]+) bits.*/\1 \2/')
    read -r tb bits <<< "$ibf"
    gb=$(awk -v t="$tb" -v b="$bits" 'BEGIN{printf "%.2f", t*b/8e9}')
    wall=$(awk '/real/{print $1; exit}' "$idx_dir/$name.index.time")
    peak=$(awk '/peak memory footprint/{printf "%.1f", $1/2^30}' "$idx_dir/$name.index.time")
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$name" "$tb" "$bits" "$gb" "$wall" "$peak" >> "$run_dir/indexes.tsv"
    echo "[$name] technical_bins=$tb ibf=${gb}GB build=${wall}s peak=${peak}GB"
    if [[ $datasets == *ont* ]]; then
        "$here/run_dehost.sh" "${name}_ont" "$bin" "$idx_dir/$name.idx" "$subsets/ont_60k.truth.tsv" "$run_dir" \
            "$subsets/ont_60k.fq.gz" | head -2
    fi
    if [[ $datasets == *ilmn* ]]; then
        "$here/run_dehost.sh" "${name}_ilmn" "$bin" "$idx_dir/$name.idx" "$subsets/ilmn_300k.truth.tsv" "$run_dir" \
            "$subsets/ilmn_300k_1.fq.gz" "$subsets/ilmn_300k_2.fq.gz" -- --min_length 50 | head -2
    fi
done
