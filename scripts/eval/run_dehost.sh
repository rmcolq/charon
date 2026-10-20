#!/usr/bin/env bash
# Run `charon dehost` on an evaluation subset, record time/peak RSS, and score against truth.
#
# run_dehost.sh <run_name> <charon_binary> <index.idx> <truth.tsv> <out_dir> <reads1> [reads2] [-- extra dehost args]
set -euo pipefail
run=$1; bin=$2; idx=$3; truth=$4; out=$5; r1=$6; shift 6
r2=""
if [[ $# -gt 0 && $1 != "--" ]]; then r2=$1; shift; fi
[[ $# -gt 0 && $1 == "--" ]] && shift
extra=("$@")
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
threads=${THREADS:-8}

/usr/bin/time -l "$bin" dehost -t "$threads" --db "$idx" --log "$out/$run.log" ${extra[@]+"${extra[@]}"} "$r1" ${r2:+"$r2"} \
    > "$out/$run.calls.tsv" 2> "$out/$run.time"

wall=$(awk '/real/ {print $1; exit}' "$out/$run.time")
rss=$(awk '/maximum resident set size/ {printf "%.2f", $1/1024/1024/1024}' "$out/$run.time")
echo "[$run] wall=${wall}s peak_rss=${rss}GB" | tee "$out/$run.resources"
printf "%s\t%s\t%s\n" "$run" "$wall" "$rss" >> "$out/resources.tsv"

strata=()
[[ -n $r2 ]] && strata=(--no-strata)
python3 -I "$here/score_dehost.py" --truth "$truth" --calls "$out/$run.calls.tsv" --run "$run" \
    --table "$out/results.tsv" ${strata[@]+"${strata[@]}"}
