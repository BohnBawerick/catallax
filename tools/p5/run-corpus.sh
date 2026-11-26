#!/usr/bin/env bash
# usage: run-corpus.sh BIN ROOT OUT LABEL THREADS [perf|callgrind] [ROWS]
# Runs each predictor-corpus row in a process of its own: no request work limit, a 180 s
# cut-off per build and per export, and a 240 s outer guard. With "perf", one-thread runs
# record user instructions and cycles.
set -uo pipefail
BIN=$1 ROOT=$2 OUT=$3 LABEL=$4 THREADS=$5 MODE=${6:-} ROWS=${7:-}
mkdir -p "$OUT"
if [ -n "$ROWS" ]; then ids=${ROWS//,/ }; else
  ids=$(python3 -c "import json,sys;print(' '.join(r['id'] for r in json.load(open(sys.argv[1]))['rows']))" "$ROOT/tools/predictor/corpus.json")
fi
extra=()
if [ -n "${P5_SOURCES:-}" ]; then extra=(--sources "$P5_SOURCES" --selection "$P5_SELECTION"); fi
for id in $ids; do
  for t in ${THREADS//,/ }; do
    f=$OUT/$id.t$t
    [ -e "$f.done" ] && continue
    rm -f "$f.jsonl"
    load=$(cut -d' ' -f1-3 /proc/loadavg)
    start=$(date +%s.%N)
    cmd=("$BIN" run --corpus "$ROOT/tools/predictor/corpus.json" --root "$ROOT" --rows "$id" --threads "$t" --out "$f.jsonl" --label "$LABEL" --max-work 4611686018427387904 --cutoff-seconds 180 "${extra[@]}")
    if [ "$MODE" = perf ] && [ "$t" = 1 ]; then
      timeout 240 perf stat -x, -e instructions:u,cycles:u -o "$f.perf" -- "${cmd[@]}" > "$f.err" 2>&1
    else
      timeout 240 "${cmd[@]}" > "$f.err" 2>&1
    fi
    code=$?
    end=$(date +%s.%N)
    echo "$id threads=$t exit=$code load=$load wall=$(echo "$end - $start" | bc -l 2>/dev/null || python3 -c "print($end-$start)")" >> "$OUT/status.txt"
    : > "$f.done"
  done
done
date -u > "$OUT/complete"
