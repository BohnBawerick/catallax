#!/usr/bin/env bash
# usage: run-census.sh catallax|occt WORKERS [LIMIT]
# P5_PY and P5_CORPUS override the work PC's Python and corpus paths.
set -euo pipefail
K=$HOME/offload/keep/cx-p5-speed-owner
PY=${P5_PY:-$K/venv/bin/python}
C=${P5_CORPUS:-$K/corpus}
case $1 in
  catallax)
    exec "$PY" -B "$K/src/tools/research/p5_census.py" --binary "$K/bin/opened_import" \
      --sources "$C/selected-step" --selection "$C/selection-private.json" \
      --out "$K/census/catallax" --workers "$2" ${3:+--limit "$3"} ;;
  occt)
    exec "$PY" -B "$K/src/tools/research/occt_census.py" collect --sources "$C/selected-step" \
      --selection "$C/selection-private.json" \
      --hashes "$K/src/docs/research/step-census/source-hashes.json" \
      --out "$K/census/occt" --workers "$2" --timeout 600 ;;
esac
