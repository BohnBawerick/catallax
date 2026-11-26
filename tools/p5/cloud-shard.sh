#!/usr/bin/env bash
# One shard of the step 5.1 census on a fresh Linux machine, both kernels, same caps.
# usage: cloud-shard.sh SHARD SHARDS [LIMIT]
# Prints P5OCCT and P5ROW lines, one per file, and P5SUMMARY lines.
set -uo pipefail
SHARD=$1 SHARDS=$2 LIMIT=${3:-0}
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
W=${P5_SCRATCH:-${TMPDIR:-/tmp}/p5-census}
mkdir -p "$W"
cd "$ROOT"
echo "P5HOST cores=$(nproc) cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs) mem_mb=$(awk '/MemTotal/{print int($2/1024)}' /proc/meminfo) kernel=$(uname -sr) revision=$(git rev-parse HEAD 2>/dev/null) shard=$SHARD/$SHARDS limit=$LIMIT"
step() { echo "P5STEP $1 at $(date -u +%H:%M:%S)"; }

step tools
if ! command -v bsdtar >/dev/null; then
  (sudo -n apt-get update -qq && sudo -n apt-get install -y -qq libarchive-tools) >/dev/null 2>&1 \
    || (apt-get update -qq && apt-get install -y -qq libarchive-tools) >/dev/null 2>&1
fi
command -v bsdtar >/dev/null || { echo "P5FAIL no bsdtar and it could not be installed"; exit 3; }
if [ ! -x "$W/uv/uv" ]; then
  mkdir -p "$W/uv"
  curl -LsSf -o "$W/uv/uv.tar.gz" https://github.com/astral-sh/uv/releases/latest/download/uv-x86_64-unknown-linux-gnu.tar.gz \
    && tar xzf "$W/uv/uv.tar.gz" -C "$W/uv" --strip-components 1 || { echo "P5FAIL uv download"; exit 3; }
fi
export UV_CACHE_DIR=$W/uv/cache UV_PYTHON_INSTALL_DIR=$W/uv/python
[ -x "$W/venv/bin/python" ] || "$W/uv/uv" venv -q --python 3.11 "$W/venv" || { echo "P5FAIL python"; exit 3; }
"$W/uv/uv" pip install -q --python "$W/venv/bin/python" cadquery-ocp==8.0.1.1.0 || { echo "P5FAIL cadquery-ocp"; exit 3; }
"$W/venv/bin/python" -c "import OCP; print('P5OCP', OCP.__version__)" || { echo "P5FAIL OCP import"; exit 3; }

step archive
# The corpus sits in a private repository the cloud can read, as 17 parts of the public ABC chunk.
if [ "$(stat -c %s "$W/abc_0000_step_v00.7z" 2>/dev/null)" != 1594129754 ]; then
  rm -rf "$W/corpus-repo"
  git clone -q --depth 1 https://github.com/ZIGRAWBOHNBAWERICK/cad-test-corpus "$W/corpus-repo" \
    || { echo "P5FAIL corpus clone"; exit 4; }
  cat "$W"/corpus-repo/abc/abc_0000_step_v00.7z.part-* > "$W/abc_0000_step_v00.7z"
  rm -rf "$W/corpus-repo"
fi
[ "$(stat -c %s "$W/abc_0000_step_v00.7z")" = 1594129754 ] || { echo "P5FAIL archive size"; exit 4; }
[ "$(md5sum "$W/abc_0000_step_v00.7z" | cut -d' ' -f1)" = 695388be7a278c7798c8c8ae239772ac ] || { echo "P5FAIL archive md5"; exit 4; }

step select
"$W/venv/bin/python" -B tools/p5/cloud_select.py --archive "$W/abc_0000_step_v00.7z" --scratch "$W" \
  --shard "$SHARD" --shards "$SHARDS" --limit "$LIMIT" --skip tools/p5/census-skip.json || { echo "P5FAIL select"; exit 5; }
rm -f "$W/abc_0000_step_v00.7z"

step build
cargo build --release -p catallax-import --example opened_import 2>&1 | tail -n 3
BIN=${CARGO_TARGET_DIR:-$ROOT/target}/release/examples/opened_import
[ -x "$BIN" ] || { echo "P5FAIL build"; exit 6; }

export P5_CPU_CAP=1 P5_MIN_AVAILABLE_MB=1500
WORKERS=${P5_WORKERS:-$(nproc)}
step occt
"$W/venv/bin/python" -B tools/research/p5_census.py --kernel occt --sources "$W/sources" \
  --selection "$W/selection.json" --out "$W/occt" --workers "$WORKERS"
step catallax
"$W/venv/bin/python" -B tools/research/p5_census.py --kernel catallax --binary "$BIN" --sources "$W/sources" \
  --selection "$W/selection.json" --out "$W/catallax" --workers "$WORKERS" ${P5_TIME_BUDGET:+--time-budget "$P5_TIME_BUDGET"}
step done
