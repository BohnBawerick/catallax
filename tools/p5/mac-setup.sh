#!/usr/bin/env bash
# One-time Mac setup for the step 5.1 census: uv, Python 3.11 with OCP 8.0.1.1, and the corpus.
set -euo pipefail
K=$HOME/offload/keep/cx-p5-speed-owner
mkdir -p "$K/uv" "$K/census" "$K/logs" "$K/bin" "$K/corpus"
cd "$K/uv"
if [ ! -x uv ]; then
  curl -LsSf -o uv.tar.gz https://github.com/astral-sh/uv/releases/latest/download/uv-aarch64-apple-darwin.tar.gz
  tar xzf uv.tar.gz --strip-components 1
  rm uv.tar.gz
fi
export UV_CACHE_DIR=$K/uv/cache UV_PYTHON_INSTALL_DIR=$K/uv/python
[ -x "$K/venv/bin/python" ] || ./uv venv --python 3.11 "$K/venv"
./uv pip install --python "$K/venv/bin/python" cadquery-ocp==8.0.1.1.0
"$K/venv/bin/python" -c "import OCP, sys; print('OCP', OCP.__version__, sys.version)"
rsync -a bm-dsk-037-wsl:cx-census/selection-private.json "$K/corpus/"
rsync -a bm-dsk-037-wsl:cx-census/selected-step "$K/corpus/"
ls "$K/corpus/selected-step" | wc -l
du -sh "$K/corpus"
echo "mac setup done"
