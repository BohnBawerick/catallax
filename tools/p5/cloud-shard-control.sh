#!/usr/bin/env bash
# Control for the step 5.1 census: today's kernel under the old protocol, 20M work per file
# and a 60 s cap, over all 3,000 models, Catallax only. It separates what kernel changes since
# the committed census gave from what removing the work limit gave.
# usage: cloud-shard-control.sh SHARD SHARDS [LIMIT]
export P5_SKIP=/dev/null P5_NO_OCCT=1 P5_WORK=20000000 P5_TIMEOUT=60
exec bash "$(dirname "$0")/cloud-shard.sh" "$@"
