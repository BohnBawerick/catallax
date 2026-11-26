#!/usr/bin/env bash
# The same shard run over the models listed in census-skip.json only: the files first
# measured on the Mac, repeated here so that every census row comes from one kind of machine.
# usage: cloud-shard-only.sh SHARD SHARDS [LIMIT]
export P5_ONLY=tools/p5/census-skip.json
exec bash "$(dirname "$0")/cloud-shard.sh" "$@"
