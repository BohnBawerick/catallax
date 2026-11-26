"""Pick one shard of the frozen 3,000 census models and find each model's source file.

The public list of model ids and source hashes is in the repository. The archive member
names are not, so this matches every STEP member under a model's folder against the
recorded SHA-256 and writes the [model_id, member] list for the shard. The list and
the extracted files stay in the scratch folder.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--shard', type=int, required=True)
    parser.add_argument('--shards', type=int, required=True)
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--skip', type=Path, help='JSON list of model ids already measured elsewhere')
    parser.add_argument('--only', type=Path, help='JSON list of model ids to take; replaces --skip')
    args = parser.parse_args()
    hashes = json.loads((ROOT / 'docs/research/step-census/source-hashes.json').read_text())
    skip = set(json.loads(args.skip.read_text())) if args.skip and args.skip.is_file() else set()
    assert skip <= set(hashes)
    remaining = [k for k in sorted(hashes) if k not in skip]
    if args.only:
        only = set(json.loads(args.only.read_text()))
        assert only <= set(hashes)
        skip = set()
        remaining = [k for k in sorted(hashes) if k in only]
    keys = [k for i, k in enumerate(remaining) if i % args.shards == args.shard]
    if args.limit:
        keys = keys[:args.limit]
    names = subprocess.run(['bsdtar', '-tf', str(args.archive)], check=True, capture_output=True,
                           text=True).stdout.splitlines()
    wanted = set(keys)
    members = [n for n in names if n.lower().endswith(('.step', '.stp')) and n.split('/')[0] in wanted]
    listing = args.scratch / 'members.txt'
    listing.write_text('\n'.join(members) + '\n')
    sources = args.scratch / 'sources'
    sources.mkdir(parents=True, exist_ok=True)
    subprocess.run(['bsdtar', '-xf', str(args.archive), '-C', str(sources), '-T', str(listing)], check=True)
    selection, missing = [], []
    for key in keys:
        found = None
        for name in members:
            if name.split('/')[0] != key:
                continue
            path = sources / name
            if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == hashes[key]:
                found = name
                break
        if found is None:
            missing.append(key)
        else:
            selection.append([key, found])
    (args.scratch / 'selection.json').write_text(json.dumps(selection))
    print(json.dumps(dict(shard=args.shard, shards=args.shards, skipped=len(skip), models=len(keys), found=len(selection),
                          missing=missing[:20])))
    if missing:
        raise SystemExit('some models have no member with the recorded hash')


if __name__ == '__main__':
    main()
