"""Build the census observer from a clean revision and bind its actual artifact."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, required=True)
    args = parser.parse_args()
    assert not args.record.resolve().is_relative_to(ROOT), 'build evidence belongs outside the checkout'
    assert not git('status', '--porcelain'), 'observer build requires a clean checkout'
    revision = git('rev-parse', 'HEAD')
    tree = git('rev-parse', 'HEAD^{tree}')
    command = ['cargo', 'build', '--release', '-p', 'catallax-import', '--example', 'r8_census',
        '--message-format=json']
    log = args.record.with_suffix('.cargo.jsonl')
    with log.open('xb') as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream, check=True, timeout=3600)
    artifacts = []
    for line in log.read_text().splitlines():
        row = json.loads(line)
        if (row['reason'] == 'compiler-artifact' and row['target']['name'] == 'r8_census'
                and row['target']['kind'] == ['example'] and row.get('executable')):
            artifacts.append(row)
    assert len(artifacts) == 1
    artifact = artifacts[0]
    assert artifact['profile']['test'] is False
    assert revision == git('rev-parse', 'HEAD') and not git('status', '--porcelain')
    binary = Path(artifact['executable']).resolve()
    record = {'format': 'catallax-census-observer-build/1', 'source_revision': revision,
        'git_tree': tree, 'build_command': command, 'artifact_path': str(binary),
        'binary_sha256': digest(binary), 'cargo_artifact': artifact,
        'cargo_json_sha256': digest(log), 'cargo_lock_sha256': digest(ROOT / 'Cargo.lock'),
        'build_tool_sha256': digest(Path(__file__)),
        'observer_sha256': digest(ROOT / 'crates/catallax-import/examples/r8_census.rs'),
        'rustc': subprocess.check_output(['rustc', '--version'], text=True).strip()}
    with args.record.open('x') as stream:
        stream.write(json.dumps(record, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'source_revision': revision, 'binary_sha256': record['binary_sha256'],
        'build_record': str(args.record)}), flush=True)


if __name__ == '__main__':
    main()
