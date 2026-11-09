#!/usr/bin/env python3
"""Build, run sequential comparisons, and remove owned runtimes and build output."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BENCH = ROOT / '.claude/skills/verify-catallax/scripts/head-to-head'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='fresh external evidence directory')
    parser.add_argument('--scratch', type=Path, default=Path('/tmp'))
    parser.add_argument('--case', action='append', default=[])
    args = parser.parse_args()
    out = args.out.resolve()
    if out == ROOT or ROOT in out.parents:
        parser.error('keep raw evidence outside the checkout')
    args.scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='catallax-comparison-', dir=args.scratch) as owned:
        work = Path(owned)
        env = dict(os.environ, CARGO_BUILD_JOBS='6', CARGO_TARGET_DIR=str(work / 'target'),
                   UV_CACHE_DIR=str(work / 'uv-cache'), UV_PYTHON_INSTALL_DIR=str(work / 'python-install'))
        def run(*command):
            subprocess.run(list(map(str, command)), cwd=ROOT, env=env, check=True)
        run('uv', 'venv', '--python', '3.11', work / 'python')
        python = work / 'python/bin/python'
        run('uv', 'pip', 'sync', '--python', python, '--require-hashes', BENCH / 'requirements.lock')
        run('nice', '-n', '10', 'cargo', 'build', '--release', '-j', '6',
            '-p', 'catallax-cli', '--bin', 'catallax-verify', '--example', 'kernel_compare',
            '-p', 'catallax-compare', '--example', 'measure')
        measure = subprocess.check_output([str(work / 'target/release/examples/measure'), '--export'], cwd=ROOT)
        (work / 'measure-cases.json').write_bytes(measure)
        command = ['taskset', '-c', str(min(os.sched_getaffinity(0))), 'nice', '-n', '10',
                   sys.executable, '-B', str(HERE / 'run.py'), '--out', str(out),
                   '--python', str(python), '--binary', str(work / 'target/release/catallax-verify'),
                   '--observer', str(work / 'target/release/examples/kernel_compare'),
                   '--measure-cases', str(work / 'measure-cases.json')]
        for case in args.case:
            command.extend(['--case', case])
        run(*command)
        run(sys.executable, '-B', HERE / 'report.py', out / 'results.json',
            '--data', out / 'data.json', '--summary', out / 'summary.md')


if __name__ == '__main__':
    main()
