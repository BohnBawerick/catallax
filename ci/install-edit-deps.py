#!/usr/bin/env python3
"""Install the edit probes' packages from cache or signed Ubuntu mirrors."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

MIRRORS = ('https://archive.ubuntu.com/ubuntu',
           'https://mirrors.edge.kernel.org/ubuntu',
           'https://mirrors.mit.edu/ubuntu')
PACKAGES = ['mesa-vulkan-drivers', 'valgrind', 'linux-tools-generic']


def command(args):
    return subprocess.run(['timeout', '--kill-after=10s', '180s', *args], check=False).returncode == 0


def install(cache, codename, run=command, sleep=time.sleep, has_perf=None):
    options = ['-o', f'Dir::Cache::archives={cache / "archives"}',
               '-o', f'Dir::State::lists={cache / "lists"}',
               '-o', 'Acquire::Retries=1', '-o', 'Acquire::http::Timeout=15',
               '-o', 'Acquire::https::Timeout=15', '-o', 'APT::Update::Error-Mode=any',
               '-o', 'APT::Keep-Downloaded-Packages=true']
    source = cache / 'sources.list'

    def source_options():
        return (['-o', f'Dir::Etc::sourcelist={source}', '-o', 'Dir::Etc::sourceparts=-']
                if source.exists() else [])

    if has_perf is None:
        has_perf = lambda: subprocess.run(
            ['apt-cache', *options, *source_options(), 'show', 'linux-perf'],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False).returncode == 0

    def packages():
        # Newer Ubuntu images ship perf separately from linux-tools-*.
        return PACKAGES + (['linux-perf'] if has_perf() else [])

    apt = ['sudo', 'apt-get', *options]
    # Matching cached indices and archives can install without working DNS.
    if run([*apt, *source_options(), '--no-download', 'install', '-y', '--no-install-recommends', *packages()]):
        return {'status': 'ready', 'source': 'cache'}
    for attempt, mirror in enumerate(MIRRORS):
        if attempt:
            sleep(10 * 2 ** (attempt - 1))
        source.write_text(''.join(
            f'deb [signed-by=/usr/share/keyrings/ubuntu-archive-keyring.gpg] {mirror} '
            f'{codename}{pocket} main restricted universe multiverse\n'
            for pocket in ('', '-updates', '-security', '-backports')))
        selected = [*apt, *source_options()]
        print(f'Package infrastructure attempt {attempt + 1}: {mirror}', flush=True)
        if (run([*selected, 'update']) and run([
                *selected, 'install', '-y', '--no-install-recommends', *packages()])):
            return {'status': 'ready', 'source': mirror}
    return {'status': 'infrastructure-unavailable',
            'reason': 'Ubuntu package fetch failed after cache and three signed mirrors',
            'packages': PACKAGES}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    lists = args.cache / 'lists'
    if not lists.exists():
        lists.mkdir()
        for path in Path('/var/lib/apt/lists').iterdir():
            if path.is_file() and path.name != 'lock':
                shutil.copyfile(path, lists / path.name)
    subprocess.run(['sudo', 'mkdir', '-p', str(lists / 'partial'),
                    str(args.cache / 'archives/partial')], check=True)
    release = dict(line.split('=', 1) for line in Path('/etc/os-release').read_text().splitlines()
                   if '=' in line)
    codename = release['VERSION_CODENAME'].strip('"')
    if release['ID'].strip('"') != 'ubuntu' or not codename.isalnum():
        raise ValueError('edit probe packages require an Ubuntu runner')
    result = install(args.cache.resolve(), codename)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + '\n')
    # apt runs as root; the cache action needs to read its result as the runner.
    subprocess.run(['sudo', 'chown', '-R', f'{os.getuid()}:{os.getgid()}', str(args.cache)], check=True)
    if result['status'] != 'ready':
        message = result['reason'] + '; no measurements ran. Rerun the entire workflow.'
        print(f'::error title=Package infrastructure unavailable::{message}', flush=True)
        if os.environ.get('GITHUB_STEP_SUMMARY'):
            with Path(os.environ['GITHUB_STEP_SUMMARY']).open('a') as output:
                output.write(message + '\n')
        raise SystemExit(75)


if __name__ == '__main__':
    main()
