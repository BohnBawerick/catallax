#!/usr/bin/env python3
"""Install checksum-pinned releases outside the checkout. No kernel code is copied in."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import urllib.request
import zipfile

DATA = Path(__file__).resolve().parent


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def download(pin, path):
    if not path.exists():
        urllib.request.urlretrieve(pin['url'], path)
    if sha(path) != pin['sha256']:
        raise ValueError(f'checksum mismatch: {path.name}')
    return path


def tree_hash(root):
    digest = hashlib.sha256()
    for path in sorted(Path(root).rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc':
            digest.update(path.relative_to(root).as_posix().encode() + b'\0')
            digest.update(bytes.fromhex(sha(path)))
    return digest.hexdigest()


def prepare(out):
    out = Path(out).resolve()
    if out == DATA or DATA.is_relative_to(out) or out.is_relative_to(DATA.parents[4]):
        raise ValueError('tool installation must be outside the repository')
    out.mkdir(parents=True, exist_ok=True)
    pin = json.loads((DATA / 'tools.json').read_text())
    def run(*args):
        subprocess.run(list(map(str, args)), check=True)
    # Check the direct release artifacts as well as the complete dependency lock.
    for name in ('cadquery', 'smlibpy'):
        item = pin['tools'][name]
        download(item, out / item['url'].rsplit('/', 1)[1])
    python = out / 'python/bin/python'
    if not python.exists():
        run('uv', 'venv', '--python', pin['python'], out / 'python')
    run('uv', 'pip', 'sync', '--python', python, '--require-hashes', DATA / 'requirements.lock')
    checker = out / 'checker/bin/python'
    if not checker.exists():
        run('uv', 'venv', '--python', pin['python'], out / 'checker')
    wheel = download(pin['checker'], out / pin['checker']['url'].rsplit('/', 1)[1])
    # Checker is not timed; pin its OCP release. Its dependencies are recorded below.
    run('uv', 'pip', 'install', '--python', checker, wheel)
    archive = download(pin['tools']['pali_ijo'], out / 'pali-source.tar.gz')
    package = out / ('PaliCAD-' + pin['tools']['pali_ijo']['revision']) / 'packages/pali_ijo'
    if not package.exists():
        with tarfile.open(archive) as bundle:
            prefix = f'PaliCAD-{pin["tools"]["pali_ijo"]["revision"]}/packages/pali_ijo/'
            members = [m for m in bundle.getmembers() if m.name.startswith(prefix)]
            bundle.extractall(out, members=members, filter='data')
    sdk = download(pin['dart'], out / 'dart-sdk.zip')
    dart = out / 'dart-sdk/bin/dart'
    if not dart.exists():
        with zipfile.ZipFile(sdk) as bundle:
            bundle.extractall(out)
        for executable in (out / 'dart-sdk/bin').rglob('*'):
            if executable.is_file() and executable.suffix not in ('.snapshot', '.sym', '.dill'):
                executable.chmod(0o755)
    subprocess.run([str(dart), 'pub', 'get', '--enforce-lockfile'], cwd=package, check=True)
    native = out / 'pali-parts'
    run(dart, 'compile', 'exe', '--packages=' + str(package / '.dart_tool/package_config.json'), DATA / 'pali_parts.dart', '-o', native)
    receipt = {'format': 'catallax-comparison-install/1', 'pins_sha256': sha(DATA / 'tools.json'),
               'requirements_sha256': sha(DATA / 'requirements.lock'), 'runner_sha256': sha(DATA / 'pali_parts.dart'),
               'python': str(python), 'python_tree_sha256': tree_hash(out / 'python/lib'),
               'checker': str(checker), 'checker_tree_sha256': tree_hash(out / 'checker/lib'),
               'pali_binary': str(native), 'pali_binary_sha256': sha(native),
               'dart_version': pin['dart']['version'], 'tools': pin['tools']}
    path = out / 'runtime.json'
    path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    prepare(parser.parse_args().out)
