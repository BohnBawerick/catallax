#!/usr/bin/env python3
"""Drive the public Catallax exchange and retain reproducible evidence."""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import statistics
import subprocess
import sys
import tempfile
import time

CHILDREN = set()

def interrupt(signum, frame):
    for child in tuple(CHILDREN):
        if child.poll() is None:
            child.terminate()
    raise KeyboardInterrupt

signal.signal(signal.SIGINT, interrupt)
signal.signal(signal.SIGTERM, interrupt)

CONTRACT = 'catallax-request/1'
TARGETS = {'stock': 100, 'boolean': 500, 'recipe': 1000, 'compare': 2000}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_hash(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        json.dump(value, f, sort_keys=True)
        f.write('\n')
        name = f.name
    os.replace(name, path)


def load(paths, source, adapter=None):
    if adapter:
        spec = importlib.util.spec_from_file_location('problem_adapter', adapter)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return [item for path in paths for item in module.load(Path(path))]
    problems = []
    for name in paths:
        path = Path(name)
        files = sorted(path.rglob('*.json')) if path.is_dir() else [path]
        for file in files:
            values = [json.loads(line) for line in file.read_text().splitlines() if line.strip()] if file.suffix == '.jsonl' else [json.loads(file.read_text())]
            for value in values:
                if isinstance(value, dict) and 'question' in value:
                    value = value['question']
                fmt = value.get('format') if isinstance(value, dict) else None
                if fmt == CONTRACT:
                    problems.append({'source': value['source'], 'problem': value['problem']})
                elif fmt == 'exact-brep-question/2':
                    problems.append({'source': 'question', 'problem': value})
                elif fmt == 'catallax-recipe/3':
                    problems.append({'source': 'compare' if source == 'compare' else 'recipe', 'problem': value})
                elif not path.is_dir():
                    raise ValueError(f'{file}: unsupported format {fmt!r}; use --adapter for another source')
    if not problems:
        raise ValueError('no questions or recipes found in the supplied paths')
    return problems


def classification(item):
    problem = item['problem']
    kind = item['source']
    if kind == 'question':
        kind = 'boolean' if problem.get('steps') else 'stock'
    operations = sorted({s.get('op', '?') for s in problem.get('steps', [])})
    shapes = sorted({o.get('shape', '?') for o in problem.get('operands', {}).values()})
    return kind, '/'.join([kind, '+'.join(operations), '+'.join(shapes)])


class Worker:
    def __init__(self, binary, core=None, expected_hash=None):
        if expected_hash is not None and file_hash(binary) != expected_hash:
            raise RuntimeError("binary changed during this run; copy the build to a stable path and resume")
        command = [str(Path(binary).resolve()), 'stream']
        if core is not None:
            command = ['taskset', '-c', str(core), *command]
        self.child = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        CHILDREN.add(self.child)
        if expected_hash is not None and file_hash(binary) != expected_hash:
            self.close()
            raise RuntimeError("binary changed while starting a worker; use a stable build path")

    def answer(self, request):
        self.child.stdin.write(canonical(request) + '\n')
        self.child.stdin.flush()
        line = self.child.stdout.readline()
        if not line:
            code = self.child.wait()
            raise RuntimeError(f'worker exited {code}: {self.child.stderr.read()}; retain the run and inspect the failing request')
        result = json.loads(line)
        if result.get('format') != 'catallax-response/1':
            raise RuntimeError('worker returned an unsupported response format')
        return result

    def close(self):
        try:
            self.child.stdin.close()
        except BrokenPipeError:
            pass
        try:
            self.child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.child.terminate()
            self.child.wait()
        self.child.stdout.close()
        self.child.stderr.close()
        CHILDREN.discard(self.child)


def valid_artifact(directory, key):
    try:
        meta = json.loads((directory / 'artifact.json').read_text())
        return meta['key'] == key and all(file_hash(directory / name) == sha for name, sha in meta['files'].items())
    except (OSError, ValueError, KeyError):
        return False


def run(args):
    started = time.perf_counter()
    binary = Path(args.binary).resolve()
    identity = file_hash(binary)
    items = load(args.inputs, args.source, args.adapter)
    out = Path(args.out).resolve()
    cache = Path(args.cache).resolve() if args.cache else None
    cores = [int(x) for x in args.cores.split(',')] if args.cores else []
    if cores and not set(cores) <= os.sched_getaffinity(0):
        raise ValueError('requested cores are outside this process affinity; inspect doctor')
    config = {'binary_sha256': identity, 'max_work': args.max_work, 'contract': CONTRACT}
    keys = [digest({**config, **item}) for item in items]
    if args.dry_run:
        return {'dry_run': True, 'cases': len(items), 'out': str(out), 'keys': keys, 'fresh': args.fresh, 'policy': 'no files or workers created'}
    out.mkdir(parents=True, exist_ok=True)
    atomic(out / 'inputs.json', items)
    atomic(out / 'config.json', {**config, 'binary': str(binary), 'workers': args.workers, 'cores': cores, 'command': sys.argv, 'load_average': os.getloadavg()})
    unique = dict(zip(keys, items))
    batches = [[] for _ in range(args.workers)]
    for i, pair in enumerate(unique.items()):
        batches[i % args.workers].append(pair)

    def batch(index):
        results = {}
        worker = None
        try:
            for key, item in batches[index]:
                destination = out / 'cases' / key
                hit = False
                if not args.fresh and valid_artifact(destination, key):
                    hit = True
                elif not args.fresh and cache and valid_artifact(cache / key, key):
                    if destination.exists():
                        shutil.rmtree(destination)
                    shutil.copytree(cache / key, destination)
                    hit = True
                else:
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with tempfile.TemporaryDirectory(prefix='.pending-', dir=destination.parent) as scratch:
                        scratch = Path(scratch)
                        request = {**item, 'format': CONTRACT, 'max_work': args.max_work, 'out': str(scratch)}
                        atomic(scratch / 'request.json', request)
                        if worker is None:
                            worker = Worker(binary, cores[index % len(cores)] if cores else None, identity)
                        response = worker.answer(request)
                        atomic(scratch / 'response.json', response)
                        files = {p.name: file_hash(p) for p in scratch.iterdir() if p.is_file()}
                        atomic(scratch / 'artifact.json', {'key': key, 'files': files})
                        if destination.exists():
                            shutil.rmtree(destination)
                        shutil.copytree(scratch, destination)
                    if cache:
                        cache.mkdir(parents=True, exist_ok=True)
                        # Complete entries publish by rename. A concurrent winner is equivalent.
                        with tempfile.TemporaryDirectory(prefix='.pending-', dir=cache) as temporary:
                            entry = Path(temporary) / 'entry'
                            shutil.copytree(destination, entry)
                            try:
                                entry.rename(cache / key)
                            except OSError:
                                if not valid_artifact(cache / key, key):
                                    raise RuntimeError(f'corrupt cache entry: {cache / key}; use a new cache directory')
                response = json.loads((destination / 'response.json').read_text())
                kind, group = classification(item)
                record = {'key': key, 'id': item['problem'].get('id', item['problem'].get('name')), 'path': str(destination), 'cache_hit': hit, 'class': kind, 'group': group, 'target_ms': TARGETS[kind], 'response': response}
                results[key] = record
                # One checkpoint per worker avoids shared writes and survives interrupted runs.
                with (out / f'checkpoint-{index}.jsonl').open('a') as log:
                    log.write(canonical({'key': key, 'path': str(destination)}) + '\n')
                    log.flush()
        finally:
            if worker:
                worker.close()
        return results

    merged = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(batch, range(args.workers)):
            merged.update(result)
    report = {'format': 'catallax-run/1', 'config': config, 'elapsed_ms': (time.perf_counter() - started) * 1000, 'cases': [merged[key] for key in keys]}
    atomic(out / 'run.json', report)
    return summarize(report)


def summarize(report):
    groups, slow, issues = {}, [], {}
    crashes = 0
    disagreements = 0
    for case in report['cases']:
        response = case['response']
        manifest = response['manifest']
        cause = (manifest.get('refusal') or {}).get('code', manifest.get('status', 'unknown'))
        crashes += manifest.get('status') == 'internal-error'
        comparison = manifest.get('comparison')
        if comparison:
            cause = 'compare.' + comparison['verdict']
            disagreements += comparison['verdict'] == 'disagree'
            for finding in comparison['findings']:
                issues.setdefault(finding['code'], []).append(case['key'])
        elif manifest.get('refusal'):
            issues.setdefault(cause, []).append(case['key'])
        groups[cause] = groups.get(cause, 0) + 1
        if response['wall_ms'] > case['target_ms']:
            steps = manifest.get('stats', {}).get('per_step')
            boolean_ms = sum(step.get('wall_ms', 0) for step in steps) if steps is not None else None
            phase = ('boolean' if boolean_ms >= response['wall_ms'] / 2 else 'outside-boolean') if boolean_ms is not None else 'unmeasured-phase'
            outside_ms = max(0, response['wall_ms'] - boolean_ms) if boolean_ms is not None else None
            slow.append({k: case[k] for k in ['id', 'key', 'group', 'target_ms']} | {'wall_ms': response['wall_ms'], 'work': response['work'], 'cause': cause, 'phase': phase, 'boolean_ms': boolean_ms, 'outside_boolean_ms': outside_ms, 'cached_timing': case['cache_hit']})
    performance = {}
    for case in slow:
        label = case['group'] + '/' + case['phase'] + '/' + case['cause']
        performance.setdefault(label, []).append(case['key'])
    return {'cases': len(report['cases']), 'elapsed_ms': report['elapsed_ms'], 'cache_hits': sum(c['cache_hit'] for c in report['cases']), 'outcomes': groups, 'issue_groups': issues, 'crashes': crashes, 'comparison_disagreements': disagreements, 'performance_groups': performance, 'slow': sorted(slow, key=lambda c: -c['wall_ms'])}


def comparable(manifest):
    return {k: v for k, v in manifest.items() if k not in ('stats', 'kernel')}


def compare(left, right):
    a, b = [json.loads(Path(path).read_text()) for path in (left, right)]
    failures = []
    if len(a['cases']) != len(b['cases']):
        failures.append({'cause': 'case-count', 'left': len(a['cases']), 'right': len(b['cases'])})
    for i, (x, y) in enumerate(zip(a['cases'], b['cases'])):
        # Inputs, including order and duplicates, are part of the proof.
        px, py = [json.loads((Path(c['path']) / 'request.json').read_text()) for c in (x, y)]
        for request in (px, py):
            request.pop('out', None)
        manifests = [json.loads((Path(c['path']) / 'manifest.json').read_text()) for c in (x, y)]
        cause = None
        if any(m.get('status') == 'internal-error' for m in manifests):
            cause = 'internal-error'
        elif px != py:
            cause = 'input-mismatch'
        elif comparable(manifests[0]) != comparable(manifests[1]):
            cause = 'manifest'
        elif a['config']['binary_sha256'] == b['config']['binary_sha256'] and x['response']['work'] != y['response']['work']:
            cause = 'work-mismatch'
        else:
            sx, sy = [(Path(c['path']) / 'answer.step') for c in (x, y)]
            if sx.exists() != sy.exists() or (sx.exists() and sx.read_bytes() != sy.read_bytes()):
                cause = 'step-bytes'
        if cause:
            failures.append({'index': i, 'id': x['id'], 'cause': cause, 'left': x['path'], 'right': y['path']})
    return {'equal': not failures, 'cases': len(a['cases']), 'failures': failures}


def pair(args):
    items = load(args.inputs, args.source, args.adapter)
    out = Path(args.out).resolve()
    if args.dry_run:
        return {'dry_run': True, 'cases': len(items), 'out': str(out), 'repeats': args.repeats, 'policy': 'no files or workers created'}
    out.mkdir(parents=True, exist_ok=True)
    samples = []
    identities = {str(i): {'path': str(Path(binary).resolve()), 'sha256': file_hash(binary)} for i, binary in enumerate((args.baseline, args.candidate))}
    atomic(out / 'config.json', {'binaries': identities, 'core': int(args.cores), 'max_work': args.max_work, 'repeats': args.repeats, 'command': sys.argv, 'load_average': os.getloadavg()})
    # AB/BA alternation reduces drift. Each process stays alive across all pairs.
    workers = []
    try:
        for side, binary in enumerate((args.baseline, args.candidate)):
            workers.append(Worker(binary, int(args.cores), identities[str(side)]['sha256']))
    except BaseException:
        for worker in workers:
            worker.close()
        raise
    try:
        for index, item in enumerate(items):
            for repeat in range(args.repeats + 1):
                row = {'index': index, 'repeat': repeat, 'warmup': repeat == 0}
                responses = {}
                for side in ((0, 1) if repeat % 2 == 0 else (1, 0)):
                    dest = out / str(index) / str(repeat) / str(side)
                    dest.mkdir(parents=True, exist_ok=True)
                    response = workers[side].answer({**item, 'format': CONTRACT, 'max_work': args.max_work, 'out': str(dest)})
                    atomic(dest / 'response.json', response)
                    row[str(side)] = response['wall_ms']
                    responses[side] = response
                step_files = [out / str(index) / str(repeat) / str(side) / 'answer.step' for side in (0, 1)]
                row['same_answer'] = comparable(responses[0]['manifest']) == comparable(responses[1]['manifest']) and all(r['manifest'].get('status') != 'internal-error' for r in responses.values())
                row['same_step'] = [p.read_bytes() if p.exists() else None for p in step_files][0] == [p.read_bytes() if p.exists() else None for p in step_files][1]
                samples.append(row)
                atomic(out / 'pairs.json', samples)
    finally:
        for worker in workers:
            worker.close()
    measured = [s for s in samples if not s['warmup']]
    cases = []
    for index in range(len(items)):
        subset = [s for s in measured if s['index'] == index]
        cases.append({'index': index, 'baseline_median_ms': statistics.median(s['0'] for s in subset), 'candidate_median_ms': statistics.median(s['1'] for s in subset), 'median_delta_ms': statistics.median(s['1'] - s['0'] for s in subset)})
    return {'pairs': len(measured), 'equal': all(s['same_answer'] and s['same_step'] for s in samples), 'median_delta_ms': statistics.median(s['1'] - s['0'] for s in measured), 'cases': cases, 'samples': str(out / 'pairs.json')}


def legacy(args):
    items = load(args.inputs, 'question')
    if args.dry_run:
        return {'dry_run': True, 'cases': len(items), 'out': str(Path(args.out).resolve()), 'policy': 'no files or workers created'}
    start = time.perf_counter()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, item in enumerate(items):
        if item['source'] != 'question':
            raise ValueError('legacy replay accepts only questions')
        dest = out / str(index)
        dest.mkdir(exist_ok=True)
        atomic(dest / 'question.json', item['problem'])
        command = [str(Path(args.binary).resolve()), 'answer', str(dest / 'question.json'), '--out', str(dest)]
        if args.cores:
            command = ['taskset', '-c', args.cores, *command]
        one = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True)
        atomic(dest / 'process.json', {'command': command, 'exit': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr})
        rows.append({'index': index, 'wall_ms': (time.perf_counter() - one) * 1000, 'exit': result.returncode})
        atomic(out / 'progress.json', rows)
    result = {'cases': len(rows), 'elapsed_ms': (time.perf_counter() - start) * 1000, 'rows': rows}
    atomic(out / 'legacy.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog='Examples: doctor --binary target/release/catallax-verify; run INPUT.jsonl --binary target/release/catallax-verify --out /tmp/proof --cache /tmp/cache --cores 2 --json; compare /tmp/main/run.json /tmp/candidate/run.json --json. All private inputs are runtime paths. No network access.')
    commands = parser.add_subparsers(dest='command', required=True)
    doctor = commands.add_parser('doctor', help='read-only build and CPU-affinity checks')
    doctor.add_argument('--binary', required=True)
    for name in ('run', 'legacy', 'pair'):
        p = commands.add_parser(name, help={'run': 'persistent workers, content cache and checkpoint resume', 'legacy': 'measure the old one-process-per-question replay', 'pair': 'uncached matched AB/BA timing on one pinned core'}[name])
        p.add_argument('inputs', nargs='+', help='JSON, JSONL corpus, or recursively discovered directory')
        p.add_argument('--out', required=True, help='owned evidence directory; retained after exit')
        p.add_argument('--dry-run', action='store_true', help='preview the operation without creating workers or writing files')
        p.add_argument('--cores', required=name == 'pair', help='CPU list for run; exactly one CPU for pair/legacy')
        if name == 'pair':
            p.add_argument('--baseline', required=True)
            p.add_argument('--candidate', required=True)
            p.add_argument('--repeats', type=int, default=3)
        else:
            p.add_argument('--binary', required=True)
        if name != 'legacy':
            p.add_argument('--max-work', type=int, default=100_000_000)
            p.add_argument('--source', choices=['auto', 'recipe', 'compare'], default='auto')
            p.add_argument('--adapter', help='trusted Python file exporting load(Path) -> source/problem records')
        if name == 'run':
            p.add_argument('--workers', type=int, choices=range(1, 65), default=1)
            p.add_argument('--cache', help='shared content-addressed cache directory')
            p.add_argument('--fresh', action='store_true', help='recompute and replace case artifacts, bypassing cache and checkpoints')
    p = commands.add_parser('compare', help='strict STEP and manifests excluding stats/build identity')
    p.add_argument('baseline')
    p.add_argument('candidate')
    p = commands.add_parser('report', help='correctness groups and separate interactive performance failures')
    p.add_argument('run')
    p = commands.add_parser('cleanup', help='confirm no persistent processes are retained; preserve all evidence')
    p.add_argument('--out', required=True)
    p.add_argument('--dry-run', action='store_true', help='show the cleanup policy without changing anything')
    for p in commands.choices.values():
        p.add_argument('--json', action='store_true', help='emit stable machine-readable JSON')
    args = parser.parse_args()
    if hasattr(args, 'max_work') and not 0 <= args.max_work <= 2**64 - 1:
        parser.error('--max-work must be an unsigned 64-bit integer')
    if hasattr(args, 'repeats') and args.repeats < 1:
        parser.error('--repeats must be positive')
    try:
        if args.command == 'doctor':
            result = {'binary': str(Path(args.binary).resolve()), 'sha256': file_hash(args.binary), 'about': json.loads(subprocess.check_output([args.binary, 'about'], text=True)), 'allowed_cores': sorted(os.sched_getaffinity(0)), 'python': sys.version.split()[0]}
            if result['about'].get('request') != CONTRACT:
                raise ValueError('binary does not implement catallax-request/1; build catallax-verify and pass its path')
        elif args.command == 'run':
            result = run(args)
        elif args.command == 'compare':
            result = compare(args.baseline, args.candidate)
        elif args.command == 'report':
            result = summarize(json.loads(Path(args.run).read_text()))
        elif args.command == 'pair':
            result = pair(args)
        elif args.command == 'legacy':
            result = legacy(args)
        else:
            result = {'policy': 'workers close on command exit; evidence is never removed', 'evidence_exists': Path(args.out).exists(), 'dry_run': args.dry_run}
        print(json.dumps(result, indent=None if args.json else 2, sort_keys=True))
        return 1 if result.get('equal') is False or result.get('crashes', 0) or result.get('comparison_disagreements', 0) else 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(json.dumps({'error': str(error), 'action': 'retain evidence, run doctor, then resume using the same --out'}), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
