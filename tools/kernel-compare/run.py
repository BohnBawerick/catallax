#!/usr/bin/env python3
"""Run public comparison inputs sequentially, retaining explicit incomplete checks."""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import sys
import time

from policy import DISTANCE_WIDTH_MM, MAX_WORK, SOLVE_LIMIT_SECONDS, TOLERANCES_MM, sweep

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BENCH = ROOT / '.claude/skills/verify-catallax/scripts/head-to-head'
sys.path.insert(0, str(BENCH.parent))
from head_to_head import edited_recipe


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def corpus(measure_cases=None):
    rows = []
    for line in (ROOT / 'parts/fillets/cases.jsonl').read_text().splitlines():
        case = json.loads(line)
        recipe = case['recipe']
        if isinstance(recipe, str):
            recipe = json.loads((ROOT / recipe).read_text())
        rows.append({'id': 'fillets/' + case['id'], 'corpus': 'fillets', 'recipe': recipe,
                     'fillet_case': case})
    for path in sorted((BENCH / 'cases').rglob('*.json')):
        case = json.loads(path.read_text())
        source = ROOT / 'parts' / case['recipe']
        if sha(source) != case['recipe_sha256']:
            raise ValueError(f'{source}: authored benchmark recipe hash changed')
        recipe = json.loads(source.read_text())
        for label in ('construction', 'edit'):
            rows.append({'id': 'five-tool/' + case['recipe'].removesuffix('.json') + '/' + label,
                         'corpus': 'five-tool', 'recipe': recipe if label == 'construction'
                         else edited_recipe(recipe, case['edit']), 'benchmark_case': case,
                         'variant': label})
    seen = set()
    for test in sorted((ROOT / 'crates/catallax-compare').rglob('*.rs')):
        for name in re.findall(r'include_str!\(\s*"([^"]+\.json)"\s*\)', test.read_text()):
            source = (test.parent / name).resolve()
            if source in seen:
                continue
            seen.add(source)
            rows.append({'id': 'comparison/' + source.relative_to(ROOT).as_posix().removesuffix('.json'),
                         'corpus': 'comparison-fixtures', 'recipe': json.loads(source.read_text()),
                         'source': source.relative_to(ROOT).as_posix()})
            if source.is_relative_to(ROOT / 'parts'):
                authored = BENCH / 'cases' / source.relative_to(ROOT / 'parts')
                if authored.exists():
                    rows[-1].update(benchmark_case=json.loads(authored.read_text()), variant='construction')
    if measure_cases is not None:
        for recipe in json.loads(measure_cases.read_text()):
            rows.append({'id': 'comparison-measure/' + recipe['name'],
                         'corpus': 'comparison-measure', 'recipe': recipe})
    return rows


def execute(command, directory, timeout, stdin=None):
    directory.mkdir(parents=True, exist_ok=True)
    save(directory / 'command.json', command)
    start = time.perf_counter()
    with (directory / 'stdout').open('w') as stdout, (directory / 'stderr').open('w') as stderr:
        child = subprocess.Popen(command, stdin=subprocess.PIPE if stdin else subprocess.DEVNULL,
                                 stdout=stdout, stderr=stderr, start_new_session=True,
                                 env=dict(os.environ, RAYON_NUM_THREADS='1', OMP_NUM_THREADS='1',
                                          OPENBLAS_NUM_THREADS='1'))
        try:
            child.communicate(stdin.encode() if stdin else None, timeout=timeout)
            status = 'completed' if child.returncode == 0 else 'crash' if child.returncode < 0 else 'process-error'
        except subprocess.TimeoutExpired:
            status = 'timeout'
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.communicate(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.communicate()
    result = {'status': status, 'exit': child.returncode,
              'wall_ms': (time.perf_counter()-start)*1000,
              'stderr': (directory / 'stderr').read_text(errors='replace')}
    save(directory / 'process.json', result)
    return result


def payload(path, last_line=False):
    try:
        text = path.read_text()
        return json.loads(text.splitlines()[-1] if last_line else text)
    except (OSError, ValueError, IndexError) as error:
        return {'status': 'protocol-error', 'detail': str(error)}


def native(args, directory, request):
    save(directory / 'input.json', request)
    result = execute([str(args.observer), str(directory / 'input.json'), str(directory / 'answer.json')],
                     directory, args.solve_limit)
    if result['status'] != 'completed':
        return result
    return payload(directory / 'answer.json')


def solve(args, case, engine, directory):
    if engine == 'catallax':
        request = {'format': 'catallax-request/1', 'source': 'recipe', 'problem': case['recipe'],
                   'out': str(directory / 'artifacts'), 'max_work': MAX_WORK}
        result = execute([str(args.binary), 'stream'], directory, args.solve_limit,
                         json.dumps(request) + '\n')
        if result['status'] == 'completed':
            response = payload(directory / 'stdout', last_line=True)
            save(directory / 'response.json', response)
            manifest = response.get('manifest')
            if isinstance(manifest, dict):
                result['status'] = manifest['status']
                result['refusal'] = manifest.get('refusal')
                result['kernel_wall_ms'] = response.get('wall_ms')
            else:
                result.update(status='protocol-error', detail=response)
        return result, directory / 'artifacts/answer.step'
    save(directory / 'case.json', case)
    result = execute([str(args.python), '-B', str(HERE / 'worker.py'), engine,
                      str(directory / 'case.json'), str(directory / 'answer.json'),
                      str(directory / 'answer.step')], directory, args.solve_limit)
    if result['status'] == 'completed':
        answer = payload(directory / 'answer.json')
        result.update(answer)
    return result, directory / 'answer.step'


def check(args, case, path, directory, samples):
    if not path.exists():
        return {'status': 'unresolved', 'reason': 'producer did not export STEP'}
    save(directory / 'case.json', case)
    save(directory / 'oracle.json', samples)
    process = execute([str(args.python), '-B', str(HERE / 'check.py'), str(path),
                       str(directory / 'case.json'), str(directory / 'oracle.json'),
                       str(directory / 'answer.json')], directory, args.solve_limit)
    if process['status'] != 'completed':
        return {'status': 'unresolved', 'process': process}
    result = json.loads((directory / 'answer.json').read_text())
    if not samples.get('points') or not result.get('points'):
        return {**result, 'complete_sample_check': False, 'oracle': samples}
    reverse = native(args, directory / 'reverse', {'action': 'distances', 'recipe': case['recipe'],
                      'distance_width_mm': DISTANCE_WIDTH_MM, 'points': result.pop('points')})
    distances = reverse.get('distances', [])
    lo = [float(Fraction(d['lower_mm'])) for d in distances]
    hi = [float(Fraction(d['upper_mm'])) for d in distances if d.get('upper_mm') is not None]
    result['reverse'] = reverse
    result['distance_lower_mm'] = max([result['forward_lower_mm'], *lo])
    complete = (reverse.get('status') == 'answered' and len(hi) == len(distances)
                and bool(distances) and samples.get('completed') is True
                and samples.get('unresolved_candidates') == 0 and samples.get('failure') is None)
    result['complete_sample_check'] = complete
    result['distance_upper_mm'] = max([result['forward_upper_mm'], *hi]) if complete else None
    result['reference_samples'] = len(samples['points'])
    return result


def summary(report):
    result = {}
    for tolerance in TOLERANCES_MM:
        groups = {engine: Counter() for engine in ('catallax', 'occt', 'smlib')}
        wins = {engine: [] for engine in groups}
        for case in report['cases']:
            statuses = {engine: next(row['outcome'] for row in item['sweep']
                                    if row['tolerance_mm'] == tolerance)
                        for engine, item in case['engines'].items()}
            for engine, status in statuses.items():
                groups[engine][status] += 1
                if status == 'done-sampled' and any(s != 'done-sampled' for s in statuses.values()):
                    wins[engine].append(case['id'])
        result[str(tolerance)] = {'outcomes': {e: dict(c) for e, c in groups.items()},
                                 'completed_where_another_did_not': wins}
    return result


def run(args):
    args.out.mkdir(parents=True, exist_ok=False)
    selected = [case for case in corpus(args.measure_cases) if not args.case or case['id'] in args.case]
    if args.case and set(args.case) != {c['id'] for c in selected}:
        raise ValueError('unknown --case')
    versions = execute([str(args.python), '-B', str(HERE / 'versions.py')],
                       args.out / 'versions', 30)
    if versions['status'] != 'completed':
        raise RuntimeError('pinned Python runtimes unavailable')
    kernel_identity = json.loads(subprocess.check_output([str(args.binary), 'about'], text=True))
    kernel_revision = json.loads(kernel_identity['kernel'])['commit']
    if args.source_tag:
        tagged_revision = subprocess.check_output(
            ['git', '-C', str(ROOT), 'rev-parse', '--verify', f'refs/tags/{args.source_tag}^{{commit}}'],
            text=True).strip()
        if tagged_revision != kernel_revision:
            raise ValueError('producer binary identity does not match the requested source tag')
    report = {'format': 'catallax-kernel-comparison/1', 'status': 'running',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_revision': kernel_revision, 'source_tag': args.source_tag,
              'harness_revision': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
              'tools': json.loads((args.out / 'versions/stdout').read_text()),
              'binary_sha256': sha(args.binary), 'observer_sha256': sha(args.observer),
              'script_sha256': {p.name: sha(p) for p in HERE.glob('*.py')},
              'upstream_pins': json.loads((BENCH / 'tools.json').read_text()),
              'requirements_sha256': sha(BENCH / 'requirements.lock'),
              'kernel_identity': kernel_identity,
              'support_script_sha256': {str(p.relative_to(ROOT)): sha(p) for p in (
                  BENCH / 'targets.py', BENCH / 'python_tools.py', ROOT / 'tools/fillet-coverage/occt.py')},
              'machine': {'system': platform.platform(), 'node': platform.node(),
                          'cpu': platform.processor(), 'affinity': sorted(os.sched_getaffinity(0)),
                          'load_start': os.getloadavg(), 'class': 'shared development host'},
              'lscpu': json.loads(subprocess.check_output(['lscpu', '--json'], text=True)),
              'policy': {'tolerances_mm': TOLERANCES_MM, 'solve_limit_seconds': args.solve_limit,
                         'distance_width_mm': DISTANCE_WIDTH_MM,
                         'scope': 'same two-way finite boundary-sample check for every producer, not global Hausdorff certification',
                         'reference': 'independently compiled Catallax field definition; shares recipe semantics and authorship with Catallax B-rep',
                         'clock': 'whole producer process including startup, construction and STEP serialization; excludes observer',
                         'sweep': 'one unchanged solve evaluated against each threshold, no tolerance retuning or retries'},
              'inventory': [{k: v for k, v in case.items() if k not in ('recipe', 'fillet_case', 'benchmark_case')}
                            for case in selected],
              'case_sha256': {case['id']: hashlib.sha256(json.dumps(case, sort_keys=True).encode()).hexdigest()
                              for case in selected}, 'cases': []}
    save(args.out / 'results.json', report)
    for case in selected:
        directory = args.out / 'cases' / case['id']
        directory.mkdir(parents=True)
        samples = native(args, directory / 'oracle', {'action': 'samples', 'recipe': case['recipe']})
        row = {'id': case['id'], 'corpus': case['corpus'], 'engines': {}}
        for engine in ('catallax', 'occt', 'smlib'):
            dest = directory / engine
            before = os.getloadavg()
            answer, step = solve(args, case, engine, dest)
            checked = check(args, case, step, dest / 'check', samples) if answer['status'] == 'answered' else {}
            item = {'solve': answer, 'check': checked, 'sweep': sweep(answer, checked, args.solve_limit),
                    'load_before': before, 'load_after': os.getloadavg(),
                    'step_sha256': sha(step) if step.exists() else None}
            row['engines'][engine] = item
            print(case['id'], engine, answer['status'], item['sweep'][2]['outcome'], flush=True)
        report['cases'].append(row)
        save(args.out / 'results.json', report)
    report.update(status='completed', completed_utc=datetime.now(timezone.utc).isoformat(), summary=summary(report))
    save(args.out / 'results.json', report)
    print(json.dumps(report['summary'], indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--python', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--observer', type=Path, required=True)
    parser.add_argument('--measure-cases', type=Path, required=True)
    parser.add_argument('--source-tag', help='require the producer identity to match this release tag')
    parser.add_argument('--solve-limit', type=float, default=SOLVE_LIMIT_SECONDS)
    parser.add_argument('--case', action='append')
    parser.add_argument('--list', action='store_true')
    args = parser.parse_args()
    import resource
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    if args.list:
        print(json.dumps(corpus(args.measure_cases), indent=2))
    else:
        run(args)


if __name__ == '__main__':
    main()
