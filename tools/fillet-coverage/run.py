#!/usr/bin/env python3
"""Run versioned fillet recipes through the unchanged CLI and independent OCCT."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import os
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CORPUS = ROOT / 'parts/fillets/cases.jsonl'
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
FILLET_OPERATIONS = ('fillet', 'fillet_variable')


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def cases(path=CORPUS):
    result = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    names = [case['id'] for case in result]
    if len(names) != len(set(names)):
        raise ValueError('duplicate fillet case id')
    for case in result:
        if not case['id'].replace('-', '').isalnum():
            raise ValueError('case ids must be safe directory names')
        if not case['tags'] or not case['target']:
            raise ValueError('each case needs tags and a target outcome')
        recipe = case.get('recipe') or case.get('reference_recipe')
        if isinstance(recipe, str):
            recipe = json.loads((ROOT / recipe).read_text())
        if not recipe or not any(s['op'] in FILLET_OPERATIONS for s in recipe['steps']):
            raise ValueError(f"{case['id']}: missing fillet recipe")
        if case.get('recipe') is None and not case.get('schema_gap'):
            raise ValueError('unexpressible cases need a schema reason')
        analytic_volume(case)
        case = dict(case)
        case['resolved_recipe'] = recipe
        yield case


def classify(response):
    manifest = response['manifest']
    status = manifest['status']
    if status == 'answered':
        return {'status': 'built', 'exact': manifest['exact'],
                'marks': manifest.get('marks', []), 'step': manifest.get('step'),
                'result': manifest.get('result'), 'wall_ms': response['wall_ms']}
    refusal = manifest.get('refusal')
    if refusal and refusal.get('code') == 'cli.invalid-recipe':
        return {'status': 'incomplete', **refusal, 'wall_ms': response['wall_ms']}
    if status == 'refused' and refusal and refusal.get('code') and refusal.get('detail'):
        operations = response.get('operation_log', {}).get('operations', [])
        return {'status': 'refused-with-reason', **refusal, 'wall_ms': response['wall_ms'],
                'reached_fillet': any(o['operation'] in FILLET_OPERATIONS for o in operations)}
    return {'status': 'incomplete', 'code': status, 'wall_ms': response.get('wall_ms')}


def command(args, log, timeout):
    with log.open('w') as output:
        child = subprocess.Popen(args, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            return child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            return None
        finally:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()


def occt(case_path, recipe_path, out, args, step=None):
    report = out.with_suffix('.json')
    cmd = [args.occt_python, '-B', str(HERE / 'occt.py'), str(case_path),
           str(recipe_path), str(report)]
    if step:
        cmd.append(str(step))
    started = time.perf_counter()
    code = command(cmd, out.with_suffix('.log'), args.timeout)
    if code != 0 or not report.exists():
        return {'status': 'incomplete', 'code': 'harness.occt-process', 'exit': code,
                'wall_ms': (time.perf_counter()-started)*1000}
    return json.loads(report.read_text())


def analytic_volume(case):
    """An independently derived exact volume, which outranks a kernel's numerical answer."""
    reference = case.get('reference_volume')
    if reference is None:
        return None
    if not reference.get('method'):
        raise ValueError(f"{case['id']}: an analytic volume needs its derivation")
    return float(reference['volume'])


def assess(catallax, reference, imported, analytic=None):
    if catallax['status'] != 'built':
        return catallax
    if imported is None or imported['status'] == 'incomplete':
        return {**catallax, 'status': 'incomplete', 'code': 'harness.raw-step-unchecked'}
    if imported['status'] != 'built':
        return {**catallax, 'status': 'wrong', 'code': 'step.invalid'}
    if analytic is not None or reference['status'] == 'built':
        actual = imported['result']['volume']
        expected = analytic if analytic is not None else reference['result']['volume']
        allowance = max(1e-5, abs(expected)*1e-7)
        catallax = {**catallax, 'volume_difference': actual-expected,
                    'volume_comparison_allowance': allowance,
                    'volume_reference': 'analytic' if analytic is not None else 'occt'}
        if abs(actual-expected) > allowance:
            return {**catallax, 'status': 'wrong', 'code': 'volume.disagreement'}
        if analytic is not None and reference['status'] == 'built':
            catallax['occt_volume_difference'] = actual-reference['result']['volume']
            if any(w['actual'] != w['expected'] for w in reference.get('target_witnesses', [])):
                catallax['reference_comparison'] = 'definition-difference'
    return catallax


def stage_result(answer, final_fillet):
    """Locate the gap without counting an answered constructor as a checked export."""
    if answer.get('origin') == 'corpus-not-kernel':
        return 'schema-gap'
    if answer['status'] == 'incomplete':
        return 'incomplete'
    if answer['status'] == 'wrong':
        return 'failed-check'
    if answer['status'] == 'built':
        return 'checked-export'
    if final_fillet['status'] == 'answered':
        stage = 'consumer'
    elif final_fillet['status'] in ('refused', 'interrupted'):
        stage = 'construction'
    else:
        stage = 'upstream'
    return stage + ('-budget' if answer.get('code') == 'numbers.operation-budget' else '-refusal')


def matches(result, expected):
    return expected is None or (result['status'] == expected['status'] and
                               (expected.get('code') is None or
                                result.get('code') == expected['code']))


def run(args):
    out = Path(args.out).resolve()
    if out == ROOT or ROOT in out.parents:
        raise ValueError('keep proof outside the repository')
    out.mkdir(parents=True, exist_ok=False)
    selected = [c for c in cases() if not args.case or c['id'] in args.case]
    if args.case and set(args.case) != {c['id'] for c in selected}:
        raise ValueError('unknown case id')
    identity = hashlib.sha256(Path(args.binary).read_bytes()).hexdigest()
    source_head = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    source_dirty = bool(subprocess.check_output(
        ['git', '-C', str(ROOT), 'status', '--porcelain'], text=True).strip())
    report = {'format': 'catallax-fillet-coverage/1', 'status': 'running',
              'binary_sha256': identity, 'source_head': source_head, 'source_dirty': source_dirty,
              'corpus_sha256': hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
              'max_work': args.max_work, 'timeout_seconds': args.timeout,
              'command': sys.argv, 'load_average': os.getloadavg(),
              'nice': os.getpriority(os.PRIO_PROCESS, 0),
              'affinity': sorted(os.sched_getaffinity(0)), 'cases': []}
    for case in selected:
        dest = out / case['id']
        dest.mkdir()
        recipe_path, case_path = dest / 'recipe.json', dest / 'case.json'
        save(recipe_path, case['resolved_recipe'])
        save(case_path, {k: v for k, v in case.items() if k != 'resolved_recipe'})
        reference = occt(case_path, recipe_path, dest / 'occt', args)
        if analytic_volume(case) is not None and reference['status'] == 'built':
            # Record, without reclassifying, how far OCCT's own solid is from the exact value.
            reference['analytic_volume_difference'] = (reference['result']['volume']
                                                       - analytic_volume(case))
        imported, construction = None, []
        if case.get('recipe') is None:
            answer = {'status': 'refused-with-reason', 'code': 'schema.unavailable',
                      'detail': case['schema_gap'], 'wall_ms': None, 'origin': 'corpus-not-kernel'}
        else:
            exit_code = command([sys.executable, '-B', str(DRIVER), 'run', str(recipe_path),
                                 '--binary', str(Path(args.binary).resolve()), '--out', str(dest / 'cli'),
                                 '--max-work', str(args.max_work), '--source', 'recipe',
                                 '--workers', '1', '--fresh', '--json'], dest / 'cli.log', args.timeout)
            if exit_code != 0:
                answer = {'status': 'incomplete', 'code': 'harness.cli-process', 'exit': exit_code}
            else:
                record = json.loads((dest / 'cli/run.json').read_text())['cases'][0]
                answer = classify(record['response'])
                answer['work'] = record['response']['work']
                construction = [o for o in record['response']['operation_log']['operations']
                                if o['operation'] in FILLET_OPERATIONS]
                if answer['status'] == 'built':
                    imported = occt(case_path, recipe_path, dest / 'raw-step', args,
                                    Path(record['path']) / 'answer.step')
                    answer = assess(answer, reference, imported, analytic_volume(case))
        last_fillet = [s['id'] for s in case['resolved_recipe']['steps']
                       if s['op'] in FILLET_OPERATIONS][-1]
        final = next((o for o in reversed(construction) if o['id'] == last_fillet),
                     {'status': 'not-reached', 'wall_ms': None})
        expected = case.get('expected')
        stage = stage_result(answer, final)
        expected_stage = case.get('expected_stage')
        row = {'id': case['id'], 'tags': case['tags'], 'target': case['target'],
               'catallax': answer, 'occt': reference, 'raw_step': imported,
               'fillet_operations': construction,
               'final_fillet': final, 'stage_result': stage, 'expected_stage': expected_stage,
               'expected': expected, 'matches_snapshot': matches(answer, expected) and
               matches(reference, case.get('expected_occt')) and
               (expected_stage is None or stage == expected_stage)}
        report['cases'].append(row)
        report['counts'] = {engine: dict(Counter(r[engine]['status'] for r in report['cases']))
                            for engine in ('catallax', 'occt')}
        report['stage_counts'] = dict(Counter(r['stage_result'] for r in report['cases']))
        save(out / 'scoreboard.json', report)
        print(f"{case['id']}: Catallax {row['stage_result']} {answer.get('code', '')}; OCCT {reference['status']}",
              flush=True)
    if hashlib.sha256(Path(args.binary).read_bytes()).hexdigest() != identity:
        raise RuntimeError('binary changed during measurement')
    report['status'] = 'completed'
    save(out / 'scoreboard.json', report)
    markdown = ['| Catallax stage result | Cases |', '| --- | ---: |']
    markdown.extend(f'| {stage} | {count} |' for stage, count in report['stage_counts'].items())
    markdown.extend(['', '| Case | Catallax stage result | Catallax export | Work used / allowance | ms | Final fillet span | ms | OCCT | ms |',
                     '| --- | --- | --- | ---: | ---: | --- | ---: | --- | ---: |'])
    for row in report['cases']:
        a, b = row['catallax'], row['occt']
        def cell(r):
            outcome = r['status'] + (f" `{r['code']}`" if r.get('code') else '')
            return outcome + (f" ({r['reference_route']})" if r.get('reference_route') else '')
        def clock(r):
            return f"{r['wall_ms']:.2f}" if r.get('wall_ms') is not None else 'n/a'
        work = a.get('work')
        work_cell = f"{work['used']} / {work['limit']}" if work else 'n/a'
        markdown.append(f"| {row['id']} | {row['stage_result']} | {cell(a)} | {work_cell} | {clock(a)} | {row['final_fillet']['status']} | "
                        f"{clock(row['final_fillet'])} | {cell(b)} | {clock(b)} |")
    (out / 'scoreboard.md').write_text('\n'.join(markdown) + '\n')
    print(json.dumps({**report['counts'], 'stage_counts': report['stage_counts']}, sort_keys=True))
    return 1 if any(r[e]['status'] in ('wrong', 'incomplete') for r in report['cases']
                   for e in ('catallax', 'occt')) or (args.check and any(
                       not r['matches_snapshot'] for r in report['cases'])) else 0


def interrupt(signum, frame):
    raise KeyboardInterrupt


def main():
    signal.signal(signal.SIGTERM, interrupt)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--occt-python', default=sys.executable)
    parser.add_argument('--max-work', type=int, default=1_000_000)
    parser.add_argument('--timeout', type=float, default=30)
    parser.add_argument('--case', action='append')
    parser.add_argument('--check', action='store_true', help='also fail on a changed baseline outcome or failure stage')
    args = parser.parse_args()
    if args.max_work <= 0 or args.timeout <= 0:
        parser.error('work and timeout must be positive')
    sys.exit(run(args))


if __name__ == '__main__':
    main()
