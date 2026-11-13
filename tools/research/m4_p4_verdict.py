#!/usr/bin/env python3
"""Collect P4 work and artifact evidence without changing its allowances."""
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
P4_STEP = 'f9a8b9a15a576a4989690f75e1bd9f697e052299b80b33c1bfb4fdc4634c7381'


def read(path):
    return json.loads(path.read_text())


def digest(data):
    return hashlib.sha256(data).hexdigest()


def without_clocks(value):
    if isinstance(value, dict):
        return {k: without_clocks(v) for k, v in value.items() if k != 'wall_ms'}
    if isinstance(value, list):
        return [without_clocks(v) for v in value]
    return value


def run(command, log, threads=1):
    with log.open('w') as output, log.with_suffix('.stderr').open('w') as errors:
        subprocess.run(command, cwd=ROOT, check=True, stdout=output,
                       stderr=errors,
                       env=dict(os.environ, CATALLAX_THREADS=str(threads)))
    result = read(log)
    if result.get('status') == 'internal-error':
        raise RuntimeError('kernel internal error; inspect the retained stderr log')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--profile', type=Path, required=True)
    parser.add_argument('--span', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    binary, profile, span = (p.resolve() for p in (args.binary, args.profile, args.span))
    recipe_file = ROOT / 'parts/P4-shelled-tray.json'
    recipe = read(recipe_file)
    edit = json.loads(json.dumps(recipe))
    next(s for s in edit['steps'] if s['id'] == 'tray')['thickness'] = '1.5'
    near = read(ROOT / 'parts/P0/near-miss.json')
    corpora = {
        'high': [(recipe, 100_000_000), (edit, 100_000_000),
                 (recipe, 20_000_000), (near, 1_000_000)],
        'low': [(recipe, 1_000_000), (edit, 1_000_000),
                (recipe, 100_000), (near, 1_000_000)],
    }
    doctor = run([sys.executable, '-B', str(DRIVER), 'doctor', '--binary', str(binary),
                  '--json'], out / 'doctor.json')
    assert doctor['about']['request'] == 'catallax-request/1'
    evidence = {
        'format': 'catallax-m4-p4-remaining/1',
        'source': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                          text=True).strip(),
        'binary_sha256': digest(binary.read_bytes()),
        'production_features': [],
        'recipe_sha256': digest(recipe_file.read_bytes()),
        'saved_answer_inventory_count': len(read(ROOT / 'tools/checks/earlier-answers.json')),
        'saved_answer_inventory_sha256': digest((ROOT / 'tools/checks/earlier-answers.json').read_bytes()),
        'production': [],
        'profiles': {},
        'raw_evidence': 'raw evidence retained privately',
    }
    for label, rows in corpora.items():
        corpus = out / f'{label}.jsonl'
        corpus.write_text(''.join(json.dumps({'format': 'catallax-question/1',
                          'source': 'recipe', 'problem': r, 'max_work': cap}) + '\n'
                          for r, cap in rows))
        baseline = None
        for threads in (1, 2, 4):
            directory = out / f'{label}-{threads}'
            run([sys.executable, '-B', str(DRIVER), 'run', str(corpus), '--binary',
                 str(binary), '--workers', str(threads), '--fresh', '--out',
                 str(directory), '--json'], out / f'{label}-{threads}.json', threads)
            report = read(directory / 'run.json')
            assert report['execution']['worker_exits'] == {str(i): 0 for i in range(threads)}
            responses, cases = [], []
            for index, case in enumerate(report['cases']):
                assert case['cache_hit'] is False
                case_dir = directory / case['path']
                response = read(case_dir / 'response.json')
                manifest = response['manifest']
                step = case_dir / 'answer.step'
                refused = label == 'low' and index < 3
                assert manifest['status'] == ('refused' if refused else 'answered')
                assert step.exists() is not refused
                if refused:
                    assert manifest['refusal']['code'] == 'numbers.operation-budget'
                elif label == 'high' and index in (0, 2):
                    assert digest(step.read_bytes()) == P4_STEP
                    assert manifest['exact'] is True
                responses.append(without_clocks(response))
                cases.append({'status': manifest['status'], 'refusal': manifest.get('refusal'),
                              'work': response['work'], 'operations': [{k: op[k] for k in ('id', 'operation', 'status', 'work', 'refusal')}
                                             for op in response['operation_log']['operations']],
                              'step_sha256': digest(step.read_bytes()) if step.exists() else None})
            if baseline is None:
                baseline = responses
            else:
                assert responses == baseline
                comparison = run([sys.executable, '-B', str(DRIVER), 'compare',
                                  str(out / f'{label}-1/run.json'), str(directory / 'run.json'),
                                  '--json'], out / f'{label}-1-v-{threads}.json')
                assert comparison['equal'] is True
            evidence['production'].append({
                'corpus': label, 'threads': threads, 'workers': threads,
                'worker_exits': report['execution']['worker_exits'],
                'responses_without_clocks_sha256': digest(json.dumps(responses, sort_keys=True).encode()),
                'cases': without_clocks(cases),
            })
    for cap in (1_000_000, 100_000_000):
        observed = run([str(profile), 'parts/P4-shelled-tray.json', '--out',
                        str(out / f'profile-{cap}'), '--max-work', str(cap)],
                       out / f'profile-{cap}.json')
        spanned = run([str(span), 'parts/P4-shelled-tray.json', '--out',
                       str(out / f'span-{cap}'), '--max-work', str(cap)],
                      out / f'span-{cap}.json')
        assert observed['summary']['unattributed'] == 0
        assert observed['summary']['work_used'] == spanned['work_used']
        assert spanned['observed_work_matches'] and spanned['paths_match_spans']
        production = evidence['production'][0 if cap == 100_000_000 else 3]['cases'][0]
        assert observed['summary']['work_used'] == production['work']['used']
        assert observed['status'] == production['status']
        callers = collections.defaultdict(lambda: [0, 0])
        for row in observed['rows']:
            key = (row['caller'], row['label'])
            callers[key][0] += row['calls']
            callers[key][1] += row['work']
        evidence['profiles'][str(cap)] = {
            'status': observed['status'], 'refusal': observed['refusal'],
            'summary': observed['summary'], 'stages': spanned['summary']['stages'],
            'all_charged_callers': [dict(caller=caller, label=label, calls=counts[0], work=counts[1])
                                    for (caller, label), counts in sorted(callers.items()) if counts[1]],
        }
    assert evidence['saved_answer_inventory_count'] == 477
    (out / 'portable.json').write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps(evidence, separators=(',', ':')))


if __name__ == '__main__':
    main()
