#!/usr/bin/env python3
"""Original M10 STEP answers, with slow top cuts at the 100M application allowance."""

import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
ALLOWANCE = 20_000_000
APPLICATION_ALLOWANCE = 100_000_000
TOP_CASES = frozenset(['M10-004', 'M10-005', 'M10-006'])


def allowance(number):
    return APPLICATION_ALLOWANCE if number in TOP_CASES else ALLOWANCE


def pose(at):
    return {'rotation': ['1', '0', '0', '0'], 'translation': at}


def loft_steps(prefix='', x='0'):
    fixture = ROOT / 'crates/catallax-features/tests/fixtures/m10-smooth-loft.json'
    steps = json.loads(fixture.read_text())['steps']
    for step in steps:
        step['id'] = prefix + step['id']
        if step['op'] == 'sketch':
            step['plane']['translation'][0] = x
        else:
            for section in step['sections']:
                section['profile'] = prefix + section['profile']
    return steps


def recipe(name, steps, result='cut'):
    return {'format': 'catallax-recipe/5', 'name': name, 'steps': steps, 'result': result}


def cases():
    """Stable case numbers identify refusals without changing kernel codes."""
    rows = []

    def add(name, problem):
        rows.append((f'M10-{len(rows) + 1:03d}', name, {'format': 'catallax-step/1', 'recipe': problem}))

    for case, size, at in [
        ('side', ['2', '10', '1'], ['1', '1', '1/2']),
        ('top', ['10', '10', '1'], ['-3', '5/2', '1/2']),
    ]:
        for kind in ['subtract', 'union', 'intersect']:
            name = f'{case}-{kind}'
            add(name, recipe(name, loft_steps() + [
                {'id': 'box', 'op': 'block', 'size': size, 'pose': pose(at)},
                {'id': 'cut', 'op': 'boolean', 'kind': kind, 'target': 'loft', 'tool': 'box'},
            ]))

    # The rectangle lies strictly inside both end profiles, so its sweep
    # removes material through both caps of the loft.
    hole = {
        'id': 'hole', 'op': 'sketch', 'plane': pose(['3/2', '1/2', '-1']),
        'loops': [[{'kind': 'line', 'id': f'h{i}', 'to': point} for i, point in enumerate(
            [['1', '0'], ['1', '1/2'], ['0', '1/2'], ['0', '0']])]],
    }
    add('through-cut', recipe('through-cut', loft_steps() + [hole, {
        'id': 'cut', 'op': 'extrude_through', 'target': 'loft', 'profile': 'hole',
    }]))

    m7 = json.loads((ROOT / 'crates/catallax-features/tests/fixtures/loft-spline-prism.json').read_text())
    m7['name'] = 'm7-cut'
    m7['steps'] += [
        {'id': 'box', 'op': 'block', 'size': ['2', '4', '1'], 'pose': pose(['1', '1', '1'])},
        {'id': 'cut', 'op': 'boolean', 'kind': 'subtract', 'target': 'loft', 'tool': 'box'},
    ]
    m7['result'] = 'cut'
    add('m7-cut', m7)

    # OpenCASCADE's independently written ruled loft has four curved
    # bilinear B-spline sides. The supplied bytes, not a path, reach import.
    source = (ROOT / 'crates/catallax-import/tests/fixtures/m27-loft-occt.step').read_bytes()
    digest = hashlib.sha256(source).hexdigest()
    for kind in ['subtract', 'union', 'intersect']:
        name = f'import-{kind}'
        part = {'id': 'part', 'op': 'import_step', 'file': {'name': 'm27-loft-occt.step', 'sha256': digest},
                'pose': pose(['0', '0', '0'])}
        problem = recipe(name, [
            part, {'id': 'box', 'op': 'block', 'size': ['4', '4', '1'], 'pose': pose(['0', '0', '1/2'])},
            {'id': 'cut', 'op': 'boolean', 'kind': kind, 'target': 'part', 'tool': 'box'},
        ])
        add(name, {'format': 'catallax-import-recipe/1', 'recipe': problem, 'files': {digest: source.decode()}})

    for kind in ['subtract', 'union', 'intersect']:
        name = f'spline-pair-{kind}'
        add(name, recipe(name, loft_steps() + loft_steps('other_', '1') + [
            {'id': 'cut', 'op': 'boolean', 'kind': kind, 'target': 'loft', 'tool': 'other_loft'},
        ]))
    return rows


def write_questions(path, rows=None):
    rows = cases() if rows is None else rows
    path.write_text(''.join(json.dumps({
        'format': 'catallax-question/1', 'source': 'recipe', 'problem': problem, 'max_work': allowance(number),
    }) + '\n' for number, _, problem in rows))
    return rows


def comparable(response):
    response = copy.deepcopy(response)
    response.pop('wall_ms')
    response['manifest'].pop('stats', None)
    for operation in response['operation_log']['operations']:
        operation.pop('wall_ms')
    return response


class M10SplineBooleanTests(unittest.TestCase):
    def test_numbered_outcomes_and_work_survive_worker_changes(self):
        self.check_outcomes([row for row in cases() if row[0] not in TOP_CASES])

    @unittest.skipUnless(os.environ.get('CATALLAX_SLOW_LOFT_TOP') == '1',
                         'slow loft top cuts; run explicitly for the phase-5 work list')
    def test_slow_top_requests_reach_step(self):
        self.check_outcomes([row for row in cases() if row[0] in TOP_CASES])

    def check_outcomes(self, rows):
        with tempfile.TemporaryDirectory(prefix='catallax-m10-') as temporary:
            root = Path(temporary)
            write_questions(root / 'questions.jsonl', rows)
            baseline = None
            for workers in [1, 4]:
                out = root / f'workers-{workers}'
                environment = dict(os.environ, CATALLAX_THREADS=str(workers))
                completed = subprocess.run([
                    str(DRIVER), 'run', str(root / 'questions.jsonl'), '--binary', str(BINARY),
                    '--out', str(out), '--workers', str(workers), '--fresh', '--json',
                ], capture_output=True, text=True, env=environment, check=True)
                self.assertEqual(json.loads(completed.stdout)['crashes'], 0)
                run = json.loads((out / 'run.json').read_text())
                self.assertEqual(len(run['cases']), len(rows))
                observed = []
                for (number, name, _), case in zip(rows, run['cases'], strict=True):
                    with self.subTest(number=number, name=name, workers=workers):
                        response = case['response']
                        manifest = response['manifest']
                        self.assertFalse(case['cache_hit'])
                        self.assertEqual(response['work']['limit'], allowance(number))
                        self.assertGreater(response['work']['used'], 0)
                        self.assertLessEqual(response['work']['used'], allowance(number))
                        operations = response['operation_log']['operations']
                        self.assert_approximate(number, case, manifest, operations)
                        observed.append(comparable(response))
                if baseline is None:
                    baseline = observed
                else:
                    self.assertEqual(observed, baseline)
                # A replay must replace a stale file with the checked STEP bytes.
                numbers = [number for number, _, _ in rows]
                if 'M10-008' not in numbers:
                    continue
                index = numbers.index('M10-008')
                first = Path(run['cases'][index]['path'])
                (first / 'answer.step').write_text('stale STEP result')
                replay = subprocess.run([
                    str(BINARY), '--threads', str(workers), 'replay', str(first / 'question.json'),
                    '--out', str(first),
                ], capture_output=True, text=True, check=True)
                self.assertEqual(comparable(json.loads(replay.stdout)), observed[index])
                self.assertEqual(hashlib.sha256((first / 'answer.step').read_bytes()).hexdigest(),
                                 observed[index]['manifest']['step']['sha256'])

    def assert_approximate(self, number, case, manifest, operations):
        """Check the actual file and the numerical stages that produced it."""
        self.assertEqual(manifest['status'], 'answered')
        self.assertIsNone(manifest['refusal'])
        self.assertFalse(manifest['exact'])
        self.assertEqual([mark['code'] for mark in manifest['marks']], ['step.numerical-approximate'])
        step = manifest['step']
        self.assertTrue(step['approximate'])
        self.assertEqual(step['category'], 'approximate')
        self.assertEqual(step['tolerance_mm'], '1/10000000')
        self.assertEqual(step['complete_body_bound'], 'unestablished')
        data = (Path(case['path']) / 'answer.step').read_bytes()
        self.assertEqual(hashlib.sha256(data).hexdigest(), step['sha256'])
        self.assertIn(b'complete_body_bound=unestablished', data)
        if number in TOP_CASES:
            self.assertEqual(manifest['proof'], 'numerical')
            self.assertEqual(manifest['category'], 'approximate')
            selected = manifest['numerical']['work_stages']
            boolean = [op for op in operations if op['operation'] == 'numerical-boolean']
            self.assertEqual([op['status'] for op in boolean], ['answered'])
            self.assertEqual(sum(selected.values()), boolean[0]['work']['used'])
            for stage in ['numerical-step-definition', 'numerical-native-admission',
                          'numerical-step-payload', 'numerical-step-reread']:
                self.assertEqual([op['status'] for op in operations if op['operation'] == stage],
                                 ['answered'])
            self.assertTrue(case['response']['work']['over_target'])
            self.assertFalse(case['response']['work']['exhausted'])
            self.assertEqual(case['response']['work']['target'], ALLOWANCE)
            self.assertFalse(any(op['operation'] == 'step' for op in operations))
            return
        construction = [op for op in operations if op['id'] == 'cut'
                        and op['operation'] in ('boolean', 'extrude_through')]
        self.assertEqual([op['status'] for op in construction], ['answered'])
        exact = [op for op in operations if op['operation'] == 'step']
        self.assertEqual(len(exact), 1)
        expected = ('numbers.geometric-gap-unproven' if 'M10-008' <= number <= 'M10-011'
                    else 'step.unsupported-geometry')
        self.assertEqual(exact[0]['refusal']['code'], expected)
        self.assertEqual([op['operation'] for op in operations[-1:]], ['numerical-step-reread'])
        numerical = [op for op in operations if op['operation'] == 'numerical-step']
        self.assertEqual([op['status'] for op in numerical], ['answered'])


if __name__ == '__main__':
    unittest.main()
