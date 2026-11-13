#!/usr/bin/env python3
"""M10 requests keep the 20M cap and publish no incomplete STEP result."""

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
REFUSALS = (
    *(['step.unsupported-geometry'] * 3),
    *(['numbers.operation-budget'] * 3),
    'step.unsupported-geometry',
    'numbers.geometric-gap-unproven',
    *(['intersect.parametric-coincidence'] * 3),
    *(['intersect.unsupported-pair'] * 3),
)


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


def write_questions(path):
    rows = cases()
    path.write_text(''.join(json.dumps({
        'format': 'catallax-question/1', 'source': 'recipe', 'problem': problem, 'max_work': ALLOWANCE,
    }) + '\n' for _, _, problem in rows))
    return rows


def comparable(response):
    response = copy.deepcopy(response)
    response.pop('wall_ms')
    response['manifest'].pop('stats', None)
    for operation in response['operation_log']['operations']:
        operation.pop('wall_ms')
    return response


class M10SplineBooleanTests(unittest.TestCase):
    def test_numbered_refusals_and_work_survive_worker_changes(self):
        with tempfile.TemporaryDirectory(prefix='catallax-m10-') as temporary:
            root = Path(temporary)
            rows = write_questions(root / 'questions.jsonl')
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
                for (number, name, _), case, expected in zip(rows, run['cases'], REFUSALS, strict=True):
                    with self.subTest(number=number, name=name, workers=workers):
                        response = case['response']
                        manifest = response['manifest']
                        self.assertFalse(case['cache_hit'])
                        self.assertEqual(manifest['status'], 'refused')
                        self.assertEqual(response['work']['limit'], ALLOWANCE)
                        self.assertGreater(response['work']['used'], 0)
                        self.assertLessEqual(response['work']['used'], ALLOWANCE)
                        code = manifest['refusal']['code']
                        self.assertEqual(code, expected)
                        self.assertTrue(manifest['refusal']['detail'])
                        self.assertEqual(response['work']['exhausted'], code == 'numbers.operation-budget')
                        self.assertFalse((Path(case['path']) / 'answer.step').exists())
                        self.assertIsNone(manifest['step'])
                        self.assertIsNone(manifest['result'])
                        operations = response['operation_log']['operations']
                        if name.startswith('side-') or name in ('through-cut', 'm7-cut'):
                            construction = [op for op in operations if op['id'] == 'cut'
                                            and op['operation'] in ('boolean', 'extrude_through')]
                            self.assertEqual(len(construction), 1)
                            self.assertEqual(construction[0]['status'], 'answered')
                            self.assertEqual(operations[-1]['operation'], 'step')
                            self.assertEqual(operations[-1]['refusal']['code'], expected)
                        if name.startswith('import-'):
                            imported = [op for op in operations if op['operation'] == 'import_step']
                            self.assertEqual(len(imported), 1)
                            self.assertEqual(imported[0]['status'], 'answered')
                        if name.startswith(('top-', 'import-', 'spline-pair-')):
                            self.assertEqual(operations[-1]['id'], 'cut')
                            self.assertEqual(operations[-1]['operation'], 'boolean')
                            self.assertEqual(operations[-1]['refusal']['code'], expected)
                        observed.append(comparable(response))
                        print(json.dumps({'case': number, 'name': name, 'workers': workers,
                                          'work': response['work'], 'refusal': manifest['refusal']}))
                if baseline is None:
                    baseline = observed
                else:
                    self.assertEqual(observed, baseline)
                # A refused replay must also remove a previous published file.
                first = Path(run['cases'][0]['path'])
                (first / 'answer.step').write_text('stale STEP result')
                replay = subprocess.run([
                    str(BINARY), '--threads', str(workers), 'replay', str(first / 'question.json'),
                    '--out', str(first),
                ], capture_output=True, text=True, check=True)
                self.assertEqual(comparable(json.loads(replay.stdout)), observed[0])
                self.assertFalse((first / 'answer.step').exists())


if __name__ == '__main__':
    unittest.main()
