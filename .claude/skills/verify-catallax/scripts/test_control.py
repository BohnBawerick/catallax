#!/usr/bin/env python3
"""Real-command regression tests for the verification driver."""
import copy
import json
import os
from fractions import Fraction
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[4]
DRIVER = Path(__file__).with_name('control-catallax.py')
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
EXAMPLES = ROOT / '.claude/skills/verify-catallax/references/examples'


def load_driver():
    import importlib.util
    spec = importlib.util.spec_from_file_location('control_catallax', DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fake_case(driver, wall_ms, status='answered', refusal=None, preview=None, source='recipe'):
    response = {'wall_ms': wall_ms, 'work': 1, 'manifest': {'status': status, 'refusal': refusal}}
    if preview is not None:
        response['preview'] = preview
    return {'key': 'k', 'id': 'part', 'cache_hit': False, 'class': source, 'group': 'g',
            'target_ms': driver.TARGETS[source], 'response': response}


class AcceptanceTests(unittest.TestCase):
    """Preview wait and proof completion are separate; no BINARY needed."""

    def setUp(self):
        self.driver = load_driver()

    def test_command_return_without_preview_evidence_is_not_a_preview_pass(self):
        verdict = self.driver.acceptance(fake_case(self.driver, 5)['response'])
        self.assertEqual(verdict['preview']['status'], 'not-measured')
        self.assertNotIn('met', verdict['preview'])
        self.assertTrue(verdict['proof']['met'])

    def test_preview_needs_unproved_label_and_delivered_display(self):
        for record in ({'state': 'proved', 'displayed': True, 'wall_ms': 1},
                       {'state': 'unproved', 'displayed': False, 'wall_ms': 1},
                       {'state': 'unproved', 'displayed': True},
                       {'state': 'unproved', 'displayed': True, 'wall_ms': True}):
            verdict = self.driver.acceptance(fake_case(self.driver, 5, preview=record)['response'])
            self.assertEqual(verdict['preview']['status'], 'invalid', record)
            self.assertFalse(verdict['preview'].get('met', False))
        good = self.driver.acceptance(fake_case(self.driver, 5000, preview={'state': 'unproved', 'displayed': True, 'wall_ms': 40})['response'])
        self.assertTrue(good['preview']['met'])
        self.assertFalse(good['proof']['met'], 'slow proof stays visible beside a fast preview')
        self.assertFalse(good['built'])

    def test_preview_alone_never_builds_a_part_and_a_refusal_is_incomplete(self):
        preview = {'state': 'unproved', 'displayed': True, 'wall_ms': 1}
        refused = fake_case(self.driver, 1, status='refused', refusal={'code': 'numbers.operation-budget'}, preview=preview)
        verdict = self.driver.acceptance(refused['response'])
        self.assertTrue(verdict['preview']['met'])
        self.assertEqual(verdict['proof']['status'], 'incomplete')
        self.assertFalse(verdict['built'])

    def test_summary_counts_the_columns_separately(self):
        preview = {'state': 'unproved', 'displayed': True, 'wall_ms': 900}
        report = {'elapsed_ms': 1, 'cases': [
            fake_case(self.driver, 1500, preview=preview),
            fake_case(self.driver, 2500),
            fake_case(self.driver, 5, status='refused', refusal={'code': 'x'}),
            fake_case(self.driver, 10, source='compare')]}
        accepted = self.driver.summarize(report)['acceptance']
        self.assertEqual(accepted['preview'], {'measured': 1, 'met': 1, 'not_measured': 3, 'invalid': 0, 'target_ms': 1000})
        self.assertEqual(accepted['proof'], {'completed': 3, 'within_target': 2, 'incomplete': 1, 'target_ms': 2000})
        self.assertEqual(accepted['built_after_proof'], 2)
        self.assertEqual(accepted['comparison_requested'], 1)

    def test_comparison_column_is_requested_not_default(self):
        import argparse
        base = dict(print_tolerance=None)
        self.assertEqual(self.driver.parts_sources(argparse.Namespace(compare=False, **base)), ('recipe', 'field'))
        self.assertEqual(self.driver.parts_sources(argparse.Namespace(compare=True, **base)), ('recipe', 'field', 'compare'))
        self.assertEqual(self.driver.parts_sources(argparse.Namespace(compare=False, print_tolerance='0.5')), ('recipe', 'field', 'print'))


class WorkerShutdownTests(unittest.TestCase):
    """Replay requires exit code 0 from every worker, so shutdown speed must not decide it."""

    def setUp(self):
        self.driver = load_driver()
        self.root = tempfile.TemporaryDirectory()
        self.addCleanup(self.root.cleanup)

    def worker(self, after_input):
        path = Path(self.root.name) / 'worker'
        path.write_text(f'#!/usr/bin/env python3\nimport sys, time\nsys.stdin.read()\n{after_input}\n')
        path.chmod(0o755)
        return self.driver.Worker(path)

    def test_slow_but_clean_exit_counts_as_clean(self):
        self.assertEqual(self.worker('time.sleep(6)').close(), 0)

    def test_worker_that_never_exits_is_still_stopped(self):
        self.driver.HANG_LIMIT_SECONDS = 0.5
        self.assertNotEqual(self.worker('time.sleep(600)').close(), 0)


@unittest.skipUnless(BINARY.exists(), 'build catallax-verify first')
class DriverTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='catallax-driver-')
        self.out = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def command(self, *args, code=0):
        completed = subprocess.run([str(DRIVER), *map(str, args), '--json'], capture_output=True, text=True)
        self.assertEqual(completed.returncode, code, completed.stderr + completed.stdout)
        return json.loads(completed.stdout)

    def run_case(self, name, *flags):
        return self.command('run', EXAMPLES / 'block.json', '--binary', BINARY, '--out', self.out / name, *flags)

    def test_cache_resume_corruption_and_step_diff(self):
        cache = self.out / 'cache'
        first = self.run_case('first', '--cache', cache)
        self.assertEqual(first['cache_hits'], 0)
        self.assertEqual(self.run_case('first')['cache_hits'], 1)
        self.assertEqual(self.run_case('second', '--cache', cache)['cache_hits'], 1)
        self.command('compare', self.out / 'first/run.json', self.out / 'second/run.json')
        second = json.loads((self.out / 'second/run.json').read_text())
        step = Path(second['cases'][0]['path']) / 'answer.step'
        step.write_bytes(step.read_bytes() + b'changed')
        diff = self.command('compare', self.out / 'first/run.json', self.out / 'second/run.json', code=1)
        self.assertEqual(diff['failures'][0]['cause'], 'step-bytes')
        self.assertEqual(self.run_case('second')['cache_hits'], 0)
        self.command('compare', self.out / 'first/run.json', self.out / 'second/run.json')

    def test_dry_run_does_not_create_artifacts(self):
        result = self.run_case('dry', '--fresh', '--dry-run')
        self.assertTrue(result['dry_run'])
        self.assertFalse((self.out / 'dry').exists())

    def test_stream_removes_stale_step_on_refusal(self):
        directory = self.out / 'same-output'
        question = json.loads((EXAMPLES / 'block.json').read_text())
        requests = [{'format': 'catallax-request/1', 'source': 'question', 'problem': question, 'out': str(directory), 'max_work': work} for work in (100_000_000, 0)]
        result = subprocess.run([str(BINARY), 'stream'], input='\n'.join(json.dumps(r) for r in requests) + '\n', capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        responses = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertEqual(responses[0]['manifest']['status'], 'answered')
        self.assertEqual(responses[1]['manifest']['refusal']['code'], 'numbers.operation-budget')
        self.assertFalse((directory / 'answer.step').exists())
        self.assertEqual(json.loads((directory / 'manifest.json').read_text())['status'], 'refused')

    def test_allowance_is_part_of_cache_key(self):
        cache = self.out / 'cache'
        self.run_case('full', '--cache', cache)
        limited = self.run_case('limited', '--cache', cache, '--max-work', 0)
        self.assertEqual(limited['cache_hits'], 0)
        self.assertEqual(limited['outcomes'], {'numbers.operation-budget': 1})
        self.assertFalse(list((self.out / 'limited').rglob('answer.step')))

    def test_worker_counts_duplicate_order_and_cleanup(self):
        corpus = self.out / 'corpus.jsonl'
        question = json.loads((EXAMPLES / 'block.json').read_text())
        rows = []
        for index in range(8):
            rows.append(dict(question, id=f'case-{index}'))
        rows.append(rows[0])
        corpus.write_text('\n'.join(json.dumps(q) for q in rows))
        for workers in (1, 2, 4):
            result = self.command('run', corpus, '--binary', BINARY, '--workers', workers, '--fresh', '--out', self.out / str(workers))
            self.assertEqual(result['cases'], 9)
        for workers in (2, 4):
            self.command('compare', self.out / '1/run.json', self.out / f'{workers}/run.json')
        self.command('cleanup', '--out', self.out, '--dry-run')
        self.command('cleanup', '--out', self.out)
        self.assertTrue((self.out / '1/run.json').exists())

    def test_certified_original_target_and_artifacts_agree_across_workers(self):
        original = json.loads((EXAMPLES / 'certified-recipe.json').read_text())
        second = copy.deepcopy(original)
        second['recipe']['steps'].append(
            {'id': 'second', 'op': 'offset', 'body': 'offset', 'distance': '-1/4'})
        second['recipe']['result'] = 'second'
        trimmed = copy.deepcopy(original)
        trimmed['recipe']['steps'].extend([
            {'id': 'trim_window', 'op': 'block', 'size': ['8', '8', '1'],
             'pose': {'rotation': ['3', '0', '1', '0'], 'translation': ['-7/5', '-4', '24/5']}},
            {'id': 'trimmed', 'op': 'boolean', 'kind': 'intersect',
             'target': 'offset', 'tool': 'trim_window'},
            {'id': 'after_trim', 'op': 'offset', 'body': 'trimmed', 'distance': '-1/4'},
        ])
        trimmed['recipe']['result'] = 'after_trim'
        trimmed['field_points'][0] = ['21/10', '0', '14/5']
        empty = copy.deepcopy(original)
        empty['recipe']['steps'][-1]['distance'] = '-2'
        inputs = self.out / 'certified-inputs'
        inputs.mkdir()
        for name, problem in [('offset', original), ('second', second), ('trimmed', trimmed), ('empty', empty)]:
            (inputs / f'{name}.json').write_text(json.dumps(problem))
        for workers in (1, 2, 4):
            output = self.out / f'certified-{workers}'
            result = self.command('run', inputs, '--binary', BINARY, '--workers', workers,
                                  '--fresh', '--out', output)
            self.assertEqual(result['outcomes'], {'answered': 4})
            run = json.loads((output / 'run.json').read_text())
            for case in run['cases']:
                response = case['response']
                manifest = response['manifest']
                self.assertFalse(response['work']['exhausted'])
                self.assertEqual(manifest['format'], 'catallax-manifest/4')
                for consumer in ('display', 'properties', 'field', 'step', 'edit'):
                    self.assertEqual(manifest[consumer]['status'], 'answered', manifest[consumer])
                point = manifest['field']['points'][0]
                target = manifest['certificate']['retained_target']
                if target['empty']:
                    self.assertEqual(manifest['result']['kind'], 'empty')
                    self.assertEqual(point['location'], 'outside')
                    self.assertEqual(point['value']['mm'], {'lo': 'infinity', 'hi': 'infinity'})
                else:
                    self.assertEqual(manifest['result']['faces_by_kind'], {'bspline': 1, 'plane': 2})
                    self.assertEqual((manifest['result']['edges'], manifest['result']['vertices']), (2, 0))
                    self.assertEqual(point['location'], 'inside')
                    expected = (Fraction(-1, 4) if Fraction(point['point'][0]) == Fraction(21, 10) else
                                Fraction(-8, 5) + Fraction(target['distance_mm']))
                    interval = point['value']['mm']
                    self.assertLessEqual(Fraction(interval['lo']), expected)
                    self.assertGreaterEqual(Fraction(interval['hi']), expected)
                    expected_refs = {('face', item['face']) for item in manifest['names']['faces']}
                    expected_refs.update(('edge', item['edge']) for item in manifest['names']['edges']
                                         if item['name'] is not None)
                    actual_refs = {(item['kind'], item['old_handle']) for item in manifest['edit']['resolved']}
                    self.assertEqual(actual_refs, expected_refs)
                    self.assertTrue(all(match['status'] == 'resolved' for match in manifest['edit']['resolved']))
            if workers != 1:
                self.command('compare', self.out / 'certified-1/run.json', output / 'run.json')

    def test_interruption_resumes_complete_artifacts(self):
        corpus = self.out / 'resume.jsonl'
        question = json.loads((EXAMPLES / 'block.json').read_text())
        corpus.write_text('\n'.join(json.dumps(dict(question, id=f'resume-{i}')) for i in range(500)))
        destination = self.out / 'resume'
        args = ['run', str(corpus), '--binary', str(BINARY), '--out', str(destination), '--json']
        child = subprocess.Popen([str(DRIVER), *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        limit = load_driver().HANG_LIMIT_SECONDS
        try:
            deadline = time.monotonic() + limit
            checkpoint = destination / 'checkpoint-0.jsonl'
            while not checkpoint.exists():
                if child.poll() is not None or time.monotonic() > deadline:
                    self.fail('worker did not publish a checkpoint')
                time.sleep(0.01)
            child.terminate()
            child.communicate(timeout=limit)
            self.assertNotEqual(child.returncode, 0)
        finally:
            if child.poll() is None:
                child.terminate()
                child.communicate(timeout=limit)
        result = self.command(*args[:-1])
        self.assertEqual(result['cases'], 500)
        self.assertGreater(result['cache_hits'], 0)

    def test_recipe_comparison_and_doctor(self):
        info = self.command('doctor', '--binary', BINARY)
        self.assertEqual(info['about']['request'], 'catallax-request/1')
        for source in ('recipe', 'compare'):
            self.command('run', EXAMPLES / 'recipe.json', '--source', source, '--binary', BINARY, '--out', self.out / source)
        report = json.loads((self.out / 'compare/run.json').read_text())
        self.assertEqual(report['cases'][0]['response']['manifest']['comparison']['verdict'], 'agree')
        self.assertFalse(list((self.out / 'compare').rglob('answer.step')))

    def test_field_display_comparison_checks_obj_bytes_but_not_clock(self):
        for name in ('field-a', 'field-b'):
            result = self.command('run', EXAMPLES / 'recipe.json', '--source', 'field',
                                  '--binary', BINARY, '--out', self.out / name)
            self.assertEqual(result['outcomes'], {'answered': 1})
        report = json.loads((self.out / 'field-b/run.json').read_text())
        case = report['cases'][0]
        directory = Path(case['path'])
        display = case['response']['manifest']['face_measurements']['display']
        self.assertTrue((directory / display['file']).exists())
        self.assertFalse((directory / 'answer.step').exists())
        self.command('compare', self.out / 'field-a/run.json', self.out / 'field-b/run.json')
        obj = directory / display['file']
        obj.write_bytes(obj.read_bytes() + b'# changed display\\n')
        diff = self.command('compare', self.out / 'field-a/run.json', self.out / 'field-b/run.json', code=1)
        self.assertEqual(diff['failures'][0]['cause'], 'display-bytes')

    def test_exported_question_keeps_its_allowance_in_batch_replay(self):
        question = json.loads((EXAMPLES / 'block.json').read_text())
        exported = self.out / 'question.json'
        exported.write_text(json.dumps({'format': 'catallax-question/1', 'source': 'question',
                                        'problem': question, 'max_work': 0}))
        result = self.command('run', exported, '--binary', BINARY, '--out', self.out / 'replay')
        self.assertEqual(result['outcomes'], {'numbers.operation-budget': 1})
        self.assertFalse(list((self.out / 'replay').rglob('answer.step')))

    def test_scoreboard_rebuilds_after_dimension_edit(self):
        targets = self.out / 'targets'
        targets.mkdir()
        recipe = json.loads((EXAMPLES / 'recipe.json').read_text())
        (targets / 'part.json').write_text(json.dumps(recipe))
        specification = {'format': 'catallax-part-edits/1', 'cases': {
            'part.json': [{'step': 'box', 'path': ['size', 2], 'value': '7'}]}}
        (targets / 'edit-checks.json').write_text(json.dumps(specification))
        result = self.command('parts', targets, '--binary', BINARY, '--out', self.out / 'edited')
        check = result['cases'][0]['edit_check']
        self.assertEqual(check['status'], 'passed', check)
        self.assertEqual(check['brep']['faces_checked'], 6)
        self.assertEqual(check['brep']['edges_checked'], 12)

    def test_scoreboard_runs_comparison_only_when_requested(self):
        targets = self.out / 'targets'
        targets.mkdir()
        (targets / 'part.json').write_text((EXAMPLES / 'recipe.json').read_text())
        plain = self.command('parts', targets, '--binary', BINARY, '--out', self.out / 'plain')
        self.assertEqual(plain['columns'], ['recipe', 'field'])
        self.assertEqual(plain['comparison'], {'requested': False, 'status': 'not-requested'})
        self.assertIsNone(plain['comparison_disagreements'])
        self.assertNotIn('compare', plain['cases'][0])
        self.assertFalse((self.out / 'plain/cases/part/compare').exists())
        row = plain['cases'][0]['acceptance']['recipe']
        self.assertEqual(row['preview']['status'], 'not-measured')
        asked = self.command('parts', targets, '--binary', BINARY, '--out', self.out / 'asked', '--compare')
        self.assertEqual(asked['columns'], ['recipe', 'field', 'compare'])
        self.assertEqual(asked['comparison_disagreements'], 0)
        self.assertEqual(asked['cases'][0]['compare']['manifest']['comparison']['verdict'], 'agree')


if __name__ == '__main__':
    unittest.main()
