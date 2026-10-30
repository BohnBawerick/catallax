#!/usr/bin/env python3
"""Exercise exported topology in OpenCASCADE, without shape repair."""
import copy
from fractions import Fraction
import json
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

# Required here, not an optional skip: CI must actually run the native checker.
import OCP  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
CAPTURE = ROOT / 'crates/catallax-cli/tests/fixtures/soliddraft-slot-cut.json'
HANG_LIMIT_SECONDS = 600


def bracket(gap, rotation=None, reverse=False):
    rotation = rotation or ['1', '0', '0', '0']
    shift = [gap, gap, Fraction(0)]
    if rotation == ['1', '0', '0', '1']:
        shift = [-gap, gap, Fraction(0)]
    elif rotation == ['1', '1', '0', '0']:
        shift = [gap, Fraction(0), gap]
    return {'format': 'catallax-question/1', 'source': 'recipe', 'max_work': 100_000_000,
            'problem': {'format': 'catallax-recipe/3', 'name': 'parallel ledge', 'steps': [
                {'id': 'base', 'op': 'block', 'size': ['12', '9', '1'],
                 'pose': {'rotation': rotation, 'translation': ['0', '0', '0']}},
                {'id': 'stand', 'op': 'block', 'size': [str(12 - 2 * gap), '1', '6'],
                 'pose': {'rotation': rotation, 'translation': list(map(str, shift))}},
                {'id': 'bracket', 'op': 'boolean', 'kind': 'union',
                 'target': 'stand' if reverse else 'base', 'tool': 'base' if reverse else 'stand'}],
                'result': 'bracket'}}


class StepImportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='catallax-step-import-')
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def export(self, question):
        source = self.directory / 'input.json'
        source.write_text(json.dumps(question))
        out = self.directory / 'run'
        result = subprocess.run([str(DRIVER), 'run', str(source), '--binary', str(BINARY),
                                 '--out', str(out), '--fresh', '--workers', '1', '--json'],
                                capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(result.stdout)['outcomes'], {'answered': 1})
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        manifest = case['response']['manifest']
        self.assertTrue(manifest['exact'], manifest)
        self.assertEqual(manifest['marks'], [])
        return Path(case['path']), manifest

    def assert_import(self, directory, faces):
        result = CHECK(str(directory / 'answer.step'), faces)
        self.assertEqual(result['status'], 'passed', result)
        return result

    def assert_regularization(self, directory, expected_gap):
        report = json.loads((directory / 'step-regularizations.json').read_text())
        self.assertEqual(len(report['regularizations']), 1)
        bound = Fraction(report['regularizations'][0]['max_deviation_mm'])
        self.assertEqual(bound, expected_gap)
        self.assertLessEqual(bound, Fraction(1, 1_000_000))
        self.assertGreaterEqual(Fraction(report['max_deviation_mm']), bound)
        return report

    def test_captured_slot_cut_and_earlier_bracket_prefixes(self):
        capture = json.loads(CAPTURE.read_text())
        gap = -Fraction(capture['problem']['steps'][2]['plane']['translation'][0])
        for result, model_faces in [('op_61', 11), ('op_89', 13), ('op_103', 14), ('op_140', 18)]:
            with self.subTest(result=result):
                question = copy.deepcopy(capture)
                steps = question['problem']['steps']
                index = next(i for i, step in enumerate(steps) if step['id'] == result)
                question['problem']['steps'] = steps[:index + 1]
                question['problem']['result'] = result
                directory, manifest = self.export(question)
                self.assertEqual(manifest['result']['faces'], model_faces)
                self.assert_import(directory, model_faces - 1)
                report = self.assert_regularization(directory, gap)
                self.assertEqual(report['faces'], model_faces - 1)
                if result == 'op_140':
                    # A rounded-away planar sliver remains after the slot cut.
                    # Partial cylinder charts must not rewrite an export whose
                    # other edges the raw translator cannot resolve.
                    self.assertEqual(manifest['step']['added_edges'], 0)

    def test_parallel_gaps_in_multiple_axes_and_operand_orders(self):
        gap = Fraction(1, 10**30)
        for rotation, reverse in [(['1', '0', '0', '0'], False),
                                  (['1', '0', '0', '1'], False),
                                  (['1', '1', '0', '0'], False),
                                  (['1', '0', '0', '0'], True)]:
            with self.subTest(rotation=rotation, reverse=reverse):
                directory, manifest = self.export(bracket(gap, rotation, reverse))
                self.assertEqual(manifest['result']['faces'], 11)
                result = self.assert_import(directory, 10)
                report = self.assert_regularization(directory, gap)
                self.assertEqual(report['faces'], 10)
                self.assertAlmostEqual(result['volume_mm3'], 168.0, delta=1e-9)

    def test_resolvable_gap_and_standalone_slot_stay_unregularized(self):
        directory, manifest = self.export(bracket(Fraction(1, 1000)))
        self.assertEqual(manifest['result']['faces'], 11)
        self.assertFalse((directory / 'step-regularizations.json').exists())
        self.assert_import(directory, 11)
        capture = json.loads(CAPTURE.read_text())
        capture['problem']['steps'] = capture['problem']['steps'][:-1]
        capture['problem']['result'] = 'op_112'
        directory, manifest = self.export(capture)
        self.assertEqual(manifest['result']['faces'], 6)
        self.assertFalse((directory / 'step-regularizations.json').exists())
        self.assert_import(directory, 6)

    def test_stream_removes_stale_regularization_evidence(self):
        directory = self.directory / 'stream'
        child = subprocess.Popen([str(BINARY), 'stream'], stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            for gap, work in [(Fraction(1, 10**30), 100_000_000),
                              (Fraction(1, 1000), 100_000_000), (Fraction(1, 1000), 0)]:
                question = bracket(gap)
                request = {'format': 'catallax-request/1', 'source': 'recipe',
                           'problem': question['problem'], 'max_work': work, 'out': str(directory)}
                child.stdin.write(json.dumps(request) + '\n')
                child.stdin.flush()
                response = json.loads(child.stdout.readline())
                if work == 0:
                    self.assertEqual(response['manifest']['refusal']['code'], 'numbers.operation-budget')
                    self.assertFalse((directory / 'answer.step').exists())
                    self.assertFalse((directory / 'step-regularizations.json').exists())
                elif gap < Fraction(1, 1_000_000):
                    self.assert_import(directory, 10)
                    self.assert_regularization(directory, gap)
                else:
                    self.assert_import(directory, 11)
                    self.assertFalse((directory / 'step-regularizations.json').exists())
            child.stdin.close()
            self.assertEqual(child.wait(timeout=HANG_LIMIT_SECONDS), 0, child.stderr.read())
        finally:
            if child.poll() is None:
                child.terminate()
                child.wait(timeout=HANG_LIMIT_SECONDS)
            child.stdout.close()
            child.stderr.close()



if __name__ == '__main__':
    unittest.main()
