"""Corpus accounting and a small end-to-end fillet coverage control."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

from test_loopless_faces import BINARY, ROOT

SPEC = importlib.util.spec_from_file_location('fillet_coverage', ROOT / 'tools/fillet-coverage/run.py')
coverage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(coverage)


class FilletCoverage(unittest.TestCase):
    def test_corpus_is_nonempty_and_keeps_schema_gaps_explicit(self):
        rows = list(coverage.cases())
        self.assertGreaterEqual(len(rows), 40)
        self.assertEqual(sum(row.get('recipe') is None for row in rows), 3)
        self.assertTrue(all(row.get('expected_stage') for row in rows))
        tags = {tag for row in rows for tag in row['tags']}
        self.assertTrue({'plane', 'cylinder', 'cone', 'sphere', 'torus', 'vertex-4',
                         'setback', 'roll-over', 'angle-limit', 'variable-radius'} <= tags)

    def test_refusal_needs_a_reason_and_errors_are_not_coverage(self):
        response = {'manifest': {'status': 'refused', 'refusal': {'code': 'fillet.unsupported-contact',
                                                               'detail': 'not constructed'}},
                    'wall_ms': 1}
        self.assertEqual(coverage.classify(response)['status'], 'refused-with-reason')
        response['manifest']['refusal']['detail'] = ''
        self.assertEqual(coverage.classify(response)['status'], 'incomplete')
        response['manifest']['status'] = 'internal-error'
        self.assertEqual(coverage.classify(response)['status'], 'incomplete')
        response['manifest'].update(status='refused', refusal={'code': 'cli.invalid-recipe',
                                                              'detail': 'unknown field'})
        self.assertEqual(coverage.classify(response)['status'], 'incomplete')

    def test_snapshot_does_not_hide_improvements_or_changed_refusals(self):
        expected = {'status': 'refused-with-reason', 'code': 'fillet.unsupported-corner'}
        self.assertTrue(coverage.matches(expected, expected))
        self.assertFalse(coverage.matches({'status': 'built'}, expected))
        self.assertFalse(coverage.matches({**expected, 'code': 'numbers.operation-budget'}, expected))

    def test_export_and_volume_failures_cannot_count_as_built(self):
        answered = {'status': 'built'}
        built = {'status': 'built', 'result': {'volume': 100}}
        self.assertEqual(coverage.assess(answered, built, None)['status'], 'incomplete')
        self.assertEqual(coverage.assess(answered, built, {'status': 'wrong'})['status'], 'wrong')
        different = {'status': 'built', 'result': {'volume': 99}}
        self.assertEqual(coverage.assess(answered, built, different)['code'], 'volume.disagreement')
        self.assertEqual(coverage.assess(answered, built, built)['status'], 'built')
        refused = {'status': 'refused-with-reason', 'code': 'occt.fillet-not-done'}
        self.assertEqual(coverage.assess(answered, refused, built)['status'], 'built')

    def test_budget_gaps_follow_the_final_requested_fillet(self):
        refusal = {'status': 'refused-with-reason', 'code': 'numbers.operation-budget'}
        for status, stage in [('answered', 'consumer'), ('refused', 'construction'),
                              ('interrupted', 'construction'), ('not-reached', 'upstream')]:
            final = {'status': status}
            with self.subTest(status=status):
                self.assertEqual(coverage.stage_result(refusal, final), stage + '-budget')
                unsupported = {**refusal, 'code': 'fillet.unsupported-contact'}
                self.assertEqual(coverage.stage_result(unsupported, final), stage + '-refusal')
        for status, expected in [('built', 'checked-export'), ('wrong', 'failed-check'),
                                  ('incomplete', 'incomplete')]:
            self.assertEqual(coverage.stage_result({'status': status}, {'status': 'answered'}), expected)
        schema = {**refusal, 'origin': 'corpus-not-kernel'}
        self.assertEqual(coverage.stage_result(schema, {'status': 'not-reached'}), 'schema-gap')

    def test_timeout_is_incomplete_not_refused(self):
        with tempfile.TemporaryDirectory(prefix='fillet-timeout-') as temporary:
            code = coverage.command([sys.executable, '-c', 'import time; time.sleep(10)'],
                                    Path(temporary) / 'log', 0.01)
            self.assertIsNone(code)

    def test_live_box_and_gap_controls(self):
        with tempfile.TemporaryDirectory(prefix='fillet-coverage-') as temporary:
            out = Path(temporary) / 'run'
            args = argparse.Namespace(binary=str(BINARY), out=str(out), occt_python=sys.executable,
                                      max_work=1_000_000, timeout=30, check=True,
                                      case=['box-single', 'torus-slice', 'variable-linear', 'narrow-cone',
                                            'roll-over-mitres'])
            self.assertEqual(coverage.run(args), 0)
            report = json.loads((out / 'scoreboard.json').read_text())
            rows = {r['id']: r for r in report['cases']}
            self.assertEqual(report['stage_counts'], {'checked-export': 2, 'consumer-budget': 1,
                                                     'schema-gap': 1, 'consumer-refusal': 1})
            markdown = (out / 'scoreboard.md').read_text()
            self.assertIn('| consumer-budget | 1 |', markdown)
            self.assertIn('| narrow-cone | consumer-refusal |', markdown)
            box = rows['box-single']
            self.assertEqual(box['catallax']['status'], 'built')
            expected = 4000 - 20*(1-math.pi/4)
            for result in (box['raw_step']['result'], box['occt']['result']):
                self.assertAlmostEqual(result['volume'], expected, delta=1e-6)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            self.assertEqual(rows['variable-linear']['catallax']['origin'], 'corpus-not-kernel')
            self.assertEqual(rows['variable-linear']['occt']['status'], 'built')
            self.assertEqual(rows['torus-slice']['occt']['status'], 'built')
            rounded = rows['roll-over-mitres']['raw_step']['result']
            # OCCT refuses this case, so check the rounded-box section integral.
            self.assertAlmostEqual(rounded['volume'], 3888 + 83*math.pi/3, delta=1e-6)
            self.assertEqual(len(rounded['witnesses']), 5)
            self.assertTrue(all(w['actual'] == w['expected'] for w in rounded['witnesses']))


if __name__ == '__main__':
    unittest.main()
