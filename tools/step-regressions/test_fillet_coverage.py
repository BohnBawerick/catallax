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

    def test_an_analytic_volume_outranks_the_reference_kernel(self):
        answered = {'status': 'built'}
        occt = {'status': 'built', 'result': {'volume': 100.01}}
        exact = {'status': 'built', 'result': {'volume': 100}}
        result = coverage.assess(answered, occt, exact, 100.0)
        self.assertEqual((result['status'], result['volume_reference']), ('built', 'analytic'))
        self.assertEqual(coverage.assess(answered, occt, occt, 100.0)['code'], 'volume.disagreement')
        refused = {'status': 'refused-with-reason', 'code': 'occt.fillet-not-done'}
        self.assertEqual(coverage.assess(answered, refused, occt, 100.0)['code'],
                         'volume.disagreement')
        with self.assertRaises(ValueError):
            coverage.analytic_volume({'id': 'x', 'reference_volume': {'volume': '1'}})

    def test_oblique_cylinder_volume_is_its_steiner_value(self):
        # The selected rim is convex and the cut cylinder is convex, so the rounded
        # rim is the cut cylinder's opening by the unit ball, apart from the
        # bottom rim, which the opening also rounds. With E the cut cylinder
        # moved one unit inward, Steiner gives V(E) + S(E) + M(E) + 4 pi / 3.
        row = next(r for r in coverage.cases() if r['id'] == 'oblique-cylinder')
        s, c = 336/625, 527/625
        rho, height = 4, 7 - 1/c
        n = 4096
        # The rim's half exterior angle times its length element, periodic and smooth.
        rim = sum(math.acos(s*math.cos(t)) * rho * math.hypot(1, s/c*math.sin(t))
                  for t in (2*math.pi*k/n for k in range(n))) * 2*math.pi/n / 2
        volume = (math.pi*rho*rho*height
                  + math.pi*rho*rho*(1 + 1/c) + 2*math.pi*rho*height
                  + math.pi*height + math.pi*math.pi*rho/2 + rim
                  + 4*math.pi/3 + 2*math.pi*(25/6 - math.pi))
        self.assertAlmostEqual(coverage.analytic_volume(row), volume, delta=1e-10)

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
                                            'roll-over-mitres', 'sphere-dome', 'unequal-ordered',
                                            'acute-5deg', 'acute-halfdeg', 'box-all',
                                            'oblique-cylinder'])
            self.assertEqual(coverage.run(args), 0)
            report = json.loads((out / 'scoreboard.json').read_text())
            rows = {r['id']: r for r in report['cases']}
            self.assertEqual(report['stage_counts'], {'checked-export': 9, 'consumer-budget': 1,
                                                     'schema-gap': 1})
            markdown = (out / 'scoreboard.md').read_text()
            self.assertIn('| consumer-budget | 1 |', markdown)
            self.assertIn('| narrow-cone | checked-export |', markdown)
            box = rows['box-single']
            self.assertEqual(box['catallax']['status'], 'built')
            expected = 4000 - 20*(1-math.pi/4)
            for result in (box['raw_step']['result'], box['occt']['result']):
                self.assertAlmostEqual(result['volume'], expected, delta=1e-6)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            self.assertEqual(rows['variable-linear']['catallax']['origin'], 'corpus-not-kernel')
            self.assertEqual(rows['variable-linear']['occt']['status'], 'built')
            self.assertEqual(rows['torus-slice']['occt']['status'], 'built')
            toroidal = rows['torus-slice']['raw_step']['result']

            def disk_integral(radius, height):
                return (height*math.sqrt(max(0, radius*radius-height*height))
                        + radius*radius*math.asin(height/radius))/2

            # Independent washers of the symmetric, rounded torus-cap section.
            volume = 24*math.pi*((4*math.sqrt(2)/5)*(16/45)
                                + disk_integral(1/5, 7/45) - disk_integral(1/5, -1/5)
                                + disk_integral(2, 2) - disk_integral(2, 14/9))
            self.assertAlmostEqual(toroidal['volume'], volume, delta=1e-6)
            self.assertEqual(len(toroidal['witnesses']), 7)
            self.assertTrue(all(w['actual'] == w['expected'] for w in toroidal['witnesses']))
            for name, height, radius in [('acute-5deg', 2, 0.1),
                                         ('acute-halfdeg', 0.2, 0.01)]:
                slope = height/20
                angle = math.atan(slope)
                cotangent = (1+math.sqrt(1+slope*slope))/slope
                volume = 10*(10*height-radius*radius*(cotangent-(math.pi-angle)/2))
                for result in (rows[name]['raw_step']['result'], rows[name]['occt']['result']):
                    self.assertEqual(result['faces'], 6)
                    self.assertAlmostEqual(result['volume'], volume, delta=1e-6)
                    self.assertEqual(len(result['witnesses']), 10)
                    self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            # Integrate washers through the tube and the remaining spherical cap.
            radius, ball, cut = 5, 0.5, 3
            centre = cut+ball
            major = math.sqrt((radius-ball)**2-centre**2)
            contact = radius*centre/(radius-ball)
            t = contact-centre
            tube = ((major*major+ball*ball)*(t+ball)-(t**3+ball**3)/3
                    + 2*major*(disk_integral(ball, t)-disk_integral(ball, -ball)))
            cap = radius*radius*(radius-contact)-(radius**3-contact**3)/3
            volume = math.pi*(tube+cap)
            sphere = rows['sphere-dome']
            self.assertAlmostEqual(sphere['raw_step']['result']['volume'], volume, delta=1e-9)
            self.assertAlmostEqual(sphere['occt']['result']['volume'], volume,
                                   delta=sphere['catallax']['volume_comparison_allowance'])
            self.assertEqual(sphere['raw_step']['result']['faces'], 5)
            self.assertEqual(sphere['catallax']['step']['approximated_curves'], 0)
            for result in (sphere['raw_step']['result'], sphere['occt']['result']):
                self.assertEqual(len(result['witnesses']), 10)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            # Integrate cone disks and the regular outer tube above its tangency.
            t = 1/math.sqrt(5)
            contact = 5+t
            major = (3-math.sqrt(5))/2
            cone = 16*contact-2*contact**2+contact**3/12
            blend = ((major*major+1)*(1-t)-(1-t**3)/3
                     + 2*major*(disk_integral(1, 1)-disk_integral(1, t)))
            volume = math.pi*(cone+blend)
            for result in (rows['narrow-cone']['raw_step']['result'],
                           rows['narrow-cone']['occt']['result']):
                self.assertEqual(result['faces'], 4)
                self.assertAlmostEqual(result['volume'], volume, delta=1e-6)
                self.assertEqual(len(result['witnesses']), 10)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            ordered = rows['unequal-ordered']
            for result in (ordered['raw_step']['result'], ordered['occt']['result']):
                self.assertEqual(result['faces'], 10)
                self.assertEqual(len(result['witnesses']), 10)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            self.assertAlmostEqual(ordered['raw_step']['result']['volume'],
                                   ordered['occt']['result']['volume'], delta=1e-6)
            oblique = rows['oblique-cylinder']
            self.assertEqual(oblique['catallax']['volume_reference'], 'analytic')
            exact = coverage.analytic_volume(
                next(c for c in coverage.cases() if c['id'] == 'oblique-cylinder'))
            self.assertAlmostEqual(oblique['raw_step']['result']['volume'], exact, delta=1e-6)
            self.assertIn('analytic_volume_difference', oblique['occt'])
            for result in (oblique['raw_step']['result'], oblique['occt']['result']):
                self.assertEqual(len(result['witnesses']), 10)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            rounded = rows['roll-over-mitres']['raw_step']['result']
            # OCCT refuses this case, so check the rounded-box section integral.
            self.assertAlmostEqual(rounded['volume'], 3888 + 83*math.pi/3, delta=1e-6)
            self.assertEqual(len(rounded['witnesses']), 5)
            self.assertTrue(all(w['actual'] == w['expected'] for w in rounded['witnesses']))


if __name__ == '__main__':
    unittest.main()
