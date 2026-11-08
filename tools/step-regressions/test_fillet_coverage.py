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
        self.assertEqual(sum(row.get('recipe') is None for row in rows), 0)
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
        self.assertEqual(coverage.assess(answered, built, None, None)['status'], 'incomplete')
        self.assertEqual(coverage.assess(answered, built, {'status': 'wrong'}, None)['status'], 'wrong')
        different = {'status': 'built', 'result': {'volume': 99}}
        self.assertEqual(coverage.assess(answered, built, different, None)['code'], 'volume.disagreement')
        self.assertEqual(coverage.assess(answered, built, built, None)['status'], 'built')
        refused = {'status': 'refused-with-reason', 'code': 'occt.fillet-not-done'}
        self.assertEqual(coverage.assess(answered, refused, built, None)['status'], 'built')

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

    def test_oblique_cone_volume_is_its_steiner_value(self):
        # The same opening argument on the cut frustum rho = 10 - z/4, with E the
        # cut frustum moved r = 1/2 inward. Its top rim is where the moved plane
        # meets E's cone, rs(t) = top / (c - k s cos t) in the axis angle t. The
        # bottom rim's torus band is added back as a Pappus moment.
        row = next(r for r in coverage.cases() if r['id'] == 'oblique-cone')
        k, zc, s, c, r = 1/4, 8.0, 336/625, 527/625, 0.5
        sk = math.sqrt(1 + k*k)
        rc0 = 10 - r*sk
        rb = rc0 - k*r
        top = c*rc0 - k*(zc*c - r)
        n = 4096
        angles = [2*math.pi*j/n for j in range(n)]
        mean = lambda f: sum(f(t) for t in angles) * 2*math.pi/n
        rs = lambda t: top/(c - k*s*math.cos(t))
        drs = lambda t: -top*k*s*math.sin(t)/(c - k*s*math.cos(t))**2

        def inner(t):
            a = rs(t)
            under_plane = ((zc*c - r)/c - r)*a*a/2 - s*math.cos(t)/c*a**3/3
            under_cone = (rc0/k - r)*(rb*rb - a*a)/2 - (rb**3 - a**3)/(3*k)
            return under_plane + under_cone
        area = (math.pi*rb*rb + mean(lambda t: rs(t)**2/2)/c
                + mean(lambda t: sk/(2*k)*(rb*rb - rs(t)**2)))
        bottom = math.pi/2 + math.atan(k)
        # The rim's half exterior angle times its length element.
        rim = mean(lambda t: math.acos((s*math.cos(t) + c*k)/sk)
                   * math.sqrt(rs(t)**2 + drs(t)**2*(1 + 1/(k*k)))/2)
        curvature = mean(lambda t: (rb - rs(t))/(2*k)) + math.pi*rb*bottom + rim
        opening = mean(inner) + area*r + curvature*r*r + 4*math.pi/3*r**3

        def moment(a, b, d):
            twice = (b[0]-a[0])*(d[1]-a[1]) - (d[0]-a[0])*(b[1]-a[1])
            return abs(twice/2*(a[0]+b[0]+d[0])/3)
        centre, corner = (rb, r), (10.0, 0.0)
        quadrilateral = (moment(centre, (rb, 0.0), corner)
                         + moment(centre, corner, (rb + r/sk, r + r*k/sk)))
        sector = rb*bottom*r*r/2 + r**3/3*(math.sin(math.atan(k)) + 1)
        volume = opening + 2*math.pi*(quadrilateral - sector)
        self.assertAlmostEqual(coverage.analytic_volume(row), volume, delta=1e-9)

    def test_unequal_mitres_volume_is_its_cliff_quadrature(self):
        # Radius-2 uprights take the cliff at the radius-1 top mitres. Each
        # corner loses a wedge of section 4 - pi over the nine units below the
        # rounds. Above that, at angle f on the round, the corner square of side
        # s = 1 + cos f lies outside the upright's disk over an area A(s), until
        # s = sqrt(2). In f the integrand is smooth, so Simpson's rule converges.
        row = next(r for r in coverage.cases() if r['id'] == 'unequal-mitres')

        def section(f):
            s = 1 + math.cos(f)
            area = s*s + math.pi - s*math.sqrt(max(0.0, 4 - s*s)) - 4*math.asin(min(1.0, s/2))
            return area*math.cos(f)

        end, n = math.acos(math.sqrt(2) - 1), 4096
        h = end/n
        top = h/3*sum((1 if k in (0, n) else 4 if k % 2 else 2)*section(k*h) for k in range(n + 1))
        volume = 4000 - 220/3 + 18*math.pi - 4*(9*(4 - math.pi) + top)
        self.assertAlmostEqual(coverage.analytic_volume(row), volume, delta=1e-9)

    def test_narrow_mitres_volume_is_its_rolling_section_integral(self):
        # Radius-1/2 uprights roll over the radius-1 top mitres. Above z = 9 the
        # ball's centre runs on c(t) = (1 - k cos t, 1 - k cos t, 9 + k sin t),
        # k = 1/2, and its contacts c + rho m1, c + rho m2 lie in the normal
        # plane. Each corner loses the integral of |c'| (A - kappa <M, nu>) over
        # the sections between the ball's arc and the rounds' traces.
        row = next(r for r in coverage.cases() if r['id'] == 'narrow-mitres')
        rho, k = 0.5, 0.5
        nodes = [(0.1834346424956498, 0.3626837833783620), (0.5255324099163290, 0.3137066458778873),
                 (0.7966664774136267, 0.2223810344533745), (0.9602898564975363, 0.1012285362903763)]

        def gauss(f, a, b, pieces=16):
            h = (b - a)/pieces
            return h/2*sum(w*(f(a + (i + 0.5)*h - x*h/2) + f(a + (i + 0.5)*h + x*h/2))
                           for i in range(pieces) for x, w in nodes)

        def dot(a, b):
            return sum(x*y for x, y in zip(a, b))

        def section(t):
            c, s = math.cos(t), math.sin(t)
            centre = (1 - k*c, 1 - k*c, 9 + k*s)
            speed = k*math.sqrt(1 + s*s)
            tangent = (s/math.sqrt(1 + s*s), s/math.sqrt(1 + s*s), c/math.sqrt(1 + s*s))
            accel = (k*c, k*c, -k*s)
            normal = [a - dot(accel, tangent)*b for a, b in zip(accel, tangent)]
            kappa = math.sqrt(dot(normal, normal))/speed**2
            nu = [v/math.sqrt(dot(normal, normal)) for v in normal]
            h = (-c/math.sqrt(2*c*c + 4*s*s), -c/math.sqrt(2*c*c + 4*s*s), 2*s/math.sqrt(2*c*c + 4*s*s))
            e = (tangent[1]*h[2] - tangent[2]*h[1], tangent[2]*h[0] - tangent[0]*h[2],
                 tangent[0]*h[1] - tangent[1]*h[0])

            def coords(p):
                q = [a - b for a, b in zip(p, centre)]
                return dot(q, h), dot(q, e)

            def trace(psi):
                # The radius-1 round along x, cut by the normal plane.
                y, z = 1 - math.cos(psi), 9 + math.sin(psi)
                x = centre[0] - ((y - centre[1])*tangent[1] + (z - centre[2])*tangent[2])/tangent[0]
                dp = (-(math.sin(psi)*tangent[1] + math.cos(psi)*tangent[2])/tangent[0],
                      math.sin(psi), math.cos(psi))
                return (x, y, z), dp

            low, high = t, math.pi/2
            for _ in range(200):
                middle = (low + high)/2
                p, _ = trace(middle)
                low, high = (middle, high) if p[0] > p[1] else (low, middle)
            mitre = coords(trace(low)[0])
            contact = coords((centre[0], centre[1] - rho*c, centre[2] + rho*s))

            def green(path, a, b):
                # Area and the first moment along h, by Green's theorem.
                return [gauss(lambda r, i=i: ((lambda u, v, du, dv: ((u*dv - v*du)/2, u*u*dv/2))
                                              (*path(r)))[i], a, b) for i in range(2)]

            def along(psi):
                p, dp = trace(psi)
                return (*coords(p), dot(dp, h), dot(dp, e))

            parts = [green(lambda a: (rho*math.cos(a), rho*math.sin(a), -rho*math.sin(a), rho*math.cos(a)),
                           math.atan2(contact[1], contact[0]), 0),
                     green(lambda r: (rho + r*(mitre[0] - rho), r*mitre[1], mitre[0] - rho, mitre[1]), 0, 1),
                     green(along, low, t)]
            area, moment = (sum(p[i] for p in parts) for i in range(2))
            if area < 0:
                area, moment = -area, -moment
            return speed*(2*area - kappa*2*moment*dot(h, nu))

        corner = gauss(section, 0, math.pi/2)
        volume = 4000 - 220/3 + 18*math.pi - 4*9*rho*rho*(1 - math.pi/4) - 4*corner
        self.assertAlmostEqual(coverage.analytic_volume(row), volume, delta=1e-9)

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
                                      max_work=1_000_000, timeout=300, check=True,
                                      case=['box-single', 'torus-slice', 'variable-linear',
                                            'variable-zero-end', 'narrow-cone',
                                            'roll-over-mitres', 'sphere-dome', 'unequal-ordered',
                                            'acute-5deg', 'acute-halfdeg', 'box-all',
                                            'oblique-cylinder'])
            self.assertEqual(coverage.run(args), 0)
            report = json.loads((out / 'scoreboard.json').read_text())
            rows = {r['id']: r for r in report['cases']}
            box = rows['box-single']
            self.assertEqual(box['catallax']['status'], 'built')
            expected = 4000 - 20*(1-math.pi/4)
            for result in (box['raw_step']['result'], box['occt']['result']):
                self.assertAlmostEqual(result['volume'], expected, delta=1e-6)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            for name, start, end in [('variable-linear', 1, 2), ('variable-zero-end', 0, 2)]:
                variable = rows[name]
                self.assertEqual(variable['catallax']['status'], 'built')
                self.assertEqual(variable['final_fillet']['status'], 'answered')
                self.assertEqual(variable['stage_result'], 'checked-export')
                volume = 4000-20*(start*start+start*end+end*end)/3*(1-math.pi/4)
                for result in (variable['raw_step']['result'], variable['occt']['result']):
                    self.assertAlmostEqual(result['volume'], volume, delta=1e-6)
                    self.assertTrue(result['analytic_volume']['matches'])
                    self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
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

    def test_four_face_vertex_has_the_inner_parallel_volume(self):
        with tempfile.TemporaryDirectory(prefix='fillet-vertex-') as temporary:
            out = Path(temporary) / 'run'
            args = argparse.Namespace(binary=str(BINARY), out=str(out), occt_python=sys.executable,
                                      max_work=20_000_000, timeout=300, check=False,
                                      case=['vertex-4-all'])
            self.assertEqual(coverage.run(args), 0)
            row = json.loads((out / 'scoreboard.json').read_text())['cases'][0]
            # Rounding every edge of this convex pyramid is its inner parallel
            # body plus the ball. The eroded pyramid is homothetic.
            radius = 0.5
            scale = 1-radius*(1+math.sqrt(2))/10
            volume = (4000*scale**3/3
                      + 400*radius*(1+math.sqrt(2))*scale**2
                      + math.pi*radius**2*(30+20*math.sqrt(3)/3)*scale
                      + 4*math.pi*radius**3/3)
            self.assertAlmostEqual(row['raw_step']['result']['volume'], volume, delta=1e-6)
            self.assertEqual(row['catallax']['result']['faces_by_kind'],
                             {'cylinder': 8, 'plane': 5, 'sphere': 5})
            for result in (row['raw_step']['result'], row['occt']['result']):
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))

    def test_one_selected_pyramid_edge_has_the_clipped_section_volume(self):
        with tempfile.TemporaryDirectory(prefix='fillet-partial-vertex-') as temporary:
            out = Path(temporary) / 'run'
            args = argparse.Namespace(binary=str(BINARY), out=str(out), occt_python=sys.executable,
                                      max_work=20_000_000, timeout=300, check=False,
                                      case=['vertex-4-one'])
            self.assertEqual(coverage.run(args), 0)
            row = json.loads((out / 'scoreboard.json').read_text())['cases'][0]
            # Integrate the removed circular segment along the selected edge.
            # Its axial extent is clipped by the base and both apex walls.
            radius = 0.5
            removed = (10*radius**2*(1-math.pi/(2*math.sqrt(3)))
                       - radius**3*(3-(math.pi+2)/math.sqrt(3))/(2*math.sqrt(2)))
            volume = 4000/3-removed
            self.assertEqual(row['catallax']['result']['faces_by_kind'],
                             {'cylinder': 1, 'plane': 5})
            for result in (row['raw_step']['result'], row['occt']['result']):
                self.assertAlmostEqual(result['volume'], volume, delta=1e-6)
                self.assertEqual(result['faces'], 6)
                self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))

    def test_partial_pyramid_corner_has_the_selected_normal_cell_and_cliff(self):
        with tempfile.TemporaryDirectory(prefix='fillet-cliff-vertex-') as temporary:
            out = Path(temporary) / 'run'
            args = argparse.Namespace(binary=str(BINARY), out=str(out), occt_python=sys.executable,
                                      max_work=100_000_000, timeout=300, check=False,
                                      case=['vertex-4-two'])
            self.assertEqual(coverage.run(args), 0)
            row = json.loads((out / 'scoreboard.json').read_text())['cases'][0]
            result = row['raw_step']['result']
            # The exact integral includes the two-strip overlap and the
            # additional selected-normal cone cell outside the common ball.
            self.assertAlmostEqual(result['volume'], 1332.8688454588128, delta=1e-6)
            self.assertEqual(row['catallax']['result']['faces_by_kind'],
                             {'cylinder': 4, 'plane': 6, 'sphere': 1})
            self.assertEqual(row['catallax']['volume_reference'], 'analytic')
            self.assertLessEqual(float(row['catallax']['step']['max_deviation']), 1e-4)
            observed = {tuple(w['point']): w['actual'] for w in result['witnesses']}
            self.assertEqual(observed[('0', '1/100', '99/10')], 'TopAbs_OUT')
            self.assertEqual(observed[('0', '-1/100', '99/10')], 'TopAbs_IN')
            self.assertEqual(observed[('0', '0', '99/10')], 'TopAbs_ON')
            self.assertEqual(observed[('0', '1/10', '97/10')], 'TopAbs_IN')
            if row['occt']['status'] == 'built':
                difference = result['volume']-row['occt']['result']['volume']
                if abs(difference) > row['catallax']['volume_comparison_allowance']:
                    self.assertEqual(row['catallax']['reference_comparison'], 'definition-difference')
                    self.assertTrue(any(w['actual'] != w['expected']
                                        for w in row['occt']['target_witnesses']))

            case_dir = out / 'vertex-4-two'
            bad_case = json.loads((case_dir / 'case.json').read_text())
            bad_case['target_witnesses'][0]['state'] = 'TopAbs_IN'
            bad_path = Path(temporary) / 'wrong-witness.json'
            bad_path.write_text(json.dumps(bad_case))
            record = json.loads((case_dir / 'cli/run.json').read_text())['cases'][0]
            wrong = coverage.occt(bad_path, case_dir / 'recipe.json',
                                  Path(temporary) / 'wrong-witness', args,
                                  Path(record['path']) / 'answer.step')
            self.assertEqual(wrong['status'], 'wrong')
            self.assertEqual(coverage.assess(row['catallax'], row['occt'], wrong,
                                            coverage.analytic_volume(bad_case))['code'],
                             'step.invalid')
            self.assertEqual(coverage.assess(row['catallax'], row['occt'], row['raw_step'],
                                            result['volume']+1)['code'],
                             'volume.disagreement')

    def test_sphere_cylinder_rings_have_the_divergence_volume(self):
        with tempfile.TemporaryDirectory(prefix='fillet-sphere-ring-') as temporary:
            out = Path(temporary) / 'run'
            args = argparse.Namespace(binary=str(BINARY), out=str(out), occt_python=sys.executable,
                                      max_work=20_000_000, timeout=300, check=False,
                                      case=['sphere-cylinder'])
            self.assertEqual(coverage.run(args), 0)
            row = json.loads((out / 'scoreboard.json').read_text())['cases'][0]
            self.assertEqual(row['stage_result'], 'checked-export')
            self.assertEqual(row['catallax']['volume_reference'], 'analytic')
            self.assertEqual(row['catallax']['result']['faces_by_kind'],
                             {'bspline': 2, 'cylinder': 2, 'plane': 2, 'sphere': 1})
            result = row['raw_step']['result']
            self.assertAlmostEqual(result['volume'], 657.67240227301200033, delta=1e-5)
            self.assertEqual(result['faces'], 7)
            self.assertTrue(all(w['actual'] == w['expected'] for w in result['witnesses']))
            # OCCT's own blend misses the independent volume by more than the
            # comparison allows.
            if row['occt']['status'] == 'built':
                self.assertGreater(abs(row['occt']['result']['volume']-657.67240227301200033),
                                   row['catallax']['volume_comparison_allowance'])


if __name__ == '__main__':
    unittest.main()
