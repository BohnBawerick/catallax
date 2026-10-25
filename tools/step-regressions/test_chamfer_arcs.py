#!/usr/bin/env python3
"""Public finite chamfers: native STEP topology, material controls and volume.

Set CATALLAX_CHAMFER_EVIDENCE to retain requests, manifests, raw STEP and import
reports. Native import and numeric volume checks are not deviation certificates.
"""
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

import OCP  # Required consumer, never a silently skipped qualification.
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/chamfer-arcs'


class ChamferArcStepTests(unittest.TestCase):
    def test_finite_cones_survive_native_import_edit_and_reflection(self):
        retained = os.environ.get('CATALLAX_CHAMFER_EVIDENCE')
        if retained:
            Path(retained).mkdir(parents=True, exist_ok=True)
        directory = Path(tempfile.mkdtemp(prefix='chamfer-step-', dir=retained))
        # Keep evidence on failures even when no persistent destination was set.
        inputs = directory / 'inputs'
        inputs.mkdir()
        expected = {}
        for name, faces in [('convex', 6), ('concave', 9), ('chord', 6), ('seam', 5), ('major', 6), ('major-seam', 6)]:
            for distance, other, mirrored in [('2', '1', False), ('2', '2', False),
                                               ('1', '1/2', False), ('2', '1', True)]:
                value = json.loads((FIXTURES / f'{name}.json').read_text())
                bevel = value['steps'][-1]
                bevel.update(distance=distance, other_distance=other)
                if mirrored:
                    value['steps'].append({'id': 'reflected', 'op': 'mirror', 'body': 'bevel',
                                          'plane': {'point': ['0', '0', '0'], 'normal': ['3', '0', '0']},
                                          'keep_original': False})
                    value['result'] = 'reflected'
                key = f'{name}-{distance}-{other.replace("/", "_")}-{mirrored}'
                value['name'] = key
                (inputs / f'{key}.json').write_text(json.dumps(value, indent=2) + '\n')
                expected[key] = (name, faces, Fraction(distance), Fraction(other), mirrored)
        run = directory / 'run'
        result = subprocess.run([str(DRIVER), 'run', str(inputs), '--source', 'recipe',
                                 '--binary', str(BINARY), '--max-work', '1000000',
                                 '--workers', '2', '--out', str(run), '--fresh', '--json'],
                                capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(result.stdout)['outcomes'], {'answered': len(expected)}, directory)
        for case in json.loads((run / 'run.json').read_text())['cases']:
            name, faces, distance, other, mirrored = expected[case['id']]
            with self.subTest(case=case['id']):
                path = Path(case['path'])
                manifest = case['response']['manifest']
                self.assertTrue(manifest['exact'], manifest)
                self.assertEqual(manifest['marks'], [])
                self.assertEqual(manifest['result']['faces'], faces)
                report = CHECK(str(path / 'answer.step'), faces)
                self.assertEqual(report['status'], 'passed', report)
                self.assertEqual(report['version'], OCP.__version__)
                if name in ('convex', 'concave', 'major', 'major-seam'):
                    direction = 1 if name == 'concave' else -1
                    base = 224 if name == 'concave' else 192
                    factor = Fraction(base) + direction * other * (8 * distance + direction * distance**2 / 3) / 4
                    if name.startswith('major'):
                        factor *= 3
                    volume = float(factor) * math.pi
                    self.assertAlmostEqual(report['volume_mm3'], volume, delta=1e-6)
                    report['analytic_volume_mm3'] = volume
                    report['analytic_volume_check'] = 'passed-numeric-readback'
                reader = STEPControl_Reader()
                self.assertEqual(reader.ReadFile(str(path / 'answer.step')), IFSelect_RetDone)
                self.assertGreater(reader.TransferRoots(), 0)
                shape = reader.OneShape()
                z = Fraction(4) + other / 2 if name == 'concave' else Fraction(12) - other / 2
                radius = Fraction(8) + distance / 2 if name == 'concave' else Fraction(8) - distance / 2
                # y=0 is a seam/control generator for all these authored domains.
                # OCCT checks it independently even if an exact B-rep ray consumer refuses it.
                checks = [(radius, Fraction(0), z, TopAbs_ON)]
                if name == 'seam':
                    checks += [(radius, Fraction(-1), z, TopAbs_OUT),
                               (radius - Fraction(1, 2), Fraction(-1), z, TopAbs_IN)]
                report['point_controls'] = []
                for x, y, height, state in checks:
                    x = -x if mirrored else x
                    classifier = BRepClass3d_SolidClassifier(shape, gp_Pnt(float(x), float(y), float(height)), 1e-7)
                    self.assertEqual(classifier.State(), state, (name, x, y, height))
                    report['point_controls'].append({'point': list(map(str, (x, y, height))), 'state': str(state)})
                (path / 'step-import.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        print(f'Retained public chamfer STEP evidence: {directory}')


if __name__ == '__main__':
    unittest.main()
