"""Independent oblique-draft polar integrals, unhealed STEP and closed displays."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_IN, TopAbs_OUT
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/checks'))
from raw_step import raw_shapes

BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
FIXTURE = ROOT / 'crates/catallax-features/tests/fixtures/draft-oblique.json'


def closed_volume(path):
    positions, triangles = [], []
    for line in path.read_text().splitlines():
        values = line.split()
        if values and values[0] == 'v':
            positions.append(tuple(map(float, values[1:4])))
        elif values and values[0] == 'f':
            triangles.append(tuple(int(v.split('/')[0]) - 1 for v in values[1:4]))
    uses, orientation = Counter(), Counter()
    volume = 0.0
    for t in triangles:
        for a, b in zip(t, t[1:] + t[:1]):
            key = tuple(sorted((a, b)))
            uses[key] += 1
            orientation[key] += 1 if a < b else -1
        a, b, c = (positions[i] for i in t)
        volume += (a[0] * (b[1]*c[2]-b[2]*c[1]) + a[1] * (b[2]*c[0]-b[0]*c[2])
                   + a[2] * (b[0]*c[1]-b[1]*c[0])) / 6
    if not triangles or set(uses.values()) != {2} or any(orientation.values()):
        raise AssertionError('display is not a closed consistently oriented triangle boundary')
    return volume


class ObliqueDraftOracles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = None
        if os.environ.get('CATALLAX_DRAFT_PROOF'):
            cls.proof = Path(os.environ['CATALLAX_DRAFT_PROOF'])
            cls.proof.mkdir(parents=True, exist_ok=True)
        else:
            cls.temporary = tempfile.TemporaryDirectory(prefix='catallax-oblique-draft-')
            cls.proof = Path(cls.temporary.name)
        cls.base = json.loads(FIXTURE.read_text())

    @classmethod
    def tearDownClass(cls):
        if cls.temporary:
            cls.temporary.cleanup()

    def run_case(self, name, recipe, source):
        input_path = self.proof / (name + '.json')
        input_path.write_text(json.dumps(recipe, indent=2) + '\n')
        directory = self.proof / (name + '-' + source)
        completed = subprocess.run([sys.executable, str(DRIVER), 'run', str(input_path),
                                    '--source', source, '--binary', str(BINARY), '--out', str(directory),
                                    '--workers', '1', '--max-work', '20000000', '--fresh', '--json'],
                                   check=True, capture_output=True, text=True)
        (self.proof / (name + '-' + source + '.driver.json')).write_text(completed.stdout)
        report = json.loads((directory / 'run.json').read_text())
        case = report['cases'][0]
        return Path(case['path']), case['response']['manifest']

    def test_raw_step_matches_independent_polar_integrals_and_material(self):
        # For rho=A(z)+B cos(u)+C sin(u), area is
        # pi*(A(z)^2+(B^2+C^2)/2). These constants are integrated separately
        # from the kernel's carrier, volume measurement and STEP control net.
        cases = []
        for name, angle, pull, neutral_z, factor in [
            ('positive', '45', '1', '1', Fraction(401, 3)),
            ('neutral-start', '45', '1', '0', Fraction(311, 3)),
            ('negative', '-45', '1', '0', Fraction(503, 3)),
            ('reverse-pull', '45', '-1', '0', Fraction(503, 3)),
        ]:
            recipe = deepcopy(self.base)
            draft = recipe['steps'][2]
            draft['angle'], draft['pull'][2], draft['neutral']['point'][2] = angle, pull, neutral_z
            cases.append((name, recipe, factor))
        rotated = deepcopy(self.base)
        rotated['steps'][0]['plane']['rotation'] = ['1', '0', '0', '2']
        rotated['steps'][2]['neutral']['normal'] = ['-11/5', '-2/5', '8']
        cases.append(('rotated', rotated, Fraction(401, 3)))
        scaled = deepcopy(self.base)
        scaled['steps'][0]['loops'][0][0]['radius'] = '4'
        scaled['steps'][1]['distance'] = '1'
        scaled['steps'][2]['neutral']['point'][2] = '1/2'
        cases.append(('half-scale', scaled, Fraction(401, 24)))
        results = []
        for name, recipe, factor in cases:
            with self.subTest(name=name):
                directory, manifest = self.run_case(name, recipe, 'recipe')
                self.assertEqual(manifest['status'], 'answered', manifest)
                self.assertEqual(manifest['format'], 'catallax-manifest/5')
                self.assertTrue(manifest['exact'])
                self.assertEqual(manifest['marks'], [])
                shapes = raw_shapes(directory / 'answer.step')
                self.assertEqual(len(shapes), 1)
                shape = shapes[0]
                self.assertTrue(BRepCheck_Analyzer(shape).IsValid(), name)
                volumes = []
                for epsilon in (1e-9, 1e-11, 1e-13):
                    props = GProp_GProps()
                    error = BRepGProp.VolumePropertiesGK_s(shape, props, epsilon, True, False, False, False, False)
                    self.assertGreaterEqual(error, 0)
                    volumes.append(props.Mass())
                oracle = float(factor) * math.pi
                self.assertLess(abs(volumes[-1]-volumes[-2]), 1e-7)
                self.assertAlmostEqual(volumes[-1], oracle, delta=1e-7)
                z = 0.5 if name == 'half-scale' else 1.0
                for xyz, expected in [((0, 0, z), TopAbs_IN), ((20, 0, z), TopAbs_OUT),
                                      ((0, 0, -1), TopAbs_OUT), ((0, 0, 3), TopAbs_OUT)]:
                    self.assertEqual(BRepClass3d_SolidClassifier(shape, gp_Pnt(*xyz), 1e-7).State(), expected)
                field_directory, field = self.run_case(name, recipe, 'field')
                self.assertEqual(field['status'], 'answered', field)
                display = field['face_measurements']['display']
                self.assertLessEqual(Fraction(display['deviation_bound_mm']), Fraction(1, 100))
                field_volume = closed_volume(field_directory / 'display.obj')
                self.assertAlmostEqual(field_volume, oracle, delta=0.5)
                results.append({'case': name, 'raw_valid': True, 'volume_mm3': volumes[-1],
                                'oracle_mm3': oracle, 'integration_spread_mm3': max(volumes)-min(volumes),
                                'field_display_volume_mm3': field_volume})
        (self.proof / 'raw-oracles.json').write_text(json.dumps(results, indent=2) + '\n')


    def test_brep_print_display_closes_without_a_field_substitute(self):
        request = {'format': 'catallax-print/1', 'recipe': self.base, 'route': 'brep', 'tolerance': '1/10'}
        directory, manifest = self.run_case('brep-display', request, 'print')
        self.assertEqual(manifest['status'], 'answered', manifest)
        self.assertEqual(manifest['printing']['route'], 'brep')
        volume = closed_volume(directory / 'display.obj')
        self.assertAlmostEqual(volume, 401 * math.pi / 3, delta=0.5)


if __name__ == '__main__':
    unittest.main()
