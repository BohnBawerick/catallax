"""Mixed circular/polygonal openings against raw STEP and section formulas."""
from copy import deepcopy
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
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures'


class MixedOpeningOracles(unittest.TestCase):
    def test_draft_and_shell_keep_both_openings_and_round_polygon_corners(self):
        with tempfile.TemporaryDirectory(prefix='catallax-mixed-openings-') as temporary:
            proof = Path(os.environ.get('CATALLAX_MIXED_PROOF', temporary))
            proof.mkdir(parents=True, exist_ok=True)
            rows = []
            source_volume = 12 * (40 * 30 - 72 - math.pi * 3**2)
            controls = []
            for angle in ['3', '5', '-3']:
                recipe = json.loads((FIXTURES / 'draft-mixed-openings.json').read_text())
                recipe['steps'][-1]['angle'] = angle
                expected = source_volume - 30 * 12**2 / 2 * math.tan(math.radians(float(angle)))
                controls.append(('draft-' + angle, recipe, expected, [
                    ((10, 15, 6), TopAbs_OUT), ((28, 15, 6), TopAbs_OUT),
                    ((20, 15, 6), TopAbs_IN), ((1, 15, 6), TopAbs_IN),
                ]))
            for family, filename, area, perimeter in [
                ('six-sided', 'shell-mixed-openings.json', 72, 32),
                ('regular', 'shell-regular-openings.json', 54*math.sqrt(3), 36),
            ]:
                for value in ['1', '5/4']:
                    t = 1 if value == '1' else 1.25
                    recipe = json.loads((FIXTURES / filename).read_text())
                    recipe['steps'][-1]['thickness'] = value
                    # Steiner's planar formula adds perimeter*t and circular
                    # corner sectors whose angles sum to a complete circle.
                    cavity_area = (40-2*t)*(30-2*t) - (area+perimeter*t+math.pi*t*t) - math.pi*(3+t)**2
                    expected = 12*(1200-area-9*math.pi) - (12-t)*cavity_area
                    controls.append(('shell-' + family + '-' + value.replace('/', '_'), recipe, expected, [
                        ((10,15,6), TopAbs_OUT), ((28,15,6), TopAbs_OUT),
                        ((20,15,6), TopAbs_OUT), ((13+0.9*t,15,6), TopAbs_IN),
                        ((13+1.1*t,15,6), TopAbs_OUT), ((22-0.9*t,15,6), TopAbs_IN),
                        ((22-1.1*t,15,6), TopAbs_OUT), ((20,15,0.5*t), TopAbs_IN),
                    ]))
            for name, recipe, expected, points in controls:
                with self.subTest(name=name):
                    path = proof / (name + '.json')
                    path.write_text(json.dumps(deepcopy(recipe), indent=2) + '\n')
                    output = proof / name
                    result = subprocess.run([sys.executable, str(DRIVER), 'run', str(path),
                        '--source', 'recipe', '--binary', str(BINARY), '--out', str(output),
                        '--workers', '1', '--max-work', '20000000', '--fresh', '--json'],
                        check=True, capture_output=True, text=True)
                    (proof / (name + '.driver.json')).write_text(result.stdout)
                    case = json.loads((output / 'run.json').read_text())['cases'][0]
                    response = case['response']
                    self.assertEqual(response['manifest']['status'], 'answered', response)
                    self.assertTrue(response['manifest']['exact'])
                    self.assertLessEqual(response['work']['used'], 20_000_000)
                    shapes = raw_shapes(Path(case['path']) / 'answer.step')
                    self.assertEqual(len(shapes), 1)
                    shape = shapes[0]
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    props = GProp_GProps()
                    error = BRepGProp.VolumePropertiesGK_s(shape, props, 1e-12, True, False, False, False, False)
                    self.assertGreaterEqual(error, 0)
                    self.assertAlmostEqual(props.Mass(), expected, delta=1e-6)
                    for point, state in points:
                        self.assertEqual(BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7).State(), state, point)
                    rows.append({'case': name, 'raw_valid': True, 'volume_mm3': props.Mass(),
                                 'oracle_mm3': expected, 'work': response['work'],
                                 'operation_log': response['operation_log']})
            (proof / 'oracles.json').write_text(json.dumps(rows, indent=2) + '\n')


if __name__ == '__main__':
    unittest.main()
