"""Three-edge bevel corners against unhealed STEP and independent volumes."""
from copy import deepcopy
from fractions import Fraction
import json
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
from test_draft_oblique import closed_volume

BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
FIXTURE = ROOT / 'crates/catallax-features/tests/fixtures/chamfer-three-edge.json'


class ChamferCornerOracles(unittest.TestCase):
    def test_equal_unequal_and_edited_corners_keep_the_material_and_closed_boundary(self):
        with tempfile.TemporaryDirectory(prefix='catallax-chamfer-corner-') as temporary:
            proof = Path(os.environ.get('CATALLAX_CHAMFER_PROOF', temporary))
            proof.mkdir(parents=True, exist_ok=True)
            base = json.loads(FIXTURE.read_text())
            rows = []
            for name, first, second, expected in [
                ('equal', '1', '1', Fraction(2837, 3)),
                ('unequal', '1', '2', Fraction(932)),
                ('edited', '5/4', '3/4', Fraction(15145, 16)),
            ]:
                with self.subTest(name=name):
                    recipe = deepcopy(base)
                    recipe['steps'][1].update(distance=first, other_distance=second)
                    results = {}
                    for source, request in [('recipe', recipe), ('print', {
                        'format': 'catallax-print/1', 'recipe': recipe,
                        'route': 'brep', 'tolerance': '1/100',
                    })]:
                        path = proof / (name + '-' + source + '.json')
                        path.write_text(json.dumps(request, indent=2) + '\n')
                        output = proof / (name + '-' + source)
                        result = subprocess.run([sys.executable, str(DRIVER), 'run', str(path),
                            '--source', source, '--binary', str(BINARY), '--out', str(output),
                            '--workers', '1', '--max-work', '20000000', '--fresh', '--json'],
                            check=True, capture_output=True, text=True)
                        (proof / (name + '-' + source + '.driver.json')).write_text(result.stdout)
                        case = json.loads((output / 'run.json').read_text())['cases'][0]
                        response = case['response']
                        self.assertEqual(response['manifest']['status'], 'answered', response)
                        results[source] = (Path(case['path']), response)
                    directory, response = results['recipe']
                    self.assertTrue(response['manifest']['exact'])
                    self.assertLessEqual(response['work']['used'], 20_000_000)
                    shapes = raw_shapes(directory / 'answer.step')
                    self.assertEqual(len(shapes), 1)
                    shape = shapes[0]
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    props = GProp_GProps()
                    error = BRepGProp.VolumePropertiesGK_s(shape, props, 1e-12, True, False, False, False, False)
                    self.assertGreaterEqual(error, 0)
                    self.assertAlmostEqual(props.Mass(), float(expected), delta=1e-7)
                    for point, state in [((0.25, 0.25, 0.25), TopAbs_OUT), ((2, 2, 2), TopAbs_IN)]:
                        self.assertEqual(BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7).State(), state)
                    display_volume = closed_volume(results['print'][0] / 'display.obj')
                    self.assertAlmostEqual(display_volume, float(expected), delta=1e-7)
                    rows.append({'case': name, 'raw_valid': True, 'volume_mm3': props.Mass(),
                                 'oracle_mm3': float(expected), 'display_volume_mm3': display_volume,
                                 'work': response['work'], 'operation_log': response['operation_log']})
            (proof / 'oracles.json').write_text(json.dumps(rows, indent=2) + '\n')


if __name__ == '__main__':
    unittest.main()
