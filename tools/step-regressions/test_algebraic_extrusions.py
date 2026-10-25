"""Independent native STEP validity and analytic volume for signed oblique sweeps."""
import copy
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/algebraic-extrusions'


class AlgebraicExtrusionImports(unittest.TestCase):
    def test_signed_line_circle_and_arc_profiles(self):
        with tempfile.TemporaryDirectory(prefix='catallax-algebraic-extrusions-') as temporary:
            directory = Path(temporary)
            for name, area, faces in [('rectangle', 8, 6), ('circle', 4 * math.pi, 3),
                                      ('slot', 8 + math.pi, 6)]:
                source = json.loads((FIXTURES / (name + '.json')).read_text())
                for distance in [5, -5]:
                    with self.subTest(profile=name, distance=distance):
                        recipe = copy.deepcopy(source)
                        recipe['steps'][1]['distance'] = str(distance)
                        path = directory / f'{name}-{distance}.json'
                        path.write_text(json.dumps(recipe))
                        out = directory / f'{name}-{distance}'
                        completed = subprocess.run([str(DRIVER), 'run', str(path), '--source', 'recipe',
                                                    '--binary', str(BINARY), '--out', str(out), '--fresh',
                                                    '--workers', '1', '--json'], check=True,
                                                   capture_output=True, text=True)
                        self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                        case = json.loads((out / 'run.json').read_text())['cases'][0]
                        manifest = case['response']['manifest']
                        self.assertTrue(manifest['exact'])
                        self.assertEqual(manifest['marks'], [])
                        checked = CHECK(str(Path(case['path']) / 'answer.step'), faces)
                        self.assertEqual(checked['status'], 'passed', checked)
                        # Import properties are a separate analytic sanity check;
                        # exact enclosure assertions live in the Rust tests.
                        self.assertAlmostEqual(checked['volume_mm3'], area * abs(distance) / math.sqrt(2), delta=1e-8)


if __name__ == '__main__':
    unittest.main()
