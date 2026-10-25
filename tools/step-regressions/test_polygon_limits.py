"""Independent STEP imports of finite convex-polygon first-hit caps."""
import copy
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/polygon-limits'


class PolygonLimitImports(unittest.TestCase):
    def test_both_modes_and_branches_have_valid_closed_solids_and_exact_volume_identity(self):
        with tempfile.TemporaryDirectory(prefix='catallax-polygon-limits-') as temporary:
            directory = Path(temporary)
            for fixture in sorted(FIXTURES.glob('*.json')):
                source = json.loads(fixture.read_text())
                points = [[Fraction(x) for x in s['to']] for s in source['steps'][0]['loops'][0]]
                area = abs(sum(a[0] * b[1] - a[1] * b[0]
                               for a, b in zip(points, points[1:] + points[:1]))) / 2
                axis = 2 if source['steps'][1]['op'] == 'sphere' else 0
                for mode in ['extrude_next', 'extrude_to']:
                    volumes = []
                    for branch in ['entry', 'exit']:
                        with self.subTest(fixture=fixture.stem, mode=mode, branch=branch):
                            recipe = copy.deepcopy(source)
                            if branch == 'exit':
                                recipe['steps'][0]['plane']['translation'][axis] = '5'
                            recipe['steps'][2]['op'] = mode
                            direction = ['0', '0', '0']
                            direction[axis] = '2'
                            recipe['steps'][2]['direction'] = direction
                            if mode == 'extrude_to':
                                recipe['steps'][2]['limit'] = {
                                    'face': {'step': 'target', 'role': {'side': {'segment': 'wall'}}}}
                            name = f'{fixture.stem}-{mode}-{branch}'
                            path = directory / (name + '.json')
                            path.write_text(json.dumps(recipe))
                            out = directory / name
                            completed = subprocess.run(
                                [str(DRIVER), 'run', str(path), '--source', 'recipe',
                                 '--binary', str(BINARY), '--out', str(out), '--fresh', '--json'],
                                check=True, capture_output=True, text=True)
                            self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                            case = json.loads((out / 'run.json').read_text())['cases'][0]
                            checked = CHECK(str(Path(case['path']) / 'answer.step'), len(points) + 2)
                            self.assertEqual(checked['status'], 'passed', checked)
                            # The shared import checker reports nonadaptive
                            # diagnostic properties. Request adaptive native
                            # integration here for the independent volume oracle.
                            reader = STEPControl_Reader()
                            self.assertEqual(reader.ReadFile(str(Path(case['path']) / 'answer.step')),
                                             IFSelect_RetDone)
                            self.assertGreater(reader.TransferRoots(), 0)
                            properties = GProp_GProps()
                            error = BRepGProp.VolumePropertiesGK_s(reader.OneShape(), properties,
                                                                 1e-12, True, False, False, False, False)
                            self.assertGreaterEqual(error, 0)
                            self.assertLessEqual(error, 1e-10)
                            volume = properties.Mass()
                            volumes.append(volume)
                            if fixture.stem == 'cylinder-triangle' and branch == 'entry':
                                exact = 2.5 - (11 * math.sqrt(15) / 4 + 4 * math.asin(0.25) - 32 / 3)
                                self.assertAlmostEqual(volume, exact, delta=1e-8)
                    # Independent import properties are a numerical sanity
                    # check; certified exact enclosures are tested in Rust.
                    with self.subTest(fixture=fixture.stem, mode=mode, check='volume_identity'):
                        self.assertAlmostEqual(sum(volumes), float(5 * area), delta=1e-7)


if __name__ == '__main__':
    unittest.main()
