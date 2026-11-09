#!/usr/bin/env python3
"""Raw STEP read-back of exact finite-face erosion of nonconvex planar solids."""

import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.gp import gp_Pnt
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_IN, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures'
CAST_SOLID = getattr(TopoDS, 'Solid_s', None) or TopoDS.Solid


class TrimmedPlanarOffsetTests(unittest.TestCase):
    def test_raw_erosions_keep_finite_trim_geometry(self):
        # Periodic chart cuts split some native cylindrical patches in the file.
        for name, faces, volume, probes in [
            ('trimmed-planar-offset', 10, 24 - math.pi / 2,
             [((2.1, 2.1, 2), True), ((2.5, 2.5, 2), False), ((1.5, 5, 2), True), ((1.5, 1.5, .5), False)]),
            ('trimmed-planar-boolean-offset', 9, 24 - math.pi / 2,
             [((2.1, 2.1, 2), True), ((2.5, 2.5, 2), False)]),
            ('trimmed-planar-hole-offset', 17, 64 - 2 * math.pi,
             [((1.5, 5, 2), True), ((2.5, 5, 2), False), ((2.2, 2.2, 2), True), ((2.4, 2.4, 2), False)]),
            ('trimmed-planar-corner-offset', 15, 44 - 5 * math.pi / 3,
             [((2.4, 2.4, 2.4), True), ((2.5, 2.5, 2.5), False), ((4, 2.1, 2.1), True), ((4, 2.5, 2.5), False)]),
        ]:
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix='catallax-trimmed-offset-') as temporary:
                out = Path(temporary) / 'run'
                result = subprocess.run([
                    str(DRIVER), 'run', str(FIXTURES / (name + '.json')),
                    '--source', 'recipe', '--binary', str(BINARY), '--out', str(out),
                    '--fresh', '--workers', '1', '--json',
                ], capture_output=True, text=True, check=True)
                self.assertEqual(json.loads(result.stdout)['outcomes'], {'answered': 1})
                case = json.loads((out / 'run.json').read_text())['cases'][0]
                directory = Path(case['path'])
                path = directory / 'answer.step'
                report = CHECK(str(path), faces)
                self.assertEqual(report['status'], 'passed', report)
                self.assertAlmostEqual(report['volume_mm3'], volume, delta=1e-8)
                reader = STEPControl_Reader()
                self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
                self.assertGreater(reader.TransferRoots(), 0)
                explorer = TopExp_Explorer(reader.OneShape(), TopAbs_SOLID)
                self.assertTrue(explorer.More())
                solid = CAST_SOLID(explorer.Current())
                explorer.Next()
                self.assertFalse(explorer.More())
                for point, inside in probes:
                    actual = BRepClass3d_SolidClassifier(solid, gp_Pnt(*point), 1e-9).State()
                    self.assertEqual(actual == TopAbs_IN, inside, point)
                if name == 'trimmed-planar-offset':
                    replay = Path(temporary) / 'replay'
                    again = subprocess.run([str(BINARY), 'replay', str(directory / 'question.json'),
                                            '--out', str(replay)], capture_output=True, text=True, check=True)
                    self.assertEqual(json.loads(again.stdout)['manifest']['status'], 'answered')
                    self.assertEqual(path.read_bytes(), (replay / 'answer.step').read_bytes())


if __name__ == '__main__':
    unittest.main()
