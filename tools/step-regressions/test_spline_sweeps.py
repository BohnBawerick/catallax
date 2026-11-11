#!/usr/bin/env python3
"""Import fixed spline sweeps without healing and check their physical faces."""

import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.TopAbs import TopAbs_FACE
from OCP.TopExp import TopExp_Explorer

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
sys.path.insert(0, str(ROOT / 'tools/checks'))
RAW = runpy.run_path(str(ROOT / 'tools/checks/raw_step.py'))['raw_shapes']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures'


class SplineSweepTests(unittest.TestCase):
    def test_planar_and_spatial_splines_import_as_valid_six_face_solids(self):
        for name in ['planar', 'spatial']:
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix='catallax-spline-sweep-') as temporary:
                out = Path(temporary) / 'run'
                result = subprocess.run([
                    str(DRIVER), 'run', str(FIXTURES / f'spline-sweep-{name}.json'),
                    '--source', 'recipe', '--binary', str(BINARY), '--out', str(out),
                    '--fresh', '--workers', '1', '--json',
                ], capture_output=True, text=True, check=True)
                self.assertEqual(json.loads(result.stdout)['outcomes'], {'answered': 1})
                case = json.loads((out / 'run.json').read_text())['cases'][0]
                path = Path(case['path']) / 'answer.step'
                report = CHECK(str(path), 6)
                self.assertEqual(report['status'], 'passed', report)
                self.assertAlmostEqual(report['volume_mm3'], 24, delta=1e-8)
                shapes = RAW(path)
                self.assertEqual(len(shapes), 1)
                self.assertTrue(BRepCheck_Analyzer(shapes[0]).IsValid())
                explorer = TopExp_Explorer(shapes[0], TopAbs_FACE)
                faces = 0
                while explorer.More():
                    faces += 1
                    explorer.Next()
                self.assertEqual(faces, 6)


if __name__ == '__main__':
    unittest.main()
