#!/usr/bin/env python3
"""Public post-Boolean planar sweeps: independent STEP read-back, without repair."""
import json
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

import OCP  # noqa: F401 -- required, never an optional passing skip
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.gp import gp_Pnt
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_IN, TopAbs_OUT, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/boolean-face-thicken'
CAST_SOLID = getattr(TopoDS, 'Solid_s', None) or TopoDS.Solid


class BooleanFaceThickenStepTests(unittest.TestCase):
    def test_polygonal_sweeps_import_with_exact_oracle_volumes_and_material(self):
        for name, faces, volume, inside, outside in [
            ('collinear-neighbours', 6, 16, (1, 1, 3.5), (1, 1, 2.5)),
            ('l-front', 8, 12, (1, 1, 3.5), (3, 3, 3.5)),
            ('l-back', 8, 12, (1, 1, 2.5), (3, 3, 2.5)),
            ('l-merged-origins', 8, 12, (1, 1, 3.5), (3, 3, 3.5)),
            ('l-tilted', 8, 12, (6, -4.2, 3.9), (8, -3, 5.5)),
            ('hole-front', 10, 32, (1, 1, 3.5), (3, 3, 3.5)),
            ('hole-back', 10, 32, (1, 1, 2.5), (3, 3, 2.5)),
            ('symmetric', 8, 24, (1, 1, 3), (3, 3, 3)),
            ('reversed-wall', 6, 6, (2.5, 3, 1), (1.5, 3, 1)),
            ('split-near', 6, 8, (1, 1, 3.5), (5, 1, 3.5)),
        ]:
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix='catallax-boolean-thicken-step-') as tmp:
                out = Path(tmp) / 'run'
                result = subprocess.run([
                    str(DRIVER), 'run', str(FIXTURES / (name + '.json')),
                    '--source', 'recipe', '--binary', str(BINARY), '--out', str(out),
                    '--fresh', '--workers', '1', '--json',
                ], capture_output=True, text=True, check=True)
                self.assertEqual(json.loads(result.stdout)['outcomes'], {'answered': 1})
                case = json.loads((out / 'run.json').read_text())['cases'][0]
                directory = Path(case['path'])
                path = directory / 'answer.step'
                self.assertTrue(path.is_file())
                self.assertFalse((directory / 'step-regularizations.json').exists())
                report = CHECK(str(path), faces)
                self.assertEqual(report['status'], 'passed', report)
                self.assertAlmostEqual(report['volume_mm3'], volume, delta=1e-9)
                reader = STEPControl_Reader()
                self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
                self.assertGreater(reader.TransferRoots(), 0)
                solids = TopExp_Explorer(reader.OneShape(), TopAbs_SOLID)
                self.assertTrue(solids.More())
                solid = CAST_SOLID(solids.Current())
                solids.Next()
                self.assertFalse(solids.More())
                for point, expected in [(inside, TopAbs_IN), (outside, TopAbs_OUT)]:
                    query = BRepClass3d_SolidClassifier(solid, gp_Pnt(*point), 1e-9)
                    self.assertEqual(query.State(), expected)


if __name__ == '__main__':
    unittest.main()
