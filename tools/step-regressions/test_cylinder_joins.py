"""Independent native STEP topology, volume and material for bounded cylinder joins."""
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

import OCP
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/cylinder-joins'


class CylinderJoinImports(unittest.TestCase):
    def test_exact_sectors_and_finite_rims_survive_native_import(self):
        directory = Path(tempfile.mkdtemp(prefix='cylinder-join-step-', dir=os.environ.get('CATALLAX_CYLINDER_JOIN_EVIDENCE')))
        cases = [
            ('boss-plate', 29, 1174 + 1761 * math.pi / 64 + math.pi**2 / 16,
             [(12.15, 8, 8.2, TopAbs_ON), (12.15, 8, 8.19, TopAbs_IN), (12.15, 8, 8.21, TopAbs_OUT),
              (12.15, 8, 3.2, TopAbs_IN), (12.3, 8, 3.3, TopAbs_OUT)]),
            ('stepped-shaft', 8, 5585 * math.pi / 64 + 5 * math.pi**2 / 16,
             [(4.15, 0, 3.2, TopAbs_ON), (4.15, 0, 3.19, TopAbs_IN), (4.15, 0, 3.21, TopAbs_OUT),
              (2.15, 0, 3.2, TopAbs_IN), (2.3, 0, 3.3, TopAbs_OUT)]),
            ('boss-plate-thicken', 4, 85 * math.pi / 16,
             [(12.15, 8, 4, TopAbs_IN), (12.15, 8, 3, TopAbs_ON), (12.15, 8, 2.99, TopAbs_OUT)]),
            ('stepped-shaft-thicken', 4, 85 * math.pi / 16,
             [(2.15, 0, 4, TopAbs_IN), (2.15, 0, 3, TopAbs_ON), (2.15, 0, 2.99, TopAbs_OUT)]),
        ]
        for name, faces, volume, points in cases:
            with self.subTest(name=name):
                out = directory / name
                completed = subprocess.run([str(DRIVER), 'run', str(FIXTURES / (name + '.json')),
                                            '--source', 'recipe', '--binary', str(BINARY), '--out', str(out),
                                            '--fresh', '--workers', '1', '--json'], check=True, capture_output=True, text=True)
                self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                case = json.loads((out / 'run.json').read_text())['cases'][0]
                manifest = case['response']['manifest']
                self.assertTrue(manifest['exact']); self.assertEqual(manifest['marks'], [])
                path = Path(case['path']) / 'answer.step'
                report = CHECK(str(path), faces)
                self.assertEqual(report['status'], 'passed', report)
                self.assertAlmostEqual(report['volume_mm3'], volume, delta=1e-7)
                reader = STEPControl_Reader()
                self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
                self.assertGreater(reader.TransferRoots(), 0)
                for x, y, z, want in points:
                    classifier = BRepClass3d_SolidClassifier(reader.OneShape(), gp_Pnt(x, y, z), 1e-8)
                    self.assertEqual(classifier.State(), want, (x, y, z))
                report['oracle_volume_mm3'] = volume
                report['ocp_version'] = OCP.__version__
                report['point_controls'] = [(x, y, z, str(want)) for x, y, z, want in points]
                (out / 'step-import.json').write_text(json.dumps(report, sort_keys=True, indent=2) + '\n')
        print(f'Cylinder join public STEP evidence: {directory}')


if __name__ == '__main__':
    unittest.main()
