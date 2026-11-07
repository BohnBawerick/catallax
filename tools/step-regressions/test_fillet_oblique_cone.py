"""The exact oblique-cone band through the CLI, read back by OCCT without healing."""
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

BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CASES = ROOT / 'parts/fillets/cases.jsonl'
# The band needs about 34M work; the corpus's 1M allowance stops it in the fillet.
WORK = '40000000'


def fraction(text):
    numerator, _, denominator = text.partition('/')
    return float(numerator) / float(denominator or 1)


class ObliqueConeBand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='catallax-oblique-cone-')
        proof = Path(cls.temporary.name)
        row = next(json.loads(line) for line in CASES.read_text().splitlines()
                   if json.loads(line)['id'] == 'oblique-cone')
        cls.row = row
        recipe = proof / 'oblique-cone.json'
        recipe.write_text(json.dumps(row['recipe'], indent=2) + '\n')
        out = proof / 'run'
        subprocess.run([sys.executable, str(DRIVER), 'run', str(recipe), '--source', 'recipe',
                        '--binary', str(BINARY), '--out', str(out), '--workers', '1',
                        '--max-work', WORK, '--fresh', '--json'],
                       check=True, capture_output=True, text=True)
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        cls.manifest = case['response']['manifest']
        cls.step = Path(case['path']) / 'answer.step'

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_the_answer_is_exact(self):
        self.assertEqual(self.manifest['status'], 'answered', self.manifest)
        self.assertEqual(self.manifest['format'], 'catallax-manifest/5')
        self.assertTrue(self.manifest['exact'])
        self.assertEqual(self.manifest['marks'], [])

    def test_raw_step_holds_the_exact_volume_and_material(self):
        shapes = raw_shapes(self.step)
        self.assertEqual(len(shapes), 1)
        shape = shapes[0]
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        props = GProp_GProps()
        BRepGProp.VolumePropertiesGK_s(shape, props, 1e-12, True, False, False, False, False)
        exact = float(self.row['reference_volume']['volume'])
        self.assertAlmostEqual(props.Mass(), exact, delta=1e-8)
        states = {'TopAbs_IN': TopAbs_IN, 'TopAbs_OUT': TopAbs_OUT}
        for witness in self.row['witnesses']:
            point = gp_Pnt(*map(fraction, witness['point']))
            with self.subTest(point=witness['point']):
                classifier = BRepClass3d_SolidClassifier(shape, point, 1e-7)
                self.assertEqual(classifier.State(), states[witness['state']])


if __name__ == '__main__':
    unittest.main()
