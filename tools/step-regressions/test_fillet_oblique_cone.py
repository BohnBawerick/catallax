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


def mirrored(recipe):
    """The recipe mirrored in x = 0 after the round, keeping only the copy."""
    steps = recipe['steps'] + [{
        'id': 'reflected', 'op': 'mirror', 'body': recipe['result'],
        'plane': {'point': ['0', '0', '0'], 'normal': ['1', '0', '0']},
        'keep_original': False}]
    return {**recipe, 'steps': steps, 'result': 'reflected'}


def negated(text):
    return text[1:] if text.startswith('-') else ('0' if text == '0' else '-' + text)


def build(recipe, proof, name):
    path = proof / f'{name}.json'
    path.write_text(json.dumps(recipe, indent=2) + '\n')
    out = proof / name
    subprocess.run([sys.executable, str(DRIVER), 'run', str(path), '--source', 'recipe',
                    '--binary', str(BINARY), '--out', str(out), '--workers', '1',
                    '--max-work', WORK, '--fresh', '--json'],
                   check=True, capture_output=True, text=True)
    case = json.loads((out / 'run.json').read_text())['cases'][0]
    return case['response']['manifest'], Path(case['path']) / 'answer.step'


class ObliqueConeBand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='catallax-oblique-cone-')
        proof = Path(cls.temporary.name)
        row = next(json.loads(line) for line in CASES.read_text().splitlines()
                   if json.loads(line)['id'] == 'oblique-cone')
        cls.row = row
        cls.manifest, cls.step = build(row['recipe'], proof, 'run')
        cls.mirrored_manifest, cls.mirrored_step = build(mirrored(row['recipe']), proof, 'mirror')

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_the_answer_is_exact(self):
        for manifest in (self.manifest, self.mirrored_manifest):
            self.assertEqual(manifest['status'], 'answered', manifest)
            self.assertEqual(manifest['format'], 'catallax-manifest/5')
            self.assertTrue(manifest['exact'])
            self.assertEqual(manifest['marks'], [])

    def test_raw_step_holds_the_exact_volume_and_material(self):
        self.check_raw(self.step, self.row['witnesses'])

    # A mirror after the round reverses the band's frame; its rails must keep
    # their sense for the raw STEP to stay valid.
    def test_a_mirrored_band_keeps_its_volume_and_material(self):
        witnesses = [{**w, 'point': [negated(w['point'][0])] + w['point'][1:]}
                     for w in self.row['witnesses']]
        self.check_raw(self.mirrored_step, witnesses)

    def check_raw(self, step, witnesses):
        shapes = raw_shapes(step)
        self.assertEqual(len(shapes), 1)
        shape = shapes[0]
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        props = GProp_GProps()
        BRepGProp.VolumePropertiesGK_s(shape, props, 1e-12, True, False, False, False, False)
        exact = float(self.row['reference_volume']['volume'])
        self.assertAlmostEqual(props.Mass(), exact, delta=1e-8)
        states = {'TopAbs_IN': TopAbs_IN, 'TopAbs_OUT': TopAbs_OUT}
        for witness in witnesses:
            point = gp_Pnt(*map(fraction, witness['point']))
            with self.subTest(point=witness['point']):
                classifier = BRepClass3d_SolidClassifier(shape, point, 1e-7)
                self.assertEqual(classifier.State(), states[witness['state']])


if __name__ == '__main__':
    unittest.main()
