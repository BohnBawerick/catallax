"""The seam contract: seamed files read raw, unseamed files refuse."""
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

from test_loopless_faces import BINARY, DRIVER, ROOT, raw_solids

FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures'
IDENTITY = ['1', '0', '0', '0']
TILTED = {'rotation': ['1', '2', '3', '4'], 'translation': ['7/3', '-11/7', '13/5']}


def pose(translation):
    return {'rotation': IDENTITY, 'translation': translation}


def primitive(op, pose_value, **fields):
    return {'format': 'catallax-recipe/3', 'name': op,
            'steps': [{'id': 'body', 'op': op, **fields, 'pose': pose_value}], 'result': 'body'}


def through_hole():
    return {'format': 'catallax-recipe/3', 'name': 'through hole', 'steps': [
        {'id': 'base', 'op': 'block', 'size': ['4', '4', '2'], 'pose': pose(['0', '0', '0'])},
        {'id': 'tool', 'op': 'cylinder', 'radius': '1', 'height': '4', 'pose': pose(['2', '2', '-1'])},
        {'id': 'hole', 'op': 'boolean', 'kind': 'subtract', 'target': 'base', 'tool': 'tool'}],
        'result': 'hole'}


def wrapped(recipe):
    return {'format': 'catallax-step/1', 'recipe': recipe}


def seamed_question(fmt):
    return {'format': fmt, 'id': 'seamed-plate', 'operands': {
        'a': {'shape': 'block', 'size': ['10', '10', '10']},
        'b': {'shape': 'cylinder', 'radius': '3', 'height': '12', 'pose': pose(['5', '5', '-1'])}},
        'steps': [{'id': 's1', 'op': 'subtract', 'args': ['a', 'b']}], 'result': 's1'}


def cut_torus(fmt):
    return {'format': fmt, 'id': 'cut-torus', 'operands': {
        'a': {'shape': 'torus', 'major': '4', 'minor': '1'},
        'b': {'shape': 'block', 'size': ['20', '20', '20'], 'pose': pose(['-10', '-10', '1/2'])}},
        'steps': [{'id': 's1', 'op': 'subtract', 'args': ['a', 'b']}], 'result': 's1'}


class SeamContract(unittest.TestCase):
    def run_problem(self, problem, root, name):
        source = root / f'{name}.json'
        source.write_text(json.dumps(problem))
        out = root / name
        completed = subprocess.run(
            [str(DRIVER), 'run', str(source), '--binary', str(BINARY), '--out', str(out),
             '--fresh', '--workers', '1', '--json'],
            capture_output=True, text=True, check=True)
        self.assertEqual(len(json.loads(completed.stdout)['outcomes']), 1)
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        return case['response']['manifest'], Path(case['path']) / 'answer.step'

    def assert_raw_valid(self, step, volume=None):
        solids = raw_solids(step)
        self.assertTrue(solids)
        for solid in solids:
            self.assertTrue(BRepCheck_Analyzer(solid).IsValid(), step)
        if volume is not None:
            total = 0.0
            for solid in solids:
                properties = GProp_GProps()
                BRepGProp.VolumeProperties_s(solid, properties)
                total += properties.Mass()
            self.assertAlmostEqual(total, volume, delta=volume * 1e-9)

    def test_seamed_parts_read_raw_and_keep_their_bytes(self):
        parts = [
            ('cylinder', primitive('cylinder', TILTED, radius='1', height='2'), 2 * math.pi, True),
            ('frustum', primitive('cone', TILTED, base_radius='3', top_radius='1', height='4'),
             52 * math.pi / 3, True),
            ('pointed-cone', primitive('cone', pose(['0', '0', '0']), base_radius='3',
                                       top_radius='0', height='4'), 12 * math.pi, True),
            ('sphere', primitive('sphere', TILTED, radius='2'), 32 * math.pi / 3, False),
            ('torus', primitive('torus', pose(['0', '0', '0']), major='3', minor='1'),
             6 * math.pi ** 2, True),
            ('through-hole', through_hole(), 32 - 2 * math.pi, True),
            ('spline-sweep', json.loads((FIXTURES / 'spline-sweep-planar.json').read_text()),
             None, False),
            ('loft-circles', json.loads((FIXTURES / 'loft-natural-circles.json').read_text()),
             None, None),
        ]
        with tempfile.TemporaryDirectory(prefix='catallax-seam-contract-') as temporary:
            root = Path(temporary)
            for name, recipe, volume, seam in parts:
                with self.subTest(part=name):
                    legacy, legacy_step = self.run_problem(recipe, root, f'{name}-legacy')
                    manifest, step = self.run_problem(wrapped(recipe), root, f'{name}-seamed')
                    self.assertEqual(manifest['status'], 'answered', manifest.get('refusal'))
                    self.assertEqual(manifest['step']['seams'], 'required')
                    self.assertNotIn('seams', legacy['step'])
                    self.assertEqual(step.read_bytes(), legacy_step.read_bytes())
                    if seam is not None:
                        self.assertEqual('SEAM_CURVE' in step.read_text(), seam)
                    self.assert_raw_valid(step, volume)

    def test_unseamed_files_refuse_under_the_contract(self):
        recipe = json.loads((FIXTURES / 'algebraic-extrusions/circle.json').read_text())
        with tempfile.TemporaryDirectory(prefix='catallax-seam-refusal-') as temporary:
            root = Path(temporary)
            legacy, legacy_step = self.run_problem(recipe, root, 'legacy')
            self.assertEqual(legacy['status'], 'answered')
            # The legacy file leaves its elliptic extrusion without pcurves.
            solids = raw_solids(legacy_step)
            self.assertFalse(solids and all(BRepCheck_Analyzer(s).IsValid() for s in solids))
            manifest, step = self.run_problem(wrapped(recipe), root, 'seamed')
            self.assertEqual(manifest['status'], 'refused')
            self.assertEqual(manifest['refusal']['code'], 'step.seam-unsupported')
            self.assertIn('SURFACE_OF_LINEAR_EXTRUSION', manifest['refusal']['detail'])
            self.assertFalse(step.exists())

    def test_seamed_questions(self):
        with tempfile.TemporaryDirectory(prefix='catallax-seam-question-') as temporary:
            root = Path(temporary)
            legacy, legacy_step = self.run_problem(
                seamed_question('exact-brep-question/2'), root, 'plate-2')
            manifest, step = self.run_problem(
                seamed_question('exact-brep-question/3'), root, 'plate-3')
            self.assertEqual(manifest['status'], 'answered')
            self.assertEqual(manifest['step']['seams'], 'required')
            self.assertEqual(step.read_bytes(), legacy_step.read_bytes())
            self.assert_raw_valid(step, 1000 - 90 * math.pi)
            manifest, step = self.run_problem(cut_torus('exact-brep-question/3'), root, 'torus-3')
            self.assertEqual(manifest['refusal']['code'], 'step.seam-unsupported')
            self.assertFalse(step.exists())


if __name__ == '__main__':
    unittest.main()
