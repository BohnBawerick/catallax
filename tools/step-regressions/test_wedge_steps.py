"""A planar wedge too narrow to read leaves the STEP file, and only then."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID

from test_loopless_faces import raw_solid, unique
import test_partial_periodic_trims as periodic

FIXTURES = Path(__file__).resolve().parents[2] / 'crates/catallax-cli/tests/fixtures'
# SHA-256 of these raw-valid exports at e7f307b, before this repair existed:
# the same brackets drawn with exact frames, and PR 212's two collapsed ribbons.
UNCHANGED = {
    'soliddraft-flush-bracket': '3e6e36dd5e86120af59038cf41d7cc587404195cb04ffacbd133d4a315225207',
    'soliddraft-flush-minimal': 'a06109a6e3aed30bb7c469d8b51278e445c833c24b31b3cdbe29f469f0691656',
    'soliddraft-narrow': 'b052c2ab1a8b8825081880b9126c21ebe7ae124d4cc2eb7253ed7b51c2dc653e',
    'soliddraft-narrow-just-under':
        'e2597e2cc18ec5c3a63ec8af0f9139b62a3c59a8f0df3346e39b13acc4ce428e',
}
# The plate's wall reconnection, then the wedge in front of the flange, then
# the 2e-28 mm strip beside it that closing the wedge leaves at its far end.
REPAIRS = [(1, 0), (1, 2), (1, 2)]


class WedgeSteps(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='catallax-wedge-step-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def answer(self, recipe):
        root = Path(tempfile.mkdtemp(dir=self.root))
        return periodic.PartialPeriodicTrimImports().answer(recipe, root)

    def recipe(self, name):
        return json.loads((FIXTURES / f'{name}.json').read_text())['problem']

    def assert_reads(self, directory, faces, volume):
        path = directory / 'answer.step'
        reader = STEPControl_Reader()
        self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
        self.assertGreater(reader.TransferRoots(), 0)
        for mode, shape in [('raw', raw_solid(path)), ('standard', reader.OneShape())]:
            with self.subTest(reader=mode):
                self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
                self.assertEqual(len(unique(shape, TopAbs_FACE)), faces)
                found, spread = periodic.converged_volume(shape)
                # The flush bracket's volume to the reported six decimals.
                self.assertLessEqual(abs(found - volume), 5e-7 + spread)

    def assert_collapse(self, recipe, model_faces, faces, volume):
        directory, manifest = self.answer(recipe)
        self.assertEqual(manifest['result']['faces'], model_faces)
        self.assertLessEqual(float(manifest['step']['max_deviation']), 1e-7)
        self.assert_reads(directory, faces, volume)
        report = json.loads((directory / 'step-regularizations.json').read_text())
        records = report['regularizations']
        self.assertEqual([(r['removed_faces'], r['removed_edges']) for r in records], REPAIRS)
        # The consumer's reconciliation with the exact body.
        self.assertEqual(report['faces'], faces)
        self.assertEqual(report['faces'] + sum(r['removed_faces'] for r in records), model_faces)
        self.assertEqual(report['edges'] + sum(r['removed_edges'] for r in records),
                         manifest['result']['edges'] + manifest['step']['added_edges'])

    def test_the_gui_bracket_reads_as_the_flush_bracket(self):
        self.assert_collapse(self.recipe('soliddraft-gui-bracket'), 18, 15, 170521.057873)

    def test_the_gui_bracket_without_its_bore_reads_as_the_flush_one(self):
        self.assert_collapse(self.recipe('soliddraft-gui-minimal'), 17, 14, 174366.367281)

    def test_the_first_union_alone_is_the_cause(self):
        recipe = copy.deepcopy(self.recipe('soliddraft-gui-minimal'))
        index = next(i for i, step in enumerate(recipe['steps']) if step['id'] == 'op_61')
        recipe['steps'] = recipe['steps'][:index + 1]
        recipe['result'] = 'op_61'
        self.assert_collapse(recipe, 11, 8, 168000.0)

    def test_raw_valid_files_keep_their_bytes(self):
        for name, digest in UNCHANGED.items():
            with self.subTest(name=name):
                directory, _ = self.answer(self.recipe(name))
                self.assertEqual(hashlib.sha256((directory / 'answer.step').read_bytes()).hexdigest(),
                                 digest)


if __name__ == '__main__':
    unittest.main()
