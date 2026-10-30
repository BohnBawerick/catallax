"""A planar step too narrow to read leaves the STEP file, and only then."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID

from test_loopless_faces import raw_solid, unique
import test_partial_periodic_trims as periodic

FIXTURES = Path(__file__).resolve().parents[2] / 'crates/catallax-cli/tests/fixtures'
# SHA-256 of these exports at e79835b, before this repair existed.
UNCHANGED = {
    ('soliddraft-narrow-just-over', '1'):
        '4a009c7fbbf63d8c7a04c6e070028cff083f35cacc3e0d4d556e462797e829ea',
    ('soliddraft-narrow', '1e-20'):
        '7001cedc9bf0690be5d75bc2a36a0c5608647a4403d492ec5dc9eed0396abd36',
    ('soliddraft-narrow', '1e-30'):
        'efb9409340f591c741de664d5858c6c04bfd0b98ff9a9dc6820be4eaea5b7fcb',
}


def scaled(recipe, factor):
    """Multiplies every length in a captured recipe by an exact factor."""
    k = Fraction(factor)
    times = lambda value: str(Fraction(value) * k)
    recipe = json.loads(json.dumps(recipe))
    for step in recipe['steps']:
        if 'distance' in step:
            step['distance'] = times(step['distance'])
        if 'plane' in step:
            step['plane']['translation'] = [times(x) for x in step['plane']['translation']]
        for loop in step.get('loops', []):
            for segment in loop:
                for key in ('to', 'centre'):
                    if key in segment:
                        segment[key] = [times(x) for x in segment[key]]
                if 'radius' in segment:
                    segment['radius'] = times(segment['radius'])
    return recipe


class NarrowSteps(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='catallax-narrow-step-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def answer(self, recipe):
        root = Path(tempfile.mkdtemp(dir=self.root))
        return periodic.PartialPeriodicTrimImports().answer(recipe, root)

    def recipe(self, name):
        return json.loads((FIXTURES / f'{name}.json').read_text())['problem']

    def assert_reads(self, directory, manifest, faces, file_volume, exact_volume):
        path = directory / 'answer.step'
        reader = STEPControl_Reader()
        self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
        self.assertGreater(reader.TransferRoots(), 0)
        for mode, shape in [('raw', raw_solid(path)), ('standard', reader.OneShape())]:
            with self.subTest(reader=mode):
                self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
                self.assertEqual(len(unique(shape, TopAbs_FACE)), faces)
                volume, spread = periodic.converged_volume(shape)
                self.assertLessEqual(abs(volume - file_volume), 1e-9 + spread)
                area = GProp_GProps()
                BRepGProp.SurfaceProperties_s(shape, area, 1e-12, False)
                deviation = float(manifest['step']['max_deviation'])
                self.assertLessEqual(abs(volume - exact_volume), area.Mass() * deviation + spread + 1e-9)

    def assert_collapse(self, name, width, file_volume, exact_volume):
        directory, manifest = self.answer(self.recipe(name))
        self.assertEqual(manifest['result']['faces'], 11)
        self.assertLessEqual(float(manifest['step']['max_deviation']), 1e-7)
        self.assert_reads(directory, manifest, 9, file_volume, exact_volume)
        report = json.loads((directory / 'step-regularizations.json').read_text())
        self.assertEqual(report['faces'], 9)
        [record] = report['regularizations']
        self.assertEqual((record['removed_faces'], record['removed_edges']), (2, 4))
        # The consumer's reconciliation with the exact body.
        self.assertEqual(report['faces'] + record['removed_faces'], manifest['result']['faces'])
        self.assertAlmostEqual(float(Fraction(record['max_deviation_mm'])), width, delta=1e-16)

    def test_width_120_collapses_a_sliver_of_noise(self):
        self.assert_collapse('soliddraft-narrow', 2.8387853417623135e-15,
                             151869.0266447077, 151869.0266447077)

    def test_the_last_width_below_the_floor_collapses(self):
        # The file loses the whole slab of 900 mm2 by the sliver's width.
        self.assert_collapse('soliddraft-narrow-just-under', 9.999999334925402e-8,
                             151869.02673470764, 151869.02682470766)

    def test_the_first_width_above_the_floor_is_untouched(self):
        directory, manifest = self.answer(self.recipe('soliddraft-narrow-just-over'))
        self.assertEqual(manifest['result']['faces'], 11)
        self.assertFalse((directory / 'step-regularizations.json').exists())
        self.assert_reads(directory, manifest, 11, 151869.0268247077, 151869.0268247077)

    def test_files_this_repair_must_not_touch_keep_their_bytes(self):
        for (name, factor), digest in UNCHANGED.items():
            with self.subTest(name=name, factor=factor):
                directory, _ = self.answer(scaled(self.recipe(name), factor))
                self.assertFalse((directory / 'step-regularizations.json').exists())
                self.assertEqual(hashlib.sha256((directory / 'answer.step').read_bytes()).hexdigest(),
                                 digest)


if __name__ == '__main__':
    unittest.main()
