"""Rigidly placed P5 countersinks preserve material and raw STEP validity."""
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_SOLID

from test_loopless_faces import raw_solid, unique
import test_partial_periodic_trims as periodic


class HolePlacementImports(unittest.TestCase):
    def test_p5_countersinks_under_a_general_rigid_placement(self):
        root = Path(__file__).resolve().parents[2]
        recipe = json.loads((root / 'parts/holes/P5-countersinks.json').read_text())
        # q=(1,2,3,4), with squared norm 30. The matrix acts on columns.
        matrix = [[Fraction(-2, 3), Fraction(2, 15), Fraction(11, 15)],
                  [Fraction(2, 3), Fraction(-1, 3), Fraction(2, 3)],
                  [Fraction(1, 3), Fraction(14, 15), Fraction(2, 15)]]

        def transform(coordinates, point):
            offset = [7, 11, -13] if point else [0, 0, 0]
            return [str(sum(a * Fraction(b) for a, b in zip(row, coordinates)) + t)
                    for row, t in zip(matrix, offset)]

        recipe['steps'][0]['plane'] = {
            'rotation': ['1', '2', '3', '4'], 'translation': ['7', '11', '-13']}
        for key in ('point', 'direction'):
            axis = recipe['steps'][1]['axis']
            axis[key] = transform(axis[key], key == 'point')
            for at in recipe['steps'][2]['at']:
                at[key] = transform(at[key], key == 'point')
        with tempfile.TemporaryDirectory(prefix='catallax-hole-placement-') as temporary:
            directory, manifest = periodic.PartialPeriodicTrimImports().answer(recipe, Path(temporary))
            self.assertLessEqual(float(manifest['step']['max_deviation']), 1e-7)
            path = directory / 'answer.step'
            raw = raw_solid(path)
            reader = STEPControl_Reader()
            self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
            self.assertGreater(reader.TransferRoots(), 0)
            for mode, shape in [('raw', raw), ('standard', reader.OneShape())]:
                with self.subTest(reader=mode):
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
                    volume, spread = periodic.converged_volume(shape)
                    area = GProp_GProps()
                    BRepGProp.SurfaceProperties_s(shape, area, 1e-12, False)
                    # Independent revolved annuli minus four drill frusta.
                    self.assertLessEqual(abs(volume - 12495.161180877803),
                                         area.Mass() * 1e-7 + spread)


if __name__ == '__main__':
    unittest.main()
