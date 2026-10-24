"""Independent geometry checks for the public comparison's analytic oracle."""
import importlib.util
import math
from pathlib import Path
import unittest

DATA = Path(__file__).with_name('head-to-head')
spec = importlib.util.spec_from_file_location('benchmark_analytic', DATA/'analytic.py')
analytic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analytic)


class AnalyticTests(unittest.TestCase):
    def test_separated_cylinders_keep_the_small_gap_and_edited_height(self):
        case = {'recipe': 'P0/near-miss.json', 'edit': {'value': '12'}}
        props, _, member = analytic.definition(case, 'edit')
        self.assertGreater(member((5.00000005, 0, 5)), 0)
        self.assertLess(member((0, 0, 11)), 0)
        self.assertGreater(member((10.0000001, 0, 11)), 0)
        self.assertAlmostEqual(props['volume'], 550*math.pi, places=9)

    def test_chamfers_remove_both_outer_rims_but_leave_the_inner_wall(self):
        case = {'recipe': 'P2-collar.json', 'edit': {'value': '0.75'}}
        props, _, member = analytic.definition(case, 'construction')
        self.assertLess(member((11.5, 6, 0)), 0)
        self.assertGreater(member((11.5, 0.25, 0)), 0)
        self.assertGreater(member((5.9, 6, 0)), 0)
        self.assertAlmostEqual(props['volume']/math.pi, 3818/3, places=10)
        self.assertEqual(props['centroid'], [0, 6, 0])

    def test_shell_cavity_and_corner_trim_are_both_required(self):
        case = {'recipe': 'P4-shelled-tray.json', 'edit': {'value': '2.5'}}
        _, _, member = analytic.definition(case, 'construction')
        self.assertLess(member((0, 0, 1)), 0)
        self.assertGreater(member((0, 0, 3)), 0)
        self.assertLess(member((19, 0, 6)), 0)
        self.assertGreater(member((19.9, 14.9, 6)), 0)

    def test_concave_boss_round_adds_material_without_changing_axis_centroid(self):
        small = analytic.properties('P3-filleted-boss.json', 2)
        large = analytic.properties('P3-filleted-boss.json', 2.5)
        self.assertGreater(large['volume'], small['volume'])
        self.assertLess(large['area'], small['area'])
        self.assertEqual(large['centroid'][:2], [20, 15])
        self.assertAlmostEqual(small['volume'], 14428.6401260479, places=7)
        self.assertAlmostEqual(small['area'], 4301.1619509751, places=7)

    def test_rotated_curved_primitives_use_the_inverse_homogeneous_rotation(self):
        cone = {'recipe': 'P0/general-cone.json', 'edit': {'value': '12'}}
        torus = {'recipe': 'P0/general-torus.json', 'edit': {'value': '9'}}
        _, _, in_cones = analytic.definition(cone, 'construction')
        _, _, in_tori = analytic.definition(torus, 'construction')
        self.assertLess(in_cones((2+10/3, 1-10/3, 1+5/3)), 0)
        self.assertLess(in_tori((1+16/3, 2+8/3, 1-16/3)), 0)
        self.assertGreater(in_cones((30, 30, 30)), 0)
        self.assertGreater(in_tori((30, 30, 30)), 0)


if __name__ == '__main__':
    unittest.main()
