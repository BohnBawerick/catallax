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

    def test_all_four_countersinks_are_empty_and_junction_round_is_filled(self):
        case = {'recipe': 'P5-flanged-housing.json', 'edit': {'value': '2.5'}}
        _, _, member = analytic.definition(case, 'construction')
        for x, z in ((17, 0), (0, -17), (-17, 0), (0, 17)):
            self.assertGreater(member((x, 2, z)), 0)
        self.assertLess(member((13, 5.1, 0)), 0)
        self.assertGreater(member((13, 7, 0)), 0)
        self.assertGreater(member((0, 10, 0)), 0)


if __name__ == '__main__':
    unittest.main()
