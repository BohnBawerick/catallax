"""Consumer-visible correctness checks for the benchmark report."""
import math
import unittest

from head_to_head import compare_properties, summarize


class CorrectnessTests(unittest.TestCase):
    def test_enclosure_must_cover_independent_value_not_just_have_small_width(self):
        truth = {'volume': 2.0, 'area': 3.0, 'centroid': [0.0, 0.0, 0.0]}
        claimed = {'volume': ['1', '1'], 'area': ['3', '3'],
                   'centroid': [['0', '0'], ['0', '0'], ['0', '0']]}
        result = compare_properties(claimed, truth, True)
        self.assertEqual(result['status'], 'failed')
        self.assertEqual([r['property'] for r in result['wrong_answers']], ['volume'])

    def test_area_centroid_and_nonfinite_values_are_not_hidden_by_correct_volume(self):
        truth = {'volume': 2.0, 'area': 3.0, 'centroid': [0.0, 0.0, 0.0]}
        claimed = {'volume': 2.0, 'area': 4.0, 'centroid': [0.0, math.nan, 0.0]}
        result = compare_properties(claimed, truth, False)
        self.assertEqual(result['status'], 'failed')
        self.assertEqual([(r['property'], r['component']) for r in result['wrong_answers']],
                         [('area', None), ('centroid', 1)])

    def test_unavailable_measurement_is_unmeasured_not_zero_wrong_answers_proof(self):
        self.assertEqual(compare_properties(None, {'volume': 1}, True)['status'], 'not-run')
        report = {'pairs': [{'part': 'P0/example.json', 'warmup': False, 'engines': {
            'catallax': {'status': 'incomplete'}, 'solidcraft': {'status': 'incomplete'}}}]}
        row = summarize(report)[0]
        self.assertIsNone(row['build_median_ms'])
        self.assertEqual(row['incomplete'], 1)
        self.assertEqual(row['properties_checked'], 0)
        self.assertEqual(row['step_passed'], 0)


if __name__ == '__main__':
    unittest.main()
