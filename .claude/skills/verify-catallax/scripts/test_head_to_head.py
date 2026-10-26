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

    def test_fast_refusal_is_not_a_successful_display_latency(self):
        engine = {'status': 'completed', 'variants': {'construction': {
            'build': {'status': 'answered', 'wall_ms': 1},
            'display': {'status': 'work-exhausted', 'wall_ms': 2, 'refusal_ms': 2}}}}
        report = {'pairs': [{'part': 'P0/example.json', 'warmup': False,
                            'engines': {'catallax': engine, 'solidcraft': engine}}]}
        row = summarize(report)[0]
        self.assertIsNone(row['display_median_ms'])
        self.assertIsNone(row['build_and_display_median_ms'])
        self.assertEqual(row['phase_refusals']['display'][0]['refusal_ms'], 2)

    def test_new_tool_refusal_does_not_become_zero_latency_or_checked_success(self):
        report = {'engines': ['cadquery', 'pali_ijo'], 'pairs': [{
            'part': 'P2-collar.json', 'warmup': False, 'engines': {
                'cadquery': {'status': 'completed', 'variants': {'construction': {
                    'build': {'status': 'answered', 'wall_ms': 2},
                    'correctness': {'membership': {'status': 'failed', 'checked': 200,
                                                  'boundary_unresolved': 16, 'wrong': [{}]}}}}},
                'pali_ijo': {'status': 'completed', 'variants': {'construction': {
                    'build': {'status': 'not-expressible'}}}}}}]}
        rows = summarize(report)
        cadquery, pali = rows[0], rows[2]
        self.assertEqual(cadquery['membership_points_checked'], 200)
        self.assertEqual(cadquery['wrong_membership'], 1)
        self.assertIsNone(pali['build_median_ms'])
        self.assertEqual(pali['not_expressible'], 1)
        self.assertEqual(pali['failed_rebuilds'], 0)
        self.assertEqual(pali['properties_checked'], 0)

    def test_one_missing_sample_suppresses_successful_median(self):
        report = {'engines': ['smlibpy'], 'pairs': [
            {'part': 'P1-bracket.json', 'warmup': False, 'engines': {'smlibpy': {
                'status': 'completed', 'variants': {'construction': {'build': build}}}}}
            for build in ({'status': 'answered', 'wall_ms': 3}, {'status': 'command-failed'})]}
        row = summarize(report)[0]
        self.assertIsNone(row['build_median_ms'])
        self.assertEqual(row['failed_rebuilds'], 1)

    def test_preview_wait_is_separate_from_proof_completion(self):
        def phase(ms, status='answered'):
            return {'status': status, 'wall_ms': ms}
        variant = {'build': phase(300), 'preview': phase(200), 'display': phase(900),
                   'measurement': phase(400), 'step': phase(500)}
        engine = {'status': 'completed', 'variants': {'construction': variant}}
        row = summarize({'pairs': [{'part': 'P0/example.json', 'warmup': False,
                                    'engines': {'catallax': engine, 'solidcraft': engine}}]})[0]
        self.assertEqual(row['preview_wait_median_ms'], 500)
        self.assertEqual(row['preview_budget_passes'], 1)
        self.assertEqual(row['proof_wait_median_ms'], 2100)
        self.assertFalse(row['built_after_proof'], 'proof over 2000 ms is not built')

    def test_missing_or_refused_preview_is_no_wait_and_preview_never_builds(self):
        quick = {'build': {'status': 'answered', 'wall_ms': 1}, 'preview': {'status': 'work-exhausted', 'wall_ms': 1},
                 'display': {'status': 'work-exhausted', 'wall_ms': 1}}
        legacy = {'build': {'status': 'answered', 'wall_ms': 1}}
        for variant in (quick, legacy):
            engine = {'status': 'completed', 'variants': {'construction': variant}}
            row = summarize({'pairs': [{'part': 'P0/example.json', 'warmup': False,
                                        'engines': {'catallax': engine, 'solidcraft': engine}}]})[0]
            self.assertIsNone(row['preview_wait_median_ms'])
            self.assertEqual(row['preview_budget_passes'], 0)
            self.assertFalse(row['built_after_proof'])


if __name__ == '__main__':
    unittest.main()
