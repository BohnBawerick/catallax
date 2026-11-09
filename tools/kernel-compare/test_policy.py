import unittest

from policy import SOLVE_LIMIT_SECONDS, outcome


class CompletionPolicy(unittest.TestCase):
    def setUp(self):
        self.solve = {'status': 'answered', 'wall_ms': 10}
        self.check = {'valid': True, 'complete_sample_check': True,
                      'distance_lower_mm': 0.0009, 'distance_upper_mm': 0.001}

    def test_inclusive_distance_threshold(self):
        self.assertEqual(outcome(self.solve, self.check, 0.001), 'done-sampled')

    def test_proven_excess_distance_is_not_done(self):
        self.assertEqual(outcome(self.solve, self.check, 0.0008), 'not-done')

    def test_straddling_distance_interval_stays_unresolved(self):
        self.assertEqual(outcome(self.solve, self.check, 0.00095), 'unresolved')

    def test_missing_or_incomplete_error_bound_never_passes(self):
        for update in ({'distance_upper_mm': None}, {'complete_sample_check': False},
                       {'distance_upper_mm': float('inf')}):
            with self.subTest(update=update):
                self.assertEqual(outcome(self.solve, {**self.check, **update}, 0.001), 'unresolved')

    def test_invalid_topology_and_failed_witnesses_never_pass(self):
        for update in ({'valid': False}, {'witness_failures': [{'actual': 'outside'}]}):
            with self.subTest(update=update):
                self.assertEqual(outcome(self.solve, {**self.check, **update}, 0.001), 'not-done')

    def test_overtime_answer_and_timeout_are_unsolved(self):
        for solve in ({'status': 'timeout', 'wall_ms': None},
                      {'status': 'answered', 'wall_ms': SOLVE_LIMIT_SECONDS*1000+1}):
            with self.subTest(solve=solve):
                self.assertEqual(outcome(solve, self.check, 0.001), 'unsolved-time-limit')

    def test_exact_time_limit_is_inclusive(self):
        solve = {**self.solve, 'wall_ms': SOLVE_LIMIT_SECONDS*1000}
        self.assertEqual(outcome(solve, self.check, 0.001), 'done-sampled')

    def test_fast_refusal_or_crash_is_not_a_solution(self):
        for status in ('refused', 'crash', 'execution-error', 'protocol-error', 'adapter-unavailable'):
            with self.subTest(status=status):
                self.assertEqual(outcome({'status': status, 'wall_ms': 1}, self.check, 0.001), status)


if __name__ == '__main__':
    unittest.main()
