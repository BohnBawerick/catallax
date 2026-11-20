#!/usr/bin/env python3
"""The three M10 loft top cuts must answer below the 20M target.

These traces take substantially longer than the regular cases. The required
`step-regressions-slow` job runs
this file beside `step-regressions`; the default `test*.py` discovery pattern
does not collect it. Run it with:

    python3 -B -m unittest discover -s tools/step-regressions -p 'slow_*.py' -v
"""

import unittest

from test_m10_spline_boolean import TOP_CASES, M10Outcomes, cases


class M10TopCutAnswers(M10Outcomes, unittest.TestCase):
    def test_top_requests_answer_with_step_across_worker_changes(self):
        self.check_outcomes([row for row in cases() if row[0] in TOP_CASES])


if __name__ == '__main__':
    unittest.main()
