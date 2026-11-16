#!/usr/bin/env python3
"""The three M10 loft top cuts must refuse at the 20M allowance.

Each request spends its whole allowance before it refuses, so these runs take
most of an hour on a CI runner. The required `step-regressions-slow` job runs
this file beside `step-regressions`; the default `test*.py` discovery pattern
does not collect it. Run it with:

    python3 -B -m unittest discover -s tools/step-regressions -p 'slow_*.py' -v
"""

import unittest

from test_m10_spline_boolean import TOP_CASES, M10Outcomes, cases


class M10TopCutRefusals(M10Outcomes, unittest.TestCase):
    def test_top_requests_refuse_with_no_step_across_worker_changes(self):
        self.check_outcomes([row for row in cases() if row[0] in TOP_CASES])


if __name__ == '__main__':
    unittest.main()
