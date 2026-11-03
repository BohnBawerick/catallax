"""Exact operation outcomes and non-increasing accepted work in hosted pairs."""
import contextlib
import io
import unittest

import edit_speed


def call(work):
    return {
        "operation": "boolean", "status": "answered", "error": None,
        "max_work": 20_000_000, "nested": False, "on_capture_thread": True,
        "work": {"used": work, "exhausted": False}, "wall_ms": 1.0,
        "operation_log": {"format": "catallax-operation-log/1", "operations": [{
            "id": "cut", "operation": "boolean", "status": "answered",
            "marks": [], "refusal": None, "local_tolerances": None,
            "work": {"before": 0, "after": work, "used": work},
            "decisions": {"filter": 20, "refined": 2}, "wall_ms": 1.0,
        }]},
    }


class WorkGateTests(unittest.TestCase):
    def compare(self, before, after):
        with contextlib.redirect_stdout(io.StringIO()):
            edit_speed.compare_calls(before, after, "paired case")

    def test_reductions_pass_but_even_one_extra_unit_fails(self):
        before, after = call(100), call(99)
        after["operation_log"]["operations"][0]["decisions"] = {"filter": 3, "refined": 0}
        self.compare([before], [after])
        with self.assertRaisesRegex(AssertionError, "accepted call work rose"):
            self.compare([before], [call(101)])

    def test_unordered_duplicates_are_matched_without_hiding_a_rise(self):
        self.compare([call(10), call(20)], [call(19), call(9)])
        with self.assertRaisesRegex(AssertionError, "accepted call work rose"):
            self.compare([call(10), call(20)], [call(21), call(8)])

    def test_lower_work_cannot_hide_changed_outcomes(self):
        for field, value in (("marks", [{"code": "changed-certificate"}]),
                             ("refusal", {"code": "changed-refusal"}),
                             ("status", "refused"),
                             ("local_tolerances", {"distance": "1/100"})):
            with self.subTest(field=field):
                after = call(99)
                after["operation_log"]["operations"][0][field] = value
                with self.assertRaisesRegex(AssertionError, "operation answers, marks or refusals changed"):
                    self.compare([call(100)], [after])

    def test_allowance_and_exhaustion_remain_part_of_the_contract(self):
        for mutate in (lambda value: value.update(max_work=20_000_001),
                       lambda value: value["work"].update(exhausted=True)):
            after = call(99)
            mutate(after)
            with self.assertRaisesRegex(AssertionError, "operation answers, marks or refusals changed"):
                self.compare([call(100)], [after])

    def test_identical_outcomes_cannot_hide_a_missing_call(self):
        with self.assertRaisesRegex(AssertionError, "operation multiplicity changed"):
            self.compare([call(100), call(90)], [call(99)])

    def test_unmeasured_calls_cannot_be_confused_with_zero_work(self):
        before = call(100)
        before.update(work=None, max_work=None, operation_log=None)
        self.compare([before], [dict(before)])
        after = dict(before, work={"used": 0, "exhausted": False})
        with self.assertRaisesRegex(AssertionError, "operation answers, marks or refusals changed"):
            self.compare([before], [after])


if __name__ == "__main__":
    unittest.main()
