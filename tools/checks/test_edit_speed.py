"""Exact operation outcomes and non-increasing accepted work in hosted pairs."""
import contextlib
import io
from pathlib import Path
import tempfile
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

    def test_malformed_json_costs_fail_on_both_sides(self):
        for invalid in (-1, True, False, 1.5, 100.0, "100", None):
            for side in ("baseline", "candidate"):
                with self.subTest(value=invalid, side=side):
                    before, after = call(100), call(0)
                    (before if side == "baseline" else after)["work"]["used"] = invalid
                    with self.assertRaises(ValueError):
                        self.compare([before], [after])

    def test_zero_and_large_integer_costs_remain_valid(self):
        self.compare([call(0)], [call(0)])
        self.compare([call((1 << 63) + 1)], [call(0)])

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


class NativeWorkGateTests(unittest.TestCase):
    def test_zero_and_large_integer_counts_are_parsed_exactly(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "work.txt"
            for count in (0, (1 << 63) + 1):
                with self.subTest(count=count):
                    path.write_text(f"{count}\n")
                    self.assertEqual(edit_speed.read_native_work(path), count)

    def test_malformed_files_cannot_seed_or_pass_a_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "work.txt"
            for reference_work in (None, 100):
                for text in ("-1\n", "true\n", "1.5\n"):
                    with self.subTest(text=text, reference_work=reference_work):
                        path.write_text(text)
                        with self.assertRaises(ValueError):
                            edit_speed.read_native_work(path, reference_work)

    def test_native_references_require_non_increasing_work_including_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "work.txt"
            for count, reference_work in ((0, 0), (99, 100), (100, 100)):
                with self.subTest(count=count, reference_work=reference_work):
                    path.write_text(f"{count}\n")
                    self.assertEqual(edit_speed.read_native_work(path, reference_work), count)
            for count, reference_work in ((1, 0), (101, 100)):
                with self.subTest(count=count, reference_work=reference_work):
                    path.write_text(f"{count}\n")
                    with self.assertRaises(AssertionError):
                        edit_speed.read_native_work(path, reference_work)


class RoundedCertificateTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name) / 'baseline'
        self.head = Path(temporary.name) / 'candidate'
        for path in (self.base, self.head):
            path.mkdir()
            (path / 'answer.step').write_bytes(b'exact STEP bytes')

    def report(self, path, numerator, denominator, faces=13, marks='[]'):
        (path / 'answer.txt').write_text(
            f'Answer {{ value: StepReport {{ max_deviation: Rational {{ '
            f'n: Integer({numerator}), d: Integer({denominator}), }}, '
            f'faces: {faces} }}, marks: {marks} }}')

    def test_identical_and_stronger_exact_certificates_pass(self):
        for numerator, denominator in ((1, 10), (2, 20), (1, 11), (0, 1)):
            with self.subTest(numerator=numerator, denominator=denominator):
                self.report(self.base, 1, 10)
                self.report(self.head, numerator, denominator)
                edit_speed.compare_rounded_step(self.base, self.head)

    def test_a_single_rational_unit_of_regression_cannot_round_away(self):
        denominator = 10 ** 200
        numerator = 6 * 10 ** 190
        self.report(self.base, numerator, denominator)
        self.report(self.head, numerator + 1, denominator)
        with self.assertRaisesRegex(AssertionError, 'deviation bound increased'):
            edit_speed.compare_rounded_step(self.base, self.head)

    def test_tighter_bounds_cannot_hide_changed_marks_or_report_fields(self):
        self.report(self.base, 1, 10)
        for changes in ({'faces': 12}, {'marks': '[changed-certificate]'}):
            with self.subTest(changes=changes):
                self.report(self.head, 1, 11, **changes)
                with self.assertRaisesRegex(AssertionError, 'report fields or marks changed'):
                    edit_speed.compare_rounded_step(self.base, self.head)

    def test_tighter_bounds_cannot_hide_changed_step_bytes(self):
        self.report(self.base, 1, 10)
        self.report(self.head, 1, 11)
        (self.head / 'answer.step').write_bytes(b'changed STEP bytes')
        with self.assertRaisesRegex(AssertionError, 'STEP bytes changed'):
            edit_speed.compare_rounded_step(self.base, self.head)

    def test_missing_duplicate_and_invalid_bounds_fail_closed_on_both_sides(self):
        for side in ('base', 'head'):
            for numerator, denominator in ((-1, 10), (1, 0), (1, -10)):
                with self.subTest(side=side, numerator=numerator, denominator=denominator):
                    self.report(self.base, 1, 10)
                    self.report(self.head, 1, 10)
                    self.report(getattr(self, side), numerator, denominator)
                    with self.assertRaises(ValueError):
                        edit_speed.compare_rounded_step(self.base, self.head)
            for malformed in ('no certificate', 'duplicate'):
                with self.subTest(side=side, malformed=malformed):
                    self.report(self.base, 1, 10)
                    self.report(self.head, 1, 10)
                    path = getattr(self, side) / 'answer.txt'
                    report = path.read_text()
                    path.write_text(report + report if malformed == 'duplicate' else malformed)
                    with self.assertRaises(ValueError):
                        edit_speed.compare_rounded_step(self.base, self.head)


class PerfLookupTests(unittest.TestCase):
    def test_kernel_tools_win_over_the_path_wrapper(self):
        with tempfile.TemporaryDirectory() as root:
            tools = Path(root)
            (tools / "6.8.0-1-azure").mkdir()
            (tools / "6.8.0-1-azure" / "perf").write_bytes(b"")
            self.assertEqual(edit_speed.find_perf(tools, lambda name: "/usr/bin/perf"),
                             tools / "6.8.0-1-azure" / "perf")

    def test_path_perf_when_kernel_tools_hold_none(self):
        # Ubuntu 26.04: linux-tools-* still create the directory, without perf.
        with tempfile.TemporaryDirectory() as root:
            (Path(root) / "7.0.0-43-generic").mkdir()
            self.assertEqual(edit_speed.find_perf(Path(root), lambda name: "/usr/bin/perf"),
                             Path("/usr/bin/perf"))

    def test_missing_perf_names_both_packages(self):
        with tempfile.TemporaryDirectory() as root, \
                self.assertRaisesRegex(SystemExit, "linux-perf"):
            edit_speed.find_perf(Path(root), lambda name: None)


if __name__ == "__main__":
    unittest.main()
