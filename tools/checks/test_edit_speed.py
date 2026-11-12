"""Exact operation outcomes and non-increasing accepted work in hosted pairs."""
import contextlib
import io
import json
import os
from fractions import Fraction
from pathlib import Path
import re
import subprocess
import tempfile
import textwrap
import unittest
from unittest import mock

import edit_speed

ABSENT = object()


def write_report(path, numerator, denominator, faces=13, marks='[]', definition=ABSENT):
    stated = ' '
    if definition is not ABSENT:
        if definition is None:
            stated = ' definition_deviation: None, '
        else:
            n, d = definition
            stated = (f' definition_deviation: Some(Rational {{ n: Integer({n}), '
                      f'd: Integer({d}), }}), ')
    (path / 'answer.txt').write_text(
        f'Answer {{ value: StepReport {{ max_deviation: Rational {{ '
        f'n: Integer({numerator}), d: Integer({denominator}), }},{stated}'
        f'faces: {faces} }}, marks: {marks} }}')


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

    def report(self, path, numerator, denominator, **fields):
        write_report(path, numerator, denominator, **fields)

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

    def test_a_revision_that_predates_the_definition_field_stays_comparable(self):
        self.report(self.base, 1, 10)
        self.report(self.head, 1, 10, definition=(1, 11))
        edit_speed.compare_rounded_step(self.base, self.head)

    def test_a_legacy_baseline_cannot_hide_a_raised_candidate_repetition(self):
        self.report(self.base, 1, 10)
        self.report(self.head, 1, 10, definition=(1, 10))
        bounds = edit_speed.DefinitionBounds()
        base_definition = edit_speed.rounded_certificate(
            self.base / 'answer.txt')[2]
        self.assertIs(base_definition, False)
        bounds.compare('rounded-step', 'baseline', base_definition, 'pair 0 baseline')
        bounds.compare('rounded-step', 'candidate',
                       edit_speed.compare_rounded_step(self.base, self.head),
                       'pair 0 candidate')
        bounds.compare('rounded-step', 'candidate',
                       edit_speed.compare_rounded_step(self.base, self.head),
                       'pair 1 candidate')
        bounds.compare('rounded-step', 'baseline', base_definition, 'pair 1 baseline')
        self.report(self.head, 1, 10, definition=(1, 5))
        with self.assertRaisesRegex(AssertionError, 'definition deviation increased'):
            bounds.compare('rounded-step', 'candidate',
                           edit_speed.compare_rounded_step(self.base, self.head),
                           'pair 2 candidate')

    def test_definition_bounds_never_increase(self):
        for definition in ((1, 11), (1, 10), (0, 1)):
            with self.subTest(definition=definition):
                self.report(self.base, 1, 10, definition=(1, 10))
                self.report(self.head, 1, 11, definition=definition)
                edit_speed.compare_rounded_step(self.base, self.head)
        self.report(self.base, 1, 10, definition=(1, 11))
        self.report(self.head, 1, 11, definition=(1, 10))
        with self.assertRaisesRegex(AssertionError, 'definition deviation increased'):
            edit_speed.compare_rounded_step(self.base, self.head)

    def test_a_lost_definition_field_or_bound_fails(self):
        self.report(self.base, 1, 10, definition=(1, 10))
        self.report(self.head, 1, 10)
        with self.assertRaisesRegex(AssertionError, 'definition deviation field removed'):
            edit_speed.compare_rounded_step(self.base, self.head)
        self.report(self.head, 1, 10, definition=None)
        with self.assertRaisesRegex(AssertionError, 'definition deviation bound removed'):
            edit_speed.compare_rounded_step(self.base, self.head)

    def test_identical_text_cannot_hide_a_raised_definition_bound(self):
        self.report(self.base, 1, 10, definition=(1, 10))
        self.report(self.head, 1, 10, definition=(2, 10))
        with self.assertRaisesRegex(AssertionError, 'definition deviation increased'):
            edit_speed.compare_rounded_step(self.base, self.head)

    def test_pretty_printed_definition_bounds_leave_the_exact_text_alone(self):
        def pretty(path, definition):
            if definition is ABSENT:
                stated = ''
            elif definition is None:
                stated = '        definition_deviation: None,\n'
            else:
                n, d = definition
                stated = ('        definition_deviation: Some(\n'
                          '            Rational {\n'
                          '                n: Integer(\n'
                          f'                    {n},\n'
                          '                ),\n'
                          '                d: Integer(\n'
                          f'                    {d},\n'
                          '                ),\n'
                          '            },\n'
                          '        ),\n')
            (path / 'answer.txt').write_text(
                'Answer {\n    value: StepReport {\n        max_deviation: Rational {\n'
                '            n: Integer(\n                1,\n            ),\n'
                '            d: Integer(\n                10,\n            ),\n'
                '        },\n' + stated +
                '        faces: 13,\n    },\n    marks: MarkLog {\n'
                '        marks: [],\n    },\n}\n')
        pretty(self.base, ABSENT)
        pretty(self.head, (2, 11))
        edit_speed.compare_rounded_step(self.base, self.head)
        pretty(self.base, (1, 11))
        pretty(self.head, (2, 11))
        with self.assertRaisesRegex(AssertionError, 'definition deviation increased'):
            edit_speed.compare_rounded_step(self.base, self.head)

    def test_invalid_and_duplicate_definition_bounds_fail_closed(self):
        for side in ('base', 'head'):
            for numerator, denominator in ((-1, 10), (1, 0)):
                with self.subTest(side=side, numerator=numerator, denominator=denominator):
                    self.report(self.base, 1, 10, definition=(1, 10))
                    self.report(self.head, 1, 10, definition=(1, 10))
                    self.report(getattr(self, side), 1, 10,
                                definition=(numerator, denominator))
                    with self.assertRaises(ValueError):
                        edit_speed.compare_rounded_step(self.base, self.head)
            with self.subTest(side=side, duplicate=True):
                self.report(self.base, 1, 10, definition=(1, 10))
                self.report(self.head, 1, 10, definition=(1, 10))
                path = getattr(self, side) / 'answer.txt'
                text = path.read_text()
                start, end = text.index('definition_deviation:'), text.index('faces:')
                path.write_text(text[:end] + text[start:end] + text[end:])
                with self.assertRaises(ValueError):
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


class NativeCertificateTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name) / 'baseline'
        self.head = Path(temporary.name) / 'candidate'
        for path in (self.base, self.head):
            path.mkdir()
            for name in edit_speed.FILES:
                (path / name).write_bytes(b'exact ' + name.encode())
            write_report(path, 1, 10)

    def certificate(self, path):
        return edit_speed.native_certificate(path)

    def test_the_reference_carries_the_definition_bound_beside_the_exact_answer(self):
        write_report(self.head, 1, 10, definition=(1, 11))
        base_answer, base_definition = self.certificate(self.base)
        head_answer, head_definition = self.certificate(self.head)
        self.assertEqual(base_answer, head_answer)
        self.assertIs(base_definition, False)
        self.assertEqual(head_definition, Fraction(1, 11))

    def test_a_stated_definition_bound_must_not_rise_or_vanish(self):
        write_report(self.base, 1, 10, definition=(1, 10))
        _, reference = self.certificate(self.base)
        for definition, message in (((2, 10), 'definition deviation increased'),
                                    (ABSENT, 'definition deviation field removed'),
                                    (None, 'definition deviation bound removed')):
            with self.subTest(definition=definition):
                write_report(self.head, 1, 10, definition=definition)
                _, candidate = self.certificate(self.head)
                with self.assertRaisesRegex(AssertionError, message):
                    edit_speed.compare_definition_bound(reference, candidate, self.head)

    def test_every_other_changed_answer_byte_still_changes_the_hash(self):
        write_report(self.head, 1, 10, definition=(1, 11))
        base_answer, _ = self.certificate(self.base)
        head_answer, _ = self.certificate(self.head)
        self.assertEqual(base_answer, head_answer)
        for name in ('answer.step', 'mesh.txt', 'display.obj', 'names.txt',
                     'operations.txt'):
            with self.subTest(name=name):
                (self.head / name).write_bytes(b'changed ' + name.encode())
                changed, _ = self.certificate(self.head)
                self.assertNotEqual(base_answer, changed)
                (self.head / name).write_bytes(b'exact ' + name.encode())
        write_report(self.head, 1, 10, faces=12, definition=(1, 11))
        changed, _ = self.certificate(self.head)
        self.assertNotEqual(base_answer, changed)


class RepeatedDefinitionBoundTests(unittest.TestCase):
    def test_a_legacy_baseline_cannot_bypass_candidate_repetitions(self):
        bounds = edit_speed.DefinitionBounds()
        bounds.compare(('bracket', 27), 'baseline', False, 'pair 0 baseline')
        bounds.compare(('bracket', 27), 'candidate', Fraction(1, 10), 'pair 0 candidate')
        bounds.compare(('bracket', 27), 'candidate', Fraction(1, 10), 'pair 1 candidate')
        bounds.compare(('bracket', 27), 'baseline', False, 'pair 1 baseline')
        bounds.compare(('bracket', 27), 'candidate', Fraction(1, 10), 'pair 2 candidate')
        for definition, message in (
                (Fraction(1, 5), 'definition deviation increased'),
                (None, 'definition deviation bound removed'),
                (False, 'definition deviation field removed')):
            with self.subTest(definition=definition):
                with self.assertRaisesRegex(AssertionError, message):
                    bounds.compare(('bracket', 27), 'candidate', definition,
                                   'pair 3 candidate')
        bounds.compare(('drawing', 24), 'candidate', None, 'another answer')


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


class ApplicationCachePathTests(unittest.TestCase):
    def test_workflow_mapping_reaches_the_build_target_before_and_after_build(self):
        workflow = (Path(__file__).resolve().parents[2] /
                    ".github/workflows/rounded-step-speed.yml").read_text()
        mapping = re.search(
            r"app -> \$\{\{ steps\.([\w-]+)\.outputs\.([\w-]+) \}\}", workflow)
        self.assertIsNotNone(mapping, "application cache must consume a computed relative target")
        step_id, output_name = mapping.groups()
        step = re.search(
            rf"        id: {re.escape(step_id)}\n        run: \|\n((?:          [^\n]*\n)+)",
            workflow)
        self.assertIsNotNone(step)
        script = textwrap.dedent(step[1])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for layout in ("work/project", "workspace with spaces/deeper/project"):
                with self.subTest(layout=layout):
                    workspace = root / layout
                    (workspace / "app").mkdir(parents=True)
                    runner_temp = root / "runner temp"
                    target = runner_temp / "edit-build"
                    for exists in (False, True):
                        with self.subTest(target_exists=exists):
                            if exists:
                                target.mkdir(parents=True)
                            output = root / "outputs"
                            output.write_text("")
                            subprocess.run(["bash", "-eu", "-c", script], check=True,
                                           env={**os.environ, "GITHUB_WORKSPACE": str(workspace),
                                                "RUNNER_TEMP": str(runner_temp),
                                                "GITHUB_OUTPUT": str(output)})
                            values = dict(line.split("=", 1) for line in output.read_text().splitlines())
                            relative = Path(values[output_name])
                            self.assertFalse(relative.is_absolute())
                            self.assertEqual((workspace / "app" / relative).resolve(), target)
                    target.rmdir()



class CpuCaptureTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.binary = self.root / "edit_answer"
        self.binary.write_bytes(b"native probe")
        self.source = self.root / "source/native.json"
        self.source.parent.mkdir()
        self.source.write_text("{}")
        for name in (*edit_speed.FILES, "stats.txt", "work.txt"):
            (self.source.parent / name).write_bytes(b"exact " + name.encode())
        self.dest = self.root / "captures"
        self.recordings = 0
        self.resolves = lambda index: True
        self.corrupt = False
        self.fail_record = False

    def run_probe(self, command, stdout, stderr, **kwargs):
        if command[1] == "record":
            self.recordings += 1
            if self.fail_record:
                raise subprocess.CalledProcessError(9, command)
            answer = Path(command[-1])
            answer.mkdir()
            for name in (*edit_speed.FILES, "stats.txt", "work.txt"):
                (answer / name).write_bytes((self.source.parent / name).read_bytes())
            if self.corrupt:
                (answer / "work.txt").write_bytes(b"changed work")
            Path(command[command.index("-o") + 1]).write_bytes(b"raw CPU capture")
        else:
            self.assertEqual(command[1], "script")
            stdout.write("catallax_boolean::graph::section\n" if self.resolves(self.recordings - 1)
                         else "unresolved native frame\n")

    def capture(self):
        with mock.patch.object(edit_speed, "find_perf", return_value=Path("/perf")), \
                mock.patch.object(edit_speed.subprocess, "check_output", return_value="perf test"), \
                mock.patch.object(edit_speed.subprocess, "run", side_effect=self.run_probe), \
                contextlib.redirect_stdout(io.StringIO()):
            edit_speed.capture_cpu_profile(self.binary, self.source, self.dest)

    def test_five_resolved_captures_are_still_required(self):
        self.capture()
        metadata = json.loads((self.dest / "metadata.json").read_text())
        self.assertEqual(self.recordings, 5)
        self.assertEqual(len(metadata["captures"]), 5)
        self.assertEqual(metadata["discarded_captures"], [])
        for index in range(5):
            self.assertIn("catallax_boolean", (self.dest / str(index) / "stacks.txt").read_text())

    def test_a_missed_phase_is_recaptured_and_its_raw_evidence_is_retained(self):
        self.resolves = lambda index: index != 0
        self.capture()
        metadata = json.loads((self.dest / "metadata.json").read_text())
        self.assertEqual(self.recordings, 6)
        self.assertEqual(len(metadata["captures"]), 5)
        self.assertEqual(metadata["captures"][0], dict(capture=0, attempts=2))
        self.assertEqual(metadata["discarded_captures"], ["unresolved/0-0"])
        discarded = self.dest / "unresolved/0-0"
        self.assertTrue((discarded / "perf.data").is_file())
        self.assertNotIn("catallax_boolean", (discarded / "stacks.txt").read_text())
        edit_speed.check_profile_answer(self.source, discarded)

    def test_missing_frames_after_the_bounded_attempts_still_fail(self):
        self.resolves = lambda index: index < 4
        with self.assertRaisesRegex(AssertionError, "after 5 attempts for capture 4"):
            self.capture()
        metadata = json.loads((self.dest / "metadata.json").read_text())
        self.assertEqual(self.recordings, 9)
        self.assertEqual(len(metadata["captures"]), 4)
        self.assertEqual(len(metadata["discarded_captures"]), 5)

    def test_a_missed_phase_cannot_hide_changed_answers_or_work(self):
        self.resolves = lambda index: False
        self.corrupt = True
        with self.assertRaisesRegex(AssertionError, "profiling changed an exact answer"):
            self.capture()
        self.assertEqual(self.recordings, 1)

    def test_perf_errors_are_not_retried_as_sampling_misses(self):
        self.fail_record = True
        with self.assertRaises(subprocess.CalledProcessError):
            self.capture()
        self.assertEqual(self.recordings, 1)

if __name__ == "__main__":
    unittest.main()
