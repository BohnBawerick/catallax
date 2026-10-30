"""Small stdlib tests. Kernel smoke runs remain outside the repository."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout
import copy
import io
import json
from pathlib import Path
import signal
import tempfile
import threading
import unittest

import edge_hunt as eh

HANG_LIMIT_SECONDS = 600


class GenerationTests(unittest.TestCase):
    def test_seed_and_index_reproduce_geometry_without_shared_rng_state(self):
        first = eh.generate_case(123, 402)
        eh.generate_case(999, 12)
        self.assertEqual(first, eh.generate_case(123, 402))
        self.assertNotEqual(first["question"], eh.generate_case(124, 402)["question"])

    def test_strata_cover_ordered_pairs_operations_chains_and_rational_rotations(self):
        cases = [eh.generate_case(4, i) for i in range(648)]
        cells = {(c["family"], tuple(c["parameters"]["shape_pair"]),
                  c["parameters"]["operation"], c["parameters"]["chain_length"]) for c in cases}
        self.assertEqual(len(cells), 6*9*3*4)
        for case in cases:
            question = case["question"]
            names = set(question["operands"])
            for step in question["steps"]:
                self.assertTrue(set(step["args"]) <= names)
                self.assertNotIn(step["id"], names)
                names.add(step["id"])
            self.assertEqual(question["result"], question["steps"][-1]["id"])
        self.assertTrue(any("/" in q for c in cases
                            for q in c["parameters"]["global_quaternion"]))
        self.assertEqual({eh.generate_case(4, 648*i)["parameters"]["scale_exponent"]
                          for i in range(7)}, set(eh.SCALES))

    def test_touching_and_near_miss_have_exact_support_separation(self):
        for family_index in (1, 5):
            for index in (family_index, 162+family_index, 324+family_index, 1296+family_index):
                case = eh.generate_case(91, index)
                q = case["question"]
                params = case["parameters"]
                scale = eh.F(10)**params["scale_exponent"]
                centers = []
                for operand in (q["operands"]["a"], q["operands"]["b"]):
                    if operand["shape"] == "block":
                        offset = [eh.F(x)/2 for x in operand["size"]]
                    elif operand["shape"] == "cylinder":
                        offset = [0, 0, eh.F(operand["height"])/2]
                    else:
                        offset = [0, 0, 0]
                    centers.append([eh.F(a)+b for a, b in zip(operand["pose"]["translation"],
                                                                           eh.rotate(operand["pose"]["rotation"], offset))])
                inverse = list(map(eh.F, params["global_quaternion"]))
                inverse[1:] = [-x for x in inverse[1:]]
                separation = eh.rotate(inverse, [b-a for a, b in zip(*centers)])
                expected = sum(map(eh.F, params["radii"]))
                if family_index == 5:
                    expected += eh.F(params["gap"])
                self.assertEqual(separation, [expected*scale, 0, 0])


class ResultsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def row(self, index=0, status="refused"):
        return {**eh.generate_case(1, index), "format": eh.RESULT_FORMAT, "status": status,
                "refusal_code": "boolean.unsupported", "stage": "imprint", "step": "s1",
                "location": None, "wall_ms": 10, "exit_code": 0, "detail": "reason",
                "manifest": {"status": status, "stats": {"wall_ms": 1}},
                "binary_sha256": "build-a", "wall_limit_s": 60, "diagnostics": False}

    def save(self, name, rows):
        path = self.root / name
        path.write_text("".join(eh.encoded(r)+"\n" for r in rows))
        return path

    def test_torn_append_repaired_but_complete_corruption_rejected(self):
        path = self.save("rows", [{"a": 1}])
        with path.open("ab") as stream:
            stream.write(b'{"a":')
        with self.assertRaises(ValueError):
            eh.read_jsonl(path)
        self.assertEqual(eh.read_jsonl(path, repair=True), [{"a": 1}])
        self.assertEqual(path.read_bytes(), b'{"a":1}\n')
        path.write_text('{broken}\n')
        with self.assertRaises(ValueError):
            eh.read_jsonl(path, repair=True)

    def test_resume_rejects_changed_binary_before_executing(self):
        binary = self.root / "binary"
        binary.write_text("not an executable")
        row = self.row()
        catalog = self.save("catalog", [row])
        results = self.save("results", [row])
        args = argparse.Namespace(binary=binary, catalog=catalog, results=results,
                                  wall_limit=60, diagnostics=False, workers=1)
        with self.assertRaisesRegex(ValueError, "resume inputs/build/options differ"):
            eh.run(args)

    def test_budget_reserves_whole_question_and_resume_skips_completed_crash(self):
        binary = self.root / "fail"
        binary.write_text("#!/usr/bin/env python3\nraise SystemExit(70)\n")
        binary.chmod(0o755)
        catalog = self.save("catalog", [eh.generate_case(1, 0)])
        results = self.root / "results"
        args = argparse.Namespace(binary=binary, catalog=catalog, results=results,
                                  wall_limit=HANG_LIMIT_SECONDS, diagnostics=False, workers=1,
                                  cpu_budget=1)
        with redirect_stdout(io.StringIO()):
            self.assertEqual(eh.run(args), 75)
        self.assertEqual(results.read_bytes(), b"")
        args.cpu_budget = 2*HANG_LIMIT_SECONDS
        with redirect_stdout(io.StringIO()):
            self.assertEqual(eh.run(args), 0)
        original = results.read_bytes()
        self.assertEqual(eh.read_jsonl(results)[0]["exit_code"], 70)
        with redirect_stdout(io.StringIO()):
            self.assertEqual(eh.run(args), 0)
        self.assertEqual(results.read_bytes(), original)
        changed = eh.read_jsonl(catalog)
        changed[0]["family"] = "other"
        self.save("catalog", changed)
        with self.assertRaisesRegex(ValueError, "resume inputs/build/options differ"):
            eh.run(args)

    def test_kernel_runs_at_absolute_nice_fifteen(self):
        binary = self.root / "priority"
        binary.write_text("#!/usr/bin/env python3\nimport os\n"
                          "raise SystemExit(os.getpriority(os.PRIO_PROCESS, 0))\n")
        binary.chmod(0o755)
        row = eh.execute(eh.generate_case(1, 0), binary, eh.file_digest(binary),
                         HANG_LIMIT_SECONDS, False, threading.Event())
        self.assertEqual(row["exit_code"], 15)

    def test_duplicate_ids_and_changed_question_hash_rejected(self):
        row = self.row()
        with self.assertRaisesRegex(ValueError, "duplicate"):
            eh.unique_rows([row, row])
        row["question"]["result"] = "a"
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            eh.unique_rows([row])

    def test_later_chain_failure_uses_actual_step_operation_and_lineage(self):
        row = self.row(162)
        row["step"] = "s3"
        cause = eh.cause(row)
        self.assertEqual(cause["operation"], row["question"]["steps"][2]["op"])
        self.assertIn("result(", cause["shape_pair"])
        row["step"] = None
        self.assertTrue(eh.cause(row)["operation"].startswith("unknown-step:"))

    def test_deadline_groups_split_pose_and_scale_and_roll_up_to_same_cause(self):
        rows = [self.row(i) for i in (1, 325, 649, 1297)]
        for row in rows:
            row["refusal_code"] = "numbers.deadline"
        keys = [eh.cause(row) for row in rows]
        self.assertEqual([k["pose"] for k in keys[:2]], ["identity", "rotated"])
        self.assertEqual([k["scale_band"] for k in keys], ["unit", "unit", "tiny", "huge"])
        rows[0]["cpu_ms"] = 4
        summaries = [eh.group_summary(key, members) for key, members in eh.groups(rows).items()]
        rolled = eh.rollups(rows, summaries)
        self.assertEqual(len(rolled), 1)
        self.assertEqual(rolled[0]["count"], 4)
        self.assertEqual(rolled[0]["cpu_observations"], 1)
        self.assertEqual(rolled[0]["total_cpu_ms"], 4)
        self.assertEqual(rolled[0]["cpu_wall_ratio"], 0.4)
        self.assertIsNone(eh.cpu_summary(rows[1:])["total_cpu_ms"])

    def test_concurrent_child_cpu_is_attributed_to_its_pid_not_other_workers(self):
        programs = []
        for name, body in (
            ("busy", "start=time.process_time()\nwhile time.process_time()-start < 0.08: pass\n"),
            ("idle", "time.sleep(0.25)\n"),
        ):
            binary = self.root / name
            binary.write_text("#!/usr/bin/env python3\nimport time\n" + body + "raise SystemExit(70)\n")
            binary.chmod(0o755)
            programs.append(binary)
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(eh.execute, eh.generate_case(1, i), binary,
                                   eh.file_digest(binary), HANG_LIMIT_SECONDS, False, threading.Event())
                       for i, binary in enumerate(programs)]
            busy, idle = [future.result() for future in futures]
        self.assertGreaterEqual(busy["cpu_ms"], 80)
        self.assertLess(idle["cpu_ms"], busy["cpu_ms"])
        for row in (busy, idle):
            self.assertAlmostEqual(row["cpu_wall_ratio"], row["cpu_ms"]/row["wall_ms"], places=5)

    def test_report_selects_three_smallest_observed_reproducers(self):
        rows = [self.row(i*4536) for i in range(4)]
        # All have the same initial cell and one-step failure cause.
        source = self.save("results", rows)
        out = self.root / "report"
        with redirect_stdout(io.StringIO()):
            eh.report(argparse.Namespace(results=source, out=out))
        summary = json.loads((out / "summary.json").read_text())
        group = summary["groups_by_count"][0]
        self.assertEqual(group["count"], 4)
        expected = sorted(rows, key=lambda r: (eh.parameter_size(r["question"]), r["question_sha256"]))[:3]
        self.assertEqual([e["question_id"] for e in group["examples"]], [r["question"]["id"] for r in expected])
        for example, row in zip(group["examples"], expected):
            self.assertEqual(json.loads((out / example["file"]).read_text()), row["question"])

    def test_comparison_detects_output_regression_and_rejects_missing_cases(self):
        before = [self.row(0), self.row(1, "answered")]
        after = copy.deepcopy(before)
        after[0]["status"] = "answered"
        after[0]["manifest"]["status"] = "answered"
        after[1]["status"] = "crash"
        after[1]["exit_code"] = 70
        left, right = self.save("before", before), self.save("after", after)
        out = self.root / "comparison"
        with redirect_stdout(io.StringIO()):
            eh.compare(argparse.Namespace(before=left, after=right, out=out))
        summary = json.loads((out / "comparison.json").read_text())
        self.assertEqual({g["change"] for g in summary["groups"]}, {"vanished", "appeared"})
        self.assertEqual({q["question_id"] for q in summary["changed_questions"]},
                         {r["question"]["id"] for r in before})
        self.save("after", after[:1])
        with self.assertRaisesRegex(ValueError, "identical question sets"):
            eh.compare(argparse.Namespace(before=left, after=right, out=out))

    def test_timing_changes_ignored_but_answered_geometry_changes_detected(self):
        a = self.row(status="answered")
        b = copy.deepcopy(a)
        b["manifest"]["stats"]["wall_ms"] = 999
        b["manifest"]["kernel"] = {"commit": "new"}
        self.assertEqual(eh.outcome(a), eh.outcome(b))
        b["manifest"]["result"] = {"faces": 99}
        self.assertNotEqual(eh.outcome(a), eh.outcome(b))

    def test_crash_missing_manifest_and_timeout_are_distinct(self):
        # Only the timeout case needs a short limit: a kill by the runner is the behavior tested.
        for script, limit, expected, exit_code in (
                ("raise SystemExit(70)", HANG_LIMIT_SECONDS, "crash", 70),
                ("pass", HANG_LIMIT_SECONDS, "crash", 0),
                ("import time; time.sleep(10)", 1, "timeout", -signal.SIGKILL)):
            path = self.root / "program"
            path.write_text("#!/usr/bin/env python3\n" + script + "\n")
            path.chmod(0o755)
            row = eh.execute(eh.generate_case(1, 0), path, eh.file_digest(path),
                             limit, False, threading.Event())
            self.assertEqual(row["status"], expected)
            self.assertEqual(row["exit_code"], exit_code)
            self.assertIsNone(row["stats"])


if __name__ == "__main__":
    unittest.main()
