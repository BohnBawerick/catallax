"""Public release evidence must reject incomplete, stale or unequal replays."""
import copy
import unittest

from release import check_summary, COMPARISONS, LABELS


class ReplayPolicyTests(unittest.TestCase):
    def summary(self):
        run = {"cases": 477, "manifest_cases": 477, "crashes": 0, "cache_hits": 0,
               "outcomes": {"answered": 396, "numbers.operation-budget": 81},
               "work_total": 1000, "elapsed_ms": 2000, "binary_sha256": "b" * 64}
        return {"format": "catallax-release-replay/1", "base_sha": "a" * 40,
                "candidate_sha": "b" * 40, "runs": {name: copy.deepcopy(run) for name in LABELS},
                "comparisons": {name: {"cases": 477, "equal": True, "failures": []} for name in COMPARISONS}}

    def test_complete_replay_accepts_ordinary_refusals(self):
        summary = self.summary()
        self.assertEqual(check_summary(summary, "b" * 40), summary)

    def test_stale_commit_cannot_release(self):
        with self.assertRaisesRegex(ValueError, "selected commit"):
            check_summary(self.summary(), "c" * 40)

    def test_reject_incomplete_cached_and_crashed_runs(self):
        for field, value in (("cases", 476), ("manifest_cases", 476), ("cache_hits", 1), ("crashes", 1)):
            with self.subTest(field=field):
                summary = self.summary()
                summary["runs"]["candidate-w4"][field] = value
                with self.assertRaises(ValueError):
                    check_summary(summary, "b" * 40)

    def test_differing_or_missing_comparison_blocks_release(self):
        for change in ("unequal", "failure", "missing"):
            with self.subTest(change=change):
                summary = self.summary()
                if change == "unequal":
                    summary["comparisons"][COMPARISONS[0]]["equal"] = False
                elif change == "failure":
                    summary["comparisons"][COMPARISONS[0]]["failures"] = [{"cause": "step-bytes"}]
                else:
                    del summary["comparisons"][COMPARISONS[2]]
                with self.assertRaises(ValueError):
                    check_summary(summary, "b" * 40)

    def test_candidate_worker_builds_must_match(self):
        summary = self.summary()
        summary["runs"]["candidate-w2"]["binary_sha256"] = "c" * 64
        with self.assertRaisesRegex(ValueError, "different candidate binaries"):
            check_summary(summary, "b" * 40)

    def test_extra_metadata_is_not_a_public_summary(self):
        summary = self.summary()
        summary["input_paths"] = ["not-public"]
        with self.assertRaisesRegex(ValueError, "unexpected replay fields"):
            check_summary(summary, "b" * 40)


if __name__ == "__main__":
    unittest.main()
