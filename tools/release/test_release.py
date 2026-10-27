"""Public release evidence rejects stale replays and unapproved changes."""
import copy
import os
import shutil
import subprocess
from types import SimpleNamespace
import unittest

from pathlib import Path
import tempfile
from unittest.mock import patch

import release
from release import check_summary, COMPARISONS, LABELS


class ReplayPolicyTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        upgrades, corrections = root / "upgrades.json", root / "corrections.json"
        release.write_json(upgrades, {"format": "catallax-replay-upgrades/2", "upgrades": []})
        release.write_json(corrections, {"format": "catallax-replay-corrections/1",
                                        "baseline_sha": "a" * 40, "corrections": []})
        patched = patch.multiple(release, UPGRADES=upgrades, CORRECTIONS=corrections)
        patched.start()
        self.addCleanup(patched.stop)

    def summary(self):
        run = {"cases": 477, "manifest_cases": 477, "crashes": 0, "cache_hits": 0,
               "outcomes": {"answered": 396, "numbers.operation-budget": 81},
               "work_total": 1000, "elapsed_ms": 2000, "binary_sha256": "b" * 64}
        summary = {"format": "catallax-release-replay/3", "base_sha": "a" * 40,
                   "candidate_sha": "b" * 40, "upgrades_sha256": release.sha256(release.UPGRADES),
                   "corrections_sha256": release.sha256(release.CORRECTIONS),
                   "runs": {name: copy.deepcopy(run) for name in LABELS},
                   "comparisons": {name: {"cases": 477, "equal": True, "failures": []} for name in COMPARISONS}}
        summary["comparisons"][COMPARISONS[0]]["upgrades"] = 0
        summary["comparisons"][COMPARISONS[0]]["corrections"] = 0
        return summary

    def test_outcome_totals_must_cover_all_manifest_cases(self):
        summary = self.summary()
        summary["runs"]["candidate-w1"]["outcomes"]["answered"] -= 1
        with self.assertRaisesRegex(ValueError, "incomplete replay outcomes"):
            check_summary(summary, "b" * 40)

    def test_stale_commit_cannot_release(self):
        with self.assertRaisesRegex(ValueError, "selected commit"):
            check_summary(self.summary(), "c" * 40)

    def test_historical_dry_run_retains_source_and_still_requires_equality(self):
        summary = self.summary()
        result = check_summary(summary, "c" * 40, allow_historical=True)
        self.assertEqual(result["candidate_sha"], "b" * 40)
        summary["comparisons"][COMPARISONS[0]]["equal"] = False
        with self.assertRaisesRegex(ValueError, "comparison differs"):
            check_summary(summary, "c" * 40, allow_historical=True)

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

    def test_summary_requires_the_committed_upgrade_digest_and_count(self):
        summary = self.summary()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "upgrades.json"
            approval = {"index": 0, "request_sha256": "1" * 64,
                        "refusal_sha256": "2" * 64, "answer_sha256": "3" * 64}
            release.write_json(path, {"format": "catallax-replay-upgrades/2", "upgrades": [approval]})
            with patch.object(release, "UPGRADES", path):
                with self.assertRaisesRegex(ValueError, "upgrades file differs"):
                    check_summary(summary, "b" * 40)
                summary["upgrades_sha256"] = release.sha256(path)
                summary["comparisons"][COMPARISONS[0]].update(equal=False, upgrades=1)
                for name in LABELS[1:]:
                    summary["runs"][name]["outcomes"] = {"answered": 397, "numbers.operation-budget": 80}
                self.assertEqual(check_summary(summary, "b" * 40)["comparisons"][COMPARISONS[0]]["upgrades"], 1)
                summary["comparisons"][COMPARISONS[0]]["upgrades"] = 2
                with self.assertRaisesRegex(ValueError, "comparison differs"):
                    check_summary(summary, "b" * 40)

    def test_upgrade_file_rejects_duplicate_or_out_of_range_cases(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "upgrades.json"
            approval = {"index": 0, "request_sha256": "1" * 64,
                        "refusal_sha256": "2" * 64, "answer_sha256": "3" * 64}
            for rows in ([approval, approval], [{**approval, "index": 477}]):
                release.write_json(path, {"format": "catallax-replay-upgrades/2", "upgrades": rows})
                with patch.object(release, "UPGRADES", path):
                    with self.assertRaisesRegex(ValueError, "invalid or duplicate"):
                        release.read_upgrades()


    def correction(self):
        return {"index": 0, "class": "invalid-to-valid-step", "request_sha256": "1" * 64,
                "before_sha256": "2" * 64, "after_sha256": "3" * 64, "manifest_fields": [],
                "sources": [{"pr": "https://github.com/BohnBawerick/catallax/pull/123", "commit": "e" * 40}],
                "evidence": {"sha256": "4" * 64, "before_step_sha256": "5" * 64, "after_step_sha256": "6" * 64}}

    def save_corrections(self, rows):
        release.write_json(release.CORRECTIONS, {"format": "catallax-replay-corrections/1",
                                                "baseline_sha": "a" * 40, "corrections": rows})

    def test_correction_digest_count_and_strict_worker_equality(self):
        summary = self.summary()
        self.save_corrections([self.correction()])
        with self.assertRaisesRegex(ValueError, "corrections file differs"):
            check_summary(summary, "b" * 40)
        summary["corrections_sha256"] = release.sha256(release.CORRECTIONS)
        summary["comparisons"][COMPARISONS[0]].update(equal=False, corrections=1)
        self.assertFalse(check_summary(summary, "b" * 40)["comparisons"][COMPARISONS[0]]["equal"])
        for comparison in COMPARISONS[1:]:
            changed = copy.deepcopy(summary)
            changed["comparisons"][comparison]["equal"] = False
            with self.assertRaisesRegex(ValueError, "comparison differs"):
                check_summary(changed, "b" * 40)
        summary["comparisons"][COMPARISONS[0]]["corrections"] = 2
        with self.assertRaisesRegex(ValueError, "comparison differs"):
            check_summary(summary, "b" * 40)

    def test_selected_baseline_cannot_change(self):
        summary = self.summary()
        summary["base_sha"] = "c" * 40
        with self.assertRaisesRegex(ValueError, "baseline differs"):
            check_summary(summary, "b" * 40)

    def test_invalid_or_unbounded_correction_approvals_are_rejected(self):
        mutations = (
            {"index": True}, {"index": -1}, {"index": 477},
            {"class": "valid-to-valid-step"}, {"request_sha256": "invalid"},
            {"manifest_fields": ["properties.volume"]},
            {"manifest_fields": ["step.max_deviation", "step.max_deviation"]},
            {"sources": []}, {"sources": [{"pr": "local-file", "commit": "e" * 40}]},
            {"evidence": {}}, {"extra": "not allowed"},
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.save_corrections([self.correction() | mutation])
                with self.assertRaises(ValueError):
                    release.read_corrections()
        self.save_corrections([self.correction(), self.correction()])
        with self.assertRaisesRegex(ValueError, "duplicate correction"):
            release.read_corrections()

    def test_upgrade_and_correction_cannot_share_a_position(self):
        self.save_corrections([self.correction()])
        release.write_json(release.UPGRADES, {"format": "catallax-replay-upgrades/2", "upgrades": [
            {"index": 0, "request_sha256": "1" * 64, "refusal_sha256": "2" * 64, "answer_sha256": "3" * 64}]})
        with self.assertRaisesRegex(ValueError, "overlapping replay approvals"):
            check_summary(self.summary(), "b" * 40)


    def test_approval_counts_and_equality_keep_their_json_types(self):
        for field, value in (("corrections", False), ("upgrades", False), ("equal", 1)):
            with self.subTest(field=field):
                summary = self.summary()
                summary["comparisons"][COMPARISONS[0]][field] = value
                with self.assertRaisesRegex(ValueError, "comparison differs"):
                    check_summary(summary, "b" * 40)

    def test_uncommitted_approval_bytes_cannot_generate_a_release_aggregate(self):
        root = release.UPGRADES.parent
        relative_driver = Path(".claude/skills/verify-catallax/scripts/control-catallax.py")
        copied_driver = root / relative_driver
        copied_driver.parent.mkdir(parents=True)
        shutil.copyfile(release.ROOT / relative_driver, copied_driver)
        for args in (("init",), ("config", "core.autocrlf", "false"), ("add", "."),
                     ("-c", "user.name=Release test", "-c", "user.email=release@example.invalid",
                      "commit", "-m", "reviewed approvals")):
            subprocess.run(["git", *args], cwd=root, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        for path in (release.UPGRADES, release.CORRECTIONS):
            with self.subTest(file=path.name):
                original = path.read_bytes()
                path.write_bytes(original + b" ")
                try:
                    with patch.object(release, "ROOT", root):
                        with self.assertRaisesRegex(ValueError, "approval files must be committed"):
                            release.replay_summary(SimpleNamespace(runs=[], out=root / "summary.json"))
                    self.assertFalse((root / "summary.json").exists())
                finally:
                    path.write_bytes(original)


class ReleaseSourceTests(unittest.TestCase):
    def test_publication_requires_matching_release_branch_ancestry(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "Cargo.toml").write_text('[workspace.package]\nversion = "0.1.1"\n')
            (root / "docs/milestones").mkdir(parents=True)
            (root / "docs/milestones/M1.md").write_text("# M1\n")

            def git(*args):
                subprocess.run(["git", *args], cwd=root, check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

            git("init")
            git("add", ".")
            git("-c", "user.name=Release test", "-c", "user.email=release@example.invalid",
                "commit", "-m", "candidate")
            git("update-ref", "refs/remotes/origin/main", "HEAD")
            git("update-ref", "refs/remotes/origin/release/0.2.0", "HEAD")
            args = SimpleNamespace(publish=True, tag="v0.1.1")
            with patch.object(release, "ROOT", root), patch.dict(
                    os.environ, GITHUB_OUTPUT=str(root / "outputs")):
                with self.assertRaises(subprocess.CalledProcessError):
                    release.prepare(args)
                git("update-ref", "refs/remotes/origin/release/0.1.1", "HEAD")
                release.prepare(args)
                (root / "Cargo.toml").write_text('[workspace.package]\nversion = "0.1.1"\n# release-only metadata\n')
                git("add", "Cargo.toml")
                git("-c", "user.name=Release test", "-c", "user.email=release@example.invalid",
                    "commit", "-m", "release-only metadata")
                git("update-ref", "refs/remotes/origin/release/0.1.1", "HEAD")
                release.prepare(args)
                with self.assertRaisesRegex(ValueError, "workspace version"):
                    release.prepare(SimpleNamespace(publish=True, tag="v0.2.0"))
                (root / "Cargo.toml").write_text('[workspace.package]\nversion = "0.1.1"\n# later work\n')
                git("add", "Cargo.toml")
                git("-c", "user.name=Release test", "-c", "user.email=release@example.invalid",
                    "commit", "-m", "unreleased work")
                git("update-ref", "refs/remotes/origin/main", "HEAD")
                with self.assertRaises(subprocess.CalledProcessError):
                    release.prepare(args)


if __name__ == "__main__":
    unittest.main()
