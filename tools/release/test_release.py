"""Public release evidence rejects stale replays and unapproved changes."""
import copy
import os
import subprocess
from types import SimpleNamespace
import unittest

from pathlib import Path
import tempfile
from unittest.mock import patch

import release
from release import check_summary, COMPARISONS, LABELS


class ReplayPolicyTests(unittest.TestCase):
    def summary(self):
        run = {"cases": 477, "manifest_cases": 477, "crashes": 0, "cache_hits": 0,
               "outcomes": {"answered": 396, "numbers.operation-budget": 81},
               "work_total": 1000, "elapsed_ms": 2000, "binary_sha256": "b" * 64}
        summary = {"format": "catallax-release-replay/2", "base_sha": "a" * 40,
                   "candidate_sha": "b" * 40, "upgrades_sha256": release.sha256(release.UPGRADES),
                   "runs": {name: copy.deepcopy(run) for name in LABELS},
                   "comparisons": {name: {"cases": 477, "equal": True, "failures": []} for name in COMPARISONS}}
        summary["comparisons"][COMPARISONS[0]]["upgrades"] = 0
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
            release.write_json(path, {"format": "catallax-replay-upgrades/1", "upgrades": [approval]})
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
                release.write_json(path, {"format": "catallax-replay-upgrades/1", "upgrades": rows})
                with patch.object(release, "UPGRADES", path):
                    with self.assertRaisesRegex(ValueError, "invalid or duplicate"):
                        release.read_upgrades()


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
