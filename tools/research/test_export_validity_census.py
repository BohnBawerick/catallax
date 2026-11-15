"""A validity headline requires the fixed population and every body of every file."""

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

from export_validity_census import WORK, observe, population, public_build, summary


class ValidityTests(unittest.TestCase):
    def setUp(self):
        self.hashes = {f"{i:08}": f"{i:064x}" for i in range(3000)}
        outcomes = ["full"] * 400 + ["partial"] * 3 + ["refused"] * 2595 + ["timeout"] * 2
        self.inputs = [
            dict(
                model_id=key,
                source_sha256=sha,
                outcome=outcome,
                coverage_complete="True" if outcome == "full" else "False",
                bodies="2" if index == 0 else "1",
            )
            for index, ((key, sha), outcome) in enumerate(zip(self.hashes.items(), outcomes))
        ]
        self.admissions, self.frozen = population(self.inputs, self.hashes)
        self.rows = [
            dict(
                model_id=key,
                source_sha256=row["source_sha256"],
                outcome="raw_valid",
                exported_body_count=int(row["bodies"]),
                raw_valid_solid_count=int(row["bodies"]),
                exported_step_hashes=["a" * 64] * int(row["bodies"]),
            )
            for key, row in self.frozen.items()
        ]

    def test_complete_population_retains_three_thousand_denominator(self):
        report = summary(self.admissions, self.rows)
        self.assertTrue(report["complete"])
        self.assertEqual(report["denominator"], 3000)
        self.assertEqual(report["frozen_admitted_inputs"], 400)
        self.assertEqual(report["raw_valid_share_percent"], "13.33")

    def test_missing_observation_cannot_publish_a_percentage(self):
        report = summary(self.admissions, self.rows[:-1])
        self.assertFalse(report["complete"])
        self.assertIsNone(report["raw_valid_share_percent"])

    def test_duplicate_changed_or_extra_sources_are_rejected(self):
        for changed in [
            self.rows + [self.rows[0]],
            self.rows
            + [
                dict(
                    model_id="00000400",
                    source_sha256=self.hashes["00000400"],
                    outcome="timeout",
                )
            ],
        ]:
            with self.assertRaises(AssertionError):
                summary(self.admissions, changed)
        changed = copy.deepcopy(self.rows)
        changed[0]["source_sha256"] = "f" * 64
        with self.assertRaisesRegex(AssertionError, "source hash"):
            summary(self.admissions, changed)

    def test_one_valid_body_does_not_qualify_a_two_body_file(self):
        self.rows[0].update(
            exported_body_count=1,
            raw_valid_solid_count=1,
            exported_step_hashes=["a" * 64],
        )
        with self.assertRaisesRegex(AssertionError, "partial export"):
            summary(self.admissions, self.rows)

    def test_invalid_refused_and_timeout_are_separate_and_not_valid(self):
        self.rows[0]["outcome"] = "raw_invalid"
        self.rows[1]["outcome"] = "export_refused"
        self.rows[2]["outcome"] = "raw_reader_refused"
        self.rows[3]["outcome"] = "timeout"
        report = summary(self.admissions, self.rows)
        self.assertTrue(report["complete"])
        self.assertEqual(report["raw_valid_files"], 396)
        self.assertEqual(report["export_outcomes"]["timeout"], 1)

    def test_process_error_and_unknown_result_cannot_publish_complete(self):
        self.rows[0]["outcome"] = "process_error"
        self.assertFalse(summary(self.admissions, self.rows)["complete"])
        self.rows[0]["outcome"] = "not_checked"
        with self.assertRaisesRegex(AssertionError, "unknown"):
            summary(self.admissions, self.rows)

    def test_observer_faults_are_process_errors_and_prevent_complete_census(self):
        for cause in ["observer.write", "observer.reread", "observer.byte-identity"]:
            with self.subTest(cause=cause), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = root / "source.step"
                source.write_bytes(b"classification fixture")
                baseline = self.admissions["00000000"] | dict(
                    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()
                )
                admission = dict(
                    source_sha256=baseline["source_sha256"],
                    outcome="admitted",
                    work_limit=WORK,
                    work_used=1,
                    bodies=2,
                    import_category="certified",
                )
                response = dict(
                    stage="export",
                    outcome="refused",
                    source_sha256=baseline["source_sha256"],
                    work_limit=WORK,
                    work_used=1,
                    files=[],
                    refusal=dict(code=cause),
                    **{"import": admission},
                )
                producer = root / "producer"
                producer.write_text(
                    f"#!{sys.executable}\nfrom pathlib import Path\nimport sys\n"
                    "out = Path(sys.argv[3])\nout.mkdir()\n"
                    f"(out / 'response.json').write_text({json.dumps(response)!r})\n"
                )
                producer.chmod(0o700)
                result = observe(producer, source, baseline, root / "capture", 5)
                self.assertEqual(result["outcome"], "process_error")
                self.assertEqual(result["cause"], cause)
                admissions = self.admissions | {"00000000": baseline}
                report = summary(admissions, [result, *self.rows[1:]])
                self.assertFalse(report["complete"])
                self.assertIsNone(report["raw_valid_share_percent"])

    def test_producer_timeout_checks_source_with_and_without_import_checkpoint(self):
        for checkpoint in [False, True]:
            for mutate in [False, True]:
                with self.subTest(
                    checkpoint=checkpoint, mutate=mutate
                ), tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    source = root / "source.step"
                    source.write_bytes(b"timeout fixture")
                    baseline = self.admissions["00000000"] | dict(
                        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()
                    )
                    admission = dict(
                        stage="import",
                        outcome="admitted",
                        source_sha256=baseline["source_sha256"],
                        work_used=1,
                        import_category="certified",
                        bodies=2,
                    )
                    producer = root / "producer"
                    producer.write_text(
                        f"#!{sys.executable}\nfrom pathlib import Path\nimport sys, time\n"
                        + (f"print({json.dumps(admission)!r}, flush=True)\n" if checkpoint else "")
                        + ("Path(sys.argv[1]).write_bytes(b'changed')\n" if mutate else "")
                        + "time.sleep(60)\n"
                    )
                    producer.chmod(0o700)
                    if mutate:
                        with self.assertRaisesRegex(
                            AssertionError, "source changed during observation"
                        ):
                            observe(producer, source, baseline, root / "capture", 2)
                    else:
                        result = observe(producer, source, baseline, root / "capture", 2)
                        self.assertEqual(result["outcome"], "timeout")
                        self.assertEqual(
                            result["timeout_stage"],
                            "export" if checkpoint else "import",
                        )

    def test_original_population_cannot_change_or_hide_missing_hashes(self):
        with self.assertRaises(AssertionError):
            population(self.inputs[:-1], self.hashes)
        self.inputs[0]["source_sha256"] = "f" * 64
        with self.assertRaisesRegex(AssertionError, "source hash"):
            population(self.inputs, self.hashes)

    def test_build_extract_omits_actual_paths_but_retains_artifact_identity(self):
        fields = "source_revision git_tree kernel_revision binary_sha256 cargo_json_sha256 cargo_lock_sha256 observer_sha256 collector_sha256 raw_helper_sha256 build_command rustc".split()
        record = {key: "a" * 64 for key in fields}
        record.update(
            artifact_path="/home/private/target/example",
            unknown="/Users/private/source",
            cargo_artifact=dict(
                executable="/home/private/target/example",
                profile={"test": False},
                features=[],
                fresh=False,
                target=dict(
                    name="export_import",
                    kind=["example"],
                    crate_types=["bin"],
                    edition="2024",
                    src_path="/home/private/example.rs",
                ),
            ),
        )
        public = public_build(record)
        self.assertEqual(public["cargo_artifact"]["target"]["name"], "export_import")
        self.assertNotIn("/home/", str(public))
        self.assertNotIn("/Users/", str(public))


if __name__ == "__main__":
    unittest.main()
