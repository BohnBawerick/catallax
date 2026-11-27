import contextlib
import hashlib
import io
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

from check import compare, migrate, parse


class DifferentialComparison(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        self.main = root / "main.jsonl"
        self.candidate = root / "candidate.jsonl"
        self.pins = root / "pins.jsonl"
        self.key = ["catallax-features", "saddle_fillets", "case", "mesh_body_impl", "options", 1]
        self.old = {"key": self.key, "sha256": "a" * 64, "carriers": "saddle", "work": None}
        self.new = {**self.old, "sha256": "b" * 64, "work": [44, 0, 0, 0]}
        self.write(self.main, self.old)
        self.write(self.candidate, self.new)

    @staticmethod
    def write(path, row):
        path.write_text(json.dumps(row) + "\n")

    def compare_quietly(self, pins=None):
        with contextlib.redirect_stdout(io.StringIO()):
            return compare(self.main, self.candidate, pins)

    def test_requires_exact_pins_for_result_and_accounting_changes(self):
        def digest(row):
            return hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()

        self.write(self.pins, {
            "key": self.key,
            "reason": "new saddle display",
            "main_sha256": digest(self.old),
            "candidate_sha256": digest(self.new),
        })
        self.assertTrue(self.compare_quietly())
        self.assertFalse(self.compare_quietly(self.pins))
        for row in ({**self.new, "work": [45, 0, 0, 0]}, self.old):
            self.write(self.candidate, row)
            self.assertTrue(self.compare_quietly(self.pins))

    def test_rejects_legacy_exception_records(self):
        self.write(self.pins, {
            "key": self.key, "reason": "new saddle display", "main": self.old, "candidate": self.new,
        })
        with self.assertRaisesRegex(ValueError, "invalid exception"):
            self.compare_quietly(self.pins)

    def test_unmeasured_work_differs_from_measured_zero(self):
        self.write(self.candidate, {**self.old, "work": [0, 0, 0, 0]})
        self.assertTrue(self.compare_quietly())

    def test_requires_an_explicit_work_field(self):
        self.write(self.candidate, {key: value for key, value in self.old.items() if key != "work"})
        with self.assertRaisesRegex(ValueError, "missing work field"):
            self.compare_quietly()

    def test_refuses_to_compare_rows_from_different_collectors(self):
        self.write(self.candidate, {**self.old, "format": 1})
        with self.assertRaisesRegex(ValueError, "different collectors"):
            self.compare_quietly()

    def test_migration_maps_each_probe_result_to_one_interface_result(self):
        def rows(path, results, **fields):
            lines = []
            for count, (options, digest) in enumerate(results, 1):
                key = [*self.key[:4], options, count]
                lines.append(json.dumps({**self.old, "key": key, "sha256": digest, **fields}))
            path.write_text("\n".join(lines) + "\n")

        probe = [("Debug", "a" * 64), ("Debug", "b" * 64)]
        rows(self.main, probe)
        rows(self.candidate, [("explicit", "c" * 64), ("explicit", "d" * 64)], format=1)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(migrate(self.main, self.candidate))
            # Two different probe results cannot become one interface result.
            rows(self.candidate, [("explicit", "c" * 64), ("explicit", "c" * 64)], format=1)
            self.assertTrue(migrate(self.main, self.candidate))


class DisplayGate(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        self.main, self.head = root / "main", root / "head"
        self.main.mkdir()
        self.head.mkdir()
        self.pins = root / "exceptions.jsonl"
        self.key = ["catallax-mesh", "native_boundaries", "triangle", "mesh_body_impl", "options", 1]
        self.old = self.mesh(1.0)
        self.write(self.main / "1.jsonl", self.old)
        self.candidate(self.old)

    def mesh(self, x):
        # Length-prefixed outcome, position count, then binary64 coordinates.
        encoded = struct.pack("<Q", 2) + b"ok" + struct.pack("<Q9d", 3, 0, 0, 0, x, 0, 0, 0, 1, 0)
        return {"key": self.key, "sha256": hashlib.sha256(encoded).hexdigest(),
                "carriers": "", "format": 2, "work": None, "refused": False}

    @staticmethod
    def write(path, row):
        path.write_text(json.dumps(row) + "\n")

    def candidate(self, row):
        for workers in (1, 4):
            self.write(self.head / f"{workers}.jsonl", row)

    def pin(self, row):
        self.write(self.pins, {"key": self.key, "reason": "move one vertex",
                               "main_sha256": self.old["sha256"], "candidate_sha256": row["sha256"]})

    @staticmethod
    def cli(*args, stdin=None):
        return subprocess.run([sys.executable, "-B", str(Path(__file__).with_name("check.py")), *args],
                              text=True, input=stdin, capture_output=True, check=False)

    def gate(self, *flags):
        return self.cli("compare", self.main, self.head, "--gate", *flags)

    def test_one_moved_vertex_needs_a_pin_for_that_digest_pair(self):
        # The seed main capture predates outcome metadata; its digest is still comparable.
        legacy = {key: value for key, value in self.old.items() if key != "refused"}
        self.write(self.main / "1.jsonl", {**legacy, "format": 1})
        moved = self.mesh(1.01)
        self.candidate(moved)
        rejected = self.gate()
        self.assertEqual(rejected.returncode, 1, rejected.stderr)
        self.assertIn("unlisted digest change", rejected.stdout)
        self.pin(moved)
        accepted = self.gate("--exceptions", self.pins)
        self.assertEqual(accepted.returncode, 0, accepted.stdout + accepted.stderr)
        self.assertIn("1 listed changes", accepted.stdout)
        self.candidate(self.mesh(1.02))
        rejected = self.gate("--exceptions", self.pins)
        self.assertEqual(rejected.returncode, 1, rejected.stderr)

    def test_pin_cannot_waive_a_new_refusal(self):
        refusal = {**self.old, "sha256": "e" * 64, "refused": True}
        self.candidate(refusal)
        self.pin(refusal)
        rejected = self.gate("--exceptions", self.pins)
        self.assertEqual(rejected.returncode, 1, rejected.stderr)
        self.assertIn("new refusal", rejected.stdout)

    def test_pin_cannot_waive_a_worker_difference(self):
        moved = self.mesh(1.01)
        self.write(self.head / "1.jsonl", moved)
        self.pin(moved)
        rejected = self.gate("--exceptions", self.pins)
        self.assertEqual(rejected.returncode, 1, rejected.stderr)
        self.assertIn("worker difference", rejected.stdout)

    def test_worker_matched_work_change_needs_no_pin(self):
        self.candidate({**self.old, "work": [30, 0, 2, 0]})
        accepted = self.gate()
        self.assertEqual(accepted.returncode, 0, accepted.stdout + accepted.stderr)

    def test_old_pins_neither_fail_unchanged_meshes_nor_authorize_new_changes(self):
        self.pin(self.mesh(1.01))
        flags = ("--exceptions", self.pins, "--baseline-exceptions", self.pins)
        accepted = self.gate(*flags)
        self.assertEqual(accepted.returncode, 0, accepted.stdout + accepted.stderr)
        self.candidate(self.mesh(1.01))
        rejected = self.gate(*flags)
        self.assertEqual(rejected.returncode, 1, rejected.stderr)
        self.assertIn("unlisted digest change", rejected.stdout)

    def test_selector_skips_non_display_paths_but_detects_display_code_and_suites(self):
        docs = "ROADMAP.md\0crates/catallax-step/src/writer.rs\0"
        self.assertEqual(self.cli("affected", stdin=docs).stdout, "false\n")
        for display_path in ("crates/catallax-mesh/src/assembly.rs",
                             "crates/catallax-field/src/patch_shell/display.rs",
                             "crates/catallax-app/tests/generator_fillet_display.rs"):
            with self.subTest(path=display_path):
                self.assertEqual(self.cli("affected", stdin=docs + display_path + "\0").stdout, "true\n")

    def test_interface_records_distinguish_refusals_from_meshes_without_rehashing(self):
        fields = f"case\tmesh_body_impl\t\toptions\tnull\t{self.old['sha256']}"
        for version, status in ((1, ""), (2, "\tok"), (2, "\trefused")):
            row = parse(f"displaycapture\t{version}\t{fields}{status}\n", "interface")
            self.assertEqual(row["sha256"], self.old["sha256"])
            self.assertEqual(row.get("refused"), status == "\trefused" if version == 2 else None)
