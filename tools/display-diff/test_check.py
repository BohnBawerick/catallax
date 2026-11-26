import contextlib
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from check import compare, migrate


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
