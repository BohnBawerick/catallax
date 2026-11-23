"""Replay source retention and missing-source coverage."""

from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import constant_walls


def source(step, text):
    digest = hashlib.sha256(text.encode()).hexdigest()
    return {"id": step, "op": "import_step", "file": {"sha256": digest}}, digest


def recipe(steps):
    return {"format": "catallax-recipe/5", "name": "wall case",
            "steps": steps, "result": steps[-1]["id"]}


class CorpusTests(unittest.TestCase):
    def replay(self, value, sources=()):
        with tempfile.TemporaryDirectory() as temporary:
            captured = Path(temporary) / "captured"
            captured.mkdir()
            (captured / "recipe-case.json").write_text(json.dumps(value))
            for digest, text in sources:
                (captured / f"source-{digest}").write_text(text)
            output = Path(temporary) / "corpus"
            errors = io.StringIO()
            with patch.object(constant_walls.subprocess, "check_output", return_value=b""), \
                    redirect_stderr(errors), redirect_stdout(io.StringIO()):
                status = constant_walls.corpus(captured, output)
            questions = [json.loads(path.read_text()) for path in output.glob("*.json")]
        return status, questions, errors.getvalue()

    def test_each_variant_embeds_only_its_imports(self):
        before, used = source("body", "used source")
        after, later = source("later", "later source")
        value = {
            "format": "catallax-import-recipe/1",
            "recipe": recipe([before, {"id": "wall", "op": "shell", "body": "body"}, after]),
            "files": {used: "used source", later: "later source", "unrelated": "private"},
        }
        status, questions, errors = self.replay(value)
        self.assertFalse(status)
        self.assertEqual(errors, "")
        self.assertEqual(len(questions), 4)
        for question in questions:
            variant = question["problem"]
            expected = {used} if variant["recipe"]["result"] == "wall" else {used, later}
            self.assertEqual(set(variant["files"]), expected)

    def test_missing_import_reports_the_variant_and_excludes_both_routes(self):
        before, used = source("body", "used source")
        after, missing = source("later", "missing source")
        value = recipe([before, {"id": "wall", "op": "shell", "body": "body"}, after])
        status, questions, errors = self.replay(value, [(used, "used source")])
        self.assertEqual(status, 1)
        self.assertIn(missing, errors)
        self.assertIn("later", errors)
        self.assertEqual({q["source"] for q in questions}, {"recipe", "field"})
        self.assertEqual(len(questions), 2)
        for question in questions:
            self.assertEqual(question["problem"]["recipe"]["result"], "wall")
            self.assertEqual(question["problem"]["files"], {used: "used source"})

    def test_a_field_only_wall_without_its_captured_source_is_not_coverage(self):
        before, missing = source("body", "missing source")
        value = recipe([before, {"id": "wall", "op": "shell", "body": "body"}])
        status, questions, errors = self.replay(value)
        self.assertEqual(status, 1)
        self.assertEqual(questions, [])
        self.assertIn(f"result wall missing replayable imports {missing}", errors)
