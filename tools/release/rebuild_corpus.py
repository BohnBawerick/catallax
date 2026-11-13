#!/usr/bin/env python3
"""Generate the wider release replay corpus from committed sources.

Every recipe in `parts/`, every question the application loop exported, and one
name-edit case per configured parameter edit become ordered portable questions.
The committed corpus file is the replay input; `--check` verifies that the
committed file still matches the sources. Never prints or writes input
identifiers outside the committed corpus.
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DRIVER = ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py"
PARTS = ROOT / "parts"
EDITS = PARTS / "edit-checks.json"
EXPORTS = ROOT / "crates/catallax-app/tests/fixtures"
CORPUS = ROOT / "tools/release/rebuild-corpus.jsonl"
# The release replay allowance. The 477 saved answers and the name gate use the
# same verification-request allowance; target recipes answer inside it.
ALLOWANCE = 100_000_000


def load_driver():
    import importlib.util
    spec = importlib.util.spec_from_file_location("catallax_rebuild_driver", DRIVER)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    return driver


def recipe_question(path):
    """A target recipe rebuilt through the ordinary recipe adapter."""
    recipe = json.loads(path.read_text())
    return {"format": "catallax-question/1", "source": "recipe",
            "problem": recipe, "max_work": ALLOWANCE}


def exported_question(path):
    """A portable question the application loop wrote, replayed with its own allowance."""
    question = json.loads(path.read_text())
    if set(question) != {"format", "source", "problem", "max_work"} or question["format"] != "catallax-question/1":
        raise ValueError(f"{path}: not a portable application question")
    if question["source"] not in ("question", "recipe", "field", "print", "iges", "compare", "name-edit"):
        raise ValueError(f"{path}: unknown exported question source")
    if type(question["max_work"]) is not int or not 0 <= question["max_work"] <= 2**64 - 1:
        raise ValueError(f"{path}: exported max_work must be an unsigned 64-bit integer")
    return question


def edited_recipe(recipe, changes, relative):
    """The driver's exact parameter-edit transformation (part_edit_check)."""
    edited = json.loads(json.dumps(recipe))
    for edit in changes:
        step = next((s for s in edited["steps"] if s["id"] == edit["step"]), None)
        if step is None:
            raise ValueError(f"{relative}: edit names unknown step {edit['step']}")
        if not edit["path"]:
            raise ValueError(f"{relative}: edit path is empty")
        target = step
        try:
            for key in edit["path"][:-1]:
                target = target[key]
            key = edit["path"][-1]
            previous = target[key]
        except (KeyError, IndexError, TypeError) as error:
            raise ValueError(f"{relative}: invalid edit path") from error
        if not isinstance(previous, str) or not isinstance(edit["value"], str):
            raise ValueError(f"{relative}: edit must replace an exact numeric parameter")
        target[key] = edit["value"]
    return {"format": "catallax-name-edit/1", "recipe": recipe, "edited_recipe": edited}


def edit_question(relative, changes):
    recipe = json.loads((PARTS / relative).read_text())
    return {"format": "catallax-question/1", "source": "name-edit",
            "problem": edited_recipe(recipe, changes, relative), "max_work": ALLOWANCE}


def rows():
    """The ordered corpus: target recipes, application exports, then name edits."""
    result = []
    recipes = sorted(p.relative_to(ROOT).as_posix() for p in PARTS.rglob("*.json")
                     if p.name != "edit-checks.json")
    if not recipes:
        raise ValueError("no target recipes found under parts/")
    for relative in recipes:
        result.append(recipe_question(ROOT / relative))
    exports = sorted(p.relative_to(ROOT).as_posix() for p in EXPORTS.rglob("*.question.json"))
    for relative in exports:
        result.append(exported_question(ROOT / relative))
    document = json.loads(EDITS.read_text())
    if document.get("format") != "catallax-part-edits/1" or not isinstance(document.get("cases"), dict):
        raise ValueError(f"{EDITS}: unsupported edit specification")
    for relative in sorted(document["cases"]):
        if (PARTS / relative).is_file():
            result.append(edit_question(relative, document["cases"][relative]))
    if not result:
        raise ValueError("empty rebuild corpus")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="verify that the committed corpus matches the sources")
    parser.add_argument("--out", default=str(CORPUS))
    args = parser.parse_args()
    driver = load_driver()
    lines = [driver.canonical(row) + "\n" for row in rows()]
    if args.check:
        committed = CORPUS.read_text().splitlines(keepends=True)
        if committed != lines:
            print("committed rebuild corpus differs from the sources; regenerate it",
                  file=__import__("sys").stderr)
            return 1
        print(f"rebuild corpus matches the sources: {len(lines)} cases")
        return 0
    Path(args.out).write_text("".join(lines))
    print(f"rebuild corpus written to {args.out}: {len(lines)} cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
