#!/usr/bin/env python3
"""Prepare public release evidence without publishing replay inputs."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[2]
LABELS = ("baseline-w1", "candidate-w1", "candidate-w2", "candidate-w4")
COMPARISONS = ("baseline-w1 vs candidate-w1", "candidate-w1 vs candidate-w2", "candidate-w1 vs candidate-w4")
COUNT = 477
UPGRADES = ROOT / "tools/release/replay-upgrades.json"
CORRECTIONS = ROOT / "tools/release/replay-corrections.json"
OUTPUTS = ("answer.step", "display.obj", "print.stl", "answer.igs")
STEP_CORRECTION_FIELDS = {"step.added_edges", "step.added_vertices", "step.max_deviation"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def read_upgrades():
    document = read_json(UPGRADES)
    require(set(document) == {"format", "upgrades"} and
            document["format"] == "catallax-replay-upgrades/2", "unsupported upgrades file")
    rows = document["upgrades"]
    require(isinstance(rows, list), "upgrades must be a list")
    seen = set()
    for row in rows:
        require(set(row) == {"index", "request_sha256", "refusal_sha256", "answer_sha256"},
                "unexpected upgrade fields")
        require(type(row["index"]) is int and 0 <= row["index"] < COUNT and row["index"] not in seen,
                "invalid or duplicate upgrade index")
        seen.add(row["index"])
        require(all(isinstance(row[key], str) and re.fullmatch(r"[0-9a-f]{64}", row[key])
                    for key in ("request_sha256", "refusal_sha256", "answer_sha256")),
                "invalid upgrade digest")
    return rows


def read_corrections():
    document = read_json(CORRECTIONS)
    require(isinstance(document, dict) and set(document) == {"format", "baseline_sha", "corrections"},
            "unexpected corrections fields")
    require(document["format"] == "catallax-replay-corrections/1", "unsupported corrections file")
    require(isinstance(document["baseline_sha"], str) and
            re.fullmatch(r"[0-9a-f]{40}", document["baseline_sha"]), "invalid corrections baseline")
    require(isinstance(document["corrections"], list), "corrections must be a list")
    seen = set()
    for row in document["corrections"]:
        require(isinstance(row, dict) and set(row) == {
            "index", "class", "request_sha256", "before_sha256", "after_sha256",
            "manifest_fields", "sources", "evidence"}, "unexpected correction fields")
        require(type(row["index"]) is int and 0 <= row["index"] < COUNT and row["index"] not in seen,
                "invalid or duplicate correction index")
        seen.add(row["index"])
        require(row["class"] == "invalid-to-valid-step", "unsupported correction class")
        require(all(isinstance(row[key], str) and re.fullmatch(r"[0-9a-f]{64}", row[key])
                    for key in ("request_sha256", "before_sha256", "after_sha256")),
                "invalid correction digest")
        fields = row["manifest_fields"]
        require(isinstance(fields, list) and all(isinstance(field, str) for field in fields) and
                len(fields) == len(set(fields)) and set(fields) <= STEP_CORRECTION_FIELDS,
                "unsupported correction manifest scope")
        sources = row["sources"]
        require(isinstance(sources, list) and sources, "missing correction source")
        for source in sources:
            require(isinstance(source, dict) and set(source) == {"pr", "commit"} and
                    isinstance(source["pr"], str) and re.fullmatch(
                        r"https://github\.com/BohnBawerick/catallax/pull/[1-9][0-9]*", source["pr"]) and
                    isinstance(source["commit"], str) and re.fullmatch(r"[0-9a-f]{40}", source["commit"]),
                    "invalid correction source")
        require(len({source["pr"] for source in sources}) == len(sources), "duplicate correction source")
        evidence = row["evidence"]
        require(isinstance(evidence, dict) and set(evidence) == {
            "sha256", "before_step_sha256", "after_step_sha256"} and all(
                isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value)
                for value in evidence.values()), "invalid correction evidence")
    return document


def answer_digest(driver, directory):
    directory = Path(directory)
    return driver.digest({
        "manifest": driver.comparable(read_json(directory / "manifest.json")),
        "files": {name: sha256(directory / name) if (directory / name).exists() else None
                  for name in OUTPUTS}})


def compare_baseline(driver, baseline, candidate, upgrades, corrections):
    result = driver.compare(baseline, candidate)
    before, after = [read_json(path)["cases"] for path in (baseline, candidate)]
    approvals = {row["index"]: row for row in upgrades}
    repairs = {row["index"]: row for row in corrections}
    require(len(approvals) == len(upgrades) and len(repairs) == len(corrections),
            "duplicate replay approval")
    require(not approvals.keys() & repairs.keys(), "overlapping replay approvals")
    used, corrected = set(), set()
    for failure in result["failures"]:
        require(failure["cause"] in ("manifest", "step-bytes"),
                "replay differs; inspect the local comparison")
        index = failure["index"]
        old, new = [(Path(report).parent / cases[index]["path"]).resolve()
                    for report, cases in ((baseline, before), (candidate, after))]
        old_manifest, new_manifest = [read_json(path / "manifest.json") for path in (old, new)]
        request = read_json(old / "request.json")
        request.pop("out", None)
        if old_manifest["status"] == "refused" and new_manifest["status"] == "answered":
            require(index in approvals, "unlisted upgrade")
            expected = {"index": index, "request_sha256": driver.digest(request),
                        "refusal_sha256": answer_digest(driver, old),
                        "answer_sha256": answer_digest(driver, new)}
            require(approvals[index] == expected, "approved upgrade differs")
            used.add(index)
        elif old_manifest["status"] == new_manifest["status"] == "answered":
            require(index in repairs, "replay differs; unlisted output correction")
            row = repairs[index]
            require(row["class"] == "invalid-to-valid-step", "unsupported correction class")
            require(row["request_sha256"] == driver.digest(request) and
                    row["before_sha256"] == answer_digest(driver, old) and
                    row["after_sha256"] == answer_digest(driver, new), "approved correction differs")
            left, right = [driver.comparable(manifest) for manifest in (old_manifest, new_manifest)]
            left_step, right_step = [manifest.get("step") for manifest in (left, right)]
            require(isinstance(left_step, dict) and isinstance(right_step, dict) and
                    left_step.keys() == right_step.keys() and
                    left_step.get("file") == right_step.get("file") == "answer.step",
                    "correction requires the same STEP publication")
            changed = {"step." + key for key in left_step if left_step[key] != right_step[key]}
            require(changed == set(row["manifest_fields"]) and changed <= STEP_CORRECTION_FIELDS and
                    {key: value for key, value in left.items() if key != "step"} ==
                    {key: value for key, value in right.items() if key != "step"},
                    "correction changes unapproved manifest scope")
            require((old / "answer.step").is_file() and (new / "answer.step").is_file(),
                    "correction requires both STEP files")
            old_step, new_step = [sha256(path / "answer.step") for path in (old, new)]
            require(old_step != new_step and old_step == row["evidence"]["before_step_sha256"] and
                    new_step == row["evidence"]["after_step_sha256"], "correction STEP evidence differs")
            require(all(
                (sha256(old / name) if (old / name).exists() else None) ==
                (sha256(new / name) if (new / name).exists() else None)
                for name in OUTPUTS if name != "answer.step"), "correction changes another output")
            corrected.add(index)
        else:
            require(False, "replay differs; unsupported status transition")
    require(used == set(approvals), "unused upgrade approval")
    require(corrected == set(repairs), "unused correction approval")
    return {"cases": result["cases"], "equal": result["equal"], "failures": [],
            "upgrades": len(used), "corrections": len(corrected)}


def check_summary(summary, commit, allow_historical=False):
    require(set(summary) == {"format", "base_sha", "candidate_sha", "runs", "comparisons",
                             "upgrades_sha256", "corrections_sha256"}, "unexpected replay fields")
    require(summary["format"] == "catallax-release-replay/3", "unsupported replay summary")
    upgrades = read_upgrades()
    require(summary["upgrades_sha256"] == sha256(UPGRADES), "replay upgrades file differs")
    correction_policy = read_corrections()
    corrections = correction_policy["corrections"]
    require(summary["corrections_sha256"] == sha256(CORRECTIONS), "replay corrections file differs")
    require(summary["base_sha"] == correction_policy["baseline_sha"], "corrections baseline differs")
    require(not {row["index"] for row in upgrades} & {row["index"] for row in corrections},
            "overlapping replay approvals")
    require(re.fullmatch(r"[0-9a-f]{40}", summary["base_sha"]), "invalid baseline commit")
    require(re.fullmatch(r"[0-9a-f]{40}", summary["candidate_sha"]), "invalid candidate commit")
    require(allow_historical or summary["candidate_sha"] == commit, "replay is not from the selected commit")
    require(set(summary["runs"]) == set(LABELS), "need baseline and candidate at 1, 2 and 4 workers")
    for name, row in summary["runs"].items():
        require(set(row) == {"cases", "cache_hits", "crashes", "manifest_cases", "outcomes", "work_total", "elapsed_ms", "binary_sha256"}, "unexpected run fields")
        require(row["cases"] == row["manifest_cases"] == COUNT, "need 477 complete manifests")
        require(row["cache_hits"] == row["crashes"] == 0, "replay must be fresh and crash-free")
        require(re.fullmatch(r"[0-9a-f]{64}", row["binary_sha256"]), "invalid binary digest")
        require(sum(row["outcomes"].values()) == COUNT, "incomplete replay outcomes")
        require(row["work_total"] >= 0 and row["elapsed_ms"] >= 0, "invalid run measurement")
    require(len({summary["runs"][name]["binary_sha256"] for name in LABELS[1:]}) == 1, "worker runs use different candidate binaries")
    require(set(summary["comparisons"]) == set(COMPARISONS), "missing replay comparisons")
    baseline = summary["comparisons"][COMPARISONS[0]]
    require(type(baseline.get("equal")) is bool and type(baseline.get("upgrades")) is int and
            type(baseline.get("corrections")) is int and baseline == {
                "cases": COUNT, "equal": not (upgrades or corrections), "failures": [],
                "upgrades": len(upgrades), "corrections": len(corrections)}, "replay comparison differs")
    for name in COMPARISONS[1:]:
        require(summary["comparisons"][name] == {"cases": COUNT, "equal": True, "failures": []},
                "replay comparison differs")
    baseline_answered = summary["runs"]["baseline-w1"]["outcomes"].get("answered", 0)
    candidate_outcomes = summary["runs"]["candidate-w1"]["outcomes"]
    require(candidate_outcomes.get("answered", 0) == baseline_answered + len(upgrades),
            "upgrade answer count differs")
    require(all(summary["runs"][name]["outcomes"] == candidate_outcomes for name in LABELS[2:]),
            "worker outcomes differ")
    return summary


def replay_summary(args):
    spec = importlib.util.spec_from_file_location("catallax_driver", ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py")
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    upgrades = read_upgrades()
    correction_policy = read_corrections()
    for path in (UPGRADES, CORRECTIONS):
        committed = subprocess.check_output(["git", "show", f"HEAD:{path.relative_to(ROOT).as_posix()}"], cwd=ROOT)
        require(committed == path.read_bytes(), "approval files must be committed")
    paths = dict(zip(LABELS, map(Path, args.runs)))
    summaries, commits = {}, {}
    for name, path in paths.items():
        report = read_json(path)
        config = read_json(path.parent / "config.json")
        require(config["workers"] == int(name[-1]), "worker count does not match run label")
        require(report["format"] == "catallax-run/1" and len(report["cases"]) == COUNT, "need 477 replay cases")
        seen_commits, outcomes, work = set(), {}, 0
        for case in report["cases"]:
            directory = Path(case["path"])
            require(driver.valid_artifact(directory, case["key"]), "replay artifact hash check failed")
            manifest = read_json(directory / "manifest.json")
            require(manifest == case["response"]["manifest"], "saved manifest differs from response")
            require(not case["cache_hit"] and manifest["status"] != "internal-error", "cached or crashed replay")
            seen_commits.add(manifest["kernel"]["commit"])
            outcome = (manifest.get("refusal") or {}).get("code", manifest["status"])
            outcomes[outcome] = outcomes.get(outcome, 0) + 1
            work += case["response"]["work"]["used"]
        require(len(seen_commits) == 1, "run mixes source commits")
        commits[name] = seen_commits.pop()
        summaries[name] = {"cases": COUNT, "manifest_cases": COUNT, "crashes": 0, "cache_hits": 0,
                           "outcomes": outcomes, "work_total": work, "elapsed_ms": report["elapsed_ms"],
                           "binary_sha256": report["config"]["binary_sha256"]}
    require(len({commits[name] for name in LABELS[1:]}) == 1, "candidate commits differ")
    require(commits[LABELS[0]] == correction_policy["baseline_sha"], "corrections baseline differs")
    comparisons = {}
    for name in COMPARISONS:
        left, right = name.split(" vs ")
        if name == COMPARISONS[0]:
            result = compare_baseline(driver, paths[left], paths[right], upgrades, correction_policy["corrections"])
        else:
            result = driver.compare(paths[left], paths[right])
            require(result["equal"], "candidate worker replay differs; inspect the local comparison")
        comparisons[name] = result
    summary = {"format": "catallax-release-replay/3", "upgrades_sha256": sha256(UPGRADES),
               "corrections_sha256": sha256(CORRECTIONS), "base_sha": commits[LABELS[0]],
               "candidate_sha": commits[LABELS[1]], "runs": summaries, "comparisons": comparisons}
    write_json(args.out, check_summary(summary, commits[LABELS[1]]))
    print(f"477 cases audited with {len(upgrades)} approved upgrades and "
          f"{len(correction_policy['corrections'])} approved output corrections; public aggregate written to {args.out}")


def check_evidence(directory, commit, allow_historical=False):
    replay = check_summary(read_json(directory / "replay-summary.json"), commit, allow_historical)
    scoreboard = read_json(directory / "target-parts.json")
    require(scoreboard["format"] == "catallax-parts/1" and scoreboard["status"] == "completed", "incomplete scoreboard")
    require(scoreboard["max_work"] == 1_000_000, "scoreboard allowance changed")
    require(all(scoreboard[key] == 0 for key in ("crashes", "comparison_disagreements", "step_failures", "edit_failures")), "scoreboard has correctness failures")
    recipes = {p.relative_to(ROOT / "parts").as_posix(): sha256(p) for p in (ROOT / "parts").rglob("*.json") if p.name != "edit-checks.json"}
    require(len(scoreboard["cases"]) == len(recipes), "incomplete target parts")
    require({row["recipe_file"]: row["sha256"] for row in scoreboard["cases"]} == recipes, "target recipes changed")
    for row in scoreboard["cases"]:
        for source in ("recipe", "field", "compare"):
            require(row[source]["manifest"]["kernel"]["commit"] == commit, "scoreboard is not from the selected commit")
        if row["recipe"]["manifest"].get("step"):
            require(row["step_check"]["status"] == "passed", "STEP import was not checked")
    print(f"Scoreboard matches {commit}; replay describes {replay['candidate_sha']}: 477 replays, {len(recipes)} target recipes")


def prepare(args):
    version = tomllib.loads((ROOT / "Cargo.toml").read_text())["workspace"]["package"]["version"]
    commit = git("rev-parse", "HEAD")
    tag = "v" + version
    if args.publish:
        require(re.fullmatch(r"v(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)", args.tag), "tag must be vX.Y.Z")
        require(args.tag == tag and version != "0.0.0", "tag must match the workspace version")
        subprocess.run(["git", "merge-base", "--is-ancestor", commit, f"origin/release/{version}"], cwd=ROOT, check=True)
        milestone = ROOT / "docs/milestones" / f"M{version.split('.')[1]}.md"
        require(milestone.is_file(), "missing milestone note under docs/milestones/MN.md")
    with Path(os.environ["GITHUB_OUTPUT"]).open("a") as stream:
        stream.write(f"tag={tag}\ncommit={commit}\nalpha={str(version.startswith('0.')).lower()}\n")
    print(f"{'Release' if args.publish else 'Dry run'} {tag} at {commit}")


def checksums(directory):
    archives = {"catallax-cli-" + target + suffix for target, suffix in (
        ("x86_64-unknown-linux-gnu", ".tar.xz"), ("aarch64-unknown-linux-gnu", ".tar.xz"),
        ("x86_64-apple-darwin", ".tar.xz"), ("aarch64-apple-darwin", ".tar.xz"),
        ("universal2-apple-darwin", ".tar.xz"), ("x86_64-pc-windows-msvc", ".zip"))}
    require(archives <= {p.name for p in directory.iterdir()}, "missing platform archives")
    require((directory / "catallax-cli.cdx.xml").stat().st_size > 0, "missing SBOM")
    require((directory / "RELEASE_NOTES.md").stat().st_size > 0, "missing generated notes")
    files = sorted(p for p in directory.iterdir() if p.is_file() and p.name != "SHA256SUMS")
    (directory / "SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in files))
    print(f"SHA256SUMS covers {len(files)} files; the checksum file is attested too")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("replay-summary", help="run on the replay machine; exports aggregates only")
    p.add_argument("runs", nargs=4, metavar="RUN_JSON", help="baseline-w1, candidate-w1, candidate-w2, candidate-w4, in order")
    p.add_argument("--out", required=True)
    p = commands.add_parser("check-replay")
    p.add_argument("summary", type=Path)
    p.add_argument("--commit", required=True)
    p.add_argument("--allow-historical-replay", action="store_true", help="non-publishing dry runs only; retain the recorded replay source commit")
    p = commands.add_parser("check-evidence")
    p.add_argument("directory", type=Path)
    p.add_argument("--commit", required=True)
    p.add_argument("--allow-historical-replay", action="store_true", help="non-publishing dry runs only; scoreboard must still match --commit")
    p = commands.add_parser("prepare")
    p.add_argument("--publish", action="store_true")
    p.add_argument("--tag", default="")
    p = commands.add_parser("checksums")
    p.add_argument("directory", type=Path)
    args = parser.parse_args()
    if args.command == "replay-summary":
        replay_summary(args)
    elif args.command == "check-replay":
        summary = check_summary(read_json(args.summary), args.commit, args.allow_historical_replay)
        print(f"Replay summary complete: source {summary['candidate_sha']}, selected commit {args.commit}")
    elif args.command == "check-evidence":
        check_evidence(args.directory, args.commit, args.allow_historical_replay)
    elif args.command == "prepare":
        prepare(args)
    else:
        checksums(args.directory)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
