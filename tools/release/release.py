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


def check_summary(summary, commit, allow_historical=False):
    require(set(summary) == {"format", "base_sha", "candidate_sha", "runs", "comparisons"}, "unexpected replay fields")
    require(summary["format"] == "catallax-release-replay/1", "unsupported replay summary")
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
    for row in summary["comparisons"].values():
        require(row == {"cases": COUNT, "equal": True, "failures": []}, "replay comparison differs")
    return summary


def replay_summary(args):
    spec = importlib.util.spec_from_file_location("catallax_driver", ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py")
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
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
    comparisons = {}
    for name in COMPARISONS:
        left, right = name.split(" vs ")
        result = driver.compare(paths[left], paths[right])
        require(result["equal"], "replay differs; inspect the private comparison on the replay machine")
        comparisons[name] = result
    summary = {"format": "catallax-release-replay/1", "base_sha": commits[LABELS[0]],
               "candidate_sha": commits[LABELS[1]], "runs": summaries, "comparisons": comparisons}
    write_json(args.out, check_summary(summary, commits[LABELS[1]]))
    print(f"477 manifests and strict STEP bytes equal; public aggregate written to {args.out}")


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
        subprocess.run(["git", "merge-base", "--is-ancestor", commit, "origin/main"], cwd=ROOT, check=True)
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
