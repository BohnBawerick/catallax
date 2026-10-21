#!/usr/bin/env python3
"""Generate, run, group and compare catallax milestone questions."""

import argparse
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import random
import re
import resource
import signal
import subprocess
import tempfile
import threading
import time

FAMILIES = ("general", "tangent", "coaxial", "equal-radius", "along-seam", "near-miss")
SHAPES = ("block", "cylinder", "sphere")
OPS = ("unite", "subtract", "intersect")
SCALES = (0, -30, 30, -15, 15, -6, 6)
ROTATIONS = ((1, 0, 0, 0), (1, 1, 1, 1), (F(1, 2), F(1, 3), F(1, 5), F(1, 7)))
RESULT_FORMAT = "catallax-edge-hunt-result/1"


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def digest(value):
    return hashlib.sha256(encoded(value).encode()).hexdigest()


def file_digest(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def rotate(q, p):
    w, x, y, z = map(F, q)
    n = w*w + x*x + y*y + z*z
    matrix = ((w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y)),
              (2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x)),
              (2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z))
    return [sum(a*b for a, b in zip(row, p))/n for row in matrix]


def multiply(a, b):
    w, x, y, z = a
    v, i, j, k = b
    return (w*v-x*i-y*j-z*k, w*i+x*v+y*k-z*j,
            w*j-x*k+y*v+z*i, w*k+x*j-y*i+z*v)


def operand(shape, radius, center, rotation, scale):
    # Centre each local primitive before carrying its pose into world space.
    offset = {"block": [-radius]*3, "cylinder": [0, 0, -radius],
              "sphere": [0, 0, 0]}[shape]
    translation = [scale*(a+b) for a, b in zip(center, rotate(rotation, offset))]
    result = {"shape": shape, "pose": {"rotation": list(map(str, rotation)),
                                         "translation": list(map(str, translation))}}
    if shape == "block":
        result["size"] = [str(2*radius*scale)]*3
    else:
        result["radius"] = str(radius*scale)
        if shape == "cylinder":
            result["height"] = str(2*radius*scale)
    return result


def generate_case(seed, index):
    rng = random.Random(int(digest([seed, index]), 16))
    family = FAMILIES[index % 6]
    pair = (SHAPES[(index//6) % 3], SHAPES[(index//18) % 3])
    operation = OPS[(index//54) % 3]
    length = (1, 3, 4, 5)[(index//162) % 4]
    exponent = SCALES[(index//648) % len(SCALES)]
    scale = F(10)**exponent
    r1, r2 = F(rng.randint(1, 3)), F(rng.randint(1, 3))
    if family == "equal-radius":
        r2 = r1
    gap = F(1, 10**rng.choice((6, 15, 30)))
    centers = {
        "general": [r1/F(2), r1/F(3), r1/F(5)],
        "tangent": [r1+r2, 0, 0],
        "coaxial": [0, 0, r1/F(2)],
        "equal-radius": [r1/F(2), 0, r1/F(3)],
        "along-seam": [r1, r2, 0],
        "near-miss": [r1+r2+gap, 0, 0],
    }
    # Also rotate within each stratum, not just once after all scale strata.
    global_q = ROTATIONS[(index//162 + index//(648*7)) % len(ROTATIONS)]
    relative_q = rng.choice(ROTATIONS[1:]) if family == "general" else ROTATIONS[0]
    center_b = rotate(global_q, centers[family])
    operands = {
        "a": operand(pair[0], r1, [0, 0, 0], global_q, scale),
        "b": operand(pair[1], r2, center_b, multiply(global_q, relative_q), scale),
    }
    steps = [{"id": "s1", "op": operation, "args": ["a", "b"]}]
    for step in range(2, length+1):
        name = f"c{step}"
        shape = SHAPES[(SHAPES.index(pair[1])+step-1) % 3]
        center = rotate(global_q, [F(step-1, 2), F((step-1) % 2, 3), F(1, 5)])
        operands[name] = operand(shape, F(rng.randint(1, 3)), center, global_q, scale)
        steps.append({"id": f"s{step}", "op": OPS[(OPS.index(operation)+step-1) % 3],
                      "args": [f"s{step-1}", name]})
    question = {"format": "exact-brep-question/2", "id": f"eh-{seed}-{index:07d}",
                "operands": operands, "steps": steps, "result": steps[-1]["id"]}
    metadata = {"seed": seed, "index": index, "family": family,
                "parameters": {"shape_pair": list(pair), "operation": operation,
                               "chain_length": length, "scale_exponent": exponent,
                               "radii": list(map(str, (r1, r2))), "gap": str(gap),
                               "global_quaternion": list(map(str, global_q)),
                               "relative_quaternion": list(map(str, relative_q))}}
    return {"question": question, "question_sha256": digest(question), **metadata}


def generate(args):
    root = args.out
    root.mkdir(parents=True, exist_ok=True)
    catalog = root / "catalog.jsonl"
    if catalog.exists():
        raise ValueError(f"catalog already exists: {catalog}")
    with catalog.open("x") as stream:
        for index in range(args.count):
            case = generate_case(args.seed, index)
            filename = case["question"]["id"] + ".json"
            dump(root / filename, case["question"])
            stream.write(encoded({**case, "file": filename}) + "\n")
    print(encoded({"questions": args.count, "catalog": str(catalog)}))


def read_jsonl(path, repair=False):
    """Only a torn final append is recoverable. Reject all other corruption."""
    rows = []
    mode = "r+b" if repair else "rb"
    with Path(path).open(mode) as stream:
        while True:
            offset = stream.tell()
            line = stream.readline()
            if not line:
                break
            if not line.endswith(b"\n"):
                if repair:
                    stream.truncate(offset)
                    break
                raise ValueError(f"unterminated JSONL record: {path}:{len(rows)+1}")
            rows.append(json.loads(line))
    return rows


def diagnostic_location(stderr):
    # The panic hook prints its origin even when catch_unwind writes a manifest.
    matches = re.findall(r"((?:crates/)?catallax-[\w-]+/src/[\w/.-]+\.rs:\d+:\d+)", stderr)
    return matches[0] if matches else None


def execute(case, binary, binary_hash, limit, diagnostics, stop):
    with tempfile.TemporaryDirectory(prefix="edge-hunt-") as tmp:
        root = Path(tmp)
        question_path = root / "question.json"
        dump(question_path, case["question"])
        env = dict(os.environ)
        env["RUST_BACKTRACE"] = "1" if diagnostics else "0"
        inherited_nice = os.getpriority(os.PRIO_PROCESS, 0)
        if inherited_nice > 15:
            raise ValueError("runner niceness exceeds 15; start it at nice 15 or lower")
        start = time.monotonic()
        timed_out = False
        with (root / "stderr").open("w+b") as stderr:
            process = subprocess.Popen(["nice", "-n", str(15-inherited_nice), str(binary), "answer",
                                        str(question_path), "--out", str(root / "out")],
                                       stdout=subprocess.DEVNULL, stderr=stderr,
                                       env=env, start_new_session=True)
            child_usage = None
            def reap(block=False):
                nonlocal child_usage
                pid, status, usage = os.wait4(process.pid, 0 if block else os.WNOHANG)
                if pid:
                    process.returncode = os.waitstatus_to_exitcode(status)
                    child_usage = usage
                return bool(pid)
            try:
                while not reap():
                    remaining = limit - (time.monotonic()-start)
                    if stop.is_set() or remaining <= 0:
                        timed_out = remaining <= 0
                        os.killpg(process.pid, signal.SIGKILL)
                        reap(block=True)
                        break
                    stop.wait(min(remaining, 0.01))
            finally:
                if process.returncode is None:
                    os.killpg(process.pid, signal.SIGKILL)
                    reap(block=True)
            elapsed = round((time.monotonic()-start)*1000, 3)
            cpu_ms = (child_usage.ru_utime + child_usage.ru_stime)*1000
            stderr.seek(0)
            error_text = stderr.read(16384).decode(errors="replace")
        if stop.is_set() and not timed_out:
            return None  # Interrupted work is retried, not recorded as a crash.
        manifest = None
        manifest_error = None
        try:
            manifest = json.loads((root / "out/manifest.json").read_text())
            if (not isinstance(manifest, dict)
                    or manifest.get("format") != "exact-brep-manifest/2"
                    or manifest.get("question_id") != case["question"]["id"]
                    or manifest.get("status") not in ("answered", "refused", "invalid-question", "internal-error")):
                raise ValueError("manifest identity, format or status is invalid")
        except (OSError, ValueError) as error:
            manifest_error = str(error)
            manifest = None
        status = ("timeout" if timed_out else "crash" if process.returncode != 0 or manifest is None
                  or manifest["status"] == "internal-error" else manifest["status"])
        refusal = (manifest or {}).get("refusal") or {}
        return {"format": RESULT_FORMAT, **case, "binary_sha256": binary_hash,
                "wall_limit_s": limit, "diagnostics": diagnostics, "status": status,
                "manifest_status": (manifest or {}).get("status"),
                "refusal_code": refusal.get("code"), "stage": refusal.get("stage"),
                "step": refusal.get("step"), "detail": refusal.get("detail"),
                "stats": (manifest or {}).get("stats"), "wall_ms": elapsed,
                "cpu_ms": round(cpu_ms, 3),
                "cpu_wall_ratio": round(cpu_ms/elapsed, 6) if elapsed else None,
                "exit_code": process.returncode, "stderr": error_text,
                "manifest_error": manifest_error,
                "location": diagnostic_location(error_text) if diagnostics else None,
                "manifest": manifest}


def unique_rows(rows):
    result = {}
    for row in rows:
        key = row["question"]["id"]
        if key in result:
            raise ValueError(f"duplicate question id: {key}")
        if digest(row["question"]) != row["question_sha256"]:
            raise ValueError(f"question hash mismatch: {key}")
        result[key] = row
    return result


def run(args):
    binary = args.binary.resolve(strict=True)
    binary_hash = file_digest(binary)
    cases = unique_rows(read_jsonl(args.catalog))
    args.results.parent.mkdir(parents=True, exist_ok=True)
    # One writer per file, including across independently launched runners.
    import fcntl
    with args.results.open("a+") as output:
        fcntl.flock(output, fcntl.LOCK_EX | fcntl.LOCK_NB)
        previous = unique_rows(read_jsonl(args.results, repair=True))
        for key, row in previous.items():
            if (key not in cases or row["question_sha256"] != cases[key]["question_sha256"]
                    or row["binary_sha256"] != binary_hash
                    or row["wall_limit_s"] != args.wall_limit
                    or row["diagnostics"] != args.diagnostics
                    or any(row[k] != cases[key][k] for k in ("seed", "index", "family", "parameters"))):
                raise ValueError(f"resume inputs/build/options differ at {key}; use a new result file")
        usage_start = resource.getrusage(resource.RUSAGE_CHILDREN)
        def cpu_seconds():
            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            return usage.ru_utime + usage.ru_stime - usage_start.ru_utime - usage_start.ru_stime
        budget_stopped = False
        pending = iter(case for key, case in cases.items() if key not in previous)
        stop = threading.Event()
        def interrupt(_signum, _frame):
            stop.set()
        old_handlers = {sig: signal.signal(sig, interrupt) for sig in (signal.SIGINT, signal.SIGTERM)}
        count = len(previous)
        try:
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                futures = set()
                def submit():
                    nonlocal budget_stopped
                    # Finished children are charged by getrusage. Reserve a full
                    # wall limit for every live single-threaded kernel process.
                    if (args.cpu_budget is not None
                            and cpu_seconds() + (len(futures)+1)*args.wall_limit > args.cpu_budget):
                        budget_stopped = True
                        return
                    case = next(pending, None)
                    if case is not None and not stop.is_set():
                        futures.add(pool.submit(execute, case, binary, binary_hash,
                                                args.wall_limit, args.diagnostics, stop))
                for _ in range(args.workers):
                    submit()
                while futures:
                    ready, _ = wait(futures, return_when=FIRST_COMPLETED)
                    for future in ready:
                        futures.remove(future)
                        row = future.result()
                        if row is not None:
                            output.write(encoded(row) + "\n")
                            output.flush()
                            os.fsync(output.fileno())
                            count += 1
                            if count % 100 == 0:
                                print(f"completed {count}/{len(cases)}", flush=True)
                        submit()
        finally:
            for sig, handler in old_handlers.items():
                signal.signal(sig, handler)
        print(encoded({"completed": count, "total": len(cases), "interrupted": stop.is_set(),
                       "cpu_seconds": round(cpu_seconds(), 3),
                       "cpu_budget_stopped": budget_stopped and count < len(cases)}))
        return 130 if stop.is_set() else 75 if count < len(cases) else 0


def cause(row):
    question = row["question"]
    step = next((s for s in question["steps"] if s["id"] == row["step"]), None)
    if step is None:
        # Do not attribute measurement/export failures or timeouts to a guessed step.
        operation = "unknown-step:" + "/".join(s["op"] for s in question["steps"])
        pair = "/".join(row["parameters"]["shape_pair"])
    else:
        operation = step["op"]
        lineage = {name: operand["shape"] for name, operand in question["operands"].items()}
        for s in question["steps"]:
            if s == step:
                break
            lineage[s["id"]] = "result(" + ",".join(sorted(set(
                re.findall(r"block|cylinder|sphere", "/".join(lineage[a] for a in s["args"]))))) + ")"
        pair = "/".join(lineage[a] for a in step["args"])
    code = row["refusal_code"] or ("runner.wall-limit" if row["status"] == "timeout"
                                    else f"runner.{row['status']}.exit-{row['exit_code']}")
    key = {"operation": operation, "shape_pair": pair, "family": row["family"],
           "code": code, "stage": row["stage"], "location": row["location"]}
    if code in ("numbers.deadline", "runner.wall-limit"):
        key["pose"] = ("rotated" if any(
            any(F(v) != 0 for v in operand.get("pose", {}).get("rotation", ["1", "0", "0", "0"])[1:])
            for operand in question["operands"].values()) else "identity")
        exponent = row["parameters"]["scale_exponent"]
        key["scale_band"] = ("tiny" if exponent <= -15 else "small" if exponent < 0 else
                             "unit" if exponent == 0 else "large" if exponent < 15 else "huge")
    return key


def failed(row):
    return row["status"] != "answered"


def parameter_size(question):
    def size(value):
        if isinstance(value, dict):
            return sum(size(v) for v in value.values())
        if isinstance(value, list):
            return sum(size(v) for v in value)
        if isinstance(value, str):
            try:
                rational = F(value)
                return abs(rational.numerator).bit_length() + rational.denominator.bit_length()
            except ValueError:
                pass
        return 0
    return size(question["operands"]) + len(question["steps"])


def groups(rows):
    result = {}
    for row in rows:
        if failed(row):
            key = encoded(cause(row))
            result.setdefault(key, []).append(row)
    return result


def group_summary(key, rows):
    return {"id": hashlib.sha256(key.encode()).hexdigest()[:16], "cause": json.loads(key),
            "count": len(rows), "total_wall_ms": round(sum(r["wall_ms"] for r in rows), 3),
            "statuses": dict(Counter(r["status"] for r in rows))}


def cpu_summary(rows):
    measured = [r for r in rows if r.get("cpu_ms") is not None]
    cpu_ms = sum(r["cpu_ms"] for r in measured)
    wall_ms = sum(r["wall_ms"] for r in measured)
    return {"cpu_observations": len(measured),
            "total_cpu_ms": round(cpu_ms, 3) if measured else None,
            "cpu_wall_ratio": round(cpu_ms/wall_ms, 6) if measured and wall_ms else None}


def rollups(rows, summaries):
    buckets = {}
    for row in rows:
        if failed(row):
            key = cause(row)
            buckets.setdefault((key["code"], key["stage"]), []).append(row)
    result = []
    for (code, stage), members in buckets.items():
        children = [g for g in summaries if g["cause"]["code"] == code and g["cause"]["stage"] == stage]
        result.append({"id": digest([code, stage])[:16],
                       "code": code, "stage": stage, "count": len(members),
                       "fraction_of_failures": len(members)/sum(failed(r) for r in rows),
                       "total_wall_ms": round(sum(r["wall_ms"] for r in members), 3),
                       **cpu_summary(members), "group_ids": [g["id"] for g in children]})
    return sorted(result, key=lambda g: (-g["count"], -g["total_wall_ms"], g["code"], g["stage"] or ""))


def cause_label(group):
    c = group["cause"]
    return " / ".join(str(c[k] or "-") for k in
                      ("operation", "shape_pair", "family", "code", "stage", "location")) + (
                          " / " + str(c["pose"]) + " / " + str(c["scale_band"]) if "pose" in c else "")


def report(args):
    rows = list(unique_rows(read_jsonl(args.results)).values())
    args.out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for key, members in groups(rows).items():
        summary = group_summary(key, members)
        summary.update(cpu_summary(members))
        directory = args.out / "examples" / summary["id"]
        directory.mkdir(parents=True, exist_ok=True)
        summary["examples"] = []
        for row in sorted(members, key=lambda r: (parameter_size(r["question"]),
                                                  r["question_sha256"]))[:3]:
            path = directory / (row["question"]["id"] + ".json")
            dump(path, row["question"])
            summary["examples"].append({"file": str(path.relative_to(args.out)),
                                         "parameter_size": parameter_size(row["question"]),
                                         "question_id": row["question"]["id"],
                                         "step": row["step"], "detail": row["detail"]})
        summaries.append(summary)
    summaries.sort(key=lambda g: (-g["count"], -g["total_wall_ms"], g["id"]))
    by_time = sorted(summaries, key=lambda g: (-g["total_wall_ms"], -g["count"], g["id"]))
    summary = {"format": "catallax-edge-hunt-summary/1", "questions": len(rows),
               "statuses": dict(Counter(r["status"] for r in rows)),
               "families": dict(Counter(r["family"] for r in rows)),
               "binary_sha256": sorted(set(r["binary_sha256"] for r in rows)),
               "total_wall_ms": round(sum(r["wall_ms"] for r in rows), 3),
               **cpu_summary(rows), "rollups": rollups(rows, summaries),
               "groups_by_count": summaries, "group_ids_by_time": [g["id"] for g in by_time]}
    dump(args.out / "summary.json", summary)
    lines = ["# Edge hunt report", "", f"Questions: {len(rows)}. Outcomes: `{encoded(summary['statuses'])}`.",
             "", "No judge was run. Answered questions have not been checked for geometric correctness.",
             "Total question wall time is not CPU time. Concurrent jobs overlap.",
             f"Child CPU observed for {summary['cpu_observations']}/{len(rows)} questions; "
             f"CPU/wall ratio over those observations: {summary['cpu_wall_ratio']}.",
             "Older records without per-question CPU are unknown, never backfilled from wall time.", "",
             "## Refusal code and stage rollups", "",
             "| Count | Share of failures | Total wall ms | Child CPU ms | CPU/wall | Code | Stage | Fine groups |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for rollup in summary["rollups"]:
        links = f"[{len(rollup['group_ids'])}](#rollup-{rollup['id']})"
        lines.append(f"| {rollup['count']} | {100*rollup['fraction_of_failures']:.1f}% | "
                     f"{rollup['total_wall_ms']} | {rollup['total_cpu_ms']} | {rollup['cpu_wall_ratio']} | "
                     f"{rollup['code']} | {rollup['stage'] or '-'} | {links} |")
    lines += ["", "## Fine groups by code and stage", ""]
    time_rank = {g["id"]: rank for rank, g in enumerate(by_time, 1)}
    count_rank = {g["id"]: rank for rank, g in enumerate(summaries, 1)}
    indexed = {g["id"]: g for g in summaries}
    for rollup in summary["rollups"]:
        lines += [f"<a id=\"rollup-{rollup['id']}\"></a>",
                  f"### {rollup['code']} / {rollup['stage'] or 'unknown stage'}", "",
                  "| Count rank | Time rank | Count | Total wall ms | Child CPU ms | CPU/wall | Cause | Reproducers |",
                  "| --- | --- | --- | --- | --- | --- | --- | --- |"]
        for key in rollup["group_ids"]:
            group = indexed[key]
            examples = ", ".join(f"[{e['question_id']}]({e['file']})" for e in group["examples"])
            lines.append(f"| {count_rank[key]} | {time_rank[key]} | {group['count']} | {group['total_wall_ms']} | "
                         f"{group['total_cpu_ms']} | {group['cpu_wall_ratio']} | {cause_label(group)} | {examples} |")
        lines.append("")
    lines += ["", "Examples are the smallest observed parameter bit-sizes, not delta-minimized cases.",
              "Run from this report directory, using an absolute path to the binary:", "",
              "```sh", "catallax answer examples/<group>/<question>.json --out /tmp/edge-hunt-repro", "```", ""]
    (args.out / "report.md").write_text("\n".join(lines))
    print(encoded({"questions": len(rows), "groups": len(summaries), "report": str(args.out / "report.md")}))


def outcome(row):
    manifest = row["manifest"]
    stable = {k: v for k, v in (manifest or {}).items() if k not in ("stats", "kernel")}
    return {"status": row["status"], "exit_code": row["exit_code"], "manifest": stable}


def compare(args):
    before = unique_rows(read_jsonl(args.before))
    after = unique_rows(read_jsonl(args.after))
    if before.keys() != after.keys():
        raise ValueError("comparison requires identical question sets; finish/resume both runs first")
    for key in before:
        if (before[key]["question_sha256"] != after[key]["question_sha256"]
                or any(before[key][k] != after[key][k] for k in
                       ("seed", "index", "family", "parameters", "wall_limit_s", "diagnostics"))):
            raise ValueError(f"comparison inputs/options differ at {key}")
    left, right = groups(before.values()), groups(after.values())
    changes = []
    for key in sorted(left.keys() | right.keys()):
        a, b = len(left.get(key, [])), len(right.get(key, []))
        state = ("appeared" if a == 0 else "vanished" if b == 0 else
                 "shrank" if b < a else "grew" if b > a else "unchanged")
        changes.append({"cause": json.loads(key), "before": a, "after": b, "change": state,
                        "before_wall_ms": round(sum(r["wall_ms"] for r in left.get(key, [])), 3),
                        "after_wall_ms": round(sum(r["wall_ms"] for r in right.get(key, [])), 3)})
    changed = [{"question_id": key, "before": outcome(before[key]), "after": outcome(after[key])}
               for key in sorted(before) if outcome(before[key]) != outcome(after[key])]
    summary = {"questions": len(before), "groups": changes, "changed_questions": changed}
    args.out.mkdir(parents=True, exist_ok=True)
    dump(args.out / "comparison.json", summary)
    lines = ["# Edge hunt comparison", "", f"Questions: {len(before)}. Changed outcomes: {len(changed)}.",
             "", "Stats and build identity are excluded. Answered output changes are listed, not judged correct.",
             "", "| Change | Before | After | Cause |", "| --- | --- | --- | --- |"]
    for group in changes:
        lines.append(f"| {group['change']} | {group['before']} | {group['after']} | {cause_label(group)} |")
    lines += ["", "## Changed questions", ""]
    for change in changed:
        lines.append(f"- `{change['question_id']}`: {change['before']['status']} -> {change['after']['status']}. "
                     "Full outcomes in comparison.json.")
    (args.out / "comparison.md").write_text("\n".join(lines) + "\n")
    print(encoded({"questions": len(before), "changed_outcomes": len(changed),
                   "group_changes": dict(Counter(g["change"] for g in changes))}))


def positive(text):
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("generate")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--count", type=positive, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(action=generate)
    p = commands.add_parser("run")
    p.add_argument("catalog", type=Path)
    p.add_argument("--binary", type=Path, default=Path("target/release/catallax"))
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--workers", type=positive, default=4)
    p.add_argument("--wall-limit", type=positive, choices=range(1, 61), default=60,
                   help="seconds, at most the judge's 60 s; lower values for smoke tests")
    p.add_argument("--diagnostics", action="store_true", help="enable panic backtraces and location grouping")
    p.add_argument("--cpu-budget", type=positive,
                   help="child CPU seconds for this invocation; reserve active workers' full limits")
    p.set_defaults(action=run)
    p = commands.add_parser("report")
    p.add_argument("results", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(action=report)
    p = commands.add_parser("compare")
    p.add_argument("before", type=Path)
    p.add_argument("after", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(action=compare)
    args = parser.parse_args()
    try:
        return args.action(args) or 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f"edge-hunt: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
