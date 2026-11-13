#!/usr/bin/env python3
"""Recheck M3 kernel exports and edit names through the public exchange."""

import argparse
from collections import Counter
import copy
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys

from m5_scoreboard import p0_copy


ROOT = Path(__file__).resolve().parents[2]
DRIVER = ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py"
APP_WORK = 20_000_000
SCOREBOARD_WORK = 1_000_000
PRINT_TOLERANCE = "0.01"
PARTS = {"P2": ("P2-collar.json", "collar", "distance", "0.75", 6),
         "P3": ("P3-filleted-boss.json", "junction", "radius", "2.5", 13)}


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def question(source, problem, allowance):
    return {"format": "catallax-question/1", "source": source,
            "problem": problem, "max_work": allowance}


def prepare(directory):
    directory.mkdir()
    near = json.loads((ROOT / "parts/P0/near-miss.json").read_text())
    for variant in ("base", "translated", "scaled", "seam-shifted"):
        recipe = near if variant == "base" else p0_copy(near, variant)
        write(directory / f"P0-{variant}.json", question("recipe", recipe, SCOREBOARD_WORK))
    for part, (file, step_id, parameter, value, _) in PARTS.items():
        original = json.loads((ROOT / "parts" / file).read_text())
        edited = copy.deepcopy(original)
        next(s for s in edited["steps"] if s["id"] == step_id)[parameter] = value
        for state, recipe in (("original", original), ("edited", edited)):
            write(directory / f"{part}-{state}-step.json", question("recipe", recipe, APP_WORK))
            write(directory / f"{part}-{state}-print.json", question("print", {
                "format": "catallax-print/1", "recipe": recipe,
                "tolerance": PRINT_TOLERANCE}, APP_WORK))
        write(directory / f"{part}-names.json", question("name-edit", {
            "format": "catallax-name-edit/1", "recipe": original,
            "edited_recipe": edited}, APP_WORK))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_stl(path, printing):
    data = path.read_bytes()
    count, = struct.unpack_from("<I", data, 80)
    assert len(data) == 84 + 50 * count and count == printing["triangles"]
    edges = Counter()
    vertices = set()
    for i in range(count):
        values = struct.unpack_from("<12fH", data, 84 + 50 * i)
        assert all(math.isfinite(value) for value in values[:12])
        triangle = [tuple(values[j:j+3]) for j in (3, 6, 9)]
        assert len(set(triangle)) == 3
        vertices.update(triangle)
        edges.update(zip(triangle, triangle[1:] + triangle[:1]))
    assert all(n == 1 and edges[(b, a)] == 1 for (a, b), n in edges.items())
    assert digest(path) == printing["sha256"]
    assert len(vertices) == printing["vertices"]
    assert printing["watertight"] and printing["route"] == "brep"
    assert printing["completeness"] == "proven-two-sided"
    assert 0 <= Fraction(printing["deviation_bound_mm"]) <= Fraction(PRINT_TOLERANCE)
    return {"status": "passed", "sha256": digest(path), "triangles": count,
            "vertices": len(vertices), "deviation_bound_mm": printing["deviation_bound_mm"]}


def audit_step(path, recipe):
    # Reuse the replay audit's direct STEP translation and numerical margin.
    sys.path.insert(0, str(ROOT / "tools/checks"))
    sys.path.insert(0, str(ROOT / "tools/step-regressions"))
    from raw_step import solid_volume
    from test_loopless_faces import raw_solids, unique
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT
    from OCP.gp import gp_Pnt

    part = recipe["name"].split()[0]
    expected_faces = 6 if part == "P0" else PARTS[part][4]
    solids = raw_solids(path)
    assert solids and all(BRepCheck_Analyzer(s).IsValid() for s in solids)
    faces = sum(len(unique(s, TopAbs_FACE)) for s in solids)
    assert faces == expected_faces, (faces, expected_faces)
    volumes = [solid_volume(s) for s in solids]
    assert all(v is not None for v in volumes)
    volume = sum(v[0] for v in volumes)
    margin = sum(v[1] for v in volumes)
    if part == "P0":
        a, b = recipe["steps"][:2]
        radius, height = float(Fraction(a["radius"])), float(Fraction(a["height"]))
        p, q = ([float(Fraction(x)) for x in s["pose"]["translation"]] for s in (a, b))
        expected_volume = 2 * math.pi * radius**2 * height
        probes = [(tuple([c[0], c[1], c[2] + height / 2]), True) for c in (p, q)]
        probes.append((((p[0] + q[0]) / 2, p[1], p[2] + height / 2), False))
    elif part == "P2":
        distance = float(Fraction(recipe["steps"][-1]["distance"]))
        expected_volume = math.pi * (1296 - 24 * distance**2 + 2 * distance**3 / 3)
        probes = [((8, 6, 0), True), ((0, 6, 0), False), ((13, 6, 0), False)]
    else:
        radius = Fraction(next(s for s in recipe["steps"] if s["id"] == "junction")["radius"])
        expected_volume = 14428.640126047918 if radius == 2 else 14455.207561827563
        probes = [((20, 15, 21), True), ((20, 15, 23), False),
                  ((28.4, 15, 11), radius == Fraction(5, 2))]
    assert abs(volume - expected_volume) <= margin, (volume, expected_volume, margin)
    for point, inside in probes:
        states = [BRepClass3d_SolidClassifier(s, gp_Pnt(*point), 1e-10).State() for s in solids]
        assert (any(s == TopAbs_IN for s in states) if inside else all(s == TopAbs_OUT for s in states)), (point, states)
    return {"status": "passed", "sha256": digest(path), "faces": faces,
            "solids": len(solids), "volume_mm3": volume, "reference_volume_mm3": expected_volume,
            "numerical_margin_mm3": margin, "membership_probes": len(probes)}


def summarize(run):
    report = json.loads(run.read_text())
    rows = []
    for case in report["cases"]:
        directory = Path(case["path"])
        request = json.loads((directory / "request.json").read_text())
        problem = request["problem"]
        source = request["source"]
        recipe = problem.get("edited_recipe", problem.get("recipe", problem))
        row = {"name": recipe["name"], "source": source, "max_work": request["max_work"],
               "work": case["response"]["work"], "request_sha256": digest(directory / "question.json")}
        part = recipe["name"].split()[0]
        if part in PARTS:
            _, step_id, parameter, edited_value, _ = PARTS[part]
            value = next(s for s in recipe["steps"] if s["id"] == step_id)[parameter]
            row["parameter"] = {"step": step_id, "name": parameter, "value": value}
            row["state"] = "edited" if Fraction(value) == Fraction(edited_value) else "original"
        manifest = case["response"]["manifest"]
        row.update(status=manifest["status"], refusal=manifest.get("refusal"))
        if manifest["status"] == "answered":
            if source == "recipe":
                assert manifest["exact"] and not manifest["marks"]
                row["step"] = audit_step(directory / "answer.step", recipe)
                row["kernel_volume"] = "not-requested"
            elif source == "print":
                row["printing"] = audit_stl(directory / "print.stl", manifest["printing"])
            elif source == "name-edit":
                names = manifest["name_edit"]
                row["names"] = names
                assert names["status"] == "passed", names
        else:
            assert manifest["status"] == "refused", manifest
            assert not any((directory / name).exists() for name in
                           ("answer.step", "print.stl", "display.obj"))
        rows.append(row)
    return {"format": "catallax-m3-kernel-acceptance/1", "cases": rows,
            "binary_sha256": report["config"]["binary_sha256"],
            "kernel": report["cases"][0]["response"]["manifest"]["kernel"],
            "kernel_passed": all(r["status"] == "answered" for r in rows),
            "application": {"status": "not-run", "code": "parts.application-loop-required"},
            "volume_comparison": "independent numerical audit, not a global deviation certificate",
            "field_comparison": "not-requested", "timing": "M26, not qualified here"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out = args.out.resolve()
    if args.out.is_relative_to(ROOT):
        parser.error("--out must be outside the checkout")
    args.out.mkdir()
    doctor = subprocess.run([sys.executable, "-B", str(DRIVER), "doctor", "--binary",
                             str(args.binary.resolve()), "--json"],
                            capture_output=True, text=True, check=True)
    identity = json.loads(doctor.stdout)
    write(args.out / "doctor.json", identity)
    assert identity["about"]["request"] == "catallax-request/1"
    prepare(args.out / "inputs")
    for workers in (1, 2, 4):
        run = args.out / f"workers-{workers}"
        subprocess.run([sys.executable, "-B", str(DRIVER), "run", str(args.out / "inputs"),
                        "--binary", str(args.binary.resolve()), "--out", str(run),
                        "--workers", str(workers), "--fresh", "--json"], check=True,
                       env={**os.environ, "CATALLAX_THREADS": str(workers)})
        if workers != 1:
            subprocess.run([sys.executable, "-B", str(DRIVER), "compare",
                            str(args.out / "workers-1/run.json"), str(run / "run.json"),
                            "--json"], check=True)
    summary = summarize(args.out / "workers-1/run.json")
    summary["workers"] = [1, 2, 4]
    write(args.out / "summary.json", summary)
    print(json.dumps({"kernel_passed": summary["kernel_passed"], "cases": len(summary["cases"])}))
    return 0 if summary["kernel_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
