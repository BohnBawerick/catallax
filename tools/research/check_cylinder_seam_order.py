"""Check authored full-cylinder STEP charts, volume and retained material.

Validity is after StepToTopoDS construction, including its internal pcurve
processing/healing. It is not unhealed validity or a complete source bound.
Generate the 224 files with CATALLAX_ANALYTIC_EVIDENCE set while running the
catallax-step numerical_cylinder_chart integration test on the work PC.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import OCP
from OCP.BRepCheck import BRepCheck_Analyzer

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/step-regressions"))
from check_numerical_step import raw_solid, state, topology, volume

RAW_HELPER_SHA = "5633cdcf16efc13f328ed85b33a738f6a9296e2afb925f4c64bbed7c91025376"


def check(folder):
    assert OCP.__version__ == "8.0.1.1"
    helper = ROOT / "tools/step-regressions/check_numerical_step.py"
    assert hashlib.sha256(helper.read_bytes()).hexdigest() == RAW_HELPER_SHA
    cases = [
        (outer, axis, phase, reverse, start)
        for outer in [False, True]
        for axis in [-1, 1]
        for phase in range(-3, 4)
        for reverse in [False, True]
        for start in range(4)
    ]
    expected = {
        f"cylinder-chart-outer-{str(o).lower()}-axis-{a}-phase-{p}"
        f"-reverse-{str(r).lower()}-start-{s}.step": o
        for o, a, p, r, s in cases
    }
    assert len(expected) == 224
    assert {p.name for p in folder.glob("*.step")} == set(expected)
    rows = []
    for name, outer in sorted(expected.items()):
        path = folder / name
        before = path.read_bytes()
        shape = raw_solid(path)
        assert BRepCheck_Analyzer(shape).IsValid(), (name, "invalid after builder")
        counts = topology(shape)
        assert counts["solids"] == 1
        assert (counts["faces"], counts["edges"], counts["vertices"]) == (
            (3, 3, 2) if outer else (7, 15, 10)
        ), (name, counts)
        measured = volume(shape)
        expected_volume = math.pi if outer else 16 - math.pi
        assert abs(measured - expected_volume) < 1e-9, (name, measured)
        probes = [
            ((0, 0, 0.5), "inside" if outer else "outside"),
            ((1.5, 0, 0.5), "outside" if outer else "inside"),
        ]
        for point, expected_state in probes:
            assert state(shape, point) == expected_state, (name, point, expected_state)
        assert path.read_bytes() == before
        rows.append(
            dict(
                file=name,
                sha256=hashlib.sha256(before).hexdigest(),
                valid_after_builder_internal_healing=True,
                volume_mm3=measured,
                expected_volume_mm3=expected_volume,
                topology=counts,
                material_probes_passed=len(probes),
            )
        )
    return dict(
        format="catallax-cylinder-seam-order-controls/1",
        ocp_version=OCP.__version__,
        raw_helper_sha256=RAW_HELPER_SHA,
        unhealed_validity_established=False,
        validity_includes_occt_builder_internal_healing=True,
        cases=rows,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("authored", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = check(args.authored)
    args.out.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps(dict(checked=len(report["cases"]), complete=True)), flush=True)


if __name__ == "__main__":
    main()
