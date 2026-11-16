"""Check finite circular STEP pcurves with small source vertex-incidence errors.

Generate fixtures outside the checkout using the numerical_cylinder_incidence
integration test. Validity includes OCCT builder internal pcurve processing;
it is not unhealed validity or a complete source-fidelity bound.
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
    expected = {
        f"incidence-oblique-{str(o).lower()}-phase-{p}-vertex-{v}-noise-{n}.step": o
        for o in [False, True]
        for p in range(4)
        for v in [1, 2, 4, 5]
        for n in range(7)
    }
    expected.update(
        {
            f"tilted-oblique-{str(o).lower()}-phase-{p}-noise-{n}.step": o
            for o in [False, True]
            for p in range(4)
            for n in range(4)
        }
    )
    assert len(expected) == 256
    assert {p.name for p in folder.glob("*.step")} == set(expected)
    rows = []
    for name, oblique in sorted(expected.items()):
        path = folder / name
        before = path.read_bytes()
        shape = raw_solid(path)
        assert BRepCheck_Analyzer(shape).IsValid(), (name, "invalid after builder")
        counts = topology(shape)
        assert counts == dict(
            solids=1, faces=5, edges=9, vertices=6, two_loop_faces=0, euler_with_face_holes=2
        ), (name, counts)
        measured = volume(shape)
        assert abs(measured - math.pi / 4) < 1e-9, (name, measured)
        probes = [
            ((0.25, 0.25, 0.5), "inside"),
            ((-0.25, 0.25, 0.5), "outside"),
            ((0.9, 0.9, 0.5), "outside"),
        ]
        for (x, y, z), expected_state in probes:
            if oblique:
                x, y = (20 * x - 21 * y) / 29, (21 * x + 20 * y) / 29
            point = (x, y, z)
            assert state(shape, point) == expected_state, (name, point, expected_state)
        assert path.read_bytes() == before
        rows.append(
            dict(
                file=name,
                sha256=hashlib.sha256(before).hexdigest(),
                valid_after_builder_internal_healing=True,
                volume_mm3=measured,
                topology=counts,
                material_probes_passed=len(probes),
            )
        )
    return dict(
        format="catallax-circular-pcurve-incidence-controls/1",
        ocp_version=OCP.__version__,
        raw_helper_sha256=RAW_HELPER_SHA,
        validity_includes_occt_builder_internal_healing=True,
        unhealed_validity_established=False,
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
