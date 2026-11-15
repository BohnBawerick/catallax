"""Check actual authored finite-circle STEP files using the frozen OCCT protocol.

Generate fixtures with CATALLAX_ANALYTIC_EVIDENCE pointing outside the checkout
and cargo test -p catallax-step --test numerical_analytic. This checks validity
including OCCT builder internal healing, volume, topology and material probes;
it is not an unhealed-validity or source-deviation proof.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

from export_validity_census import RAW_HELPER_SHA

ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / "tools/step-regressions/check_numerical_step.py"
assert hashlib.sha256(HELPER.read_bytes()).hexdigest() == RAW_HELPER_SHA
sys.path.insert(0, str(HELPER.parent))
import OCP
from OCP.BRepCheck import BRepCheck_Analyzer
from check_numerical_step import raw_solid, state, topology, volume


def check(folder):
    assert OCP.__version__ == "8.0.1.1"
    names = [f"quarter-chart-{phase}-{turn}.step" for phase in range(4) for turn in [-1, 0, 1]]
    names += ["quarter-cylinder.step", "quarter-cylinder-moved.step", "quarter-oblique-anchor.step"]
    rows = []
    for name in names:
        path = folder / name
        before = path.read_bytes()
        assert b"catallax numerical approximate" in before
        assert b"complete_body_bound=unestablished" in before
        shape = raw_solid(path)
        assert BRepCheck_Analyzer(shape).IsValid(), (name, "invalid after builder processing")
        counts = topology(shape)
        assert counts == dict(
            solids=1, faces=5, edges=9, vertices=6, two_loop_faces=0, euler_with_face_holes=2
        ), (name, counts)
        measured = volume(shape)
        assert abs(measured - math.pi / 4) < 1e-10, (name, measured)
        shift = [3, 5, 7] if name == "quarter-cylinder-moved.step" else [0, 0, 0]
        probes = [
            ((0.25, 0.25, 0.5), "inside"),
            ((-0.25, 0.25, 0.5), "outside"),
            ((0.9, 0.9, 0.5), "outside"),
        ]
        for (x, y, z), expected in probes:
            if name == "quarter-oblique-anchor.step":
                x, y = (20 * x - 21 * y) / 29, (21 * x + 20 * y) / 29
            point = [a + b for a, b in zip([x, y, z], shift)]
            assert state(shape, point) == expected, (name, expected)
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
        format="catallax-finite-circle-exchange-controls/1",
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
