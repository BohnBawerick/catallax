"""Check the fixed periodic import qualification using raw OCCT construction.

Inputs are actual authored exports, the three fixed real import/pose captures,
and the pinned 00002345 source. Its own metre declaration, planar halfspaces
and circle/cylinder carriers independently determine material and volume checks.
No ShapeFix, sewing, SameParameter or standard root transfer is called.
The cut is an OCCT reference operation; it does not qualify a Catallax Boolean.
"""

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from fractions import Fraction
import OCP
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source-revision", required=True)
parser.add_argument("--capture", type=Path, required=True)
parser.add_argument("--imports", type=Path, required=True)
parser.add_argument("--real-source", type=Path, required=True)
parser.add_argument("--example-binary", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument(
    "--raw-helper",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "tools/step-regressions/check_numerical_step.py",
)
args = parser.parse_args()
sys.path.insert(0, str(args.raw_helper.parent))
from check_numerical_step import raw_solid, topology, volume, state

head = args.source_revision
capture = args.capture
probe = args.imports


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def checked(p):
    b = p.read_bytes()
    assert b"catallax numerical approximate" in b and b"complete_body_bound=unestablished" in b
    assert b"UNCERTAINTY_MEASURE_WITH_UNIT" not in b
    s = raw_solid(p)
    assert BRepCheck_Analyzer(s).IsValid(), (p.parent.name, p.name, "raw invalid")
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(s, props)
    c = props.CentreOfMass()
    assert props.Mass() > 0 and p.read_bytes() == b
    return s, props.Mass(), [c.X(), c.Y(), c.Z()], sha(p)


authored = capture / "authored"
s, v, c, h = checked(authored / "periodic-bore.step")
assert abs(v - (16 - math.pi)) < 1e-9
assert topology(s) == dict(
    solids=1, faces=7, edges=15, vertices=10, two_loop_faces=2, euler_with_face_holes=0
)
assert state(s, (0, 0, 0.5)) == "outside" and state(s, (1.5, 0, 0.5)) == "inside"
cut = BRepAlgoAPI_Cut(s, BRepPrimAPI_MakeBox(gp_Pnt(0, -3, -1), 4, 6, 3).Shape())
cut.Build()
assert cut.IsDone() and BRepCheck_Analyzer(cut.Shape()).IsValid()
assert abs(volume(cut.Shape()) - (8 - math.pi / 2)) < 1e-9
assert (
    state(cut.Shape(), (-1.5, 0, 0.5)) == "inside"
    and state(cut.Shape(), (1.5, 0, 0.5)) == "outside"
)
# Independent source carrier inspection and material expectations for the named real part.
source = args.real_source
assert sha(source) == "8b138ab5f82850b19b017f86283b8074a88c7a271c4c4a01576a610bbfa7b40e"
# The pinned source declares METRE via a unit conversion of 1 metre.
assert re.search(r"CONVERSION_BASED_UNIT\(\s*'METRE'", source.read_text())
assert re.search(r"LENGTH_MEASURE\(\s*1\.0+\s*\)", source.read_text())
assert re.search(r"SI_UNIT\(\s*\$,\s*\.METRE\.\s*\)", source.read_text())
entities = {
    int(n): body for n, body in re.findall(r"#(\d+)\s*=\s*(.*?);", source.read_text(), re.S)
}


def refs(b):
    return [int(n) for n in re.findall(r"#(\d+)", b)]


def xyz(n):
    return [
        float(Fraction(x.strip()))
        for x in entities[n].rsplit("(", 1)[1].split(")", 1)[0].split(",")
    ]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    return [x / norm(a) for x in a]


planes = []
for b in entities.values():
    if not b.startswith("ADVANCED_FACE("):
        continue
    surface = entities[refs(b)[-1]]
    if not surface.startswith("PLANE("):
        continue
    placement = entities[refs(surface)[0]]
    ids = refs(placement)
    p = [1000 * x for x in xyz(ids[0])]
    n = unit(xyz(ids[1]))
    sense = 1 if re.sub(r"\s+", "", b).endswith(".T.)") else -1
    planes.append((p, [sense * x for x in n]))
assert len(planes) == 6
groups = []
for p, n in planes:
    group = next((g for g in groups if abs(abs(dot(n, g[0][1])) - 1) < 1e-12), None)
    if group is None:
        groups.append([(p, n)])
    else:
        group.append((p, n))
assert len(groups) == 3 and all(len(g) == 2 for g in groups)
assert all(abs(dot(groups[i][0][1], groups[j][0][1])) < 1e-12 for i in range(3) for j in range(i))
widths = [abs(dot(sub(g[1][0], g[0][0]), g[0][1])) for g in groups]


def frame(surface):
    placement = entities[refs(surface)[0]]
    ids = refs(placement)
    return [1000 * x for x in xyz(ids[0])], unit(xyz(ids[1])), unit(xyz(ids[2]))


circles = [b for b in entities.values() if b.startswith("CIRCLE(")]
cylinders = [b for b in entities.values() if b.startswith("CYLINDRICAL_SURFACE(")]
assert len(circles) == 2 and len(cylinders) == 1
centres = [frame(b)[0] for b in circles]
centre = [(a + b) / 2 for a, b in zip(*centres)]
r = 1000 * float(Fraction(cylinders[0].rsplit(",", 1)[1].strip().rstrip(" )\n")))
radial = frame(cylinders[0])[2]
inside_outer = lambda p: all(dot(sub(p, o), n) < -1e-6 for o, n in planes)
assert inside_outer(centre)
axis = unit(sub(*centres))
side = [
    axis[1] * radial[2] - axis[2] * radial[1],
    axis[2] * radial[0] - axis[0] * radial[2],
    axis[0] * radial[1] - axis[1] * radial[0],
]
stock_candidates = [
    add(
        centre,
        [factor * r * (math.cos(angle) * u + math.sin(angle) * v) for u, v in zip(radial, side)],
    )
    for factor in [1.05, 1.1, 1.25, 1.5]
    for angle in [i * math.pi / 8 for i in range(16)]
]
stock = next((p for p in stock_candidates if inside_outer(p)), None)
assert stock is not None, "no independent stock probe selected"
outside = add(centre, [100 * r * x for x in radial])
assert not inside_outer(outside)
expected = math.prod(widths) - math.pi * r * r * norm(sub(*centres))
rows = []
for key in ["00001818", "00009734", "00002345"]:
    folder = probe / key
    response = json.loads((folder / "response.json").read_text())
    assert response["status"] == "numerical" and response["exact"] is False
    assert response["source_admission"]["method"] == "retained_shared_analytic_boundaries/3"
    assert response["source_admission"]["tolerance_mm"] == "1/100"
    assert response["work_used"] <= response["work_limit"] == 20_000_000
    a, av, ac, ah = checked(folder / "solid-0.step")
    b, bv, bc, bh = checked(folder / "solid-0-edited.step")
    assert topology(a) == topology(b) and abs(av - bv) < 1e-9 * max(1, av)
    assert all(abs(y - x - shift) < 1e-8 for x, y, shift in zip(ac, bc, [3, 5, 7]))
    assert response["solids"][0]["edited_step"]["source_history_stages"] == 1
    repairs = response["source_admission"]["topology_repairs"]
    if key == "00002345":
        assert len(repairs) == 1 and len(response["source_admission"]["intersection_checks"]) == 21
        assert b"=SEAM_CURVE(" in (folder / "solid-0.step").read_bytes()
        assert abs(av - expected) < 1e-9 * max(1, expected)
        assert topology(a) == dict(
            solids=1, faces=7, edges=15, vertices=10, two_loop_faces=2, euler_with_face_holes=0
        )
        for p, wanted in [(centre, "outside"), (stock, "inside"), (outside, "outside")]:
            assert state(a, p) == wanted and state(b, add(p, [3, 5, 7])) == wanted
    else:
        assert not repairs
    rows.append(
        dict(
            model_id=key,
            source_sha256=response["source_sha256"],
            source_step_sha256=ah,
            edited_step_sha256=bh,
            raw_before_and_pose_valid=True,
            pose_volume_and_centroid_checked=True,
            material_probes_checked=key == "00002345",
            topology=topology(a),
            topology_repairs=repairs,
            work_used=response["work_used"],
            work_limit=response["work_limit"],
            import_work=response["work"]["import"],
            stages=response["solids"][0]["work"],
        )
    )
report = dict(
    format="catallax-periodic-numerical-workflow/1",
    source_revision=head,
    ocp_version=OCP.__version__,
    healing=False,
    same_parameter=False,
    real_parts=rows,
    authored=dict(
        raw_valid=True,
        source_step_sha256=h,
        volume_mm3="16-pi",
        bore_void_and_stock_probes_checked=True,
        reference_boolean_cut_valid=True,
        reference_cut_volume_mm3="8-pi/2",
        catallax_boolean_qualification=False,
    ),
    real_2345_independent_carrier_volume_and_material_checked=True,
    source_unit_to_mm=1000,
)
report["checker_sha256"] = sha(Path(__file__))
report["example_binary_sha256"] = sha(args.example_binary)
report["raw_builder_sha256"] = sha(args.raw_helper)
args.out.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
print(json.dumps(report, sort_keys=True), flush=True)
