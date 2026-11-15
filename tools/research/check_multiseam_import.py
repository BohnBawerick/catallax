"""Qualify the fixed two-seam census source 00003390 and its actual pose exports.

Export validity uses raw StepToTopoDS construction without healing. The source
material/volume reference uses normal OCCT root transfer and is labelled
separately. Source-local names and carrier/incidence comparisons check identity;
they are not cross-file replacement correspondence or a global deviation proof.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import OCP
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader

SOURCE_SHA = "9ee5c0862a70c2f7eef35ce8bf61fb1a748fc173fa6c6b606745a9842e36e4a7"
QUOTED = r"'(?:''|[^'])*'"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Entities:
    def __init__(self, text):
        self.bodies = {
            int(n): body.strip()
            for n, body in re.findall(r"#(\d+)\s*=\s*((?:" + QUOTED + r"|[^';])*);", text)
        }

    def kind(self, n):
        return self.bodies[n].split("(", 1)[0].strip()

    def ids(self, kind):
        return [n for n in self.bodies if self.kind(n) == kind]

    def refs(self, n):
        return [int(v) for v in re.findall(r"#(\d+)", re.sub(QUOTED, "''", self.bodies[n]))]

    def name(self, n):
        return re.search(QUOTED, self.bodies[n]).group()[1:-1].replace("''", "'")

    def named(self, kind):
        ids = self.ids(kind)
        result = {self.name(n): n for n in ids}
        assert len(result) == len(ids), (kind, "duplicate source-local identity")
        return result

    def xyz(self, n):
        return [
            Fraction(x.strip()) for x in self.bodies[n].rsplit("(", 1)[1].split(")")[0].split(",")
        ]

    def frame(self, n):
        return [self.xyz(i) for i in self.refs(self.refs(n)[0])]

    def face_edges(self, n):
        return [
            self.refs(use)[0]
            for bound in self.refs(n)[:-1]
            for use in self.refs(self.refs(bound)[0])
        ]


def normalized(v):
    length = math.sqrt(sum(float(x) ** 2 for x in v))
    assert length > 0
    return [float(x) / length for x in v]


def check_identity(source, emitted, response, shift):
    vertices = emitted.named("VERTEX_POINT")
    edges = emitted.named("EDGE_CURVE")
    faces = emitted.named("ADVANCED_FACE")
    repairs = response["source_admission"]["topology_repairs"]
    source_vertices = {f"#{n}" for n in source.ids("VERTEX_POINT")}
    source_edges = {f"#{n}" for n in source.ids("EDGE_CURVE")}
    derived = {r["derived_seam_edge"] for r in repairs}
    assert set(vertices) == source_vertices, "vertex identity differs"
    assert set(edges) == source_edges | derived, "edge identity differs"
    assert set(faces) == {f"#{n}" for n in source.ids("ADVANCED_FACE")}, "face identity differs"
    step = response["solids"][0]["edited_step" if any(shift) else "source_step"]
    rounding = Fraction(step["max_scalar_rounding_error_mm"])

    def point(a, b):
        assert all(
            abs(1000 * x + s - y) <= rounding for x, y, s in zip(a, b, shift)
        ), "carrier/vertex position differs"

    def direction(a, b):
        assert (
            max(abs(x - y) for x, y in zip(normalized(a), normalized(b))) < 1e-12
        ), "carrier direction differs"

    def frame(a, b, reference=True):
        af, bf = source.frame(a), emitted.frame(b)
        point(af[0], bf[0])
        direction(af[1], bf[1])
        if reference:
            direction(af[2], bf[2])

    for n in source.ids("VERTEX_POINT"):
        point(source.xyz(source.refs(n)[0]), emitted.xyz(emitted.refs(vertices[f"#{n}"])[0]))
    for n in source.ids("EDGE_CURVE"):
        written = edges[f"#{n}"]
        original_refs, written_refs = source.refs(n), emitted.refs(written)
        assert [f"#{v}" for v in original_refs[:2]] == [
            emitted.name(v) for v in written_refs[:2]
        ], "edge endpoint identity differs"
        a, b = original_refs[2], written_refs[2]
        if emitted.kind(b) in ["SURFACE_CURVE", "SEAM_CURVE"]:
            b = emitted.refs(b)[0]
        assert source.kind(a) == emitted.kind(b), "source carrier type differs"
        if source.kind(a) == "LINE":
            ar, br = source.refs(a), emitted.refs(b)
            point(source.xyz(ar[0]), emitted.xyz(br[0]))
            direction(source.xyz(source.refs(ar[1])[0]), emitted.xyz(emitted.refs(br[1])[0]))
        else:
            assert source.kind(a) == "CIRCLE"
            # Complete circles retain the geometric carrier; their reference
            # phase may rotate to the retained shared endpoint for STEP ranges.
            frame(a, b, reference=original_refs[0] != original_refs[1])
            ra = Fraction(source.bodies[a].rsplit(",", 1)[1].strip(" )"))
            rb = Fraction(emitted.bodies[b].rsplit(",", 1)[1].strip(" )"))
            assert abs(1000 * ra - rb) <= rounding, "circle radius differs"
    for n in source.ids("ADVANCED_FACE"):
        written = faces[f"#{n}"]
        a, b = source.refs(n)[-1], emitted.refs(written)[-1]
        assert source.kind(a) == emitted.kind(b)
        frame(a, b)
        if source.kind(a) == "CYLINDRICAL_SURFACE":
            ra = Fraction(source.bodies[a].rsplit(",", 1)[1].strip(" )"))
            rb = Fraction(emitted.bodies[b].rsplit(",", 1)[1].strip(" )"))
            assert abs(1000 * ra - rb) <= rounding
        assert source.bodies[n].replace(" ", "").endswith(".T.)") == emitted.bodies[
            written
        ].replace(" ", "").endswith(".T.)"), "face orientation differs"
        expected = Counter(f"#{edge}" for edge in source.face_edges(n))
        for repair in repairs:
            if repair["face"] == f"#{n}":
                expected[repair["derived_seam_edge"]] += 2
        assert (
            Counter(emitted.name(edge) for edge in emitted.face_edges(written)) == expected
        ), "face-edge incidence differs"
    for repair in repairs:
        assert repair["solid"] == "#73" and repair["method"] == "paired_cylinder_generatrix/1"
        edge = edges[repair["derived_seam_edge"]]
        assert [emitted.name(v) for v in emitted.refs(edge)[:2]] == repair["source_vertices"]
        seam = emitted.refs(edge)[2]
        assert emitted.kind(seam) == "SEAM_CURVE" and len(emitted.refs(seam)) == 3
        face = faces[repair["face"]]
        assert (
            Counter(emitted.name(e) for e in emitted.face_edges(face))[repair["derived_seam_edge"]]
            == 2
        )
        uses = [u for u in emitted.ids("ORIENTED_EDGE") if emitted.refs(u) == [edge]]
        assert len(uses) == 2 and sorted(emitted.bodies[u].endswith(".T.)") for u in uses) == [
            False,
            True,
        ]
    return {
        "faces": len(faces),
        "retained_source_edges": len(source_edges),
        "derived_seams": len(derived),
        "vertices": len(vertices),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--example-binary", type=Path, required=True)
    parser.add_argument("--example-sha256", required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--phased-authored", type=Path, required=True)
    parser.add_argument(
        "--precision-authored",
        nargs=2,
        type=Path,
        required=True,
        metavar=("SHIFTED_CENTER", "OBLIQUE_ANCHOR"),
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--raw-helper",
        type=Path,
        default=Path(__file__).resolve().parents[2]
        / "tools/step-regressions/check_numerical_step.py",
    )
    args = parser.parse_args()
    assert (
        sha(args.raw_helper) == "5633cdcf16efc13f328ed85b33a738f6a9296e2afb925f4c64bbed7c91025376"
    )
    assert OCP.__version__ == "8.0.1.1"
    sys.path.insert(0, str(args.raw_helper.parent))
    from check_numerical_step import raw_solid, state, topology, volume

    assert re.fullmatch("[0-9a-f]{40}", args.source_revision)
    assert re.fullmatch("[0-9a-f]{64}", args.example_sha256)
    assert sha(args.source) == SOURCE_SHA and sha(args.example_binary) == args.example_sha256
    source_text = args.source.read_text()
    assert re.search(r"CONVERSION_BASED_UNIT\(\s*'METRE'", source_text)
    source = Entities(source_text)
    r = json.loads((args.capture / "response.json").read_text())
    assert r["source_sha256"] == SOURCE_SHA and r["source_record"]["source"]["sha256"] == SOURCE_SHA
    assert r["status"] == "numerical" and r["exact"] is False and len(r["solids"]) == 1
    assert r["source_admission"]["tolerance_mm"] == "1/100"
    assert len(r["source_admission"]["topology_repairs"]) == 2
    assert len(r["source_admission"]["intersection_checks"]) == 16 * 15 // 2
    assert r["work_used"] <= r["work_limit"] == 20_000_000
    assert r["work_used"] == r["work"]["import"] + sum(r["solids"][0]["work"].values())
    phase_checks = {}
    for field in ["source_step", "edited_step"]:
        check = r["solids"][0][field]["closed_circle_phase_checks"]
        assert check["count"] == 4 and check["angular_limit_rad"] == "1/1000000000000"
        assert check["distance_limit_mm"] == "1/1000000000"
        assert 0 <= Fraction(check["max_tangent_bound_mm"]) <= Fraction(1, 1_000_000_000)
        phase_checks[field] = check

    artifacts, shapes = [], []
    identity = []
    for filename, shift, field in [
        ("solid-0.step", [0, 0, 0], "source_step"),
        ("solid-0-edited.step", [3, 5, 7], "edited_step"),
    ]:
        path = args.capture / filename
        original_bytes = path.read_bytes()
        assert sha(path) == r["solids"][0][field]["sha256"]
        assert (
            b"catallax numerical approximate" in original_bytes
            and b"complete_body_bound=unestablished" in original_bytes
        )
        assert b"UNCERTAINTY_MEASURE_WITH_UNIT" not in original_bytes
        shape = raw_solid(path)
        assert BRepCheck_Analyzer(shape).IsValid(), (filename, "raw invalid")
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(shape, props)
        assert props.Mass() > 0
        centre = props.CentreOfMass()
        identity.append(check_identity(source, Entities(original_bytes.decode()), r, shift))
        assert path.read_bytes() == original_bytes
        shapes.append(shape)
        artifacts.append(
            dict(
                sha256=sha(path),
                raw_valid=True,
                topology=topology(shape),
                volume_mm3=props.Mass(),
                centroid_mm=[centre.X(), centre.Y(), centre.Z()],
            )
        )
    assert (
        identity[0]
        == identity[1]
        == dict(faces=16, retained_source_edges=40, derived_seams=2, vertices=28)
    )
    assert artifacts[0]["topology"] == artifacts[1]["topology"]
    assert artifacts[0]["topology"]["solids"] == 1
    assert (
        abs(artifacts[0]["volume_mm3"] - artifacts[1]["volume_mm3"])
        < 1e-9 * artifacts[0]["volume_mm3"]
    )
    assert all(
        abs(b - a - delta) < 1e-8
        for a, b, delta in zip(artifacts[0]["centroid_mm"], artifacts[1]["centroid_mm"], [3, 5, 7])
    )

    # This is a separately labelled normal source transfer, with OCCT defaults.
    reader = STEPControl_Reader()
    assert Interface_Static.CVal_s("xstep.cascade.unit") == "MM"
    assert reader.ReadFile(str(args.source)) == IFSelect_RetDone
    roots = reader.NbRootsForTransfer()
    transferred = reader.TransferRoots()
    reference = reader.OneShape()
    assert transferred == roots == 1 and BRepCheck_Analyzer(reference).IsValid()
    assert topology(reference)["solids"] == 1
    assert abs(volume(reference) - artifacts[0]["volume_mm3"]) < 1e-8 * volume(reference)
    probes = []
    for repair in r["source_admission"]["topology_repairs"]:
        circles = [source.refs(int(edge[1:]))[2] for edge in repair["source_rim_edges"]]
        centres = [source.frame(c)[0] for c in circles]
        centre = [500 * float(a + b) for a, b in zip(*centres)]
        radial = normalized(source.frame(circles[0])[2])
        radius = 1000 * float(Fraction(source.bodies[circles[0]].rsplit(",", 1)[1].strip(" )")))
        assert state(reference, centre) == "outside"
        stock = next(
            (
                [x + factor * radius * d for x, d in zip(centre, radial)]
                for factor in [1.25, -1.25, 2, -2, 3, -3]
                if state(reference, [x + factor * radius * d for x, d in zip(centre, radial)])
                == "inside"
            ),
            None,
        )
        assert stock is not None
        for point, expected in [(centre, "outside"), (stock, "inside")]:
            assert (
                state(shapes[0], point) == expected
                and state(shapes[1], [p + d for p, d in zip(point, [3, 5, 7])]) == expected
            )
            probes.append(dict(face=repair["face"], expected=expected, before_and_pose_match=True))
    # Negative identity control: geometrically identical STEP with exchanged names.
    text = (args.capture / "solid-0.step").read_text()
    names = list(Entities(text).named("VERTEX_POINT"))[:2]
    changed = (
        text.replace("'" + names[0] + "'", "'identity-swap'")
        .replace("'" + names[1] + "'", "'" + names[0] + "'")
        .replace("'identity-swap'", "'" + names[1] + "'")
    )
    try:
        check_identity(source, Entities(changed), r, [0, 0, 0])
    except AssertionError:
        pass
    else:
        raise AssertionError("exchanged source-local vertex identities were accepted")
    phased = raw_solid(args.phased_authored)
    assert BRepCheck_Analyzer(phased).IsValid()
    assert abs(volume(phased) - (16 - math.pi)) < 1e-9
    assert topology(phased) == dict(
        solids=1, faces=7, edges=15, vertices=10, two_loop_faces=2, euler_with_face_holes=0
    )
    assert state(phased, (0, 0, 0.5)) == "outside" and state(phased, (1.5, 0, 0.5)) == "inside"
    precision_authored = []
    for path, shift, role in zip(
        args.precision_authored, [-1.0004, 0], ["shifted_center", "oblique_anchor"]
    ):
        shape = raw_solid(path)
        assert BRepCheck_Analyzer(shape).IsValid(), role
        assert abs(volume(shape) - (16 - math.pi)) < 1e-9
        assert topology(shape) == topology(phased)
        assert (
            state(shape, (0, shift, 0.5)) == "outside"
            and state(shape, (1.5, shift, 0.5)) == "inside"
        )
        precision_authored.append(
            dict(
                role=role,
                sha256=sha(path),
                raw_valid=True,
                volume_and_material_checked=True,
                significant_digits=17,
            )
        )
    # Keep the actual shape properties private. Public artifacts retain checked
    # relationships, source-local identities and hashes, rather than coordinates.
    public_artifacts = [
        {k: v for k, v in artifact.items() if k not in ["volume_mm3", "centroid_mm"]}
        for artifact in artifacts
    ]
    report = dict(
        format="catallax-multiseam-qualification/1",
        model_id="00003390",
        source_sha256=SOURCE_SHA,
        source_revision=args.source_revision,
        revision_binary_association="author_attested_frozen_offload_artifact",
        example_binary_sha256=args.example_sha256,
        checker_sha256=sha(Path(__file__)),
        raw_builder_sha256=sha(args.raw_helper),
        response_sha256=sha(args.capture / "response.json"),
        ocp_version=OCP.__version__,
        tolerance_mm="1/100",
        work_limit=r["work_limit"],
        work_used=r["work_used"],
        import_work=r["work"]["import"],
        stages=r["solids"][0]["work"],
        topology_repairs=r["source_admission"]["topology_repairs"],
        identity=identity[0],
        closed_circle_identity="retained_center_axis_radius_and_endpoint; reference_phase_may_rotate",
        source_identity_swap_rejected=True,
        closed_circle_phase_checks=phase_checks,
        source_unit_to_mm=1000,
        exported_step_healing=False,
        exported_step_same_parameter=False,
        source_reference=dict(
            method="normal_STEPControl_Reader_TransferRoots_defaults",
            healing_may_run=True,
            roots=roots,
            transferred=transferred,
            valid=True,
            exported_volume_matches_within_relative="1e-8",
        ),
        source_step=public_artifacts[0],
        pose_step=public_artifacts[1],
        pose_volume_matches_within_relative="1e-9",
        centroid_shift_3_5_7_mm_checked=True,
        phased_authored=dict(
            sha256=sha(args.phased_authored), raw_valid=True, volume_and_material_checked=True
        ),
        precision_authored=precision_authored,
        material_probes=probes,
        complete_body_bound="unestablished",
        catallax_boolean_qualified=False,
    )
    args.out.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps(report, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
