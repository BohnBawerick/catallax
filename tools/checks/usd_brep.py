#!/usr/bin/env python3
"""Decode the golden export with OpenUSD and check analytic obligations.

Requires usd-core. This does not load the NVIDIA consumer or qualify its radial
convention. Run DIRECTORY [--edit DIRECTORY] [--negative-controls].
"""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from pxr import Gf, Sdf, Usd, UsdGeom, Vt


SCHEMA_COMMIT = "c9c5979f04cea8c38213004d4f04df0dcfb40d05"
APIS = [
    "BrepSurfacePlaneAPI",
    "BrepPointAPI:vertexPoint",
    "BrepCurve3dLineAPI:edge3dLine",
    "BrepCurveUvNurbAPI",
]
POINTS = "brep:vertexPoint:point:position"
ORIGINS = "brep:surface:plane:origin"
NORMALS = "brep:surface:plane:axis"
REFS = "brep:surface:plane:refDirection"
LINE = "brep:edge3dLine:curve3d:line:"
UV = "brep:curveUv:nurb:"


def need(ok, message):
    if not ok:
        raise ValueError(message)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(a, s):
    return tuple(x * s for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vec(v):
    return tuple(Q(float(x)) for x in v)


def cyclic(a, b):
    return len(a) == len(b) and any(a[i:] + a[:i] == b for i in range(len(a)))


def fixture(height):
    bottom = [
        (-5, -4, 0),
        (5, -4, 0),
        (5, 4, 0),
        (-5, 4, 0),
        (-1, -1, 0),
        (1, -1, 0),
        (1, 1, 0),
        (-1, 1, 0),
    ]
    points = bottom + [(x, y, height) for x, y, _ in bottom]
    loops = [[[3, 2, 1, 0], [4, 5, 6, 7]], [[8, 9, 10, 11], [12, 15, 14, 13]]]
    loops += [[[i, (i + 1) % 4, (i + 1) % 4 + 8, i + 8]] for i in range(4)]
    loops += [[[i, i + 8, 4 + (i - 3) % 4 + 8, 4 + (i - 3) % 4]] for i in range(4, 8)]
    return points, loops


def attributes(prim, expected):
    authored = {a.GetName(): a for a in prim.GetAuthoredAttributes()}
    need(
        set(authored) == set(expected),
        f"unexpected/missing attributes: {set(authored) ^ set(expected)}",
    )
    for name, ty in expected.items():
        need(str(authored[name].GetTypeName()) == ty, f"wrong declared type for {name}")
        need(not authored[name].GetTimeSamples(), f"animated attribute {name}")
    return {name: attr.Get() for name, attr in authored.items()}


def properties(points, polygons):
    volume = Q(0)
    area = Q(0)
    moment = (Q(0),) * 3
    for polygon in polygons:
        a = points[polygon[0]]
        # Signed fans retain inner-loop subtraction without filling the hole.
        normal = cross(sub(points[polygon[1]], a), sub(points[polygon[2]], a))
        length = sum(abs(x) for x in normal)
        need(length > 0 and sum(x != 0 for x in normal) == 1, "non-axis planar polygon")
        normal = mul(normal, Q(1, length))
        for i in range(1, len(polygon) - 1):
            b, c = points[polygon[i]], points[polygon[i + 1]]
            dv = Q(dot(a, cross(b, c)), 6)
            volume += dv
            moment = add(moment, mul(add(add(a, b), c), dv / 4))
            area += dot(cross(sub(b, a), sub(c, a)), normal) / 2
    return volume, area, mul(moment, 1 / volume)


def check_geometry(stage, height, placed):
    expected_points, expected_loops = fixture(height)
    need(stage.HasAuthoredMetadata("metersPerUnit"), "absent units")
    need(stage.GetMetadata("metersPerUnit") == 0.001, "wrong units")
    need(stage.GetMetadata("upAxis") == "Z", "wrong up axis")
    need(str(stage.GetDefaultPrim().GetPath()) == "/Part", "wrong default prim")
    need(
        [(str(p.GetPath()), p.GetTypeName()) for p in stage.Traverse()]
        == [
            ("/Part", "Xform"),
            ("/Part/Solid", "BrepArray"),
            ("/Part/Display", "Mesh"),
        ],
        "unexpected scene content",
    )
    part = stage.GetPrimAtPath("/Part")
    pa = attributes(part, {"xformOp:transform": "matrix4d", "xformOpOrder": "token[]"})
    need(list(pa["xformOpOrder"]) == ["xformOp:transform"], "transform op order")
    m = pa["xformOp:transform"]
    transform = lambda p: (20 - p[1], 30 + p[0], 40 + p[2]) if placed else p
    for p in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]:
        need(
            tuple(m.Transform(Gf.Vec3d(*p))) == transform(p),
            "transform origin/basis action",
        )
    solid = stage.GetPrimAtPath("/Part/Solid")
    need(solid.GetMetadata("apiSchemas").GetAppliedItems() == APIS, "wrong APIs")
    need(solid.GetCustomData() == {"source": "Catallax"}, "producer customData")
    types = {
        name: "uint[]"
        for name in [
            "brep:regionCount",
            "region:shellCount",
            "shell:faceuseCount",
            "shell:wireEdgeCount",
            "faceuse:faceIndex",
            "face:loopCount",
            "loop:edgeuseCount",
            "loop:vertexIndex",
            "edgeuse:edgeIndex",
            "edgeuse:nextRadialEUIndex",
            UV + "vertexCount",
            UV + "order",
        ]
    }
    types.update(
        {
            name: "token[]"
            for name in [
                "region:type",
                "shell:pointType",
                "faceuse:orientationType",
                "face:surfaceType",
                "face:trimType",
                "edgeuse:orientationType",
                "edgeuse:thisRadialEntryType",
                "edge:curveType",
                "vertex:pointType",
            ]
        }
    )
    types.update(
        {
            "purpose": "token",
            "brep:intersectTol3d": "double[]",
            "brep:extent": "double3[]",
            "face:range": "double2[]",
            "edge:range": "double[]",
            "edge:vertexIndices": "int2[]",
            ORIGINS: "point3d[]",
            NORMALS: "vector3d[]",
            REFS: "vector3d[]",
            POINTS: "point3d[]",
            LINE + "origin": "point3d[]",
            LINE + "direction": "vector3d[]",
            UV + "controlVertices": "double2[]",
            UV + "knots": "double[]",
            UV + "weights": "double[]",
        }
    )
    a = attributes(solid, types)
    need(
        all(
            solid.GetAttribute(n).GetVariability() == Sdf.VariabilityUniform
            for n in types
        ),
        "nonuniform BRep data",
    )
    need(a["purpose"] == "guide", "BRep must not duplicate the display")
    expected = {
        "brep:regionCount": [2],
        "region:type": ["voidRegion", "solidRegion"],
        "region:shellCount": [1, 1],
        "shell:faceuseCount": [10, 10],
        "shell:wireEdgeCount": [0, 0],
        "shell:pointType": ["none", "none"],
        "faceuse:faceIndex": list(range(10)) * 2,
        "faceuse:orientationType": ["same"] * 10 + ["opposite"] * 10,
        "face:loopCount": [2, 2] + [1] * 8,
        "loop:edgeuseCount": [4] * 12,
        "loop:vertexIndex": [0] * 12,
        "face:surfaceType": ["BrepSurfacePlaneAPI"] * 10,
        "face:trimType": ["rectangular"] * 10,
        "edge:curveType": ["BrepCurve3dLineAPI"] * 24,
        "vertex:pointType": ["BrepPointAPI"] * 16,
        UV + "vertexCount": [2] * 48,
        UV + "order": [2] * 48,
        UV + "weights": [1] * 96,
        "brep:intersectTol3d": [0.000001],
    }
    for name, value in expected.items():
        need(list(a[name]) == value, f"invalid {name}")
    points = list(map(vec, a[POINTS]))
    need(
        points == expected_points and len(set(points)) == 16,
        "vertex coordinates/duplicates",
    )
    need(
        list(map(vec, a["brep:extent"])) == [(-5, -4, 0), (5, 4, height)],
        "local extent",
    )
    edges = [tuple(e) for e in a["edge:vertexIndices"]]
    need(
        len(edges) == 24 and len(set(edges)) == 24 and edges == sorted(edges),
        "duplicate/unsorted edges",
    )
    need(all(0 <= x < y < 16 for x, y in edges), "edge vertex index")
    need(
        len(a["edge:range"]) == 48
        and len(a[LINE + "origin"]) == len(a[LINE + "direction"]) == 24,
        "line counts",
    )
    for e, (x, y) in enumerate(edges):
        lo, hi = a["edge:range"][(2 * e) : (2 * e + 2)]
        origin, direction = vec(a[LINE + "origin"][e]), vec(a[LINE + "direction"][e])
        need(
            lo == 0 and hi > 0 and dot(direction, direction) == 1,
            "line range/unit direction",
        )
        need(
            origin == points[x] and add(origin, mul(direction, Q(hi))) == points[y],
            "line endpoint incidence",
        )
    need(
        all(
            len(a[name]) == 48
            for name in [
                "edgeuse:edgeIndex",
                "edgeuse:orientationType",
                "edgeuse:nextRadialEUIndex",
                "edgeuse:thisRadialEntryType",
            ]
        ),
        "edgeuse counts",
    )
    need(
        len(a[UV + "controlVertices"]) == 96 and len(a[UV + "knots"]) == 192,
        "missing UV records",
    )
    need(
        len(a[ORIGINS]) == len(a[NORMALS]) == len(a[REFS]) == 10
        and len(a["face:range"]) == 20,
        "plane counts",
    )
    users = defaultdict(list)
    paths = []
    decoded_faces = []
    cursor = 0
    for f, nloops in enumerate(a["face:loopCount"]):
        origin, n, u = vec(a[ORIGINS][f]), vec(a[NORMALS][f]), vec(a[REFS][f])
        need(dot(n, n) == dot(u, u) == 1 and dot(n, u) == 0, "plane frame")
        v = cross(n, u)
        face_paths = []
        projected = []
        for _ in range(nloops):
            path, ends = [], []
            for _ in range(a["loop:edgeuseCount"][len(paths)]):
                k = cursor
                e = a["edgeuse:edgeIndex"][k]
                need(0 <= e < 24, "edgeuse index")
                sense = a["edgeuse:orientationType"][k]
                need(sense in ["same", "opposite"], "edgeuse sense")
                x, y = edges[e] if sense == "same" else edges[e][::-1]
                path.append(x)
                ends.append(y)
                users[e].append((k, f, x, y))
                controls = list(map(vec, a[UV + "controlVertices"][2 * k : 2 * k + 2]))
                for p, uv in zip([points[x], points[y]], controls):
                    need(
                        add(origin, add(mul(u, uv[0]), mul(v, uv[1]))) == p,
                        "UV lifting/plane incidence",
                    )
                projected.extend(controls)
                length = a["edge:range"][2 * e + 1]
                need(
                    list(a[UV + "knots"][4 * k : 4 * k + 4]) == [0, 0, length, length],
                    "UV parameter range",
                )
                # Material lies left of traversal. Derive the entry side from
                # the right-hand angular derivative around the directed edge.
                tangent = vec(a[LINE + "direction"][e])
                inward = cross(n, sub(points[y], points[x]))
                entry = (
                    "bottomEntry" if dot(cross(tangent, inward), n) > 0 else "topEntry"
                )
                need(
                    a["edgeuse:thisRadialEntryType"][k] == entry,
                    "radial entry geometry",
                )
                cursor += 1
            need(ends == path[1:] + path[:1], "loop not closed")
            face_paths.append(path)
            paths.append(path)
        need(
            all(
                cyclic(path, expected_loops[f][j]) for j, path in enumerate(face_paths)
            ),
            "wrong winding/boundary or overlapping patch",
        )
        outer = face_paths[0]
        normal = cross(
            sub(points[outer[1]], points[outer[0]]),
            sub(points[outer[2]], points[outer[1]]),
        )
        need(dot(normal, n) > 0, "outward chart normal")
        need(origin == points[outer[0]], "chart origin")
        need(
            list(map(vec, a["face:range"][2 * f : 2 * f + 2]))
            == [
                tuple(min(p[i] for p in projected) for i in range(2)),
                tuple(max(p[i] for p in projected) for i in range(2)),
            ],
            "face UV range",
        )
        decoded_faces.append(face_paths)
    need(cursor == 48 and len(paths) == 12, "prefix sums")
    for e in range(24):
        uses = users[e]
        need(len(uses) == 2, "edge is not shared twice")
        (k, f, x, y), (j, g, z, w) = uses
        need((x, y) == (w, z) and f != g, "edge traversal/face incidence")
        need(
            a["edgeuse:nextRadialEUIndex"][k] == j
            and a["edgeuse:nextRadialEUIndex"][j] == k,
            "broken radial cycle",
        )
    # Each vertex link is one connected cycle of incident faces.
    for vertex in range(16):
        adjacency = defaultdict(list)
        for e, (x, y) in enumerate(edges):
            if vertex in (x, y):
                f, g = [use[1] for use in users[e]]
                adjacency[f].append(g)
                adjacency[g].append(f)
        need(
            all(len(neighbors) == 2 for neighbors in adjacency.values()),
            "noncircular vertex link",
        )
        visited, todo = set(), [next(iter(adjacency))]
        while todo:
            f = todo.pop()
            if f not in visited:
                visited.add(f)
                todo.extend(adjacency[f])
        need(len(visited) == len(adjacency), "disconnected vertex link")
    need(
        16 - 24 + sum(2 - n for n in a["face:loopCount"]) == 0,
        "genus-one Euler characteristic",
    )
    # Integrate signed loop fans, including the negative inner-loop areas.
    volume, _, centroid = properties(points, paths)
    area = Q(0)
    for f, face in enumerate(expected_loops):
        n = vec(a[NORMALS][f])
        for path in face:
            area += (
                sum(
                    dot(cross(points[path[i]], points[path[(i + 1) % 4]]), n)
                    for i in range(4)
                )
                / 2
            )
    need(
        (volume, area, centroid)
        == (76 * height, 152 + 44 * height, (0, 0, Q(height, 2))),
        "analytic mass properties",
    )

    # The exact boundary comparison above proves nonoverlap and no extra voids.
    # These witnesses additionally exercise material/through-hole classification.
    def in_loop(p, path):
        crossings = 0
        for i, index in enumerate(path):
            x, y, _ = points[index]
            xx, yy, _ = points[path[(i + 1) % len(path)]]
            if (y > p[1]) != (yy > p[1]):
                crossings += p[0] < x + (p[1] - y) * (xx - x) / (yy - y)
        return crossings % 2 == 1

    def inside(p):
        winding = 0
        for f, face in enumerate(decoded_faces):
            normal = vec(a[NORMALS][f])
            if normal[2] and points[face[0][0]][2] > p[2]:
                if in_loop(p, face[0]) and not any(
                    in_loop(p, hole) for hole in face[1:]
                ):
                    winding += normal[2]
        return winding == 1

    need(
        inside((3, 0, 1))
        and not any(inside(p) for p in [(0, 0, 1), (0, 0, -1), (0, 0, height + 1)]),
        "material/void witnesses",
    )
    mesh = stage.GetPrimAtPath("/Part/Display")
    ma = attributes(
        mesh,
        {
            "subdivisionScheme": "token",
            "orientation": "token",
            "points": "point3f[]",
            "faceVertexCounts": "int[]",
            "faceVertexIndices": "int[]",
        },
    )
    need(
        ma["subdivisionScheme"] == "none" and ma["orientation"] == "rightHanded",
        "mesh orientation",
    )
    need(
        list(map(vec, ma["points"])) == points
        and list(ma["faceVertexCounts"]) == [4] * 16,
        "display mesh points/counts",
    )
    indices = list(ma["faceVertexIndices"])
    need(len(indices) == 64 and all(0 <= i < 16 for i in indices), "mesh indices")
    quads = [indices[i : i + 4] for i in range(0, 64, 4)]
    expected_quads = []
    for i in range(4):
        j = (i + 1) % 4
        expected_quads += [[j, i, i + 4, j + 4], [i + 8, j + 8, j + 12, i + 12]]
    expected_quads += [f[0] for f in expected_loops[2:]]
    need(quads == expected_quads, "display filled hole/overlap")
    mv, ma_area, mc = properties(points, quads)
    need((mv, ma_area, mc) == (volume, area, centroid), "display geometry differs")
    mesh_edges = Counter(
        tuple(sorted((q[i], q[(i + 1) % 4]))) for q in quads for i in range(4)
    )
    need(set(mesh_edges.values()) == {2}, "display cracks")
    world = [
        vec(
            UsdGeom.XformCache()
            .GetLocalToWorldTransform(mesh)
            .Transform(Gf.Vec3d(*map(float, p)))
        )
        for p in points
    ]
    need(world == [transform(p) for p in points], "transform applied more than once")
    bounds = [tuple(fn(p[i] for p in world) for i in range(3)) for fn in [min, max]]
    need(
        bounds
        == (
            [(16, 25, 40), (24, 35, 40 + height)]
            if placed
            else [(-5, -4, 0), (5, 4, height)]
        ),
        "world bounds",
    )
    wv, wa, wc = properties(world, quads)
    need((wv, wa, wc) == (volume, area, transform(centroid)), "world mass/orientation")
    return (
        {
            "volume_mm3": str(volume),
            "area_mm2": str(area),
            "centroid_mm": list(map(str, centroid)),
            "world_bounds_mm": [list(map(str, p)) for p in bounds],
        },
        users,
        edges,
    )


def identity_check(metadata, users, edges):
    origin = lambda role: {"step": "solid", "role": role}
    faces = [[origin("start")], [origin("end")]] + [
        [origin({"side": {"segment": prefix + str(i)}})]
        for prefix in ["outer", "inner"]
        for i in range(4)
    ]
    key = lambda o: json.dumps(o, separators=(",", ":"))
    names = {
        "faces": faces,
        "edges": [
            {"between": sorted([faces[u[1]][0] for u in users[e]], key=key)}
            for e in range(24)
        ],
        "vertices": [
            {
                "faces": sorted(
                    [
                        faces[f][0]
                        for f in sorted(
                            {
                                u[1]
                                for e, vs in enumerate(edges)
                                if v in vs
                                for u in users[e]
                            }
                        )
                    ],
                    key=key,
                )
            }
            for v in range(16)
        ],
    }
    need(metadata["identities"] == names, "altered persistent identity mapping")


def check_package(directory):
    metadata = json.loads((directory / "metadata.json").read_bytes())
    need(
        metadata["format"] == "catallax-usd-golden/1"
        and metadata["fixture"] == "analytic-square-through-hole/1",
        "metadata format",
    )
    height, placement = metadata["height_mm"], metadata["placement"]
    need(height in [2, 3] and placement in ["identity", "placed"], "fixture parameters")
    for file, expected in metadata["files"].items():
        need(
            file == "part.usda"
            and hashlib.sha256((directory / file).read_bytes()).hexdigest() == expected,
            "artifact digest",
        )
    need(set(metadata["files"]) == {"part.usda"}, "artifact inventory")
    recipe_bytes = (directory / "recipe.json").read_bytes()
    need(
        metadata["recipe"]["file"] == "recipe.json"
        and hashlib.sha256(recipe_bytes).hexdigest() == metadata["recipe"]["sha256"],
        "recipe digest",
    )
    recipe = json.loads(recipe_bytes)
    need(
        recipe["format"] == "catallax-recipe/3"
        and recipe["steps"][1]
        == {
            "id": "solid",
            "op": "extrude",
            "profile": "profile",
            "distance": str(height),
        },
        "recipe revision",
    )
    sketch = recipe["steps"][0]
    need(
        sketch["id"] == "profile"
        and sketch["op"] == "sketch"
        and sketch["plane"]
        == {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]},
        "recipe chart",
    )
    expected_loops = []
    for prefix, points in [
        ("outer", [(5, -4), (5, 4), (-5, 4), (-5, -4)]),
        ("inner", [(1, -1), (1, 1), (-1, 1), (-1, -1)]),
    ]:
        expected_loops.append(
            [
                {"kind": "line", "id": prefix + str(i), "to": list(map(str, p))}
                for i, p in enumerate(points)
            ]
        )
    need(
        sketch["loops"] == expected_loops
        and recipe["result"] == "solid"
        and len(recipe["steps"]) == 2,
        "recipe geometry/identity",
    )
    need(
        metadata["source"]["schema_commit"] == SCHEMA_COMMIT
        and metadata["source"]["commit_authority"] == "caller-supplied",
        "source provenance",
    )
    commit = metadata["source"]["commit"]
    need(
        len(commit) == 40 and all(c in "0123456789abcdefABCDEF" for c in commit),
        "source commit",
    )
    need(
        metadata["units"]["nominal_meters_per_unit"] == "1/1000"
        and metadata["units"]["encoded_meters_per_unit"] == 0.001,
        "nominal units",
    )
    need(
        metadata["intersect_tolerance"]["nominal_stage_units"] == "1/1000000"
        and metadata["intersect_tolerance"]["encoded_stage_units"] == 0.000001,
        "tolerance policy",
    )
    need(
        metadata["proof"]
        == {
            "source_topology": "checked",
            "source_geometry": "exact analytic construction",
            "acceptance": "requires independent checker",
            "consumer_interoperability": "unproved",
            "healing": "none",
        },
        "proof downgrade",
    )
    need(
        len(metadata["charts"]) == 10
        and all(
            c["export_chart"] == "outer first vertex, U first segment, V=N cross U"
            for c in metadata["charts"]
        ),
        "chart provenance",
    )
    for kind, count in [("faces", 10), ("edges", 24), ("vertices", 16), ("uv", 48)]:
        rows = metadata["local_bounds"][kind]
        need(
            len(rows) == count
            and all(
                row
                == {
                    "index": i,
                    "approximation_mm": "0",
                    "encoding_mm": "0",
                    "basis": "exact integer endpoints and coordinate-axis carriers",
                }
                for i, row in enumerate(rows)
            ),
            "per-entity bounds",
        )
    need(
        set(metadata["preservation"])
        == {
            "geometry",
            "topology",
            "trims",
            "names",
            "recipe",
            "local_certificates",
            "units",
            "proof_provenance",
        }
        and len(metadata["information_loss"]) == 2,
        "loss accounting",
    )
    stage = Usd.Stage.Open(str(directory / "part.usda"))
    need(stage is not None, "OpenUSD decode failed")
    result, users, edges = check_geometry(stage, height, placement == "placed")
    identity_check(metadata, users, edges)
    points, faces = fixture(height)
    for chart, face in zip(metadata["charts"], faces):
        origin = tuple(chart["source_origin"])
        normal = tuple(chart["source_normal"])
        reference = tuple(chart["source_ref_direction"])
        need(chart["source_sense"] in ["Same", "Reversed"], "source chart sense")
        need(
            dot(normal, normal) > 0
            and dot(reference, reference) > 0
            and dot(normal, reference) == 0,
            "source chart frame",
        )
        need(
            all(
                dot(sub(points[i], origin), normal) == 0 for loop in face for i in loop
            ),
            "source chart incidence",
        )
        outer = face[0]
        outward = cross(
            sub(points[outer[1]], points[outer[0]]),
            sub(points[outer[2]], points[outer[1]]),
        )
        sign = 1 if chart["source_sense"] == "Same" else -1
        need(dot(outward, normal) * sign > 0, "source chart orientation")
    unit_rounding = Q(stage.GetMetadata("metersPerUnit")) - Q(1, 1000)
    scale_bound = abs(unit_rounding) * (40 + height if placement == "placed" else 5)
    need(
        Q(metadata["units"]["rounding_m_per_unit"]) == unit_rounding
        and Q(metadata["units"]["max_world_coordinate_error_m"]) == scale_bound,
        "physical unit-scale error accounting",
    )
    result["unit_ratio_rounding_m_per_unit"] = str(unit_rounding)
    result["max_world_coordinate_error_m"] = str(scale_bound)
    tol = (
        stage.GetPrimAtPath("/Part/Solid").GetAttribute("brep:intersectTol3d").Get()[0]
    )
    result["tolerance_rounding_stage_units"] = str(Q(tol) - Q(1, 1000000))
    need(
        Q(metadata["intersect_tolerance"]["rounding_stage_units"])
        == Q(tol) - Q(1, 1000000),
        "interchange tolerance rounding",
    )
    result["consumer_interoperability"] = "unproved"
    return stage, metadata, result, users, edges


def negative_controls(stage, metadata, users, edges):
    rejected = []

    def reject(name, mutate):
        copy = Usd.Stage.CreateInMemory()
        need(
            copy.GetRootLayer().ImportFromString(stage.GetRootLayer().ExportToString()),
            "copy stage",
        )
        mutate(copy)
        try:
            check_geometry(
                copy, metadata["height_mm"], metadata["placement"] == "placed"
            )
        except ValueError:
            rejected.append(name)
        else:
            raise ValueError(f"negative control accepted: {name}")

    def edit(name, fn):
        def mutate(stage):
            attr = stage.GetPrimAtPath("/Part/Solid").GetAttribute(name)
            attr.Set(fn(attr.Get()))

        return mutate

    def reverse_inner(stage):
        prim = stage.GetPrimAtPath("/Part/Solid")
        remap = lambda k: 11 - k if 4 <= k < 8 else k
        for name in [
            "edgeuse:edgeIndex",
            "edgeuse:orientationType",
            "edgeuse:thisRadialEntryType",
            "edgeuse:nextRadialEUIndex",
        ]:
            attr = prim.GetAttribute(name)
            old = attr.Get()
            values = type(old)(old)
            for k in range(48):
                value = old[remap(k)]
                if name == "edgeuse:nextRadialEUIndex":
                    value = remap(value)
                elif 4 <= k < 8 and name == "edgeuse:orientationType":
                    value = "opposite" if value == "same" else "same"
                elif 4 <= k < 8 and name == "edgeuse:thisRadialEntryType":
                    value = "bottomEntry" if value == "topEntry" else "topEntry"
                values[k] = value
            attr.Set(values)
        attr = prim.GetAttribute(UV + "controlVertices")
        old = attr.Get()
        values = Vt.Vec2dArray(old)
        for k in range(4, 8):
            values[2 * k], values[2 * k + 1] = old[2 * remap(k) + 1], old[2 * remap(k)]
        attr.Set(values)

    reject("reversed inner loop", reverse_inner)
    reject(
        "duplicate edge",
        edit("edge:vertexIndices", lambda a: Vt.Vec2iArray([a[0], a[0]] + list(a)[2:])),
    )
    reject(
        "omitted UV record",
        edit(UV + "controlVertices", lambda a: Vt.Vec2dArray(list(a)[:-2])),
    )
    reject(
        "swapped faceuse side",
        edit(
            "faceuse:orientationType",
            lambda a: Vt.TokenArray(["opposite"] + list(a)[1:]),
        ),
    )
    reject(
        "broken radial link",
        edit("edgeuse:nextRadialEUIndex", lambda a: Vt.UIntArray([0] + list(a)[1:])),
    )
    reject("absent units", lambda s: s.ClearMetadata("metersPerUnit"))

    def twice(s):
        attr = s.GetPrimAtPath("/Part").GetAttribute("xformOp:transform")
        m = Gf.Matrix4d(0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 20, 30, 40, 1)
        attr.Set(m * m)

    reject("twice-applied transform", twice)
    changed = json.loads(json.dumps(metadata))
    changed["identities"]["faces"][0], changed["identities"]["faces"][1] = (
        changed["identities"]["faces"][1],
        changed["identities"]["faces"][0],
    )
    try:
        identity_check(changed, users, edges)
    except ValueError:
        rejected.append("altered identity map")
    else:
        raise ValueError("negative control accepted: altered identity map")
    return rejected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--edit", type=Path)
    parser.add_argument("--negative-controls", action="store_true")
    args = parser.parse_args()
    stage, metadata, result, users, edges = check_package(args.directory)
    if args.edit:
        _, edited, edit_result, _, _ = check_package(args.edit)
        need(
            metadata["height_mm"] == 2 and edited["height_mm"] == 3,
            "height edit must be 2 to 3",
        )
        need(
            metadata["identities"] == edited["identities"],
            "identities did not survive edit",
        )
        need(
            metadata["recipe"]["sha256"] != edited["recipe"]["sha256"],
            "recipe edit missing",
        )
        result["height_edit"] = {"identities": "preserved", "geometry": edit_result}
    if args.negative_controls:
        result["rejected_controls"] = negative_controls(stage, metadata, users, edges)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
