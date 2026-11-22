#!/usr/bin/env python3
"""Check a Catallax binary STL with libigl and CGAL, without mesh repair."""

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import struct
import sys

import igl
from igl.copyleft import cgal
import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_mesh(data, manifest):
    require(len(data) >= 84 and data[:8] == b"CXPRINT1", "missing print header")
    route, unit, encoding, field_cap, num, den, cells, cap, *bounds = struct.unpack_from(
        "<BBBBqQQQ4d", data, 8
    )
    count, = struct.unpack_from("<I", data, 80)
    require(count > 0 and len(data) == 84 + 50 * count, "invalid binary STL length")
    require(den > 0 and num > 0, "invalid accuracy")
    accuracy = Fraction(num, den)
    require(
        all(np.isfinite(value) and value >= 0 for value in bounds),
        "invalid error bound",
    )
    model, tessellation, rounding, total = map(Fraction, bounds)
    require(model + tessellation + rounding <= total <= accuracy, "error budget exceeded")

    require(manifest["status"] == "answered", "manifest did not approve STL")
    record = manifest["printing"]
    require((route, unit, encoding, field_cap) == (
        {"field": 2, "brep": 1}[record["route"]], 1, 1, int(record["cell_cap_applicable"])
    ), "route or encoding disagrees with the manifest")
    require(data[76:80] == bytes(4), "unknown print header flags")
    require((cells, cap) == (record["max_cells"], record["max_vertices"]), "print caps disagree")
    require(accuracy == Fraction(record["tolerance_mm"]), "requested accuracy disagrees")
    names = ("model", "tessellation", "encoding_rounding", "deviation")
    for name, bound in zip(names, (model, tessellation, rounding, total)):
        require(Fraction(record[f"{name}_bound_mm"]) == bound, f"{name} bound disagrees")
    digest = sha256(data).hexdigest()
    require(digest == record["sha256"] == manifest["files"]["print.stl"]["sha256"], "STL hash disagrees")
    require(count == record["triangles"], "facet count disagrees")
    require(not field_cap or 0 < record["cells"] <= cells, "field cell cap exceeded")

    layout = np.dtype([("normal", "<f4", (3,)), ("vertices", "<f4", (3, 3)), ("flags", "<u2")])
    facets = np.frombuffer(data, dtype=layout, count=count, offset=84)
    require(np.isfinite(facets["vertices"]).all(), "nonfinite vertex")
    require(np.isfinite(facets["normal"]).all(), "nonfinite normal")
    points, indices = np.unique(facets["vertices"].reshape(-1, 3), axis=0, return_inverse=True)
    require(len(points) == record["vertices"] <= cap, "vertex count disagrees or exceeds cap")
    points = points.astype(np.float64)
    faces = indices.reshape(-1, 3).astype(np.int64)
    triangles = points[faces]
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    require(np.all(np.einsum("ij,ij->i", cross, facets["normal"]) > 0), "degenerate or reversed facet")
    return points, faces, digest, accuracy, total


def check_shells(points, faces):
    require(igl.is_edge_manifold(faces)[0], "nonmanifold or inconsistently oriented edge")
    require(igl.is_vertex_manifold(faces).all(), "nonmanifold vertex")
    require(igl.boundary_facets(faces)[0].size == 0, "open boundary")
    _, _, pairs, _, _ = cgal.remesh_self_intersections(points, faces, detect_only=True, first_only=True)
    require(pairs.size == 0, f"intersecting facets {pairs.tolist()}")

    count, component = igl.facet_components(faces)
    shells = [faces[component == i] for i in range(count)]
    # A vertex of one disjoint shell can be tested against every other shell.
    probes = points[[shell[0, 0] for shell in shells]]
    depth = np.zeros(count, dtype=np.int64)
    details = []
    for i, shell in enumerate(shells):
        other = np.arange(count) != i
        windings = igl.winding_number(points, shell, probes[other]) if count > 1 else []
        depth[other] += np.abs(windings) > 0.5
        triangles = points[shell] - probes[i]
        signed_volume = float(np.einsum(
            "ij,ij->i", triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])
        ).sum() / 6)
        # Each shell is closed and triangulated, so E = 3F/2.
        euler = int(np.unique(shell).size) - len(shell) // 2
        details.append({"euler": euler, "triangles": len(shell), "signed_volume_mm3": signed_volume})
    for shell, nesting in zip(details, depth.tolist()):
        volume = shell["signed_volume_mm3"]
        require(volume < 0 if nesting % 2 else volume > 0, "shell orientation disagrees with nesting")
        shell["nesting_depth"] = nesting
    return details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stl", type=Path)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        points, faces, digest, accuracy, bound = read_mesh(
            args.stl.read_bytes(), json.loads(args.manifest.read_text())
        )
        result = {
            "status": "passed",
            "libigl": version("libigl"),
            "sha256": digest,
            "accuracy_mm": str(accuracy),
            "recorded_bound_mm": str(bound),
            "vertices": len(points),
            "shells": check_shells(points, faces),
        }
    except (OSError, KeyError, ValueError, RuntimeError) as error:
        print(json.dumps({"status": "failed", "detail": str(error)}))
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
