#!/usr/bin/env python3
"""Run authored target scripts against unmodified released Python packages."""
import json
import math
from pathlib import Path
import sys
import time

from targets import build


class CadQuery:
    def __init__(self):
        import cadquery as cq
        self.cq = cq

    def box(self, size, corner=(0, 0, 0)):
        return self.cq.Solid.makeBox(*size, self.cq.Vector(*corner))

    def cylinder(self, radius, height, base=(0, 0, 0), axis=(0, 0, 1)):
        return self.cq.Solid.makeCylinder(radius, height, self.cq.Vector(*base), self.cq.Vector(*axis))

    def sphere(self, radius, center=(0, 0, 0)):
        return self.cq.Solid.makeSphere(radius, self.cq.Vector(*center), angleDegrees1=-90)

    def cone(self, radius, top, height, base=(0, 0, 0), axis=(0, 0, 1)):
        return self.cq.Solid.makeCone(radius, top, height, self.cq.Vector(*base), self.cq.Vector(*axis))

    def torus(self, major, minor):
        return self.cq.Solid.makeTorus(major, minor)

    def rotate_translate(self, part, quaternion, offset):
        w, x, y, z = quaternion
        angle = math.degrees(2 * math.atan2(math.sqrt(x*x+y*y+z*z), w))
        return part.rotate((0, 0, 0), (x, y, z), angle).translate(offset)

    def union(self, a, b):
        return a.fuse(b).clean()

    def cut(self, a, b):
        return a.cut(b).clean()

    def common(self, a, b):
        return a.intersect(b).clean()

    def blend(self, part, points, radius, chamfer=False, curve=None):
        edges = []
        for point in points:
            vertex = self.cq.Vertex.makeVertex(*point)
            candidates = [e for e in part.Edges() if curve is None or e.geomType().lower() == curve]
            distance, edge = min(((e.distance(vertex), e) for e in candidates), key=lambda pair: pair[0])
            if distance > 1e-7:
                raise ValueError(f'edge pick misses target: {point}, {distance}')
            if not any(edge.isSame(old) for old in edges):
                edges.append(edge)
        # CadQuery Boolean results are compounds even when they contain one solid.
        solids = part.Solids()
        if len(solids) != 1:
            raise ValueError('blend requires the complete single target solid')
        return solids[0].chamfer(radius, radius, edges) if chamfer else solids[0].fillet(radius, edges)

    def shell_top(self, part, thickness):
        face = max(part.Faces(), key=lambda f: f.Center().z)
        return part.shell([face], -thickness)

    def coil(self, radius, z, angle, rise, tube):
        # Circle in the axial plane y=0, carried along a right-handed helix about z.
        # isFrenet=False selects OCCT's corrected Frenet frame, its torsion-free frame.
        cq = self.cq
        path = cq.Wire.makeHelix(rise * 360 / angle, rise, radius, cq.Vector(0, 0, z))
        profile = cq.Wire.makeCircle(tube, cq.Vector(radius, 0, z), cq.Vector(0, 1, 0))
        return cq.Solid.sweep(profile, [], path, makeSolid=True, isFrenet=False)

    def loft(self, sections):
        cq = self.cq
        def wire(kind, shape, z):
            if kind == 'circle':
                return cq.Wire.makeCircle(shape, cq.Vector(0, 0, z), cq.Vector(0, 0, 1))
            points = [[cq.Vector(x, y, z) for x, y in segment[1:]] for segment in shape]
            return cq.Wire.assembleEdges([cq.Edge.makeLine(*p) if len(p) == 2 else cq.Edge.makeThreePointArc(*p)
                                          for p in points])
        return cq.Solid.makeLoft([wire(*section) for section in sections], ruled=False)

    def bezier_prism(self, controls, top, y, length):
        cq = self.cq
        (x0, z0), (x3, z3) = controls[0], controls[-1]
        arch = cq.Edge.makeBezier([cq.Vector(x, y, z) for x, z in controls])
        rest = cq.Wire.makePolygon([(x3, y, z3), (x3, y, top), (x0, y, top), (x0, y, z0)])
        wire = cq.Wire.assembleEdges([arch, *rest.Edges()])
        return cq.Solid.extrudeLinear(wire, [], cq.Vector(0, length, 0))

    def blend_above(self, part, z, radius):
        solids = part.Solids()
        if len(solids) != 1:
            raise ValueError('blend requires the complete single target solid')
        edges = [e for e in solids[0].Edges() if e.BoundingBox().zmin > z]
        return solids[0].fillet(radius, edges)

    def import_step(self, path):
        solids = self.cq.importers.importStep(str(path)).solids().vals()
        if len(solids) != 1:
            raise ValueError(f'expected one source solid, read {len(solids)}')
        return solids[0]

    def measure(self, part):
        center = part.Center()
        return {'volume': part.Volume(), 'area': part.Area(), 'centroid': list(center.toTuple())}

    def faces(self, part):
        return len(part.Faces())

    def export(self, part, path):
        self.cq.exporters.export(part, str(path))


class SMLib:
    def __init__(self):
        import smlibpy as sm
        self.sm = sm

    def box(self, size, corner=(0, 0, 0)):
        return self.sm.Box(*size, align=(self.sm.Align.MIN,)*3).translate(*corner)

    def cylinder(self, radius, height, base=(0, 0, 0), axis=(0, 0, 1)):
        part = self.sm.Cylinder(radius, height, align=(self.sm.Align.CENTER, self.sm.Align.CENTER, self.sm.Align.MIN))
        return self.orient(part, axis).translate(*base)

    def cone(self, radius, top, height, base=(0, 0, 0), axis=(0, 0, 1)):
        part = self.sm.Cone(radius, top, height, align=(self.sm.Align.CENTER, self.sm.Align.CENTER, self.sm.Align.MIN))
        return self.orient(part, axis).translate(*base)

    def orient(self, part, axis):
        if axis == (0, 0, 1):
            return part
        if axis == (1, 0, 0):
            return part.rotate(self.sm.Y, 90)
        if axis == (0, 1, 0):
            return part.rotate(self.sm.X, -90)
        if axis == (0, -1, 0):
            return part.rotate(self.sm.X, 90)
        raise NotImplementedError(f'axis {axis}')

    def sphere(self, radius, center=(0, 0, 0)):
        return self.sm.Sphere(radius).translate(*center)

    def torus(self, major, minor):
        return self.sm.Torus(major, minor)

    def rotate_translate(self, part, quaternion, offset):
        w, x, y, z = quaternion
        angle = math.degrees(2 * math.atan2(math.sqrt(x*x+y*y+z*z), w))
        return part.rotate(self.sm.Axis((0, 0, 0), (x, y, z)), angle).translate(*offset)

    def union(self, a, b):
        return a + b

    def cut(self, a, b):
        return a - b

    def common(self, a, b):
        return a & b

    def blend(self, part, points, radius, chamfer=False, curve=None):
        edges = []
        for point in points:
            candidates = part.edges().of_type(curve) if curve else part.edges()
            selection = candidates.at(point)
            if len(selection) != 1:
                raise ValueError(f'edge pick not unique: {point}, {len(selection)}')
            edges.extend(selection)
        return part.chamfer(edges, radius) if chamfer else part.fillet(edges, radius)

    def shell_top(self, part, thickness):
        return part.shell(thickness, part.faces().top())

    def coil(self, radius, z, angle, rise, tube):
        sm = self.sm
        plane = sm.Plane(origin=(radius, 0, z), x_dir=(1, 0, 0), z_dir=(0, 1, 0))
        path = sm.Helix(rise * 360 / angle, rise, radius, center=(0, 0, z))
        return sm.sweep(plane * sm.Circle(tube), path, orient=True)

    def loft(self, sections):
        sm = self.sm
        def wire(kind, shape, z):
            if kind == 'circle':
                return sm.Circle(shape).translate(0, 0, z)
            points = [[(x, y, z) for x, y in segment[1:]] for segment in shape]
            return sm.Wire([sm.Line(*p) if len(p) == 2 else sm.ThreePointArc(*p) for p in points])
        return sm.loft([wire(*section) for section in sections])

    def bezier_prism(self, controls, top, y, length):
        sm = self.sm
        (x0, z0), (x3, z3) = controls[0], controls[-1]
        arch = sm.Bezier(*[(x, y, z) for x, z in controls])
        rest = sm.Polyline((x3, y, z3), (x3, y, top), (x0, y, top), (x0, y, z0))
        return sm.extrude(sm.make_face(sm.Wire([arch, rest])), abs(length),
                          direction=(0, 1 if length > 0 else -1, 0))

    def blend_above(self, part, z, radius):
        edges = [e for e in part.edges() if e.bbox().min.z > z]
        return part.fillet(edges, radius)

    def import_step(self, path):
        part = self.sm.import_step(str(path))
        solids = part.solids() if hasattr(part, 'solids') else [part]
        if len(solids) != 1:
            raise ValueError(f'expected one source solid, read {len(solids)}')
        return solids[0]

    def measure(self, part):
        center = part.center()
        return {'volume': part.volume, 'area': part.area, 'centroid': [center.x, center.y, center.z]}

    def faces(self, part):
        return len(part.faces())

    def export(self, part, path):
        part.export(str(path))


def run(engine, case, directory):
    api = {'cadquery': CadQuery, 'smlibpy': SMLib}[engine]()
    result = {'status': 'completed', 'variants': {}, 'edit_mode': 'complete target script rebuild, no incremental reuse',
              'edit_survival': {'status': 'not-run', 'reason': 'released API has no persistent per-face/edge naming contract'}}
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    def checkpoint():
        (directory / 'phases.json').write_text(json.dumps(result, allow_nan=False))
    def phase(row, key, fn):
        started = time.perf_counter()
        try:
            value = fn()
            row[key] = {'status': 'answered', 'wall_ms': (time.perf_counter()-started)*1000}
            return value
        except Exception as error:
            row[key] = {'status': 'not-expressible' if isinstance(error, NotImplementedError) else 'command-failed',
                        'wall_ms': None, 'refusal_upper_bound_ms': (time.perf_counter()-started)*1000,
                        'detail': f'{type(error).__name__}: {error}'}
            return None
        finally:
            checkpoint()
    for label in ('construction', 'edit'):
        row = result['variants'][label] = {}
        edit = float(case['edit']['value']) if label == 'edit' else None
        part = phase(row, 'build', lambda: build(api, case['recipe'], edit, case))
        if part is None:
            continue
        phase(row, 'display', lambda: api.export(part, directory / f'{label}.stl'))
        measured = phase(row, 'measurement', lambda: api.measure(part))
        if measured is not None:
            row['measurement']['claimed'] = measured
        phase(row, 'step', lambda: api.export(part, directory / f'{label}.step'))
        row['step']['expected_faces'] = api.faces(part)
        checkpoint()
    return result


if __name__ == '__main__':
    run(sys.argv[1], json.loads(Path(sys.argv[2]).read_text()), sys.argv[3])
