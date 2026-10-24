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
        part = phase(row, 'build', lambda: build(api, case['recipe'], edit))
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
