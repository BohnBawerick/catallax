"""Construct complete cases with released public APIs. Never export a prefix."""
import json
import math
from fractions import Fraction
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / '.claude/skills/verify-catallax/scripts/head-to-head'
sys.path.insert(0, str(BENCH))
from python_tools import CadQuery, SMLib
from targets import build as benchmark_build



class NativeRefusal(RuntimeError):
    pass

def number(value):
    return float(Fraction(value))


def coordinates(value):
    if hasattr(value, 'toTuple'):
        return tuple(value.toTuple())
    return (value.x, value.y, value.z)


def placed(api, part, pose):
    q = tuple(map(number, pose['rotation']))
    t = tuple(map(number, pose['translation']))
    if q[1:] == (0, 0, 0):
        return part.translate(t) if isinstance(api, CadQuery) else part.translate(*t)
    return api.rotate_translate(part, q, t)


def faces(api, part):
    return part.Faces() if isinstance(api, CadQuery) else list(part.faces())


def center(api, face):
    return coordinates(face.Center() if isinstance(api, CadQuery) else face.center())


def normal(api, face):
    try:
        return coordinates(face.normalAt() if isinstance(api, CadQuery) else face.normal_at())
    except Exception:
        return None


def face_type(api, face):
    return face.geomType().lower() if isinstance(api, CadQuery) else face.geom_type.lower()


def distance(api, face, point):
    return face.distance(api.cq.Vertex.makeVertex(*point)) if isinstance(api, CadQuery) else face.distance_to(point)


def references(api, shape, step):
    result = {}
    fs = faces(api, shape)
    result['all'] = fs
    curved = [f for f in fs if face_type(api, f) != 'plane']
    if curved:
        result['wall'] = curved
    for label, axis, value in (('start', 2, -1), ('end', 2, 1),
                               ('x_min', 0, -1), ('x_max', 0, 1),
                               ('y_min', 1, -1), ('y_max', 1, 1)):
        selected = []
        for face in fs:
            n = normal(api, face)
            if face_type(api, face) == 'plane' and n and n[axis]*value > 1-1e-7:
                selected.append(face)
        if selected:
            result[label] = selected
    return result


def selected_faces(api, part, selectors, origins):
    result = []
    for selector in selectors:
        origin = selector['face']
        role = origin.get('role')
        label = 'all' if role is None else role if isinstance(role, str) else role['side']['segment']
        originals = origins[origin['step']].get(label)
        if not originals:
            raise NotImplementedError(f'no geometric witness for face origin {origin}')
        for face in faces(api, part):
            if any(face_type(api, face) == face_type(api, old)
                   and distance(api, old, surface_point(api, face)) <= 1e-6 for old in originals):
                if face not in result:
                    result.append(face)
    if not result:
        raise ValueError('authored face selection has no surviving geometric match')
    return result


def surface_point(api, face):
    if not isinstance(api, CadQuery):
        return coordinates(face.interior_point())
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.TopAbs import TopAbs_IN
    from OCP.gp import gp_Pnt2d
    surface = BRepAdaptor_Surface(face.wrapped)
    u0, u1 = surface.FirstUParameter(), surface.LastUParameter()
    v0, v1 = surface.FirstVParameter(), surface.LastVParameter()
    fractions = (0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875)
    for a in fractions:
        for b in fractions:
            u, v = u0+(u1-u0)*a, v0+(v1-v0)*b
            if BRepClass_FaceClassifier(face.wrapped, gp_Pnt2d(u, v), 1e-9).State() == TopAbs_IN:
                p = surface.Value(u, v)
                return (p.X(), p.Y(), p.Z())
    raise NotImplementedError('no trim-checked interior face witness found')


def profile(api, sketch):
    loops = []
    for loop in sketch['loops']:
        if any(s['kind'] != 'line' for s in loop):
            raise NotImplementedError('non-polygon sketch adapter not implemented')
        points = [tuple(map(number, s['to'])) for s in loop]
        if isinstance(api, CadQuery):
            wire = api.cq.Wire.makePolygon([api.cq.Vector(x, y, 0) for x, y in points], close=True)
            loops.append(wire)
        else:
            loops.append(api.sm.Polygon(*points))
    return loops


def extrude(api, sketch, height):
    loops = profile(api, sketch)
    if isinstance(api, CadQuery):
        shape = api.cq.Solid.extrudeLinear(loops[0], loops[1:], api.cq.Vector(0, 0, height))
    else:
        shape = api.sm.extrude(loops[0], height)
        for loop in loops[1:]:
            shape = api.cut(shape, api.sm.extrude(loop, height))
    return placed(api, shape, sketch['plane'])


def generic(api, recipe, fillet_case=None):
    bodies, sketches, origins = {}, {}, {}
    for step in recipe['steps']:
        name, op = step['id'], step['op']
        if op == 'sketch':
            sketches[name] = step
            continue
        if op == 'block':
            shape = api.box(tuple(map(number, step['size'])))
        elif op == 'cylinder':
            shape = api.cylinder(number(step['radius']), number(step['height']))
        elif op == 'sphere':
            shape = api.sphere(number(step['radius']))
        elif op == 'cone':
            shape = api.cone(number(step['base_radius']), number(step['top_radius']), number(step['height']))
        elif op == 'torus':
            shape = api.torus(number(step['major']), number(step['minor']))
        elif op == 'extrude':
            sketch = sketches[step['profile']]
            raw = dict(sketch, plane={'rotation': ['1', '0', '0', '0'], 'translation': ['0', '0', '0']})
            height = number(step['distance'])
            shape = extrude(api, raw, height)
            origins[name] = references(api, shape, step)
            for loop in sketch['loops']:
                for index, segment in enumerate(loop):
                    a, b = loop[index-1]['to'], segment['to']
                    point = ((number(a[0])+number(b[0]))/2, (number(a[1])+number(b[1]))/2, height/2)
                    origins[name][segment['id']] = [f for f in faces(api, shape) if distance(api, f, point) <= 1e-6]
            shape = placed(api, shape, sketch['plane'])
            origins[name] = {role: [placed(api, f, sketch['plane']) for f in fs]
                             for role, fs in origins[name].items()}
        elif op == 'boolean':
            fn = {'union': api.union, 'subtract': api.cut, 'intersect': api.common}[step['kind']]
            shape = fn(bodies[step['target']], bodies[step['tool']])
        elif op == 'revolve':
            sketch = sketches[step['profile']]
            loops = [placed(api, loop, sketch['plane']) for loop in profile(api, sketch)]
            start = tuple(map(number, step['axis']['point']))
            direction = tuple(map(number, step['axis']['direction']))
            angle = number(step['angle'])
            if isinstance(api, CadQuery):
                end = tuple(a+b for a, b in zip(start, direction))
                shape = api.cq.Solid.revolve(loops[0], loops[1:], angle, start, end)
            else:
                axis = api.sm.Axis(start, direction)
                shape = api.sm.revolve(loops[0], axis, angle)
                for loop in loops[1:]:
                    shape = api.cut(shape, api.sm.revolve(loop, axis, angle))
        elif op == 'extrude_through':
            target, sketch = bodies[step['target']], sketches[step['profile']]
            bounds = target.BoundingBox() if isinstance(api, CadQuery) else target.bbox()
            lower = bounds.zmin if isinstance(api, CadQuery) else bounds.min.z
            upper = bounds.zmax if isinstance(api, CadQuery) else bounds.max.z
            shifted = dict(sketch, plane=dict(sketch['plane'], translation=[
                sketch['plane']['translation'][0], sketch['plane']['translation'][1], str(lower-1)]))
            shape = api.cut(target, extrude(api, shifted, upper-lower+2))
        elif op == 'hole':
            shape = bodies[step['target']]
            bounds = shape.BoundingBox() if isinstance(api, CadQuery) else shape.bbox()
            length = bounds.DiagonalLength if isinstance(api, CadQuery) else bounds.diagonal
            for placement in step['at']:
                direction = tuple(map(number, placement['direction']))
                start = tuple(number(x)-d for x, d in zip(placement['point'], direction))
                depth = length+2 if step['depth'] == 'through' else number(step['depth'])+1
                if isinstance(api, CadQuery):
                    tool = api.cylinder(number(step['diameter'])/2, depth, start, direction)
                else:
                    tool = api.cylinder(number(step['diameter'])/2, depth)
                    if direction == (0, 0, -1):
                        tool = tool.rotate(api.sm.X, 180)
                    elif direction != (0, 0, 1):
                        raise NotImplementedError('hole direction is outside the authored adapter')
                    tool = tool.translate(*start)
                shape = api.cut(shape, tool)
        elif op in ('fillet', 'fillet_variable'):
            base = bodies[step['body']]
            if fillet_case is None:
                raise NotImplementedError('fillet case needs authored geometric edge witnesses')
            picks = fillet_case['occt']['edges'][name]
            if picks == 'all':
                edges = list(base.edges())
            else:
                edges = []
                for point in picks:
                    selection = base.edges().at(tuple(map(number, point)))
                    if not selection:
                        raise ValueError('authored edge witness selects no edge')
                    if len(selection) > 1:
                        first = selection[0]
                        same_circle = all(
                            e.geom_type == first.geom_type == 'circle'
                            and math.dist(coordinates(e.arc_center), coordinates(first.arc_center)) <= 1e-7
                            and abs(e.radius-first.radius) <= 1e-7
                            and abs(abs(sum(a*b for a, b in zip(coordinates(e.normal), coordinates(first.normal))))-1) <= 1e-7
                            for e in selection)
                        if not same_circle:
                            raise ValueError('authored edge witness matches distinct geometric contours')
                    edges.extend(selection)
            if op == 'fillet_variable':
                record, = step['edges']
                radius = tuple(number(record['radius'][key]) for key in ('start', 'end'))
            else:
                radius = number(step['radius'])
            shape = base.fillet(edges, radius)
            window = fillet_case['occt'].get('setback_window', {}).get(name)
            if window:
                region = api.box(tuple(map(number, window['size'])), tuple(map(number, window['min'])))
                shape = api.cut(base, api.common(api.cut(base, shape), region))
        elif op == 'offset':
            base, d = bodies[step['body']], number(step['distance'])
            if isinstance(api, CadQuery):
                from OCP.BRepOffsetAPI import BRepOffsetAPI_MakeOffsetShape
                from OCP.BRepOffset import BRepOffset_Skin
                from OCP.GeomAbs import GeomAbs_Arc
                builder = BRepOffsetAPI_MakeOffsetShape()
                builder.PerformByJoin(base.wrapped, d, 1e-7, BRepOffset_Skin, False, False, GeomAbs_Arc)
                if not builder.IsDone():
                    raise NativeRefusal('OCCT offset not done')
                shape = api.cq.Shape.cast(builder.Shape())
            else:
                shape = base.offset(d)
        elif op == 'shell':
            base = bodies[step['body']]
            opening = selected_faces(api, base, step['open'], origins) if step.get('open') else []
            if isinstance(api, CadQuery):
                solids = base.Solids()
                if len(solids) != 1:
                    raise ValueError('shell requires the complete single target solid')
                base = solids[0]
            d = number(step['thickness'])
            shape = base.shell(opening, -d) if isinstance(api, CadQuery) else base.shell(d, opening)
        elif op == 'thicken':
            base = bodies[step['body']]
            selection = selected_faces(api, base, step['faces'], origins)
            d = number(step['thickness'])
            side = step.get('side', 'front')
            if isinstance(api, CadQuery):
                from OCP.BRepOffsetAPI import BRepOffsetAPI_MakeThickSolid
                from OCP.TopTools import TopTools_ListOfShape
                from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing
                sewing = BRepBuilderAPI_Sewing(1e-7)
                for face in selection:
                    sewing.Add(face.wrapped)
                sewing.Perform()
                sheet = api.cq.Shape.cast(sewing.SewedShape())
                def thick(amount):
                    builder = BRepOffsetAPI_MakeThickSolid()
                    builder.MakeThickSolidBySimple(sheet.wrapped, amount)
                    if not builder.IsDone():
                        raise NativeRefusal('OCCT thicken not done')
                    return api.cq.Shape.cast(builder.Shape())
                shape = api.union(thick(d/2), thick(-d/2)) if side == 'symmetric' else thick(-d if side == 'back' else d)
            else:
                sheet = selection[0].to_shell()
                for face in selection[1:]:
                    sheet = sheet + face.to_shell()
                shape = sheet.thicken(-d if side == 'back' else d, both=side == 'symmetric')
        elif op == 'extrude_next':
            sketch = sketches[step['profile']]
            target = bodies[step['target']]
            bounds = target.BoundingBox() if isinstance(api, CadQuery) else target.bbox()
            top = bounds.zmax if isinstance(api, CadQuery) else bounds.max.z
            source_z = number(sketch['plane']['translation'][2])
            prism = extrude(api, sketch, top-source_z+1)
            remainder = api.cut(prism, target)
            solids = remainder.Solids() if isinstance(api, CadQuery) else list(remainder.solids())
            shape = min(solids, key=lambda s: center(api, s)[2])
        else:
            raise NotImplementedError(f'no equivalent public recipe adapter for {op}')
        if op in ('block', 'cylinder', 'sphere', 'cone', 'torus'):
            origins[name] = references(api, shape, step)
        if 'pose' in step:
            shape = placed(api, shape, step['pose'])
            if name in origins:
                origins[name] = {role: [placed(api, f, step['pose']) for f in fs]
                                 for role, fs in origins[name].items()}
        bodies[name] = shape
    return bodies[recipe['result']]


def run(engine, case, output):
    if engine == 'occt' and case.get('fillet_case'):
        sys.path.insert(0, str(ROOT / 'tools/fillet-coverage'))
        import occt
        report = occt.build(case['fillet_case'], case['recipe'], output)
        return {'status': 'answered' if output.exists() else 'refused', 'native': report}
    api = CadQuery() if engine == 'occt' else SMLib()
    if case.get('benchmark_case'):
        source = case['benchmark_case']
        edit = number(source['edit']['value']) if case['variant'] == 'edit' else None
        part = benchmark_build(api, source['recipe'], edit)
    else:
        part = generic(api, case['recipe'], case.get('fillet_case'))
    api.export(part, output)
    return {'status': 'answered'}


if __name__ == '__main__':
    engine, case_path, report_path, output_path = sys.argv[1:]
    started = time.perf_counter()
    try:
        result = run(engine, json.loads(Path(case_path).read_text()), Path(output_path))
    except NotImplementedError as error:
        result = {'status': 'adapter-unavailable', 'detail': str(error), 'kernel_refusal': False}
    except Exception as error:
        frames = error.__traceback__
        released_api = isinstance(error, NativeRefusal)
        while frames is not None:
            released_api |= frames.tb_frame.f_globals.get('__name__', '').startswith(('smlibpy', 'OCP', 'cadquery'))
            frames = frames.tb_next
        result = {'status': 'refused' if released_api else 'execution-error',
                  'detail': f'{type(error).__module__}.{type(error).__name__}: {error}',
                  'refusal_scope': 'released public API; not a geometric impossibility claim' if released_api else None}
    result['kernel_wall_ms'] = (time.perf_counter()-started)*1000
    Path(report_path).write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
