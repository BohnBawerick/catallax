"""Independent OCCT construction for the public fillet corpus, not STEP repair."""
import json
import math
from fractions import Fraction
from pathlib import Path
import sys
import time

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,
                                BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_Transform)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone,
                           BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism,
                           BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeTorus)
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.GC import GC_MakeArcOfCircle
from OCP.GProp import GProp_GProps
from OCP.GeomAbs import GeomAbs_C0
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID, TopAbs_VERTEX
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt, gp_Trsf, gp_Vec

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'step-regressions'))
from test_loopless_faces import raw_solid, unique


def number(value):
    return float(Fraction(value))


def placed(shape, pose):
    w, x, y, z = map(number, pose['rotation'])
    n = w*w + x*x + y*y + z*z
    tx, ty, tz = map(number, pose['translation'])
    transform = gp_Trsf()
    transform.SetValues(1-2*(y*y+z*z)/n, 2*(x*y-z*w)/n, 2*(x*z+y*w)/n, tx,
                        2*(x*y+z*w)/n, 1-2*(x*x+z*z)/n, 2*(y*z-x*w)/n, ty,
                        2*(x*z-y*w)/n, 2*(y*z+x*w)/n, 1-2*(x*x+y*y)/n, tz)
    return BRepBuilderAPI_Transform(shape, transform, True).Shape()


def witness_observations(shape, witnesses):
    return [{'point': item['point'], 'expected': item['state'],
             'actual': str(BRepClass3d_SolidClassifier(
                 shape, gp_Pnt(*map(number, item['point'])), 1e-7).State()).split('.')[-1]}
            for item in witnesses]


def summary(shape, witnesses=(), analytic_volume=None):
    properties = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-8, False, False)
    volume = properties.Mass()
    tolerances = [BRep_Tool.Tolerance_s(cast(item))
                  for kind, cast in [(TopAbs_FACE, TopoDS.Face), (TopAbs_EDGE, TopoDS.Edge),
                                     (TopAbs_VERTEX, TopoDS.Vertex)]
                  for item in unique(shape, kind)]
    result = {'valid': BRepCheck_Analyzer(shape).IsValid(),
              'witnesses': witness_observations(shape, witnesses),
              'max_tolerance_mm': max(tolerances, default=0),
              'solids': len(unique(shape, TopAbs_SOLID)),
              'faces': len(unique(shape, TopAbs_FACE)), 'volume': volume,
              'positive_volume': math.isfinite(volume) and volume > 0}
    if analytic_volume:
        if analytic_volume['kind'] != 'right-angle-linear-fillet':
            raise ValueError('unknown independent analytic volume control')
        r0, r1 = map(number, analytic_volume['radii'])
        expected = (number(analytic_volume['base'])
                    - number(analytic_volume['length'])*(r0*r0+r0*r1+r1*r1)/3*(1-math.pi/4))
        allowance = max(1e-5, abs(expected)*1e-7)
        result['analytic_volume'] = {'expected': expected, 'difference': volume-expected,
                                     'allowance': allowance,
                                     'matches': math.isfinite(volume) and abs(volume-expected) <= allowance}
    return result


def selected(shape, points):
    edges = [TopoDS.Edge(e) for e in unique(shape, TopAbs_EDGE)]
    if points == 'all':
        faces = unique(shape, TopAbs_FACE)
        result = []
        for edge in edges:
            adjacent = [f for f in faces if any(edge.IsSame(e) for e in unique(f, TopAbs_EDGE))]
            if len(adjacent) == 2 and BRep_Tool.Continuity_s(
                    edge, TopoDS.Face(adjacent[0]), TopoDS.Face(adjacent[1])) == GeomAbs_C0:
                result.append(edge)
        return result
    result = []
    for point in points:
        vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*map(number, point))).Vertex()
        found = []
        for edge in edges:
            distance = BRepExtrema_DistShapeShape(vertex, edge)
            distance.Perform()
            if distance.IsDone() and distance.Value() < 1e-6:
                found.append(edge)
        if len(found) != 1:
            raise ValueError(f'edge witness {point} matches {len(found)} edges, expected one')
        if not any(found[0].IsSame(e) for e in result):
            result.append(found[0])
    if not result:
        raise ValueError('empty OCCT edge selection')
    return result


def native_linear_observation(base, edge, radii, start, witnesses):
    """Record the native endpoint overload without making it the linear oracle."""
    observation = {'route': 'BRepFilletAPI_MakeFillet endpoint overload'}
    try:
        fillet = BRepFilletAPI_MakeFillet(base)
        fillet.Add(*radii, edge)
        fillet.Build()
        if not fillet.IsDone():
            return {**observation, 'status': 'refused-with-reason',
                    'code': 'occt.native-fillet-not-done'}
        result = summary(fillet.Shape(), witnesses)
        curve = BRepAdaptor_Curve(edge)
        a, b = curve.FirstParameter(), curve.LastParameter()
        forward = curve.Value(a).Distance(gp_Pnt(*start)) <= 1e-7
        if not forward and curve.Value(b).Distance(gp_Pnt(*start)) > 1e-7:
            raise ValueError('the linear-law origin is not a native edge endpoint')
        law = fillet.GetLaw(fillet.Contour(edge), edge)
        samples = []
        for fraction in (0, 0.25, 0.5, 0.75, 1):
            parameter = a+(b-a)*(fraction if forward else 1-fraction)
            samples.append({'fraction': fraction, 'native_radius': law.Value(parameter),
                            'requested_linear_radius': radii[0]+fraction*(radii[1]-radii[0])})
        same_law = all(abs(s['native_radius']-s['requested_linear_radius']) <= 1e-9
                       for s in samples)
        return {**observation, 'status': 'built' if same_law and acceptable(result)
                else 'different-radius-law' if not same_law else 'wrong',
                'law_type': type(law).__name__, 'same_requested_law': same_law,
                'law_samples': samples, 'result': result}
    except Exception as error:
        return {**observation, 'status': 'refused-with-reason',
                'code': 'occt.native-variable-exception', 'detail': str(error)}


def ruled_linear_reference(base, step, case):
    """Public OCCT loft/cut construction, independent of Catallax's control net."""
    name = step['id']
    spec = case['occt']['ruled_linear_sections'][name]
    [edge] = selected(base, case['occt']['edges'][name])
    [record] = step['edges']
    if record['radius']['kind'] != 'linear':
        raise ValueError('the ruled-section oracle requires an authored linear law')
    radii = [number(record['radius'][end]) for end in ('start', 'end')]
    if radii != list(map(number, case['occt']['linear_radius'][name])):
        raise ValueError('the native observation and authored endpoint radii differ')
    if min(radii) < 0 or max(radii) <= 0:
        raise ValueError('a linear reference needs non-negative, not both zero radii')
    points = [tuple(map(number, spec[end])) for end in ('start', 'end')]
    along = [b-a for a, b in zip(*points)]
    length = math.sqrt(sum(x*x for x in along))
    normals = [tuple(map(number, n)) for n in spec['normals']]
    normals = [tuple(x/math.sqrt(sum(y*y for y in n)) for x in n) for n in normals]
    dot = sum(a*b for a, b in zip(*normals))
    if length <= 0 or not -1 < dot < 1 or any(
            abs(sum(a*b for a, b in zip(n, along))) > 1e-10 for n in normals):
        raise ValueError('reference section normals must be transverse to a straight edge')
    curve = BRepAdaptor_Curve(edge)
    actual = [curve.Value(curve.FirstParameter()), curve.Value(curve.LastParameter())]
    if not all(any(p.Distance(gp_Pnt(*point)) <= 1e-7 for p in actual) for point in points):
        raise ValueError('reference section endpoints must be the selected edge endpoints')
    shift = [-(a+b)/(1+dot) for a, b in zip(*normals)]
    bisector = [a+b for a, b in zip(*normals)]
    magnitude = math.sqrt(sum(x*x for x in bisector))
    bisector = [x/magnitude for x in bisector]
    loft = BRepOffsetAPI_ThruSections(True, True, 1e-7)
    loft.CheckCompatibility(False)
    for corner, radius in zip(points, radii):
        if radius == 0:
            loft.AddVertex(BRepBuilderAPI_MakeVertex(gp_Pnt(*corner)).Vertex())
            continue
        centre = [p+radius*d for p, d in zip(corner, shift)]
        first, last = [gp_Pnt(*(c+radius*n for c, n in zip(centre, normal)))
                       for normal in normals]
        middle = gp_Pnt(*(c+radius*n for c, n in zip(centre, bisector)))
        wire = BRepBuilderAPI_MakeWire()
        wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(*corner), first).Edge())
        wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(first, middle, last).Value()).Edge())
        wire.Add(BRepBuilderAPI_MakeEdge(last, gp_Pnt(*corner)).Edge())
        loft.AddWire(wire.Wire())
    loft.Build()
    if not loft.IsDone():
        raise ValueError('OCCT ruled contact-section construction did not complete')
    cut = BRepAlgoAPI_Cut(base, loft.Shape())
    cut.Build()
    if not cut.IsDone():
        raise ValueError('OCCT ruled contact-section subtraction did not complete')
    observation = native_linear_observation(base, edge, radii, points[0], case.get('witnesses', []))
    return cut.Shape(), observation


def build(case, recipe):
    bodies, sketches, operations, native_observations = {}, {}, [], []
    for step in recipe['steps']:
        name, op = step['id'], step['op']
        started = time.perf_counter()
        if op == 'sketch':
            sketches[name] = step
            continue
        if op == 'block':
            shape = BRepPrimAPI_MakeBox(*map(number, step['size'])).Shape()
        elif op == 'cylinder':
            shape = BRepPrimAPI_MakeCylinder(number(step['radius']), number(step['height'])).Shape()
        elif op == 'cone':
            shape = BRepPrimAPI_MakeCone(number(step['base_radius']), number(step['top_radius']),
                                        number(step['height'])).Shape()
        elif op == 'sphere':
            shape = BRepPrimAPI_MakeSphere(number(step['radius'])).Shape()
        elif op == 'torus':
            shape = BRepPrimAPI_MakeTorus(number(step['major']), number(step['minor'])).Shape()
        elif op == 'extrude':
            sketch = sketches[step['profile']]
            if len(sketch['loops']) != 1 or any(s['kind'] != 'line' for s in sketch['loops'][0]):
                raise ValueError('OCCT corpus adapter requires a single polygon')
            points = [gp_Pnt(*map(number, s['to']), 0) for s in sketch['loops'][0]]
            wire = BRepBuilderAPI_MakeWire()
            for a, b in zip(points, points[1:] + points[:1]):
                wire.Add(BRepBuilderAPI_MakeEdge(a, b).Edge())
            shape = BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(wire.Wire()).Face(),
                                         gp_Vec(0, 0, number(step['distance']))).Shape()
            shape = placed(shape, sketch['plane'])
        elif op == 'boolean':
            cls = {'union': BRepAlgoAPI_Fuse, 'subtract': BRepAlgoAPI_Cut,
                   'intersect': BRepAlgoAPI_Common}[step['kind']]
            boolean = cls(bodies[step['target']], bodies[step['tool']])
            boolean.Build()
            if not boolean.IsDone():
                return {'status': 'refused-with-reason', 'code': 'occt.boolean', 'step': name}
            unify = ShapeUpgrade_UnifySameDomain(boolean.Shape(), True, True, False)
            unify.Build()
            shape = unify.Shape()
        elif op == 'fillet_variable':
            shape, observation = ruled_linear_reference(bodies[step['body']], step, case)
            native_observations.append({'step': name, **observation})
        elif op == 'fillet':
            base = bodies[step['body']]
            edges = selected(base, case['occt']['edges'][name])
            fillet = BRepFilletAPI_MakeFillet(base)
            for edge in edges:
                if fillet.Contour(edge) != 0:
                    continue
                fillet.Add(number(step['radius']), edge)
            try:
                fillet.Build()
            except Exception as error:
                return {'status': 'refused-with-reason', 'code': 'occt.fillet-exception',
                        'detail': str(error), 'step': name}
            if not fillet.IsDone():
                return {'status': 'refused-with-reason', 'code': 'occt.fillet-not-done',
                        'step': name, 'faulty_contours': fillet.NbFaultyContours(),
                        'stripe_status': [str(fillet.StripeStatus(fillet.FaultyContour(i)))
                                          for i in range(1, fillet.NbFaultyContours()+1)]}
            shape = fillet.Shape()
            window = case['occt'].get('setback_window', {}).get(name)
            if window:
                box = BRepPrimAPI_MakeBox(gp_Pnt(*map(number, window['min'])),
                                           *map(number, window['size'])).Shape()
                removed = BRepAlgoAPI_Common(BRepAlgoAPI_Cut(base, shape).Shape(), box).Shape()
                shape = BRepAlgoAPI_Cut(base, removed).Shape()
        else:
            raise ValueError(f'no independent OCCT adapter for {op}')
        if 'pose' in step:
            shape = placed(shape, step['pose'])
        bodies[name] = shape
        operations.append({'step': name, 'op': op,
                           'wall_ms': (time.perf_counter()-started)*1000})
    result = summary(bodies[recipe['result']], case.get('witnesses', []), case.get('analytic_volume'))
    status = 'built' if acceptable(result) else 'wrong'
    report = {'status': status, 'result': result, 'operations': operations}
    if case.get('target_witnesses'):
        report['target_witnesses'] = witness_observations(
            bodies[recipe['result']], case['target_witnesses'])
    if native_observations:
        report.update(reference_route='OCCT ruled linear contact sections',
                      native_fillet_observations=native_observations)
    return report


def acceptable(result):
    return (result['valid'] and result['solids'] == 1 and result['positive_volume']
            and result['max_tolerance_mm'] <= 0.001
            and all(w['actual'] == w['expected'] for w in result['witnesses'])
            and result.get('analytic_volume', {}).get('matches', True))


def main():
    case_path, recipe_path, output = map(Path, sys.argv[1:4])
    started = time.perf_counter()
    try:
        case = json.loads(case_path.read_text())
        if len(sys.argv) == 5:
            witnesses = case.get('witnesses', []) + case.get('target_witnesses', [])
            result = summary(raw_solid(Path(sys.argv[4])), witnesses,
                             case.get('analytic_volume'))
            report = {'status': 'built' if acceptable(result) else 'wrong',
                      'result': result}
        else:
            report = build(case, json.loads(recipe_path.read_text()))
    except Exception as error:
        report = {'status': 'incomplete', 'code': 'harness.occt', 'detail': str(error)}
    import OCP
    report.update(wall_ms=(time.perf_counter()-started)*1000, ocp_version=OCP.__version__)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
