"""Independent OCCT construction for the public fillet corpus, not STEP repair."""
import json
import math
from fractions import Fraction
from pathlib import Path
import sys
import time

from OCP.BRep import BRep_Tool
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


def summary(shape, witnesses=()):
    properties = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-8, False, False)
    volume = properties.Mass()
    tolerances = [BRep_Tool.Tolerance_s(cast(item))
                  for kind, cast in [(TopAbs_FACE, TopoDS.Face), (TopAbs_EDGE, TopoDS.Edge),
                                     (TopAbs_VERTEX, TopoDS.Vertex)]
                  for item in unique(shape, kind)]
    observations = [{'point': item['point'], 'expected': item['state'],
                     'actual': str(BRepClass3d_SolidClassifier(
                         shape, gp_Pnt(*map(number, item['point'])), 1e-7).State()).split('.')[-1]}
                    for item in witnesses]
    return {'valid': BRepCheck_Analyzer(shape).IsValid(), 'witnesses': observations,
            'max_tolerance_mm': max(tolerances, default=0),
            'solids': len(unique(shape, TopAbs_SOLID)),
            'faces': len(unique(shape, TopAbs_FACE)), 'volume': volume,
            'positive_volume': math.isfinite(volume) and volume > 0}


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


def build(case, recipe):
    bodies, sketches, operations = {}, {}, []
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
        elif op == 'fillet':
            base = bodies[step['body']]
            edges = selected(base, case['occt']['edges'][name])
            fillet = BRepFilletAPI_MakeFillet(base)
            for edge in edges:
                if fillet.Contour(edge) != 0:
                    continue
                law = case['occt'].get('linear_radius', {}).get(name)
                if law:
                    fillet.Add(*map(number, law), edge)
                else:
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
    result = summary(bodies[recipe['result']], case.get('witnesses', []))
    status = 'built' if acceptable(result) else 'wrong'
    return {'status': status, 'result': result, 'operations': operations}


def acceptable(result):
    return (result['valid'] and result['solids'] == 1 and result['positive_volume']
            and result['max_tolerance_mm'] <= 0.001
            and all(w['actual'] == w['expected'] for w in result['witnesses']))


def main():
    case_path, recipe_path, output = map(Path, sys.argv[1:4])
    started = time.perf_counter()
    try:
        case = json.loads(case_path.read_text())
        if len(sys.argv) == 5:
            result = summary(raw_solid(Path(sys.argv[4])), case.get('witnesses', []))
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
