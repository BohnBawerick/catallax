"""Observe every producer's raw STEP without healing or changing tolerances."""
from fractions import Fraction
import json
import math
from pathlib import Path
import sys

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass import BRepClass_FaceClassifier
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_IN, TopAbs_ON
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.gp import gp_Pnt, gp_Pnt2d


def unique(shape, kind):
    result = []
    explorer = TopExp_Explorer(shape, kind)
    while explorer.More():
        item = explorer.Current()
        if not any(item.IsSame(old) for old in result):
            result.append(item)
        explorer.Next()
    return result


def number(value):
    return float(Fraction(value))


def point(p):
    return [str(Fraction.from_float(v)) for v in (p.X(), p.Y(), p.Z())]


def observe(step, case, oracle):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(step)) != IFSelect_RetDone or reader.TransferRoots() == 0:
        return {'valid': False, 'reason': 'raw STEP read or root transfer failed'}
    shape = reader.OneShape()
    fs = [TopoDS.Face_s(f) for f in unique(shape, TopAbs_FACE)]
    solids = unique(shape, TopAbs_SOLID)
    properties = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-8, False, False)
    volume = properties.Mass()
    result = {'valid': not shape.IsNull() and BRepCheck_Analyzer(shape).IsValid()
              and bool(solids) and math.isfinite(volume) and volume > 0,
              'solids': len(solids), 'faces': len(fs), 'volume_mm3': volume,
              'checker': 'OpenCASCADE 7.8.1 via cadquery-ocp 7.8.1.1.post1',
              'checker_limitation': 'OCCT also competes; raw validity is not an independent geometric accuracy proof',
              'max_topological_tolerance_mm': max(
                  [BRep_Tool.Tolerance_s(cast(item)) for kind, cast in (
                      (TopAbs_FACE, TopoDS.Face_s), (TopAbs_EDGE, TopoDS.Edge_s),
                      (TopAbs_VERTEX, TopoDS.Vertex_s)) for item in unique(shape, kind)], default=0)}
    analytic = case.get('fillet_case', {}).get('reference_volume')
    if analytic is not None:
        expected = float(analytic['volume'])
        result['analytic_volume_check'] = {'expected_mm3': expected,
            'error_mm3': volume-expected, 'method': analytic['method'],
            'scope': 'independent scalar diagnostic; not a geometric distance bound'}
    classifier = [BRepClass3d_SolidClassifier(s) for s in solids]
    witnesses = case.get('fillet_case', {}).get('witnesses', []) + case.get('fillet_case', {}).get('target_witnesses', [])
    failures = []
    for witness in witnesses:
        states = []
        for solid in classifier:
            solid.Perform(gp_Pnt(*map(number, witness['point'])), 1e-7)
            states.append(str(solid.State()).split('.')[-1])
        state = 'TopAbs_IN' if 'TopAbs_IN' in states else 'TopAbs_ON' if 'TopAbs_ON' in states else 'TopAbs_OUT'
        if state != witness['state']:
            failures.append({**witness, 'actual': state})
    result['witness_failures'] = failures
    boundary = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(boundary)
    for face in fs:
        builder.Add(boundary, face)
    forward = []
    for p in oracle.get('points', []):
        distance = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*map(number, p))).Vertex(), boundary)
        distance.SetDeflection(1e-7)
        distance.Perform()
        if not distance.IsDone():
            return {**result, 'reason': 'point-to-boundary observer did not finish'}
        forward.append(distance.Value())
    result['forward_max_mm'] = max(forward, default=0)
    result['forward_numerical_allowance_mm'] = 1e-7
    result['forward_lower_mm'] = max(0, result['forward_max_mm']-1e-7)
    result['forward_upper_mm'] = result['forward_max_mm']+1e-7
    samples = []
    for face in fs:
        surface = BRepAdaptor_Surface(face)
        u0, u1 = surface.FirstUParameter(), surface.LastUParameter()
        v0, v1 = surface.FirstVParameter(), surface.LastVParameter()
        if not all(math.isfinite(x) and abs(x) < 1e50 for x in (u0, u1, v0, v1)):
            return {**result, 'reason': 'unbounded face parameter domain'}
        for i in range(5):
            for j in range(5):
                u, v = u0+(u1-u0)*i/4, v0+(v1-v0)*j/4
                state = BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-9).State()
                if state in (TopAbs_IN, TopAbs_ON):
                    samples.append(point(surface.Value(u, v)))
    for edge in unique(shape, TopAbs_EDGE):
        curve = BRepAdaptor_Curve(TopoDS.Edge_s(edge))
        a, b = curve.FirstParameter(), curve.LastParameter()
        if all(math.isfinite(x) and abs(x) < 1e50 for x in (a, b)):
            for i in range(5):
                samples.append(point(curve.Value(a+(b-a)*i/4)))
    for vertex in unique(shape, TopAbs_VERTEX):
        samples.append(point(BRep_Tool.Pnt_s(TopoDS.Vertex_s(vertex))))
    result['points'] = [list(p) for p in dict.fromkeys(tuple(p) for p in samples)]
    return result


if __name__ == '__main__':
    step, case, oracle, output = map(Path, sys.argv[1:])
    result = observe(step, json.loads(case.read_text()), json.loads(oracle.read_text()))
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
