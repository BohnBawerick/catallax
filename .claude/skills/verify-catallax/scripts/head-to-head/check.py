#!/usr/bin/env python3
"""Independent OpenCASCADE constructions from the target dimensions, not imported answers."""
import json
import math
from pathlib import Path
import sys


def check(case, variant, step_path, expected_faces):
    import OCP
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut, BRepAlgoAPI_Common
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform, BRepBuilderAPI_MakeVertex
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet, BRepFilletAPI_MakeChamfer
    from OCP.BRepGProp import BRepGProp
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeCone, BRepPrimAPI_MakeTorus
    from OCP.GProp import GProp_GProps
    from OCP.gp import gp_Pnt, gp_Dir, gp_Ax2, gp_Trsf, gp_Quaternion, gp_Vec
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SOLID, TopAbs_IN, TopAbs_ON
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    def entries(shape, kind):
        result = []
        ex = TopExp_Explorer(shape, kind)
        while ex.More():
            obj = ex.Current()
            if not any(obj.IsSame(old) for old in result):
                result.append(obj)
            ex.Next()
        return result

    def props(shape):
        v, a = GProp_GProps(), GProp_GProps()
        volume_error = BRepGProp.VolumeProperties_s(shape, v, 1e-10, False, False)
        area_error = BRepGProp.SurfaceProperties_s(shape, a, 1e-10, False)
        if max(volume_error, area_error) > 1e-9:
            raise ValueError('independent integration accuracy not reached')
        c = v.CentreOfMass()
        return {'volume': v.Mass(), 'area': a.Mass(), 'centroid': [c.X(), c.Y(), c.Z()],
                'relative_integration_estimates': {'volume': volume_error, 'area': area_error}}

    def box(size, corner=(0, 0, 0)):
        return BRepPrimAPI_MakeBox(gp_Pnt(*corner), *size).Shape()

    def cylinder(r, h, base=(0, 0, 0), axis=(0, 0, 1)):
        return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*base), gp_Dir(*axis)), r, h).Shape()

    def cone(r, t, h, base=(0, 0, 0), axis=(0, 0, 1)):
        return BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(*base), gp_Dir(*axis)), r, t, h).Shape()

    def boolean(a, b, kind='union'):
        builder = {'union': BRepAlgoAPI_Fuse, 'subtract': BRepAlgoAPI_Cut, 'intersect': BRepAlgoAPI_Common}[kind](a, b)
        builder.Build()
        if not builder.IsDone():
            raise ValueError('independent Boolean failed')
        return builder.Shape()

    def blend(shape, points, radius, chamfer=False):
        builder = BRepFilletAPI_MakeChamfer(shape) if chamfer else BRepFilletAPI_MakeFillet(shape)
        edges = entries(shape, TopAbs_EDGE)
        chosen = []
        for point in points:
            vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*point)).Shape()
            distances = [(BRepExtrema_DistShapeShape(vertex, edge).Value(), edge) for edge in edges]
            distance, edge = min(distances, key=lambda item: item[0])
            if distance > 1e-7:
                raise ValueError(f'independent blend pick misses edge: {point}, {distance}')
            if not any(edge.IsSame(old) for old in chosen):
                chosen.append(edge)
                builder.Add(radius, TopoDS.Edge_s(edge))
        builder.Build()
        if not builder.IsDone():
            raise ValueError('independent blend failed')
        return builder.Shape()

    name = case['recipe']
    edit = float(case['edit']['value']) if variant == 'edit' else None
    if name.startswith('P0/'):
        family = Path(name).stem
        h = edit if edit is not None and family not in ('block-sphere', 'general-torus') else None
        if family == 'block-sphere':
            oracle = boolean(box((edit or 10, 10, 10)), BRepPrimAPI_MakeSphere(gp_Pnt(10, 5, 5), 6).Shape())
        elif family == 'crossing-pipes':
            a = boolean(cylinder(5, h or 20, (0, 0, -10)), cylinder(3, 22, (0, 0, -11)), 'subtract')
            b = boolean(cylinder(5, 20, (-10, 0, 0), (1, 0, 0)), cylinder(3, 22, (-11, 0, 0), (1, 0, 0)), 'subtract')
            oracle = boolean(a, b)
        elif family == 'equal-cylinder-sphere':
            oracle = boolean(cylinder(5, h or 20, (0, 0, -10)), BRepPrimAPI_MakeSphere(5).Shape())
        elif family in ('line-tangent', 'near-miss', 'seam'):
            x = {'line-tangent': 10, 'near-miss': 10.0000001, 'seam': 4}[family]
            oracle = boolean(cylinder(5, h or 10), cylinder(5, 10, (x, 0, 0)))
        elif family in ('general-cone', 'general-torus'):
            a = cone(5, 2, h or 10) if family == 'general-cone' else BRepPrimAPI_MakeTorus(edit or 8, 2).Shape()
            b = cone(5, 2, 10) if family == 'general-cone' else BRepPrimAPI_MakeTorus(8, 2).Shape()
            tr = gp_Trsf()
            norm = math.sqrt(6)
            tr.SetRotation(gp_Quaternion(1 / norm, 1 / norm, 0, 2 / norm))
            tr.SetTranslationPart(gp_Vec(*( (2, 1, 1) if family == 'general-cone' else (1, 2, 1))))
            oracle = boolean(a, BRepBuilderAPI_Transform(b, tr, True).Shape())
        else:
            raise ValueError(f'unknown P0 {family}')
    elif name.startswith('P1-'):
        oracle = boolean(box((60, 40, 5)), box((60, 5, 25), (0, 0, 5)))
        oracle = boolean(oracle, cylinder(edit or 8, 10, (30, 25, 5)))
        oracle = boolean(oracle, cylinder(4, 17, (30, 25, -1)), 'subtract')
        for x in (10, 50):
            for y in (12, 32):
                oracle = boolean(oracle, cylinder(2, 7, (x, y, -1)), 'subtract')
        # This capsule equals the eight overlapping round tools plus their bridge.
        slot = boolean(cylinder(3, 7, (23, 6, 18), (0, -1, 0)), cylinder(3, 7, (37, 6, 18), (0, -1, 0)))
        slot = boolean(slot, box((14, 7, 6), (23, -1, 15)))
        oracle = boolean(oracle, slot, 'subtract')
    elif name.startswith('P2-'):
        oracle = boolean(cylinder(12, 12, axis=(0, 1, 0)), cylinder(6, 12, axis=(0, 1, 0)), 'subtract')
        oracle = blend(oracle, [(12, 0, 0), (12, 12, 0)], edit or 1, chamfer=True)
    elif name.startswith('P3-'):
        oracle = boolean(box((40, 30, 10)), cylinder(8, 12, (20, 15, 10)))
        oracle = blend(oracle, [(28, 15, 10)], edit or 2)
        oracle = blend(oracle, [(20, 0, 10), (40, 15, 10), (20, 30, 10), (0, 15, 10)], 1)
    elif name.startswith('P4-'):
        t = edit or 2
        oracle = boolean(box((40, 30, 12), (-20, -15, 0)), box((40 - 2*t, 30 - 2*t, 13), (-20+t, -15+t, t)), 'subtract')
        oracle = boolean(oracle, cylinder(23, 14, (0, 0, -1)), 'intersect')
    elif name.startswith('P5-'):
        oracle = boolean(cylinder(22, 5, axis=(0, 1, 0)), cylinder(12, 30, axis=(0, 1, 0)))
        oracle = boolean(oracle, cylinder(8, 30, axis=(0, 1, 0)), 'subtract')
        for x, z in [(17, 0), (0, -17), (-17, 0), (0, 17)]:
            tool = boolean(cylinder(2, 7, (x, 6, z), (0, -1, 0)), cone(4, 2, 2, (x, 5, z), (0, -1, 0)))
            oracle = boolean(oracle, tool, 'subtract')
        oracle = blend(oracle, [(12, 5, 0)], edit or 2)
    else:
        raise ValueError(f'no independent definition for {name}')
    if not BRepCheck_Analyzer(oracle).IsValid():
        raise ValueError('independent construction has invalid topology')
    reference = props(oracle)
    result = {'binding': 'OCP', 'version': OCP.__version__, 'oracle': reference,
              'oracle_source': 'independent OpenCASCADE construction from fixed target dimensions',
              'numeric_comparison': {'relative': 1e-8, 'absolute': 1e-6,
                 'scope': 'diagnostic double-precision comparison, not a kernel acceptance tolerance'}}
    if step_path is None:
        result['step'] = {'status': 'not-run', 'reason': 'no exported STEP'}
        return result
    reader = STEPControl_Reader()
    read_ok = reader.ReadFile(str(step_path)) == IFSelect_RetDone
    roots = reader.TransferRoots() if read_ok else 0
    if not read_ok or roots == 0:
        result['step'] = {'status': 'failed', 'opened': False, 'roots': roots, 'reason': 'read or root transfer failed'}
        return result
    actual = reader.OneShape()
    valid = not actual.IsNull() and BRepCheck_Analyzer(actual).IsValid()
    faces = len(entries(actual, TopAbs_FACE))
    result['step'] = {'status': 'passed' if valid and faces == expected_faces else 'failed',
                      'opened': True, 'roots': roots, 'valid': valid, 'faces': faces,
                      'expected_faces': expected_faces, 'face_count_matches': faces == expected_faces}
    # Fixed off-boundary sample grid spans the independent part's bounds.
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    bounds = Bnd_Box()
    BRepBndLib.Add_s(oracle, bounds)
    lower, upper = bounds.CornerMin(), bounds.CornerMax()
    lohi = (lower.X(), lower.Y(), lower.Z(), upper.X(), upper.Y(), upper.Z())
    axes = [[lohi[i] + (lohi[i+3]-lohi[i]) * k / 3 for k in range(-1, 5)] for i in range(3)]
    expected_solids = [BRepClass3d_SolidClassifier(s) for s in entries(oracle, TopAbs_SOLID)]
    actual_solids = [BRepClass3d_SolidClassifier(s) for s in entries(actual, TopAbs_SOLID)]
    def member(solids, point):
        states = []
        for classifier in solids:
            classifier.Perform(gp_Pnt(*point), 1e-9)
            states.append(classifier.State())
        return None if TopAbs_ON in states else TopAbs_IN in states
    checked, boundary, wrong = 0, 0, []
    for x in axes[0]:
        for y in axes[1]:
            for z in axes[2]:
                point = (x, y, z)
                expected, got = member(expected_solids, point), member(actual_solids, point)
                if expected is None or got is None:
                    boundary += 1
                    continue
                checked += 1
                if expected != got:
                    wrong.append({'point': point, 'expected': expected, 'actual': got})
    result['membership'] = {'status': 'failed' if wrong else 'passed', 'checked': checked,
                            'boundary_unresolved': boundary, 'wrong': wrong,
                            'coverage': 'fixed 216-point grid, not a global set equality proof'}
    return result

if __name__ == '__main__':
    try:
        case = json.loads(Path(sys.argv[1]).read_text())
        value = check(case, sys.argv[2], None if sys.argv[3] == '-' else Path(sys.argv[3]), int(sys.argv[4]))
    except ImportError as error:
        value = {'status': 'not-run', 'reason': 'OpenCASCADE binding unavailable', 'detail': str(error)}
    except Exception as error:
        value = {'status': 'failed', 'reason': 'independent checker failed', 'detail': str(error)}
    print(json.dumps(value, sort_keys=True))
