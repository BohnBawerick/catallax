#!/usr/bin/env python3
"""Check analytic values and membership; OpenCASCADE is only the STEP consumer."""
import json
from pathlib import Path
import sys

from analytic import definition


def check(case, variant, step_path, expected_faces):
    reference, bounds, implicit = definition(case, variant)
    result = {'oracle': reference, 'oracle_source': 'independent elementary formulas and analytic solid predicates; no CAD-kernel reference construction',
              'numeric_comparison': {'relative': 1e-8, 'absolute': 1e-6,
                 'scope': 'diagnostic double-precision comparison, not a kernel acceptance tolerance'}}
    if reference is None:
        result['property_coverage'] = {'status': 'not-run', 'reason': 'no analytic property formula implemented for this intersecting curved P0 family'}
    if step_path is None:
        result['step'] = {'status': 'not-run', 'reason': 'no exported STEP'}
        return result
    try:
        import OCP
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.gp import gp_Pnt
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.STEPControl import STEPControl_Reader
        from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_IN, TopAbs_ON
        from OCP.TopExp import TopExp_Explorer
    except ImportError as error:
        result['step'] = {'status': 'not-run', 'reason': 'OpenCASCADE binding unavailable', 'detail': str(error)}
        return result
    def entries(shape, kind):
        values = []
        iterator = TopExp_Explorer(shape, kind)
        while iterator.More():
            item = iterator.Current()
            if not any(item.IsSame(old) for old in values):
                values.append(item)
            iterator.Next()
        return values
    result.update(binding='OCP', version=OCP.__version__)
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
    # Same 6^3 sampling scheme as the original report, with analytic target bounds.
    axes = [[bounds[i]+(bounds[i+3]-bounds[i])*k/3 for k in range(-1, 5)] for i in range(3)]
    solids = [BRepClass3d_SolidClassifier(s) for s in entries(actual, TopAbs_SOLID)]
    checked, boundary, wrong = 0, 0, []
    for x in axes[0]:
        for y in axes[1]:
            for z in axes[2]:
                point = (x, y, z)
                distance = implicit(point)
                states = []
                for classifier in solids:
                    classifier.Perform(gp_Pnt(*point), 1e-9)
                    states.append(classifier.State())
                if abs(distance) <= 1e-9 or TopAbs_ON in states:
                    boundary += 1
                    continue
                expected, got = distance < 0, TopAbs_IN in states
                checked += 1
                if expected != got:
                    wrong.append({'point': point, 'expected': expected, 'actual': got})
    result['membership'] = {'status': 'failed' if wrong else 'passed', 'checked': checked,
                            'boundary_unresolved': boundary, 'wrong': wrong,
                            'coverage': '216-point grid against analytic predicates, not a global set equality proof'}
    return result


if __name__ == '__main__':
    try:
        case = json.loads(Path(sys.argv[1]).read_text())
        value = check(case, sys.argv[2], None if sys.argv[3] == '-' else Path(sys.argv[3]), int(sys.argv[4]))
    except Exception as error:
        value = {'status': 'failed', 'reason': 'independent checker failed', 'detail': str(error)}
    print(json.dumps(value, sort_keys=True, allow_nan=False))
