#!/usr/bin/env python3
"""Check analytic values and membership; OpenCASCADE is only the STEP consumer."""
from itertools import chain, product
import json
from pathlib import Path
import sys

from analytic import definition


def healing_off_read(step_path):
    """Report the raw StepToTopoDS read at 1e-7 mm without ShapeFix, sewing or SameParameter."""
    try:
        helpers = str(Path(__file__).resolve().parents[5] / 'tools/step-regressions')
        if helpers not in sys.path:
            sys.path.insert(0, helpers)
        from check_numerical_step import raw_solid
        shape = raw_solid(step_path)
        from OCP.BRepCheck import BRepCheck_Analyzer
        return {'opened': True, 'valid': BRepCheck_Analyzer(shape).IsValid()}
    except Exception as error:
        return {'opened': False, 'valid': False, 'reason': f'{type(error).__name__}: {error}'}


def check(case, variant, step_path, expected_faces):
    reference, bounds, implicit = definition(case, variant)
    result = {'oracle': reference, 'oracle_source': 'independent elementary formulas and analytic solid predicates; no CAD-kernel reference construction',
              'numeric_comparison': {'relative': 1e-8, 'absolute': 1e-6,
                 'scope': 'diagnostic double-precision comparison, not a kernel acceptance tolerance'}}
    if reference is None:
        result['property_coverage'] = {'status': 'not-run', 'reason': 'no independent analytic property formula implemented for this complete target'}
    if step_path is None:
        result['step'] = {'status': 'not-run', 'reason': 'no exported STEP'}
        return result
    try:
        import OCP
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.gp import gp_Dir, gp_Lin, gp_Pnt
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
        from OCP.Precision import Precision
        from OCP.STEPControl import STEPControl_Reader
        from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_IN, TopAbs_ON, TopAbs_OUT
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
    # A valid open shell is not a part: require at least one solid, every shell
    # of every solid closed, and no face outside a solid. Separate solids are
    # allowed because line contact and near misses legitimately yield two.
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_SHELL
    solid_list = entries(actual, TopAbs_SOLID)
    faces_by_solid = [entries(solid, TopAbs_FACE) for solid in solid_list]
    # Union of distinct faces owned by any solid; a face shared by two solids counts once.
    owned = []
    for solid_faces in faces_by_solid:
        owned += [face for face in solid_faces if not any(face.IsSame(o) for o in owned)]
    closed = (bool(solid_list) and len(owned) == faces
              and all(BRep_Tool.IsClosed_s(shell) for solid in solid_list
                      for shell in entries(solid, TopAbs_SHELL)))
    result['step'] = {'status': 'passed' if valid and closed and faces == expected_faces else 'failed',
                      'opened': True, 'roots': roots, 'valid': valid, 'closed_solid': closed,
                      'solids': len(solid_list), 'faces': faces,
                      'expected_faces': expected_faces, 'face_count_matches': faces == expected_faces}
    step_read = case.get('step_read')
    if step_read in ('healing-off', 'healing-off-required'):
        raw = healing_off_read(step_path)
        result['step']['healing_off'] = raw
        if step_read == 'healing-off-required' and not (raw['opened'] and raw['valid']):
            result['step'].update(status='failed', reason='raw read invalid' if raw['opened'] else 'raw read failed')
    if 'census' in case and closed:
        # Material removal is judged from each kernel's exported solid, measured in OCCT.
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(actual, props)
        result['material'] = {'volume_mm3': abs(props.Mass()),
                              'source_volume_mm3': case['census']['cutter']['source_volume_mm3']}
    if implicit is None:
        result['membership'] = {'status': 'not-run', 'reason': 'no independent analytic predicate for this complete target (P5 trimmed blend, census parts, phase-2 targets)'}
        return result
    # Same 6^3 sampling scheme as the original report, with analytic target bounds.
    axes = [[bounds[i]+(bounds[i+3]-bounds[i])*k/3 for k in range(-1, 5)] for i in range(3)]
    probes = case.get('membership_probes_mm', [])
    tolerance = 1e-9
    solids = []
    for solid, solid_faces in zip(solid_list, faces_by_solid):
        face_boxes = []
        if valid and closed:
            for face in solid_faces:
                box = Bnd_Box()
                BRepBndLib.Add_s(face, box, False)
                face_boxes.append(box)
        solids.append((solid, BRepClass3d_SolidClassifier(solid), face_boxes))
    checked, boundary, wrong = 0, 0, []
    for point in chain(product(*axes), map(tuple, probes)):
        distance = implicit(point)
        origin = gp_Pnt(*point)
        upward_ray = Bnd_Box()
        upward_ray.Update(*point)
        upward_ray.OpenZmax()
        states = []
        for solid, classifier, face_boxes in solids:
            classifier.Perform(origin, tolerance)
            state = classifier.State()
            # Use face boxes to select candidates, then check the trimmed faces.
            if state == TopAbs_IN and face_boxes and all(
                    not box.IsVoid() and box.IsOut(upward_ray) for box in face_boxes):
                intersections = IntCurvesFace_ShapeIntersector()
                intersections.Load(solid, tolerance)
                intersections.Perform(gp_Lin(origin, gp_Dir(0, 0, 1)), 0., Precision.Infinite_s())
                if intersections.IsDone() and intersections.NbPnt() == 0:
                    state = TopAbs_OUT
            states.append(state)
        if abs(distance) <= tolerance or TopAbs_ON in states:
            boundary += 1
            continue
        expected, got = distance < 0, TopAbs_IN in states
        checked += 1
        if expected != got:
            wrong.append({'point': point, 'expected': expected, 'actual': got})
    coverage = '216-point grid'
    if probes:
        coverage += f' and {len(probes)} case-specific points'
    result['membership'] = {'status': 'failed' if wrong else 'passed', 'checked': checked,
                            'boundary_unresolved': boundary, 'wrong': wrong,
                            'coverage': f'{coverage} against analytic predicates, not a global set equality proof'}
    return result


if __name__ == '__main__':
    try:
        case = json.loads(Path(sys.argv[1]).read_text())
        value = check(case, sys.argv[2], None if sys.argv[3] == '-' else Path(sys.argv[3]), int(sys.argv[4]))
    except Exception as error:
        value = {'status': 'failed', 'reason': 'independent checker failed', 'detail': str(error)}
    print(json.dumps(value, sort_keys=True, allow_nan=False))
