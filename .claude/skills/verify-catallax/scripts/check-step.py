#!/usr/bin/env python3
"""Open one exported STEP file with OpenCASCADE, without repairing its topology."""
import json
import sys


def check(path, expected_faces):
    try:
        import OCP
        from OCP.STEPControl import STEPControl_Reader
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        from OCP.TopAbs import TopAbs_FACE
        from OCP.TopExp import TopExp_Explorer
    except ImportError as error:
        return {'status': 'not-run', 'code': 'step.binding-unavailable', 'detail': str(error)}
    reader = STEPControl_Reader()
    if reader.ReadFile(path) != IFSelect_RetDone:
        return {'status': 'failed', 'code': 'step.read-failed'}
    roots = reader.TransferRoots()
    if roots == 0:
        return {'status': 'failed', 'code': 'step.no-roots'}
    shape = reader.OneShape()
    if shape.IsNull():
        return {'status': 'failed', 'code': 'step.null-shape'}
    faces = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = explorer.Current()
        if not any(face.IsSame(previous) for previous in faces):
            faces.append(face)
        explorer.Next()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    valid = BRepCheck_Analyzer(shape).IsValid()
    same_faces = len(faces) == expected_faces
    return {'status': 'passed' if valid and same_faces else 'failed',
            'code': None if valid and same_faces else 'step.invalid-topology' if not valid else 'step.face-count',
            'binding': 'OCP', 'version': OCP.__version__, 'valid': valid,
            'faces': len(faces), 'expected_faces': expected_faces,
            'volume_mm3': props.Mass(),
            'volume_comparison': {'status': 'not-run', 'code': 'step.kernel-volume-unmeasured'}}


if __name__ == '__main__':
    try:
        result = check(sys.argv[1], int(sys.argv[2]))
    except Exception as error:
        result = {'status': 'failed', 'code': 'step.check-error', 'detail': str(error)}
    print(json.dumps(result, sort_keys=True))
