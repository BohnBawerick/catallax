"""Check the two numerical loft exports raw, including their physical holes.

Uses the repository's healing-off StepToTopoDS helper with OCP 8. It does not
call TransferRoots, ShapeFix, sewing or SameParameter before validation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import OCP
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Pnt
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_ON,
                       TopAbs_OUT, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE)
from test_loopless_faces import raw_solid, unique


def topology(shape):
    faces = unique(shape, TopAbs_FACE)
    return {'solids': len(unique(shape, TopAbs_SOLID)), 'faces': len(faces),
            'edges': len(unique(shape, TopAbs_EDGE)),
            'vertices': len(unique(shape, TopAbs_VERTEX)),
            'two_loop_faces': sum(len(unique(face, TopAbs_WIRE)) == 2 for face in faces)}


def volume(shape):
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return props.Mass()


def state(shape, point):
    result = BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7).State()
    return {TopAbs_IN: 'inside', TopAbs_OUT: 'outside', TopAbs_ON: 'boundary'}.get(result, 'unknown')


def check(root, name, counts, expected_volume, material, void, cut_x):
    out = root / name
    response = json.loads((out / 'response.json').read_text())
    assert response['status'] == 'numerical' and response['category'] == 'approximate'
    assert response['exact'] is False
    assert response['counts'] == counts and response['two_loop_faces'] == 1
    assert response['work']['limit'] == response['work']['target'] == 20_000_000
    assert not response['work']['exhausted'] and response['work']['used'] <= 20_000_000
    assert sum(response['work_stages'].values()) == response['work']['used']
    assert response['source_evidence']['tolerance_mm'] == '1/10000000'
    assert response['emitted_evidence']['tolerance_mm'] == '1/10000000'
    assert response['complete_body_bound'] == 'unestablished'
    path = out / 'answer.step'
    before = path.read_bytes()
    assert b'catallax numerical approximate' in before
    assert b'complete_body_bound=unestablished' in before
    assert b'UNCERTAINTY_MEASURE_WITH_UNIT' not in before
    shape = raw_solid(path)
    assert BRepCheck_Analyzer(shape).IsValid(), (name, 'raw import invalid')
    actual = topology(shape)
    assert actual == dict(solids=1, faces=counts[0], edges=counts[1],
                          vertices=counts[2], two_loop_faces=1), (name, actual)
    before_volume = volume(shape)
    assert abs(before_volume - expected_volume) < 1e-7, (name, before_volume)
    assert state(shape, material) == 'inside', (name, 'material lost')
    assert state(shape, void) == 'outside', (name, 'physical hole filled')
    edit = BRepAlgoAPI_Cut(shape, BRepPrimAPI_MakeBox(
        gp_Pnt(-1., -1., -1.), cut_x + 1., 10., 5.).Shape())
    assert edit.IsDone() and BRepCheck_Analyzer(edit.Shape()).IsValid(), (name, 'later cut invalid')
    after_volume = volume(edit.Shape())
    assert 0 < after_volume < before_volume, (name, after_volume)
    assert state(edit.Shape(), material) == 'inside', (name, 'later cut lost material')
    assert state(edit.Shape(), void) == 'outside', (name, 'later cut filled physical hole')
    assert path.read_bytes() == before
    return {'case': name, 'sha256': hashlib.sha256(before).hexdigest(),
            'ocp': OCP.__version__, 'healing': False, 'raw_valid': True,
            'topology': actual, 'volume_mm3': before_volume,
            'material': 'inside', 'physical_void': 'outside',
            'later_cut_valid': True, 'later_cut_volume_mm3': after_volume,
            'later_cut_material': 'inside', 'later_cut_physical_void': 'outside',
            'input_bytes_unchanged': True, 'work': response['work'],
            'work_stages': response['work_stages']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    cases = [check(args.root, 'tiny', [7, 12, 8], 16.,
                   [2., 1., 1.5], [2., 1.9999995, 1.5], 1.),
             check(args.root, 'm10-008', [9, 18, 12], 43. / 3.,
                   [2., .5, 1.5], [2., 1.5, 1.5], .5)]
    controls = []
    for name in ['reverse-hole', 'missing-hole']:
        directory = args.root / name
        response = json.loads((directory / 'response.json').read_text())
        assert response['status'] == 'refused' and not response['work']['exhausted']
        assert not (directory / 'answer.step').exists(), (name, 'stale STEP remains')
        controls.append({'control': name, 'refused': True,
                         'no_stale_step': True, 'reason': response['reason'],
                         'work': response['work']})
    report = {'cases': cases, 'controls': controls}
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
