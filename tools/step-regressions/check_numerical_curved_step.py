"""Check the two numerical loft exports raw, including their physical holes.

Uses the repository's healing-off StepToTopoDS helper with OCP 8. It does not
call TransferRoots, ShapeFix, sewing or SameParameter before validation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import OCP
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GeomAbs import GeomAbs_BSplineCurve
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Pnt
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_ON,
                       TopAbs_OUT, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE)
from OCP.TopoDS import TopoDS
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


def check_native_uv(root):
    result = check(root, 'native', [9, 18, 12], 43. / 3.,
                   [2., .5, 1.5], [2., 1.5, 1.5], .5)
    response = json.loads((root / 'native/response.json').read_text())
    bindings = response['native_trims']
    assert bindings['original'] == bindings['source']
    assert len(bindings['source']) == len(bindings['emitted']) == 2
    for source, emitted in zip(bindings['source'], bindings['emitted']):
        for key in ['face', 'boundary', 'use_index', 'edge', 'degree']:
            assert source[key] == emitted[key]
        assert source['knots'][0] == '-1' and source['knots'][-1] == '2'
        assert emitted['knots'][0] == '0' and emitted['knots'][-1] == '1'
    shape = raw_solid(root / 'native/answer.step')
    residual, mapping_error, uses = 0., 0., 0
    for item in unique(shape, TopAbs_FACE):
        face = TopoDS.Face(item)
        surface = BRepAdaptor_Surface(face)
        for item in unique(face, TopAbs_EDGE):
            edge = TopoDS.Edge(item)
            curve = BRepAdaptor_Curve(edge)
            if curve.GetType() != GeomAbs_BSplineCurve or curve.Degree() != 3:
                continue
            uses += 1
            lo, hi = curve.FirstParameter(), curve.LastParameter()
            assert abs(lo) < 1e-12 and abs(hi - 1.) < 1e-12
            pcurve = BRep_Tool.CurveOnSurface_s(edge, face, lo, hi)
            assert pcurve is not None, 'native pcurve missing after raw import'
            first, last = curve.Value(lo), curve.Value(hi)
            for i in range(33):
                t = i / 32.
                point = curve.Value(t)
                share = t / 2. + 1.5 * t * t - t * t * t
                expected = gp_Pnt(*(a + share * (b - a) for a, b in
                                    zip(first.Coord(), last.Coord())))
                mapping_error = max(mapping_error, point.Distance(expected))
                uv = pcurve.Value(t)
                residual = max(residual, point.Distance(surface.Value(uv.X(), uv.Y())))
    assert uses == 2 and residual <= 1e-7 and mapping_error <= 1e-7
    controls = []
    for name in ['missing', 'reversed', 'partial', 'wrong', 'stale', 'duplicate']:
        out = root / name
        response = json.loads((out / 'response.json').read_text())
        assert response['status'] == 'refused' and not response['work']['exhausted']
        assert not (out / 'answer.step').exists(), (name, 'stale STEP remains')
        controls.append({'control': name, 'refused': True, 'no_stale_step': True,
                         'reason': response['reason'], 'work': response['work']})
    result['native_uv'] = {'uses': uses, 'sampled_same_parameter_residual_mm': residual,
                           'sampled_cubic_mapping_error_mm': mapping_error,
                           'controls': controls}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--repeated-endpoint', type=Path)
    parser.add_argument('--native-uv-root', type=Path)
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
    if args.native_uv_root:
        report['native_uv'] = check_native_uv(args.native_uv_root)
    if args.repeated_endpoint:
        shape = raw_solid(args.repeated_endpoint)
        assert BRepCheck_Analyzer(shape).IsValid(), 'repeated-endpoint raw import invalid'
        assert topology(shape) == dict(solids=1, faces=10, edges=24, vertices=16, two_loop_faces=2)
        assert abs(volume(shape) - 12.) < 1e-7, 'unadmitted traversal changed material'
        for edge in unique(shape, TopAbs_EDGE):
            curve = BRepAdaptor_Curve(TopoDS.Edge(edge))
            lo, hi = curve.FirstParameter(), curve.LastParameter()
            for i in range(9):
                point = curve.Value(lo + (hi - lo) * i / 8)
                assert -1e-7 <= point.X() <= 4. + 1e-7
                assert -1e-7 <= point.Y() <= 4. + 1e-7
                assert -1e-7 <= point.Z() <= 1. + 1e-7, 'unadmitted off-plane carrier span used'
        report['repeated_endpoint'] = {'raw_valid': True, 'volume_mm3': volume(shape),
                                       'finite_edge_samples_in_body': True, 'healing': False}
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
