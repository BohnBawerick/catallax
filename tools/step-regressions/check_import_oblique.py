"""Check the two fixed real workflows and authored oblique STEP controls raw.

Requires the existing OCCT Python dependency. Uses N03's healing-off raw builder;
does not call TransferRoots, ShapeFix, sewing or SameParameter. Reports metadata
only; real-file coordinates and volume values are not published.
"""
import argparse
import hashlib
import json
from pathlib import Path

import OCP
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from check_numerical_step import raw_solid, state, topology


def checked(path, expected_hash, expected_topology):
    before = path.read_bytes()
    digest = hashlib.sha256(before).hexdigest()
    assert digest == expected_hash, (path.name, 'file hash differs')
    assert b'catallax numerical approximate' in before
    assert b'complete_body_bound=unestablished' in before
    assert b'UNCERTAINTY_MEASURE_WITH_UNIT' not in before
    shape = raw_solid(path)
    assert BRepCheck_Analyzer(shape).IsValid(), (path.name, 'raw import invalid')
    assert topology(shape) == expected_topology, (path.name, 'topology differs')
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    centre = props.CentreOfMass()
    assert props.Mass() > 0, (path.name, 'nonpositive volume')
    assert path.read_bytes() == before, (path.name, 'input changed')
    return shape, props.Mass(), [centre.X(), centre.Y(), centre.Z()]


def check_real(root, expected):
    key = expected['model_id']
    response = json.loads((root / key / 'response.json').read_text())
    assert response['status'] == 'numerical'
    assert response['source_sha256'] == expected['source_sha256']
    record = response['source_record']
    assert record['source']['sha256'] == expected['source_sha256']
    assert record['derivative']['sha256'] == expected['repaired_definition_sha256']
    assert len(response['solids']) == 1
    solid = response['solids'][0]
    assert solid['source_step']['source_history_stages'] == 0
    assert solid['edited_step']['source_history_stages'] == 1
    a, av, ac = checked(root / key / 'solid-0.step',
                        expected['source_step']['sha256'],
                        expected['source_step']['topology'])
    b, bv, bc = checked(root / key / 'solid-0-edited.step',
                        expected['edited_step']['sha256'],
                        expected['edited_step']['topology'])
    assert state(a, ac) == 'inside' and state(b, bc) == 'inside'
    assert abs(av - bv) <= 1e-9 * max(1., av), (key, 'pose volume differs')
    assert all(abs(y - x - shift) < 1e-8
               for x, y, shift in zip(ac, bc, [3., 5., 7.])), (key, 'pose differs')
    assert response['work_used'] <= response['work_limit'] == 20_000_000
    return {'model_id': key, 'source_sha256': expected['source_sha256'],
            'repaired_definition_sha256': expected['repaired_definition_sha256'],
            'same_source_and_repaired_definition': True,
            'raw_source_and_edited_valid': True, 'pose_checked': True,
            'source_history_stages': 1, 'input_bytes_unchanged': True}


def check_control(root, expected):
    three = expected['control'] == 'three_component_tube'
    filename = 'oblique-3d-tube.step' if three else 'oblique-tube.step'
    shape, volume, _ = checked(root / 'oblique' / filename,
                               expected['file_sha256'], expected['topology'])
    assert abs(volume - 12.) < 1e-9
    probes = [((.75, 2., .5), 'inside'), ((2., 2., .5), 'outside'),
              ((-.25, 2., .5), 'outside')]
    for (x, y, z), answer in probes:
        a, b, c = .28 * x + .96 * z, y, -.96 * x + .28 * z
        point = (.6 * a - .8 * b, .8 * a + .6 * b, c) if three else (a, b, c)
        assert state(shape, point) == answer, (filename, 'material probe differs')
    return {'control': expected['control'], 'sha256': expected['file_sha256'],
            'raw_import_valid': True, 'material_probes_checked': True,
            'input_bytes_unchanged': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--real-results', type=Path, required=True,
                        help='model-ID directories produced by numerical_import')
    parser.add_argument('--controls', type=Path, required=True,
                        help='capture root containing oblique/ and n03-axis/')
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[2]
                        / 'docs/phase1/import-oblique-step/real-controls.json')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    expected = json.loads(args.manifest.read_text())
    assert [row['model_id'] for row in expected['rows']] == ['00000022', '00000134']
    rows = [check_real(args.real_results, row) for row in expected['rows']]
    controls = [check_control(args.controls, row) for row in expected['public_controls']]
    for row in expected['unchanged_planar_controls']:
        path = args.controls / 'n03-axis' / (row['control'] + '.step')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256']
    report = {'format': 'catallax-oblique-step-raw-check/1',
              'ocp_version': OCP.__version__, 'healing': False,
              'same_parameter': False, 'builder_precision_mm': '1e-7',
              'builder_max_tolerance_mm': '1e-7', 'real_workflows': rows,
              'public_controls': controls, 'planar_file_hashes_unchanged': True}
    text = json.dumps(report, indent=2) + '\n'
    if args.out:
        args.out.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
