"""Check a fixed real material-edit refusal and a supported material change raw.

Requires the existing OCCT Python dependency. Does not heal or rewrite STEP.
Real coordinates and volume values remain private; only authored volumes and
reference outcomes are printed.
"""
import argparse
import json
from pathlib import Path

import OCP
from check_import_oblique import checked
from check_numerical_step import state


def check(root, manifest):
    report = json.loads((root / 'report.json').read_text())
    expected = next(row for row in manifest['rows']
                    if row['model_id'] == '00000022')
    assert report['source_sha256'] == expected['source_sha256']
    assert report['repaired_definition_sha256'] == expected['repaired_definition_sha256']
    assert report['work_used'] <= report['work_limit'] == 20_000_000
    assert report['requested_edit']['kind'] == 'resize'
    assert report['requested_edit']['axis'] == 0
    assert report['requested_edit']['expected_volume_ratio'] == '2'
    assert report['requested_edit']['expected_delta_volume_ratio'] == '1'
    real = report['real']
    assert real['status'] == 'refused'
    assert real['refusal']['code'] == 'features.numerical-unsupported'
    assert real['prior_geometry_and_history_unchanged']
    assert not real['material_delta_applied']
    assert report['protected_history_stages'] == 1
    a = real['before']
    b = real['preserved_after_refusal']
    assert a['sha256'] == b['sha256'] == expected['edited_step']['sha256']
    assert a['history_stages'] == b['history_stages'] == 1
    assert a['source_history_stages'] == b['source_history_stages'] == 1
    topology = expected['edited_step']['topology']
    before, av, ac = checked(root / a['file'], a['sha256'], topology)
    preserved, bv, bc = checked(root / b['file'], b['sha256'], topology)
    assert av == bv and ac == bc
    assert state(before, ac) == state(preserved, bc) == 'inside'

    control = report['positive_control']
    assert control['status'] == 'numerical'
    assert control['history_retained'] and control['replay_checked']
    assert control['expected_before_volume_mm3'] == '12'
    assert control['expected_after_volume_mm3'] == '24'
    assert control['expected_delta_mm3'] == '12'
    a, b = control['before'], control['after']
    assert b['history_stages'] == a['history_stages'] + 1
    assert b['source_history_stages'] == a['source_history_stages'] + 1
    topology = {'solids': 1, 'faces': 10, 'edges': 24, 'vertices': 16,
                'two_loop_faces': 2, 'euler_with_face_holes': 0}
    before, av, _ = checked(root / a['file'], a['sha256'], topology)
    after, bv, _ = checked(root / b['file'], b['sha256'], topology)
    assert abs(av - 12.) < 1e-9 and abs(bv - 24.) < 1e-9
    assert abs(bv - av - 12.) < 1e-9
    assert state(before, (.75, 2., 1.5)) == 'outside'
    assert state(after, (.75, 2., 1.5)) == 'inside'
    for shape in (before, after):
        assert state(shape, (.75, 2., .5)) == 'inside'
        assert state(shape, (2., 2., .5)) == 'outside'
    return {
        'format': 'catallax-import-material-check/1',
        'ocp_version': OCP.__version__, 'healing': False, 'same_parameter': False,
        'builder_precision_mm': '1e-7', 'builder_max_tolerance_mm': '1e-7',
        'model_id': '00000022', 'source_sha256': report['source_sha256'],
        'repaired_definition_sha256': report['repaired_definition_sha256'],
        'real_status': 'refused', 'real_refusal': real['refusal'],
        'real_prior_geometry_and_history_unchanged': True,
        'real_raw_before_and_preserved_valid': True,
        'real_material_delta_applied': False,
        'authored_control': {
            'raw_before_and_after_valid': True, 'volume_delta_mm3': '12',
            'volume_before_mm3': '12', 'volume_after_mm3': '24',
            'history_retained': True, 'replay_checked': True,
            'material_probes_checked': True,
            'before_sha256': a['sha256'], 'after_sha256': b['sha256']},
        'input_bytes_unchanged': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', required=True, type=Path)
    parser.add_argument('--manifest', type=Path,
                        default=Path(__file__).resolve().parents[2]
                        / 'docs/phase1/import-oblique-step/real-controls.json')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    report = check(args.results, json.loads(args.manifest.read_text()))
    text = json.dumps(report, indent=2) + '\n'
    if args.out:
        args.out.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
