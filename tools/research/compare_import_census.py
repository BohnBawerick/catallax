"""Validate and publish a same-source Catallax/standard-OCCT census comparison."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
csv.field_size_limit(64 * 1024 * 1024)
COUNT = 3000
HASHES = '3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7'
SELECTION = '1b11410b220d5e8a3156794f6c8f9cc09091a43dcac9af801215fb3a38d571a5'
OCCT_OUTCOMES = {'healed_valid_solid', 'open_shell', 'opened_other', 'failed', 'timeout'}
OCCT_FIELDS = ['outcome', 'opened', 'stage', 'all_roots_transferred', 'roots_available',
    'roots_transferred', 'result_shapes', 'transferred_whole_shape_valid',
    'healed_whole_shape_valid', 'solid_count', 'valid_solid_count', 'shell_count',
    'open_shell_count', 'free_shell_count', 'free_face_count', 'free_wire_count',
    'free_edge_count', 'free_vertex_count', 'face_occurrences']


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def indexed(rows, hashes):
    result = {row['model_id']: row for row in rows}
    assert len(result) == len(rows) == len(hashes) == COUNT, 'incomplete or duplicate model IDs'
    assert set(result) == set(hashes), 'different model IDs'
    assert all(result[key]['source_sha256'] == value for key, value in hashes.items()), 'source bytes differ'
    return result


def compare(native, native_manifest, reference, reference_manifest, hashes):
    for manifest in [native_manifest, reference_manifest]:
        assert manifest['complete'] is True and manifest['denominator'] == COUNT
        assert manifest['seed'] == 20261111
        assert manifest['source_hash_manifest_sha256'] == HASHES
        assert manifest['selection_sha256'] == SELECTION
    build = native_manifest['build_record']
    for field in ['source_revision', 'binary_sha256', 'observer_sha256']:
        assert build[field] == native_manifest[field], 'observer build record differs from capture'
    assert build['cargo_artifact']['target']['name'] == 'r8_census'
    assert build['cargo_artifact']['target']['kind'] == ['example']
    native = indexed(native, hashes)
    reference = indexed(reference, hashes)
    assert all(r['outcome'] in {'full', 'partial', 'refused', 'timeout'} for r in native.values())
    assert all(r['outcome'] in OCCT_OUTCOMES and r['source_unchanged'] is True for r in reference.values())
    for row in reference.values():
        if row['outcome'] == 'healed_valid_solid':
            assert row['opened'] is True and row['stage'] == 'complete'
            assert row['healed_whole_shape_valid'] is True
            assert row['solid_count'] > 0 and row['valid_solid_count'] == row['solid_count']
            assert all(row['free_' + name + '_count'] == 0 for name in ['shell', 'face', 'wire', 'edge', 'vertex'])
        elif row['outcome'] == 'open_shell':
            assert row['opened'] is True and row['solid_count'] == 0 and row['open_shell_count'] > 0
        elif row['outcome'] == 'opened_other':
            assert row['opened'] is True
    assert dict(Counter(r['outcome'] for r in native.values())) == native_manifest['outcomes']
    assert dict(Counter(r['outcome'] for r in reference.values())) == reference_manifest['outcomes']
    rows = []
    for key in sorted(hashes):
        n, r = native[key], reference[key]
        rows.append({'model_id': key, 'source_sha256': hashes[key],
            'catallax_outcome': n['outcome'], 'catallax_category': n['import_category'],
            'catallax_first_cause': n['after_fix_cause'],
            **{'occt_' + field: r.get(field) for field in OCCT_FIELDS}})
    report = {'format': 'catallax-occt-import-comparison/1', 'complete': True,
        'denominator': COUNT, 'seed': 20261111, 'source_hash_manifest_sha256': HASHES,
        'selection_sha256': SELECTION, 'catallax_source_revision': native_manifest['source_revision'],
        'catallax_binary_sha256': native_manifest['binary_sha256'],
        'occt_protocol_sha256': reference_manifest['protocol_sha256'],
        'ocp_versions': reference_manifest['ocp_versions'],
        'catallax_outcomes': native_manifest['outcomes'], 'occt_outcomes': reference_manifest['outcomes'],
        'occt_opened_observed': sum(bool(r.get('opened')) for r in reference.values()),
        'occt_all_roots_transferred_observed': sum(bool(r.get('all_roots_transferred'))
            and r.get('roots_available', 0) > 0 for r in reference.values()),
        'outcome_pairs': dict(Counter(row['catallax_outcome'] + ' / ' + row['occt_outcome'] for row in rows)),
        'occt_opened_catallax_not_full_model_ids': [row['model_id'] for row in rows
            if row['occt_opened'] and row['catallax_outcome'] != 'full'],
        'catallax_full_occt_opening_unobserved_model_ids': [row['model_id'] for row in rows
            if row['catallax_outcome'] == 'full' and not row['occt_opened']],
        'occt_opened_catallax_first_causes': dict(Counter(row['catallax_first_cause'] for row in rows
            if row['occt_opened'] and row['catallax_outcome'] == 'refused')),
        'limits': [
            'OCCT opened means a non-null transferred shape, not complete source coverage or validity.',
            'healed_valid_solid requires a valid whole healed shape, at least one valid solid and no free lower-dimensional topology.',
            'OCCT root transfer counts are recorded separately; they do not prove source face/item coverage.',
            'OCCT uses default STEP transfer and ShapeFix tolerances, not the Catallax 0.01 mm cap.',
            'An opened checkpoint before a timeout establishes only that opening occurred; final healing/validity may be unknown.',
            'The external 60-second deadline is comparable; OCCT has no Catallax 20M work accounting.']}
    return rows, report


def public_reference_manifest(manifest, rows):
    fields = '''format complete denominator seed selection_sha256 source_hash_manifest_sha256
        protocol_sha256 workers external_timeout_seconds reader_policy healing_policy
        catallax_work_limit_comparable rows_sha256 outcomes opened transferred_whole_shape_valid
        healed_whole_shape_valid ocp_versions'''.split()
    public = {key: manifest[key] for key in fields if key in manifest}
    for field in ['reader_integer_defaults', 'reader_real_defaults', 'shape_fix_defaults']:
        profiles = Counter(json.dumps(row[field], sort_keys=True) for row in rows if field in row)
        public[field + '_observed_profiles'] = [
            {'settings': json.loads(settings), 'files_observed': count}
            for settings, count in sorted(profiles.items())]
        public[field + '_unobserved_files'] = sum(field not in row for row in rows)
    public['validity_transitions'] = dict(Counter(
        str(row.get('transferred_whole_shape_valid')) + ' -> ' + str(row.get('healed_whole_shape_valid'))
        for row in rows))
    public['scope'] = 'Standard transferred results after default ShapeFix; not complete source-item coverage or a common tolerance bound.'
    return public


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['native_csv', 'native_manifest', 'occt_rows', 'occt_manifest', 'out']:
        parser.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    args = parser.parse_args()
    native_manifest = json.loads(args.native_manifest.read_text())
    reference_manifest = json.loads(args.occt_manifest.read_text())
    assert digest(args.native_csv) == native_manifest['per_model_csv_sha256']
    assert digest(args.occt_rows) == reference_manifest['rows_sha256']
    hashes_path = ROOT / 'docs/research/step-census/source-hashes.json'
    assert digest(hashes_path) == HASHES
    reference_rows = json.loads(args.occt_rows.read_text())
    rows, report = compare(list(csv.DictReader(args.native_csv.open())), native_manifest,
        reference_rows, reference_manifest, json.loads(hashes_path.read_text()))
    retained = json.loads((ROOT / 'docs/research/step-census/retained-controls.json').read_text())
    by_key = {row['model_id']: row for row in rows}
    assert retained['count'] == len(retained['rows']) == 10
    assert all(by_key[row['model_id']]['catallax_outcome'] == 'full'
        and by_key[row['model_id']]['source_sha256'] == row['source_sha256'] for row in retained['rows'])
    args.out.mkdir()
    table = args.out / 'comparison.csv'
    with table.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report.update(comparison_csv_sha256=digest(table), retained_controls_count=10,
        catallax_input_csv_sha256=digest(args.native_csv), occt_input_rows_sha256=digest(args.occt_rows))
    (args.out / 'comparison.json').write_text(json.dumps(report, sort_keys=True, indent=2) + '\n')
    (args.out / 'occt-reference.json').write_text(json.dumps(
        public_reference_manifest(reference_manifest, reference_rows), sort_keys=True, indent=2) + '\n')
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
