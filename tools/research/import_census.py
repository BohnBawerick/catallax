"""Rerun the fixed 3,000-source R8 sample through opt-in numerical admission.

Use a release r8_census binary with --numerical support. Raw diagnostics and
archive member names remain private. A complete report requires one terminal
observation for every original model ID, with the pinned individual byte hashes.
External observation timeouts remain a separate outcome.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess

from r8_step_census import qualify_coverage

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261111
COUNT = 3000
HASH_MANIFEST = '3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7'
PRIOR_CSV = '3da02551580d917b4b15475222abb453d7bf23f8a1beb25c76248df54a38ed9a'
SELECTION_SHA = '1b11410b220d5e8a3156794f6c8f9cc09091a43dcac9af801215fb3a38d571a5'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def baseline_outcome(row):
    if row['outcome'] == 'admitted':
        return 'full' if row['coverage_complete'] == 'True' else 'partial'
    if row['outcome'] == 'observation_timeout':
        return 'timeout'
    assert row['outcome'] == 'refused'
    return 'refused'


def observe(binary, path, key, expected_hash, raw, timeout, baseline):
    assert digest(path) == expected_hash, (key, 'source bytes differ')
    stdout, stderr = raw / (key + '.jsonl'), raw / (key + '.stderr')
    with stdout.open('wb') as out, stderr.open('wb') as err:
        try:
            process = subprocess.run([str(binary), '--numerical', str(path)],
                                     stdout=out, stderr=err, timeout=timeout, check=False,
                                     env={**os.environ, 'CATALLAX_THREADS': '1'})
            outcome = 'process_error' if process.returncode else 'missing_result'
        except subprocess.TimeoutExpired:
            outcome = 'timeout'
    assert digest(path) == expected_hash, (key, 'source changed during observation')
    checkpoints = []
    for line in stdout.read_text().splitlines():
        try:
            checkpoints.append(json.loads(line))
        except json.JSONDecodeError:
            # Killing the observer can interrupt the current checkpoint write.
            if outcome != 'timeout':
                raise
    sources = [row for row in checkpoints if row.get('stage') == 'faces']
    imports = [row for row in checkpoints if row.get('stage') == 'import']
    assert len(sources) <= 1 and len(imports) <= 1, (key, 'duplicate checkpoint')
    source = qualify_coverage(sources[0]) if sources else {}
    imported = imports[0] if imports else {}
    if imports:
        assert imported['work_limit'] == 20_000_000
    if outcome == 'missing_result' and imports:
        if imported['outcome'] == 'admitted':
            assert imported['bodies'] > 0
            assert imported['import_category'] in ('certified', 'numerical approximate')
            assert imported['work_used'] <= imported['work_limit']
            if imported['import_category'] == 'numerical approximate':
                assert imported['exact'] is False
                assert imported['complete_body_bound'] == 'unestablished'
                assert set(imported['tolerances_mm']) == {'1/100'}
                assert Fraction(imported['max_checked_residual_mm']) <= Fraction(1, 100)
            outcome = 'full' if source.get('coverage_complete') is True else 'partial'
        else:
            assert imported['outcome'] == 'refused'
            outcome = 'refused'
    return {
        'model_id': key, 'source_sha256': expected_hash,
        'baseline_outcome': baseline_outcome(baseline),
        'baseline_first_refusal': baseline['refusal_kind'],
        'outcome': outcome, 'after_fix_cause': imported.get('refusal_kind', ''),
        'import_category': imported.get('import_category', ''),
        'coverage_complete': source.get('coverage_complete') is True,
        'source_parse': source.get('source_parse', 'unobserved'),
        'bodies': imported.get('bodies'),
        'work_used': imported.get('work_used'),
        'work_limit': imported.get('work_limit'),
        'tolerances_mm': imported.get('tolerances_mm', []),
        'source_to_result_max_mm': imported.get('source_to_result_max_mm'),
        'source_change_count': imported.get('source_change_count'),
        'source_admission_method': imported.get('source_admission_method'),
        'source_pair_disjoint_count': (None if imported.get('source_intersection_checks') is None else
                                      sum(check['status'] == 'disjoint' for check in
                                          imported['source_intersection_checks'])),
        'source_pair_unresolved': (None if imported.get('source_intersection_checks') is None else
                                  [{'solid': check['solid'], 'faces': check['faces'],
                                   'code': check['refusal']['code']} for check in
                                  imported['source_intersection_checks']
                                  if check['status'] == 'unresolved']),
        'max_checked_residual_mm': imported.get('max_checked_residual_mm'),
        'referenced_faces': source.get('referenced_faces'),
        'face_kinds': source.get('face_kinds', {}),
        'solid_entities': source.get('solid_entities'),
        'unknown_representation_items': source.get('unknown_representation_items')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--build-record', type=Path, required=True,
                        help='actual clean-checkout Cargo artifact record from build_import_observer.py')
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--selection', type=Path, required=True,
                        help='private original [model_id, archive_member] selection')
    parser.add_argument('--source-hashes', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True,
                        help='private new capture directory, or existing with --resume')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--timeout', type=int, default=60)
    parser.add_argument('--baseline', type=Path,
                        default=ROOT / 'docs/milestones/R8-public-step-census.csv')
    parser.add_argument('--baseline-manifest', type=Path,
                        default=ROOT / 'docs/milestones/R8-public-step-census.json')
    parser.add_argument('--prior', type=Path,
                        default=ROOT / 'docs/research/step-census/results.csv',
                        help='the complete 9.20 percent checkpoint, with identical inputs')
    args = parser.parse_args()
    assert args.workers > 0 and args.timeout == 60
    baseline_manifest = json.loads(args.baseline_manifest.read_text())
    assert baseline_manifest['seed'] == SEED
    assert baseline_manifest['evaluated_model_ids'] == COUNT
    assert digest(args.source_hashes) == HASH_MANIFEST
    assert digest(args.prior) == PRIOR_CSV
    assert digest(args.selection) == SELECTION_SHA
    hashes = json.loads(args.source_hashes.read_text())
    selection = json.loads(args.selection.read_text())
    baseline_rows = list(csv.DictReader(args.baseline.open()))
    baseline = {row['model_id']: row for row in baseline_rows}
    prior_rows = list(csv.DictReader(args.prior.open()))
    prior = {row['model_id']: row for row in prior_rows}
    assert len(baseline_rows) == len(baseline) == len(selection) == len(hashes) == COUNT
    assert len(prior_rows) == len(prior) == COUNT and set(prior) == set(hashes)
    assert all(prior[key]['source_sha256'] == value for key, value in hashes.items())
    assert [key for key, _ in selection] == sorted(baseline) == sorted(hashes)
    files = {}
    for key, member in selection:
        path = PurePosixPath(member)
        assert not path.is_absolute() and '..' not in path.parts
        candidate = (args.sources / member).resolve()
        assert candidate.is_relative_to(args.sources.resolve())
        files[key] = candidate
    subprocess.run(['git', '-C', str(ROOT), 'diff', '--quiet', 'HEAD'], check=True)
    source_revision = subprocess.check_output(
        ['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    build_record = json.loads(args.build_record.read_text())
    assert build_record['format'] == 'catallax-census-observer-build/1'
    assert build_record['source_revision'] == source_revision
    assert build_record['binary_sha256'] == digest(args.binary)
    assert Path(build_record['artifact_path']).resolve() == args.binary.resolve()
    assert build_record['observer_sha256'] == digest(ROOT / 'crates/catallax-import/examples/r8_census.rs')
    assert build_record['cargo_lock_sha256'] == digest(ROOT / 'Cargo.lock')
    assert build_record['git_tree'] == subprocess.check_output(
        ['git', '-C', str(ROOT), 'rev-parse', 'HEAD^{tree}'], text=True).strip()
    identity = {
        'schema': 'catallax-numerical-import-census/2', 'seed': SEED,
        'denominator': COUNT, 'binary_sha256': digest(args.binary),
        'collector_sha256': digest(Path(__file__)),
        'source_revision': source_revision,
        'build_record': build_record, 'build_record_sha256': digest(args.build_record),
        'source_root': str(args.sources.resolve()),
        'observer_sha256': digest(ROOT / 'crates/catallax-import/examples/r8_census.rs'),
        'rustc': subprocess.check_output(['rustc', '--version'], text=True).strip(),
        'source_hash_manifest_sha256': HASH_MANIFEST,
        'selection_sha256': digest(args.selection),
        'baseline_csv_sha256': digest(args.baseline),
        'prior_csv_sha256': digest(args.prior),
        'archive_sha256': baseline_manifest['archive_sha256'],
        'archive_md5': baseline_manifest['archive_md5'],
        'archive_bytes': baseline_manifest['archive_bytes'],
        'observation_timeout_seconds': args.timeout, 'workers': args.workers,
        'per_import_workers': 1,
        'per_import_work_limit': 20_000_000,
        'verification_byte_work': 'charged before the initial kernel verification hash',
        'full_definition': 'admitted with complete source coverage',
        'partial_definition': 'admitted with unverified source coverage',
        'scope': 'import admission; does not assert an edit or STEP export for every admission'}
    rows = {}
    if args.resume:
        assert json.loads((args.out / 'identity.json').read_text()) == identity
        for line in (args.out / 'rows.jsonl').read_text().splitlines():
            row = json.loads(line)
            key = row['model_id']
            assert key not in rows and row['source_sha256'] == hashes[key]
            rows[key] = row
    else:
        args.out.mkdir()
        (args.out / 'raw').mkdir()
        save(args.out / 'identity.json', identity)
    with (args.out / 'rows.jsonl').open('a') as checkpoint:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            pending = {pool.submit(observe, args.binary, files[key], key, hashes[key],
                                   args.out / 'raw', args.timeout, baseline[key]): key
                       for key in sorted(files) if key not in rows}
            for future in as_completed(pending):
                row = future.result()
                key = row['model_id']
                row['prior_outcome'] = prior[key]['outcome']
                row['prior_first_cause'] = prior[key]['after_fix_cause']
                rows[key] = row
                checkpoint.write(json.dumps(row, sort_keys=True) + '\n')
                checkpoint.flush()
                if len(rows) % 25 == 0:
                    print(json.dumps({'observed': len(rows), 'denominator': COUNT,
                                      'outcomes': dict(Counter(r['outcome'] for r in rows.values()))}),
                          flush=True)
    assert digest(args.binary) == identity['binary_sha256'], 'binary changed during run'
    assert set(rows) == set(baseline)
    ordered = [rows[key] for key in sorted(rows)]
    with (args.out / 'census.csv').open('w', newline='') as table:
        writer = csv.DictWriter(table, fieldnames=list(ordered[0]))
        writer.writeheader()
        for row in ordered:
            writer.writerow({key: json.dumps(value, sort_keys=True) if isinstance(value, (list, dict))
                             else value for key, value in row.items()})
    outcomes = Counter(row['outcome'] for row in ordered)
    complete = set(outcomes).issubset({'full', 'partial', 'refused', 'timeout'})
    report = {**identity, 'complete': complete, 'evaluated': len(rows),
              'per_model_csv_sha256': digest(args.out / 'census.csv'),
              'baseline_outcomes': dict(Counter(row['baseline_outcome'] for row in ordered)),
              'outcomes': dict(outcomes),
              'baseline_first_refusals': dict(Counter(row['baseline_first_refusal'] for row in ordered
                                                     if row['baseline_outcome'] == 'refused')),
              'after_fix_causes': dict(Counter(row['after_fix_cause'] for row in ordered
                                              if row['outcome'] == 'refused')),
              'import_categories': dict(Counter(row['import_category'] for row in ordered
                                               if row['outcome'] in ('full', 'partial'))),
              'outcome_transitions': dict(Counter(row['baseline_outcome'] + ' -> ' + row['outcome']
                                                 for row in ordered)),
              'prior_outcomes': dict(Counter(row['prior_outcome'] for row in ordered)),
              'prior_first_causes': dict(Counter(row['prior_first_cause'] for row in ordered
                                                if row['prior_outcome'] == 'refused')),
              'prior_outcome_transitions': dict(Counter(row['prior_outcome'] + ' -> ' + row['outcome']
                                                       for row in ordered)),
              'recovered_since_prior_model_ids': [row['model_id'] for row in ordered
                                                  if row['outcome'] == 'full' and row['prior_outcome'] != 'full'],
              'lost_prior_full_model_ids': [row['model_id'] for row in ordered
                                            if row['prior_outcome'] == 'full' and row['outcome'] != 'full'],
              'numerical_source_methods': dict(Counter(row['source_admission_method'] for row in ordered
                                                       if row['source_admission_method'])),
              'recovered_full_model_ids': [row['model_id'] for row in ordered
                                           if row['outcome'] == 'full' and row['baseline_outcome'] != 'full'],
              'lost_full_model_ids': [row['model_id'] for row in ordered
                                      if row['baseline_outcome'] == 'full' and row['outcome'] != 'full'],
              'full_import_share_percent': f'{100 * outcomes["full"] / COUNT:.6f}'}
    save(args.out / 'summary.json', report)
    save(args.out / 'ordered-rows.json', ordered)
    print(json.dumps(report, sort_keys=True), flush=True)
    assert complete, 'process errors or missing results; not a complete census'


if __name__ == '__main__':
    main()
