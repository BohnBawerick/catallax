"""Phase 5 step 5.1 census driver: no request work limit, one time cap per file.

Runs the opened-import observer, the OpenCASCADE probe, or both over the frozen
selection or a shard of it, with a resumable row store. Each row keeps the
collector's fields plus wall seconds and the observer's own open and export
seconds and work. Compact rows also go to standard output, one line each,
prefixed P5ROW (Catallax) or P5OCCT, so a run whose only record is its log
still yields every finished row.
"""
import argparse
import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import opened_import_census as oic  # noqa: E402
import rss_guard  # noqa: E402

COMPACT = ['model_id', 'outcome', 'export_validity', 'timeout_stage', 'wall_seconds',
           'observer_open_seconds', 'observer_export_seconds', 'open_work_used', 'export_work_used',
           'open_work_exhausted', 'observed_body_count', 'all_solids_export_valid',
           'has_valid_exported_solid', 'any_valid_exported_solid_count',
           'any_valid_exported_solid_count_is_lower_bound', 'prototype_solids_export_valid',
           'source_document_exported', 'valid_exported_solid_count', 'invalid_exported_body_count',
           'reader_refused_body_count', 'first_refusal', 'export_refusal', 'proofs',
           'valid_solid_body_count', 'repaired_solid_body_count', 'tolerant_solid_body_count',
           'surface_body_count', 'observer_peak_rss_mb', 'cause', 'error']


def compact(row):
    out = {k: row[k] for k in COMPACT if row.get(k) is not None}
    blockers = {}
    for blocker in row.get('body_first_blockers') or []:
        key = '|'.join(str(blocker.get(k)) for k in ('stage', 'refusal', 'resource_limited'))
        blockers[key] = blockers.get(key, 0) + 1
    if blockers:
        out['blockers'] = blockers
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--selection', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--kernel', choices=['catallax', 'occt'], default='catallax')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=600)
    parser.add_argument('--raw-timeout', type=float, default=600)
    parser.add_argument('--work', type=int, default=1 << 62)
    parser.add_argument('--max-rss-mb', type=int, default=6144)
    parser.add_argument('--only', type=Path, help='JSON list of model ids')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--time-budget', type=float, default=0,
                        help='seconds after which no further file starts; 0 for none')
    args = parser.parse_args()
    csv.field_size_limit(64 * 1024 * 1024)
    with oic.BASELINE.open(newline='') as stream:
        baseline = {r['model_id']: r for r in csv.DictReader(stream)}
    hashes = json.loads((oic.ROOT / 'docs/research/step-census/source-hashes.json').read_text())
    selection = dict(json.loads(args.selection.read_text()))
    assert set(hashes) == set(baseline) and len(hashes) == 3000 and set(selection) <= set(hashes)
    keys = sorted(selection)
    if args.only:
        wanted = set(json.loads(args.only.read_text()))
        keys = [k for k in keys if k in wanted]
    if args.limit:
        keys = keys[:args.limit]
    files = {k: (args.sources / selection[k]).resolve() for k in keys}
    rows_dir = args.out / 'rows'
    rows_dir.mkdir(parents=True, exist_ok=True)
    done = {p.stem for p in rows_dir.glob('*.json')}
    todo = [k for k in keys if k not in done]
    started = time.monotonic()
    print(json.dumps(dict(kernel=args.kernel, selected=len(keys), done=len(done & set(keys)), todo=len(todo),
                          cpu_cap=rss_guard.CPU_CAP)), flush=True)

    def over_budget():
        return args.time_budget and time.monotonic() - started > args.time_budget

    def catallax(key):
        if over_budget():
            return key, None
        began = time.monotonic()
        try:
            row = oic.observe(args.binary, files[key], baseline[key], args.out / 'parts' / key, False,
                              args.max_rss_mb, timeout=args.timeout, work=args.work,
                              raw_timeout=args.raw_timeout)
        except Exception as error:  # one bad row must not end a long run
            row = dict(model_id=key, outcome='driver_error', error=type(error).__name__ + ': ' + str(error)[:500])
        row['wall_seconds'] = time.monotonic() - began
        attempts = sorted((args.out / 'parts' / key).glob('attempt-*'))
        response = attempts[-1] / 'capture/response.json' if attempts else None
        opened = attempts[-1] / 'capture/opened.json' if attempts else None
        try:
            if response is not None and response.is_file():
                value = json.loads(response.read_text())
                row['observer_export_seconds'] = value.get('export_seconds')
                row['observer_open_seconds'] = (value.get('open') or {}).get('open_seconds')
                row['observer_stage'] = value.get('stage')
            elif opened is not None and opened.is_file():
                row['observer_open_seconds'] = json.loads(opened.read_text()).get('open_seconds')
        except Exception as error:
            row['timing_error'] = str(error)[:200]
        (rows_dir / (key + '.json')).write_text(json.dumps(row, sort_keys=True) + '\n')
        print('P5ROW ' + json.dumps(compact(row), sort_keys=True, separators=(',', ':')), flush=True)
        return key, row

    def occt(key):
        if over_budget():
            return key, None
        began = time.monotonic()
        folder = args.out / 'raw'
        folder.mkdir(parents=True, exist_ok=True)
        output = folder / (key + '.json')
        row = {}
        try:
            assert oic.digest(files[key]) == hashes[key], 'source hash'
            with (folder / (key + '.stdout')).open('wb') as out, (folder / (key + '.stderr')).open('wb') as err:
                process = rss_guard.run([sys.executable, '-B', str(Path(__file__).with_name('occt_census.py')),
                                         'probe', str(files[key]), str(output)], stdout=out, stderr=err,
                                        timeout=args.timeout, max_rss_mb=args.max_rss_mb)
            row = json.loads(output.read_text()) if output.exists() else {}
            if process.returncode:
                row.update(outcome='failed', failure='process_exit_' + str(process.returncode))
        except subprocess.TimeoutExpired:
            row = json.loads(output.read_text()) if output.exists() else {}
            row.update(outcome='timeout', failure='external_deadline')
        except rss_guard.MemoryExceeded as error:
            row = json.loads(output.read_text()) if output.exists() else {}
            row.update(outcome='out_of_memory', failure='memory_cap', peak_mb=error.peak_mb)
        except Exception as error:
            row.update(outcome='driver_error', failure=type(error).__name__ + ': ' + str(error)[:300])
        row.update(model_id=key, wall_seconds=time.monotonic() - began)
        (rows_dir / (key + '.json')).write_text(json.dumps(row, sort_keys=True) + '\n')
        keep = ['model_id', 'outcome', 'stage', 'failure', 'wall_seconds', 'elapsed_seconds', 'opened',
                'solid_count', 'valid_solid_count', 'healed_whole_shape_valid',
                'transferred_whole_shape_valid', 'open_shell_count', 'ocp_version', 'peak_mb']
        print('P5OCCT ' + json.dumps({k: row[k] for k in keep if row.get(k) is not None}, sort_keys=True,
                                     separators=(',', ':')), flush=True)
        return key, row

    one = catallax if args.kernel == 'catallax' else occt
    count = skipped = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(one, k) for k in todo]):
            key, row = future.result()
            if row is None:
                skipped += 1
                continue
            count += 1
    print('P5SUMMARY ' + json.dumps(dict(kernel=args.kernel, complete=skipped == 0, ran=count, skipped=skipped,
                                         rows=len(list(rows_dir.glob('*.json'))),
                                         seconds=round(time.monotonic() - started))), flush=True)


if __name__ == '__main__':
    main()
