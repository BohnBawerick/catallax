#!/usr/bin/env python3
"""Project a completed raw run into reviewable data and a short outcome summary."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

from policy import TOLERANCES_MM
from run import save, summary


def project(path, data, markdown, supplements=()):
    raw = json.loads(path.read_text())
    if raw['status'] != 'completed' or len(raw['cases']) != len(raw['inventory']):
        raise ValueError('only a complete inventory can be projected')
    result = copy.deepcopy(raw)
    result['raw_result_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    result['projection_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    sources = [path, *supplements]
    records = [raw, *(json.loads(p.read_text()) for p in supplements)]
    result['observation_runs'] = []
    result['superseded_cases'] = []
    positions = {case['id']: index for index, case in enumerate(result['cases'])}
    for index, (source, record) in enumerate(zip(sources, records)):
        if record['status'] != 'completed' or len(record['cases']) != len(record['inventory']):
            raise ValueError('supplement has an incomplete inventory')
        if any(record[key] != raw[key] for key in ('binary_sha256', 'observer_sha256', 'policy')):
            raise ValueError('supplements must use the same kernel binaries and completion policy')
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        header = {key: value for key, value in record.items() if key not in ('cases', 'summary', 'inventory')}
        header['raw_result_sha256'] = digest
        result['observation_runs'].append(header)
        for case in record['cases']:
            name = case['id']
            if name not in positions or record['case_sha256'][name] != raw['case_sha256'][name]:
                raise ValueError('supplement changes a case definition')
            replacement = copy.deepcopy(case)
            replacement['observation_run'] = index
            if index:
                result['superseded_cases'].append(result['cases'][positions[name]])
            result['cases'][positions[name]] = replacement
    result['executed_producer_attempts'] = sum(len(record['cases'])*3 for record in records)
    result['completed_utc'] = records[-1]['completed_utc']
    for case in [*result['cases'], *result['superseded_cases']]:
        if set(case['engines']) != {'catallax', 'occt', 'smlib'}:
            raise ValueError('every case must retain all three attempted engines')
        for engine in case['engines'].values():
            checked = engine['check']
            checked.pop('points', None)
            reverse = checked.get('reverse', {})
            distances = reverse.pop('distances', None)
            if distances is not None:
                reverse['sample_count'] = len(distances)
                reverse['distances_sha256'] = hashlib.sha256(
                    json.dumps(distances, sort_keys=True).encode()).hexdigest()
            process = checked.get('process', {})
            if 'stderr' in process and not process['stderr']:
                process.pop('stderr')
    result['summary'] = summary(result)
    save(data, result)
    source_label = f"`{raw['source_tag']}` at `{raw['source_revision']}`" if raw.get('source_tag') else f"`{raw['source_revision']}`"
    lines = ['# Kernel comparison run', '',
             f"Kernel source {source_label}. Started {raw['started_utc']}. Shared development host, not quiet-host timing evidence.", '',
             f"{len(raw['cases'])} case variants, each attempted sequentially with Catallax, OpenCASCADE and SMLib. {result['executed_producer_attempts']} producer attempts are retained across {len(records)} runs, including superseded attempts. No site pages were regenerated.", '',
             f"The acceptance policy is `tools/kernel-compare/policy.py`. Every producer has the same {raw['policy']['solve_limit_seconds']:g}-second whole-process solve limit and the same geometric thresholds. A successful build is not enough. `done-sampled` requires valid raw STEP, passing authored witnesses, and both finite boundary-sample distance checks within the threshold. Unknown error bounds remain unresolved. This is not a global Hausdorff or correctness certificate.", '',
             'The reference is an independently compiled Catallax field definition. It shares recipe semantics and authorship with Catallax B-rep. OpenCASCADE both competes and imports all STEP files. Unsupported reference definitions, absent bounds and observer timeouts cannot become passes.', '',
             '| Tolerance mm | Kernel | Sampled done | Not done | Refused | Solve time limit | Other or unresolved |',
             '| ---: | --- | ---: | ---: | ---: | ---: | ---: |']
    for tolerance in TOLERANCES_MM:
        counts = result['summary'][str(tolerance)]['outcomes']
        for engine, rows in counts.items():
            named = ('done-sampled', 'not-done', 'refused', 'unsolved-time-limit')
            other = sum(value for key, value in rows.items() if key not in named)
            lines.append(f'| {tolerance:g} | {engine} | '+ ' | '.join(str(rows.get(k, 0)) for k in named)+f' | {other} |')
    lines.extend(['', '## Differences at 0.001 mm', ''])
    middle = result['summary'][str(TOLERANCES_MM[2])]
    for engine, cases in middle['completed_where_another_did_not'].items():
        examples = ', '.join(f'`{c}`' for c in cases[:4]) or 'none'
        lines.append(f"- {engine}: {len(cases)} cases sampled done where another engine was not. Examples: {examples}. Complete lists are in the JSON.")
    control = next((case for case in result['cases'] if case['id'] == 'fillets/oblique-cylinder'), None)
    if control is not None:
        lines.extend(['', '## Analytic attribution control', '',
                      'The oblique-cylinder corpus carries an independently derived Steiner-formula volume. This scalar diagnostic does not replace the missing global distance certificate.', '',
                      '| Kernel | Volume error mm³ | Outcome |', '| --- | ---: | --- |'])
        for engine, row in control['engines'].items():
            scalar = row['check'].get('analytic_volume_check')
            error = f"{scalar['error_mm3']:.12g}" if scalar is not None else 'unmeasured'
            lines.append(f"| {engine} | {error} | {row['solve']['status']} |")
    lines.extend(['', 'Refusals, crashes, protocol errors, adapter gaps and timeouts remain distinct in the JSON. An unresolved observer is not a kernel loss. The sweep reclassifies one unchanged result at five thresholds; it does not silently retune or retry a kernel.', '',
                  'The committed JSON keeps per-case solve times, error intervals, threshold outcomes, source and binary identities, native library hashes, package versions, machine identity and dates. Raw commands, STEP files, standard streams and per-sample intervals remain in the external run directory, bound by the raw-result digest.'])
    markdown.parent.mkdir(parents=True, exist_ok=True)
    markdown.write_text('\n'.join(lines) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('raw', type=Path)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--supplement', type=Path, action='append', default=[])
    args = parser.parse_args()
    project(args.raw, args.data, args.summary, args.supplement)
