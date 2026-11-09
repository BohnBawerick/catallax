#!/usr/bin/env python3
"""Project a completed raw run into reviewable data and a short outcome summary."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

from policy import TOLERANCES_MM
from run import save, summary


def project(path, data, markdown):
    raw = json.loads(path.read_text())
    if raw['status'] != 'completed' or len(raw['cases']) != len(raw['inventory']):
        raise ValueError('only a complete inventory can be projected')
    result = copy.deepcopy(raw)
    result['raw_result_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    for case in result['cases']:
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
    lines = ['# Kernel comparison run', '',
             f"Kernel source `{raw['source_revision']}`. Started {raw['started_utc']}. Shared development host, not quiet-host timing evidence.", '',
             f"{len(raw['cases'])} case variants, each attempted sequentially with Catallax, OpenCASCADE and SMLib. No site pages were regenerated.", '',
             'The acceptance policy is `tools/kernel-compare/policy.py`. Every producer has the same 120-second whole-process solve limit and the same geometric thresholds. A successful build is not enough. `done-sampled` requires valid raw STEP, passing authored witnesses, and both finite boundary-sample distance checks within the threshold. Unknown error bounds remain unresolved. This is not a global Hausdorff or correctness certificate.', '',
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
        lines.append(f"- {engine}: sampled done where another engine was not: " + (', '.join(f'`{c}`' for c in cases) or 'none') + '.')
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
    args = parser.parse_args()
    project(args.raw, args.data, args.summary)
