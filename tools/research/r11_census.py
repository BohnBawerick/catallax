#!/usr/bin/env python3
"""Collect public R11 logs using the existing verifier driver; never read private replay."""
import argparse
import collections
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()


def inventory(root):
    # Deliberately exclude tools/checks inventories, data summaries and all
    # external checkouts. Public saved replay means exported question/1 here.
    names = git(root, 'ls-files', 'crates', '.claude', 'parts').splitlines()
    parts, replay, excluded = [], [], []
    for name in names:
        if not name.endswith(('.json', '.jsonl')):
            continue
        if name.startswith('parts/') and re.fullmatch(r'parts/(P0/[^/]+|P[1-5][^/]*)\.json', name):
            parts.append(name)
            continue
        if name.endswith('.jsonl'):
            excluded.append({'path': name, 'reason': 'not an exported question/1 JSON document'})
            continue
        value = json.loads((root / name).read_text())
        fmt = value.get('format') if isinstance(value, dict) else None
        if fmt != 'catallax-question/1':
            excluded.append({'path': name, 'reason': 'outside complete P0-P5 and exported question/1 scope', 'format': fmt})
            continue
        # No file-dependent or oracle adapter is admitted by this collection.
        if value.get('source') != 'recipe' or any(
            s.get('op', '').startswith('import') for s in value.get('problem', {}).get('steps', [])
        ):
            excluded.append({'path': name, 'reason': 'external-file or non-recipe replay requires separate eligibility review'})
            continue
        if value.get('max_work') not in (1_000_000, 20_000_000, 100_000_000):
            excluded.append({'path': name, 'reason': 'saved allowance is outside the three authorized profiles'})
            continue
        replay.append(name)
    return parts, replay, excluded


def count(path):
    evidence = json.loads(path.read_text())
    selected, refusals, roots, seen = [], [], collections.Counter(), {}
    for record in evidence['records']:
        roots[(record['group'], record['status'])] += 1
        for event in record['operation_log']['operations']:
            # Artifact key identifies the driver execution; root/manifest copies
            # never enter this loop. Different adapters are different executions.
            key = (record['group'], record['artifact_key'], event['sequence'])
            if key in seen:
                assert seen[key] == event, 'conflicting duplicate event'
                continue
            seen[key] = event
            item = {'input': record['input'], 'group': record['group'],
                    'max_work': record['max_work'], **event}
            if event['operation'] in ('shell', 'offset', 'thicken'):
                selected.append(item)
            if event['status'] != 'answered':
                refusals.append(item)
    totals = collections.Counter((e['operation'], e['status']) for e in selected)
    result = {'requests': len(evidence['records']), 'unique_inputs': len({r['input'] for r in evidence['records']}),
              'unique_executed_spans': len(seen), 'selected_events': selected,
              'selected_totals': [{'operation': op, 'status': status, 'count': n} for (op, status), n in sorted(totals.items())],
              'request_statuses': [{'group': group, 'status': status, 'count': n} for (group, status), n in sorted(roots.items())],
              'nonanswered_spans': refusals,
              'incomplete_groups': [g for g in evidence['groups'] if g['state'] != 'completed']}
    print(json.dumps(result, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count', type=Path, help='count a preserved evidence.json without executing geometry')
    parser.add_argument('--source-root', type=Path)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--workers', type=int, choices=range(1, 5), default=4)
    args = parser.parse_args()
    if args.count:
        count(args.count)
        return
    if not all((args.source_root, args.binary, args.out)):
        parser.error('collection requires --source-root, --binary and --out')
    root, binary, out = args.source_root.resolve(), args.binary.resolve(), args.out.resolve()
    if git(root, 'status', '--porcelain', '--untracked-files=no'):
        raise SystemExit('source checkout has tracked modifications')
    out.mkdir(parents=True, exist_ok=False)
    revision = git(root, 'rev-parse', 'HEAD')
    parts, replay, excluded = inventory(root)
    driver_path = root / '.claude/skills/verify-catallax/scripts/control-catallax.py'
    spec = importlib.util.spec_from_file_location('r11_driver', driver_path)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    doctor = subprocess.run(['python3', '-B', str(driver_path), 'doctor', '--binary', str(binary), '--json'], capture_output=True, text=True)
    (out / 'doctor.json').write_text(doctor.stdout)
    (out / 'doctor.stderr').write_text(doctor.stderr)
    if doctor.returncode:
        raise SystemExit('doctor failed; retained diagnostics')
    about = json.loads(doctor.stdout)
    assert json.loads(about['about']['kernel'])['commit'] == revision, 'binary/source commit mismatch'
    evidence = {'format': 'catallax-r11-census/1', 'source_revision': revision,
                'binary_sha256': sha(binary), 'driver_sha256': sha(driver_path),
                'collector_sha256': sha(Path(__file__)), 'started_utc': utc(),
                'workers': args.workers, 'eligible_parts': parts, 'eligible_replays': replay,
                'excluded_public_files': excluded,
                'excluded_sources': [{'source': 'protected 477-case replay', 'reason': 'private inputs per docs/RELEASING.md; not opened or enumerated'},
                                     {'source': 'other public files outside crates/.claude/parts', 'reason': 'not included in exported-question discovery; no complete-repository replay claim'},
                                     {'source': 'C++/external oracle checks', 'reason': 'not authorized; no oracle run'}],
                'groups': [], 'records': [], 'observation_timeout_seconds': None}
    # Doctor may contain machine paths; preserve raw privately and publish only
    # its protocol/identity subsection, plus its raw digest.
    evidence['doctor_sha256'] = sha(out / 'doctor.json')
    evidence['about'] = about.get('about')
    groups = [(f'parts-{source}', parts, source, 1_000_000) for source in ('recipe', 'field', 'compare')]
    groups.append(('public-replay', replay, 'auto', 100_000_000))
    for group, names, source, cap in groups:
        info = {'name': group, 'started_utc': utc(), 'inputs': len(names), 'default_max_work': cap}
        evidence['groups'].append(info)
        write(out / 'evidence.json', evidence)
        try:
            driver.run(argparse.Namespace(binary=str(binary), inputs=[str(root / n) for n in names],
                       source=source, adapter=None, out=str(out / group), cache=None, cores=None,
                       workers=args.workers, max_work=cap, fresh=True, dry_run=False))
            run = json.loads((out / group / 'run.json').read_text())
            items = json.loads((out / group / 'inputs.json').read_text())
            assert len(run['cases']) == len(names) == len(items)
            info['execution'] = run['execution']
            for name, item, case in zip(names, items, run['cases']):
                folder = Path(case['path'])
                response = case['response']
                op_log = json.loads((folder / 'operations.json').read_text())
                assert op_log == response['operation_log']
                assert not case['cache_hit']
                assert response['manifest']['kernel']['commit'] == revision
                assert response['work']['limit'] == item.get('max_work', cap)
                evidence['records'].append({'input': name, 'input_sha256': sha(root / name),
                    'group': group, 'source': item['source'], 'max_work': item.get('max_work', cap),
                    'artifact_key': case['key'], 'operation_log_sha256': sha(folder / 'operations.json'),
                    'response_sha256': sha(folder / 'response.json'),
                    'status': response['manifest']['status'], 'refusal': response['manifest'].get('refusal'),
                    'comparison': response['manifest'].get('comparison'),
                    'kernel': response['manifest'].get('kernel'), 'work': response['work'], 'operation_log': op_log})
            info['state'] = 'completed'
        except Exception as error:
            # A failed group is not a set of geometric refusals. Checkpoints and
            # pending artifacts remain in the raw evidence directory for review.
            info['state'] = 'incomplete'
            info['error_type'] = type(error).__name__
            (out / f'{group}-error.txt').write_text(str(error))
        info['finished_utc'] = utc()
        write(out / 'evidence.json', evidence)
        print(group, info['state'], flush=True)
    assert sha(binary) == evidence['binary_sha256']
    evidence['finished_utc'] = utc()
    write(out / 'evidence.json', evidence)
    if any(g['state'] != 'completed' for g in evidence['groups']):
        raise SystemExit('incomplete collection; do not report full coverage')


if __name__ == '__main__':
    main()
