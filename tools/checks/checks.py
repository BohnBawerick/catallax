#!/usr/bin/env python3
"""Conservative merge selection and exact-run evidence. Never prints replay inputs."""
import argparse
from collections import Counter
import fnmatch
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import tarfile
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
RULES = 'tools/checks/rules.json'
INVENTORY = 'tools/checks/earlier-answers.json'
DRIVER = '.claude/skills/verify-catallax/scripts/control-catallax.py'
RAW_POLICY = 'tools/checks/known-raw-invalid.json'
CORRECTIONS = 'tools/checks/replay-corrections.json'
STEP_CORRECTION_FIELDS = {
    'step.added_edges', 'step.added_vertices',
    'step.max_deviation', 'step.approximated_curves'}
# Captain's 29 October 2026 option-A ruling. Data cannot grant itself authority.
RAW_STEP_RULING_DIGESTS = {
    '7755b0bd5ef52fd26b40197d73a8e9d767683b78644b24d496795dc0b4217e9e'}


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args])


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def matches(path, patterns):
    return any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns)


def cargo_graph(directory):
    # No target filter: metadata retains optional, build, dev and target-specific
    # dependencies. All features are requested; only workspace edges select jobs.
    metadata = json.loads(subprocess.check_output([
        'cargo', 'metadata', '--format-version', '1', '--no-deps', '--all-features'],
        cwd=directory))
    names = {package['name'] for package in metadata['packages']}
    all_edges, build_edges = {}, {}
    for package in metadata['packages']:
        deps = package['dependencies']
        all_edges[package['name']] = {dep['name'] for dep in deps if dep['name'] in names}
        build_edges[package['name']] = {dep['name'] for dep in deps
                                        if dep['name'] in names and dep['kind'] != 'dev'}
    return all_edges, build_edges


def closure(seeds, edges):
    result = set(seeds)
    pending = list(seeds)
    while pending:
        for target in edges.get(pending.pop(), ()):
            if target not in result:
                result.add(target)
                pending.append(target)
    return result


def select(paths, rules, graph, build_graph, event='pull_request', ref=''):
    names = set(graph)
    reverse = {name: {consumer for consumer, deps in graph.items() if name in deps}
               for name in names}
    binary = closure({'catallax-cli'}, build_graph)
    packages, reasons = set(), []
    replay = scoreboard = build = False
    full_event = event == 'workflow_dispatch' or ref.startswith('refs/heads/release/')
    for path in paths:
        if full_event or matches(path, rules['full_paths']):
            packages |= names
            replay = scoreboard = build = True
            reasons.append(f'{path}: check policy or shared build input')
        elif path.startswith('crates/') and len(path.split('/')) > 2 and path.split('/')[1] in names:
            owner = path.split('/')[1]
            packages |= closure({owner}, reverse)
            build = build or owner in binary
            scoreboard = scoreboard or owner in binary
            replay = replay or (owner in binary and matches(path, rules['replay_paths']))
            reasons.append(f'{path}: crate and reverse consumers, unknown corpus owners')
        elif path.startswith('parts/'):
            packages |= closure({'catallax-app', 'catallax-cli', 'catallax-compare',
                                 'catallax-features', 'catallax-field', 'catallax-iges'}, reverse)
            scoreboard = build = True
            reasons.append(f'{path}: target-part fixtures and consumers')
        elif path.startswith('.claude/skills/verify-catallax/'):
            packages |= closure({'catallax-app', 'catallax-cli', 'catallax-measure'}, reverse)
            build = True
            replay = replay or matches(path, rules['replay_paths'])
            scoreboard = True
            reasons.append(f'{path}: verification driver or shared example')
        elif matches(path, rules['documentation_paths']):
            reasons.append(f'{path}: documentation only')
        else:
            packages |= names
            replay = scoreboard = build = True
            reasons.append(f'{path}: unknown mapping, full single-configuration checks')
    if full_event:
        packages = names
        replay = scoreboard = build = True
        reasons.append('release push or manual dispatch: every check')
    jobs = {'fmt': 1}
    if packages:
        jobs.update(clippy=1, tests=rules['jobs']['tests'])
    if build:
        jobs['build'] = 1
    if replay:
        jobs['replay'] = 1
    if scoreboard:
        jobs['scoreboard'] = 1
    # A newly declared job is mandatory until its own reviewed routing exists.
    for name, partitions in rules['jobs'].items():
        if name not in {'fmt', 'clippy', 'tests', 'build', 'replay', 'scoreboard'}:
            jobs[name] = partitions
    return {'packages': sorted(packages), 'jobs': jobs, 'selected_cases': 477 if replay else 0,
            'reasons': reasons}


def union_selection(base, candidate):
    return {'packages': sorted(set(base['packages']) | set(candidate['packages'])),
            'jobs': {name: max(base['jobs'].get(name, 0), candidate['jobs'].get(name, 0))
                     for name in base['jobs'].keys() | candidate['jobs'].keys()},
            'selected_cases': max(base['selected_cases'], candidate['selected_cases']),
            'base_selected': base['selected_cases'], 'candidate_selected': candidate['selected_cases'],
            'reasons': sorted(set(base['reasons'] + candidate['reasons']))}


def metadata_errors(base, candidate, corrections=None, base_policy=None):
    if base == candidate:
        return []
    failure = ['accepted-answer metadata changed without base approval']
    if not corrections or corrections.get('format') != 'catallax-replay-corrections/1':
        return failure
    if not base_policy or raw_step_errors([], base_policy):
        return failure
    rows = corrections.get('corrections', [])
    if len({row['index'] for row in rows}) != len(rows) or len(base) != len(candidate):
        return failure
    approved = {row['index']: row['step_sha256']
                for ruling in base_policy['rulings'] for row in ruling['approved']}
    by_index = {row['index']: row for row in rows}
    for before, after in zip(base, candidate):
        if before == after:
            continue
        index = before['index']
        correction = by_index.get(index)
        if not correction or index not in base_policy['known_raw_invalid']:
            return failure
        if (correction.get('class') != 'invalid-to-valid-step'
                or correction.get('request_sha256') != before['request_sha256']
                or correction.get('before_sha256') != before['accepted_answer_sha256']
                or correction.get('before_manifest_sha256') != before['accepted_manifest_sha256']
                or correction.get('before_step_sha256') != before['accepted_files']['answer.step']
                or correction.get('before_step_sha256') != approved.get(index)
                or not set(correction.get('manifest_fields', [])) <= STEP_CORRECTION_FIELDS
                or not valid_digest(correction.get('unchanged_manifest_sha256'))):
            return failure
        expected = json.loads(json.dumps(before))
        expected['accepted_answer_sha256'] = correction.get('after_sha256')
        expected['accepted_manifest_sha256'] = correction.get('after_manifest_sha256')
        expected['accepted_files']['answer.step'] = correction.get('after_step_sha256')
        if (expected != after or any(not valid_digest(value) for value in (
                expected['accepted_answer_sha256'], expected['accepted_manifest_sha256'],
                expected['accepted_files']['answer.step']))):
            return failure
    return []


def valid_digest(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def protected_manifest_digest(manifest, driver):
    comparable = driver.comparable(manifest)
    comparable['step'] = dict(comparable['step'])
    for field in STEP_CORRECTION_FIELDS:
        comparable['step'].pop(field.split('.')[1], None)
    return driver.digest(comparable)


def correction_manifest_errors(manifest, correction, driver, base_manifest):
    if driver.digest(driver.comparable(base_manifest)) != correction['before_manifest_sha256']:
        return ['base manifest does not match the accepted before manifest']
    protected = protected_manifest_digest(base_manifest, driver)
    if protected != correction['unchanged_manifest_sha256']:
        return ['correction commitment does not match the authenticated base manifest']
    if protected_manifest_digest(manifest, driver) != protected:
        return ['correction changes a manifest field outside the four STEP fields']
    return []


def raw_step_errors(observations, policy):
    if policy.get('format') != 'catallax-known-raw-invalid/1':
        return ['unsupported known-raw-invalid policy']
    rulings = policy.get('rulings')
    ruling_digest = hashlib.sha256(json.dumps(
        rulings, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
    if ruling_digest not in RAW_STEP_RULING_DIGESTS:
        return ['known-raw-invalid expansion requires an explicit trusted captain ruling']
    approved = {row['index']: row['step_sha256']
                for ruling in rulings for row in ruling['approved']}
    known = policy.get('known_raw_invalid')
    if (not isinstance(known, list) or any(type(index) is not int for index in known)
            or len(known) != len(set(known))):
        return ['invalid or duplicate known-raw-invalid indices']
    errors = [f'index {index}: no captain ruling for new raw-invalid entry'
              for index in set(known) - approved.keys()]
    for row in observations:
        index, fingerprint = row['index'], row['step_sha256']
        if fingerprint is None:
            if index in known:
                errors.append(f'index {index}: listed raw-invalid STEP is missing')
        elif row['valid']:
            if index in known:
                errors.append(f'index {index}: raw-valid upgrade must leave the known-raw-invalid list')
        elif index not in known:
            errors.append(f'index {index}: unlisted raw-invalid STEP')
        elif fingerprint != approved.get(index):
            errors.append(f'index {index}: STEP bytes changed while still raw-invalid')
    return errors


def decide(plan, reported, records, binary_sha256=None, actual_jobs=()):
    errors = []
    required, declared = plan['selection']['jobs'], plan['declared']
    allowed_names = set(declared) | {'select', 'check'}
    for name, count in required.items():
        if count > 1:
            allowed_names.discard(name)
            allowed_names.update(f'{name} ({part})' for part in range(1, count + 1))
    for actual in actual_jobs:
        if actual.split(' / ', 1)[0] not in allowed_names:
            errors.append(f'{actual}: unexpected job')
    for name in reported.keys() - set(declared):
        errors.append(f'{name}: unexpected job')
    for name in set(declared) | set(required):
        result = reported.get(name, {}).get('result')
        if name not in required:
            if result not in (None, 'skipped'):
                errors.append(f'{name}: unexpected selected job ({result})')
            continue
        if result != 'success':
            errors.append(f'{name}: {result or "missing"}')
    seen = set()
    for record in records:
        name, part = record.get('job'), record.get('partition')
        key = (name, part)
        if name not in required:
            errors.append(f'{name}: unexpected evidence')
            continue
        if type(part) is not int or not 1 <= part <= required[name]:
            errors.append(f'{name}: unexpected partition {part}')
            continue
        if key in seen:
            errors.append(f'{name}: duplicate partition {part}')
        seen.add(key)
        expected = dict(plan['identity'])
        expected['result'] = 'success'
        if name in {'build', 'replay', 'scoreboard'}:
            expected['binary_sha256'] = binary_sha256
        else:
            expected['binary_sha256'] = None
        if name == 'replay':
            expected.update(workers=4, fresh=True, cache_hits=0,
                            selected_cases=plan['selection']['selected_cases'],
                            verified_corrections=[row['index'] for row in plan['correction_base_answers']])
        for field, value in expected.items():
            if type(record.get(field)) is not type(value) or record.get(field) != value:
                errors.append(f'{name}: stale {field}: evidence={record.get(field)!r}, candidate={value!r}')
        if name in {'build', 'replay', 'scoreboard'} and not binary_sha256:
            errors.append(f'{name}: missing binary_sha256 from current build')
    for name, count in required.items():
        if name in {'release-evidence', 'release-dry-run'}:
            continue  # Reusable release workflows retain their own artifact gates.
        for part in range(1, count + 1):
            if (name, part) not in seen:
                errors.append(f'{name}: partition {part}/{count} missing')
    return errors


def prepare(args):
    base_rules = json.loads(git('show', f'{args.base}:{RULES}')) if git_has(args.base, RULES) else read(ROOT / RULES)
    base_metadata = json.loads(git('show', f'{args.base}:{INVENTORY}')) if git_has(args.base, INVENTORY) else read(ROOT / INVENTORY)
    candidate_rules, candidate_metadata = read(ROOT / RULES), read(ROOT / INVENTORY)
    base_policy = (json.loads(git('show', f'{args.base}:{RAW_POLICY}'))
                   if git_has(args.base, RAW_POLICY) else read(ROOT / RAW_POLICY))
    errors = metadata_errors(base_metadata, candidate_metadata, read(ROOT / CORRECTIONS), base_policy)
    errors += raw_step_errors([], read(ROOT / RAW_POLICY))
    errors += raw_step_errors([], base_policy)
    if errors:
        raise ValueError('; '.join(errors))
    if base_rules['workers'] != 4 or candidate_rules['workers'] != 4:
        raise ValueError('fixed workers must be 4 under both rules')
    paths = [p.decode() for p in git('diff', '--name-only', '--no-renames', '-z', args.base, 'HEAD').split(b'\0') if p]
    with tempfile.TemporaryDirectory(prefix='catallax-base-metadata-') as temporary:
        with tarfile.open(fileobj=io.BytesIO(git('archive', args.base))) as archive:
            archive.extractall(temporary, filter='data')
        base_graph, base_build = cargo_graph(temporary)
    graph, build = cargo_graph(ROOT)
    before = select(paths, base_rules, base_graph, base_build, args.event, args.ref)
    after = select(paths, candidate_rules, graph, build, args.event, args.ref)
    selection = union_selection(before, after)
    # A removed crate contributes its old consumers, not an impossible -p flag.
    selection['packages'] = sorted(set(selection['packages']) & set(graph))
    if args.release_evidence:
        selection['jobs']['release-evidence'] = 1
    if args.release_dry_run:
        selection['jobs']['release-dry-run'] = 1
    plan = {'format': 'catallax-check-plan/1', 'base': args.base, 'head': args.head,
            'correction_base_answers': [before for before, after in zip(base_metadata, candidate_metadata)
                                        if before != after],
            'declared': sorted(set(base_rules['jobs']) | set(candidate_rules['jobs']) |
                               {'select', 'release-evidence', 'release-dry-run', 'timings'}),
            'identity': {'commit': git('rev-parse', 'HEAD').decode().strip(),
                         'head': args.head,
                         'run_id': str(args.run_id), 'run_attempt': str(args.run_attempt),
                         'corpus_sha256': base_rules['corpus_sha256'],
                         'raw_policy_sha256': sha(ROOT / RAW_POLICY),
                         'corrections_sha256': sha(ROOT / CORRECTIONS),
                         'metadata_sha256': sha(ROOT / INVENTORY), 'driver_sha256': sha(ROOT / DRIVER)},
            'selection': selection}
    plan['decisions'] = {
        name: {'selected': name in selection['jobs'],
               'partitions': selection['jobs'].get(name, 0),
               'reasons': ['always: formatting and check policy'] if name == 'fmt'
               else selection['reasons'] if name in selection['jobs']
               else ['no affected owner or consumer under either rule set']}
        for name in plan['declared'] if name != 'select'}
    plan['decisions']['timings'] = {
        'selected': True, 'partitions': 1, 'required_by_gate': False,
        'reasons': ['observational timing artifact after the final gate']}
    write(args.out, plan)
    if args.github_output:
        with Path(args.github_output).open('a') as out:
            for name in candidate_rules['jobs']:
                out.write(f'{name}={str(name in selection["jobs"]).lower()}\n')
            out.write('packages=' + ' '.join('-p ' + name for name in selection['packages']) + '\n')
            out.write('base_corrections=' + str(bool(plan['correction_base_answers'])).lower() + '\n')
    print(json.dumps(selection, sort_keys=True))


def git_has(commit, path):
    return subprocess.run(['git', '-C', str(ROOT), 'cat-file', '-e', f'{commit}:{path}'],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def evidence(args):
    plan = read(args.plan)
    record = {**plan['identity'], 'job': args.job, 'partition': args.partition,
              'binary_sha256': sha(args.binary) if args.binary else None, 'result': args.result}
    if args.job == 'replay':
        summary = read(args.summary)
        record.update(summary)
    write(args.out, record)


def gate(args):
    plan = read(args.plan)
    current = {'commit': git('rev-parse', 'HEAD').decode().strip(),
               'run_id': os.environ.get('GITHUB_RUN_ID', 'local'),
               'run_attempt': os.environ.get('GITHUB_RUN_ATTEMPT', '1'),
               'raw_policy_sha256': sha(ROOT / RAW_POLICY),
               'corrections_sha256': sha(ROOT / CORRECTIONS),
               'metadata_sha256': sha(ROOT / INVENTORY), 'driver_sha256': sha(ROOT / DRIVER)}
    for field, value in current.items():
        if plan['identity'].get(field) != value:
            raise ValueError(f'select: stale {field}')
    # The selector is mandatory even if a candidate deletes its own dependency.
    reported = json.loads(os.environ['RESULTS'])
    if reported.get('select', {}).get('result') != 'success':
        raise ValueError('select: missing or failed')
    reported.pop('select', None)
    records = [read(path) for path in Path(args.evidence).glob('*.json')]
    builds = [record for record in records if record.get('job') == 'build']
    binary_hash = builds[0].get('binary_sha256') if len(builds) == 1 else None
    errors = decide(plan, reported, records, binary_hash, Path(args.jobs).read_text().splitlines())
    if errors:
        raise ValueError('\n'.join(errors))
    print('Every selected check and partition passed with exact-candidate evidence')


def build_base(args):
    plan = read(args.plan)
    if not plan['correction_base_answers']:
        raise ValueError('no new STEP corrections require a base verifier')
    target = ROOT / 'target' / 'check-base'
    with tempfile.TemporaryDirectory(prefix='catallax-base-verifier-') as temporary:
        with tarfile.open(fileobj=io.BytesIO(git('archive', plan['base']))) as archive:
            archive.extractall(temporary, filter='data')
        subprocess.run([
            'cargo', 'build', '--release', '--locked', '-p', 'catallax-cli',
            '--bin', 'catallax-verify', '--manifest-path', str(Path(temporary) / 'Cargo.toml'),
            '--target-dir', str(target)], cwd=temporary, check=True)
    destination = Path(args.out)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(target / 'release' / 'catallax-verify', destination)
    print(json.dumps({'base_commit': plan['base'], 'binary_sha256': sha(destination)}))


def run_fresh(driver, binary, requests, out, workers):
    run_args = argparse.Namespace(binary=binary, inputs=[str(requests)], source='auto',
                                 adapter=None, out=str(out), cache=None, cores=None,
                                 max_work=100_000_000, workers=workers, fresh=True, dry_run=False)
    try:
        driver.run(run_args)
    except Exception:
        raise ValueError('replay execution failed; private diagnostics retained on runner') from None
    return read(Path(out) / 'run.json')


def correction_bases(args, plan, driver, release, inputs):
    accepted = plan['correction_base_answers']
    if not accepted:
        return {}
    if not args.base_binary:
        raise ValueError('new STEP corrections require the base verifier')
    source = Path(args.out) / 'base-requests.jsonl'
    source.parent.mkdir(parents=True, exist_ok=True)
    with source.open('w') as stream:
        for row in accepted:
            stream.write(driver.canonical({'max_work': 100_000_000, **inputs[row['index']],
                                           'format': driver.CONTRACT}) + '\n')
    report = run_fresh(driver, args.base_binary, source, Path(args.out) / 'base', 1)
    if (len(report['cases']) != len(accepted)
            or report['execution']['worker_exits'] != {'0': 0}
            or report['config']['binary_sha256'] != sha(args.base_binary)):
        raise ValueError('incomplete base correction replay')
    manifests = {}
    for row, case in zip(accepted, report['cases']):
        directory = Path(case['path'])
        request = read(directory / 'request.json')
        request.pop('out', None)
        manifest = case['response']['manifest']
        if (driver.digest(request) != row['request_sha256']
                or driver.digest(driver.comparable(manifest)) != row['accepted_manifest_sha256']
                or release.answer_digest(driver, directory) != row['accepted_answer_sha256']
                or sha(directory / 'answer.step') != row['accepted_files']['answer.step']
                or case['cache_hit'] or case['response']['work']['used'] > 100_000_000):
            raise ValueError(f"index {row['index']}: base replay does not match the accepted answer")
        manifests[row['index']] = manifest
    return manifests


def replay(args):
    started = time.monotonic()
    rules, metadata = read(ROOT / RULES), read(ROOT / INVENTORY)
    corrections = {row['index']: row for row in read(ROOT / CORRECTIONS)['corrections']}
    plan = read(args.plan)
    if (plan['identity']['metadata_sha256'] != sha(ROOT / INVENTORY)
            or plan['identity']['corrections_sha256'] != sha(ROOT / CORRECTIONS)
            or plan['identity']['commit'] != git('rev-parse', 'HEAD').decode().strip()):
        raise ValueError('stale correction replay plan')
    if sha(args.inputs) != rules['corpus_sha256']:
        raise ValueError('corpus_sha256 mismatch')
    if len(os.sched_getaffinity(0)) < 4:
        raise ValueError('fixed replay requires at least four available cores')
    driver = load_module('catallax_merge_driver', ROOT / DRIVER)
    release = load_module('catallax_merge_release', ROOT / 'tools/release/release.py')
    inputs = read(args.inputs)
    base_manifests = correction_bases(args, plan, driver, release, inputs)
    # Do not invoke the CLI: its ordinary JSON summary contains private ids.
    # The pinned file is an array; the public driver accepts request JSONL.
    normalized = Path(args.out) / 'replay-requests.jsonl'
    normalized.parent.mkdir(parents=True, exist_ok=True)
    with normalized.open('w') as stream:
        for row in inputs:
            stream.write(driver.canonical({'max_work': 100_000_000, **row,
                                           'format': driver.CONTRACT}) + '\n')
    report = run_fresh(driver, args.binary, normalized, args.out, 4)
    execution = report['execution']
    if execution['workers'] != 4 or execution['fresh'] is not True:
        raise ValueError('replay requires four fresh workers')
    if execution['worker_exits'] != {str(i): 0 for i in range(4)}:
        raise ValueError('replay has missing or failed workers')
    if report['config']['binary_sha256'] != sha(args.binary):
        raise ValueError('binary_sha256 changed during replay')
    failures, outcomes, work = [], Counter(), 0
    verified_corrections = []
    cases = report['cases']
    if len(cases) != len(metadata) or len(cases) != 477:
        raise ValueError('replay needs 477 complete cases')
    for index, (case, accepted) in enumerate(zip(cases, metadata)):
        directory = Path(case['path'])
        request = read(directory / 'request.json')
        request.pop('out', None)
        fingerprint = driver.digest(request)
        response = case['response']
        manifest = response['manifest']
        cause = None
        if fingerprint != accepted['request_sha256']:
            cause = 'request identity'
        answer_sha256 = release.answer_digest(driver, directory)
        if cause is None and answer_sha256 != accepted['accepted_answer_sha256']:
            manifest_sha256 = driver.digest(driver.comparable(manifest))
            if manifest_sha256 != accepted['accepted_manifest_sha256']:
                cause = 'manifest'
            else:
                changed = [name for name in release.OUTPUTS
                           if (sha(directory / name) if (directory / name).exists() else None)
                           != accepted['accepted_files'][name]]
                cause = 'STEP bytes' if changed == ['answer.step'] else 'output bytes'
        if index in base_manifests and correction_manifest_errors(
                manifest, corrections[index], driver, base_manifests[index]):
            cause = 'unapproved correction manifest field'
        if response['work']['used'] > 100_000_000:
            cause = 'above 100M work'
        if case['cache_hit'] or manifest.get('status') == 'internal-error':
            cause = 'cache or internal error'
        if cause:
            failures.append({'index': index, 'request_sha256': fingerprint, 'cause': cause,
                             'answer_sha256': answer_sha256})
        elif index in base_manifests:
            verified_corrections.append(index)
        refusal = manifest.get('refusal')
        outcomes[refusal['code'] if refusal else manifest['status']] += 1
        work += response['work']['used']
    summary = {'workers': execution['workers'], 'fresh': execution['fresh'],
               'cache_hits': sum(c['cache_hit'] for c in cases),
               'selected_cases': len(cases), 'elapsed_ms': (time.monotonic() - started) * 1000,
               'outcomes': dict(outcomes), 'work_total': work,
               'verified_corrections': verified_corrections,
               'equal_cases': len(cases) - len(failures), 'failures': failures}
    write(args.summary, summary)
    print(json.dumps(summary, sort_keys=True))
    if failures:
        raise ValueError('replay: changed earlier answers')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('select')
    p.add_argument('--base', required=True)
    p.add_argument('--head', default='')
    p.add_argument('--event', default='pull_request')
    p.add_argument('--ref', default='')
    p.add_argument('--run-id', default=os.environ.get('GITHUB_RUN_ID', 'local'))
    p.add_argument('--run-attempt', default=os.environ.get('GITHUB_RUN_ATTEMPT', '1'))
    p.add_argument('--release-evidence', action='store_true')
    p.add_argument('--release-dry-run', action='store_true')
    p.add_argument('--github-output')
    p.add_argument('--out', required=True)
    p.set_defaults(function=prepare)
    p = commands.add_parser('build-base')
    p.add_argument('--plan', required=True)
    p.add_argument('--out', required=True)
    p.set_defaults(function=build_base)
    p = commands.add_parser('evidence')
    p.add_argument('--plan', required=True)
    p.add_argument('--job', required=True)
    p.add_argument('--partition', type=int, default=1)
    p.add_argument('--result', default='success')
    p.add_argument('--binary')
    p.add_argument('--summary')
    p.add_argument('--out', required=True)
    p.set_defaults(function=evidence)
    p = commands.add_parser('gate')
    p.add_argument('--plan', required=True)
    p.add_argument('--evidence', required=True)
    p.add_argument('--jobs', required=True)
    p.set_defaults(function=gate)
    p = commands.add_parser('replay')
    p.add_argument('--plan', required=True)
    p.add_argument('--base-binary')
    p.add_argument('--inputs', required=True)
    p.add_argument('--binary', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--summary', required=True)
    p.set_defaults(function=replay)
    args = parser.parse_args()
    try:
        args.function(args)
    except (ValueError, OSError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
