"""Explicit baseline evidence reuse for the driver's strict six-run gate.

No discovery, global cache, retries, acceptance exceptions or geometry changes.
An archive is an integrity-checked evidence directory, not a trusted verdict.
"""
import argparse
from decimal import Decimal
import json
import os
from pathlib import Path
import platform
import re
import resource
import shutil
import subprocess
import sys
import time


class EvidenceError(ValueError):
    """Evidence is incomplete, inconsistent or outside this strict gate."""


class BaselineDefect(EvidenceError):
    """Intact prior evidence fails the current validator; diagnosis is required."""


def require(condition, message):
    if not condition:
        raise EvidenceError(message)


def read(path):
    # Keep JSON numbers exact too: never round a published decimal interval.
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), parse_float=Decimal,
                      parse_constant=lambda _: require(False, 'nonfinite JSON'),
                      object_pairs_hook=pairs)


def known(value):
    if value is None or (isinstance(value, str) and
                         value.strip().lower() in ('', '?', 'unknown', 'unspecified', 'unavailable')):
        return False
    if isinstance(value, dict):
        return all(known(v) for v in value.values())
    if isinstance(value, list):
        return all(known(v) for v in value)
    return True


def build_identity(path, binary, driver):
    """Require a build-time attestation; an executable cannot reveal its flags."""
    data = json.loads(Path(path).read_text())
    require(known(data), 'unknown build provenance')
    require(set(data) == {'format', 'source', 'executable_sha256', 'build'}, 'build provenance fields')
    require(data['format'] == 'catallax-build-evidence/1', 'build provenance version')
    source = data['source']
    require(set(source) == {'commit', 'tree', 'dirty', 'submodules'}, 'source identity fields')
    require(source['dirty'] is False and isinstance(source['submodules'], list), 'unclean source identity')
    require(all(re.fullmatch('[0-9a-f]{40}', source[k]) for k in ('commit', 'tree')), 'source hashes')
    require(all(isinstance(s, dict) and set(s) == {'path', 'commit'} and
                isinstance(s['path'], str) and re.fullmatch('[0-9a-f]{40}', s['commit'])
                for s in source['submodules']), 'submodule identities')
    require(data['executable_sha256'] == driver.file_hash(binary), 'executable identity')
    build = data['build']
    require(set(build) == {'toolchain', 'dependencies', 'options'}, 'build identity fields')
    require(set(build['toolchain']) == {'rustc', 'cargo', 'linker'}, 'toolchain identity fields')
    require(all(isinstance(v, str) and v for v in build['toolchain'].values()), 'toolchain identity')
    require(set(build['dependencies']) == {'cargo_lock_sha256', 'resolved_sha256', 'native_sha256'}, 'dependency identity fields')
    options = build['options']
    require(set(options) == {'command', 'env_sha256', 'config_sha256', 'target', 'profile', 'features'}, 'build options fields')
    require(isinstance(options['command'], list) and options['command'] and
            all(isinstance(v, str) for v in options['command']) and
            isinstance(options['features'], list) and
            all(isinstance(v, str) for v in options['features']) and
            all(isinstance(options[k], str) for k in ('target', 'profile')), 'build options')
    hashes = [*build['dependencies'].values(), options['env_sha256'], options['config_sha256']]
    require(all(isinstance(v, str) and re.fullmatch('[0-9a-f]{64}', v) for v in hashes), 'build hashes')
    return data


def effective_inputs(args, driver):
    return [dict(item, max_work=item.get('max_work', args.max_work))
            for item in driver.load(args.inputs, args.source, args.adapter)]


def reuse_key(args, driver):
    identity = build_identity(args.baseline_identity, args.baseline, driver)
    about = json.loads(subprocess.check_output([args.baseline, 'about'], text=True))
    kernel = json.loads(about['kernel']) if isinstance(about['kernel'], str) else about['kernel']
    require(about['request'] == driver.CONTRACT and kernel['commit'] == identity['source']['commit'],
            'executable/source provenance disagreement')
    linked = subprocess.run(['ldd', args.baseline], capture_output=True, text=True)
    libraries = sorted(set(re.findall(r'(?:=>\s+|^\s*)(/\S+)', linked.stdout, re.MULTILINE)))
    require(libraries and linked.returncode == 0 or
            'statically linked' in linked.stdout + linked.stderr, 'unknown runtime dependencies')
    here = Path(__file__).resolve()
    cores = [int(c) for c in args.cores.split(',')] if args.cores else []
    return {
        'format': 'catallax-baseline-key/1', 'baseline': identity,
        'ordered_requests': effective_inputs(args, driver),
        'interpretation': {'source': args.source, 'default_max_work': args.max_work,
                           'adapter_sha256': driver.file_hash(args.adapter) if args.adapter else 'builtin'},
        'semantics': {'driver': driver.file_hash(driver.__file__), 'comparator': driver.file_hash(here),
                      'request': driver.CONTRACT, 'accounting': 'complete-request; unchanged per-request allowance and embedded caps'},
        'worker_policy': {'counts': [1, 2, 4], 'cores': cores, 'dispatch': 'ordered unique keys round-robin; duplicate occurrences retained'},
        'environment': {'variables_sha256': driver.digest(dict(os.environ)),
                        'runtime_libraries': {name: driver.file_hash(name) for name in libraries},
                        'uname': list(platform.uname()), 'python': sys.version,
                        'python_sha256': driver.file_hash(sys.executable),
                        'cpu_sha256': driver.digest([line for line in Path('/proc/cpuinfo').read_text().splitlines()
                                                   if line.partition(':')[0].strip() not in ('cpu MHz', 'bogomips')]),
                        'affinity': sorted(os.sched_getaffinity(0)),
                        'scheduler': [os.sched_getscheduler(0), os.sched_getparam(0).sched_priority,
                                      os.getpriority(os.PRIO_PROCESS, 0)],
                        'limits': [list(resource.getrlimit(k)) for k in
                                   (resource.RLIMIT_AS, resource.RLIMIT_STACK, resource.RLIMIT_CPU, resource.RLIMIT_NOFILE)]}}


def inventory(root, driver):
    root = Path(root)
    result = {}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'archive contains a symlink')
        if path.is_file() and path != root / 'archive.json':
            result[path.relative_to(root).as_posix()] = driver.file_hash(path)
    return result


def without_time(value):
    if isinstance(value, dict):
        return {k: without_time(v) for k, v in value.items() if k != 'wall_ms'}
    if isinstance(value, list):
        return [without_time(v) for v in value]
    return value


def exact_manifest(manifest, driver):
    return driver.comparable(manifest)


def audit_run(path, items, binary_hash, workers, cores, driver, evidence_kind='fresh'):
    """Check every occurrence and every file, including non-time log content."""
    path = Path(path)
    report = read(path)
    require(report['format'] == 'catallax-run/1', 'run version')
    require(report['evidence']['kind'] == evidence_kind, 'fresh/reused provenance')
    config = report['config']
    require(config['binary_sha256'] == binary_hash and config['contract'] == driver.CONTRACT, 'run identity')
    execution = report['execution']
    require(execution['fresh'] is True and execution['workers'] == workers and execution['cores'] == cores, 'fresh worker policy')
    # driver hashes the loaded item before resolving a default allowance.
    loaded = json.loads((path.parent / 'inputs.json').read_text())
    require([dict(x, max_work=x.get('max_work', config['max_work'])) for x in loaded] == items, 'ordered inputs')
    expected_keys = [driver.digest({**config, **item}) for item in loaded]
    expected_workers = min(workers, len(set(expected_keys)))
    require(execution['worker_exits'] == {str(i): 0 for i in range(expected_workers)}, 'worker exits')
    saved_config = read(path.parent / 'config.json')
    require(all(saved_config[k] == v for k, v in config.items()) and
            saved_config['workers'] == workers and saved_config['cores'] == cores, 'saved config')
    require([c['key'] for c in report['cases']] == expected_keys, 'case completeness/order/duplicates')
    deterministic = []
    for case, item, key in zip(report['cases'], items, expected_keys):
        require(case['cache_hit'] is False, 'cached candidate or baseline evidence')
        directory = path.parent / 'cases' / key
        require((path.parent / case['path']).resolve() == directory.resolve(), 'case outside owned evidence')
        metadata = read(directory / 'artifact.json')
        require(metadata['key'] == key and
                set(metadata['files']) == {p.name for p in directory.iterdir() if p.name != 'artifact.json'}, 'artifact completeness')
        require(all(re.fullmatch('[a-zA-Z0-9_.-]+', name) and
                    driver.file_hash(directory / name) == sha for name, sha in metadata['files'].items()), 'artifact integrity')
        require({'request.json', 'response.json', 'manifest.json', 'operations.json', 'question.json'} <= set(metadata['files']), 'missing exchange evidence')
        request = read(directory / 'request.json')
        request.pop('out')
        require(request == dict(item, format=driver.CONTRACT), 'request interpretation/allowance')
        portable = read(directory / 'question.json')
        require(portable == dict(item, format='catallax-question/1'), 'portable question')
        response = read(directory / 'response.json')
        manifest = read(directory / 'manifest.json')
        log = read(directory / 'operations.json')
        require(response == case['response'] and response['manifest'] == manifest and
                response['operation_log'] == log, 'response/manifest/log disagreement')
        require(response['format'] == 'catallax-response/1' and
                manifest['status'] in ('answered', 'refused') and
                log['format'] == 'catallax-operation-log/1' and log['status'] == manifest['status'], 'invalid outcome')
        require((manifest.get('comparison') or {}).get('verdict') != 'disagree' and
                (manifest.get('name_edit') or {}).get('status') not in ('failed', 'invalid'), 'failed comparison/name-edit gate')
        require(all(log[k] == manifest[k] for k in ('refusal', 'marks')) and
                isinstance(log['operations'], list), 'operation log completeness')
        work = response['work']
        require(type(work['used']) is int and 0 <= work['used'] <= item['max_work'] and
                work['limit'] == item['max_work'] and type(work['exhausted']) is bool, 'work accounting')
        for sequence, operation in enumerate(log['operations']):
            counted = operation['work']
            require(operation['sequence'] == sequence and
                    all(type(counted[k]) is int for k in ('before', 'after', 'used')) and
                    0 <= counted['before'] <= counted['after'] <= work['used'] and
                    counted['used'] == counted['after'] - counted['before'], 'operation work accounting')
        artifacts = {name: driver.file_hash(directory / name) for name in
                     ('answer.step', 'answer.igs', 'display.obj', 'print.stl') if (directory / name).exists()}
        require(manifest['status'] == 'answered' or not artifacts, 'stale refused artifacts')
        require(bool(manifest.get('step')) == ('answer.step' in artifacts), 'STEP presence')
        if manifest.get('step'):
            require(manifest['step']['file'] == 'answer.step', 'STEP artifact name')
        for field, filename in (('iges', 'answer.igs'), ('printing', 'print.stl')):
            require(bool(manifest.get(field)) == (filename in artifacts), 'export presence')
            if manifest.get(field):
                require(manifest[field]['file'] == filename, 'export artifact name')
        if manifest.get('iges'):
            require(manifest['iges']['sha256'] == artifacts['answer.igs'] and
                    manifest['iges']['bytes'] == (directory / 'answer.igs').stat().st_size, 'IGES artifact identity')
        display = manifest.get('face_measurements', {}).get('display', {})
        require(('display.obj' in artifacts) == (bool(display) or 'print.stl' in artifacts), 'display presence')
        if display:
            require(display['file'] == 'display.obj', 'display artifact name')
        deterministic.append({'manifest': exact_manifest(manifest, driver), 'artifacts': artifacts,
                              'work': work, 'log': without_time(log)})
    return deterministic


def gates(root, items, hashes, cores, driver, reused_baseline=False):
    observed = {}
    for role in ('baseline', 'candidate'):
        for workers in (1, 2, 4):
            observed[role, workers] = audit_run(Path(root) / f'{role}{workers}/run.json', items,
                hashes[role], workers, cores, driver,
                'reused' if role == 'baseline' and reused_baseline else 'fresh')
        require(observed[role, 1] == observed[role, 2] == observed[role, 4], 'within-build manifest/artifact/work/non-time-log determinism')
    # Conservative preservation: no changed export or changed certificate is
    # automatically accepted. Such a change needs the parent's separate audit.
    for workers in (1, 2, 4):
        left, right = observed['baseline', workers], observed['candidate', workers]
        require([{k: c[k] for k in ('manifest', 'artifacts')} for c in left] ==
                [{k: c[k] for k in ('manifest', 'artifacts')} for c in right], 'preservation/changed-export audit required')
        result = driver.compare(Path(root) / f'baseline{workers}/run.json', Path(root) / f'candidate{workers}/run.json')
        require(result['equal'], 'current strict comparator failed')


def validate_archive(root, key, driver):
    record = read(Path(root) / 'archive.json')
    require(record['format'] == 'catallax-baseline-archive/1' and record['origin'] == 'fresh-six', 'archive provenance')
    require(record['key'] == key, 'complete key changed or unknown')
    require(record['files'] == inventory(root, driver), 'archive inventory/integrity')
    require(record['hashes']['baseline'] == key['baseline']['executable_sha256'], 'archive executable')
    try:
        gates(root, key['ordered_requests'], record['hashes'], key['worker_policy']['cores'], driver)
    except (ValueError, TypeError, KeyError, OSError) as error:
        raise BaselineDefect(f'intact archive fails current validator: {error}; diagnose before requalification') from error
    return record


def qualify(args, driver):
    started = time.perf_counter()
    out = Path(args.out).resolve()
    require(not out.exists(), 'qualification output already exists; retain and diagnose, never resume a qualification')
    items = effective_inputs(args, driver)
    hashes = {role: driver.file_hash(getattr(args, role)) for role in ('baseline', 'candidate')}
    key, key_error = None, None
    try:
        key = reuse_key(args, driver)
    except (OSError, ValueError, TypeError, KeyError, subprocess.SubprocessError) as error:
        key_error = f'unknown/incomplete key: {type(error).__name__}'
    cores = [int(c) for c in args.cores.split(',')] if args.cores else []
    reuse, decision, defect = False, 'no archive selected', None
    audit_started = time.perf_counter()
    if args.fresh_audit:
        decision = 'explicit fresh audit'
    elif key_error:
        decision = key_error
    elif args.reuse_baseline:
        try:
            validate_archive(args.reuse_baseline, key, driver)
            reuse, decision = True, 'current comparator revalidated earlier fresh qualified six'
        except BaselineDefect as error:
            defect = str(error)
            decision = 'baseline/validator defect; stopped before executing jobs'
        except (OSError, ValueError, TypeError, KeyError) as error:
            decision = f'archive ineligible: {type(error).__name__}: {error}'
    audit_ms = (time.perf_counter() - audit_started) * 1000
    out.mkdir(parents=True)
    result = {'format': 'catallax-qualification/1', 'mode': 'REUSED-BASELINE' if reuse else 'FULL-SIX',
              'decision': decision, 'status': 'running', 'equal': False,
              'archive_audit_ms': audit_ms, 'runs': [], 'hashes': hashes,
              'provenance': {'archive_sha256': driver.file_hash(Path(args.reuse_baseline) / 'archive.json') if reuse else None}}
    driver.atomic(out / 'qualification.json', result)
    if defect:
        result.update(status='failed', error=defect,
                      action='retain evidence; diagnose and explicitly requalify; do not retry until green',
                      elapsed_ms=(time.perf_counter() - started) * 1000)
        driver.atomic(out / 'qualification.json', result)
        return result
    try:
        for role in ('baseline', 'candidate'):
            for workers in (1, 2, 4):
                name = f'{role}{workers}'
                destination = out / name
                run_started = time.perf_counter()
                if role == 'baseline' and reuse:
                    shutil.copytree(Path(args.reuse_baseline) / name, destination)
                else:
                    # Candidate has no cache/resume/reuse branch, including when
                    # both executable hashes happen to be identical.
                    driver.run(argparse.Namespace(inputs=args.inputs, source=args.source, adapter=args.adapter,
                               binary=getattr(args, role), out=str(destination), cache=None, cores=args.cores,
                               max_work=args.max_work, workers=workers, dry_run=False, fresh=True))
                report = json.loads((destination / 'run.json').read_text())
                report['evidence'] = {'kind': 'reused' if role == 'baseline' and reuse else 'fresh',
                                      'archive_sha256': result['provenance']['archive_sha256'] if role == 'baseline' and reuse else None}
                for case in report['cases']:
                    case['path'] = f"cases/{case['key']}"
                driver.atomic(destination / 'run.json', report)
                result['runs'].append({'role': role, 'workers': workers,
                    'evidence': 'reused' if role == 'baseline' and reuse else 'fresh',
                    'historical_elapsed_ms': report['elapsed_ms'] if role == 'baseline' and reuse else None,
                    'observed_elapsed_ms': (time.perf_counter() - run_started) * 1000})
                driver.atomic(out / 'qualification.json', result)
        require(hashes == {role: driver.file_hash(getattr(args, role)) for role in hashes}, 'executable changed during qualification')
        if key is not None:
            require(reuse_key(args, driver) == key, 'key changed during qualification')
        gate_started = time.perf_counter()
        gates(out, items, hashes, cores, driver, reused_baseline=reuse)
        result.update(status='qualified', equal=True, gate_audit_ms=(time.perf_counter() - gate_started) * 1000)
        if not reuse and key is not None:
            # Only fresh six-run evidence can seed reuse. Reused descendants
            # cannot launder historical baseline evidence into a fresh archive.
            driver.atomic(out / 'archive.json', {'format': 'catallax-baseline-archive/1',
                         'origin': 'fresh-six', 'key': key, 'hashes': hashes, 'files': inventory(out, driver)})
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
        result.update(status='failed', equal=False, error=str(error),
                      action='retain evidence; diagnose and requalify; do not retry until green')
    result['elapsed_ms'] = (time.perf_counter() - started) * 1000
    driver.atomic(out / 'qualification.json', result)
    # qualification.json is a mutable progress report, so seal it only after
    # the final write. No stored equal flag is used by validate_archive.
    if (out / 'archive.json').exists():
        record = json.loads((out / 'archive.json').read_text())
        record['files'] = inventory(out, driver)
        driver.atomic(out / 'archive.json', record)
    return result
