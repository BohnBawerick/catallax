"""Development-only released-tool comparison for the fixed target scoreboard."""
import copy
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tarfile
import time
import urllib.request

DATA = Path(__file__).with_name('head-to-head')


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    temporary.replace(path)


def fetch_reference(out):
    pin = json.loads((DATA / 'reference.json').read_text())
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    archive = out / Path(pin['url']).name
    if not archive.exists():
        with urllib.request.urlopen(pin['url']) as source, archive.open('xb') as dest:
            while block := source.read(1024 * 1024):
                dest.write(block)
    if sha(archive) != pin['archive_sha256']:
        raise ValueError('reference archive checksum mismatch; refusing extraction')
    binary = out / 'solidcraft-cli'
    with tarfile.open(archive) as bundle:
        member = bundle.getmember(pin['member'])
        if not member.isfile():
            raise ValueError('reference executable is not a regular archive member')
        content = bundle.extractfile(member).read()
    if hashlib.sha256(content).hexdigest() != pin['binary_sha256']:
        raise ValueError('reference binary checksum mismatch; refusing execution')
    binary.write_bytes(content)
    binary.chmod(0o755)
    return {'binary': str(binary), **pin}


def process(command, directory, core, timeout):
    directory.mkdir(parents=True, exist_ok=True)
    command = ['taskset', '-c', str(core), *map(str, command)]
    save(directory / 'command.json', command)
    started = time.perf_counter()
    env = dict(os.environ, RAYON_NUM_THREADS='1')
    with (directory / 'stdout').open('w') as stdout, (directory / 'stderr').open('w') as stderr:
        try:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, env=env, timeout=timeout)
            status, code = 'completed', result.returncode
        except subprocess.TimeoutExpired:
            status, code = 'timeout', None
    row = {'status': status, 'exit': code, 'process_ms': (time.perf_counter()-started)*1000}
    save(directory / 'process.json', row)
    return row


def edited_recipe(recipe, edit):
    result = copy.deepcopy(recipe)
    target = next(step for step in result['steps'] if step['id'] == edit['step'])
    for key in edit['path'][:-1]:
        target = target[key]
    if not isinstance(target[edit['path'][-1]], str):
        raise ValueError('edit does not replace an exact scalar')
    target[edit['path'][-1]] = edit['value']
    return result


def instruction(command, **params):
    return {'command': command, 'params': params}


def solidcraft(case, binary, directory, args):
    body = case['reference']['body']
    result = {'status': 'completed', 'variants': {}, 'work': None,
              'work_reason': 'released CLI has no deterministic work counter',
              'consumer_mode': 'separate commands reopening the saved design; internal command times exclude reopen and save'}
    previous = None
    for label in ('construction', 'edit'):
        variant = {}
        result['variants'][label] = variant
        if label == 'edit' and previous is None:
            variant['build'] = {'status': 'not-run', 'reason': 'successful original design required'}
            continue
        design = directory / f'{label}.solidcraft'
        def execute(phase, commands, source=None, save_design=None):
            dest = directory / label / phase
            dest.mkdir(parents=True, exist_ok=True)
            save(dest / 'script.json', {'commands': commands})
            command = [binary, 'run', dest / 'script.json', '--results', '--timings', '--profile']
            if source is not None:
                command.extend(['--in', source])
            if save_design is not None:
                command.extend(['--save', save_design])
            execution = process(command, dest, args.core, args.timeout)
            if execution['exit'] != 0:
                return {'status': 'incomplete' if execution['status'] == 'timeout' else 'command-failed',
                        'process': execution, 'wall_ms': None, 'work': None,
                        'refusal_upper_bound_ms': execution['process_ms'] if execution['status'] != 'timeout' else None,
                        'refusal_scope': 'whole process including reopen and preceding commands; per-command timer unavailable on error',
                        'stderr': (dest / 'stderr').read_text()}
            payload = json.loads((dest / 'stdout').read_text())
            return {'status': 'answered', 'process': execution,
                    'wall_ms': sum(t['ms'] for t in payload['timings']),
                    'embedded_display_ms': sum(t['profile']['tessellate_ms'] for t in payload['timings']),
                    'values': payload['results'], 'timings': payload['timings'], 'work': None, 'refusal_ms': None}
        commands = case['reference']['commands'] if label == 'construction' else [case['reference']['edit']]
        variant['build'] = execute('build', commands, previous, design)
        build = variant['build']
        if build['status'] != 'answered':
            previous = None
            continue
        build['operations_rebuilt'] = (build['values'][0].get('recomputed') if label == 'edit'
                                        else sum('feature' in value for value in build['values']))
        previous = design
        for phase, command in (
            ('display', instruction('file.export', path=str(directory / f'{label}.obj'), format='obj', bodies=[body])),
            ('measurement', instruction('inspect.measure', bodies=[body])),
            ('names_faces', instruction('model.faces', body=body)),
            ('names_edges', instruction('model.edges', body=body)),
            ('step', instruction('file.export', path=str(directory / f'{label}.step'), format='step', bodies=[body]))):
            variant[phase] = execute(phase, [command], design)
        if all(variant[p]['status'] == 'answered' for p in ('names_faces', 'names_edges')):
            faces = variant['names_faces']['values'][0]['faces']
            edges = variant['names_edges']['values'][0]['edges']
            variant['names'] = {'faces': [f.get('name') for f in faces], 'edges': [e.get('name') for e in edges]}
            variant['step']['expected_faces'] = len(faces)
        if variant['measurement']['status'] == 'answered':
            measured = variant['measurement']['values'][0]['bodies']
            if len(measured) != 1:
                raise ValueError('reference must measure exactly the complete target body')
            m = measured[0]
            variant['measurement']['claimed'] = {'volume': m['volume_mm3'], 'area': m['area_mm2'], 'centroid': m['center_of_mass']}
    old, new = (result['variants'][v].get('names') for v in ('construction', 'edit'))
    if old is not None and new is not None:
        missing = {kind: sorted(set(old[kind]) - set(new[kind]), key=str) for kind in ('faces', 'edges')}
        unnamed = sum(n is None or n == '' for names in (old, new) for kind in names for n in names[kind])
        result['edit_survival'] = {'status': 'failed' if unnamed or any(missing.values()) else 'passed',
                                 'missing': missing, 'unnamed': unnamed,
                                 'coverage': 'every original result face/edge string name remains present; not a semantic correspondence proof'}
    else:
        result['edit_survival'] = {'status': 'not-run', 'reason': 'both complete result name lists required'}
    return result


def catallax(case, recipe, binary, directory, args):
    save(directory / 'request.json', {'recipe': recipe, 'edited_recipe': edited_recipe(recipe, case['edit']), 'max_work': args.max_work})
    execution = process([binary, directory / 'request.json', directory], directory, args.core, args.timeout)
    partial = directory / 'phases.json'
    result = json.loads(partial.read_text()) if partial.exists() else {}
    result['process'] = execution
    result['status'] = 'completed' if execution['exit'] == 0 else 'incomplete' if execution['status'] == 'timeout' else 'failed-rebuild'
    for variant in result.get('variants', {}).values():
        if 'values' in variant.get('measurement', {}):
            variant['measurement']['claimed'] = variant['measurement']['values']
    return result


def compare_properties(claimed, oracle, bounded):
    if claimed is None or oracle is None:
        return {'status': 'not-run', 'reason': 'both requested properties and independent values required'}
    wrong = []
    for property_name in ('volume', 'area', 'centroid'):
        actual = claimed.get(property_name)
        expected = oracle[property_name]
        actuals = actual if property_name == 'centroid' else [actual]
        expecteds = expected if property_name == 'centroid' else [expected]
        if actuals is None:
            wrong.append({'property': property_name, 'reason': 'missing claimed value'})
            continue
        if len(actuals) != len(expecteds):
            raise ValueError('malformed centroid')
        for index, (value, truth) in enumerate(zip(actuals, expecteds)):
            if bounded:
                lo, hi = map(lambda n: float(Fraction(n)), value)
            else:
                lo = hi = value
            allowance = max(1e-6, abs(truth)*1e-8)
            if not all(map(__import__('math').isfinite, (lo, hi, truth))) or lo > hi or truth < lo-allowance or truth > hi+allowance:
                wrong.append({'property': property_name, 'component': index if property_name == 'centroid' else None,
                              'claimed': value, 'independent': truth, 'diagnostic_allowance': allowance})
    return {'status': 'failed' if wrong else 'passed', 'wrong_answers': wrong,
            'scope': 'numeric independent-oracle check; passing is not a certificate'}


def correctness(case_path, case, engine, result, directory, args):
    for label, variant in result.get('variants', {}).items():
        step = directory / label / 'answer.step' if engine == 'catallax' else directory / f'{label}.step'
        checker_dir = directory / f'{label}-check'
        faces = variant.get('step', {}).get('expected_faces', -1)
        p = process([args.occt_python, DATA / 'check.py', case_path, label, step if step.exists() else '-', faces], checker_dir, args.core, args.timeout)
        if p['exit'] == 0:
            check = json.loads((checker_dir / 'stdout').read_text().splitlines()[-1])
        else:
            check = {'status': 'incomplete', 'process': p}
        check['claimed_properties'] = compare_properties(variant.get('measurement', {}).get('claimed'), check.get('oracle'), engine == 'catallax')
        if step.exists():
            check['step_sha256'] = sha(step)
        variant['correctness'] = check
    save(directory / 'row.json', result)


def benchmark(args):
    if args.core not in os.sched_getaffinity(0):
        raise ValueError('core is outside allowed affinity')
    pin = json.loads((DATA / 'reference.json').read_text())
    if sha(args.reference) != pin['binary_sha256']:
        raise ValueError('reference executable does not match the pinned released binary')
    root = Path(args.inputs).resolve()
    cases = sorted((DATA / 'cases').rglob('*.json'))
    if args.part:
        cases = [p for p in cases if p.relative_to(DATA / 'cases').as_posix() in args.part]
        if len(cases) != len(set(args.part)):
            raise ValueError('--part must name exact committed case filenames')
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'head-to-head.json').exists():
        raise ValueError('use a fresh output directory; never mix measured runs')
    report = {'format': 'catallax-head-to-head/1', 'status': 'running', 'reference': pin,
              'binary_sha256': sha(args.binary), 'host': {'system': platform.platform(), 'cpu': platform.processor(),
              'affinity': sorted(os.sched_getaffinity(0))},
              'timing': {'core': args.core, 'order': 'alternating AB/BA', 'warmups': 1, 'measured_pairs': args.repeats,
                 'threads': 1, 'max_work_per_catallax_phase': args.max_work, 'timeout_seconds_per_engine': args.timeout,
                 'claim': 'shared-machine pinned-core observations; startup/checker/export costs not silently subtracted'}, 'pairs': []}
    save(out / 'head-to-head.json', report)
    try:
        for path in cases:
            case = json.loads(path.read_text())
            recipe_path = root / case['recipe']
            if sha(recipe_path) != case['recipe_sha256']:
                raise ValueError(f'{case["recipe"]}: recipe changed since reference script was authored')
            recipe = json.loads(recipe_path.read_text())
            for index in range(args.repeats + 1):
                order = ['catallax', 'solidcraft'] if index % 2 == 0 else ['solidcraft', 'catallax']
                pair = {'part': case['recipe'], 'sample': index, 'warmup': index == 0,
                        'order': order, 'load_average_before': os.getloadavg(), 'engines': {}}
                for engine in order:
                    directory = out / 'cases' / case['recipe'].removesuffix('.json') / str(index) / engine
                    directory.mkdir(parents=True, exist_ok=True)
                    if engine == 'catallax':
                        result = catallax(case, recipe, Path(args.binary).resolve(), directory, args)
                    else:
                        result = solidcraft(case, Path(args.reference).resolve(), directory, args)
                    pair['engines'][engine] = result
                pair['load_average_after'] = os.getloadavg()
                # Checks are outside the matched timing window and retained even for refusals.
                for engine, result in pair['engines'].items():
                    directory = out / 'cases' / case['recipe'].removesuffix('.json') / str(index) / engine
                    correctness(path.resolve(), case, engine, result, directory, args)
                report['pairs'].append(pair)
                save(out / 'head-to-head.json', report)
    except BaseException:
        report['status'] = 'interrupted'
        save(out / 'head-to-head.json', report)
        raise
    report['status'] = 'completed'
    report['summary'] = summarize(report)
    measured = [r for p in report['pairs'] if not p['warmup'] for r in p['engines'].values()]
    report['correctness_failures'] = sum(
        check.get('status') == 'failed'
        for result in measured for variant in result.get('variants', {}).values()
        for check in (variant.get('correctness', {}),
                      variant.get('correctness', {}).get('claimed_properties', {}),
                      variant.get('correctness', {}).get('step', {}),
                      variant.get('correctness', {}).get('membership', {})))
    report['correctness_failures'] += sum(r.get('edit_survival', {}).get('status') == 'failed' for r in measured)
    save(out / 'head-to-head.json', report)
    (out / 'report.md').write_text(markdown(report))
    return report


def summarize(report):
    rows = []
    for part in sorted({p['part'] for p in report['pairs']}):
        for engine in ('catallax', 'solidcraft'):
            pairs = [p['engines'][engine] for p in report['pairs'] if p['part'] == part and not p['warmup']]
            for label in ('construction', 'edit'):
                variants = [p.get('variants', {}).get(label, {}) for p in pairs]
                row = {'part': part, 'engine': engine, 'phase': label, 'samples': len(pairs)}
                for phase in ('build', 'display', 'measurement'):
                    values = [v.get(phase, {}).get('wall_ms') for v in variants]
                    row[phase+'_median_ms'] = statistics.median(values) if values and all(v is not None for v in values) else None
                waits = []
                queries = []
                for variant in variants:
                    build, display, measure = (variant.get(k, {}) for k in ('build', 'display', 'measurement'))
                    ready = all(p.get('status') == 'answered' and p.get('wall_ms') is not None for p in (build, display))
                    waits.append(build['wall_ms'] + display['wall_ms'] if ready else None)
                    queries.append(waits[-1] + measure['wall_ms'] if ready and measure.get('status') == 'answered' else None)
                row['build_and_display_budget_ms'] = 1000
                row['build_and_display_median_ms'] = statistics.median(waits) if waits and all(v is not None for v in waits) else None
                row['build_and_display_budget_passes'] = sum(v is not None and v <= 1000 for v in waits)
                row['with_requested_measurement_median_ms'] = statistics.median(queries) if queries and all(v is not None for v in queries) else None
                row['work'] = [v.get('build', {}).get('work') for v in variants]
                row['operations_rebuilt'] = [v.get('build', {}).get('operations_rebuilt') for v in variants]
                row['refusal_ms'] = [v.get('build', {}).get('refusal_ms') for v in variants]
                row['refusal_upper_bound_ms'] = [v.get('build', {}).get('refusal_upper_bound_ms') for v in variants]
                row['failed_rebuilds'] = sum(v.get('build', {}).get('status') in ('refused', 'work-exhausted', 'command-failed') for v in variants)
                row['incomplete'] = sum(p['status'] == 'incomplete' or v.get('build', {}).get('status') in ('incomplete', 'not-run') for p, v in zip(pairs, variants))
                row['wrong_properties'] = sum(len(v.get('correctness', {}).get('claimed_properties', {}).get('wrong_answers', [])) for v in variants)
                row['wrong_membership'] = sum(len(v.get('correctness', {}).get('membership', {}).get('wrong', [])) for v in variants)
                row['properties_checked'] = sum(v.get('correctness', {}).get('claimed_properties', {}).get('status') in ('passed','failed') for v in variants)
                row['membership_checked'] = sum(v.get('correctness', {}).get('membership', {}).get('status') in ('passed','failed') for v in variants)
                row['step_passed'] = sum(v.get('correctness', {}).get('step', {}).get('status') == 'passed' for v in variants)
                row['names_passed'] = sum(p.get('edit_survival', {}).get('status') == 'passed' for p in pairs) if label == 'edit' else None
                rows.append(row)
    return rows


def markdown(report):
    lines = ['# Target-part head-to-head', '',
             f'Reference release {report["reference"]["version"]}. One warmup and {report["timing"]["measured_pairs"]} alternating matched pairs per part on core {report["timing"]["core"]}.', '',
             'Shared-machine observations. See raw JSON for load averages, work, operations rebuilt, failures, STEP checks and every property discrepancy.', '',
             '| Part | Tool | Variant | Build ms | Display/export ms | Measurement ms | Failed rebuilds | Incomplete | Property errors / checked | Membership errors | STEP passes | Name passes |',
             '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |']
    def number(value):
        return 'unmeasured' if value is None else f'{value:.3f}'
    for row in report['summary']:
        lines.append('| ' + ' | '.join([row['part'], row['engine'], row['phase'],
            *[number(row[k+'_median_ms']) for k in ('build','display','measurement')],
            str(row['failed_rebuilds']), str(row['incomplete']), f'{row["wrong_properties"]} / {row["properties_checked"]}',
            f'{row["wrong_membership"]} / {row["membership_checked"]}', str(row['step_passed']), str(row['names_passed'])]) + ' |')
    return '\n'.join(lines) + '\n'
