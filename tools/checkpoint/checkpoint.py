"""Checkpoint runs in one command each: the edit-chain scoreboard and the census rerun.

    python3 -B tools/checkpoint/checkpoint.py scoreboard --tools T --sources S --selection L --out O
    python3 -B tools/checkpoint/checkpoint.py census --tools T --sources S --selection L --out O

`--tools` is the comparison install from head-to-head/prepare.py; it is prepared on
first use. `--sources` and `--selection` are the private census corpus and its
[model_id, member] list. Every output goes to `--out`, outside the checkout.
Both commands build their own binaries from the clean current checkout.

The scoreboard runs P0 to P5 and ten real census parts through Catallax,
OpenCASCADE (CadQuery) and SMLib (smlibpy) with the existing head-to-head
driver, then reports one definition of done for all three kernels. The census
reruns import admission on all 3,000 files, then exports every full admission
and checks it in OCCT, and reports post-builder validity beside OCCT's own
healed-valid count.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / '.claude/skills/verify-catallax/scripts'
RESEARCH = ROOT / 'tools/research'
CONTROLS = ROOT / 'docs/research/step-census/retained-controls.json'
HASHES = ROOT / 'docs/research/step-census/source-hashes.json'
PREVIOUS = ROOT / 'docs/research/step-census-circular-pcurves/results.csv'
ENGINES = ['catallax', 'cadquery', 'smlibpy']
KERNEL_NAMES = {'catallax': 'Catallax', 'cadquery': 'OpenCASCADE (CadQuery 2.6.1)',
                'smlibpy': 'SMLib (smlibpy 0.4.0)'}
# Published reference capture: STEPControl_Reader, default ShapeFix, BRepCheck.
OCCT_HEALED_VALID = 2701
DENOMINATOR = 3000
IDENTITY = {'rotation': ['1', '0', '0', '0'], 'translation': ['0', '0', '0']}
# Exact homogeneous quaternions taking the stock cylinder's z axis to x, y and z.
AXIS_ROTATION = {0: ['1', '0', '1', '0'], 1: ['1', '-1', '0', '0'], 2: ['1', '0', '0', '0']}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args], text=True).strip()


def fresh(out):
    out = Path(out).resolve()
    if out.is_relative_to(ROOT):
        raise ValueError('checkpoint output belongs outside the checkout')
    out.mkdir(parents=True)
    return out


def runtime(tools):
    """Prepare the pinned comparison tools once; later runs verify the receipt."""
    path = Path(tools).resolve() / 'runtime.json'
    if not path.exists():
        subprocess.run([sys.executable, '-B', str(SCRIPTS / 'head-to-head/prepare.py'),
                        '--out', str(path.parent)], check=True)
    return path


def cargo_example(package, example, log):
    command = ['cargo', 'build', '--release', '-p', package, '--example', example,
               '--message-format=json']
    with Path(log).open('xb') as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream, check=True, timeout=3600)
    found = [row['executable'] for line in Path(log).read_text().splitlines()
             if (row := json.loads(line))['reason'] == 'compiler-artifact'
             and row['target']['name'] == example and row.get('executable')]
    assert len(found) == 1, found
    return Path(found[0])


def selected_sources(sources, selection):
    members = dict(json.loads(Path(selection).read_text()))
    root = Path(sources).resolve()
    paths = {}
    for key, member in members.items():
        path = (root / member).resolve()
        assert path.is_relative_to(root), member
        paths[key] = path
    return paths


def census_cases(paths, checker, directory):
    """Write one private case per retained real part: import, drill, edit the radius."""
    directory.mkdir()
    for row in json.loads(CONTROLS.read_text())['rows']:
        key, expected = row['model_id'], row['source_sha256']
        path = paths[key]
        if digest(path) != expected:
            raise ValueError(f'{key}: source bytes differ from the retained control')
        cut = json.loads(subprocess.check_output(
            [checker, '-B', str(ROOT / 'tools/checkpoint/cutter.py'), str(path)], text=True))
        axis = cut['axis'].index(1)
        recipe = {'format': 'catallax-recipe/3', 'name': f'ABC {key} drilled',
                  'steps': [{'id': 'imported', 'op': 'import_step',
                             'file': {'name': 'part.step', 'sha256': expected}, 'pose': IDENTITY},
                            {'id': 'cutter', 'op': 'cylinder', 'radius': cut['radius'],
                             'height': cut['height'],
                             'pose': {'rotation': AXIS_ROTATION[axis], 'translation': cut['base']}},
                            {'id': 'drilled', 'op': 'boolean', 'kind': 'subtract',
                             'target': 'imported', 'tool': 'cutter'}],
                  'result': 'drilled'}
        save(directory / f'{key}.json', {
            'format': 'catallax-head-to-head-case/1', 'recipe': f'census/{key}.json',
            'edit': {'step': 'cutter', 'path': ['radius'], 'value': cut['edited_radius']},
            'census': {'model_id': key, 'source': str(path), 'source_sha256': expected,
                       'cutter': cut, 'recipe': recipe}})


# A drill or radius change counts only if the OCCT-measured volume drops by more
# than this fraction of the source volume; equal volumes mean nothing was cut.
REMOVAL = 1e-6


def outcome(result, census=False):
    """One definition of done for every kernel, from the head-to-head row."""
    for label in ('construction', 'edit'):
        variant = result.get('variants', {}).get(label, {})
        build = variant.get('build', {})
        if build.get('status') != 'answered':
            detail = build.get('error') or build.get('detail') or build.get('reason') or ''
            return False, f'{label} build {build.get("status", "not-run")}', detail
        if result.get('process', {}).get('status') == 'timeout':
            for phase in ('display', 'measurement', 'step'):
                if 'status' not in variant.get(phase, {}):
                    return False, f'{label} {phase} timed out', 'external per-engine timeout; not a refusal'
        check = variant.get('correctness', {})
        step = check.get('step', {})
        if variant.get('step', {}).get('status') not in (None, 'answered'):
            return False, f'{label} STEP export {variant["step"]["status"]}', variant['step'].get('error') or ''
        if step.get('status') != 'passed':
            reason = ('not valid' if step.get('valid') is False else
                      'face count differs' if step.get('face_count_matches') is False else
                      step.get('reason') or step.get('status', 'not-run'))
            return False, f'{label} STEP {reason}', ''
        wrong = (len(check.get('claimed_properties', {}).get('wrong_answers', []))
                 + len(check.get('membership', {}).get('wrong', [])))
        if wrong:
            return False, f'{label} {wrong} wrong checked values', ''
    if census:
        before = result['variants']['construction']['correctness']['material']
        after = result['variants']['edit']['correctness']['material']
        source = before['source_volume_mm3']
        steps = (('drill', source, before['volume_mm3']),
                 ('radius change', before['volume_mm3'], after['volume_mm3']))
        for step, old, new in steps:
            if new - old > REMOVAL * source:
                return False, f'{step} raised the volume', \
                    f'a subtraction cannot add material; OCCT volume {old!r} -> {new!r} mm3'
            if not old - new > REMOVAL * source:
                return False, f'attempted subtraction: {step} removed no material', \
                    f'OCCT volume {old!r} -> {new!r} mm3'
    return True, 'done', ''


def scoreboard_tables(report):
    measured = [pair for pair in report['pairs'] if not pair['warmup']]
    parts = sorted({pair['part'] for pair in measured}, key=lambda p: (p.startswith('census/'), p))
    rows = []
    for part in parts:
        row = {'part': part, 'engines': {}}
        for engine in report['engines']:
            samples = [pair['engines'][engine] for pair in measured if pair['part'] == part]
            verdicts = [outcome(sample, part.startswith('census/')) for sample in samples]
            done = bool(verdicts) and all(v[0] for v in verdicts)
            first = next((v for v in verdicts if not v[0]), verdicts[0] if verdicts else (False, 'not run', ''))
            entry = {'done': done, 'summary': 'done' if done else first[1], 'detail': first[2]}
            if engine == 'catallax' and samples:
                variants = samples[-1].get('variants', {})
                entry['work'] = {label: variants.get(label, {}).get('build', {}).get('work')
                                 for label in ('construction', 'edit')}
                entry['names'] = samples[-1].get('edit_survival', {}).get('status')
            row['engines'][engine] = entry
        rows.append(row)
    return rows


def scoreboard_markdown(report, rows, seconds):
    engines = report['engines']
    lines = ['# Edit-chain scoreboard', '',
             f'Source `{report["source_revision"]}`, {report["host"]["cpu_model"]}, '
             f'{report["timing"]["measured_pairs"]} measured sample(s) after one warmup, '
             f'{report["timing"]["max_work_per_catallax_phase"]:,} work per Catallax phase, '
             f'{report["timing"]["timeout_seconds_per_engine"]:.0f} s per engine. '
             f'Wall time {seconds / 60:.1f} min.', '',
             'Done: the original and the edited target both build, both STEP exports open '
             'in OCCT 8.0.1.1 as one valid closed solid with the kernel\'s own face count, and no '
             'independent check found a wrong value. On a census part the drill and the radius '
             'change must each lower the OCCT-measured volume.', '',
             '| Target | ' + ' | '.join(KERNEL_NAMES[e] for e in engines) + ' | Catallax work (build, edit) | Catallax names |',
             '| --- |' + ' --- |' * len(engines) + ' ---: | --- |']
    for row in rows:
        work = row['engines'].get('catallax', {}).get('work') or {}
        def fmt(value):
            return 'n/a' if value is None else f'{value:,}'
        lines.append('| ' + ' | '.join([row['part'].removesuffix('.json'),
            *[row['engines'][e]['summary'] for e in engines],
            f'{fmt(work.get("construction"))}, {fmt(work.get("edit"))}',
            str(row['engines'].get('catallax', {}).get('names'))]) + ' |')
    lines += ['', '| Kernel | Done, P0 to P5 (13) | Done, census parts (10) |', '| --- | ---: | ---: |']
    for engine in engines:
        synthetic = sum(r['engines'][engine]['done'] for r in rows if not r['part'].startswith('census/'))
        real = sum(r['engines'][engine]['done'] for r in rows if r['part'].startswith('census/'))
        lines.append(f'| {KERNEL_NAMES[engine]} | {synthetic} | {real} |')
    details = [(r['part'], e, r['engines'][e]['detail']) for r in rows for e in engines
               if r['engines'][e]['detail']]
    if details:
        lines += ['', '## First refusal per target', '']
        lines += [f'- {part} / {engine}: {detail.splitlines()[0][:300]}' for part, engine, detail in details]
    return '\n'.join(lines) + '\n'


def scoreboard(args):
    out = fresh(args.out)
    started = time.monotonic()
    path = runtime(args.tools)
    receipt = json.loads(path.read_text())
    binary = cargo_example('catallax-cli', 'head_to_head', out / 'head-to-head-build.cargo.jsonl')
    cases = out / 'census-cases'
    census_cases(selected_sources(args.sources, args.selection), receipt['checker'], cases)
    sys.path.insert(0, str(SCRIPTS))
    from head_to_head import benchmark
    report = benchmark(argparse.Namespace(
        inputs=str(ROOT / 'parts'), binary=str(binary), reference=None, out=str(out / 'head-to-head'),
        core=args.core, max_work=args.max_work, repeats=args.repeats, timeout=args.timeout, part=args.part,
        occt_python=receipt['checker'], tool_runtime=str(path), machine_class=args.machine_class,
        engines=ENGINES, census_cases=str(cases)))
    rows = scoreboard_tables(report)
    seconds = time.monotonic() - started
    save(out / 'scoreboard.json', {'format': 'catallax-checkpoint-scoreboard/1',
                                   'source_revision': report['source_revision'],
                                   'wall_seconds': seconds, 'rows': rows})
    (out / 'scoreboard.md').write_text(scoreboard_markdown(report, rows, seconds))
    print((out / 'scoreboard.md').read_text())


def census(args):
    try:
        import OCP  # noqa: F401  The export check needs OCP in this interpreter.
    except ImportError:
        checker = json.loads(runtime(args.tools).read_text())['checker']
        os.execv(checker, [checker, '-B', str(Path(__file__).resolve()), *sys.argv[1:]])
    assert not git('status', '--porcelain'), 'the census binds a clean checkout'
    out = fresh(args.out)
    revision = git('rev-parse', 'HEAD')
    clock = {}
    started = time.monotonic()
    subprocess.run([sys.executable, '-B', str(RESEARCH / 'build_import_observer.py'),
                    '--record', str(out / 'import-build.json')], check=True)
    importer = json.loads((out / 'import-build.json').read_text())['artifact_path']
    clock['import_build'] = time.monotonic() - started
    started = time.monotonic()
    subprocess.run([sys.executable, '-B', str(RESEARCH / 'import_census.py'), '--binary', importer,
                    '--build-record', str(out / 'import-build.json'), '--sources', str(args.sources),
                    '--selection', str(args.selection), '--source-hashes', str(HASHES),
                    '--out', str(out / 'import'), '--workers', str(args.import_workers),
                    '--max-rss-mb', str(args.max_rss_mb)], check=False)
    clock['import_census'] = time.monotonic() - started
    admission = json.loads((out / 'import/summary.json').read_text())
    rows = {row['model_id']: row for row in json.loads((out / 'import/ordered-rows.json').read_text())}
    started = time.monotonic()
    subprocess.run([sys.executable, '-B', str(RESEARCH / 'export_validity_census.py'), 'build',
                    '--record', str(out / 'export-build.json'), '--kernel-revision', revision], check=True)
    exporter = Path(json.loads((out / 'export-build.json').read_text())['artifact_path'])
    clock['export_build'] = time.monotonic() - started
    sys.path.insert(0, str(RESEARCH))
    import export_validity_census as validity
    paths = selected_sources(args.sources, args.selection)
    admitted = sorted(key for key, row in rows.items() if row['outcome'] == 'full')
    (out / 'export/rows').mkdir(parents=True)
    (out / 'export/parts').mkdir()
    observed = {}
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.export_workers) as pool:
        jobs = {pool.submit(validity.observe, exporter, paths[key],
                            {'model_id': key, 'source_sha256': rows[key]['source_sha256'],
                             'bodies': rows[key]['bodies']},
                            out / 'export/parts' / key, 60, args.max_rss_mb): key for key in admitted}
        for job in as_completed(jobs):
            row = job.result()
            observed[row['model_id']] = row
            save(out / 'export/rows' / f'{row["model_id"]}.json', row)
            if len(observed) % 50 == 0:
                print(json.dumps({'exported': len(observed), 'of': len(admitted)}), flush=True)
    clock['export_census'] = time.monotonic() - started
    with PREVIOUS.open(newline='') as stream:
        csv.field_size_limit(64 * 1024 * 1024)
        before = {r['model_id'] for r in csv.DictReader(stream) if r['validity_outcome'] == 'raw_valid'}
    valid = {key for key, row in observed.items() if row['outcome'] == 'raw_valid'}
    # Runs with no verdict: aborted for memory or time, or a process error. A
    # previous pass among them was not measured, so it is neither kept nor lost.
    unmeasured = {'timeout': 'timed out', 'out_of_memory': 'out of memory',
                  'process_error': 'process error', 'missing_result': 'process error, no result'}
    not_covered = [{'model_id': k, 'stage': stage, 'kind': unmeasured[row['outcome']],
                    'detail': row.get('cause') or row.get('timeout_stage') or ''}
                   for stage, table in (('import', rows), ('export', observed))
                   for k, row in sorted(table.items()) if row['outcome'] in unmeasured]
    aborted = {entry['model_id'] for entry in not_covered}
    report = {
        'format': 'catallax-checkpoint-census/1', 'source_revision': revision,
        'import_binary_sha256': admission['binary_sha256'], 'export_binary_sha256': digest(exporter),
        'denominator': DENOMINATOR, 'complete': admission['complete'] and len(observed) == len(admitted)
        and 'process_error' not in Counter(r['outcome'] for r in observed.values()),
        'valid_after_builder': len(valid), 'valid_share_percent': f'{100 * len(valid) / DENOMINATOR:.2f}',
        'occt_healed_valid': OCCT_HEALED_VALID,
        'occt_healed_valid_share_percent': f'{100 * OCCT_HEALED_VALID / DENOMINATOR:.2f}',
        'previous_valid_after_builder': len(before),
        'gained_since_previous': sorted(valid - before),
        'lost_since_previous': sorted(before - valid - aborted),
        'previous_not_measured': sorted((before - valid) & aborted),
        'admission_outcomes': admission['outcomes'], 'import_categories': admission['import_categories'],
        'refusal_causes': admission['after_fix_causes'],
        'not_covered': not_covered,
        'max_rss_mb_per_observer': args.max_rss_mb,
        'export_outcomes': dict(Counter(r['outcome'] for r in observed.values())),
        'export_refusal_causes': dict(Counter(r.get('cause') for r in observed.values()
                                              if r['outcome'] == 'export_refused')),
        'wall_seconds': clock, 'import_workers': args.import_workers,
        'export_workers': args.export_workers, 'per_file_timeout_seconds': 60,
        'scope': 'a file counts when every imported body exports and passes BRepCheck after '
                 'OCCT builder construction, which can heal pcurves internally'}
    save(out / 'census.json', report)
    (out / 'census.md').write_text(census_markdown(report))
    print((out / 'census.md').read_text())


def census_markdown(r):
    total = sum(r['wall_seconds'].values())
    lines = ['# Census rerun', '',
             f'Source `{r["source_revision"]}`. Complete: {r["complete"]}. Wall time {total / 60:.1f} min '
             f'({", ".join(f"{k} {v / 60:.1f}" for k, v in r["wall_seconds"].items())}).', '',
             '| Measure | Files | Share of 3,000 |', '| --- | ---: | ---: |',
             f'| Catallax, valid after OCCT builder construction | {r["valid_after_builder"]} | {r["valid_share_percent"]}% |',
             f'| OpenCASCADE, healed-valid (published reference) | {r["occt_healed_valid"]} | {r["occt_healed_valid_share_percent"]}% |',
             f'| Catallax, previous capture | {r["previous_valid_after_builder"]} | {100 * r["previous_valid_after_builder"] / 3000:.2f}% |',
             '', f'Gained since the previous capture: {len(r["gained_since_previous"])}. '
             f'Lost, measured and no longer valid: {len(r["lost_since_previous"])} {r["lost_since_previous"]}. '
             f'Previously valid but without a verdict this run (see Not covered): '
             f'{len(r["previous_not_measured"])} {r["previous_not_measured"]}.', '',
             '| Admission | Files |', '| --- | ---: |']
    lines += [f'| {k} | {v} |' for k, v in sorted(r['admission_outcomes'].items())]
    lines += ['', '| Import category | Files |', '| --- | ---: |']
    lines += [f'| {k} | {v} |' for k, v in sorted(r['import_categories'].items())]
    lines += ['', '| First refusal cause | Files |', '| --- | ---: |']
    lines += [f'| `{k}` | {v} |' for k, v in sorted(r['refusal_causes'].items(), key=lambda kv: -kv[1])]
    lines += ['', '| Export outcome of full admissions | Files |', '| --- | ---: |']
    lines += [f'| {k} | {v} |' for k, v in sorted(r['export_outcomes'].items())]
    lines += ['', '## Not covered', '',
              'Runs with no verdict: aborted for memory or time, or ended by a process error. They are '
              'not passes, not losses and not refusals of geometry. The memory cap is sampled per '
              f'observer ({r["max_rss_mb_per_observer"]} MiB), not a host-wide bound.', '']
    lines += [f'- {e["model_id"]}: {e["stage"]} {e["kind"]}' + (f' ({e["detail"]})' if e['detail'] else '')
              for e in r['not_covered']] or ['- None.']
    if r['export_refusal_causes']:
        lines += ['', 'Export refusals: ' + ', '.join(f'`{k}` {v}' for k, v in r['export_refusal_causes'].items()) + '.']
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('scoreboard', 'census'):
        command = commands.add_parser(name)
        command.add_argument('--tools', required=True, help='comparison install directory')
        command.add_argument('--sources', required=True, help='private census source root')
        command.add_argument('--selection', required=True, help='private [model_id, member] list')
        command.add_argument('--out', required=True, help='fresh directory outside the checkout')
    scoreboard_command = commands.choices['scoreboard']
    scoreboard_command.add_argument('--core', type=int, default=min(os.sched_getaffinity(0)))
    scoreboard_command.add_argument('--repeats', type=int, default=1)
    scoreboard_command.add_argument('--max-work', type=int, default=100_000_000)
    scoreboard_command.add_argument('--timeout', type=float, default=120)
    scoreboard_command.add_argument('--machine-class', default='work PC, shared')
    scoreboard_command.add_argument('--part', action='append',
                                    help='smoke runs only: exact case name such as P2-collar.json or census/00000131.json')
    census_command = commands.choices['census']
    # Two workers and a per-observer cap keep one huge source (00005642 peaks at
    # 4.6 GB) from exhausting the shared 20 GB work PC.
    census_command.add_argument('--import-workers', type=int, default=2)
    census_command.add_argument('--max-rss-mb', type=int, default=6144,
                                help='kill one observer above this resident memory and record it as out of memory')
    census_command.add_argument('--export-workers', type=int, default=2)
    args = parser.parse_args()
    {'scoreboard': scoreboard, 'census': census}[args.command](args)


if __name__ == '__main__':
    main()
