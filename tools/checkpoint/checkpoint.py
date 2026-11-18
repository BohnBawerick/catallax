"""Checkpoint runs in one command each: the edit-chain scoreboard and the census rerun.

    python3 -B tools/checkpoint/checkpoint.py scoreboard --tools T --sources S --selection L --out O
    python3 -B tools/checkpoint/checkpoint.py census --tools T --sources S --selection L --out O

`--tools` is the comparison install from head-to-head/prepare.py; it is prepared on
first use. `--sources` and `--selection` are the private census corpus and its
[model_id, member] list. Every output goes to `--out`, outside the checkout.
Both commands build their own binaries from the clean current checkout.

The scoreboard runs P0 to P5, the phase-2 targets in `phase2/` (P8, P9 and the
curved-housing chain) and ten real census parts through Catallax,
OpenCASCADE (CadQuery) and SMLib (smlibpy) with the existing head-to-head
driver, then reports one definition of done for all three kernels. The census
uses the existing default per-body tier observer and selected-document export,
then reports certified, numerical and sampled-native admissions separately
beside occurrence-preserving export validity and OCCT's published reference.
"""
import argparse
from collections import Counter
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
PHASE2 = ROOT / 'tools/checkpoint/phase2'
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


def group(part):
    return 'census' if part.startswith('census/') else 'phase2' if part.startswith('phase2/') else 'phase1'


def scoreboard_tables(report):
    measured = [pair for pair in report['pairs'] if not pair['warmup']]
    order = {'phase1': 0, 'phase2': 1, 'census': 2}
    parts = sorted({pair['part'] for pair in measured}, key=lambda p: (order[group(p)], p))
    rows = []
    for part in parts:
        row = {'part': part, 'engines': {}}
        for engine in report['engines']:
            samples = [pair['engines'][engine] for pair in measured if pair['part'] == part]
            verdicts = [outcome(sample, part.startswith('census/')) for sample in samples]
            done = bool(verdicts) and all(v[0] for v in verdicts)
            first = next((v for v in verdicts if not v[0]), verdicts[0] if verdicts else (False, 'not run', ''))
            entry = {'done': done, 'summary': 'done' if done else first[1], 'detail': first[2]}
            # The curved-housing chain also reads its STEP with healing off; reported, not gating.
            raw = [variant.get('correctness', {}).get('step', {}).get('healing_off')
                   for sample in samples for variant in sample.get('variants', {}).values()]
            if any(r is not None for r in raw):
                entry['healing_off_valid'] = all(r is not None and r['valid'] for r in raw)
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
             'in OCCT 8.0.1.1 as valid closed solids with the kernel\'s own face count, and no '
             'independent check found a wrong value. On a census part the drill and the radius '
             'change must each lower the OCCT-measured volume.', '',
             'Groups: P0 to P5 (13), phase 2 (P8, P9 and the curved-housing chain), census parts (10).', '',
             '| Target | ' + ' | '.join(KERNEL_NAMES[e] for e in engines) + ' | Catallax work (build, edit) | Catallax names |',
             '| --- |' + ' --- |' * len(engines) + ' ---: | --- |']
    for row in rows:
        work = row['engines'].get('catallax', {}).get('work') or {}
        def fmt(value):
            return 'n/a' if value is None else f'{value:,}'
        lines.append('| ' + ' | '.join([row['part'].removesuffix('.json'),
            *[row['engines'][e]['summary'] for e in engines],
            f'{fmt(work.get("construction"))}, {fmt(work.get("edit"))}',
            str(row['engines'].get('catallax', {}).get('names') or 'not-run')]) + ' |')
    sizes = {g: sum(group(r['part']) == g for r in rows) for g in ('phase1', 'phase2', 'census')}
    lines += ['', f'| Kernel | Done, P0 to P5 ({sizes["phase1"]}) | Done, phase 2 ({sizes["phase2"]}) '
              f'| Done, census parts ({sizes["census"]}) |', '| --- | ---: | ---: | ---: |']
    for engine in engines:
        done = {g: sum(r['engines'][engine]['done'] for r in rows if group(r['part']) == g) for g in sizes}
        lines.append(f'| {KERNEL_NAMES[engine]} | {done["phase1"]} | {done["phase2"]} | {done["census"]} |')
    raw = [(r['part'], e, r['engines'][e]['healing_off_valid']) for r in rows for e in engines
           if 'healing_off_valid' in r['engines'][e]]
    if raw:
        lines += ['', 'Healing-off read (raw StepToTopoDS at 1e-7 mm, no ShapeFix), reported beside '
                  'the verdict: ' + '; '.join(f'{part.removesuffix(".json")} / {KERNEL_NAMES[e]}: '
                                              f'{"valid" if ok else "not valid"}' for part, e, ok in raw) + '.']
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
        engines=ENGINES, census_cases=str(cases),
        extra_cases=[(str(PHASE2 / 'cases'), str(PHASE2.parent), 'phase2')]))
    rows = scoreboard_tables(report)
    seconds = time.monotonic() - started
    save(out / 'scoreboard.json', {'format': 'catallax-checkpoint-scoreboard/1',
                                   'source_revision': report['source_revision'],
                                   'wall_seconds': seconds, 'rows': rows})
    (out / 'scoreboard.md').write_text(scoreboard_markdown(report, rows, seconds))
    print((out / 'scoreboard.md').read_text())


def census(args):
    try:
        import OCP  # noqa: F401  The tier collector needs OCP in this interpreter.
    except ImportError:
        checker = json.loads(runtime(args.tools).read_text())['checker']
        os.execv(checker, [checker, '-B', str(Path(__file__).resolve()), *sys.argv[1:]])
    assert not git('status', '--porcelain'), 'the census binds a clean checkout'
    out = fresh(args.out)
    collector = RESEARCH / 'opened_import_census.py'
    record = out / 'opened-build.json'
    started = time.monotonic()
    subprocess.run([sys.executable, '-B', str(collector), 'build',
                    '--record', str(record)], check=True)
    build_seconds = time.monotonic() - started
    started = time.monotonic()
    subprocess.run([sys.executable, '-B', str(collector), 'collect',
                    '--build-record', str(record), '--sources', str(args.sources),
                    '--selection', str(args.selection), '--source-hashes', str(HASHES),
                    '--out', str(out / 'opened'), '--workers', '2',
                    '--max-rss-mb', str(args.max_rss_mb)], check=True)
    collect_seconds = time.monotonic() - started
    report = json.loads((out / 'opened/results.json').read_text())
    proofs = Counter()
    full_files = Counter()
    valid_files = Counter()
    with (out / 'opened/results.csv').open(newline='') as stream:
        csv.field_size_limit(64 * 1024 * 1024)
        for row in csv.DictReader(stream):
            inventory = json.loads(row.get('proofs') or '{}')
            proofs.update(inventory)
            if row['outcome'] != 'full':
                continue
            # File groups are disjoint, labelled by their least-checked body.
            category = next((name for name in ('native_sampled', 'numerical', 'certified')
                             if inventory.get(name)), 'unknown')
            full_files[category] += 1
            if row.get('all_solids_export_valid') == 'True':
                valid_files[category] += 1
    report.update(
        format='catallax-checkpoint-census/2',
        body_proofs=dict(proofs), full_file_proofs=dict(full_files),
        valid_file_proofs=dict(valid_files),
        wall_seconds={'build': build_seconds, 'collect': collect_seconds},
        occt_healed_valid=OCCT_HEALED_VALID,
        checkpoint_1_legacy_valid=374)
    save(out / 'census.json', report)
    (out / 'census.md').write_text(census_markdown(report))
    print((out / 'census.md').read_text())


def census_markdown(r):
    total = sum(r['wall_seconds'].values())
    lines = ['# Census rerun', '',
             f'Source `{r["source_revision"]}`. Complete: {r["complete"]}. '
             f'Evaluated: {r["evaluated"]}/{r["denominator"]}. '
             f'Wall time {total / 60:.1f} min, including build.', '',
             '| Measure | Files | Share of 3,000 |', '| --- | ---: | ---: |']
    for label, count in [
        ('Catallax, complete-source occurrence-preserving valid solid export',
         r['all_solids_export_valid_files']),
        ('Catallax, at least one observed valid exported solid',
         r['files_with_any_valid_exported_solid']),
        ('Catallax, all isolated prototype exports valid',
         r['prototype_solids_export_valid_files']),
        ('OpenCASCADE, healed-valid, published reference', r['occt_healed_valid']),
        ('Checkpoint 1 legacy import/export, historical', r['checkpoint_1_legacy_valid']),
    ]:
        lines.append(f'| {label} | {count} | {100 * count / DENOMINATOR:.2f}% |')
    lines += ['', 'File groups below use the least-checked body in each full admission. '
              'Numerical and sampled-native bodies are approximate, not certified. '
              'Valid export requires the complete occurrence-preserving document, '
              'not just isolated prototypes.', '',
              '| Proof group | Full solid admissions | Complete-source valid exports | Bodies, including partial files |',
              '| --- | ---: | ---: | ---: |']
    for proof in ('certified', 'numerical', 'native_sampled', 'unknown'):
        lines.append(f'| {proof} | {r["full_file_proofs"].get(proof, 0)} | '
                     f'{r["valid_file_proofs"].get(proof, 0)} | '
                     f'{r["body_proofs"].get(proof, 0)} |')
    lines += ['', '| Opening outcome | Files |', '| --- | ---: |']
    lines += [f'| {k} | {v} |' for k, v in sorted(r['outcomes'].items())]
    lines += ['', f'Retained source-surface bodies: {r["body_proofs"].get("retained_source_surfaces", 0)}.',
              '', '| Export outcome | Files |', '| --- | ---: |']
    lines += [f'| {k} | {v} |' for k, v in sorted(r['export_outcomes'].items())]
    lines += ['', '| First body blocker | Refusal | Resource-limited | Bodies | Files, overlapping |',
              '| --- | --- | --- | ---: | ---: |']
    lines += [f'| {b["stage"]} | `{b["refusal"]}` | {b["resource_limited"]} | '
              f'{b["bodies"]} | {b["files"]} |' for b in r['body_first_blockers']]
    lines += ['', '## Observation limits', '',
              f'Any-valid verdict unknown: {len(r["any_valid_export_unmeasured_files"])} files. '
              f'Valid-solid counts are lower bounds in {len(r["any_valid_export_lower_bound_files"])} files. '
              'Unobserved siblings are not invalid solids. Time, memory and work limits '
              'are not geometric refusals.', '',
              'Unknown any-valid file IDs: ' + ', '.join(r['any_valid_export_unmeasured_files']) + '.',
              '', f'No common-tolerance source-fidelity or unhealed-validity claim. '
              f'Observer memory cap: {r["max_rss_mb_per_observer"]} MiB. '
              'The raw collector receipts and per-file rows are in `opened/`.']
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
                                    help='smoke runs only: exact case name such as P2-collar.json, '
                                         'phase2/P9-lofted-duct.json or census/00000131.json')
    census_command = commands.choices['census']
    # Two workers and a per-observer cap keep one huge source (00005642 peaks at
    # 4.6 GB) from exhausting the shared 20 GB work PC.
    census_command.add_argument('--max-rss-mb', type=int, default=6144,
                                help='kill one observer above this resident memory and record it as out of memory')
    args = parser.parse_args()
    {'scoreboard': scoreboard, 'census': census}[args.command](args)


if __name__ == '__main__':
    main()
