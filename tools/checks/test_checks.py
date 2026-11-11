"""R1 acceptance cases G1-G6, S1-S3 and C3-C4 from the independent designer."""
import copy
import contextlib
import io
import json
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import metrics
from pathlib import Path
import re
import subprocess
import unittest

import checks


class TimingTests(unittest.TestCase):
    def test_jobs_502_records_unavailable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stub = root / 'gh'
            stub.write_text("#!/usr/bin/env python3\n"
                            "import base64, json, sys\n"
                            f"with open({str(root / 'calls')!r}, 'a') as log: log.write(sys.argv[2] + '\\n')\n"
                            "if '/jobs?' in sys.argv[2]:\n"
                            "    print('gh: Bad Gateway (HTTP 502)', file=sys.stderr)\n"
                            "    sys.exit(1)\n"
                            "print(base64.b64encode(json.dumps({'run_attempt': 1}).encode()).decode())\n")
            stub.chmod(0o755)
            output = root / 'timings.json'
            console = io.StringIO()
            with patch.object(metrics, 'GITHUB_CLI', str(stub)), \
                    patch.object(metrics, 'sleep') as sleep, \
                    patch('sys.argv', ['metrics.py', '--run-id', '123', '--label', 'test',
                                       '--selection', str(root / 'selection.json'),
                                       '--out', str(output)]), contextlib.redirect_stdout(console):
                metrics.main()
            self.assertEqual(json.loads(output.read_text())['status'], 'timings unavailable')
            self.assertIn('::warning::', console.getvalue())
            self.assertIn('HTTP 502', console.getvalue())
            calls = (root / 'calls').read_text().splitlines()
            self.assertEqual([url.split('?')[1] for url in calls if '/jobs?' in url],
                             ['per_page=100&page=1'] * 2 + ['per_page=30&page=1'] * 2)
            self.assertEqual([call.args[0] for call in sleep.call_args_list], [20, 40, 60])


class GateTests(unittest.TestCase):
    def setUp(self):
        self.plan = {'identity': {'commit': 'candidate', 'run_id': '123', 'run_attempt': '1',
                                 'corpus_sha256': 'corpus', 'metadata_sha256': 'metadata',
                                 'driver_sha256': 'driver', 'raw_policy_sha256': 'raw-policy',
                                 'corrections_sha256': 'corrections'},
                     'correction_base_answers': [{'index': 262}],
                     'selection': {'jobs': {'fmt': 1, 'clippy': 1, 'test-build': 1, 'tests': 4,
                                            'build': 1, 'step-regressions': 1, 'replay': 4},
                                   'packages': ['catallax-cli'],
                                   'replay_partitions': checks.replay_catalog(checks.read(checks.ROOT / checks.INVENTORY)),
                                   'selected_cases': 477},
                     'declared': ['fmt', 'clippy', 'test-build', 'tests', 'build',
                                  'step-regressions', 'replay', 'scoreboard']}
        self.reported = {name: {'result': 'success'} for name in self.plan['selection']['jobs']}
        self.reported['scoreboard'] = {'result': 'skipped'}
        self.records = []
        for name, partitions in self.plan['selection']['jobs'].items():
            for part in range(1, partitions + 1):
                record = dict(self.plan['identity'], job=name, partition=part,
                              binary_sha256='binary' if name in {'build', 'step-regressions', 'replay'} else None,
                              result='success')
                if name == 'replay':
                    selected = checks.replay_partition(self.plan, part)
                    record.update(workers=4, fresh=True, cache_hits=0,
                                  family=selected['family'], indices=selected['indices'],
                                  selected_cases=len(selected['indices']),
                                  verified_corrections=[262] if 262 in selected['indices'] else [])
                self.records.append(record)

    def errors(self, reported=None, records=None, plan=None):
        return checks.decide(plan or self.plan, reported or self.reported,
                             records if records is not None else self.records, 'binary')

    def test_g1_g2_g3_g5_g6(self):
        self.assertEqual(self.errors(), [])
        for result in ('cancelled', 'skipped'):
            with self.subTest(result=result):
                reported = copy.deepcopy(self.reported)
                reported['clippy']['result'] = result
                self.assertIn(f'clippy: {result}', self.errors(reported=reported))
        reported = copy.deepcopy(self.reported)
        reported['extra'] = {'result': 'success'}
        self.assertIn('extra: unexpected job', self.errors(reported=reported))
        for unexpected in ('extra', 'tests (5)'):
            self.assertIn(f'{unexpected}: unexpected job', checks.decide(
                self.plan, self.reported, self.records, 'binary', actual_jobs=[unexpected]))
        records = [r for r in self.records if (r['job'], r['partition']) != ('tests', 4)]
        self.assertIn('tests: partition 4/4 missing', self.errors(records=records))
        # G6's legitimately unselected scoreboard has no evidence and passes G1.

    def test_docs_only_skipped_replay_placeholder(self):
        plan = copy.deepcopy(self.plan)
        plan['selection'] = checks.select(
            ['docs/PLAN.md'], checks.read(checks.ROOT / checks.RULES), {}, {})
        plan['declared'] += ['release-evidence', 'release-dry-run', 'timings']
        reported = {name: {'result': 'success' if name == 'fmt' else 'skipped'}
                    for name in plan['declared']}
        records = [record for record in self.records if record['job'] == 'fmt']
        actual = ['select', 'release-evidence', 'release-dry-run', 'fmt', 'clippy',
                  'tests', 'test-build', 'build', 'scoreboard',
                  'replay (${{ matrix.partition }})', 'step-regressions', 'check', 'timings']
        self.assertEqual(checks.decide(plan, reported, records, actual_jobs=actual), [])
        for unexpected in ('extra', 'replay (5)', 'replay (${{ matrix.other }})'):
            with self.subTest(job=unexpected):
                self.assertIn(f'{unexpected}: unexpected job', checks.decide(
                    plan, reported, records, actual_jobs=actual + [unexpected]))
        for result in (None, 'success', 'failure', 'cancelled'):
            with self.subTest(result=result):
                reported['replay'] = {'result': result}
                self.assertIn('replay (${{ matrix.partition }}): unexpected job', checks.decide(
                    plan, reported, records, actual_jobs=actual))

    def test_selected_replay_rejects_unexpanded_placeholder(self):
        actual = 'replay (${{ matrix.partition }})'
        self.assertIn(f'{actual}: unexpected job', checks.decide(
            self.plan, self.reported, self.records, 'binary', actual_jobs=[actual]))

    def test_g4_c4_identity_mismatches(self):
        for field in ('commit', 'run_id', 'binary_sha256', 'corpus_sha256',
                      'metadata_sha256', 'driver_sha256', 'raw_policy_sha256',
                      'corrections_sha256', 'cache_hits', 'workers', 'fresh', 'verified_corrections',
                      'family', 'indices', 'selected_cases'):
            with self.subTest(field=field):
                records = copy.deepcopy(self.records)
                replay = next(r for r in records if r['job'] == 'replay')
                replay[field] = {'cache_hits': 1, 'workers': 1, 'fresh': False}.get(field, 'stale')
                self.assertTrue(any(f'replay: stale {field}:' in e for e in self.errors(records=records)))

    def test_missing_job_ignores_candidate_needs(self):
        reported = copy.deepcopy(self.reported)
        del reported['clippy']
        self.assertIn('clippy: missing', self.errors(reported=reported))

    def test_every_family_is_required_and_cannot_borrow_other_family_evidence(self):
        for part in range(1, 5):
            with self.subTest(partition=part):
                records = [r for r in self.records if (r['job'], r['partition']) != ('replay', part)]
                self.assertIn(f'replay: partition {part}/4 missing', self.errors(records=records))
                records = copy.deepcopy(self.records)
                replay = [r for r in records if r['job'] == 'replay']
                replay[part - 1]['indices'] = replay[part % 4]['indices']
                self.assertTrue(any('stale indices' in e for e in self.errors(records=records)))

    def test_empty_or_malformed_plan_fails_closed(self):
        for fault in ('jobs', 'packages', 'replay_partitions', 'selected_cases'):
            with self.subTest(field=fault):
                plan = copy.deepcopy(self.plan)
                plan['selection'][fault] = {} if fault == 'jobs' else [] if fault != 'selected_cases' else 0
                self.assertTrue(self.errors(plan=plan))
        for result in ('skipped', 'failure', 'cancelled'):
            reported = copy.deepcopy(self.reported)
            reported['replay']['result'] = result
            self.assertIn(f'replay: {result}', self.errors(reported=reported))

    def test_archive_and_step_jobs_cannot_be_skipped_or_lose_evidence(self):
        for name in ('test-build', 'step-regressions'):
            with self.subTest(job=name):
                for result in ('failure', 'cancelled', 'skipped'):
                    reported = copy.deepcopy(self.reported)
                    reported[name]['result'] = result
                    self.assertIn(f'{name}: {result}', self.errors(reported=reported))
                records = [record for record in self.records if record['job'] != name]
                self.assertIn(f'{name}: partition 1/1 missing', self.errors(records=records))
        records = copy.deepcopy(self.records)
        next(record for record in records if record['job'] == 'step-regressions')['binary_sha256'] = 'stale'
        self.assertTrue(any('step-regressions: stale binary_sha256:' in error
                            for error in self.errors(records=records)))

    def test_c3_metadata_cannot_hide_changed_answer(self):
        base = [{'index': 7, 'accepted_answer_sha256':
                 '0203b8b177ca59e519192e31b591c51d3668e031ed6a83e0cae0539f84a8bdcd'}]
        changed = [{'index': 7, 'accepted_answer_sha256': 'changed'}]
        self.assertEqual(checks.metadata_errors(base, base), [])
        # C3(a) and C3(b) cannot grant approval without an admissible correction
        # backed by the base revision's trusted raw-invalid ruling.
        self.assertIn('without base approval', checks.metadata_errors(base, changed)[0])


class SelectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules = checks.read(checks.ROOT / checks.RULES)
        cls.graph, cls.build = checks.cargo_graph(checks.ROOT)

    def route(self, paths, rules=None, event='pull_request', ref='', changed_jobs=()):
        return checks.select(paths, rules or self.rules, self.graph, self.build, event, ref, changed_jobs)

    def test_s1_routing(self):
        def packages(*names):
            return {'catallax-' + name for name in names}
        step = packages('app', 'cli', 'compare', 'features', 'field', 'iges', 'import', 'measure', 'step')
        rows = [
            (['crates/catallax-numbers/src/lib.rs'], set(self.graph), True),
            (['crates/catallax-step/src/export.rs'], step, True),
            (['crates/catallax-iges/src/lib.rs'], packages('iges', 'cli'), True),
            (['crates/catallax-cli/src/verification.rs'], packages('cli'), True),
            (['crates/catallax-import/src/lib.rs'], packages('import'), None),
            (['parts/P2-collar.json'], packages('app', 'cli', 'compare', 'features', 'field', 'iges'), None),
            (['.claude/skills/verify-catallax/references/examples/certified-recipe.json'],
             packages('app', 'cli', 'measure'), None),
            (['.claude/skills/verify-catallax/scripts/control-catallax.py'], set(), True),
            (['Cargo.lock'], set(self.graph), True),
            (['docs/PLAN.md'], set(), False),
            (['newdir/unmapped.txt'], set(self.graph), True),
            (['crates/catallax-iges/src/format.rs', 'crates/catallax-step/src/format.rs'], step, True),
            (['.github/workflows/check.yml'], set(), False),
            (['tools/checks/rules.json'], set(), False),
            (['tools/checks/earlier-answers.json'], set(self.graph), True),
        ]
        for paths, minimum, replay in rows:
            with self.subTest(paths=paths):
                selected = self.route(paths)
                self.assertTrue(minimum <= set(selected['packages']), selected)
                if replay is not None:
                    self.assertEqual(selected['selected_cases'], 477 if replay else 0)
                if paths == ['docs/PLAN.md']:
                    self.assertEqual(selected['packages'], [])
                if paths == ['parts/P2-collar.json']:
                    self.assertIn('scoreboard', selected['jobs'])
                if paths == ['.claude/skills/verify-catallax/scripts/control-catallax.py']:
                    self.assertIn('build', selected['jobs'])
                self.assertEqual('test-build' in selected['jobs'], 'tests' in selected['jobs'])
                self.assertEqual('step-regressions' in selected['jobs'], 'build' in selected['jobs'])

    def test_s2_base_and_candidate_union(self):
        changed = copy.deepcopy(self.rules)
        changed['replay_paths'] = []
        paths = ['crates/catallax-boolean/src/lib.rs']
        base, candidate = self.route(paths), self.route(paths, changed)
        union = checks.union_selection(base, candidate)
        self.assertEqual((union['base_selected'], union['candidate_selected'], union['selected_cases']),
                         (477, 0, 477))
        changed['jobs']['extra'] = 1
        union = checks.union_selection(base, self.route(paths, changed))
        plan = {'identity': {}, 'declared': list(changed['jobs']), 'selection': union}
        plan['selection']['replay_partitions'] = checks.replay_catalog(checks.read(checks.ROOT / checks.INVENTORY))
        errors = checks.decide(plan, {}, [])
        self.assertIn('extra: missing', errors)

    def test_base_only_selection_still_requires_split_checks(self):
        old_rules = copy.deepcopy(self.rules)
        del old_rules['jobs']['test-build']
        del old_rules['jobs']['step-regressions']
        base = self.route(['crates/catallax-cli/src/verification.rs'], old_rules)
        candidate = self.route(['docs/PLAN.md'])
        jobs = checks.union_selection(base, candidate)['jobs']
        self.assertEqual(jobs['test-build'], 1)
        self.assertEqual(jobs['step-regressions'], 1)

    def test_ci_only_changes_run_tooling_and_changed_jobs_without_shards_or_replay(self):
        before = {'jobs': {
            'fmt': {'steps': [{'run': 'ci/run fmt'}]},
            'clippy': {'steps': [{'run': 'ci/run clippy --workspace'}]},
            'build': {'steps': [{'run': 'cargo build --release'}]},
            'tests': {'needs': 'test-build', 'steps': [{'run': 'old tests'}]},
            'test-build': {'steps': [{'run': 'archive tests'}]},
            'step-regressions': {'needs': ['select', 'build'], 'steps': [{'run': 'old STEP check'}]},
            'replay': {'needs': 'build', 'steps': [{'run': 'old replay'}]}}}
        after = copy.deepcopy(before)
        after['jobs']['clippy']['steps'] = [{'run': 'ci/run clippy'}]
        changed = checks.changed_check_jobs(before, after)
        paths = ['.github/workflows/check.yml', 'ci/run', 'tools/checks/rules.json', 'docs/RELEASING.md']
        selected = self.route(paths, changed_jobs=changed)
        self.assertEqual(selected['jobs'], {'fmt': 1, 'clippy': 1})
        self.assertEqual(selected['selected_cases'], 0)
        self.assertEqual(set(selected['packages']), set(self.graph))
        unchanged = self.route(paths, changed_jobs=checks.changed_check_jobs(after, after))
        self.assertEqual(unchanged['jobs'], {'fmt': 1})
        self.assertEqual(unchanged['packages'], [])
        after['jobs']['step-regressions']['steps'] = [{'run': 'changed STEP check'}]
        selected = self.route(paths, changed_jobs=checks.changed_check_jobs(before, after))
        self.assertEqual(set(selected['jobs']), {'fmt', 'clippy', 'build', 'step-regressions'})
        source = self.route(paths + ['crates/catallax-cli/src/verification.rs'], changed_jobs=changed)
        self.assertEqual(source['selected_cases'], 477)
        self.assertEqual(source['jobs']['tests'], 4)
        for path in ('tools/checks/earlier-answers.json', 'tools/checks/replay-corrections.json',
                     'tools/checks/raw_step.py'):
            self.assertEqual(self.route([path])['selected_cases'], 477)

    def test_ci_rule_migration_runs_full_once_then_obeys_changed_job_selection(self):
        old_rules = copy.deepcopy(self.rules)
        del old_rules['ci_paths']
        before = self.route(['tools/checks/rules.json'], old_rules)
        after = self.route(['tools/checks/rules.json'])
        migration = checks.union_selection(before, after)
        self.assertEqual(migration['selected_cases'], 477)
        self.assertEqual(migration['jobs']['tests'], 4)
        future = checks.union_selection(after, after)
        self.assertEqual(future['jobs'], {'fmt': 1})
        self.assertEqual(future['selected_cases'], 0)
        build_only = self.route(['.github/workflows/check.yml'], changed_jobs={'build'})
        future = checks.union_selection(build_only, build_only)
        self.assertEqual(future['jobs'], {'fmt': 1, 'build': 1})

    def test_workflow_environment_changes_effective_jobs_and_their_dependencies(self):
        before = {'env': {'CARGO_BUILD_JOBS': '2'}, 'jobs': {
            'build': {}, 'tests': {}, 'replay': {}, 'step-regressions': {'needs': 'build'}}}
        after = copy.deepcopy(before)
        after['env']['CARGO_BUILD_JOBS'] = '4'
        self.assertEqual(checks.changed_check_jobs(before, after), {'build', 'step-regressions', 'tests', 'replay'})

    def test_s3_full_events(self):
        for event, ref in [('push', 'refs/heads/release/0.2.0'), ('workflow_dispatch', 'refs/heads/main'),
                           ('schedule', 'refs/heads/main')]:
            with self.subTest(event=event, ref=ref):
                selected = self.route(['docs/PLAN.md'], event=event, ref=ref)
                self.assertEqual(set(selected['packages']), set(self.graph))
                self.assertEqual(selected['jobs'], self.rules['jobs'])
                self.assertEqual(selected['selected_cases'], 477)

    def test_changed_test_or_replay_jobs_cannot_take_ci_only_shortcut(self):
        for name in ('tests', 'replay'):
            with self.subTest(job=name):
                before = {'jobs': {name: {'steps': [{'run': 'before'}]}}}
                after = {'jobs': {name: {'steps': [{'run': 'after'}]}}}
                changed = checks.changed_check_jobs(before, after)
                self.assertIn(name, changed)
                selected = self.route(['.github/workflows/check.yml'], changed_jobs=changed)
                self.assertEqual(selected['jobs'], self.rules['jobs'])
                self.assertEqual(selected['selected_cases'], 477)

    def test_mixed_fixtures_keep_consumers_and_kernel_coverage(self):
        fixtures = ['parts/P2-collar.json',
                    '.claude/skills/verify-catallax/references/examples/certified-recipe.json']
        sources = ['crates/catallax-import/src/lib.rs', 'crates/catallax-cli/tests/exchange.rs',
                   'crates/catallax-boolean/src/lib.rs']
        for fixture in fixtures:
            for source in sources:
                with self.subTest(fixture=fixture, source=source):
                    baseline = self.route([fixture])
                    mixed = self.route([fixture, source])
                    self.assertTrue(set(baseline['packages']) <= set(mixed['packages']))
                    self.assertTrue(set(baseline['jobs']) <= set(mixed['jobs']))
                    if self.route([source])['selected_cases']:
                        self.assertEqual(mixed['selected_cases'], 477)


class ReplayFamilyTests(unittest.TestCase):
    def test_public_catalog_and_unknown_owner_fallback(self):
        metadata = checks.read(checks.ROOT / checks.INVENTORY)
        catalog = checks.replay_catalog(metadata)
        self.assertEqual({r['family']: len(r['indices']) for r in catalog},
                         {'unite': 85, 'subtract': 89, 'intersect': 86, 'mixed': 217})
        self.assertEqual({r['owner'] for r in metadata}, {'unknown'})
        self.assertEqual(sorted(i for r in catalog for i in r['indices']), list(range(len(metadata))))
        selection = {'jobs': {'fmt': 1, 'build': 1, 'replay': len(catalog)},
                     'selected_cases': len(metadata), 'replay_partitions': catalog}
        self.assertEqual(checks.selection_errors(selection, metadata), [])
        narrowed = copy.deepcopy(selection)
        narrowed['replay_partitions'] = catalog[:1]
        narrowed['jobs']['replay'] = 1
        narrowed['selected_cases'] = len(catalog[0]['indices'])
        self.assertIn('unknown replay owners require the full public catalog',
                      checks.selection_errors(narrowed, metadata))
        for rows in ([], [metadata[1]], metadata + [metadata[0]]):
            with self.assertRaises(ValueError):
                checks.replay_catalog(rows)

    def test_raw_step_positions_are_original_corpus_indices(self):
        cases = [{'path': 'first'}, {'path': 'second'}]
        self.assertEqual(checks.indexed_replay_cases({'cases': cases, 'replay_indices': [262, 431]}),
                         [(262, cases[0]), (431, cases[1])])
        for indices in ([], [262], [262, 262], [-1, 431], [True, 431]):
            with self.subTest(indices=indices), self.assertRaises(ValueError):
                checks.indexed_replay_cases({'cases': cases, 'replay_indices': indices})
        with self.assertRaises(ValueError):
            checks.indexed_replay_cases({'cases': []})

    def test_family_timings_require_each_complete_unique_summary(self):
        catalog = checks.replay_catalog(checks.read(checks.ROOT / checks.INVENTORY))
        selection = {'identity': {'commit': 'candidate'},
                     'selection': {'replay_partitions': catalog}}
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for part in catalog:
                path = Path(directory) / f"{part['family']}.json"
                checks.write(path, dict(part, selected_cases=len(part['indices']),
                                        elapsed_ms=part['partition'] * 100, commit='candidate'))
                paths.append(path)
            self.assertEqual(metrics.family_aggregate(selection, paths)['selected_cases'], 477)
            self.assertEqual(metrics.family_aggregate(selection, paths)['elapsed_ms'], 400)
            for bad in (paths[:-1], paths + paths[:1], [paths[0]] * 4):
                with self.assertRaises(ValueError):
                    metrics.family_aggregate(selection, bad)
            row = checks.read(paths[0])
            row['indices'] = row['indices'][:-1]
            checks.write(paths[0], row)
            with self.assertRaises(ValueError):
                metrics.family_aggregate(selection, paths)

    def test_replay_executes_only_family_and_retains_failures(self):
        # Tiny synthetic requests exercise execution without private inputs or Rust.
        driver = checks.load_module('family_test_driver', checks.ROOT / checks.DRIVER)
        for fault in (None, 'answer', 'cache', 'work', 'internal-error', 'missing', 'worker'):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                inputs = [{'source': 'question', 'problem': {'format': 'test', 'position': i}}
                          for i in range(6)]
                metadata = [dict(index=i, family='first' if i % 2 else 'second', owner='unknown',
                                 request_sha256=driver.digest(dict(max_work=100_000_000, **row,
                                                                   format=driver.CONTRACT)),
                                 accepted_answer_sha256='accepted', accepted_manifest_sha256='accepted',
                                 accepted_files={}) for i, row in enumerate(inputs)]
                checks.write(root / checks.INVENTORY, metadata)
                checks.write(root / checks.CORRECTIONS, {'corrections': []})
                checks.write(root / checks.UPGRADES, {'format': checks.UPGRADE_FORMAT, 'upgrades': []})
                checks.write(root / 'inputs.json', inputs)
                (root / 'binary').write_text('test binary')
                checks.write(root / checks.RULES, {'corpus_sha256': checks.sha(root / 'inputs.json')})
                plan = {'identity': {'commit': 'candidate',
                                     'metadata_sha256': checks.sha(root / checks.INVENTORY),
                                     'corrections_sha256': checks.sha(root / checks.CORRECTIONS),
                                     'upgrades_sha256': checks.sha(root / checks.UPGRADES)},
                        'correction_base_answers': [],
                        'selection': {'jobs': {'fmt': 1, 'build': 1, 'replay': 2}, 'selected_cases': 6,
                                      'replay_partitions': checks.replay_catalog(metadata)}}
                checks.write(root / 'plan.json', plan)
                args = SimpleNamespace(plan=root / 'plan.json', partition=1, inputs=root / 'inputs.json',
                                       binary=root / 'binary', out=root / 'out', summary=root / 'summary.json')
                seen = []

                def run(_driver, binary, requests, out, workers):
                    cases = []
                    for line in Path(requests).read_text().splitlines():
                        request = json.loads(line)
                        seen.append(request['problem']['position'])
                        target = root / f"case-{len(cases)}"
                        checks.write(target / 'request.json', dict(request, out=str(target)))
                        manifest = {'status': 'answered'}
                        if not cases and fault == 'internal-error':
                            manifest['status'] = 'internal-error'
                        cases.append({'path': str(target), 'cache_hit': not cases and fault == 'cache',
                                      'response': {'manifest': manifest,
                                                   'work': {'used': 100_000_001 if fault == 'work' else 1}}})
                    if fault == 'missing':
                        cases.pop()
                    return {'cases': cases, 'config': {'binary_sha256': checks.sha(binary)},
                            'execution': {'workers': workers, 'fresh': True,
                                          'worker_exits': {str(i): int(fault == 'worker') for i in range(4)}}}

                release = SimpleNamespace(OUTPUTS=[], answer_digest=lambda *_: 'changed' if fault == 'answer' else 'accepted')
                with patch.object(checks, 'ROOT', root), patch.object(checks, 'git', return_value=b'candidate\n'), \
                        patch.object(checks.os, 'sched_getaffinity', return_value=set(range(4))), \
                        patch.object(checks, 'run_fresh', side_effect=run), \
                        patch.object(checks, 'load_module', side_effect=[driver, release]), \
                        contextlib.redirect_stdout(io.StringIO()):
                    if fault:
                        with self.assertRaises(ValueError):
                            checks.replay(args)
                    else:
                        checks.replay(args)
                self.assertEqual(seen, [1, 3, 5])
                if fault not in {'missing', 'worker'}:
                    summary = checks.read(args.summary)
                    self.assertEqual(summary['indices'], [1, 3, 5])
                    self.assertEqual(summary['selected_cases'], 3)
                    expected_work = 100_000_001 if fault == 'work' else 1
                    self.assertEqual(summary['work_cases'],
                                     [{'index': i, 'used': expected_work} for i in [1, 3, 5]])
                    self.assertEqual(summary['work_total'],
                                     sum(row['used'] for row in summary['work_cases']))
                    self.assertEqual(bool(summary['failures']), bool(fault))
                    self.assertEqual(checks.read(root / 'out/run.json')['replay_indices'], [1, 3, 5])


class RawStepPolicyTests(unittest.TestCase):
    def test_captain_ruling_table(self):
        policy = checks.read(checks.ROOT / checks.RAW_POLICY)
        policy['known_raw_invalid'] = [262, 431]
        observed = [
            {'index': 262, 'step_sha256': '770dfffc0237a182539afff13d14b0307dbbf928be5a089bc43a665ee06319e1', 'valid': False},
            {'index': 431, 'step_sha256': '4229c04a36d92d059398573423059238c34732ba8d220933094c260c7ef4a9a3', 'valid': False}]
        rows = [('unchanged approved invalid bytes', observed, policy, None)]
        extra_invalid = copy.deepcopy(observed)
        extra_invalid.append({'index': 7, 'step_sha256': 'other', 'valid': False})
        rows.append(('other invalid output', extra_invalid, policy, 'index 7: unlisted'))
        changed_invalid = copy.deepcopy(observed)
        changed_invalid[0]['step_sha256'] = 'changed'
        rows.append(('listed bytes changed, still invalid', changed_invalid, policy, 'bytes changed'))
        unapproved = copy.deepcopy(policy)
        unapproved['known_raw_invalid'].append(7)
        rows.append(('unapproved list addition', observed, unapproved, 'no captain ruling'))
        forged = copy.deepcopy(unapproved)
        forged['rulings'][0]['approved'].append({'index': 7, 'step_sha256': 'other'})
        rows.append(('self-declared authority expansion', observed, forged, 'trusted captain ruling'))
        valid = copy.deepcopy(observed)
        valid[0].update(step_sha256='upgraded-valid-bytes', valid=True)
        rows.append(('valid upgrade still listed', valid, policy, 'must leave'))
        removed = copy.deepcopy(policy)
        removed['known_raw_invalid'].remove(262)
        rows.append(('valid upgrade removed', valid, removed, None))
        rows.append(('invalid entry removed', observed, removed, 'index 262: unlisted'))
        for name, observations, candidate, expected_error in rows:
            with self.subTest(name=name):
                errors = checks.raw_step_errors(observations, candidate)
                if expected_error is None:
                    self.assertEqual(errors, [])
                else:
                    self.assertTrue(any(expected_error in error for error in errors), errors)


class StepCorrectionTests(unittest.TestCase):
    def test_invalid_to_valid_correction_table(self):
        driver = checks.load_module('correction_test_driver', checks.ROOT / checks.DRIVER)
        policy = checks.read(checks.ROOT / checks.RAW_POLICY)
        policy['known_raw_invalid'] = [262, 431]
        approved = policy['rulings'][0]['approved'][0]['step_sha256']
        before = {
            'index': 262, 'request_sha256': 'a' * 64, 'family': 'mixed', 'owner': 'unknown',
            'accepted_answer_sha256': 'b' * 64,
            'accepted_files': {'answer.step': approved, 'display.obj': None,
                               'print.stl': None, 'answer.igs': None}}
        manifest = {
            'status': 'answered', 'exact': True, 'marks': [],
            'properties': {'volume': {'lo': '1', 'hi': '1'}},
            'step': {'added_edges': 2, 'added_vertices': 1,
                     'max_deviation': '1e-9', 'approximated_curves': 1}}
        before['accepted_manifest_sha256'] = driver.digest(manifest)
        unchanged = copy.deepcopy(manifest)
        for field in checks.STEP_CORRECTION_FIELDS:
            unchanged['step'].pop(field.split('.')[1])
        correction = {
            'index': 262, 'class': 'invalid-to-valid-step', 'request_sha256': 'a' * 64,
            'before_sha256': 'b' * 64, 'after_sha256': 'd' * 64,
            'before_manifest_sha256': before['accepted_manifest_sha256'],
            'after_manifest_sha256': driver.digest(manifest),
            'before_step_sha256': approved, 'after_step_sha256': 'e' * 64,
            'manifest_fields': sorted(checks.STEP_CORRECTION_FIELDS),
            'unchanged_manifest_sha256': driver.digest(unchanged)}
        after = copy.deepcopy(before)
        after.update(accepted_answer_sha256=correction['after_sha256'],
                     accepted_manifest_sha256=correction['after_manifest_sha256'])
        after['accepted_files']['answer.step'] = correction['after_step_sha256']
        cases = [
            ('exact approved correction', None, False),
            ('candidate still raw-invalid', 'invalid', True),
            ('not on base known list', 'unlisted', True),
            ('before bytes differ from captain ruling', 'before-hash', True),
            ('properties changed', 'properties', True),
            ('inventory changed without a correction', 'missing', True),
            ('protected field and self-supplied commitment changed', 'self-commitment', True)]
        for name, fault, fails in cases:
            with self.subTest(case=name):
                base, candidate, row, base_policy, output = copy.deepcopy(
                    (before, after, correction, policy, manifest))
                valid = fault != 'invalid'
                if fault == 'unlisted':
                    base_policy['known_raw_invalid'].remove(262)
                elif fault == 'before-hash':
                    base['accepted_files']['answer.step'] = 'f' * 64
                    row['before_step_sha256'] = 'f' * 64
                elif fault == 'properties':
                    output['properties']['volume']['hi'] = '2'
                    row['after_manifest_sha256'] = driver.digest(output)
                    candidate['accepted_manifest_sha256'] = row['after_manifest_sha256']
                elif fault == 'self-commitment':
                    output['exact'] = False
                    row['after_manifest_sha256'] = driver.digest(output)
                    candidate['accepted_manifest_sha256'] = row['after_manifest_sha256']
                    row['after_sha256'] = candidate['accepted_answer_sha256'] = 'f' * 64
                    forged = copy.deepcopy(output)
                    for field in checks.STEP_CORRECTION_FIELDS:
                        forged['step'].pop(field.split('.')[1], None)
                    row['unchanged_manifest_sha256'] = driver.digest(forged)
                rows = [] if fault == 'missing' else [row]
                errors = checks.metadata_errors(
                    [base], [candidate],
                    {'format': 'catallax-replay-corrections/1', 'corrections': rows}, base_policy)
                errors += checks.correction_manifest_errors(output, row, driver, manifest)
                candidate_policy = copy.deepcopy(base_policy)
                candidate_policy['known_raw_invalid'] = []
                errors += checks.raw_step_errors(
                    [{'index': 262, 'step_sha256': row['after_step_sha256'], 'valid': valid}],
                    candidate_policy)
                self.assertEqual(bool(errors), fails, errors)


class BudgetUpgradeTests(unittest.TestCase):
    """A listed budget refusal may answer, audited by CI; nothing else may change."""

    def setUp(self):
        self.driver = checks.load_module('upgrade_test_driver', checks.ROOT / checks.DRIVER)
        self.refusal = {'status': 'refused', 'refusal': {'code': checks.BUDGET_REFUSAL,
                                                         'detail': 'operation work allowance exhausted'}}
        self.answer = {'status': 'answered', 'exact': True, 'marks': [],
                       'properties': {'volume': {'lo': '1', 'hi': '1'}}}
        files = {'answer.step': None, 'display.obj': None, 'print.stl': None, 'answer.igs': None}
        self.before = {'index': 194, 'request_sha256': 'a' * 64, 'family': 'mixed', 'owner': 'unknown',
                       'accepted_answer_sha256': 'b' * 64,
                       'accepted_manifest_sha256': self.driver.digest(self.refusal),
                       'accepted_files': files}
        self.after = copy.deepcopy(self.before)
        self.after.update(accepted_answer_sha256='c' * 64,
                          accepted_manifest_sha256=self.driver.digest(self.answer))
        self.after['accepted_files']['answer.step'] = 'd' * 64
        self.row = {'index': 194, 'class': checks.UPGRADE_CLASS, 'request_sha256': 'a' * 64,
                    'before_sha256': 'b' * 64, 'before_manifest_sha256': self.before['accepted_manifest_sha256'],
                    'after_sha256': 'c' * 64, 'after_manifest_sha256': self.after['accepted_manifest_sha256'],
                    'after_files': dict(self.after['accepted_files'])}
        self.other = {'index': 5, 'request_sha256': 'f' * 64, 'family': 'unite', 'owner': 'unknown',
                      'accepted_answer_sha256': '1' * 64, 'accepted_manifest_sha256': '2' * 64,
                      'accepted_files': {'answer.step': '3' * 64, 'display.obj': None,
                                         'print.stl': None, 'answer.igs': None}}

    def errors(self, base, candidate, rows, upgrade_format=None):
        upgrades = {'format': upgrade_format or checks.UPGRADE_FORMAT, 'upgrades': rows}
        corrections = {'format': 'catallax-replay-corrections/1', 'corrections': []}
        policy = checks.read(checks.ROOT / checks.RAW_POLICY)
        return checks.metadata_errors(base, candidate, corrections, policy, upgrades)

    def test_upgrade_metadata_table(self):
        listed = (([self.other, self.before], [self.other, self.after]), [self.row])
        self.assertEqual(self.errors(*listed[0], listed[1]), [])
        self.assertEqual(self.errors([self.other], [self.other], []), [])
        cases = []
        cases.append(('unlisted upgrade', [self.before], [self.after], []))
        proven = copy.deepcopy(self.before)
        proven['accepted_files']['answer.step'] = '4' * 64
        row = copy.deepcopy(self.row)
        changed = copy.deepcopy(self.after)
        cases.append(('listed change to a proven answer', [proven], [changed], [row]))
        for field in ('before_sha256', 'before_manifest_sha256', 'request_sha256', 'after_sha256'):
            row = copy.deepcopy(self.row)
            row[field] = '9' * 64
            cases.append((f'row {field} disagrees', [self.before], [self.after], [row]))
        moved = copy.deepcopy(self.after)
        moved['family'] = 'unite'
        cases.append(('family changed', [self.before], [moved], [self.row]))
        # Without a STEP file there is nothing for CI to audit independently.
        no_step = copy.deepcopy(self.after)
        no_step['accepted_files']['answer.step'] = None
        row = copy.deepcopy(self.row)
        row['after_files']['answer.step'] = None
        cases.append(('upgrade without a STEP file', [self.before], [no_step], [row]))
        row = copy.deepcopy(self.row)
        row['class'] = 'answer-to-answer'
        cases.append(('other class', [self.before], [self.after], [row]))
        other_changed = copy.deepcopy(self.other)
        other_changed['accepted_answer_sha256'] = '8' * 64
        cases.append(('unlisted proven answer beside a listed upgrade',
                      [self.other, self.before], [other_changed, self.after], [self.row]))
        cases.append(('duplicate rows', [self.before], [self.after], [self.row, self.row]))
        for name, base, candidate, rows in cases:
            with self.subTest(case=name):
                self.assertTrue(self.errors(base, candidate, rows), name)
        with self.subTest(case='wrong list format'):
            self.assertTrue(self.errors([self.before], [self.after], [self.row], 'catallax-replay-upgrades/2'))

    def test_upgrade_audit_requires_the_listed_step_and_a_passing_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            listed = root / 'listed' / 'answer.step'
            other = root / 'other' / 'answer.step'
            listed.parent.mkdir()
            other.parent.mkdir()
            listed.write_text('listed step')
            other.write_text('another case step')
            row = copy.deepcopy(self.row)
            row['after_files']['answer.step'] = checks.sha(listed)
            upgrades = {194: row}
            case = {'path': str(listed.parent), 'response': {'manifest': self.answer}}
            report = {'cases': [case], 'replay_indices': [194], 'upgrade_indices': [194]}
            seen = []

            def audit(path, manifest):
                seen.append((path, manifest))
                return True
            self.assertEqual(checks.upgrade_audit_errors(report, upgrades, audit), [])
            self.assertEqual(seen, [(listed, self.answer)])
            failures = {}
            failures['audit fails'] = (report, upgrades, lambda *_: False)
            failures['audit raises'] = (report, upgrades, lambda *_: 1 / 0)
            stale = copy.deepcopy(row)
            stale['after_files']['answer.step'] = '9' * 64
            failures['stale listed STEP'] = (report, {194: stale}, audit)
            moved = dict(report, cases=[dict(case, path=str(other.parent))])
            failures['another case STEP'] = (moved, upgrades, audit)
            listed_missing = root / 'missing'
            listed_missing.mkdir()
            absent = dict(report, cases=[dict(case, path=str(listed_missing))])
            failures['absent STEP'] = (absent, upgrades, audit)
            failures['unlisted audited index'] = (report, {}, audit)
            unreplayed = dict(report, upgrade_indices=[194, 7])
            failures['upgrade without a replayed case'] = (unreplayed, upgrades, audit)
            for name, (candidate, rows, check) in failures.items():
                with self.subTest(case=name):
                    self.assertTrue(checks.upgrade_audit_errors(candidate, rows, check), name)

    def test_upgrade_replay_manifests(self):
        self.assertEqual(checks.upgrade_manifest_errors(self.answer, self.refusal), [])
        other_refusal = {'status': 'refused', 'refusal': {'code': 'boolean.self-touching-shell'}}
        self.assertTrue(checks.upgrade_manifest_errors(self.answer, other_refusal))
        self.assertTrue(checks.upgrade_manifest_errors(self.answer, self.answer))
        self.assertTrue(checks.upgrade_manifest_errors(self.refusal, self.refusal))

    def test_replay_runs_the_base_verifier_for_a_listed_upgrade(self):
        driver = self.driver
        for fault in (None, 'base-answers', 'candidate-refuses'):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                inputs = [{'source': 'question', 'problem': {'format': 'test', 'position': 0}}]
                request = dict(max_work=100_000_000, **inputs[0], format=driver.CONTRACT)
                before = dict(self.before, index=0, request_sha256=driver.digest(request), family='mixed')
                after = dict(self.after, index=0, request_sha256=before['request_sha256'])
                row = dict(self.row, index=0, request_sha256=before['request_sha256'])
                checks.write(root / checks.INVENTORY, [after])
                checks.write(root / checks.CORRECTIONS, {'corrections': []})
                checks.write(root / checks.UPGRADES, {'format': checks.UPGRADE_FORMAT, 'upgrades': [row]})
                checks.write(root / 'inputs.json', inputs)
                (root / 'binary').write_text('candidate binary')
                (root / 'base-binary').write_text('base binary')
                checks.write(root / checks.RULES, {'corpus_sha256': checks.sha(root / 'inputs.json')})
                plan = {'identity': {'commit': 'candidate',
                                     'metadata_sha256': checks.sha(root / checks.INVENTORY),
                                     'corrections_sha256': checks.sha(root / checks.CORRECTIONS),
                                     'upgrades_sha256': checks.sha(root / checks.UPGRADES)},
                        'correction_base_answers': [before],
                        'selection': {'jobs': {'fmt': 1, 'build': 1, 'replay': 1}, 'selected_cases': 1,
                                      'replay_partitions': checks.replay_catalog([after])}}
                checks.write(root / 'plan.json', plan)
                args = SimpleNamespace(plan=root / 'plan.json', partition=1, inputs=root / 'inputs.json',
                                       binary=root / 'binary', base_binary=root / 'base-binary',
                                       out=root / 'out', summary=root / 'summary.json')

                def run(_driver, binary, requests, out, workers):
                    base = Path(binary).name == 'base-binary'
                    manifest = self.refusal if base != (fault == 'base-answers') else self.answer
                    if not base and fault == 'candidate-refuses':
                        manifest = self.refusal
                    target = Path(out) / 'case-0'
                    line = Path(requests).read_text().splitlines()[0]
                    checks.write(target / 'request.json', dict(json.loads(line), out=str(target)))
                    return {'cases': [{'path': str(target), 'cache_hit': False,
                                       'response': {'manifest': manifest, 'work': {'used': 1}}}],
                            'config': {'binary_sha256': checks.sha(binary)},
                            'execution': {'workers': workers, 'fresh': True,
                                          'worker_exits': {str(i): 0 for i in range(workers)}}}

                def answer_digest(_driver, case):
                    return 'b' * 64 if Path(case).parent.name == 'base' else 'c' * 64
                release = SimpleNamespace(OUTPUTS=[], answer_digest=answer_digest)
                with patch.object(checks, 'ROOT', root), patch.object(checks, 'git', return_value=b'candidate\n'), \
                        patch.object(checks.os, 'sched_getaffinity', return_value=set(range(4))), \
                        patch.object(checks, 'run_fresh', side_effect=run), \
                        patch.object(checks, 'load_module', side_effect=[driver, release]), \
                        contextlib.redirect_stdout(io.StringIO()):
                    if fault:
                        with self.assertRaises(ValueError):
                            checks.replay(args)
                    else:
                        checks.replay(args)
                        self.assertEqual(checks.read(args.summary)['verified_corrections'], [0])
                        self.assertEqual(checks.read(root / 'out/run.json')['upgrade_indices'], [0])


def workflow_step(name):
    """(shell, run script) of a named check.yml step, read as text."""
    lines = (checks.ROOT / '.github/workflows/check.yml').read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == f'- name: {name}')
    keys = len(lines[start]) - len(lines[start].lstrip()) + 2
    shell, run, i = None, None, start + 1

    def indent(line):
        return len(line) - len(line.lstrip())
    while i < len(lines) and (not lines[i].strip() or indent(lines[i]) >= keys):
        stripped = lines[i].strip()
        if stripped.startswith('shell:'):
            shell = stripped.split(':', 1)[1].strip()
        if stripped == 'run: |':
            body, i = [], i + 1
            while i < len(lines) and (not lines[i].strip() or indent(lines[i]) > keys):
                body.append(lines[i][keys + 2:])
                i += 1
            run = '\n'.join(body) + '\n'
            continue
        i += 1
    return shell, run


class RawStepWorkflowTests(unittest.TestCase):
    def test_direct_step_audit_failure_fails_the_job(self):
        # Run the step exactly as the hosted runner does, with a stand-in
        # audit interpreter, so a masked exit status cannot pass.
        shell, run = workflow_step('Direct no-healing STEP validity')
        self.assertIn('tools/checks/raw_step.py', run)
        command = (['bash', '--noprofile', '--norc', '-eo', 'pipefail'] if shell == 'bash'
                   else ['bash', '-e'])
        audits = {'passes': ("print('{\"raw_step_failures\": []}')", 0),
                  'returns nonzero': ("print('{\"raw_step_failures\": [\"x\"]}'); raise SystemExit(1)", 1),
                  'raises': ("raise RuntimeError('audit crashed')", 1)}
        for name, (body, expected) in audits.items():
            with self.subTest(audit=name), tempfile.TemporaryDirectory() as directory:
                temp = Path(directory)
                python = temp / 'occt' / 'bin' / 'python'
                python.parent.mkdir(parents=True)
                (temp / 'audit.py').write_text(body + '\n')
                python.write_text(f'#!/usr/bin/env bash\nexec python3 -B {temp / "audit.py"}\n')
                python.chmod(0o755)
                (temp / 'step.sh').write_text(run)
                result = subprocess.run([*command, str(temp / 'step.sh')], cwd=checks.ROOT,
                                        env={'PATH': '/usr/bin:/bin', 'RUNNER_TEMP': str(temp)},
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode != 0, bool(expected), result.stderr)
                if not expected:
                    self.assertIn('raw_step_failures', (temp / 'raw-step-summary.json').read_text())


class HangLimitTests(unittest.TestCase):
    def test_every_hang_guard_uses_the_one_ten_minute_limit(self):
        python = r'^HANG_LIMIT_SECONDS = (\d+)$'
        for path, pattern in [('.claude/skills/verify-catallax/scripts/control-catallax.py', python),
                              ('tools/step-regressions/test_step_import.py', python),
                              ('tools/edge-hunt/test_edge_hunt.py', python),
                              ('crates/catallax-cli/tests/exchange.rs', r'let limit = Duration::from_secs\((\d+)\);')]:
            with self.subTest(path=path):
                found = re.findall(pattern, (checks.ROOT / path).read_text(), re.M)
                self.assertEqual(found, ['600'])
        nextest = (checks.ROOT / '.config/nextest.toml').read_text()
        period, count = re.search(r'period = "(\d+)s", terminate-after = (\d+)', nextest).groups()
        self.assertEqual(int(period) * int(count), 600)
        self.assertIn('--profile ci', (checks.ROOT / '.github/workflows/check.yml').read_text())


if __name__ == '__main__':
    unittest.main()
