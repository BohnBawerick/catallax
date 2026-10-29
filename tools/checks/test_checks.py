"""R1 acceptance cases G1-G6, S1-S3 and C3-C4 from the independent designer."""
import copy
from pathlib import Path
import unittest

import checks


class GateTests(unittest.TestCase):
    def setUp(self):
        self.plan = {'identity': {'commit': 'candidate', 'run_id': '123', 'run_attempt': '1',
                                 'corpus_sha256': 'corpus', 'metadata_sha256': 'metadata',
                                 'driver_sha256': 'driver', 'raw_policy_sha256': 'raw-policy'},
                     'selection': {'jobs': {'fmt': 1, 'clippy': 1, 'tests': 4, 'build': 1, 'replay': 1},
                                   'selected_cases': 477},
                     'declared': ['fmt', 'clippy', 'tests', 'build', 'replay', 'scoreboard']}
        self.reported = {name: {'result': 'success'} for name in self.plan['selection']['jobs']}
        self.reported['scoreboard'] = {'result': 'skipped'}
        self.records = []
        for name, partitions in self.plan['selection']['jobs'].items():
            for part in range(1, partitions + 1):
                record = dict(self.plan['identity'], job=name, partition=part,
                              binary_sha256='binary' if name in {'build', 'replay'} else None,
                              result='success')
                if name == 'replay':
                    record.update(workers=4, fresh=True, cache_hits=0, selected_cases=477)
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

    def test_g4_c4_identity_mismatches(self):
        for field in ('commit', 'run_id', 'binary_sha256', 'corpus_sha256',
                      'metadata_sha256', 'driver_sha256', 'raw_policy_sha256',
                      'cache_hits', 'workers', 'fresh'):
            with self.subTest(field=field):
                records = copy.deepcopy(self.records)
                replay = next(r for r in records if r['job'] == 'replay')
                replay[field] = {'cache_hits': 1, 'workers': 1, 'fresh': False}.get(field, 'stale')
                self.assertTrue(any(f'replay: stale {field}:' in e for e in self.errors(records=records)))

    def test_missing_job_ignores_candidate_needs(self):
        reported = copy.deepcopy(self.reported)
        del reported['clippy']
        self.assertIn('clippy: missing', self.errors(reported=reported))

    def test_c3_metadata_cannot_hide_changed_answer(self):
        base = [{'index': 7, 'accepted_answer_sha256':
                 '0203b8b177ca59e519192e31b591c51d3668e031ed6a83e0cae0539f84a8bdcd'}]
        changed = [{'index': 7, 'accepted_answer_sha256': 'changed'}]
        self.assertEqual(checks.metadata_errors(base, base), [])
        # C3(a) and C3(b) share this preflight: kernel changes cannot grant an
        # approval because the comparison takes only base and candidate metadata.
        self.assertIn('without base approval', checks.metadata_errors(base, changed)[0])


class SelectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules = checks.read(checks.ROOT / checks.RULES)
        cls.graph, cls.build = checks.cargo_graph(checks.ROOT)

    def route(self, paths, rules=None, event='pull_request', ref=''):
        return checks.select(paths, rules or self.rules, self.graph, self.build, event, ref)

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
            (['.github/workflows/check.yml'], set(self.graph), True),
            (['tools/checks/rules.json'], set(self.graph), True),
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
        errors = checks.decide(plan, {}, [])
        self.assertIn('extra: missing', errors)

    def test_s3_full_events(self):
        for event, ref in [('push', 'refs/heads/release/0.2.0'), ('workflow_dispatch', 'refs/heads/main')]:
            with self.subTest(event=event, ref=ref):
                selected = self.route(['docs/PLAN.md'], event=event, ref=ref)
                self.assertEqual(set(selected['packages']), set(self.graph))
                self.assertEqual(selected['jobs'], self.rules['jobs'])
                self.assertEqual(selected['selected_cases'], 477)


class RawStepPolicyTests(unittest.TestCase):
    def test_captain_ruling_table(self):
        policy = checks.read(checks.ROOT / checks.RAW_POLICY)
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


if __name__ == '__main__':
    unittest.main()
