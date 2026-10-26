"""Public synthetic evidence mutation tests; never read a qualification corpus."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from test_control import BINARY, EXAMPLES, load_driver
import qualification as q


def write(path, value):
    Path(path).write_text(json.dumps(value))


@unittest.skipUnless(BINARY.exists(), 'build catallax-verify first')
class QualificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='catallax-public-reuse-')
        cls.root = Path(cls.temp.name)
        cls.driver = load_driver()
        about = json.loads(subprocess.check_output([str(BINARY), 'about'], text=True))
        # These are explicitly synthetic build attestations, not measurements
        # of the compiler provenance of the test executable.
        cls.identity = {'format': 'catallax-build-evidence/1',
            'source': {'commit': json.loads(about['kernel'])['commit'], 'tree': 'b' * 40, 'dirty': False, 'submodules': []},
            'executable_sha256': cls.driver.file_hash(BINARY),
            'build': {'toolchain': {'rustc': 'synthetic rustc', 'cargo': 'synthetic cargo', 'linker': 'synthetic linker'},
                'dependencies': {k: 'c' * 64 for k in ('cargo_lock_sha256', 'resolved_sha256', 'native_sha256')},
                'options': {'command': ['synthetic-build'], 'env_sha256': 'd' * 64, 'config_sha256': 'e' * 64,
                            'target': 'synthetic-target', 'profile': 'release', 'features': []}}}
        write(cls.root / 'identity.json', cls.identity)
        block = json.loads((EXAMPLES / 'block.json').read_text())
        # Include duplicates, ordering and two different per-request allowances.
        full = {'format': 'catallax-question/1', 'source': 'question', 'problem': block, 'max_work': 100_000_000}
        limited = dict(full, max_work=0)
        (cls.root / 'corpus.jsonl').write_text('\n'.join(json.dumps(v) for v in (full, limited, full)))
        cls.args = argparse.Namespace(inputs=[str(cls.root / 'corpus.jsonl')], source='auto', adapter=None,
            baseline=str(BINARY), candidate=str(BINARY), out=str(cls.root / 'fresh'),
            baseline_identity=str(cls.root / 'identity.json'), reuse_baseline=None, fresh_audit=False,
            cores=None, max_work=100_000_000)
        result = q.qualify(cls.args, cls.driver)
        if not result['equal']:
            raise AssertionError(result)
        cls.archive = cls.root / 'fresh'
        cls.key = q.reuse_key(cls.args, cls.driver)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def setUp(self):
        self.case = tempfile.TemporaryDirectory(dir=self.root)
        self.work = Path(self.case.name)

    def tearDown(self):
        self.case.cleanup()

    def args_for(self, **updates):
        return argparse.Namespace(**{**vars(self.args), 'out': str(self.work / 'qualification'), **updates})

    def archive_copy(self):
        return Path(shutil.copytree(self.archive, self.work / 'archive'))

    def reseal(self, archive):
        record = json.loads((archive / 'archive.json').read_text())
        record['files'] = q.inventory(archive, self.driver)
        write(archive / 'archive.json', record)

    def mutate_case(self, archive, edit):
        run = archive / 'baseline2/run.json'
        report = json.loads(run.read_text())
        index = 0
        directory = run.parent / report['cases'][index]['path']
        edit(directory, report)
        # Update duplicate response records as an ordinary writer would do.
        response = json.loads((directory / 'response.json').read_text())
        for case in report['cases']:
            if case['key'] == report['cases'][index]['key']:
                case['response'] = response
        write(run, report)
        meta = json.loads((directory / 'artifact.json').read_text())
        meta['files'] = {p.name: self.driver.file_hash(p) for p in directory.iterdir() if p.name != 'artifact.json'}
        write(directory / 'artifact.json', meta)
        self.reseal(archive)

    def test_valid_archive_recomputed_and_candidates_always_execute_fresh(self):
        # Stored verdicts have no authority. Identical executable identities
        # still require three real fresh candidate jobs.
        archive = self.archive_copy()
        report = json.loads((archive / 'qualification.json').read_text())
        report['equal'] = False
        write(archive / 'qualification.json', report)
        self.reseal(archive)
        with patch.object(self.driver, 'run', wraps=self.driver.run) as calls:
            result = q.qualify(self.args_for(reuse_baseline=str(archive)), self.driver)
        self.assertEqual(result['status'], 'qualified', result)
        self.assertEqual(result['mode'], 'REUSED-BASELINE')
        self.assertEqual([c.args[0].workers for c in calls.call_args_list], [1, 2, 4])
        self.assertTrue(all(c.args[0].fresh and c.args[0].cache is None and
                            Path(c.args[0].out).name.startswith('candidate') for c in calls.call_args_list))
        self.assertEqual([r['evidence'] for r in result['runs']], ['reused'] * 3 + ['fresh'] * 3)
        self.assertTrue(all(r['historical_elapsed_ms'] is not None for r in result['runs'][:3]))
        self.assertFalse((Path(self.args_for().out) / 'archive.json').exists(), 'reused lineage never becomes fresh')
        baseline_report = json.loads((Path(self.args_for().out) / 'baseline1/run.json').read_text())
        self.assertEqual(self.driver.summarize(baseline_report)['timing'], 'historical')
        self.assertEqual(self.driver.summarize(baseline_report)['acceptance']['preview']['not_measured'], 3)
        candidate_report = json.loads((Path(self.args_for().out) / 'candidate1/run.json').read_text())
        self.assertEqual(self.driver.summarize(candidate_report)['timing'], 'observed')

    def test_each_key_leaf_change_rejects_reuse(self):
        def leaves(value, path=()):
            if isinstance(value, dict):
                for k, v in value.items():
                    yield from leaves(v, (*path, k))
            elif isinstance(value, list) and value:
                for i, v in enumerate(value):
                    yield from leaves(v, (*path, i))
            else:
                yield path
        paths = list(leaves(self.key))
        self.assertGreater(len(paths), 50)
        for path in paths:
            with self.subTest(path=path):
                changed = copy.deepcopy(self.key)
                parent = changed
                for k in path[:-1]:
                    parent = parent[k]
                parent[path[-1]] = 'changed'
                with self.assertRaises(q.EvidenceError):
                    q.validate_archive(self.archive, changed, self.driver)

    def test_unknown_or_incomplete_build_fields_disable_key(self):
        for path in [('source',), ('source', 'commit'), ('executable_sha256',),
                     ('build', 'toolchain'), ('build', 'dependencies'), ('build', 'options')]:
            for bad in (None, 'unknown', {}):
                with self.subTest(path=path, bad=bad):
                    data = copy.deepcopy(self.identity)
                    target = data
                    for k in path[:-1]:
                        target = target[k]
                    target[path[-1]] = bad
                    write(self.work / 'identity.json', data)
                    with self.assertRaises((ValueError, TypeError, KeyError)):
                        q.reuse_key(self.args_for(baseline_identity=str(self.work / 'identity.json')), self.driver)

    def test_real_key_construction_tracks_build_controls_corpus_and_environment(self):
        paths = [('source', 'tree'), ('source', 'commit'), ('executable_sha256',),
                 *[('build', 'toolchain', k) for k in self.identity['build']['toolchain']],
                 *[('build', 'dependencies', k) for k in self.identity['build']['dependencies']],
                 *[('build', 'options', k) for k in self.identity['build']['options']]]
        for path in paths:
            for unknown in (False, True):
                with self.subTest(path=path, unknown=unknown):
                    data = copy.deepcopy(self.identity)
                    parent = data
                    for k in path[:-1]:
                        parent = parent[k]
                    original = parent[path[-1]]
                    parent[path[-1]] = (None if unknown else ['changed'] if isinstance(original, list)
                                          else 'f' * len(original) if q.re.fullmatch('[0-9a-f]{40}|[0-9a-f]{64}', original) else 'changed')
                    write(self.work / 'changed-identity.json', data)
                    args = self.args_for(baseline_identity=str(self.work / 'changed-identity.json'))
                    try:
                        key = q.reuse_key(args, self.driver)
                    except (ValueError, TypeError, KeyError):
                        continue  # An unknown key also cannot authorize reuse.
                    self.assertNotEqual(key, self.key)
        for updates in ({'max_work': 1}, {'source': 'recipe'},
                        {'cores': str(min(os.sched_getaffinity(0)))}):
            self.assertNotEqual(q.reuse_key(self.args_for(**updates), self.driver), self.key)
        rows = (self.root / 'corpus.jsonl').read_text().splitlines()
        for index, changed in enumerate(([rows[1], rows[0], rows[2]], rows[:-1])):
            corpus = self.work / f'changed-{index}.jsonl'
            corpus.write_text('\n'.join(changed))
            self.assertNotEqual(q.reuse_key(self.args_for(inputs=[str(corpus)]), self.driver), self.key)
        with patch.dict(os.environ, {'CATALLAX_SYNTHETIC_DETERMINISM_CHANGE': '1'}):
            self.assertNotEqual(q.reuse_key(self.args_for(), self.driver), self.key)
        adapter = self.work / 'adapter.py'
        adapter.write_text('import json\ndef load(path):\n for line in path.read_text().splitlines():\n  v=json.loads(line)\n  yield {k:v[k] for k in ("source","problem","max_work")}\n')
        args = self.args_for(adapter=str(adapter))
        first = q.reuse_key(args, self.driver)
        adapter.write_text(adapter.read_text() + '# changed interpretation implementation\n')
        self.assertNotEqual(q.reuse_key(args, self.driver), first)

    def test_fresh_audit_unknown_key_missing_and_corrupt_archive_choose_six(self):
        for index, updates in enumerate(({'fresh_audit': True, 'reuse_baseline': str(self.archive)},
                {'baseline_identity': None, 'reuse_baseline': str(self.archive)},
                {'reuse_baseline': str(self.work / 'missing')},
                {'reuse_baseline': str(self.work)})):
            with self.subTest(updates=updates), patch.object(self.driver, 'run', wraps=self.driver.run) as calls:
                result = q.qualify(self.args_for(out=str(self.work / str(index)), **updates), self.driver)
                self.assertTrue(result['equal'], result)
                self.assertEqual(result['mode'], 'FULL-SIX')
                self.assertEqual(calls.call_count, 6)

    def test_run_policy_exit_and_completeness_mutations(self):
        mutations = [lambda r: r['execution'].update(fresh=False),
                     lambda r: r['execution'].update(workers=4),
                     lambda r: r['execution'].update(cores=[0]),
                     lambda r: r['execution']['worker_exits'].update({'0': 1}),
                     lambda r: r['execution'].update(worker_exits={}),
                     lambda r: r['cases'].pop(),
                     lambda r: r['cases'].reverse(),
                     lambda r: r['cases'][0].update(cache_hit=True),
                     lambda r: r['cases'][0].update(path='../escape'),
                     lambda r: r['config'].update(binary_sha256='f' * 64)]
        for i, mutate in enumerate(mutations):
            with self.subTest(mutation=i):
                archive = Path(shutil.copytree(self.archive, self.work / str(i)))
                run = archive / 'baseline2/run.json'
                data = json.loads(run.read_text())
                # Swap unequal allowance occurrences; reversing [A,B,A]
                # alone would be a no-op.
                if i == 6:
                    data['cases'][0], data['cases'][1] = data['cases'][1], data['cases'][0]
                else:
                    mutate(data)
                write(run, data)
                self.reseal(archive)
                with self.assertRaises((ValueError, OSError, KeyError)):
                    q.validate_archive(archive, self.key, self.driver)

    def test_missing_corrupt_extra_and_symlink_files(self):
        for i, mutation in enumerate(('missing', 'corrupt', 'extra', 'symlink')):
            with self.subTest(mutation=mutation):
                archive = Path(shutil.copytree(self.archive, self.work / str(i)))
                path = next((archive / 'baseline1/cases').glob('*/manifest.json'))
                if mutation == 'missing':
                    path.unlink()
                elif mutation == 'corrupt':
                    path.write_text('{')
                elif mutation == 'extra':
                    (path.parent / 'unexpected').write_text('extra')
                else:
                    path.unlink()
                    path.symlink_to(self.work / 'absent')
                with self.assertRaises((ValueError, OSError)):
                    q.validate_archive(archive, self.key, self.driver)

    def test_semantic_mutations_fail_even_with_resealed_hashes(self):
        def change_json(directory, filename, edit):
            path = directory / filename
            value = json.loads(path.read_text())
            edit(value)
            write(path, value)

        def coherent_response(directory, edit):
            response = json.loads((directory / 'response.json').read_text())
            edit(response)
            for filename, value in (('response.json', response), ('manifest.json', response['manifest']),
                                    ('operations.json', response['operation_log'])):
                write(directory / filename, value)

        mutations = [
            lambda d, r: change_json(d, 'question.json', lambda v: v.update(max_work=1)),
            lambda d, r: change_json(d, 'request.json', lambda v: v.update(max_work=1)),
            lambda d, r: change_json(d, 'manifest.json', lambda v: v.update(exact=False)),
            lambda d, r: (d / 'answer.step').write_bytes(b'changed STEP'),
            lambda d, r: (d / 'answer.step').unlink(),
            lambda d, r: coherent_response(d, lambda v: v['work'].update(used=v['work']['used'] + 1)),
            lambda d, r: coherent_response(d, lambda v: v['operation_log'].update(refusal={'code': 'changed'})),
            lambda d, r: coherent_response(d, lambda v: v['manifest'].update(properties={'interval': ['1.00000000000000000000001', '1.00000000000000000000002']})),
            lambda d, r: coherent_response(d, lambda v: v['manifest'].update(status='internal-error')),
            lambda d, r: (d / 'operations.json').unlink(),
        ]
        for i, edit in enumerate(mutations):
            with self.subTest(mutation=i):
                archive = Path(shutil.copytree(self.archive, self.work / str(i)))
                self.mutate_case(archive, edit)
                with self.assertRaises((ValueError, OSError, KeyError)):
                    q.validate_archive(archive, self.key, self.driver)

    def test_validation_defect_stops_without_retry_or_archive(self):
        with patch.object(self.driver, 'run', wraps=self.driver.run) as runs, \
             patch.object(q, 'gates', side_effect=q.EvidenceError('synthetic comparator defect')):
            result = q.qualify(self.args_for(), self.driver)
        self.assertFalse(result['equal'])
        self.assertEqual(result['status'], 'failed')
        self.assertEqual(runs.call_count, 6)
        self.assertIn('diagnose', result['action'])
        self.assertFalse((Path(self.args_for().out) / 'archive.json').exists())

    def test_intact_archive_with_validator_or_baseline_defect_is_not_retried(self):
        archive = self.archive_copy()
        self.mutate_case(archive, lambda d, r: (d / 'answer.step').write_bytes(b'changed'))
        with patch.object(self.driver, 'run', wraps=self.driver.run) as runs:
            result = q.qualify(self.args_for(reuse_baseline=str(archive)), self.driver)
        self.assertEqual(result['status'], 'failed')
        self.assertIn('diagnose', result['action'])
        self.assertEqual(runs.call_count, 0)

    def test_no_archive_chaining_and_no_resume(self):
        archive = self.archive_copy()
        record = json.loads((archive / 'archive.json').read_text())
        record['origin'] = 'reused-baseline'
        write(archive / 'archive.json', record)
        with self.assertRaises(q.EvidenceError):
            q.validate_archive(archive, self.key, self.driver)
        with self.assertRaises(q.EvidenceError):
            q.qualify(self.args_for(out=str(self.archive)), self.driver)

    def test_all_export_types_use_same_strict_gate(self):
        recipe = json.loads((EXAMPLES / 'recipe.json').read_text())
        rows = [{'format': 'catallax-question/1', 'source': source, 'max_work': 100_000_000, 'problem': problem}
                for source, problem in [('recipe', recipe), ('field', recipe),
                    ('iges', {'format': 'catallax-iges/1', 'recipe': recipe}),
                    ('print', {'format': 'catallax-print/1', 'recipe': recipe, 'tolerance': '1'})]]
        corpus = self.work / 'exports.jsonl'
        corpus.write_text('\n'.join(json.dumps(v) for v in rows))
        args = self.args_for(inputs=[str(corpus)])
        result = q.qualify(args, self.driver)
        self.assertTrue(result['equal'], result)
        out = Path(args.out)
        artifacts = list((out / 'baseline1/cases').glob('*/*'))
        for filename in ('answer.step', 'answer.igs', 'display.obj', 'print.stl'):
            self.assertTrue(any(p.name == filename for p in artifacts), filename)
        key = q.reuse_key(args, self.driver)
        q.validate_archive(out, key, self.driver)
        for filename in ('answer.igs', 'display.obj', 'print.stl'):
            archive = Path(shutil.copytree(out, self.work / filename))
            artifact = next((archive / 'baseline2/cases').glob('*/' + filename))
            artifact.write_bytes(artifact.read_bytes() + b'changed')
            metadata = json.loads((artifact.parent / 'artifact.json').read_text())
            metadata['files'][filename] = self.driver.file_hash(artifact)
            write(artifact.parent / 'artifact.json', metadata)
            self.reseal(archive)
            with self.assertRaises(q.EvidenceError):
                q.validate_archive(archive, key, self.driver)

    def test_full_decimal_numeric_intervals_are_not_rounded(self):
        first, second = self.work / 'first.json', self.work / 'second.json'
        first.write_text('{"lo":1.00000000000000000000001}')
        second.write_text('{"lo":1.00000000000000000000002}')
        self.assertEqual(json.loads(first.read_text()), json.loads(second.read_text()))
        self.assertNotEqual(q.read(first), q.read(second))


if __name__ == '__main__':
    unittest.main()
