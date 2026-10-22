#!/usr/bin/env python3
"""Real-command regression tests for the verification driver."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[4]
DRIVER = Path(__file__).with_name('control-catallax.py')
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
EXAMPLES = ROOT / '.claude/skills/verify-catallax/references/examples'


@unittest.skipUnless(BINARY.exists(), 'build catallax-verify first')
class DriverTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='catallax-driver-')
        self.out = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def command(self, *args, code=0):
        completed = subprocess.run([str(DRIVER), *map(str, args), '--json'], capture_output=True, text=True)
        self.assertEqual(completed.returncode, code, completed.stderr + completed.stdout)
        return json.loads(completed.stdout)

    def run_case(self, name, *flags):
        return self.command('run', EXAMPLES / 'block.json', '--binary', BINARY, '--out', self.out / name, *flags)

    def test_cache_resume_corruption_and_step_diff(self):
        cache = self.out / 'cache'
        first = self.run_case('first', '--cache', cache)
        self.assertEqual(first['cache_hits'], 0)
        self.assertEqual(self.run_case('first')['cache_hits'], 1)
        self.assertEqual(self.run_case('second', '--cache', cache)['cache_hits'], 1)
        self.command('compare', self.out / 'first/run.json', self.out / 'second/run.json')
        second = json.loads((self.out / 'second/run.json').read_text())
        step = Path(second['cases'][0]['path']) / 'answer.step'
        step.write_bytes(step.read_bytes() + b'changed')
        diff = self.command('compare', self.out / 'first/run.json', self.out / 'second/run.json', code=1)
        self.assertEqual(diff['failures'][0]['cause'], 'step-bytes')
        self.assertEqual(self.run_case('second')['cache_hits'], 0)
        self.command('compare', self.out / 'first/run.json', self.out / 'second/run.json')

    def test_dry_run_does_not_create_artifacts(self):
        result = self.run_case('dry', '--fresh', '--dry-run')
        self.assertTrue(result['dry_run'])
        self.assertFalse((self.out / 'dry').exists())

    def test_stream_removes_stale_step_on_refusal(self):
        directory = self.out / 'same-output'
        question = json.loads((EXAMPLES / 'block.json').read_text())
        requests = [{'format': 'catallax-request/1', 'source': 'question', 'problem': question, 'out': str(directory), 'max_work': work} for work in (100_000_000, 0)]
        result = subprocess.run([str(BINARY), 'stream'], input='\n'.join(json.dumps(r) for r in requests) + '\n', capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        responses = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertEqual(responses[0]['manifest']['status'], 'answered')
        self.assertEqual(responses[1]['manifest']['refusal']['code'], 'numbers.operation-budget')
        self.assertFalse((directory / 'answer.step').exists())
        self.assertEqual(json.loads((directory / 'manifest.json').read_text())['status'], 'refused')

    def test_allowance_is_part_of_cache_key(self):
        cache = self.out / 'cache'
        self.run_case('full', '--cache', cache)
        limited = self.run_case('limited', '--cache', cache, '--max-work', 0)
        self.assertEqual(limited['cache_hits'], 0)
        self.assertEqual(limited['outcomes'], {'numbers.operation-budget': 1})
        self.assertFalse(list((self.out / 'limited').rglob('answer.step')))

    def test_worker_counts_duplicate_order_and_cleanup(self):
        corpus = self.out / 'corpus.jsonl'
        question = json.loads((EXAMPLES / 'block.json').read_text())
        rows = []
        for index in range(8):
            rows.append(dict(question, id=f'case-{index}'))
        rows.append(rows[0])
        corpus.write_text('\n'.join(json.dumps(q) for q in rows))
        for workers in (1, 2, 4):
            result = self.command('run', corpus, '--binary', BINARY, '--workers', workers, '--fresh', '--out', self.out / str(workers))
            self.assertEqual(result['cases'], 9)
        for workers in (2, 4):
            self.command('compare', self.out / '1/run.json', self.out / f'{workers}/run.json')
        self.command('cleanup', '--out', self.out, '--dry-run')
        self.command('cleanup', '--out', self.out)
        self.assertTrue((self.out / '1/run.json').exists())

    def test_interruption_resumes_complete_artifacts(self):
        corpus = self.out / 'resume.jsonl'
        question = json.loads((EXAMPLES / 'block.json').read_text())
        corpus.write_text('\n'.join(json.dumps(dict(question, id=f'resume-{i}')) for i in range(500)))
        destination = self.out / 'resume'
        args = ['run', str(corpus), '--binary', str(BINARY), '--out', str(destination), '--json']
        child = subprocess.Popen([str(DRIVER), *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            deadline = time.monotonic() + 30
            checkpoint = destination / 'checkpoint-0.jsonl'
            while not checkpoint.exists():
                if child.poll() is not None or time.monotonic() > deadline:
                    self.fail('worker did not publish a checkpoint')
                time.sleep(0.01)
            child.terminate()
            child.communicate(timeout=30)
            self.assertNotEqual(child.returncode, 0)
        finally:
            if child.poll() is None:
                child.terminate()
                child.communicate(timeout=30)
        result = self.command(*args[:-1])
        self.assertEqual(result['cases'], 500)
        self.assertGreater(result['cache_hits'], 0)

    def test_recipe_comparison_and_doctor(self):
        info = self.command('doctor', '--binary', BINARY)
        self.assertEqual(info['about']['request'], 'catallax-request/1')
        for source in ('recipe', 'compare'):
            self.command('run', EXAMPLES / 'recipe.json', '--source', source, '--binary', BINARY, '--out', self.out / source)
        report = json.loads((self.out / 'compare/run.json').read_text())
        self.assertEqual(report['cases'][0]['response']['manifest']['comparison']['verdict'], 'agree')
        self.assertFalse(list((self.out / 'compare').rglob('answer.step')))


if __name__ == '__main__':
    unittest.main()
