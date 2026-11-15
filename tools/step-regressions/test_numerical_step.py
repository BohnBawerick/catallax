"""Numerical STEP stays approximate across real driver workers and replay."""
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from test_m10_spline_boolean import cases

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'


def numerical_cases():
    m10 = next(problem['recipe'] for key, _, problem in cases() if key == 'M10-008')
    tiny = copy.deepcopy(m10)
    tiny['name'] = 'loft crossing'
    tiny['steps'][-2]['size'] = ['10', '10', '1']
    tiny['steps'][-2]['pose']['translation'] = ['-3', '1999999/1000000', '1']
    return [{'format': 'catallax-numerical-step/1', 'recipe': recipe,
             'tolerance_mm': '1/10000000'} for recipe in (m10, tiny)]


class NumericalStepTests(unittest.TestCase):
    def test_named_holes_keep_manifests_work_and_step_at_all_worker_counts(self):
        with tempfile.TemporaryDirectory(prefix='catallax-numerical-step-') as temporary:
            root = Path(temporary)
            questions = root / 'questions.jsonl'
            questions.write_text(''.join(json.dumps(problem) + '\n' for problem in numerical_cases()))
            baseline = None
            for workers in (1, 2, 4):
                out = root / str(workers)
                subprocess.run(['python3', str(DRIVER), 'run', str(questions), '--binary', str(BINARY),
                                '--out', str(out), '--workers', str(workers), '--max-work', '20000000',
                                '--fresh', '--json'], check=True, capture_output=True, text=True)
                run = json.loads((out / 'run.json').read_text())
                observed = []
                for case, counts in zip(run['cases'], [(9, 18, 12), (7, 12, 8)]):
                    response = case['response']
                    manifest = response['manifest']
                    self.assertEqual(manifest['format'], 'catallax-numerical-manifest/1')
                    self.assertEqual(manifest['status'], 'answered')
                    self.assertEqual(manifest['proof'], 'numerical')
                    self.assertIs(manifest['exact'], False)
                    self.assertEqual(manifest['category'], 'approximate')
                    self.assertEqual(tuple(manifest['result'][key] for key in ('faces', 'edges', 'vertices')), counts)
                    self.assertEqual(manifest['numerical']['complete_body_bound'], 'unestablished')
                    self.assertLessEqual(response['work']['used'], 20_000_000)
                    directory = Path(case['path'])
                    data = (directory / 'answer.step').read_bytes()
                    self.assertEqual(manifest['step']['sha256'], hashlib.sha256(data).hexdigest())
                    observed.append((manifest, response['work'], data))
                if baseline is None:
                    baseline = observed
                self.assertEqual(observed, baseline)
                if workers != 1:
                    subprocess.run(['python3', str(DRIVER), 'compare', str(root / '1/run.json'),
                                    str(out / 'run.json'), '--json'], check=True, capture_output=True, text=True)
            original = Path(run['cases'][0]['path'])
            replay = root / 'replay'
            response = json.loads(subprocess.check_output([
                str(BINARY), 'replay', str(original / 'question.json'), '--out', str(replay)], text=True))
            self.assertEqual(response['manifest'], baseline[0][0])
            self.assertEqual(response['work'], baseline[0][1])
            self.assertEqual((replay / 'answer.step').read_bytes(), baseline[0][2])


if __name__ == '__main__':
    unittest.main()
