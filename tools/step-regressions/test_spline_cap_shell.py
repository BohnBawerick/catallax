#!/usr/bin/env python3
"""Spline cap shells preserve the source family's explicit STEP refusal."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures'


class SplineCapShellStepTests(unittest.TestCase):
    def test_unsupported_polynomial_export_publishes_no_partial_step(self):
        for name in ['spline-cap-shell', 'spline-two-cap-shell', 'spline-cut-cap-shell']:
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix='catallax-spline-cap-shell-') as temporary:
                out = Path(temporary) / 'run'
                result = subprocess.run([
                    str(DRIVER), 'run', str(FIXTURES / (name + '.json')),
                    '--source', 'recipe', '--binary', str(BINARY), '--out', str(out),
                    '--fresh', '--workers', '1', '--json',
                ], capture_output=True, text=True, check=True)
                self.assertEqual(json.loads(result.stdout)['outcomes'], {'step.native-bspline-curve-unproven': 1})
                case = json.loads((out / 'run.json').read_text())['cases'][0]
                self.assertEqual(case['response']['manifest']['refusal']['code'], 'step.native-bspline-curve-unproven')
                self.assertFalse((Path(case['path']) / 'answer.step').exists())


if __name__ == '__main__':
    unittest.main()
