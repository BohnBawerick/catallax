"""Evidence controls: timeouts/errors cannot become admissions from checkpoints."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from import_census import observe


class ObservationTests(unittest.TestCase):
    def observe(self, source, imported, *, timeout=False, returncode=0, truncated=False):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / 'source.step'
            path.write_bytes(b'fixed observation input')
            raw = root / 'raw'
            raw.mkdir()

            def execute(_args, **kwargs):
                for value in (source, imported):
                    if value is not None:
                        kwargs['stdout'].write((json.dumps(value) + '\n').encode())
                if truncated:
                    kwargs['stdout'].write(b'{"stage":"import"')
                if timeout:
                    raise subprocess.TimeoutExpired(_args, 60)
                return subprocess.CompletedProcess(_args, returncode)

            baseline = {'outcome': 'refused', 'refusal_kind': 'not_solid',
                        'coverage_complete': 'True'}
            with patch('import_census.subprocess.run', execute):
                return observe(Path('/unused-binary'), path, '00000022',
                               hashlib.sha256(path.read_bytes()).hexdigest(),
                               raw, 60, baseline)

    @staticmethod
    def source(complete=True):
        return {'stage': 'faces', 'coverage_complete': complete,
                'source_parse': 'parsed', 'data_section_scan_complete': True,
                'data_sections': 1, 'unsupported_data_sections': 0,
                'all_face_kinds_buildable': True, 'needs_spline_faces': False}

    @staticmethod
    def admitted():
        return {'stage': 'import', 'outcome': 'admitted', 'bodies': 1,
                'import_category': 'certified', 'work_used': 100,
                'work_limit': 20_000_000}

    def test_unverified_source_coverage_is_partial(self):
        row = self.observe(self.source(False), self.admitted())
        self.assertEqual(row['outcome'], 'partial')
        self.assertFalse(row['coverage_complete'])

    def test_admission_checkpoint_cannot_erase_an_external_timeout(self):
        row = self.observe(self.source(), self.admitted(), timeout=True, truncated=True)
        self.assertEqual(row['outcome'], 'timeout')

    def test_admission_checkpoint_cannot_erase_a_process_failure(self):
        row = self.observe(self.source(), self.admitted(), returncode=1)
        self.assertEqual(row['outcome'], 'process_error')

    def test_numerical_zero_residual_cannot_be_recorded_exact(self):
        imported = {**self.admitted(), 'import_category': 'numerical approximate',
                    'exact': True, 'complete_body_bound': 'unestablished'}
        with self.assertRaises(AssertionError):
            self.observe(self.source(), imported)

    def test_numerical_method_and_pointwise_residual_keep_their_scope(self):
        imported = {**self.admitted(), 'import_category': 'numerical approximate',
                    'exact': False, 'complete_body_bound': 'unestablished',
                    'tolerances_mm': ['1/100'], 'max_checked_residual_mm': '0',
                    'source_admission_method': 'retained_shared_analytic_boundaries/1'}
        imported['source_intersection_checks'] = [
            {'status': 'disjoint', 'solid': '#1', 'faces': ['#2', '#3']},
            {'status': 'unresolved', 'solid': '#1', 'faces': ['#2', '#4'],
             'refusal': {'code': 'not_solid', 'detail': 'private source diagnostic'}}]
        row = self.observe(self.source(), imported)
        self.assertEqual(row['outcome'], 'full')
        self.assertEqual(row['import_category'], 'numerical approximate')
        self.assertIsNone(row['source_to_result_max_mm'])
        self.assertEqual(row['max_checked_residual_mm'], '0')
        self.assertEqual(row['source_admission_method'], imported['source_admission_method'])
        self.assertEqual(row['source_pair_disjoint_count'], 1)
        self.assertEqual(row['source_pair_unresolved'], [
            {'solid': '#1', 'faces': ['#2', '#4'], 'code': 'not_solid'}])
        imported['max_checked_residual_mm'] = '0.0101'
        with self.assertRaises(AssertionError):
            self.observe(self.source(), imported)


if __name__ == '__main__':
    unittest.main()
