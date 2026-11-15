"""Evidence controls: timeouts/errors cannot become admissions from checkpoints."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from import_census import observe, public_manifest


class ObservationTests(unittest.TestCase):
    def test_public_manifest_omits_machine_paths_and_unknown_private_fields(self):
        report = {'complete': True, 'source_revision': 'a' * 40,
            'source_root': '/home/private/corpus', 'future_private_path': '/Users/private/scratch',
            'build_record': {'binary_sha256': 'b' * 64, 'artifact_path': '/home/private/r8_census',
                'cargo_artifact': {'executable': '/home/private/r8_census',
                    'filenames': ['/home/private/r8_census'], 'package_id': 'path+file:///home/private/crate',
                    'target': {'name': 'r8_census', 'kind': ['example'], 'src_path': '/Users/private/crate.rs'}}}}
        public = public_manifest(report)
        self.assertTrue(public['complete'])
        self.assertEqual(public['build_record']['binary_sha256'], 'b' * 64)
        self.assertEqual(public['build_record']['cargo_artifact']['target']['name'], 'r8_census')
        for private in ['/home/', '/Users/', 'file://', 'future_private_path']:
            self.assertNotIn(private, json.dumps(public))

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

    def test_admitted_body_kinds_survive_collection(self):
        for kinds in [{'surface': 1}, {'solid': 1, 'surface': 2}]:
            imported = {**self.admitted(), 'body_kinds': kinds,
                        'bodies': sum(kinds.values())}
            row = self.observe(self.source(), imported)
            self.assertEqual(row.get('body_kinds'), kinds)

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
            {'status': 'checked_no_crossing', 'solid': '#1', 'faces': ['#2', '#3']},
            {'status': 'unresolved', 'solid': '#1', 'faces': ['#2', '#4'],
             'refusal': {'code': 'not_solid', 'detail': 'private source diagnostic'}}]
        imported['source_topology_repairs'] = [{'solid': '#1', 'face': '#2',
            'source_rim_edges': ['#3', '#4'], 'source_vertices': ['#5', '#6'],
            'derived_seam_edge': 'derived:#2:periodic-seam',
            'method': 'paired_cylinder_generatrix/1'}]
        row = self.observe(self.source(), imported)
        self.assertEqual(row['source_topology_repair_count'], 1)
        self.assertEqual(json.loads(row['source_topology_repairs']), imported['source_topology_repairs'])
        self.assertEqual(row['outcome'], 'full')
        self.assertEqual(row['import_category'], 'numerical approximate')
        self.assertIsNone(row['source_to_result_max_mm'])
        self.assertEqual(row['max_checked_residual_mm'], '0')
        self.assertEqual(row['source_admission_method'], imported['source_admission_method'])
        self.assertEqual(row['source_pair_checked_count'], 1)
        self.assertEqual(row['source_pair_unresolved'], [
            {'solid': '#1', 'faces': ['#2', '#4'], 'code': 'not_solid'}])
        imported['max_checked_residual_mm'] = '0.0101'
        with self.assertRaises(AssertionError):
            self.observe(self.source(), imported)


if __name__ == '__main__':
    unittest.main()
