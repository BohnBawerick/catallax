"""Publication must reject incomplete/different inputs and retain unknown outcomes."""
from collections import Counter
import unittest

from compare_import_census import compare, public_reference_manifest, HASHES, SELECTION


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.hashes = {f'{i:08}': f'{i:064x}' for i in range(3000)}
        self.native = [{'model_id': key, 'source_sha256': value, 'outcome': 'refused',
            'import_category': '', 'after_fix_cause': 'incidence_gap'} for key, value in self.hashes.items()]
        self.reference = [{'model_id': key, 'source_sha256': value, 'outcome': 'opened_other',
            'source_unchanged': True, 'opened': True} for key, value in self.hashes.items()]
        self.identity = {'complete': True, 'denominator': 3000, 'seed': 20261111,
            'source_hash_manifest_sha256': HASHES, 'selection_sha256': SELECTION,
            'source_revision': 'a' * 40, 'binary_sha256': 'b' * 64,
            'observer_sha256': 'd' * 64,
            'build_record': {'source_revision': 'a' * 40, 'binary_sha256': 'b' * 64,
                'observer_sha256': 'd' * 64, 'cargo_artifact': {'target': {'name': 'r8_census', 'kind': ['example']}}},
            'protocol_sha256': 'c' * 64, 'ocp_versions': ['8.0.1.1']}

    def run_comparison(self):
        native = {**self.identity, 'outcomes': dict(Counter(r['outcome'] for r in self.native))}
        reference = {**self.identity, 'outcomes': dict(Counter(r['outcome'] for r in self.reference))}
        return compare(self.native, native, self.reference, reference, self.hashes)

    def test_partial_and_timed_out_shapes_stay_distinct(self):
        self.native[0]['outcome'] = 'partial'
        self.reference[0].update(outcome='timeout', stage='shape_fix')
        rows, report = self.run_comparison()
        self.assertEqual(rows[0]['catallax_outcome'], 'partial')
        self.assertEqual(rows[0]['occt_outcome'], 'timeout')
        self.assertIsNone(rows[0]['occt_healed_whole_shape_valid'])
        self.assertEqual(report['outcome_pairs']['partial / timeout'], 1)

    def test_missing_or_duplicate_terminal_observation_prevents_publication(self):
        self.reference[-1] = dict(self.reference[0])
        with self.assertRaisesRegex(AssertionError, 'duplicate'):
            self.run_comparison()
        self.reference.pop()
        with self.assertRaises(AssertionError):
            self.run_comparison()

    def test_changed_source_bytes_prevent_comparison(self):
        self.reference[0]['source_sha256'] = 'd' * 64
        with self.assertRaisesRegex(AssertionError, 'source bytes'):
            self.run_comparison()

    def test_valid_surface_cannot_be_published_as_a_solid(self):
        self.reference[0].update(outcome='healed_valid_solid', stage='complete',
            healed_whole_shape_valid=True, solid_count=0, valid_solid_count=0)
        with self.assertRaises(AssertionError):
            self.run_comparison()

    def test_different_observer_build_prevents_attribution(self):
        self.identity['build_record']['source_revision'] = 'e' * 40
        with self.assertRaisesRegex(AssertionError, 'build record'):
            self.run_comparison()

    def test_reference_extract_omits_paths_and_preserves_unobserved_profiles(self):
        manifest = {**self.identity, 'source_root': '/home/private/sources',
            'selection_path': '/Users/private/selection.json', 'unknown_path': '/private/cache'}
        public = public_reference_manifest(manifest, self.reference)
        self.assertNotIn('/home/', str(public))
        self.assertNotIn('/Users/', str(public))
        self.assertNotIn('/private/', str(public))
        self.assertEqual(public['shape_fix_defaults_unobserved_files'], 3000)


if __name__ == '__main__':
    unittest.main()
