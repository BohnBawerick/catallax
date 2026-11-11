"""Runner routing and package failures must preserve every check."""
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'ci' / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


routing = load('runner_route', 'runner-route.py')
deps = load('install_edit_deps', 'install-edit-deps.py')


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.env = dict(TRUSTED_EVENT='true', BOAT_ALLOWED='true',
                        USE_LOCAL_RUNNERS='true', USE_BOAT='true',
                        LOCAL_OK_AT='1000', BOAT_OK_AT='1000', GITHUB_RUN_ATTEMPT='1')

    def test_fresh_leases_route_only_light_jobs_locally(self):
        result = routing.route(self.env, 1100)
        self.assertIn('ci-linux-local', result['light'])
        self.assertEqual(result['boat'], 'true')
        workflow = (ROOT / '.github/workflows/check.yml').read_text()
        scoreboard = workflow.split('  scoreboard:\n', 1)[1].split('\n  # Preserve', 1)[0]
        self.assertIn("needs.select.outputs.boat == 'true'", scoreboard)
        self.assertNotIn('outputs.light', scoreboard)

    def test_stale_missing_malformed_and_future_leases_fall_back(self):
        for stamp in ('0', '979', '1101', '', None, 'true', 'NaN'):
            with self.subTest(stamp=stamp):
                result = routing.route(dict(self.env, LOCAL_OK_AT=stamp, BOAT_OK_AT=stamp), 1100)
                self.assertEqual(result, {'light': ['ubuntu-latest'], 'boat': 'false'})

    def test_forks_disabled_switches_and_reruns_use_github(self):
        for change in (dict(TRUSTED_EVENT='false'), dict(GITHUB_RUN_ATTEMPT='2'),
                       dict(USE_LOCAL_RUNNERS='false', USE_BOAT='false')):
            self.assertEqual(routing.route(dict(self.env, **change), 1100),
                             {'light': ['ubuntu-latest'], 'boat': 'false'})
        self.assertEqual(routing.route(dict(self.env, BOAT_ALLOWED='false'), 1100)['boat'], 'false')

    def test_watchdog_metadata_and_full_rerun_paths(self):
        subprocess.run(['node', '--test', str(ROOT / 'ci/runner-fallback.test.cjs')], check=True)


class PackageTests(unittest.TestCase):
    def attempt(self, outcomes, perf=True):
        calls, delays = [], []

        def run(command):
            calls.append(command)
            return next(outcomes)

        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            result = deps.install(cache, 'resolute', run=run, sleep=delays.append, has_perf=lambda: perf)
            source = (cache / 'sources.list').read_text() if (cache / 'sources.list').exists() else ''
        return result, calls, delays, source

    def test_cache_needs_no_network_or_update(self):
        result, calls, delays, source = self.attempt(iter([True]))
        self.assertEqual(result, {'status': 'ready', 'source': 'cache'})
        self.assertEqual(len(calls), 1)
        self.assertIn('--no-download', calls[0])
        self.assertIn('linux-perf', calls[0])
        self.assertEqual((delays, source), ([], ''))

    def test_dns_failure_moves_to_an_independent_signed_mirror(self):
        result, calls, delays, source = self.attempt(iter([False, False, True, True]))
        self.assertEqual(result['source'], deps.MIRRORS[1])
        self.assertEqual(delays, [10])
        self.assertIn('signed-by=/usr/share/keyrings/ubuntu-archive-keyring.gpg', source)
        self.assertIn('resolute-security', source)
        self.assertIn('APT::Update::Error-Mode=any', calls[1])
        self.assertTrue(all('--allow-unauthenticated' not in call for call in calls))

    def test_restored_cache_uses_the_mirror_that_supplied_its_indices(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            source = cache / 'sources.list'
            source.write_text('deb https://mirrors.edge.kernel.org/ubuntu resolute main\n')
            calls = []
            result = deps.install(cache, 'resolute', run=lambda call: calls.append(call) or True,
                                  has_perf=lambda: True)
            self.assertEqual(result['source'], 'cache')
            self.assertIn(f'Dir::Etc::sourcelist={source}', calls[0])

    def test_archive_fetch_failure_also_changes_mirror(self):
        result, calls, delays, _ = self.attempt(iter([False, True, False, True, True]), perf=False)
        self.assertEqual(result['source'], deps.MIRRORS[1])
        self.assertNotIn('linux-perf', calls[-1])
        self.assertEqual(delays, [10])

    def test_exhaustion_is_infrastructure_and_never_success(self):
        result, calls, delays, _ = self.attempt(iter([False] * 4))
        self.assertEqual(result['status'], 'infrastructure-unavailable')
        self.assertEqual(delays, [10, 20])
        self.assertEqual(len(calls), 4)


if __name__ == '__main__':
    unittest.main()
