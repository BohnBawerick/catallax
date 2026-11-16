"""The memory-capped observer runner never turns a late exit into a result."""
import subprocess
import sys
import time
import unittest
from unittest.mock import patch

import rss_guard


class DeadlineTests(unittest.TestCase):
    def test_exit_after_the_deadline_is_a_timeout_even_inside_one_poll(self):
        # The poll interval is longer than the timeout, so an unbounded wait would
        # see this child exit successfully at 0.5 s, after the 0.2 s deadline.
        command = [sys.executable, '-c', 'import time; time.sleep(0.5)']
        started = time.monotonic()
        with self.assertRaises(subprocess.TimeoutExpired):
            rss_guard.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                          timeout=0.2, max_rss_mb=1024, poll=1.0)
        self.assertLess(time.monotonic() - started, 0.45)

    def test_exit_observed_after_the_deadline_is_a_timeout(self):
        # The child has exited before the wait, so the wait succeeds, but the
        # clock read after it is past the deadline.
        command = [sys.executable, '-c', 'pass']
        clock = iter([0.0, 0.5, 2.0])
        with patch('rss_guard.time.monotonic', lambda: next(clock, 2.0)):
            with self.assertRaises(subprocess.TimeoutExpired):
                rss_guard.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                              timeout=1.0, max_rss_mb=1024, poll=5.0)

    def test_exit_before_the_deadline_is_a_result(self):
        result = rss_guard.run([sys.executable, '-c', 'pass'], stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, timeout=30, max_rss_mb=1024)
        self.assertEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
