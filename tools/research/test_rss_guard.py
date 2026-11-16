"""The memory-capped observer runner never turns a late exit into a result."""
import subprocess
import sys
import time
import unittest

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


if __name__ == '__main__':
    unittest.main()
