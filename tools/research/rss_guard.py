"""Run one census observer under a resident-memory cap.

The work PC is shared and its out-of-memory guard picks victims by score, not
size, so one huge source can take down unrelated processes before the observer
itself. This runner polls the child's VmRSS and kills only that child when it
passes the cap. The caller records the file as out of memory and moves on.
A cap of 0 disables the check.
"""
import subprocess
import time


class MemoryExceeded(Exception):
    def __init__(self, peak_mb):
        super().__init__(f'observer exceeded its memory cap at {peak_mb} MiB')
        self.peak_mb = peak_mb


def rss_mb(pid):
    try:
        with open(f'/proc/{pid}/status') as status:
            for line in status:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) // 1024
    except (FileNotFoundError, ProcessLookupError):
        pass
    return 0


def run(command, *, stdout, stderr, timeout, env=None, max_rss_mb=0, poll=0.05):
    """Like subprocess.run(check=False), raising TimeoutExpired or MemoryExceeded."""
    process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=env)
    started, peak = time.monotonic(), 0
    while True:
        try:
            return subprocess.CompletedProcess(command, process.wait(timeout=poll))
        except subprocess.TimeoutExpired:
            pass
        peak = max(peak, rss_mb(process.pid))
        if max_rss_mb and peak > max_rss_mb:
            process.kill()
            process.wait()
            raise MemoryExceeded(peak)
        if time.monotonic() - started > timeout:
            process.kill()
            process.wait()
            raise subprocess.TimeoutExpired(command, timeout)
