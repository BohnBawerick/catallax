"""Run one census observer under a resident-memory cap.

The work PC is shared and its out-of-memory guard picks victims by score, not
size, so one huge source can take down unrelated processes before the observer
itself. This runner polls the child's VmRSS and kills only that child when it
passes the cap. The caller records the file as out of memory and moves on.
A cap of 0 disables the check.
"""
import os
import subprocess
import time

# Shared-machine floor: when the system's available memory falls under this many
# MiB, an observer above 1 GiB is stopped and recorded as out of memory.
MIN_AVAILABLE_MB = int(os.environ.get('P5_MIN_AVAILABLE_MB', '0'))


def available_mb():
    with open('/proc/meminfo') as info:
        for line in info:
            if line.startswith('MemAvailable:'):
                return int(line.split()[1]) // 1024
    return 1 << 30


class MemoryExceeded(Exception):
    def __init__(self, peak_mb):
        super().__init__(f'observer exceeded its memory cap at {peak_mb} MiB')
        self.peak_mb = peak_mb


_SAMPLES = {}
# Shared machines: with P5_CPU_CAP=1 the time cap counts the child's own CPU seconds, so a
# loaded host does not turn a slow wall clock into a timeout. A wall backstop of six times
# the cap still ends a child that hangs without using CPU.
CPU_CAP = os.environ.get('P5_CPU_CAP') == '1'
WALL_FACTOR = 6


def cpu_seconds(pid):
    try:
        if os.path.isdir('/proc'):
            with open(f'/proc/{pid}/stat') as stat:
                fields = stat.read().rsplit(')', 1)[1].split()
            return (int(fields[11]) + int(fields[12])) / os.sysconf('SC_CLK_TCK')
        text = subprocess.run(['ps', '-o', 'time=', '-p', str(pid)], capture_output=True,
                              text=True, timeout=5).stdout.strip()
        if not text:
            return 0.0
        days, _, clock = text.rpartition('-')
        total = 0.0
        for part in clock.split(':'):
            total = total * 60 + float(part)
        return total + (int(days) * 86400 if days else 0)
    except (ValueError, OSError, IndexError, subprocess.SubprocessError):
        return 0.0


def rss_mb(pid):
    if not os.path.isdir('/proc'):
        # macOS: ask ps, at most once a second for each child.
        now = time.monotonic()
        at, value = _SAMPLES.get(pid, (0.0, 0))
        if now - at < 1.0:
            return value
        try:
            text = subprocess.run(['ps', '-o', 'rss=', '-p', str(pid)], capture_output=True,
                                  text=True, timeout=5).stdout.strip()
            value = int(text) // 1024 if text else 0
        except (ValueError, OSError, subprocess.SubprocessError):
            value = 0
        _SAMPLES[pid] = (now, value)
        return value
    try:
        with open(f'/proc/{pid}/status') as status:
            for line in status:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) // 1024
    except (FileNotFoundError, ProcessLookupError):
        pass
    return 0


def run(command, *, stdout, stderr, timeout, env=None, max_rss_mb=0, poll=0.05):
    """Like subprocess.run(check=False), raising TimeoutExpired or MemoryExceeded.

    Every wait is bounded by the absolute deadline, so a child that exits after
    the timeout is a timeout, never a result. The child is killed and reaped on
    any error, including an interrupt in the caller.
    """
    deadline = time.monotonic() + timeout * (WALL_FACTOR if CPU_CAP else 1)
    cpu_checked = time.monotonic()
    process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=env)
    try:
        peak = 0
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise subprocess.TimeoutExpired(command, timeout)
            try:
                code = process.wait(timeout=min(poll, remaining))
            except subprocess.TimeoutExpired:
                pass
            else:
                # The wait can return an exit only seen after the deadline passed;
                # count any run observed past its deadline as timed out.
                if time.monotonic() > deadline:
                    raise subprocess.TimeoutExpired(command, timeout)
                return subprocess.CompletedProcess(command, code)
            peak = max(peak, rss_mb(process.pid))
            if CPU_CAP and time.monotonic() - cpu_checked >= 1.0:
                cpu_checked = time.monotonic()
                if cpu_seconds(process.pid) > timeout:
                    raise subprocess.TimeoutExpired(command, timeout)
            if max_rss_mb and peak > max_rss_mb:
                raise MemoryExceeded(peak)
            if MIN_AVAILABLE_MB and peak > 1024 and available_mb() < MIN_AVAILABLE_MB:
                raise MemoryExceeded(peak)
    except BaseException:
        if process.poll() is None:
            process.kill()
        process.wait()
        raise
