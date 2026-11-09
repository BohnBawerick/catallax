"""One acceptance rule for every kernel. Distances are sampled, not global."""
import math

TOLERANCES_MM = (0.1, 0.01, 0.001, 0.0001, 0.00001)
SOLVE_LIMIT_SECONDS = 120.0
DISTANCE_WIDTH_MM = "0.000001"
MAX_WORK = 100_000_000


def outcome(solve, check, tolerance, limit=SOLVE_LIMIT_SECONDS):
    elapsed = solve.get('wall_ms')
    if solve['status'] == 'timeout' or (elapsed is not None and elapsed > limit * 1000):
        return 'unsolved-time-limit'
    if solve['status'] != 'answered':
        return solve['status']
    if elapsed is None or not math.isfinite(elapsed):
        return 'unresolved'
    if check.get('valid') is False or check.get('witness_failures'):
        return 'not-done'
    lower, upper = check.get('distance_lower_mm'), check.get('distance_upper_mm')
    if lower is not None and math.isfinite(lower) and lower > tolerance:
        return 'not-done'
    if (check.get('valid') is True and check.get('complete_sample_check') is True
            and upper is not None and math.isfinite(upper) and upper <= tolerance):
        return 'done-sampled'
    return 'unresolved'


def sweep(solve, check, limit=SOLVE_LIMIT_SECONDS):
    return [{'tolerance_mm': t, 'solve_limit_seconds': limit,
             'outcome': outcome(solve, check, t, limit),
             'solve_wall_ms': solve.get('wall_ms'),
             'distance_lower_mm': check.get('distance_lower_mm'),
             'distance_upper_mm': check.get('distance_upper_mm')}
            for t in TOLERANCES_MM]
