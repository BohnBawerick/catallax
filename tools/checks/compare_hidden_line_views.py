#!/usr/bin/env python3
"""Compare hidden-line drawing answers from two builds of the same tests.

Build `hidden_line_views` and `drawing_curves` from each revision on the same
host with `cargo test -p catallax-app --test hidden_line_views --test drawing_curves --no-run`.
The baseline needs the test recorder applied to its unmodified view code.
Pass the four executable paths as printed by Cargo.
"""

import argparse
from itertools import zip_longest
import os
import subprocess


def records(binary, *filters):
    env = {**os.environ, "CATALLAX_RECORD_VIEWS": "1", "CATALLAX_THREADS": "1"}
    run = subprocess.run(
        [binary, "--nocapture", "--test-threads=1", *filters],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if run.returncode:
        raise SystemExit(f"{binary} failed:\n{run.stdout.decode(errors='replace')}")
    marker = b"VIEW_RECORD "
    rows = [line.split(marker, 1)[1] for line in run.stdout.splitlines() if marker in line]
    if not rows:
        raise SystemExit(f"{binary} emitted no view records")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_hidden")
    parser.add_argument("baseline_drawing")
    parser.add_argument("candidate_hidden")
    parser.add_argument("candidate_drawing")
    args = parser.parse_args()

    drawing_case = "face_loops_chain_and_strokes_keep_their_exact_ends_on_the_p1_bracket"
    baseline = records(args.baseline_hidden) + records(args.baseline_drawing, "--exact", drawing_case)
    candidate = records(args.candidate_hidden) + records(args.candidate_drawing, "--exact", drawing_case)
    for index, (old, new) in enumerate(zip_longest(baseline, candidate), 1):
        if old != new:
            test = (old or new).split(maxsplit=1)[0].decode()
            raise SystemExit(f"view record {index} differs in {test}")
    print(f"{len(baseline)} hidden-line view records match byte for byte")


if __name__ == "__main__":
    main()
