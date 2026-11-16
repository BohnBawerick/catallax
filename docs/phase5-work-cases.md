# Phase-5 work cases

The application allows 100M until phase 5; the target remains 20M. These completed
edits exceed the target and retain their measured stage work for later speed work.
No optimisation is part of their capability qualification.

| Case | Result | Work used | Evidence |
|---|---|---:|---|
| M10-004 | builds at 100M, over the 20M target | 26,777,693 | [Stage work and raw STEP checks](phase1/spline-boolean-top-cuts/REPORT.md) |
| M10-005 | builds at 100M, over the 20M target | 27,289,073 | [Stage work and raw STEP checks](phase1/spline-boolean-top-cuts/REPORT.md) |
| M10-006 | builds at 100M, over the 20M target | 25,945,726 | [Stage work and raw STEP checks](phase1/spline-boolean-top-cuts/REPORT.md) |

The slow loft-source admission test timed out at 1,200 seconds in CI. It stays
available with `--ignored`. The complete 100M top requests also have an explicit
slow regression. Required Python CI still checks all three 20M refusals and artifact
absence, in the `step-regressions-slow` job. The listed full-request counts come from the recorded source before the
review's narrower predispatch checks; those long requests were not repeated.
