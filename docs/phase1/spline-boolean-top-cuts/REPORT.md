# Loft top cuts at the application allowance

All three original `catallax-step/1` requests build at 100M, over the 20M target.
Together with the eleven earlier STEP answers, this closes the recorded 14-request
set for milestone 1.3. Each result and its STEP are numerical and approximate at
1e-7 mm. No complete body approximation bound is established.

The caller accepts an unchanged loft constructor and a block spanning its width.
Every block vertex must be above the floor and strictly between the loft caps,
clear of the stated tolerance. A sampled original roof point must be strictly
inside the block. Disjoint and unresolved routing cases keep ordinary evaluation.
It reuses the exact finite intersection and imprint, then the separate numerical
material selection and retained-source STEP embedding. The selected collection has
no exact solid authority. Numerical admission, decimal-payload admission and disk
reread retain their existing checks. No geometry algorithm or optimisation changed.

The numerical classifier follows OCCT's
[nearest regular finite-ray transition](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/BRepClass3d/BRepClass3d_SClassifier.cxx)
and the bounded native curve/surface search in the
[SMLib port map](../../M25-INTERSECTION-PORT-SCOPE.md).

## Actual requests and files

Measured source commit `dd6767990347520ae3064886fab6b513afad1e73`. The [receipt](acceptance.json) records the source, binary hash, work and file checks. The later review change narrows predispatch only. These full-request measurements were not repeated after that routing change.

| Request | Whole-request work | Raw faces/edges/vertices | Physical holes | Volume mm3 |
|---|---:|---|---:|---:|
| M10-004 | 26,777,693 | 7/12/8 | 1 | 13.610476086914 |
| M10-005 | 27,289,073 | 11/24/16 | 1 | 113.610476086914 |
| M10-006 | 25,945,726 | 4/6/4 | 0 | 0.389523913086 |

OCP 8.0.1.1 opens each actual file through the raw StepToTopoDS path with healing
off. Each file passes BRepCheck, expected topology and physical-loop counts,
independent volume within 1e-6 mm3, six material witnesses and a later upper-slab
cut. The later edit is valid and matches its independent volume and material
witnesses. Candidate validation uses no TransferRoots, ShapeFix, sewing or
SameParameter. Source and emitted numerical evidence remain separate.

The work PC stopped for the checkpoint before any request completed. The Blade
then completed each request once with the same binary. The three concurrent
requests took about 136 to 140 minutes each on that run. These are observed proof
waits, not a latency qualification. No worker-count matrix or optimisation was run.

## Work per stage

The stages partition the complete request without counting parent spans twice.

| Stage | M10-004 | M10-005 | M10-006 |
|---|---:|---:|---:|
| face collection | 107 | 171 | 49 |
| imprint | 8,978,883 | 8,978,883 | 8,978,883 |
| intersection | 1,995,012 | 1,995,012 | 1,995,012 |
| material classification | 3,184,107 | 3,184,107 | 3,184,107 |
| source admission | 120 | 120 | 120 |
| numerical-step-definition | 201,966 | 203,831 | 196,438 |
| numerical-native-admission | 5,822,825 | 5,839,022 | 5,515,982 |
| numerical-step-payload | 6,489,909 | 6,978,479 | 5,976,637 |
| numerical-step-reread | 64,878 | 69,562 | 58,612 |
| prefix and dispatch | 39,886 | 39,886 | 39,886 |

All three cases join the [phase-5 work list](../../phase5-work-cases.md). The
application allowance is 100M until phase 5. The target remains 20M.

The review follow-up measured 9,187 routing work on the original stock, compared
with 2,540 in the recorded runs. That focused test checks only the read-only
predispatch predicate. It is not a new full-request measurement; the receipt
continues to identify the source of its 25.9M to 27.3M totals.

## Regression and reproduction

Default STEP CI runs all fourteen original requests at 20M, including worker
comparisons and stale-file replay. The three top requests must exhaust that
allowance, refuse with no result or STEP metadata, and leave no STEP artifact.
The STEP job has a 360-minute proof-wait guard for these required bounded checks;
the 20M request limit and target are unchanged.

The three 100M success requests have an additional explicit slow test. Run it with
`CATALLAX_SLOW_LOFT_TOP=1` and
`CATALLAX_VERIFY_BINARY=/path/to/catallax-verify`:

```sh
python3 -B -m unittest discover -s tools/step-regressions -p test_m10_spline_boolean.py -v
```

That opt-in also compares the top requests across the existing worker settings.
It was not repeated for this change. The qualification above uses the three
original requests completed once, their actual STEP files and the raw checks.
The separately ignored Rust source-admission test remains runnable with `--ignored`.
Only the additional 100M success matrix is opt-in. The required 20M refusal and
artifact checks remain in default CI. The existing Rust source-test exception
and the measured full-request times remain listed for phase 5.
