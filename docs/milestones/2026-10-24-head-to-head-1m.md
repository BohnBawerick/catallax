# P0–P5 head-to-head at 1,000,000 work units

This first run measures commit `74bf31a8fae4e2d814cc9a4224dee401621401de`, based on `94ee7ebe72c986cf8944646b55d2643079962583`, against released SolidCraft v0.21.0. Each Catallax phase has a fixed allowance of 1,000,000 deterministic work units. This is the documented scoreboard setting, not the exchange's 100,000,000-unit default. The allowance is not increased to turn a refusal into a passing row.

Construction answers are not complete interactive results. Catallax constructs four of thirteen targets and their edits at this allowance. Every attempted Catallax display exhausts its own allowance, so none of these rows proves a completed displayed edit. Five target constructions exhaust work; the general cone and torus, P2 chamfer and P4 shell return typed coverage refusals. The reference constructs eleven targets; P1 and P5 fail in native Boolean commands, and its equal-cylinder/sphere edit also fails. Failed targets are retained rather than replaced by prefixes.

One warmup and three measured alternating pairs run on pinned core 0, with one thread. The host exposes 24 allowed cores and reports Linux under WSL2. Pair-start one-minute load averages range from 0.291 to 1.700; all before/after load triples remain in the raw data. These are pinned shared-machine observations, not an idle-host speed claim. Correctness checks occur outside the timing window.

## Measured rows

Times are medians in milliseconds. A successful-phase median requires every measured sample to answer that phase. Refusal time is separate. Counts cover the three measured samples, not warmups. Work counts, executed/recomputed operations, embedded reference tessellation and complete per-sample outcomes remain in the raw rows.

Property errors count discrepant quantities or centroid components, not distinct parts. They use an independent OpenCASCADE construction and the documented numeric diagnostic allowance, not a kernel acceptance tolerance. Membership is a fixed 216-point grid, not a global set equality proof. Name checks resolve every original result name but do not prove semantic face correspondence.

STEP columns distinguish files that open from checks that also pass topology and expected face count. Most reference STEP failures are face-count changes, not failed reads. Crossing-pipes exports actually fail OpenCASCADE topology validation in both variants. P2/P3 STEP and membership checks were unavailable in this original run because of the checker cast issue described below. A denominator of zero is unmeasured, never a pass.

| Part | Tool | Variant | Build ms | Display/export ms | Measurement ms | Refusal ms | Failed rebuilds | Incomplete | Property errors / checked | Membership errors / checked | STEP opened / checked / full pass | Names passed / checked |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: | --- | --- | --- | --- |
| P0/block-sphere.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 14.898 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/block-sphere.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 13.664 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/block-sphere.json | solidcraft | construction | 619.283 | 31.074 | 133.075 | none observed | 0 | 0 | 9 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P0/block-sphere.json | solidcraft | edit | 647.706 | 29.600 | 130.476 | none observed | 0 | 0 | 9 / 3 | 0 / 3 | 3 / 3 / 0 | 3 / 3 |
| P0/crossing-pipes.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 71.815 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/crossing-pipes.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 29.618 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/crossing-pipes.json | solidcraft | construction | 14025.278 | 4476.490 | 30809.803 | none observed | 0 | 0 | 6 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P0/crossing-pipes.json | solidcraft | edit | 7727.364 | 2478.897 | 16193.898 | none observed | 0 | 0 | 12 / 3 | 0 / 3 | 3 / 3 / 0 | 0 / 3 |
| P0/equal-cylinder-sphere.json | catallax | construction | 2.889 | unmeasured | 1.014 | display 30.344 | 0 | 0 | 0 / 3 | 0 / 3 | 3 / 3 / 3 | not applicable |
| P0/equal-cylinder-sphere.json | catallax | edit | 11.180 | unmeasured | 0.993 | display 36.259 | 0 | 0 | 0 / 3 | 0 / 3 | 3 / 3 / 3 | 3 / 3 |
| P0/equal-cylinder-sphere.json | solidcraft | construction | 144.130 | 12.712 | 38.547 | none observed | 0 | 0 | 6 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P0/equal-cylinder-sphere.json | solidcraft | edit | unmeasured | unmeasured | unmeasured | build upper bound 5082.142 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/general-cone.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 0.640 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/general-cone.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 0.439 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/general-cone.json | solidcraft | construction | 1314.354 | 791.114 | 3339.442 | none observed | 0 | 0 | 15 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P0/general-cone.json | solidcraft | edit | 1498.447 | 778.569 | 3188.006 | none observed | 0 | 0 | 15 / 3 | 0 / 3 | 3 / 3 / 0 | 0 / 3 |
| P0/general-torus.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 9.530 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/general-torus.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 0.147 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/general-torus.json | solidcraft | construction | 8335.539 | 3682.457 | 20611.575 | none observed | 0 | 0 | 15 / 3 | 0 / 3 | 3 / 3 / 3 | not applicable |
| P0/general-torus.json | solidcraft | edit | 8333.433 | 3733.029 | 22742.086 | none observed | 0 | 0 | 15 / 3 | 0 / 3 | 3 / 3 / 3 | 0 / 3 |
| P0/line-tangent.json | catallax | construction | 14.759 | unmeasured | 1.316 | display 47.604 | 0 | 0 | 0 / 3 | 0 / 3 | 3 / 3 / 3 | not applicable |
| P0/line-tangent.json | catallax | edit | 16.932 | unmeasured | 1.352 | display 45.313 | 0 | 0 | 0 / 3 | 0 / 3 | 3 / 3 / 3 | 3 / 3 |
| P0/line-tangent.json | solidcraft | construction | 62.498 | 16.208 | 59.728 | none observed | 0 | 0 | 6 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P0/line-tangent.json | solidcraft | edit | 579.735 | 15.714 | 59.417 | none observed | 0 | 0 | 6 / 3 | 0 / 3 | 3 / 3 / 0 | 0 / 3 |
| P0/near-miss.json | catallax | construction | 4.868 | unmeasured | 1.495 | display 51.784 | 0 | 0 | 0 / 3 | 0 / 3 | 3 / 3 / 3 | not applicable |
| P0/near-miss.json | catallax | edit | 4.388 | unmeasured | 1.400 | display 45.780 | 0 | 0 | 0 / 3 | 0 / 3 | 3 / 3 / 3 | 3 / 3 |
| P0/near-miss.json | solidcraft | construction | 62.380 | 17.268 | 62.846 | none observed | 0 | 0 | 12 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P0/near-miss.json | solidcraft | edit | 579.960 | 16.168 | 61.511 | none observed | 0 | 0 | 6 / 3 | 0 / 3 | 3 / 3 / 0 | 0 / 3 |
| P0/seam.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 14.000 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/seam.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 12.752 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/seam.json | solidcraft | construction | 50.862 | 14.103 | 57.118 | none observed | 0 | 0 | 9 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P0/seam.json | solidcraft | edit | 58.180 | 17.437 | 70.294 | none observed | 0 | 0 | 15 / 3 | 0 / 3 | 3 / 3 / 0 | 0 / 3 |
| P1-bracket.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 361.874 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P1-bracket.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 356.341 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P1-bracket.json | solidcraft | construction | unmeasured | unmeasured | unmeasured | build upper bound 11386.652 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P1-bracket.json | solidcraft | edit | unmeasured | unmeasured | unmeasured | none observed | 0 | 3 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P2-collar.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 0.301 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P2-collar.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 0.134 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P2-collar.json | solidcraft | construction | 75.562 | 17.733 | 67.852 | none observed | 0 | 0 | 6 / 3 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P2-collar.json | solidcraft | edit | 57.464 | 17.278 | 63.965 | none observed | 0 | 0 | 6 / 3 | 0 / 0 | 0 / 0 / 0 | 3 / 3 |
| P3-filleted-boss.json | catallax | construction | 10.377 | unmeasured | 4.410 | display 45.946 | 0 | 0 | 0 / 3 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P3-filleted-boss.json | catallax | edit | 9.778 | unmeasured | 4.625 | display 42.377 | 0 | 0 | 0 / 3 | 0 / 0 | 0 / 0 / 0 | 3 / 3 |
| P3-filleted-boss.json | solidcraft | construction | 939.210 | 17.987 | 511.197 | none observed | 0 | 0 | 9 / 3 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P3-filleted-boss.json | solidcraft | edit | 900.636 | 24.746 | 543.660 | none observed | 0 | 0 | 9 / 3 | 0 / 0 | 0 / 0 / 0 | 3 / 3 |
| P4-shelled-tray.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 0.207 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P4-shelled-tray.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 0.089 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P4-shelled-tray.json | solidcraft | construction | 150.024 | 1.938 | 7.772 | none observed | 0 | 0 | 9 / 3 | 0 / 3 | 3 / 3 / 0 | not applicable |
| P4-shelled-tray.json | solidcraft | edit | 147.065 | 1.949 | 7.974 | none observed | 0 | 0 | 9 / 3 | 0 / 3 | 3 / 3 / 0 | 3 / 3 |
| P5-flanged-housing.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 97.725 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P5-flanged-housing.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 96.968 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P5-flanged-housing.json | solidcraft | construction | unmeasured | unmeasured | unmeasured | build upper bound 6024.767 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P5-flanged-housing.json | solidcraft | edit | unmeasured | unmeasured | unmeasured | none observed | 0 | 3 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |

## Data and checker correction

- [Original raw rows](2026-10-24-head-to-head-1m-raw.json) retain all 52 pairs, including warmups, failures, phase logs and the original checker errors. Output paths and machine names were normalized in transfer. The old raw `summary` also retains refused-phase timing medians; the table above is regenerated from `pairs`, never from that summary.
- [Property rechecks](2026-10-24-head-to-head-1m-property-rechecks.json) contain fresh independent P2/P3 constructions with the corrected OpenCASCADE 8 cast. Every original claimed property is compared against these values, without changing any measured time or work count. Catallax's P3 properties agree; the reference P2/P3 requested properties differ at the diagnostic allowance. These rechecks do not manufacture STEP or membership evidence.
- [Derived table rows](2026-10-24-head-to-head-1m-summary.json) pin the raw and recheck SHA-256 values. They retain phase outcomes, refusal latency, work and operations rebuilt alongside the table fields.

The original checker used `TopoDS.Edge_s`, which the pinned OpenCASCADE 8.0.1.1 binding does not export. The current checker uses `TopoDS.Edge`. Focused independent P2 and P3 original/edit constructions now complete. P5's independent junction fillet still refuses in OpenCASCADE; neither tool constructed that full target in this run, so no P5 property agreement is asserted. The original aggregate `correctness_failures` includes checker errors and must not be read as a count of wrong geometric answers.

All checked successful Catallax properties agree with the independent values. Reference numeric discrepancies concern the requested measurement results, not a proof that every corresponding exact body is wrong. No overall speed ranking follows when displays or correctness gates do not complete.

## Earlier-answer proof

The [replay aggregate](2026-10-24-head-to-head-replay.json) records 477 regression inputs. Baseline versus candidate and candidate at one, two and four workers have identical manifests and strict STEP bytes, no differing IDs and no crashes. Each run answers 396 and refuses 81, including seven self-touching shells and 74 work-budget refusals. Work is unchanged. This proof compares the named commits above; it is not a timing ratio for the CAD targets.

| Build | Workers | Wall seconds | Deterministic work |
| --- | ---: | ---: | ---: |
| Baseline | 1 | 776.09 | 2,246,696,840 |
| Candidate | 1 | 774.28 | 2,246,696,840 |
| Candidate | 2 | 393.04 | 2,246,696,840 |
| Candidate | 4 | 249.00 | 2,246,696,840 |

## Rerun

Use the [benchmark README](../../.claude/skills/verify-catallax/scripts/head-to-head/README.md). It pins the released archive and executable by SHA-256, lists every authored reference script and separates construction, native edit, display export and requested measurement. Reproduce this table with `--max-work 1000000`; record a different allowance as a separate experiment. Nothing in this measurement adds work to a kernel or application build path.
