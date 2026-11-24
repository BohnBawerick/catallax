# Checkpoint 4

This report shares the [v0.5.0 cut and rerun commands](checkpoint-3.md#rerun-the-cut) with Checkpoint 3. It records the default-tier STEP census, phase 4 delivery and the work carried into phase 5. The [coverage recount](checkpoint-3-coverage.md) uses the 61 reference-kernel lines. Independent field capability sits outside that count.

## STEP census

The headline counts complete-source occurrence-preserving exports that a separate OpenCASCADE reader validates as solids. Full solid admission, valid exports of isolated prototypes and at least one observed valid solid are separate measures. The published OpenCASCADE reference is 2,701 healed-valid files of the same 3,000 (90.03%), using normal transfer and default ShapeFix. A direct parity comparison still needs a common protocol and a target.

The tagged collector uses the same frozen source hashes and selection as [Checkpoint 2](checkpoint-2.md#step-census), with two workers, 20M opening work, 20M aggregate export work, 60 seconds per observation and 6,144 MiB per observer. The published reference and the older acceptance baseline serve different purposes. The latter compares against 405 full admissions and 374 valid exports, with no per-body tier inventory for a same-body comparison.

The cut collected all 3,000 rows in 1,461.96 seconds (24.4 minutes), including 2.01 seconds to build. The observer binary SHA-256 is `5de86e97a52a5c6b070e2fcd2440505ff33287f2e4a9fca9606fc8b022078aa2`. The per-file CSV SHA-256 is `eb2c458d17d161a04ec92823ad0117b080ba9f0836d246aa4b4bdfce92a210e4`. Compared with the v0.4.0 cut, complete-source valid exports changed by +8 files and full solid admissions by +5 files. These are aggregate differences.

### Baseline preservation

The collector reports `acceptance_passed = false`. The frozen baseline admitted one solid in ABC 00002204. At v0.5.0, the same file exhausts opening at 19,999,998 charged units and retains one source-surface body. The recorded cause is `numbers.operation-budget`, with a body blocker at `request_limit`. The collector counts this loss of solid coverage under tier losses. Every frozen valid-export baseline file had a measured verdict. The run lost none of those exports. [Frozen baseline](../research/step-census-periodic/results.csv)

| Baseline comparison | Files |
| --- | --- |
| Tier losses | 1 |
| Frozen valid-export losses | 0 |
| Baseline items without a measurement | 0 |
| Unknown historical tier comparisons | 3,000 |

### File and body results

| Measure | Files | Share of 3,000 |
| --- | --- | --- |
| Catallax complete-source valid export | 767 | 25.57% |
| Catallax full solid admission | 889 | 29.63% |
| Catallax at least one observed valid exported solid | 1,148 | 38.27% |
| Catallax all isolated prototype exports valid | 767 | 25.57% |
| OpenCASCADE healed-valid reference | 2,701 | 90.03% |
| Checkpoint 1 legacy import/export, historical | 374 | 12.47% |

A full file takes the group of its least-checked admitted body. Checked numerical and sampled-native bodies are approximate. Body totals include partial files.

| Proof group | Full solid admissions | Complete-source valid exports | Bodies |
| --- | --- | --- | --- |
| Certified only | 187 | 177 | 572 |
| Checked numerical | 625 | 586 | 2,361 |
| Sampled-native tolerant | 77 | 4 | 154 |

The importer retained 13,630 source-surface bodies.

| Opening outcome | Files |
| --- | --- |
| Full solid admission | 889 |
| Partial or mixed | 409 |
| Refused | 40 |
| Source surfaces only | 1,662 |

| Export outcome | Files |
| --- | --- |
| Complete-source all-solids-valid | 767 |
| Export refused | 6 |
| Contains reader-invalid output | 202 |
| Not opened | 40 |
| Source surfaces without checked B-rep | 1,662 |
| Export timed out | 1 |
| Valid exported subset | 322 |

### First blockers and observation limits

| First body blocker | Cause | Resource-limited | Bodies | Files, overlapping |
| --- | --- | --- | --- | --- |
| body build | `not_solid` | no | 74 | 11 |
| body build | `numbers.operation-budget` | yes | 1 | 1 |
| body build | `resource_limit` | yes | 5 | 5 |
| curve fit | `import.tolerant-broken-geometry` | no | 52 | 52 |
| curve fit | `numbers.operation-budget` | yes | 604 | 604 |
| curve fit | `resource_limit` | yes | 4 | 4 |
| curve fit | `shapes.spline-exchange-data` | no | 8 | 8 |
| request limit | `resource_limit` | yes | 11,339 | 915 |
| sewing | `import.tolerant-broken-geometry` | no | 1 | 1 |
| sewing | `numbers.operation-budget` | yes | 4 | 4 |
| shell incidence | `import.tolerant-non-manifold` | no | 697 | 27 |
| surface fit | `numbers.operation-budget` | yes | 841 | 841 |

The any-valid verdict is unknown for 7 files. Valid-solid counts are lower bounds in 7 files. First-body-blocker inventories are absent for 40 files.

ABC 00004777 timed out during export. The six export refusals recorded one `numbers.operation-budget` and five `import.selected-document` at or near the export limit. No observer stopped for memory or reported a process error.

Unknown any-valid IDs are 00004777, 00005466, 00005469, 00006077, 00006301, 00006487, 00007518.

Body blocker file counts overlap. Work limits, timeouts, memory stops and unobserved siblings cannot establish geometric invalidity. Reader validity gives no common-tolerance source-fidelity result.

## Phase 4 delivery

| Roadmap item | Delivered | Open at the cut |
|---|---|---|
| 4.1 field preview | SolidDraft limits field previews to the field-only lattice and gyroid routes. It removed the proved/bounded lights and kept a refusal and working indicator. [Preview](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/344) [Status indicators](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/345) [Preview integration](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/358) | Preview and export keep separate checks. |
| 4.2 field builds P5 | P5 radius 2 builds its own field fillet. The recorded compare checked 3,192 vertices within 9.2e-8 mm. [P5 proof and controls](../M42A-FIELD-FILLETS.md#step-4-part-1-rim-contact) | The P5 result closes the named example. It does not close the full 4.2a comparison. |
| 4.2a independent fillets | Steps 1 to 6 have independent rolling-ball cases and typed refusals. Step 7 proves the spline-pair tube and includes dense checks, query times and a near-tangent control. The recorded 20M field proof and worker agreement predate basis-table charging. [Field fillet evidence](../M42A-FIELD-FILLETS.md) | The frozen compare table remains open. Tangent boss and narrow cone stop before the field check. The sphere rings lack B-rep display and the spline-plane row lacks a B-rep fillet. [Remaining comparisons](../M42A-FIELD-FILLETS.md#open-for-the-done-when-line) |
| 4.3 printing | P6 writes a closed STL with a 0.009997891320 mm bound. The independent reader checked its cavity paths and mesh. The reference takes 80,742,140 work, inside 100M and above the 20M target. [P6 print receipt](../milestones/4.3-printing.md#p6-at-001-mm) | The canonical fine P7 composite still lacks an exact grid-corner decision. Passing app gyroid scenarios cover their own inputs. [Gyroid scope](../milestones/4.3-printing.md#gyroid-scope) |
| 4.4 imports and scans | An OpenCASCADE export of the M27 loft exercises a proper rectangular spline trim in both kernels. A measured sphere scan becomes an analytic sphere with a maximum point residual of 0.051628 mm at a requested 0.06 mm. [STEP trim and scan](../milestones/4.4-imports-and-scans.md) | Independently sourced production trims, general UV, segmentation and other fitted scan families remain open. |
| 4.5 FEA | Parked on 13 Nov. [FEA note](../research/fea-first-beam.md) | Restart needs an owner decision. |
| SolidDraft field work | The release pin includes supported imported-body lattices and gyroids, planar graded walls and field grades. On the Mac, all 131 application scenarios passed their functional checks, 118 met every strict timing budget, and none timed out on STL. [Kernel release scope](../milestones/M5.md#050-phase-4-alpha-release) [SolidDraft validation](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/releases/tag/v0.8.0) | Nonplanar graded walls and wider import and scan domains remain open. The 13 app timing misses belong to phase 5. |

## Work for phase 5

The 20M target and 100M application allowance still apply. [Phase 5.1](../../ROADMAP.md#speed-and-simplification) carries the display work and size lines. [The work ledger](../phase5-work-cases.md) records successful requests and refused stages above 20M, including the P6 reference print and the thread sweep.

### Scoreboard work above 20M

The table uses the measured sample, one row per affected construction or edit. Each listed phase has its own allowance. A refused or exhausted phase records work spent. These counts belong to the speed phase.

| Target | Variant | Phase work and result |
| --- | --- | --- |
| P0 crossing-pipes | construction | measurement 88,701,624 (answered) |
| P0 crossing-pipes | edit | measurement 45,059,508 (answered) |
| ABC 00000405 | construction | measurement 27,007,529 (answered) |
| ABC 00000405 | edit | measurement 54,661,467 (answered) |
| P8-spring-threaded-boss | construction | build 55,626,285 (refused) |
| P8-spring-threaded-boss | edit | build 55,626,285 (refused) |
