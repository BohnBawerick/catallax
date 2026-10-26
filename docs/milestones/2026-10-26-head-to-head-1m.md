# P0-P5 head-to-head refresh at 1,000,000 work units

Compared with [24 Oct](2026-10-24-head-to-head-1m.md), Catallax now constructs **five of thirteen** complete targets and their scalar edits, up from four. P2 collar is the new construction, edit and STEP-export success; its original displays, but its edited display still exhausts work. Equal-cylinder/sphere, line-tangent, near-miss and P3 now display both original and edit, where every display previously exhausted its allowance. Their STEP exports and name checks still pass. No other previously refused target now constructs. These are kernel CLI outcomes, not completed SolidDraft acceptance.

P1 still refuses at `slot-cut-1` with `numbers.operation-budget`, using 999,999 construction work and 999,997 edit work. The single first blocker is completing the full overlapping-slot Boolean subtraction chain within the fixed allowance, without replacing the target by a prefix or a different recipe. Until the full original and boss-radius edit construct, their display, requested properties, STEP and name survival cannot run. That is the next target-part lane supported by this run. It is a dependency diagnosis, not evidence that removing this blocker alone will make the later phases interactive.

## Snapshot and measurement conditions

- Source is main snapshot `866c9a30c4730ebe906cba352f61bf350b251652`. Benchmark executable SHA-256 is `f649eb931e004eb1750b41e5b3092ea661d39fbcb028cdd91be1d4d694741cc2`.
- Reference is released SolidCraft v0.21.0, executable SHA-256 `b334c87706833af67fe3b32e86c0c41938baa07f367a78e0a14ef6afc2de04be`. The existing downloader verified the pinned release archive and executable; no source build was substituted.
- Each Catallax phase has its own fixed 1,000,000-work allowance. Targets, authored edits, diagnostic tolerances and benchmark code were unchanged. One warmup and three measured alternating AB/BA pairs per target give 52 complete pairs. Each worker uses one thread pinned to allowed core 0.
- Timings come from a **shared loaded host**, an Intel Core i7-9750H running Linux. Twelve cores were allowed. Pair-start one-minute load ranged from 2.667 to 12.642; pair-end load ranged from 2.937 to 11.284. Every before/after one-, five- and fifteen-minute load triple is in the raw data. Work counts and outcomes are the main result; timings are secondary and not comparable as speed ratios against the earlier WSL2 host.
- The complete run took 5,325.24 seconds. An earlier attempt was interrupted by the outer command's one-hour deadline after 17 checkpointed pairs. It was discarded from the published results and the whole benchmark restarted in a fresh directory without that outer deadline. The published 52 pairs all belong to the complete second run. The driver's unchanged per-command timeout was 120 seconds. No 100M variant or 477-input replay was run.

## Remaining refusals

All three measured samples have identical Catallax phase work and outcomes. Counts below are consumed work, not an increased allowance. A charge may be rejected before the remaining units reach zero. Original and edit constructions were attempted independently.

| Target and phase | Original work | Edit work | Typed reason and blocked operation |
| --- | ---: | ---: | --- |
| P0/block-sphere build | 999,996 | 1,000,000 | `numbers.operation-budget`, `joined` |
| P0/crossing-pipes build | 999,979 | 1,000,000 | `numbers.operation-budget`, `cross` |
| P0/general-cone build | 845 | 845 | `intersect.traced-pair-required`, non-coaxial cone/cone at `joined` |
| P0/general-torus build | 4 | 4 | `intersect.traced-pair-required`, non-coaxial torus/torus at `joined` |
| P0/seam build | 999,998 | 1,000,000 | `numbers.operation-budget`, `joined` |
| P1 bracket build | 999,999 | 999,997 | `numbers.operation-budget`, `slot-cut-1` |
| P4 shelled tray build | 1,000,000 | 999,998 | `numbers.operation-budget`, `rounded` |
| P5 flanged housing build | 483,916 | 483,916 | `intersect.traced-pair-required`, non-coaxial cylinder/torus at `filleted` |
| P2 collar display | answered, 804,822 | 1,000,000 | Edit only, `numbers.operation-budget`, face 8 |

P4 has progressed from an unsupported shell to executing the final rounding intersection before work exhaustion. P5 has progressed from a work refusal to the typed junction-fillet tracing gap. Neither change is whole-target success. Cone and torus branch tracing remains unsupported. No display, measurement or STEP success is inferred from a refused construction.

## Successful construction and consumer work

Each cell is original / edit work. All five targets export STEP and pass all original face/edge name-resolution checks after the edit. P2's edit export and properties complete independently of its refused display.

| Target | Construction | Display | Requested properties | STEP export |
| --- | ---: | ---: | ---: | ---: |
| P0/equal-cylinder-sphere | 3,246 / 137,145 | 224,431 / 370,359 | 580 / 695 | 1,170 / 1,960 |
| P0/line-tangent | 243,114 / 371,668 | 449,546 / 449,546 | 784 / 784 | 2,340 / 2,340 |
| P0/near-miss | 4,202 / 4,167 | 453,453 / 453,453 | 784 / 784 | 4,432 / 4,432 |
| P2 collar | 2,751 / 2,751 | 804,822 / exhausted at 1,000,000 | 1,052 / 1,064 | 4,268 / 3,514 |
| P3 filleted boss | 25,081 / 26,523 | 477,581 / 483,011 | 3,849 / 3,863 | 10,931 / 10,931 |

The four fully displayed original/edit pairs finish build plus display below 1,000 ms in all three measured samples. This does not establish the whole requested click budget, the 500 ms Boolean target, the background comparison target or application interaction. The benchmark separately measures properties and STEP, and does not drive SolidDraft. P2's edit is explicitly not a completed displayed edit.

For P1, the original's first slot cut consumes 307,350 work; the next cut consumes 427,876 before refusing. The boss-radius edit uses 303,297 and 431,927 respectively. Both have 21 logged operations, with 20 completed before the refusal. The first slot cut already costs more than all preceding construction work combined. The [interactive contract](../PLAN.md#73-useful-latency-kept) requires the whole requested wait, not a quick refusal or a successful prefix. After complete P1 construction and edit become affordable, display, requested output, name survival and the actual application loop still need end-to-end proof.

## Correctness and historical comparability

The current benchmark checker differs from the 24 Oct workflow: it now compares properties and membership against analytic definitions. To retain the same historical property, membership and STEP checks without changing current tooling, the saved outputs were also checked with the unchanged corrected OpenCASCADE checker from commit `01c4326`. That is the historical `TopoDS.Edge` cast correction already documented on 24 Oct, not a tolerance change. All 208 original/edit, warmup/measured checker records are retained under `historical_correctness_rechecks`; the producing STEP SHA-256 is matched to the current-run record. These checks ran after all timing windows. Name checks are the original benchmark's checks, not a replacement.

The table below uses those historical checks for direct comparison. Every successful Catallax original and edit passes requested properties, sampled membership, STEP topology and expected face count. Every available edit name check passes. The fixed 216-point grid is sampled evidence, not global set equality; successful name lookup is not semantic face correspondence. A zero checked denominator means unmeasured.

The current analytic checks are also retained unchanged under each raw variant's `correctness` and in the summary's `current_analytic_checker_rows`. They report P3 membership discrepancies in every measured sample: one point for the original and two for the edit, totalling nine discrepancies. One is the outside point `[40, -10, 29.333333333333332]`, classified inside the imported STEP; the edit also disagrees at `[26.666666666666668, 0, 22]`. The historical grid passes on the same STEP digests, and both checkers pass P3 properties and STEP validity/face count. These different-grid outcomes remain unresolved, not erased by choosing the historical table. No unrestricted correctness or speed claim follows.

SolidCraft still constructs eleven originals. P1 and P5 fail native construction commands; their edits are not run. Equal-cylinder/sphere's native edit still fails. Reference numeric discrepancies, name losses and STEP face-count failures remain visible below. Crossing-pipes STEP topology remains invalid. The historical independent P5 junction blend also fails, so its property agreement is unmeasured, not a tool geometry failure. The complete current driver exits 1 with `correctness_failures: 126`; that aggregate mixes numeric, membership, STEP and name failures and is not a count of wrong parts.

## Data and reproduction

- [Raw run and historical checks](2026-10-26-head-to-head-1m-raw.json) retain every pair, load triple, phase log, refusal, interval endpoint, artifact path and STEP digest.
- [Derived summary](2026-10-26-head-to-head-1m-summary.json) pins the raw SHA-256, source and tooling hashes, verifier doctor identity, both checker summaries and all measured Catallax phase work samples.
- The existing [rerun guide](../../.claude/skills/verify-catallax/scripts/head-to-head/README.md) pins the reference release and environment. Build the unchanged `head_to_head` example at the named source, prepare `cadquery-ocp==8.0.1.1.0`, then run the command below. `$PROOF` is an external fresh directory. The additional historical checks execute `check.py` from `01c4326` against each saved original/edit STEP with the producer's expected face count, then use the current driver's unchanged `compare_properties` and `summarize` functions. No kernel, benchmark, target, site or old milestone file changed.

```sh
nice -n 10 python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py \
  head-to-head parts --binary "$PROOF/target/release/examples/head_to_head" \
  --reference "$PROOF/reference/solidcraft-cli" \
  --occt-python "$PROOF/occt/bin/python" --core 0 --max-work 1000000 \
  --repeats 3 --machine-class 'shared loaded build host' \
  --out "$PROOF/run-1m-complete" --json
```

## Phase timings and historical correctness checks

Medians are milliseconds over the three measured samples, excluding warmups. Successful-phase medians require all samples to answer. Refusal latency is separate. Property errors count discrepant quantities/components, membership errors count sampled points, and STEP full passes require opening, valid topology and the expected face count. Checks are outside the timing window. Timings are secondary shared-host observations.

| Part | Tool | Variant | Build ms | Display/export ms | Measurement ms | Refusal ms | Failed rebuilds | Incomplete | Not expressible | Property errors / checked sets | Membership errors / checked points | STEP opened / checked / full pass | Names passed / checked |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- | --- | --- | --- |
| P0/block-sphere.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 72.849 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/block-sphere.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 88.689 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/block-sphere.json | solidcraft | construction | 1858.762 | 98.527 | 463.194 | none observed | 0 | 0 | 0 | 9 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P0/block-sphere.json | solidcraft | edit | 1752.228 | 96.441 | 307.940 | none observed | 0 | 0 | 0 | 9 / 3 | 0 / 648 | 3 / 3 / 0 | 3 / 3 |
| P0/crossing-pipes.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 129.632 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/crossing-pipes.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 130.525 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/crossing-pipes.json | solidcraft | construction | 29036.981 | 8250.947 | 54977.419 | none observed | 0 | 0 | 0 | 6 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P0/crossing-pipes.json | solidcraft | edit | 13391.755 | 5629.540 | 32522.861 | none observed | 0 | 0 | 0 | 12 / 3 | 0 / 648 | 3 / 3 / 0 | 0 / 3 |
| P0/equal-cylinder-sphere.json | catallax | construction | 12.710 | 123.285 | 3.224 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | not applicable |
| P0/equal-cylinder-sphere.json | catallax | edit | 37.909 | 498.158 | 3.778 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | 3 / 3 |
| P0/equal-cylinder-sphere.json | solidcraft | construction | 525.079 | 60.568 | 131.695 | none observed | 0 | 0 | 0 | 6 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P0/equal-cylinder-sphere.json | solidcraft | edit | unmeasured | unmeasured | unmeasured | build upper bound 17776.224 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/general-cone.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 2.291 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/general-cone.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 2.061 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/general-cone.json | solidcraft | construction | 2737.042 | 2028.001 | 8590.018 | none observed | 0 | 0 | 0 | 15 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P0/general-cone.json | solidcraft | edit | 3226.810 | 1453.290 | 6626.881 | none observed | 0 | 0 | 0 | 15 / 3 | 0 / 648 | 3 / 3 / 0 | 0 / 3 |
| P0/general-torus.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 0.745 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/general-torus.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 0.422 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/general-torus.json | solidcraft | construction | 19511.332 | 8380.637 | 50079.003 | none observed | 0 | 0 | 0 | 15 / 3 | 0 / 648 | 3 / 3 / 3 | not applicable |
| P0/general-torus.json | solidcraft | edit | 15718.497 | 8115.159 | 51883.037 | none observed | 0 | 0 | 0 | 15 / 3 | 0 / 648 | 3 / 3 / 3 | 0 / 3 |
| P0/line-tangent.json | catallax | construction | 32.577 | 146.377 | 3.979 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | not applicable |
| P0/line-tangent.json | catallax | edit | 44.350 | 166.914 | 4.494 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | 3 / 3 |
| P0/line-tangent.json | solidcraft | construction | 198.600 | 65.754 | 168.377 | none observed | 0 | 0 | 0 | 6 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P0/line-tangent.json | solidcraft | edit | 1584.219 | 56.711 | 202.859 | none observed | 0 | 0 | 0 | 6 / 3 | 0 / 648 | 3 / 3 / 0 | 0 / 3 |
| P0/near-miss.json | catallax | construction | 10.451 | 139.285 | 3.192 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | not applicable |
| P0/near-miss.json | catallax | edit | 7.657 | 116.816 | 2.797 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | 3 / 3 |
| P0/near-miss.json | solidcraft | construction | 166.082 | 39.537 | 129.962 | none observed | 0 | 0 | 0 | 12 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P0/near-miss.json | solidcraft | edit | 1262.782 | 42.427 | 135.056 | none observed | 0 | 0 | 0 | 6 / 3 | 0 / 648 | 3 / 3 / 0 | 0 / 3 |
| P0/seam.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 66.007 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P0/seam.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 55.430 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P0/seam.json | solidcraft | construction | 86.987 | 28.461 | 121.954 | none observed | 0 | 0 | 0 | 9 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P0/seam.json | solidcraft | edit | 113.017 | 39.455 | 142.156 | none observed | 0 | 0 | 0 | 15 / 3 | 0 / 648 | 3 / 3 / 0 | 0 / 3 |
| P1-bracket.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 699.165 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P1-bracket.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 648.232 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P1-bracket.json | solidcraft | construction | unmeasured | unmeasured | unmeasured | build upper bound 18901.767 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P1-bracket.json | solidcraft | edit | unmeasured | unmeasured | unmeasured | none observed | 0 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P2-collar.json | catallax | construction | 0.809 | 174.539 | 3.755 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | not applicable |
| P2-collar.json | catallax | edit | 0.512 | unmeasured | 3.069 | display 220.384 | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | 3 / 3 |
| P2-collar.json | solidcraft | construction | 152.382 | 37.476 | 152.027 | none observed | 0 | 0 | 0 | 6 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P2-collar.json | solidcraft | edit | 109.567 | 38.063 | 127.495 | none observed | 0 | 0 | 0 | 6 / 3 | 0 / 648 | 3 / 3 / 0 | 3 / 3 |
| P3-filleted-boss.json | catallax | construction | 33.390 | 242.112 | 9.693 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | not applicable |
| P3-filleted-boss.json | catallax | edit | 26.914 | 284.147 | 13.583 | none observed | 0 | 0 | 0 | 0 / 3 | 0 / 648 | 3 / 3 / 3 | 3 / 3 |
| P3-filleted-boss.json | solidcraft | construction | 2373.193 | 46.563 | 1849.545 | none observed | 0 | 0 | 0 | 9 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P3-filleted-boss.json | solidcraft | edit | 2675.905 | 73.800 | 1510.952 | none observed | 0 | 0 | 0 | 9 / 3 | 0 / 648 | 3 / 3 / 0 | 3 / 3 |
| P4-shelled-tray.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 261.337 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P4-shelled-tray.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 228.487 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P4-shelled-tray.json | solidcraft | construction | 315.861 | 3.921 | 19.192 | none observed | 0 | 0 | 0 | 9 / 3 | 0 / 648 | 3 / 3 / 0 | not applicable |
| P4-shelled-tray.json | solidcraft | edit | 318.760 | 4.231 | 16.768 | none observed | 0 | 0 | 0 | 9 / 3 | 0 / 648 | 3 / 3 / 0 | 3 / 3 |
| P5-flanged-housing.json | catallax | construction | unmeasured | unmeasured | unmeasured | build 364.079 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P5-flanged-housing.json | catallax | edit | unmeasured | unmeasured | unmeasured | build 345.522 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
| P5-flanged-housing.json | solidcraft | construction | unmeasured | unmeasured | unmeasured | build upper bound 14393.861 | 3 | 0 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | not applicable |
| P5-flanged-housing.json | solidcraft | edit | unmeasured | unmeasured | unmeasured | none observed | 0 | 3 | 0 | 0 / 0 | 0 / 0 | 0 / 0 / 0 | 0 / 0 |
