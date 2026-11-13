# M12 no-growth evidence

The no-growth claim covers the supported later-feature chains below. The
[M12 exit rerun](milestones/M12.md) at main
`f112dbd853b5ecc3a589df9870e784a86f02d576` refreshes these measurements after
M10. It does not qualify arbitrary cuts through fitted faces or finish M12.

For these supported operations, a rebuild measures each fitted face against
its retained definition. It never adds a prior build's deviation to the new
bound. An edit that changes the exact target can change the measured error;
the final bound must match a direct fit of that same target and meet the
unchanged local cap.

## What main already proves

- [Certified inward elliptic offsets](https://github.com/BohnBawerick/catallax/pull/53)
  retain the original exact law. In `certified_prism.rs`,
  `a_second_offset_refits_the_original_law_instead_of_offsetting_old_controls`
  checks a second offset, and
  `axial_planar_trim_and_later_offset_keep_the_original_exact_profile` checks
  a trim followed by another offset.
- [Fit marks on later Boolean results](https://github.com/BohnBawerick/catallax/pull/488)
  state the worst measured bound over the fitted faces the result keeps.
  `later_cuts.rs` checks that a disjoint cut keeps and remeasures all four P5
  rolling-envelope faces. It also checks that discarded fitted faces leave no
  fit mark on the result.

These implementations remain unchanged. The test makes their no-growth
claim explicit across longer chains and worker counts.

## Repeatable check

```sh
cargo test --release -p catallax-features --test no_growth -- --nocapture
```

`recorded_error_does_not_accumulate_at_one_two_or_four_workers` rebuilds the
canonical P1 bracket, P2 collar, P3 filleted boss, P4 shelled tray and P5 flanged
housing. Each part then subtracts a disjoint tool, reflects in an oblique
plane, reflects back, and subtracts the tool again. Each build retains the
20,000,000-unit application allowance.

At each supported stage the test checks topology and entity counts. It
remeasures every fitted P5 carrier against its retained definition, checks its
recorded bound and the body's fit mark, and requires the same bound as before
the later features. P1-P4 have no fitted carriers in these chains. The test
requires that fact explicitly, so they cannot pass by losing a fitted face.

The offset case compares four successive inward offsets of 1/16 mm with one
inward offset of 1/4 mm on the checked oblique-circle prism. It verifies the
complete final certificate, independently measures every intermediate carrier
against its exact offset law, and requires each region's recorded bound to
equal that measurement. The final bound must equal the direct fit's bound.
Intermediate fits need not have equal errors because their exact targets
differ. The claim is that the final fit does not add its ancestors' errors.

All cases run freshly at 1, 2 and 4 kernel threads. The test compares measured
bounds, work counts, entity counts, face origins, edge names and senses,
128-bit vertex enclosures, and fitted carrier knots, weights and control-point
enclosures across those runs. It prints the attained rational bounds and build
work for review. Raw evidence retained privately.

The measured values below use nanometres for readability. Assertions compare
the exact rational bounds, without rounding them.

| Supported case | Stages checked | Fitted faces per stage | Measured upper bound | Build work units |
| --- | ---: | ---: | --- | ---: |
| P1 bracket | 5 | 0 | No fitted deviation | 13,744,536 |
| P2 collar | 5 | 0 | No fitted deviation | 95,969 |
| P3 filleted boss | 5 | 0 | No fitted deviation | 226,065 |
| P4 shelled tray | 5 | 0 | No fitted deviation | 1,836,197 |
| P5 flanged housing | 5 | 4 | About 4.170328 nm on each face, unchanged | 15,749,769 |
| Direct 1/4 mm offset | 1 | 3 certificate regions | About 7.322819 nm | 798,205 |
| Four 1/16 mm offsets | 4 | 3 certificate regions | About 3.819635, 4.952060, 6.137439, 7.322819 nm | 1,655,989 |

## Earlier refusals resolved at the exit source

The earlier run recorded P1's base refusal at `slot-cut-6`, 19,999,999 units,
and P5's second disjoint cut at `cut_twice`, 20,000,000 units. Both now build
inside the same 20M allowance. The existing test follows its success branch
and checks all five stages, including the final cut, at 1, 2 and 4 threads.
No refusal expectation, tolerance or allowance changed for this rerun.

## Coverage boundary

The target-part chain preserves the fitted domains. It proves no growth for
disjoint cuts and rigid reflections, including the return to the original
placement. The existing certified-prism test covers an axial planar trim and
a later refit. Arbitrary cuts through fitted faces, new blend families, import
uncertainty, the release replay, area and STEP bounds, and the rest of the M12
exit run are separate checks. A refusal on one of those routes is not evidence of
zero error or a passing no-growth result.
