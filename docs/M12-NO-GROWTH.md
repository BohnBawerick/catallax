# M12 no-growth evidence

The no-growth claim covers the supported later-feature chains below. It does
not finish M12 or qualify cuts through fitted faces. Those cuts belong to M10
and the M12 exit run after M10.

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

These implementations remain unchanged. The new test makes their no-growth
claim explicit across longer chains and worker counts.

## Repeatable check

```sh
cargo test --release -p catallax-features --test no_growth -- --nocapture
```

`recorded_error_does_not_accumulate_at_one_two_or_four_workers` rebuilds the
canonical P1 bracket, P2 collar, P3 filleted boss, P4 shelled tray and P5 flanged
housing. P2-P4 then subtract a disjoint tool, reflect in an oblique plane,
reflect back, and subtract the tool again. P5 supports the chain through the
second reflection. Each build retains the 20,000,000-unit application allowance.

At each supported stage the test checks topology and entity counts. It
remeasures every fitted P5 carrier against its retained definition, checks its
recorded bound and the body's fit mark, and requires the same bound as before
the later features. P2-P4 have no fitted carriers in these chains. The test
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
| P2 collar | 5 | 0 | No fitted deviation | 100,733 |
| P3 filleted boss | 5 | 0 | No fitted deviation | 245,753 |
| P4 shelled tray | 5 | 0 | No fitted deviation | 3,433,733 |
| P5 through the second reflection | 4 | 4 | About 4.170328 nm on each face, unchanged | 12,341,919 |
| Direct 1/4 mm offset | 1 | 3 certificate regions | About 7.322819 nm | 1,601,665 |
| Four 1/16 mm offsets | 4 | 3 certificate regions | About 3.819635, 4.952060, 6.137439, 7.322819 nm | 3,708,172 |

## Recorded gaps

| Case | Typed refusal | Step | Work units | Owner |
| --- | --- | --- | ---: | --- |
| P1 base recipe | `numbers.operation-budget` | `slot-cut-6` | 19,999,999 | M2, `cx-m2-p1-allowance` |
| P5 after a disjoint cut and two reflections, then another cut | `boolean.parametric-query` | `cut_twice` | 20,000,000 | M10/M18 |

P1's refusal reports `operation work allowance exhausted`; no spline query
spends that allowance, so it keeps the plain budget code. P5's longer chain
spends it inside a fitted spline tangent query while imprinting, so the refusal
names the fitted spline edge with a trimmed face (the observed face kind is a
cylinder) and keeps the allowance exhaustion as its reason. The test compares
the full refusal, including that text, and the work at 1, 2 and 4 workers.
P1's extended recipe must refuse identically to the original recipe, before any
added feature. P5's longer chain is separate from the supported prefix checked
for no growth. Neither refusal counts as a no-growth pass. The claim extends to
those cases in a follow-up after their owning lanes land.
If the longer P5 chain builds, the test checks all five stages for no growth
instead of requiring the old refusal.

## Coverage boundary

The target-part chain preserves the fitted domains. It proves no growth for
disjoint cuts and rigid reflections, including the return to the original
placement. The existing certified-prism test covers an axial planar trim and
a later refit. Arbitrary cuts through fitted faces, new blend families, import
uncertainty, the release replay, area and STEP bounds, and the full M12 exit
run are separate checks. A refusal on one of those routes is not evidence of
zero error or a passing no-growth result.
