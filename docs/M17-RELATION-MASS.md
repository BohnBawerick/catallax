# M17 retained-body relation and mass agreement

## Existing work reused

The exact overlap-volume, clipped contact-face and regularised empty-body
queries landed in [the application-query change](https://github.com/BohnBawerick/catallax/pull/537).
The field's bounded assembly search landed in
[the relation-search change](https://github.com/BohnBawerick/catallax/pull/506).
This change composes those answers with existing trimmed-face distances. It
does not add a Boolean engine or a distance approximation.

## B-rep relation

`catallax_app::queries::body_relation` takes two retained body handles in one
model. Their geometry already includes their placements. It returns:

- `Clash` only when the overlap-volume enclosure has a positive lower bound.
- `Clear` only after proving zero overlap, then enclosing the minimum distance
  over every pair of closed boundary faces. Its lower bound must be positive
  and its width must meet the requested `Accuracy`.
- `QueryRefusal` for contact, empty operands, unsupported distances, uncertain
  overlap, a wider gap enclosure or exhausted arithmetic work.

Zero overlap establishes the noncrossing precondition of the existing
trimmed-entity distance query. A minimum over every boundary pair covers all
lumps and cavity walls. Bounding-box separation alone supplies no answer.
Face, edge and point contact all refuse `assembly.contact` with a proven zero
gap. A positive gap never becomes contact through a tolerance threshold.

Nonzero distance or angle construction marks refuse
`app.relation-construction-bound`, reporting the mark's code, kind and numeric
bound. Distance on a fitted carrier alone cannot establish the relation of its
original definition. Other query refusals retain their underlying code,
undecided number and both operands' topology counts.

The query leaves its source model unchanged. It shares the caller's active
work allowance or uses the unchanged 20,000,000-unit application default.
`ApplicationBody` exposes the matching method.

## One reported number

`catallax_app::assembly::relation_agreement` reconciles B-rep and field relation
answers for the same placed definitions. Two clashes agree on positive-volume
interference. Two clear answers must have overlapping gap enclosures, and the
width of their hull must meet the requested accuracy. The result carries that
single hull. A disagreement, wider hull or field refusal remains unresolved.

`assembly::uniform_mass` measures B-rep volume and centroid without requesting
inertia. Density is an exact positive value in grams per cubic millimetre.
Volume and centroid retain their separate requested accuracies. A wider
enclosure refuses with the number reached. A field-only region envelope refuses
`assembly.mass-envelope`; it is not the lattice's material volume.

`assembly::mass_agreement` takes independent mass answers for the same recipe,
placement and density. It applies the same hull rule to grams and each world
centroid coordinate. Absolute gram accuracy is density times volume accuracy;
relative accuracy is unchanged. It retains B-rep construction marks. Zero mass
has no centre. A missing centre on either positive-mass answer refuses.

The hull encloses both measured values. Intersecting two intervals alone cannot
prove equal geometry, and the API does not claim that it does. Callers retain
both query outcomes alongside the reported value. No refused route silently
becomes independent agreement.

`assembly::resolve_mass` reports one checked result with `MassSource::Agreed`,
`Brep` or `Field`. A single-source value retains the other query's original
refusal in `ResolvedMass::unavailable`. Two proven but incompatible results
refuse; neither is selected. Both unavailable routes remain a numeric typed
refusal. This selection performs no build and creates no work allowance.

## Reproducible checks

`assembly_queries` exercises the public queries through complete recipe
evaluation and independent field compilation. Controls cover overlap,
containment, disjoint spheres, rational and irrational gaps, face/edge/point
contact, clearance inside a cavity, a positive 10^-12 mm gap, empty operands,
exhausted work, invalid accuracy, density units, one-number agreement,
disagreement, wide hulls and missing centres.

The explicit release-mode `m17_mass_gate` probe evaluates the unchanged public
P1-P6 recipes with default mass accuracies, the default 2,000,000-cell field
limit and the 20,000,000-unit application allowance. It records independent
construction and measurement work, enclosures or typed refusals, and agreement.
Fresh runs at 1, 2 and 4 kernel threads must give the same records, work and
exhaustion flags. Every successful construction or measurement must remain
non-exhausted; a failed route may legitimately exhaust its allowance.
This probe reports acceptance and refusal separately; a passing determinism
assertion does not close a refused mass gate.

```sh
cargo test -p catallax-app --test assembly_queries
cargo test -p catallax-app --release --test m17_mass_gate -- --ignored --nocapture --test-threads=1
cargo clippy -p catallax-app --lib --test assembly_queries --test m17_mass_gate -- -D warnings
```

Raw evidence retained privately. The saved replay answers, measurement
tolerances, arithmetic policies and field integration allowances are unchanged.

## Target mass verdicts

The complete probe passed at source
`b390d5f75f2d282f2e55364667495d035b03cfb0`, with identical exact enclosures,
selected sources, typed refusals, undecided numeric records and operation work
at 1, 2 and 4 workers. Density is the explicit test input `0.00785` g/mm3.
The table reports
outward-rounded bounds; their displayed width is not the finer width used to
check default accuracy. Each answer carries its source. A single-source answer
retains the other kernel's typed refusal.

| Part | Proven mass and centre, if available | Remaining refusal |
| --- | --- | --- |
| P1 | Field mass `[156.5599131409369,156.5599131409370]` g; centre near `(30,14.201472877614,8.271602017314)` mm. | B-rep `numbers.operation-budget` at `slot-cut-6`, 19,999,999 construction work. The construction cost belongs to M2. |
| P2 | Agreed mass `[31.38587196617859,31.38587196617860]` g; centre `(0,6,0)` mm. | None on the opt-in meridian query. The original field integrator retains `field.print-coverage`. |
| P3 | B-rep mass `[113.2648249,113.2648250]` g; centre x=20, y=15, z=`[6.846759374,6.846759375]` mm. | Independent field mass refuses `field.print-coverage`, 3,628 work including compilation, placement and the opt-in coverage probes. |
| P4 | Field mass `[37.92534104783266,37.92534104783267]` g; centre x/y enclose zero and z=`[4.062614653686839,4.062614653686840]` mm. B-rep independently proves the same mass and centre. | The new prismatic route closes the field integration gap. The original integrator retains its default-cap `field.mass-budget` refusal. |
| P5 | Both routes refuse at default accuracy. | B-rep `assembly.mass-accuracy`: volume enclosure `[12558.68134,12558.68205]` mm3 versus allowed width `1.255868135e-8`; measurement work 6,204,328. Field `fillet.unsupported-overflow`, 3,367 work including compilation. |
| P6 | Independent field mass `[3.025218110,3.025218111]` g; centre x=`[5.980637620,5.980637621]`, y=`[3.060448818,3.060448819]`, z=`[3.076344223,3.076344224]` mm. Field route work 2,438,426 including compilation, placement and the opt-in coverage probes. | B-rep `assembly.mass-envelope` correctly excludes its single field-region envelope. A missing B-rep authority does not invalidate the proven field mass. |

P1's new opt-in prismatic CSG query closes its field mass without waiting for
the M2 B-rep construction cost. At the same density, default accuracy and
allowance it returns `[156.5599131409369,156.5599131409370]` g and centre
x=`[29.99999999999999,30.00000000000001]`,
y=`[14.20147287761387,14.20147287761388]`,
z=`[8.271602017314395,8.271602017314396]` mm. Query work is 2,795,348 and exact
enclosures/work agree at 1, 2 and 4 threads. The completed recipe is unchanged.
[The CSG moment record](interface-changes/2026-11-13-prismatic-field-mass.md)
states the proof and independent component-formula check. The original field
mass route retains its `field.print-coverage` refusal.

P2 now has an independent exact field route for Boolean combinations of coaxial
polygonal revolutions. Its mass and centre agree with B-rep at default accuracy,
with identical exact enclosures and 389 measurement work at 1, 2 and 4 threads.
[The meridian integration record](interface-changes/2026-11-13-coaxial-field-mass.md)
states the proof and coverage. The full probe now attempts this route explicitly
before the original field integrator; it also tries the new prismatic route
for its covered definitions. Its exact thread comparison retains unrounded
results and undecided bounds, separate from the printed decimal reports.

P4's complete shell definition also qualifies for the prismatic query. Every
halfspace must prove an axis-box definition; a finite enclosing box alone never
qualifies. A separate circular-segment formula checks the clipped outer rectangle
and the complete inner cavity. At 1, 2 and 4 threads, exact enclosures and work
agree. Compilation, placement and mass use 40,411 work units.

P1 and P6 have field-only proven values. P2's independent agreement is closed.
P4's independent mass and centre agree with B-rep at default accuracy. P3 has a B-rep-only
proven value and a visible field refusal. The full gate
requires those five numeric answers and authored P5's named refusal, with the
same selected sources, exact enclosures, refusals and operation work at 1, 2 and
4 threads. It claims independent agreement only where both queries succeeded.
Independent field integration on P3 remains an open coverage item. No build
requirement is parked by this record.

## Authored P5 refusal

Authored P5 keeps the named `fillet.unsupported-overflow` field refusal.
[The radius-2 verdict](M18-P5-VERDICT.md) records the 1 mm rail/rim overlap,
two crossings per countersink mouth and missing trimmed continuation. A smaller
fillet, a B-rep copy or its untrimmed torus is not an independent field answer.
