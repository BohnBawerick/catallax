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
Fresh runs at 1, 2 and 4 kernel threads must give the same records and work.
This probe reports acceptance and refusal separately; a passing determinism
assertion does not close a refused mass gate.

```sh
cargo test -p catallax-app --test assembly_queries
cargo test -p catallax-app --release --test m17_mass_gate -- --ignored --nocapture --test-threads=1
cargo clippy -p catallax-app --lib --test assembly_queries --test m17_mass_gate -- -D warnings
```

Raw evidence retained privately. The saved replay answers, measurement
tolerances, arithmetic policies and field integration allowances are unchanged.

## Remaining mass coverage

Authored P5 keeps the named `fillet.unsupported-overflow` field refusal.
[The radius-2 verdict](M18-P5-VERDICT.md) records the 1 mm rail/rim overlap,
two crossings per countersink mouth and missing trimmed continuation. A smaller
fillet, a B-rep copy or its untrimmed torus is not an independent field answer.
The remaining target verdicts will be recorded from the full probe.
