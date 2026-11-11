# Bounded field assembly queries (M17 first cut)

`catallax_field::assembly` adds an opt-in Rust API. Existing recipe, field,
measurement, proof and export routes do not call it. M17 is not complete.

## Placement and material

`PlacedPart::new` applies an exact checked rigid `Pose` after all recipe
placements. Queries use world millimetres. Scaling is not a placement.
`uniform_mass` requires a caller-supplied, positive exact uniform density in
**g/mm³**. It returns a mass enclosure in grams and centre-of-gravity enclosures
in world mm. The caller must supply the material model; there is no implicit
steel, per-region density, inferred material or gravity-field model.

The existing field `mass::properties` integrator supplies volume and centroid.
Its `Accuracy` still refers to unit-density volume, centroid and inertia, and
its cell cap and numeric policy are unchanged. Mass width is density times
volume width. Centroid is absent for exactly zero volume. Integration refusals
(including inability to certify the requested inertia) are preserved. No print
mesh, preview volume or STEP stand-in is substituted. This first cut does not
reconcile independent B-rep measurements into one user-facing value.

## Clash and global smallest gap

`relation` returns one of:

- `Clash`: a rational world point certified strictly inside both complete
  solids. By interior openness this establishes positive-volume overlap.
- `Clear`: a strictly positive enclosure of the global smallest distance
  between the closed solids, with width no greater than requested `gap_width`.
- `Unresolved(Refusal)`: a typed unsupported-definition, contact, empty-input,
  numeric-budget or pair-budget outcome. Touch is not positive-volume clash
  and never becomes clear by a tolerance rule.

The root boxes are conservative bounds of the complete field definitions.
Top-level unions retain separate component bounds instead of filling their
voids with one outer box; the initial Cartesian cover must fit the pair cap.
Each pending pair covers a subset of the Cartesian product of the solids.
The rational squared distance between its axis-aligned boxes is a lower bound
on every point pair it covers. Subdivision retains **both** children, and a
pair is removed only when the existing interval classifier proves one entire
closed box strictly outside its solid. Unknown cells remain in the cover.
The minimum of all pending lower bounds is consequently global.

Closest box points and box centres are only proposals. The existing proven
point query must certify membership in the actual closed solid before any
proposal supplies a distance upper bound. A coincident proposal proves contact
but does not by itself prove absence of overlap elsewhere. Such a case returns
unresolved. No difference of field values, nearest-owner vector, or nearest
witness alone certifies clearance. Square roots have outward rational bounds;
rational roots remain exact. The fixed 192-bit root enclosure can itself limit
very tight requests; these remain unresolved rather than widening accuracy.

`Options` has no silent default: callers choose a positive gap width and a
positive maximum number of processed box pairs. Each iteration charges work
and checks the existing numeric policy. The first cut uses deterministic serial
subdivision; independent calls may run in any number of workers. It neither
adds a thread pool nor consumes the unproved preview batch API.

## Validation and limits

Analytic controls cover axis boxes, spheres, exact diagonal distance, irrational
distance enclosure, tiny positive gaps, overlap, contact, invalid density and
options, and overlapping outer envelopes with disjoint material. A common
90-degree rotation plus translation checks clearance and centre of gravity.
The tests run complete public P1–P6 recipes with unchanged default property
accuracy, a 128-cell integration cap and a 16-pair assembly cap. Their independent
runs at 1, 2 and 4 workers must produce identical outcomes. These small caps are
explicit test controls for this new API, not changes to existing allowances.

The algorithm is a bounded proof search, not a completeness claim. Curved and
nonconvex closest pairs may exhaust subdivision. Existing interval or mass
coverage may refuse field features. Empty operands have no finite smallest
gap. There is no assembly recipe syntax, CLI/UI integration, heterogeneous
material model, contact classification, penetration depth or B-rep reconciliation
in this first cut. No other lane's proposed APIs are required.

Observed complete public-target results at those controls (second copy translated
100 mm in world x; density 0.00785 g/mm³ is an explicit test input):

| Target | Mass and centre of gravity | Relation to translated copy |
| --- | --- | --- |
| P1 bracket | `field.print-coverage` | `assembly.pair-budget` |
| P2 collar | `field.print-coverage` | `field.print-coverage` |
| P3 filleted boss | `field.print-coverage` | `field.print-coverage` |
| P4 shelled tray | `field.mass-budget` | clear, exact 60 mm |
| P5 flanged housing | compilation: `fillet.unsupported-overflow` | same compilation refusal |
| P6 lattice bracket | `field.mass-budget` | clear, exact 88 mm |

Thus **no complete target's mass/centre-of-gravity gate is closed** by this
change. The bounded mass API is established on analytic controls; P4/P6's finite
cell run does not establish impossibility at a larger caller-selected budget.
Existing mass/field algorithms, accuracy defaults and work allowances are
unchanged. The only shared helper change exposes the existing exact axis-box
bounds to this sibling module without altering its implementation.
