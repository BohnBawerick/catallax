# N02: numerical hole trims and explicit endpoint sewing

Piece: topology, opt-in numerical admission only. Existing exact consumers,
`Answer`, `BoundedEvidence`, proof/export boundaries and saved answers keep their
contracts. No caller migration is required.

`catallax_topology::numerical::build` admits one simple outer linear u/v loop
and strictly interior simple holes per face. The outer orientation follows the
face sense; holes have the opposite orientation. Hole/outer touching, crossing,
outside holes, nested holes, self-crossing, open or reversed loops refuse.
Numerical admission also rejects shared walked endpoint IDs across distinct
boundaries after endpoint remapping, independently of their unchanged UV trims.
The completed `trim.boundary-endpoints-disjoint` record names this local check;
no stronger guarantee is inferred from the generic topology checker.
Every loop retains the existing pointwise distance checks and carrier-specific
admission rules. `trim.hole-contained-disjoint` reports the finite linear
parameter-space check; it is not a whole-body geometric certificate.

Planar complete-circle boundaries may be outer rims or holes. Projected disks
enclose circular holes with the stated tolerance added to each radius.
An outer circle's projected ellipse contains a disk of radius
`r * abs(cos(angle between normals))`; the stated tolerance is subtracted from
that interior radius. Hole disks must fit inside it or remain disjoint from
other hole disks. Finite-line neighbours receive complete segment clearance
or interior endpoint checks. Unsupported curved neighbours refuse.
Simple closed-loop checks bind every segment end to the next start.
Circular pairs use these continuous checks rather than intersecting inscribed
chord polygons, which can cross for a valid thin annulus with different rim phases.
The mixed record is `trim.sampled-boundaries-with-continuous-circle-clearance`,
with separate `trim.continuous-circle-hole-clearance` and
`trim.continuous-circle-outer-clearance` entries where applicable.
Other projected/native curved walks retain their sampled scope.
The continuous-circle records are emitted only for the planar family that
actually ran projected-disk clearance. A rational patch with native UV circles
retains sampled-only evidence even when its shared 3D edges are complete circles.
No tolerance, default work allowance or exact-proof authority is increased.

New API:

```rust,ignore
numerical::sew(
    input: &NumericalDefinition,
    pairs: &[[usize; 2]], // [keeper, removed] in the original vertex table
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal>
```

This bounded family unifies caller-selected disjoint endpoint pairs on an already
shared edge table. It does not match or merge independent edge carriers. Each
removed vertex is relocated to its keeper only if its 128-bit outward movement
enclosure is within the explicit local tolerance. Edge endpoint references are
remapped and the removed table entries are compacted. Carriers and finite trims
remain unchanged. The entire new shell is rebuilt and checked before return.
The borrowed input is unchanged on success or refusal. No search, tolerance
increase or iterative repair occurs.

Each nonzero movement retains N01's existing `NumericalRepair::VertexMoved`, with
the removed key, before/after positions and measured upper distance. That key is
historical; the keeper key identifies the retained vertex. A completed
`sewing.explicit-vertex-pair` check records the removed-to-keeper mapping, including
zero-gap joins. Zero movement adds no `VertexMoved` record. Later N01 continuations
retain movement history and recheck the resulting definition; operation-specific
sewing check records are not themselves carried forward as new admission checks.

The return value is exactly N01's separate `NumericalBody`/`NumericalEvidence`,
including when all residuals are zero. There is no new numerical wire schema,
implicit conversion, exact incidence fact, bounded result or STEP authority.
N01's generic-constructor trust caveat still applies. Existing structural,
`numerical.invalid-definition`, `numerical.invalid-trim`,
`numerical.repair-tolerance-exceeded`, sample, precision and budget refusal codes
are reused. Sewing and rebuilding share the active request allowance, or N01's
unchanged 20M default when none is active.

Executable accept/refuse examples and public evidence are in
[the N02 report](../phase1/N02/REPORT.md).
