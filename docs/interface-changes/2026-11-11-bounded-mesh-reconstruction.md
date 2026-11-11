# M23 first cut: exact-ring box/bore reconstruction

`catallax_field::reconstruction` is an opt-in Rust operation. It genuinely reads
source bytes through M22 `ClosedMesh::read_with_limit`, retaining the source
SHA-256 and exact scale separately from the inferred definition. There is no
input recipe, feature tree, ground-truth dimension, or fixture lookup. Existing
import, proof and export routes do not invoke it. Baseline:
`439aa03061fef1078bc82688463fe7f627fde4a0`.

## Accepted family and inference

The first cut accepts one outward connected, closed, noncrossing, manifold mesh
shell representing an axis-aligned rectangular box with **one through-bore along
z**. The bore may be cylindrical or a truncated cone with two positive rational
end radii. It must be strictly interior to the four side planes. The input has
exact planar patches and exactly two matching polygonal bore rings, with each
bore strip triangulated as a complete quad. No noisy-scan acceptance is claimed.

The recognizer derives the six outer planes from the mesh's exact coordinate
extrema. This is a constrained exact plane fit, not general least-squares plane
fitting. Triangle-adjacency region growing must yield exactly one connected
region on each plane and one connected residual bore region. The bore's two
circles are solved algebraically from three vertices each and **every** remaining
ring vertex must satisfy the resulting exact circle equation. Equal radii select
a cylinder; unequal radii and the observed z levels define a cone's linear
radius law. Both rings must have the same inferred centre and matching normalized
rays. Exact rational fitting has no residual tolerance or snapping threshold.
Irrational radii, approximate circles and noisy planes explicitly refuse.

This is more than carrier-surface fitting. The recognizer checks:

- M22's exact closedness, orientations, vertex fans, noncrossing triangles and
  winding/material validity before fitting;
- all seven connected regions, with counts accounting for every source triangle;
- one simple convex cyclic bore boundary, matching upper/lower edges, two
  triangles covering every four-vertex strip, and no unmatched bore triangle;
- cap boundary edges belonging only to the outer rectangular perimeter or the
  recognized bore ring; no disconnected cap patch or additional component;
- positive outer dimensions and both complete circles strictly within the sides.

These restrictions retain the through-hole and its adjacency. They refuse
unsupported topology instead of silently repairing it. The inferred definition
is a box minus a cylinder/cone, built by the existing exact feature/Boolean
kernel and independently compiled to the field. Structural B-rep validation is
required before return. Each control produces **seven analytic faces** rather
than retaining 192 mesh triangles. Ordinary build or field refusals propagate.
There is no STEP or UI integration in this operation.

## Proof and fidelity boundary

The output fields are immutable through this API. The returned recipe, B-rep and
field define the fitted mathematical solid. Exact queries on that solid do not
certify the original polyhedral mesh, let alone an unknown scanned object.
For example, in the cylinder control `(1.98, 0.28, 3)` mm is material in the
source polygonal bore but void in the reconstructed circular bore; the regression
checks these different answers explicitly.

Every successful reconstruction includes **both unproved finite-sample deviation
numbers**, at fixed 4/8/128 sampling described below. If the diagnostic refuses,
the operation does not return a partial reconstruction. A separate diagnostic
can resample with explicit bounded controls. Source bytes must match the retained
hash; the retained scale converts file coordinates to millimetres.

Both distances compare **boundaries**, not distances to solid interiors:

- Scan to part: barycentric grid on every source triangle, with integer weights
  summing to `triangle_steps`. At 4 subdivisions this is 15 samples per triangle,
  including vertices and edge midpoints: 2880 samples for each control. Distance
  to the analytic boundary is computed in binary64 as the minimum over the four
  bounded side rectangles, two rectangular caps minus their circular openings,
  and the bounded cylinder/cone. Cone distance uses projection onto its radial/z
  line segment, not an infinite supporting cone.
- Part to scan: each outer face uses an inclusive uniform parameter grid with
  `surface_steps` subdivisions; cap samples inside the hole are excluded. Bore
  samples use uniform angle and axial fraction, including both end rings.
  Distance to the entire source triangle surface uses M22's unproved fast
  nearest-triangle distance, not vertex distance.

Default `surface_steps=8`, `angles=128` gives 1620 cylinder and 1608 cone boundary
samples. Samples at shared edges may repeat. Explicit resampling ranges are
triangle subdivisions 1–16, surface subdivisions 2–64, angular samples 8–512.
These are sampling controls, **not geometric accuracy tolerances**. Nonfinite
samples refuse. No finite sample maximum is a global Hausdorff bound; unseen
extrema, binary64 rounding and between-sample error are not enclosed.

The operation caps source triangles at 4096 without raising M22's existing
limit. M22 and exact construction retain the caller's Decider. Region visits and
strip scans charge work; each diagnostic part sample conservatively charges the
whole triangle count because the fast nearest-triangle method has no Decider.
No per-stage allowance reset or fallback to approximate membership occurs.

## Public controls and results

The original [procedural generator](../../crates/catallax-field/examples/reconstruction/fixture.rs)
and its [cylinder](../../crates/catallax-field/examples/reconstruction/box-hole.obj)
and [cone](../../crates/catallax-field/examples/reconstruction/cone-hole.obj)
OBJ outputs are contributed under this repository's MIT licence. No external
asset or proprietary implementation is used. Each contains 96 vertices and 192
triangles. The generator's analytic ground truth is a 10 × 10 × 6 mm plate;
the cylindrical bore has radius 2 mm, the conical bore radii 2 and 3 mm. Exact
terminating-decimal ring coordinates come from rational rotations. Ground truth
and fixture generation are separate from the production recognizer.

Observed **unproved sampled maxima**, in mm:

| Control | Scan → part | Part → scan |
| --- | ---: | ---: |
| Cylindrical bore | 0.0201010126776671 | 0.0201003703140566 |
| Conical bore | 0.0297412751521858 | 0.0285088506565897 |
| Cylinder source versus deliberately filled box | 3.0000000000000000 | 1.9798989873223327 |
| Cone source versus deliberately filled box | 3.0000000000000000 | 2.9302325581395343 |

The filled-envelope diagnostic is deliberately adverse: it removes the hole
from the candidate boundary, supplies only diagnostic numbers, and cannot replace
the reconstructed part. The differing directions catch missing-feature coverage
that a single reported number would obscure. It samples 486 box-boundary points.

B-rep volumes agree numerically with `600 - 24*pi` and `600 - 38*pi` mm³:
524.601776313845 and 480.61947916358787. This is a floating analytic regression
comparison; the exact B-rep measurement retains its own enclosure. Source mesh
volume is intentionally different. Tests also check B-rep/field classification
of hole, material and exterior, translated/scaled inference, missing-triangle
refusal, an unmatched ring vertex, changed source hash and an extra disconnected
component. The latter is refused, never dropped.

Independent calls on exactly 1, 2 and 4 workers compare inferred recipe text,
both sample-max bit patterns, sample counts and work. The algorithm is serial
and uses deterministic traversal/reductions; this is same-platform evidence,
not a cross-platform transcendental reproducibility claim.

## Reproduce and remaining work

```sh
cargo test --release -p catallax-field --test reconstruction -- --nocapture
cargo run --release -p catallax-field --example reconstruct_mesh
cargo run --release -p catallax-field --example reconstruct_mesh -- crates/catallax-field/examples/reconstruction/box-hole.obj
cargo run --release -p catallax-field --example reconstruct_mesh -- --write-fixture box-hole.obj
cargo run --release -p catallax-field --example reconstruct_mesh -- --write-cone-fixture cone-hole.obj
cargo clippy -p catallax-field --all-targets -- -D warnings
cargo fmt --all --check
```

M23 is **not complete**: no 30,000 noisy-point P3 round trip, open-scan inference,
free-form fitting, general axis discovery, arbitrary triangulated multi-ring
bores, several holes or geometric fidelity certification. The six focused tests
and controls support only the admitted family. Classic exact constrained fitting
suffices for these controls; they provide no evidence that learned segmentation
is necessary or sufficient for general scans. No learned model or weights were
used. Deciding that for noisy P3 remains open.

No existing tolerances, saved answers, allowances or production routes change.
The only shared wiring is the new module export; there is no dependency on M17,
M24 or M9 changes. Full workspace, private replay and C++ gates are not run;
existing hosted CI must supply final-head integration evidence. No overall PASS
or general reconstruction-fidelity claim follows from these controls.
