## 2026-11-07: whole tori from inscribed coaxial ring pairs

Pieces: shapes, features, recipe CLI
Change: `catallax_shapes` exports `RollingTorus::inscribed` and the
`InscribedTorus` type. The first returns the half tube a ball sweeps when
the radius is exactly the limit where two coaxial supports' offsets touch.
The B-rep fillet now closes both rings of a two-face shell of revolution
into a whole ring torus at that radius, where it refused
`fillet.radius-too-large`. The pair adds both consumed input faces to
`Built::consumed`, crediting its blend origins. `edge_blend::Replacement`
gains a set of consumed faces that leave their shell.
Why: `torus-intersection` is the plane z = 1 through the torus R = 6,
r = 2 with radius 1/2, exactly the inscribed radius of its meridian section.
RECIPE section 4.5.1 defines the result as the whole torus R = 6, r = 1/2.
OpenCASCADE refuses it. At 20M work the unchanged recipe exports one
raw-valid solid of volume 3 pi^2. Its input Boolean still needs about 1.94M
work, so the 1M snapshot is unchanged.
Migration: none for existing callers. `RollingTorus::new` keeps its strict
contract and refusals. In the pair's `Built::consumed` records, `by` holds
the blend origins. Only complete pairs build: both
rings selected, equal convexity, one centre circle, opposite halves and
exact shared contact circles. Concave pairs round a lens-shaped void to a
toroidal void. A convex shell must be the only shell of its lump and a
concave shell a void of the body's only lump; otherwise the pair refuses
`fillet.unsupported-overflow`. Other configurations at the limit keep the
radius refusal. The
independent field evaluator is unchanged and does not build this case.
