## 2026-11-07: whole tori from inscribed coaxial ring pairs

Pieces: shapes, features, recipe CLI
Change: `catallax_shapes` exports `RollingTorus::inscribed` and the
`InscribedTorus` type. The first returns the half tube a ball sweeps when
the radius is exactly the limit where two coaxial supports' offsets touch.
The B-rep fillet now closes both rings of a two-face shell of revolution
into a whole ring torus at that radius, where it refused
`fillet.radius-too-large`. Version 4 and 5 fillets fill
`Built::consumed` with the input faces they remove. `ConsumedFace` now
documents fillet records as well as draft records.
Why: `torus-intersection` is the plane z = 1 through the torus R = 6,
r = 2 with radius 1/2, exactly the inscribed radius of its meridian section.
RECIPE section 4.5.1 defines the result as the whole torus R = 6, r = 1/2.
OpenCASCADE refuses it. At 20M work the unchanged recipe exports one
raw-valid solid of volume 3 pi^2. Its input Boolean still needs about 1.94M
work, so the 1M snapshot is unchanged.
Migration: none for existing callers. `RollingTorus::new` keeps its strict
contract and refusals. Readers of `Built::consumed` see fillet records:
`by` holds the blend origins of the pair. Only complete pairs build: both
rings selected, equal convexity, one centre circle, opposite halves and
exact shared contact circles. Concave pairs round a lens-shaped void to a
toroidal void. Other configurations at the limit keep the refusal. The
independent field evaluator is unchanged and does not build this case.
