## 2026-10-18: ND11 seams and poles in the operand's frame

Pieces: 2, 3, 7, judge interface
Change: every periodic surface of an operand has its seam on the local half-plane `y = 0, x > 0`, a sphere's poles lie on the local `z` axis, and a torus's tube seam is its outer equator. Booleans keep operand surfaces, so answers keep those seams, and STEP placements show them in the same place. This confirms the judge's provisional behaviour.
Why: this was already the builders' and piece 2's parametrisation; it is now stated in the judge contract so "along the seam" poses can aim at it.
Migration: none; piece 3's builders already use `ref_dir = +x` and `axis = +z` in the local frame.

