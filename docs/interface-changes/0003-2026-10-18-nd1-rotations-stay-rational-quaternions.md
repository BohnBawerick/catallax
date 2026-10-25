## 2026-10-18: ND1 rotations stay rational quaternions

Pieces: judge interface, cli
Change: a pose rotation is a quaternion `[w, x, y, z]` of exact numbers; integer quaternions are the usual form. Rotation matrices are not accepted.
Why: every rotation matrix with rational entries comes from an integer quaternion, so one form reaches every exact pose and the kernel never has to check a matrix for orthogonality.
Migration: none.

