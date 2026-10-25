## 2026-10-24: direct draft face proofs

Pieces: field

Change: named-face queries on a sole result leaf no longer compute the field's full distance. Exact leaf membership and named trim already prove result ownership. A point strictly inside a published polyhedral draft boundary tile also uses its stored outward plane normal directly. Points on tile edges retain the existing ray proof. Boolean results retain their global ownership and coactive-support checks.

Why: comparison asks for named boundary membership and direction at many points. Repeating trimmed-boundary distance and whole-cell ray calculations for a certified facet interior does not strengthen either answer.

Migration: no public signature, answer format, tolerance or work-allowance changes.

