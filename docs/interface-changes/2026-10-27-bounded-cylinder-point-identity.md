## 2026-10-27: bounded certificates for cylinder-point identity

Pieces: 1, 5

Change: `Decider::certify_root_at` attempts a root-box certificate at one
working precision, within the current policy and request allowance. Unlike
`certify_root`, it does not escalate precision after an inconclusive attempt.
Boolean point canonicalisation uses metered coordinate bounds to separate
candidates, or a unique-root certificate to equate candidates with recorded
incidence on the same three rational cylinders. The first point and both
original coordinate witnesses are retained.

Why: independently constructed algebraic witnesses can denote the same point.
A nearly singular box can also contain two distinct points, so shared
incidence and overlapping coordinate bounds are not an equality proof.

Migration: existing root-certification callers are unchanged. An inconclusive
optional certificate falls back to the existing exact point comparison;
deadline and exhausted request allowances still propagate. Padding a candidate
box is not a geometric tolerance. This shortcut proves neither trim membership
nor orientation. Budgets, tolerances, file formats and refusal codes are
unchanged. Boundary regressions cover two roots separated by less than a
nanometre, incomplete incidence, closed-box roots, and caller accounting.
Full canonical P1 construction and SolidDraft U2 acceptance remain unproved.
