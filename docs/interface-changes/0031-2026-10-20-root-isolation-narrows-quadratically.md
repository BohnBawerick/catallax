## 2026-10-20: root isolation narrows quadratically

Pieces: 1 (and every piece that calls `real_roots`)
Change: `real_roots` narrows isolating intervals, including the narrowing that decides whether a root is rational, by quadratically convergent steps sized relative to the root, kept only when exact signs confirm the bracket, with bisection only as a fallback. Piece 1's case 40 gains a quartic with a 2,800-bit leading coefficient that must be isolated under the default budget. No signature, outcome, code or budget default changes.
Why: on main, when the rational candidates cannot be enumerated, each root is bisected down to width `1 / (4 q^2)`, with `q` the leading integer coefficient. For piece 5 case 23 at scale `10^-30`, a probe-plane quartic with `q` of 1,428 bits needed about 2,760 bisections per root, exhausted the 10^9 work budget and refused the boolean. A scale-relative safeguarded Newton narrowing isolated the same quartic, rationality check included, with 1.3e7 work, and the boolean answered with the expected counts.
Migration: piece 1 replaces bisection-only narrowing in `algebraic::isolate`. Callers need no change.

