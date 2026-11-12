## 2026-11-12: root-list consumers refuse instead of panicking

Pieces: 1, 4
Change: `RealPoly::real_roots` and the rational curve-surface solve in catallax-intersect (`rur_points`) refuse when `Decider::real_roots` answers with anything other than a root list or a refusal: an undecided sign keeps its `numbers.*` code and bound, and a coincidence refuses with the new code `numbers.root-isolation`. Before, both arms were `unreachable!` and would panic.
Why: a follow-up of [NOT-EXACT-INVENTORY.md](../NOT-EXACT-INVENTORY.md), section 2.4 (dead arms). `Decider::real_roots` answers only with a list or a refusal today, and a numbers test now holds it to that across allowances, so no answer changes; the arms are kept as refusals so a future change to root isolation cannot crash a request.
Migration: none. `numbers.root-isolation` is unreachable today.
