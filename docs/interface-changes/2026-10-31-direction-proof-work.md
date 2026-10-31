## 2026-10-31: STEP direction proofs without the radical

Pieces: step
Change: `write_step` proves that a written unit-direction component equals its decimal through `v^2 = q^2 n` instead of comparing `v / sqrt(n)` with `q`. The 3- and 5-degree planar draft recipes now answer with exact analytic STEP inside a shared 20M request, at about 4.2M and 2.0M work, where they refused with `numbers.operation-budget` at both 20M and 100M.
Why: proving an edge direction written as `1.` or `-1.` through the radical cost 4M to 14M work each on the drafted planes' nested angle towers. The squared form costs about 5,000 and yields the same exact result, so the writer's certified errors do not change.
Migration: a caller that matched the budget refusal for these drafts now sees an exact manifest/3 record with a STEP file. No format, allowance, significant-digit or tolerance changes. Recipes that already answered keep their files and bounds.
