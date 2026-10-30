## 2026-10-30: finite cone bands exclude a tube inside their smaller ring

Pieces: 4, 5
Change: The Boolean fallback now excludes a certified finite cone band and a parallel cylinder or cone band when the second tube lies strictly inside the first band's smaller ring. This builds a plain hole inside a bevel ring, such as the P2 collar's annulus between radius 6 and the bevels at radius 11. A band that reaches its apex has a smaller ring of radius zero, so it never takes this exclusion.
Why: The hole's support meets the cone's unbounded support beyond the finite band, and the band's axis-aligned box contains the hole, so the box and external tube proofs could not separate them. No trimmed curve exists, so the refusal was a missed exclusion and not a missing trace.
Migration: None. The refusal code is unchanged for every pair that is not proven apart. Exact contact, crossing and skew pairs still refuse `intersect.traced-pair-required`. This does not add cone/cylinder tracing or change STEP output code. The SolidDraft collar scenario and pin remain a separate companion change. See the prerequisite [countersink placement](2026-10-30-countersink-placement.md) change.
