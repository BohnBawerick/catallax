## 2026-11-11: unproved processor preview evaluator

Pieces: catallax-field
Change: `catallax_field::Plan::unproved_preview` returns a new
`UnprovedValue` (`value`, `true_distance_formula`, `hardware`), where
`hardware` is the new non-exhaustive `PreviewHardware` with the single
variant `Processor`. It evaluates the same tape and leaf formulas as
`Plan::fast` with the processor's square root for each point, and refuses the
same inputs with the same codes. `Plan::fast`, `Approximate` and `fast_batch`
keep their signatures and return the same bits as before.
Why: PLAN M14 F4, under the [unproved hardware preview decision of
2026-11-10](../ARCHITECTURE.md#unproved-hardware-preview-evaluator-decision-of-2026-11-10).
The deterministic software square root is most of `Plan::fast`'s cost. The
new value is a separate type with no conversion into `Approximate`, `Query`,
`Value`, `DistanceBound` or `FaceQuery`, and an integration test fails if any
workspace library source outside the plan module names it.
Migration: none for existing callers. A display caller that adopts the
preview adds its module to the allowed list in
`crates/catallax-field/tests/unproved_preview.rs` in its own reviewed change;
export, print, measurement, mass, distance and comparison modules stay on the
forbidden list.
