## 2026-11-07: explicit variable-fillet radius laws

Pieces: recipe language, fillet corpus and independent OCCT observer.

Change: `catallax-recipe/5` accepts `fillet_variable {body, edges}`. Each edge record stores an `EdgeRef`, a `from: VertexName` endpoint and `radius: RadiusLaw`. The initial tagged law is `{"kind":"linear","start":"1","end":"2"}`. Its parameter is normalized arc length of the original input edge, oriented from the named vertex. Either endpoint may be exactly zero; negative radii and an identically zero law fail the structural check. The public `RadiusLaw` enum is non-exhaustive so later law kinds do not need a different operation.

Why: the scalar `fillet` cannot express interpolation, orientation or a collapsed end. Approximating a law with ordered constant-radius steps would change its definition. The contact-section and continuation contract is in [RECIPE section 4.5.2](../RECIPE.md#452-variable-radius-fillets).

Migration: existing `fillet` operations and saved versions are unchanged. Producers of variable fillets must explicitly select version 5, store one endpoint name per selected open edge and write exact string radii. An edge reference must resolve uniquely; a ring has no `from` endpoint under this contract. Selected junctions require equal endpoint radii, checked by geometry rather than by list order. Empty lists and repeated identical references fail the structural check. Both `blend.edge` and `corner.vertex` accept origins from the new operation.

This is the separate schema change required before construction. Neither evaluator constructs the operation yet. The actual CLI returns its existing `features.unsupported-operation` refusal, not a constant-radius fallback. Both variable corpus rows now execute that path and retain their original independent OCCT endpoint laws. They count as construction refusals rather than synthetic schema gaps. No geometry, whole-face certificate, STEP, mesh or SolidDraft acceptance follows from this change.

SolidDraft producers should preserve the named orientation and tagged law through edits and persistence, and retain the typed geometry refusal until the construction change lands. Do not enable a working variable-fillet control solely because this schema parses.
