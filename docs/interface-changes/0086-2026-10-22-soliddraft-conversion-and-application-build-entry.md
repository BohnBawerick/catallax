## 2026-10-22: SolidDraft conversion and application build entry

Pieces: recipe, features, mesh, application
Change: new `catallax-app` crate exposes `soliddraft::convert`, `build` and `build_soliddraft`. Results carry the exact `Built`, marked display mesh with face and edge picking handles, complete face-origin sets and optional persistent edge/vertex names. SolidDraft conversion returns the original-tag to canonical-step table. Typed errors separate invalid source, not-yet conversion or operations, step-local build refusals, topology faults and display refusals; `Error::code()` returns the stable code. No legacy CLI or exchange format changes.
Why: SolidDraft needs one entry from recipes to named, displayable bodies without inventing topology identities or taking mesh geometry as the model. Fillet, chamfer, shell and unresolved conversion constructions remain explicit not-yet results.
Migration: SolidDraft's backend consumes this entry in its own repository. Map returned handles explicitly to backend ids, preserve merged origin sets, and persist names rather than handle indices. `docs/RECIPE.md` section 6.3 describes the output contract and current conversion limits. Recipe booleans also apply exact repeated-operand set identities before invoking the distinct-operand intersection graph, so a body combined with itself no longer panics.

