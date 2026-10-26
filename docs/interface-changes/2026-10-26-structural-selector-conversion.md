## 2026-10-26: structural selector conversion in the SolidDraft converter

Pieces: application converter

Change: `soliddraft::convert` now converts `fillet` and `chamfer` with `all` edges and `shell` with an `all` or `created_by` open-face selector into the native `fillet`, `chamfer` and `shell` operations instead of returning `features.unsupported-operation`. `nearest`, `within`, `all_of`, and `all` open faces on a body built through a boolean or shell (whose faces a mandatory reference could lose) still return that `NotYet` code. A `created_by` open face whose step makes no faces also returns it.

Why: these selectors have an exact structural mapping to edge and face references, so no built body is needed to resolve them.

Migration: callers that expected `NotYet` for these recipes now receive a converted recipe and, through `build_soliddraft`, a built body or the native kernel's own refusal. Callers matching `features.unsupported-operation` for other selectors are unaffected. No recipe format, budget, tolerance or STEP option changes.
