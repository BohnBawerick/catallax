# Product scope

Catallax is the MIT-licensed geometry kernel for SolidDraft. It owns editable part geometry and the answers computed from it. SolidDraft owns the feature workflow, drawings, documents and user interface. See the [vision](../VISION.md) for the product principles, the [architecture](ARCHITECTURE.md) for implementation rules and the [roadmap](../ROADMAP.md) for scope and order of work.

## What a user must be able to do

- Build a part, edit an earlier feature, reopen it and deliver the result. Preserve references by construction name across an edit. A reference that cannot be resolved must fail explicitly, not attach to a different face.
- Open a bought-in part, see what imported, and edit an admitted solid. Keep a retained surface preview distinct from an editable solid. When reconstructing a scan, report deviation in both directions. A scan result need not recover a feature tree.
- Measure and export the requested geometry. When a result is numerical, report its actual tolerance, checks and repairs. Label approximate STEP output. An unproved preview must never look like a delivered result.
- Print parts at a stated accuracy. Let the caller choose a different output accuracy explicitly. For supported field-only parts, a mass stand-in is separate from a print mesh and must report which properties it preserves.

## How to tell whether it works

Use [real edit chains](TARGET-PARTS.md) and repeatable public census cases. Compare like with like, including the same file limits and a distinct count for checked exports. A library-only example, preview or faster refusal does not complete the app workflow. See [capabilities](CAPABILITIES.md) for released evidence and known limits.

The kernel keeps exact answers deterministic and does not link another kernel as a fallback. Use reference kernels as comparisons, not as hidden substitutes. A possible input that still refuses is a gap to fix in the [owning roadmap step](../ROADMAP.md).
