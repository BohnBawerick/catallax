Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Explicit numerical STEP requests

The verification recipe source now reads `catallax-numerical-step/1` with a required positive rational `tolerance_mm`. It builds the retained recipe geometry, converts its finite definition, admits it numerically at that preset tolerance, and writes approximate native STEP. The whole operation shares the caller's allowance. It does not retry exact requests or convert a numerical body into an exact model.

The separate `catallax-numerical-manifest/1` record retains source/emitted samples and repairs, file identity/checks, separate rounding scopes and measured work per stage. Success is answered with numerical proof, approximate category and `exact: false`. No whole-body bound is established. Refusals remove stale output, including work exhaustion after a successful write. The ordinary exact request adapter is unchanged.

The near-contact loft and recorded M10-008 use the native curved writer introduced in [its interface note](2026-11-15-numerical-curved-step.md). This adapter reuses its checks and the numerical evidence serializer; it introduces no new intersection or export algorithm. It follows the same OCCT native curve/surface representation and P04 consumer boundary. Multi-chart and fitted carrier admission remain later slice work.

`tools/step-regressions/test_numerical_step.py` drives both named cases with one, two and four workers and replays a portable question. `crates/catallax-cli/tests/numerical_step.rs` also checks malformed tolerances, unknown controls, wrong source, zero-work refusal, stale-file removal, unchanged exact refusal and a successful request after refusal. See [the wire contract](../../VERIFY-INTERFACE.md#numerical-step-output) for field meanings.
