# Typed spline-face offset limit

A nonzero body offset that retains a B-spline face requires M10's certified
metric offset and finite trim arrangement. The ordinary B-rep offset builder now
returns `distance.spline-offset-unsupported` with that cause at the offset step.
The represented spline itself remains valid geometry.

Zero distance still returns the input body by identity before the guard. Empty
input also retains its identity rule. The new check reads carrier kinds without
asking a numeric question. Other analytic carrier kinds continue through their
existing dispatch branches.

The guard belongs to `offset`. Existing cap-only spline shell operations retain
their proved planar normal construction and are covered by their regression tests.

The public fixture is `crates/catallax-features/tests/fixtures/spline-face-offset.json`.
Its stock is an exact extrusion of a cubic Bezier profile. Both offset signs have
a typed M10 refusal. A separate test checks zero-offset body and face identity.

The complete CLI's zero-offset control retains its existing
`step.unsupported-geometry` result because that represented spline has no specified
STEP form. The kernel body identity test succeeds before export. This exporter
limit remains distinct from the nonzero offset refusal.

All verification uses the existing 1M work allowance. The saved 477-case inputs
contain no spline-offset request; their unchanged hosted replay remains required.
Raw evidence retained privately.

## Measured checks

[Portable evidence](evidence/m4-spline-offset-refusal.json) records the exact
requests, source hashes and fresh outcomes. Positive and negative offsets each
use 47 work units before and after the guard. Their reason changes from the
generic envelope code to `distance.spline-offset-unsupported`.

Zero-offset CLI control uses 450 units and retains its existing export refusal.
The cap-only shell control retains its 1M budget refusal at 999,998 units; its
B-rep construction remains covered by the six passing spline-cap-shell tests.
Both new kernel tests pass. Clippy passes with warnings denied.

All four requests agree at 1, 2, 4 and all 12 available kernel threads with
1, 2, 4 and 4 persistent workers. Work, manifests, operation logs without clocks
and absence of STEP output agree. No tolerance or work allowance changed.
