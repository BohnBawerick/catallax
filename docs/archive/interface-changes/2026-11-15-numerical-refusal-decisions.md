Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-15: retain numerical refusal decisions

Pieces: CLI verification.

Change: a refused `catallax-numerical-result/1` now includes optional
`refusal.undecided` when the kernel retained a `NotExact` decision. It preserves
the decision's `code`, `detail` and `bound`. The bound uses the exact rational
string representation of the numerical result format, including nonterminating
fractions such as `1/3`. The field is absent when no decision was retained.

Why: the initial [numerical build result](2026-11-15-checked-numerical-build-results.md)
serialized only the outer refusal code and detail. Distinct retained residual
bounds could therefore produce identical wire records.

Migration: numerical refusal consumers may read the optional decision. Its bound
is the magnitude of the stopped numerical quantity, whose value lies in
`[-bound, bound]`. It has no implied millimetre unit and grants no accepted tolerance or
complete bound on a body. Refused records still publish no geometry, numerical
acceptance evidence or STEP result. Existing exact manifests and bounded-result
adapters remain unchanged.
