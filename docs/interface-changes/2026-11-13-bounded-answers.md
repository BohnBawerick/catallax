# Bounded operation answers

Date: 13 November 2026. Affected pieces: numbers, STEP and CLI verification.

Implements the [owner's decision](../ARCHITECTURE.md#named-case-levels-and-bounded-answers-decision-of-2026-11-13). This adds a public result contract without changing any existing operation's answer or tolerance.

- `BoundedRule` names a case kind and versioned rule, with an exact rational level and units fixed in code before evaluation. There are no production rules in this change. A request cannot supply a rule or tune its level after work starts.
- `BoundedEvidence::new` accepts only a nonnegative bound at or below that level, in matching known units. It retains `CertificateIdentity<String>`, including the kept definition, recipe and representation revisions and complete measured scope. Producers must establish the bound over that scope. Structural validation does not prove the measurement method.
- `BoundedAnswer<T>` stores the value and checked evidence separately from `Answer<T>`. Neither it nor `OperationResult<T>` converts into a proven answer. Even a zero measured bound retains bounded status.
- `OperationResult<T>` distinguishes `Proven`, `Bounded` and `Refused`. Existing `Answered<T>` values convert only to their corresponding proven or refused variant. Decision-level `Outcome<T>` is unchanged.
- `StepResult` names `OperationResult<StepReport>`. `StepSummary::from_result` preserves bounded evidence and refuses to summarize a refused export as a successful file.
- CLI `Status::Bounded`, `BoundedSummary` and optional `bounded` fields on `Manifest` and `StepSummary` carry the same record. Existing manifests omit absent records, preserving their previous byte layout. `Manifest::with_bounded_evidence` sets bounded status and clears `exact` and any refusal.
- `verification::operation_result_summary` serializes `catallax-operation-result/1`. The record has `status` and `proof`, a result value or refusal, and `exact: false` for bounded answers. Proven results keep diagnostic marks. Bounded results include `bounded`, with `case_kind`, `rule`, `level`, `measured_bound`, `unit` and `definition`. Rational values are exact fractions, without decimal rounding. Definition scope is `whole-solid`, `entities`, `parameter-box` or `query`. Units are `mm`, `rad`, `scalar`, `mm2` or `mm3`; raw residuals are rejected.

New stable refusal codes are `numbers.invalid-bounded-answer` for invalid rule, units, bounds or missing definition identities, and `numbers.bounded-level-exceeded` when the measured bound exceeds the preset level. Scope structure errors retain `numbers.geometric-gap-unproven` from the existing certificate validator. Existing refusal codes and text do not change.

Tangent touch, helix and general sweeps, surface fitting and freeform Boolean classification remain follow-ups. Each must choose its rule in advance and establish its own kept-definition bound. No preview, sampled agreement, weakened proof, raised allowance or tolerance retry earns bounded status by itself.

Validation includes separate wire and display records, compile-fail coverage for promotion to `Answer<T>`, invalid bound and scope checks, legacy STEP and refusal serialization checks, and the unchanged fresh 477-case replay in hosted CI.
