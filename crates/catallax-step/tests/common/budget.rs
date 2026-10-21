// The local budget refusal plus stable codes from pieces 1, 2 and 4. Trace can
// promote a numbers NotExact code to a refusal, including precision-limit.
pub fn assert_budget_refusal(code: &str, before: &[u8], after: &[u8]) {
    assert!(
        matches!(
            code,
            "step.budget-exhausted"
                | "numbers.budget-exhausted"
                | "numbers.precision-limit"
                | "numbers.division-by-zero"
                | "numbers.divisor-not-separated"
                | "numbers.sqrt-of-negative"
                | "numbers.sqrt-not-separated"
                | "numbers.zero-polynomial"
                | "numbers.root-not-certified"
                | "numbers.deadline"
                | "shapes.zero-quaternion"
                | "shapes.zero-vector"
                | "shapes.not-perpendicular"
                | "shapes.non-positive-radius"
                | "shapes.non-positive-slope"
                | "shapes.torus-self-intersecting"
                | "shapes.ellipse-axes"
                | "intersect.unsupported-pair"
                | "intersect.unsupported-curve"
        ),
        "only the local budget code or unchanged upstream codes are allowed, got {code}"
    );
    assert_eq!(after, before, "refused export changed out");
}
