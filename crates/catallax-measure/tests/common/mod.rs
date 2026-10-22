#![allow(dead_code, unused_imports, unused_macros)]
pub mod boolean_cases;
pub mod fixtures;
pub mod oracle;
pub use catallax_measure::{self as measure, Accuracy, Enclosure, MassProperties};
pub use catallax_numbers::{
    Answer, Answered, BoundKind, BudgetPolicy, Decider, Mark, MarkLog, Rational,
};
pub use catallax_shapes::{Pose, Rotation};
pub use fixtures::*;
pub use oracle::{Expr, contains, exact, meets, pi, rat, sqrt, width};

// Completed comparison bodies run by default; unresolved bodies name their failure.
macro_rules! blind_test {
    ($name:ident, $body:block) => {
        #[test]
        fn $name() $body
    };
    ($name:ident, $reason:literal, $body:block) => {
        #[test]
        #[ignore = $reason]
        fn $name() $body
    };
}
pub(crate) use blind_test;

pub fn accuracy(digits: i32) -> Accuracy {
    Accuracy {
        relative: ten(-digits),
        absolute: ten(-40),
    }
}
pub fn default_accuracy() -> Accuracy {
    Accuracy {
        relative: ten(-12),
        absolute: ten(-30),
    }
}

pub fn assert_scalar(
    a: &Answer<Enclosure>,
    value: &Expr,
    requested: &Accuracy,
    require_accuracy: bool,
    require_exact: bool,
) {
    contains(&a.value, value);
    if require_exact {
        exact(&a.value, value);
    }
    if require_accuracy {
        assert!(
            meets(&a.value, requested),
            "requested width not met: {:?}",
            a.value
        );
    }
    let missed = if meets(&a.value, requested) {
        vec![]
    } else {
        vec![width(&a.value)]
    };
    assert_marks(&a.marks, &missed);
}
pub fn assert_marks(marks: &MarkLog, missed: &[Rational]) {
    assert_eq!(
        marks.len(),
        missed.len(),
        "one mark per missed property, no inherited body marks"
    );
    for (mark, w) in marks.iter().zip(missed) {
        assert_eq!(mark.code, "measure.accuracy-not-reached");
        assert_eq!(mark.kind, BoundKind::Value);
        assert_eq!(&mark.bound, w);
    }
}
pub fn assert_mass(
    a: &Answer<MassProperties>,
    e: &Expected,
    requested: &Accuracy,
    require_accuracy: bool,
) -> usize {
    let mut missed = vec![];
    for (actual, value) in [(&a.value.volume, &e.volume), (&a.value.area, &e.area)] {
        contains(actual, value);
        if e.exact {
            exact(actual, value);
        }
        if !meets(actual, requested) {
            missed.push(width(actual));
        }
    }
    let count = match (&a.value.centroid, &e.centroid) {
        (None, None) => 2,
        (Some(actual), Some(values)) => {
            // A Solid centroid must never divide an interval containing zero,
            // even when the requested absolute floor is larger than its volume.
            assert!(
                a.value.volume.lo > r(0),
                "solid volume not separated from zero"
            );
            for (c, value) in actual.iter().zip(values) {
                contains(c, value);
                if e.exact {
                    exact(c, value);
                }
            }
            if actual.iter().any(|c| !meets(c, requested)) {
                missed.push(actual.iter().map(width).max().unwrap());
            }
            5
        }
        _ => panic!("centroid must be None exactly for Empty"),
    };
    match (&a.value.inertia, &e.centroid) {
        (None, None) => {}
        (Some(tensor), Some(_)) => {
            if tensor
                .iter()
                .flatten()
                .any(|value| !meets(value, requested))
            {
                missed.push(tensor.iter().flatten().map(width).max().unwrap());
            }
        }
        _ => panic!("inertia must be None exactly for Empty"),
    }
    if require_accuracy {
        assert!(missed.is_empty(), "requested accuracy not reached: {a:?}");
    }
    assert_marks(&a.marks, &missed);
    count
}
pub fn run_mass(f: &Fixture, a: &Accuracy) -> Answer<MassProperties> {
    measure::mass_properties(&f.model, f.body, a, &mut Decider::default()).unwrap()
}
pub fn assert_faces(f: &Fixture, a: &Accuracy) {
    assert_eq!(
        f.model.faces_of_body(f.body),
        f.faces.iter().map(|(id, _)| *id).collect::<Vec<_>>()
    );
    let mut sum = Enclosure { lo: r(0), hi: r(0) };
    for (id, expected) in &f.faces {
        let value = measure::face_area(&f.model, *id, a, &mut Decider::default()).unwrap();
        assert_scalar(&value, expected, a, true, f.expected.exact);
        sum.lo = &sum.lo + &value.value.lo;
        sum.hi = &sum.hi + &value.value.hi;
    }
    // Independent closed form and interval-sum identity, not equality of widths.
    contains(&sum, &f.expected.area);
    let total = measure::area(&f.model, f.body, a, &mut Decider::default()).unwrap();
    assert_scalar(&total, &f.expected.area, a, true, f.expected.exact);
    assert!(sum.lo <= total.value.hi && total.value.lo <= sum.hi);
}
pub fn upstream_budget_refusal(code: &str) {
    // Supported, valid fixtures cannot legitimately hit constructor/domain errors.
    // trace propagates NotExact codes as refusals, including precision-limit.
    assert!(
        matches!(
            code,
            "numbers.budget-exhausted"
                | "numbers.precision-limit"
                | "numbers.divisor-not-separated"
                | "numbers.sqrt-not-separated"
                | "numbers.root-not-certified"
        ),
        "unexpected upstream refusal: {code}"
    );
}
