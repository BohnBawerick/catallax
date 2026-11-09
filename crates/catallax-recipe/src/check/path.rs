//! Rational spline invariants that need no geometric decisions.

use catallax_numbers::Sign;

use crate::{PathSegment, Point3, RecipeError, StepId};

pub(super) fn check(
    step: &StepId,
    start: &Point3,
    segments: &[PathSegment],
) -> Result<(), RecipeError> {
    let mut previous = start;
    for segment in segments {
        match segment {
            PathSegment::Line { to, .. } | PathSegment::Arc { to, .. } => previous = to,
            PathSegment::Spline {
                id,
                degree,
                control_points,
                weights,
                knots,
                ..
            } => {
                let invalid = |detail: &str| RecipeError::PathSpline {
                    step: step.clone(),
                    segment: id.clone(),
                    detail: detail.into(),
                };
                let width = degree
                    .checked_add(1)
                    .filter(|&width| *degree > 0 && width <= control_points.len())
                    .ok_or_else(|| {
                        invalid("degree must be positive and less than the control count")
                    })?;
                let count = control_points.len();
                if count.checked_add(width) != Some(knots.len()) {
                    return Err(invalid(
                        "knot count must equal controls plus degree plus one",
                    ));
                }
                if weights.len() != count
                    || weights.iter().any(|w| w.value().sign() != Sign::Positive)
                {
                    return Err(invalid(
                        "one strictly positive weight is required per control",
                    ));
                }
                if knots.windows(2).any(|pair| pair[0] > pair[1]) {
                    return Err(invalid("knots must be nondecreasing"));
                }
                let one = crate::Number::from(1);
                if !knots[..width].iter().all(|k| k.value().is_zero())
                    || !knots[count..].iter().all(|k| k == &one)
                {
                    return Err(invalid("knots must be clamped to zero and one"));
                }
                let mut previous_knot = None;
                let mut multiplicity = 0;
                for knot in &knots[width..count] {
                    if knot.value().sign() != Sign::Positive || knot >= &one {
                        return Err(invalid(
                            "interior knots must lie strictly between zero and one",
                        ));
                    }
                    multiplicity = if previous_knot == Some(knot) {
                        multiplicity + 1
                    } else {
                        1
                    };
                    if multiplicity > *degree {
                        return Err(invalid("interior knot multiplicity exceeds the degree"));
                    }
                    previous_knot = Some(knot);
                }
                if control_points.first() != Some(previous) {
                    return Err(invalid(
                        "first control must equal the preceding endpoint exactly",
                    ));
                }
                previous = control_points
                    .last()
                    .expect("validated positive control count");
            }
        }
    }
    Ok(())
}
