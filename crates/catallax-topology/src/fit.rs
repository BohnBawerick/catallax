//! The fit mark a body's fitted rolling-ball faces call for.
//!
//! A fitted carrier keeps the bound its fit's own two-way check measured, and
//! a trimmed or rigidly moved face of it keeps that bound. A body's fit mark is
//! therefore the worst bound over the fitted faces the body has, not a bound
//! copied from the operands of a later feature.

use crate::{BodyId, Model};
use catallax_numbers::{BoundKind, Mark, MarkLog, Rational};
use catallax_shapes::{SplineSupport, Surface};

/// Code of the mark for fitted rolling-envelope, rolling-saddle and
/// rolling-sphere faces. Its bound is the worst fitted carrier's measured
/// distance to its rolling-ball definition, in model units.
pub const FIT_MARK: &str = "features.rolling-envelope-fit";

/// The worst measured fit bound over the body's fitted rolling-ball faces, or
/// `None` when it has none.
pub fn fitted_face_bound(model: &Model, body: BodyId) -> Option<Rational> {
    model
        .faces_of_body(body)
        .into_iter()
        .filter_map(
            |face| match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RollingEnvelope(fit)) => {
                    Some(fit.bound().upper.clone())
                }
                Surface::BSpline(SplineSupport::RollingSaddle(fit)) => {
                    Some(fit.bound().upper.clone())
                }
                Surface::BSpline(SplineSupport::RollingSphere(fit)) => {
                    Some(fit.bound().upper.clone())
                }
                _ => None,
            },
        )
        .max()
}

/// `marks` with its fit marks replaced by the one the body's own fitted faces
/// call for. A fit mark that already states that bound is kept as it is,
/// with its place; other marks keep their order.
pub fn with_fit_mark(model: &Model, body: BodyId, marks: &MarkLog) -> MarkLog {
    fit_mark_for(fitted_face_bound(model, body), marks)
}

/// `marks` without fit marks, for a result with no fitted face, such as an
/// empty body.
pub fn without_fit_mark(marks: &MarkLog) -> MarkLog {
    fit_mark_for(None, marks)
}

fn fit_mark_for(bound: Option<Rational>, marks: &MarkLog) -> MarkLog {
    let mut kept = MarkLog::new();
    let mut stated = false;
    for mark in marks.iter() {
        if mark.code != FIT_MARK {
            kept.push(mark.clone());
        } else if !stated && bound.as_ref() == Some(&mark.bound) && mark.kind == BoundKind::Distance
        {
            kept.push(mark.clone());
            stated = true;
        }
    }
    if let (Some(bound), false) = (bound, stated) {
        kept.push(Mark::new(
            FIT_MARK,
            BoundKind::Distance,
            bound,
            "fitted rolling-ball faces; the mark alone is not authority",
        ));
    }
    kept
}
