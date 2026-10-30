//! Boundary rows of a rolling-envelope face.

use catallax_numbers::{Decider, Outcome, Refusal, Sign};
use catallax_shapes::{BoundaryRow, Curve, Point3, SplineSupport, Surface};

use crate::{EdgeId, FaceId, Loop, Model, Sense};

fn unsupported(detail: &'static str) -> Refusal {
    Refusal::new("topology.rolling-envelope-boundary", detail)
}

fn same_point(a: &Point3, c: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(&a.minus(c).norm_squared()) {
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::Proven(_) => Ok(false),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(refusal) => Err(refusal),
    }
}

/// Each edge of an envelope face with its chart row and its use by the face.
/// Meridians end on the exact corners of one transition. The rail is the only
/// spline edge and the mouth is the remaining arc.
pub fn envelope_rows(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Vec<(EdgeId, BoundaryRow, Sense)>, Refusal> {
    let Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) =
        model.geom().surface(model.face(face).surface)
    else {
        return Err(unsupported("the face is not a rolling envelope"));
    };
    let negative = envelope.definition().corners(Sign::Negative);
    let positive = envelope.definition().corners(Sign::Positive);
    let [boundary] = model.face(face).loops.as_slice() else {
        return Err(unsupported("an envelope face has one outer boundary"));
    };
    let Loop::Coedges(coedges) = model.loop_(*boundary) else {
        return Err(unsupported("an envelope face has one outer boundary"));
    };
    let mut rows = Vec::with_capacity(4);
    for coedge in coedges {
        let use_ = model.coedge(*coedge);
        let edge = model.edge(use_.edge);
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Err(unsupported("an envelope edge has two ends"));
        };
        let points = [start, end].map(|v| model.geom().position(model.vertex(v).point));
        let mut on = |corners: &[Point3; 2]| -> Result<bool, Refusal> {
            Ok((same_point(&points[0], &corners[0], dec)?
                && same_point(&points[1], &corners[1], dec)?)
                || (same_point(&points[0], &corners[1], dec)?
                    && same_point(&points[1], &corners[0], dec)?))
        };
        let row = if matches!(model.geom().curve(edge.curve), Curve::BSpline(_)) {
            BoundaryRow::Rail
        } else if on(&negative)? {
            BoundaryRow::Negative
        } else if on(&positive)? {
            BoundaryRow::Positive
        } else {
            BoundaryRow::Mouth
        };
        rows.push((use_.edge, row, use_.sense));
    }
    Ok(rows)
}
