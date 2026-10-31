//! Heights of a finite face along one direction, in any placement. A world box
//! of a tilted face has corners the face never reaches; these values depend
//! only on the face and the direction, so a rigid motion leaves them unchanged.

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::{Point3, Surface, Vector3};
use catallax_topology::{EdgeId, FaceId, Loop, Model};
use std::collections::BTreeSet;

/// Exact reals whose least and greatest enclose `direction . (p - origin)` over
/// the face, or `None` when the face has no such closed form.
///
/// Height is affine along every ruling of a plane, cylinder or cone, so its
/// extremes over a compact face lie on the boundary or at an apex vertex loop.
/// An arc is bounded by its full circle, so an arc's values may be outward.
pub(crate) fn face(
    model: &Model,
    id: FaceId,
    origin: &Point3,
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Vec<Real>>, Refusal> {
    let face = model.face(id);
    if !matches!(
        model.geom().surface(face.surface),
        Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_)
    ) {
        return Ok(None);
    }
    let height = |p: &Point3| p.minus(origin).dot(direction);
    let position = |v| model.geom().position(model.vertex(v).point);
    let mut heights = Vec::new();
    let mut edges = BTreeSet::new();
    for lp in &face.loops {
        match model.loop_(*lp) {
            Loop::Vertex(v) => heights.push(height(&position(*v))),
            Loop::Coedges(cs) => edges.extend(cs.iter().map(|c| model.coedge(*c).edge)),
        }
    }
    for edge in edges {
        let Some(values) = self::edge(model, edge, origin, direction, dec)? else {
            return Ok(None);
        };
        heights.extend(values);
    }
    Ok((!heights.is_empty()).then_some(heights))
}

/// Exact reals whose least and greatest enclose `direction . (p - origin)` over
/// an edge, or `None` when the edge has no such closed form. An arc is bounded
/// by its full circle, as on a face.
pub(crate) fn edge(
    model: &Model,
    id: EdgeId,
    origin: &Point3,
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Vec<Real>>, Refusal> {
    let height = |p: &Point3| p.minus(origin).dot(direction);
    let position = |v| model.geom().position(model.vertex(v).point);
    let edge = model.edge(id);
    Ok(Some(match curve_form(model.geom(), edge.curve, dec)? {
        CurveForm::Line { .. } => {
            let (Some(start), Some(end)) = (edge.start, edge.end) else {
                return Ok(None);
            };
            vec![height(&position(start)), height(&position(end))]
        }
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } => {
            let along = &direction.dot(&normal) / &normal.norm_squared();
            let squared = direction.minus(&normal.scaled(&along)).norm_squared();
            let middle = height(&center);
            if super::nonzero(&squared, dec)? {
                let reach = (&radius_sq * &squared).sqrt();
                vec![&middle + &reach, &middle - &reach]
            } else {
                vec![middle]
            }
        }
        _ => return Ok(None),
    }))
}
