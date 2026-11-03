//! Diagnose a failed equal-cylinder strip whose input arc reaches a pole.
//!
//! At a crossing pole the two cylinder-axis feet coincide. Equal source and
//! offset radii then give the same ball contact twice. Its spherical contact
//! cone has zero area, so a three-independent-contact SphereSector cannot close
//! it. This is missing collapsed-boundary topology, not impossible geometry.

use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign, ZeroProof};
use catallax_shapes::{Curve, Cylinder, Point3, Surface};
use catallax_topology::{EdgeId, FaceId, Model, Sense};

use crate::{fillet::UNSUPPORTED_CORNER, profile::sign};

fn exact_zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::ProvenCoincident(proof) => Ok(proof.proof != ZeroProof::SeparationBound),
        Outcome::Refused(error) => Err(error),
        _ => Ok(false),
    }
}

fn axis_foot(cylinder: &Cylinder, point: &Point3) -> Point3 {
    cylinder.origin.plus(&cylinder.axis.scaled(
        &(point.minus(&cylinder.origin).dot(&cylinder.axis) / cylinder.axis.norm_squared()),
    ))
}

fn reaches_collapsed_pole(
    model: &Model,
    edge: EdgeId,
    faces: [FaceId; 2],
    radius: &Real,
    convex: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [first, second] = faces.map(|face| model.face(face));
    let (Surface::Cylinder(a), Surface::Cylinder(b)) = (
        model.geom().surface(first.surface),
        model.geom().surface(second.surface),
    ) else {
        return Ok(false);
    };
    let edge = model.edge(edge);
    if first.sense != second.sense
        || !matches!(
            model.geom().curve(edge.curve),
            Curve::Ellipse(_) | Curve::Intersection(_)
        )
        || !exact_zero(&(&a.radius_sq - &b.radius_sq), dec)?
    {
        return Ok(false);
    }
    // Equal face senses give equal signed offsets. A zero or negative inward
    // offset has a different failure; no positive-radius pole is claimed there.
    if (first.sense == Sense::Same) == convex
        && sign(&(&a.radius_sq - &radius.square()), dec)? != Sign::Positive
    {
        return Ok(false);
    }
    let normal = a.axis.cross(&b.axis);
    if sign(&normal.norm_squared(), dec)? != Sign::Positive
        || !exact_zero(&b.origin.minus(&a.origin).dot(&normal), dec)?
    {
        return Ok(false);
    }
    if edge.start.is_none() && edge.end.is_none() {
        // A whole equal-cylinder intersection branch includes both poles.
        return Ok(true);
    }
    for vertex in [edge.start, edge.end].into_iter().flatten() {
        let point = model.geom().position(model.vertex(vertex).point);
        if exact_zero(
            &axis_foot(a, &point)
                .minus(&axis_foot(b, &point))
                .norm_squared(),
            dec,
        )? {
            return Ok(true);
        }
    }
    // A trimmed arc that stops before either pole is not this corner case.
    Ok(false)
}

pub(super) fn corner_refusal(
    model: &Model,
    edge: EdgeId,
    faces: [FaceId; 2],
    radius: &Real,
    convex: bool,
    refusal: Refusal,
    dec: &mut Decider,
) -> Refusal {
    if refusal.code == "fillet.unsupported-contact"
        && matches!(
            reaches_collapsed_pole(model, edge, faces, radius, convex, dec),
            Ok(true)
        )
    {
        Refusal::new(
            UNSUPPORTED_CORNER,
            "equal-radius cylinders on crossing axes have coincident pole contacts; their zero-area spherical sector needs collapsed-boundary topology, not a finite SphereSector patch",
        )
    } else {
        refusal
    }
}
