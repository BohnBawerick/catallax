//! Normal-section construction for two parallel circular cylinders.

use catallax_numbers::{Decider, OperationSpan, Real, Refusal, Sign};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Cylinder, Line, Point3, Surface, Vector3};

use super::{FilletDefinition, bad};
use crate::profile::sign;

fn sqrt(value: &Real) -> Real {
    value
        .as_rational()
        .and_then(|r| crate::extrude::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| value.sqrt())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn cylinders(
    a: &Cylinder,
    b: &Cylinder,
    edge: &Line,
    same: [bool; 2],
    radius: &Real,
    k: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Surface, [Curve; 2]), Refusal> {
    let contact_span = OperationSpan::begin("brep", step, "fillet.contact", dec);
    for c in [a, b] {
        if sign(&c.axis.norm_squared(), dec)? == Sign::Zero
            || sign(&edge.dir.norm_squared(), dec)? == Sign::Zero
            || sign(&c.radius_sq, dec)? != Sign::Positive
        {
            return Err(Refusal::new(
                "fillet.invalid-edge",
                "degenerate cylinder or line",
            ));
        }
        if sign(&c.axis.cross(&edge.dir).norm_squared(), dec)? != Sign::Zero {
            return Err(bad("cylinder generator fillets require parallel axes"));
        }
    }
    // Both axis feet use the input edge's station. Axial placement and direction
    // scaling therefore cannot choose a different normal section or branch.
    let axes = [a, b].map(|c| {
        c.origin.plus(
            &c.axis
                .scaled(&(edge.origin.minus(&c.origin).dot(&c.axis) / c.axis.norm_squared())),
        )
    });
    for (c, axis) in [a, b].into_iter().zip(&axes) {
        if sign(
            &(&edge.origin.minus(axis).norm_squared() - &c.radius_sq),
            dec,
        )? != Sign::Zero
        {
            return Err(Refusal::new(
                "fillet.invalid-edge",
                "line is not on both cylinders",
            ));
        }
    }
    let delta = axes[1].minus(&axes[0]);
    let distance_sq = delta.norm_squared();
    if sign(&distance_sq, dec)? != Sign::Positive {
        return Err(bad("coaxial cylinders have no regular sharp generator"));
    }
    let radial = edge.origin.minus(&axes[0]);
    let side = radial.minus(&delta.scaled(&(&radial.dot(&delta) / &distance_sq)));
    let side_sq = side.norm_squared();
    if sign(&side_sq, dec)? != Sign::Positive {
        return Err(bad("tangent cylinders need singular contact continuation"));
    }
    let source_radii = [a.radius(), b.radius()];
    let offset_radii = [0, 1].map(|i| {
        &source_radii[i] - &(&(k * &Real::from_i64(if same[i] { 1 } else { -1 })) * radius)
    });
    for value in &offset_radii {
        if sign(value, dec)? != Sign::Positive {
            return Err(bad(
                "a cylinder offset reaches its axis and needs contact continuation",
            ));
        }
    }
    let along = (&distance_sq + &offset_radii[0].square() - offset_radii[1].square())
        / (&Real::from_i64(2) * &distance_sq);
    let height_sq = offset_radii[0].square() - &along.square() * &distance_sq;
    if sign(&height_sq, dec)? != Sign::Positive {
        return Err(bad(
            "offset cylinders have no regular contact branch at this radius",
        ));
    }
    // Keep the half-plane of the selected input generator. Both signed offset
    // radii are affine in the homotopy radius. Strict circle inequalities at
    // both ends prove this branch stays regular throughout that homotopy.
    let centre = axes[0]
        .plus(&delta.scaled(&along))
        .plus(&side.scaled(&sqrt(&(height_sq / side_sq))));
    let contacts = [0, 1].map(|i| {
        axes[i].plus(
            &centre
                .minus(&axes[i])
                .scaled(&(&source_radii[i] / &offset_radii[i])),
        )
    });
    contact_span.finish(dec, None, None);

    let section_span = OperationSpan::begin("brep", step, "fillet.sections", dec);
    let surface = Surface::Cylinder(Cylinder {
        ref_dir: contacts[0].minus(&centre),
        origin: centre,
        axis: edge.dir.clone(),
        radius_sq: radius.square(),
    });
    let rails = contacts.map(|origin| {
        Curve::Line(Line {
            origin,
            dir: edge.dir.clone(),
        })
    });
    section_span.finish(dec, None, None);
    Ok((surface, rails))
}

/// Whether `p` lies on the arc from `from` to `to`, a turn of less than half a
/// circle about `axis`. All three rays are normal to the axis.
fn on_arc(
    axis: &Vector3,
    [from, to]: [&Vector3; 2],
    p: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let turn = sign(&from.cross(to).dot(axis), dec)?;
    for pair in [from.cross(p), p.cross(to)] {
        let value = sign(&pair.dot(axis), dec)?;
        if value != Sign::Zero && value != turn {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether two arcs of a common circle, each under half a turn, share a point.
fn arcs_meet(
    axis: &Vector3,
    [first, second]: [[&Vector3; 2]; 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for (arc, ends) in [(first, second), (second, first)] {
        for end in ends {
            if on_arc(axis, arc, end, dec)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// The arc of a cylinder support between a strip's input edge and its contact
/// rail. The contact and the edge share a half-plane of the two axes, so the
/// arc is under half a turn.
fn consumed(strip: &FilletDefinition, support: usize) -> Option<(&Cylinder, [Vector3; 2])> {
    let (Surface::Cylinder(c), Curve::Line(rail), Curve::Line(edge)) = (
        &strip.supports[support],
        &strip.contacts[support],
        &strip.edge,
    ) else {
        return None;
    };
    let ray = |p: &Point3| {
        let d = p.minus(&c.origin);
        d.minus(&c.axis.scaled(&(d.dot(&c.axis) / c.axis.norm_squared())))
    };
    Some((c, [ray(&edge.origin), ray(&rail.origin)]))
}

pub(super) fn wall_arc_meets(
    strips: [(&FilletDefinition, usize); 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [first, second] = strips.map(|(strip, support)| consumed(strip, support));
    let (Some((wall, first)), Some((_, second))) = (first, second) else {
        return Ok(false);
    };
    arcs_meet(
        &wall.axis,
        [[&first[0], &first[1]], [&second[0], &second[1]]],
        dec,
    )
}

#[cfg(test)]
mod tests {
    use catallax_numbers::{Decider, Real};
    use catallax_shapes::Vector3;

    use super::arcs_meet;

    fn ray(x: i64, y: i64) -> Vector3 {
        Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::zero())
    }

    #[test]
    fn wall_arcs_meet_only_when_they_overlap_in_either_turn() {
        let up = Vector3::new(Real::zero(), Real::zero(), Real::one());
        let rays = [
            ray(5, 0),
            ray(4, 3),
            ray(3, 4),
            ray(0, 5),
            ray(-3, 4),
            ray(-5, 1),
        ];
        let mut dec = Decider::default();
        for (first, second, want) in [
            ([0, 1], [2, 3], false),
            ([1, 0], [3, 2], false),
            ([0, 2], [1, 3], true),
            ([2, 0], [1, 3], true),
            ([0, 3], [1, 2], true),
            ([3, 4], [5, 4], true),
            ([0, 1], [4, 5], false),
        ] {
            let [a, b] = first.map(|i| &rays[i]);
            let [c, d] = second.map(|i| &rays[i]);
            assert_eq!(
                arcs_meet(&up, [[a, b], [c, d]], &mut dec).unwrap(),
                want,
                "{first:?} {second:?}"
            );
        }
    }
}
