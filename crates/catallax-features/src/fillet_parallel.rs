//! Normal-section construction for straight generators: two parallel circular
//! cylinders, or a cylinder and a plane parallel to its axis.

use catallax_numbers::{Decider, OperationSpan, Real, Refusal, Sign};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Cylinder, Line, Plane, Point3, Surface, Vector3};

use super::{FilletDefinition, bad, normal};
use crate::profile::sign;

fn sqrt(value: &Real) -> Real {
    value
        .as_rational()
        .and_then(|r| crate::extrude::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| value.sqrt())
}

/// The foot on a cylinder's axis of the normal section through the edge's
/// station. Axial placement and direction scaling then cannot choose a
/// different section or branch.
fn axis_foot(cylinder: &Cylinder, edge: &Line) -> Point3 {
    cylinder.origin.plus(&cylinder.axis.scaled(
        &(edge.origin.minus(&cylinder.origin).dot(&cylinder.axis) / cylinder.axis.norm_squared()),
    ))
}

/// Checks a cylinder has an exact generator along the edge, and returns its
/// axis foot in the edge's normal section.
fn generator(cylinder: &Cylinder, edge: &Line, dec: &mut Decider) -> Result<Point3, Refusal> {
    if sign(&cylinder.axis.norm_squared(), dec)? == Sign::Zero
        || sign(&edge.dir.norm_squared(), dec)? == Sign::Zero
        || sign(&cylinder.radius_sq, dec)? != Sign::Positive
    {
        return Err(Refusal::new(
            "fillet.invalid-edge",
            "degenerate cylinder or line",
        ));
    }
    if sign(&cylinder.axis.cross(&edge.dir).norm_squared(), dec)? != Sign::Zero {
        return Err(bad("cylinder generator fillets require parallel axes"));
    }
    let foot = axis_foot(cylinder, edge);
    if sign(
        &(&edge.origin.minus(&foot).norm_squared() - &cylinder.radius_sq),
        dec,
    )? != Sign::Zero
    {
        return Err(Refusal::new(
            "fillet.invalid-edge",
            "line is not on both supports",
        ));
    }
    Ok(foot)
}

/// The signed radius of a cylinder's offset at the rolling radius, which must
/// stay on the contact's side of the axis.
fn offset_radius(
    cylinder: &Cylinder,
    same: bool,
    radius: &Real,
    k: &Real,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let offset = &cylinder.radius() - &(&(k * &Real::from_i64(if same { 1 } else { -1 })) * radius);
    if sign(&offset, dec)? != Sign::Positive {
        return Err(bad(
            "a cylinder offset reaches its axis and needs contact continuation",
        ));
    }
    Ok(offset)
}

fn strip(
    centre: Point3,
    contacts: [Point3; 2],
    edge: &Line,
    radius: &Real,
) -> (Surface, [Curve; 2]) {
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
    (surface, rails)
}

/// The contact of a ball centred on a plane's inward offset and on a cylinder's
/// offset, for a plane parallel to the cylinder's axis. Contacts follow the
/// `[plane, cylinder]` order.
#[allow(clippy::too_many_arguments)]
pub(super) fn plane_cylinder(
    plane: &Plane,
    cylinder: &Cylinder,
    edge: &Line,
    same: [bool; 2],
    radius: &Real,
    k: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Surface, [Curve; 2]), Refusal> {
    let contact_span = OperationSpan::begin("brep", step, "fillet.contact", dec);
    if sign(&plane.normal.norm_squared(), dec)? == Sign::Zero {
        return Err(Refusal::new("fillet.invalid-edge", "degenerate plane"));
    }
    let axis = generator(cylinder, edge, dec)?;
    if sign(&plane.normal.dot(&edge.dir), dec)? != Sign::Zero
        || sign(&edge.origin.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
    {
        return Err(Refusal::new(
            "fillet.invalid-edge",
            "line is not on both supports",
        ));
    }
    let outward = normal(&plane.normal, same[0]);
    let radial = edge.origin.minus(&axis);
    let height = radial.dot(&outward);
    let side = radial.minus(&outward.scaled(&height));
    let side_sq = side.norm_squared();
    if sign(&side_sq, dec)? != Sign::Positive {
        return Err(bad(
            "a tangent cylinder and plane need singular contact continuation",
        ));
    }
    let offset = offset_radius(cylinder, same[1], radius, k, dec)?;
    let across = &height - &(k * radius);
    let chord_sq = offset.square() - across.square();
    if sign(&chord_sq, dec)? != Sign::Positive {
        return Err(bad(
            "the offset line misses the offset circle at this radius",
        ));
    }
    // Keep the chord end on the selected generator's side. `chord_sq` is the
    // product of `offset - across` and `offset + across`, both affine in the
    // homotopy radius, so strict signs at both ends hold throughout.
    let centre = axis
        .plus(&outward.scaled(&across))
        .plus(&side.scaled(&sqrt(&(chord_sq / side_sq))));
    let contacts = [
        centre.plus(&outward.scaled(&(k * radius))),
        axis.plus(&centre.minus(&axis).scaled(&(&cylinder.radius() / &offset))),
    ];
    contact_span.finish(dec, None, None);

    let section_span = OperationSpan::begin("brep", step, "fillet.sections", dec);
    let answer = strip(centre, contacts, edge, radius);
    section_span.finish(dec, None, None);
    Ok(answer)
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
    let axes = [generator(a, edge, dec)?, generator(b, edge, dec)?];
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
    let offset_radii = [
        offset_radius(a, same[0], radius, k, dec)?,
        offset_radius(b, same[1], radius, k, dec)?,
    ];
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
    let answer = strip(centre, contacts, edge, radius);
    section_span.finish(dec, None, None);
    Ok(answer)
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

/// The band of one support between a strip's input edge and its contact rail,
/// as rays from the axis of a cylinder or positions across a plane.
enum Band {
    Arc(Vector3, [Vector3; 2]),
    Across(Vector3, [Real; 2]),
}

/// The contact and the edge share a half-plane of the section, whose boundary
/// passes through the axis, so an arc band is under half a turn.
fn band(strip: &FilletDefinition, support: usize) -> Option<Band> {
    let (Curve::Line(rail), Curve::Line(edge)) = (&strip.contacts[support], &strip.edge) else {
        return None;
    };
    match &strip.supports[support] {
        Surface::Cylinder(c) => {
            let ray = |p: &Point3| {
                let d = p.minus(&c.origin);
                d.minus(&c.axis.scaled(&(d.dot(&c.axis) / c.axis.norm_squared())))
            };
            Some(Band::Arc(
                c.axis.clone(),
                [ray(&edge.origin), ray(&rail.origin)],
            ))
        }
        Surface::Plane(p) => {
            let across = p.normal.cross(&rail.dir);
            let at = |q: &Point3| q.minus(&p.origin).dot(&across);
            Some(Band::Across(
                rail.dir.clone(),
                [at(&edge.origin), at(&rail.origin)],
            ))
        }
        _ => None,
    }
}

fn ordered(pair: &[Real; 2], dec: &mut Decider) -> Result<[Real; 2], Refusal> {
    Ok(if sign(&(&pair[0] - &pair[1]), dec)? == Sign::Positive {
        [pair[1].clone(), pair[0].clone()]
    } else {
        pair.clone()
    })
}

pub(super) fn bands_meet(
    strips: [(&FilletDefinition, usize); 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    // Plane-plane strips keep their carrier and kite proofs.
    if strips
        .iter()
        .all(|(strip, _)| matches!(strip.supports, [Surface::Plane(_), Surface::Plane(_)]))
    {
        return Ok(false);
    }
    let [first, second] = strips.map(|(strip, support)| band(strip, support));
    match (first, second) {
        (Some(Band::Arc(axis, first)), Some(Band::Arc(_, second))) => arcs_meet(
            &axis,
            [[&first[0], &first[1]], [&second[0], &second[1]]],
            dec,
        ),
        (Some(Band::Across(dir, first)), Some(Band::Across(other, second))) => {
            if sign(&dir.cross(&other).norm_squared(), dec)? != Sign::Zero {
                return Ok(false);
            }
            // Each position follows its own rail's direction, so antiparallel
            // rails measure across the plane in opposite senses.
            let flip = sign(&dir.dot(&other), dec)? == Sign::Negative;
            let second = if flip {
                [-&second[0], -&second[1]]
            } else {
                second
            };
            let [a_low, a_high] = ordered(&first, dec)?;
            let [b_low, b_high] = ordered(&second, dec)?;
            Ok(sign(&(&a_high - &b_low), dec)? != Sign::Negative
                && sign(&(&b_high - &a_low), dec)? != Sign::Negative)
        }
        _ => Ok(false),
    }
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
