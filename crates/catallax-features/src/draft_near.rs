//! Exact squared distances to the trimmed faces of a plane-capped axial cylinder
//! or of a coaxial annulus between two such cylinders.
//!
//! Nearest-face references compare closed trimmed domains, never carriers: a
//! point past the rim is nearest the rim, whichever carrier it is nearer.

use super::sign;
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Cylinder, Plane, Point3};

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-near", detail)
}

/// Axial coordinate in axis lengths, and squared radial distance, of a point.
fn chart(cylinder: &Cylinder, point: &Point3) -> (Real, Real, Real) {
    let axis_sq = cylinder.axis.norm_squared();
    let offset = point.minus(&cylinder.origin);
    let t = offset.dot(&cylinder.axis) / axis_sq.clone();
    let radial = offset.minus(&cylinder.axis.scaled(&t));
    (t, radial.norm_squared(), axis_sq)
}

/// Axial coordinate of a cap perpendicular to the cylinder axis.
fn cap_coordinate(cylinder: &Cylinder, cap: &Plane, dec: &mut Decider) -> Result<Real, Refusal> {
    if sign(&cap.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero {
        return Err(unsupported(
            "nearest-face selection needs caps perpendicular to the cylinder axis",
        ));
    }
    Ok(chart(cylinder, &cap.origin).0)
}

/// `(rho - radius)^2` from the squared radial distance, avoiding a branch on sign.
fn radial_gap_squared(cylinder: &Cylinder, rho_sq: &Real) -> Real {
    let product = (rho_sq * &cylinder.radius_sq).sqrt();
    rho_sq + &cylinder.radius_sq - Real::from_i64(2) * product
}

/// Squared distance to the closed disc a perpendicular cap cuts from the cylinder.
pub fn cap_distance_squared(
    cylinder: &Cylinder,
    cap: &Plane,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let height = cap_coordinate(cylinder, cap, dec)?;
    let (t, rho_sq, axis_sq) = chart(cylinder, point);
    let axial = (&t - &height).square() * axis_sq;
    if sign(&(&rho_sq - &cylinder.radius_sq), dec)? != Sign::Positive {
        return Ok(axial);
    }
    Ok(axial + radial_gap_squared(cylinder, &rho_sq))
}

/// The two walls of an annulus, inner first. They must share one axis line and
/// differ in radius; a coincident pair is one physical cylinder, not an annulus.
pub fn annulus_walls<'a>(
    a: &'a Cylinder,
    b: &'a Cylinder,
    dec: &mut Decider,
) -> Result<[&'a Cylinder; 2], Refusal> {
    let offset = b.origin.minus(&a.origin);
    if sign(&a.axis.cross(&b.axis).norm_squared(), dec)? != Sign::Zero
        || sign(&offset.cross(&a.axis).norm_squared(), dec)? != Sign::Zero
    {
        return Err(unsupported(
            "nearest-face selection needs coaxial annular walls",
        ));
    }
    match sign(&(&a.radius_sq - &b.radius_sq), dec)? {
        Sign::Negative => Ok([a, b]),
        Sign::Positive => Ok([b, a]),
        Sign::Zero => Err(unsupported(
            "nearest-face selection needs annular walls of different radii",
        )),
    }
}

/// Squared distance to the closed annulus a perpendicular cap cuts between two
/// coaxial walls, inner wall first. The hole is not part of the face.
pub fn annular_cap_distance_squared(
    walls: [&Cylinder; 2],
    cap: &Plane,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let [inner, outer] = walls;
    let height = cap_coordinate(outer, cap, dec)?;
    let (t, rho_sq, axis_sq) = chart(outer, point);
    let axial = (&t - &height).square() * axis_sq;
    if sign(&(&rho_sq - &inner.radius_sq), dec)? == Sign::Negative {
        return Ok(axial + radial_gap_squared(inner, &rho_sq));
    }
    if sign(&(&rho_sq - &outer.radius_sq), dec)? == Sign::Positive {
        return Ok(axial + radial_gap_squared(outer, &rho_sq));
    }
    Ok(axial)
}

/// Squared distance to the cylindrical wall between two perpendicular caps.
pub fn wall_distance_squared(
    cylinder: &Cylinder,
    caps: [&Plane; 2],
    point: &Point3,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let heights = [
        cap_coordinate(cylinder, caps[0], dec)?,
        cap_coordinate(cylinder, caps[1], dec)?,
    ];
    let (low, high) = if sign(&(&heights[0] - &heights[1]), dec)? == Sign::Positive {
        (&heights[1], &heights[0])
    } else {
        (&heights[0], &heights[1])
    };
    let (t, rho_sq, axis_sq) = chart(cylinder, point);
    let excess = if sign(&(&t - low), dec)? == Sign::Negative {
        low - &t
    } else if sign(&(&t - high), dec)? == Sign::Positive {
        &t - high
    } else {
        Real::zero()
    };
    Ok(excess.square() * axis_sq + radial_gap_squared(cylinder, &rho_sq))
}

/// Index of the strictly smallest distance, or the existing ambiguous-reference error.
pub fn unique_nearest(distances: &[Real], dec: &mut Decider) -> Result<usize, Refusal> {
    let mut best = 0;
    let mut tied = false;
    for index in 1..distances.len() {
        match sign(&(&distances[index] - &distances[best]), dec)? {
            Sign::Negative => {
                best = index;
                tied = false;
            }
            Sign::Zero => tied = true,
            Sign::Positive => (),
        }
    }
    if tied {
        return Err(Refusal::new(
            "recipe.ambiguous-reference",
            "near point has equal minimum distance to several trimmed faces",
        ));
    }
    Ok(best)
}
