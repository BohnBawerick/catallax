//! Finite fillet-region bounds, never full-carrier clearance substitutes.
//!
//! OCCT's BRepBndLib bounds topological faces. SMLib's SmFace bounds the face's
//! UV domain or its finite edges. Here neighbour bounds come from Catallax's
//! proven trimmed-face extrema; blend bounds enclose their retained regions.

use super::*;
use crate::fillet_corner::SphereSector;
use catallax_measure::{Accuracy, BoundingBox};

/// The full sweep of the meridian triangle used by `arc_hull`. This applies
/// only to straight meridian supports. Curved meridians need their own bounds.
pub(super) fn ring_hull(
    definition: &FilletDefinition,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    let (Curve::Circle(edge), [Curve::Circle(a), Curve::Circle(b)]) =
        (&definition.edge, &definition.contacts)
    else {
        return Ok(None);
    };
    if definition.supports.iter().any(|surface| {
        !matches!(
            surface,
            Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_)
        )
    }) {
        return Ok(None);
    }
    let mut low = [Real::zero(), Real::zero(), Real::zero()];
    let mut high = low.clone();
    for (index, circle) in [edge, a, b].into_iter().enumerate() {
        let squared = circle.normal.norm_squared();
        for axis in 0..3 {
            let normal = [&circle.normal.x, &circle.normal.y, &circle.normal.z][axis];
            let centre = [&circle.center.x, &circle.center.y, &circle.center.z][axis];
            let extent = (&circle.radius_sq * &(Real::one() - normal.square() / &squared)).sqrt();
            let lo = centre - &extent;
            let hi = centre + &extent;
            if index == 0 || sign(&(&lo - &low[axis]), dec)? == Sign::Negative {
                low[axis] = lo;
            }
            if index == 0 || sign(&(&hi - &high[axis]), dec)? == Sign::Positive {
                high[axis] = hi;
            }
        }
    }
    Ok(Some(
        (0..8)
            .map(|mask| {
                Point3::new(
                    if mask & 1 == 0 {
                        low[0].clone()
                    } else {
                        high[0].clone()
                    },
                    if mask & 2 == 0 {
                        low[1].clone()
                    } else {
                        high[1].clone()
                    },
                    if mask & 4 == 0 {
                        low[2].clone()
                    } else {
                        high[2].clone()
                    },
                )
            })
            .collect(),
    ))
}

/// The sector's rays meet the plane `dual_sum . (p-centre) = r*|dual_sum|`
/// in a triangle. Cauchy-Schwarz places the spherical sector inside the
/// pyramid joining that triangle to the centre. No sampling enters this bound.
pub(super) fn corner_hull(sector: &SphereSector) -> [Point3; 4] {
    let [a, b, c] = sector.dual_basis();
    let dual_sum = a.plus(&b).plus(&c);
    let reach = &sector.radius * &dual_sum.norm_squared().sqrt();
    [
        sector.centre.clone(),
        sector.centre.plus(&sector.directions[0].scaled(&reach)),
        sector.centre.plus(&sector.directions[1].scaled(&reach)),
        sector.centre.plus(&sector.directions[2].scaled(&reach)),
    ]
}

pub(super) fn hull_box(points: &[Point3], dec: &Decider) -> Result<BoundingBox, Refusal> {
    let first = points.first().ok_or_else(|| {
        Refusal::new(
            "fillet.invalid-trim",
            "a retained fillet region has an empty hull",
        )
    })?;
    let values = [
        dec.enclose_at(&first.x, 128)?,
        dec.enclose_at(&first.y, 128)?,
        dec.enclose_at(&first.z, 128)?,
    ];
    let mut bounds = BoundingBox {
        min: values.each_ref().map(|value| value.lo()),
        max: values.each_ref().map(|value| value.hi()),
    };
    for point in &points[1..] {
        for (axis, coordinate) in [&point.x, &point.y, &point.z].into_iter().enumerate() {
            let value = dec.enclose_at(coordinate, 128)?;
            let (lo, hi) = (value.lo(), value.hi());
            if lo < bounds.min[axis] {
                bounds.min[axis] = lo;
            }
            if hi > bounds.max[axis] {
                bounds.max[axis] = hi;
            }
        }
    }
    Ok(bounds)
}

pub(super) fn outside_face(
    model: &Model,
    face: FaceId,
    region: &BoundingBox,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let bounds = match catallax_measure::face_bounding_box(model, face, &Accuracy::default(), dec) {
        Ok(answer) => answer.value,
        Err(refusal)
            if matches!(
                refusal.code.as_str(),
                "measure.unsupported-surface"
                    | "measure.unsupported-curve"
                    | "measure.rational-patch-domain"
            ) =>
        {
            return Ok(false);
        }
        Err(refusal) => return Err(refusal),
    };
    Ok((0..3)
        .any(|axis| region.max[axis] < bounds.min[axis] || bounds.max[axis] < region.min[axis]))
}
