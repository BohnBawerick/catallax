//! Finite fillet-region bounds, never full-carrier clearance substitutes.
//!
//! OCCT's BRepBndLib bounds topological faces. SMLib's SmFace bounds the face's
//! UV domain or its finite edges. Here neighbour bounds come from Catallax's
//! proven trimmed-face extrema; blend bounds enclose their retained regions.

use super::*;
use catallax_measure::{Accuracy, BoundingBox};

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
                refusal.code,
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
