//! Shared exact model-space placements for the two recipe kernels.

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{Origin, PatternLayout, Plane, Role, StepId};
use catallax_shapes::{Isometry, RPoint3, RVector3};

/// Exact placement of one body-pattern copy, including copy zero.
pub fn pattern_placement(
    layout: &PatternLayout,
    instance: u32,
    dec: &mut Decider,
) -> Result<Isometry, Refusal> {
    if layout.count() == 0 {
        return Err(Refusal::new(
            "pattern.zero-count",
            "a pattern requires at least one copy",
        ));
    }
    if instance >= layout.count() {
        return Err(Refusal::new(
            "pattern.instance",
            "copy index is outside the pattern",
        ));
    }
    let n = Rational::from_i64(i64::from(instance));
    match layout {
        PatternLayout::Linear { vector, .. } => {
            let [x, y, z] = vector.0.each_ref().map(|v| v.value() * &n);
            Ok(Isometry::translation(RVector3::new(x, y, z)))
        }
        PatternLayout::Circular { axis, angle, .. } => {
            let [x, y, z] = axis.point.0.each_ref().map(|v| v.value().clone());
            let [dx, dy, dz] = axis.direction.0.each_ref().map(|v| v.value().clone());
            Isometry::around_axis(
                &RPoint3::new(x, y, z),
                &RVector3::new(dx, dy, dz),
                &(angle.value() * &n),
                dec,
            )
        }
    }
}

/// Exact reflection in a recipe plane with a normal of any nonzero length.
pub fn mirror_placement(plane: &Plane) -> Result<Isometry, Refusal> {
    let [x, y, z] = plane.point.0.each_ref().map(|v| v.value().clone());
    let [dx, dy, dz] = plane.normal.0.each_ref().map(|v| v.value().clone());
    Isometry::reflection(&RPoint3::new(x, y, z), &RVector3::new(dx, dy, dz))
}

/// Origin of a copied completed body. Copy zero retains its seed origins.
pub fn copy_origin(step: &StepId, instance: u32, source: &Origin) -> Origin {
    if instance == 0 {
        source.clone()
    } else {
        Origin {
            step: step.clone(),
            role: Role::Copy {
                instance,
                of: Box::new(source.clone()),
            },
        }
    }
}
