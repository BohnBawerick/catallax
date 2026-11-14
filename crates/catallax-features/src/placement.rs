//! Shared exact model-space placements for the two recipe kernels.

use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{Origin, PatternLayout, Plane, Role, StepId};
use catallax_shapes::{Isometry, Pose, RPoint3, RVector3, Rotation, Vector3};

/// Places local XY on one supporting plane, oriented toward `normal_direction`.
/// Split faces may supply the same plane; differing planes are an ambiguous reference.
pub fn face_placement(
    planes: &[catallax_shapes::Plane],
    normal_direction: &RVector3,
    x_direction: &RVector3,
    offset: &Rational,
    placement: &catallax_recipe::Pose,
    dec: &mut Decider,
) -> Result<Pose, Refusal> {
    let Some(first) = planes.first() else {
        return Err(Refusal::new(
            "reference.lost-face",
            "no surviving face supplies the sketch plane",
        ));
    };
    for other in &planes[1..] {
        let parallel = first.normal.cross(&other.normal).norm_squared();
        let separation = first.normal.dot(&other.origin.minus(&first.origin));
        if crate::profile::sign(&parallel, dec)? != Sign::Zero
            || crate::profile::sign(&separation, dec)? != Sign::Zero
        {
            return Err(Refusal::new(
                "reference.ambiguous-face",
                "matching faces do not share one plane",
            ));
        }
    }
    let sense = crate::profile::sign(&first.normal.dot(&normal_direction.to_real()), dec)?;
    let normal = match sense {
        Sign::Positive => first.normal.clone(),
        Sign::Negative => first.normal.scaled(&Real::from_i64(-1)),
        Sign::Zero => {
            return Err(Refusal::new(
                "extrude.face-normal",
                "normal direction lies in the supporting plane",
            ));
        }
    };
    let rotation = Rotation::from_frame(&normal, &x_direction.to_real(), dec)?;
    let norm = first.normal.norm_squared();
    let origin = first.normal.scaled(
        &(&first
            .normal
            .dot(&first.origin.minus(&catallax_shapes::Point3::origin()))
            / &norm),
    );
    let normal = rotation.apply(&Vector3::new(Real::zero(), Real::zero(), Real::one()));
    let pose = Pose::new(
        rotation,
        origin.plus(&normal.scaled(&offset.clone().into())),
        dec,
    )?;
    let placement = crate::profile::placement(placement, dec)?;
    // Avoid algebraic zero terms when a copy does not move the tool.
    if placement.same_motion(&Pose::identity(), dec)? {
        return Ok(pose);
    }
    Ok(pose.then(&placement))
}

/// Face extrusion coordinates are local XY, not another world-space frame.
pub fn check_face_profile(pose: &Pose, dec: &mut Decider) -> Result<(), Refusal> {
    if pose.same_motion(&Pose::identity(), dec)? {
        Ok(())
    } else {
        Err(Refusal::new(
            "extrude.face-profile",
            "a face-mounted extrusion requires a local XY profile",
        ))
    }
}

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
