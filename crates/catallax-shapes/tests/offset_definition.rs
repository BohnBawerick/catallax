#![allow(missing_docs)]

use catallax_numbers::{Decider, Real};
use catallax_shapes::{EllipticCylinder, EllipticOffsetDefinition, OffsetSide, Point3, Vector3};

fn original() -> EllipticCylinder {
    EllipticCylinder {
        origin: Point3::origin(),
        axis: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        major_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        major_sq: Real::from_i64(4),
        minor_sq: Real::one(),
    }
}

#[test]
fn signed_algebraic_offsets_keep_their_branch_and_unsigned_distance() {
    let mut dec = Decider::default();
    let length = Real::from_i64(2).sqrt();
    for (signed, side, magnitude) in [
        (-&length, OffsetSide::Inward, length.clone()),
        (length.clone(), OffsetSide::Outward, length),
        (Real::zero(), OffsetSide::Outward, Real::zero()),
    ] {
        let definition = EllipticOffsetDefinition::new(original(), signed, &mut dec).unwrap();
        assert_eq!(definition.side(), side);
        assert!(
            dec.sign(&(definition.distance() - &magnitude))
                .is_coincident()
        );
    }
}

#[test]
fn invalid_supports_and_undefined_distances_cannot_become_definition_records() {
    let mut dec = Decider::default();
    let mut equal_axes = original();
    equal_axes.minor_sq = equal_axes.major_sq.clone();
    assert_eq!(
        EllipticOffsetDefinition::new(equal_axes, Real::one(), &mut dec)
            .unwrap_err()
            .code,
        "shapes.ellipse-axes"
    );
    let mut invalid_frame = original();
    invalid_frame.major_dir = invalid_frame.axis.clone();
    assert_eq!(
        EllipticOffsetDefinition::new(invalid_frame, Real::one(), &mut dec)
            .unwrap_err()
            .code,
        "shapes.not-perpendicular"
    );
    let undefined = Real::one() / Real::zero();
    assert!(EllipticOffsetDefinition::new(original(), undefined, &mut dec).is_err());
}
