//! Legacy compatibility and Required torus lifts at the outer equator.

use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_step::{Seams, StepOptions, check_seams, write_step};
use catallax_topology::{Model, block, torus};

fn near_equator_step(seams: Seams) -> Vec<u8> {
    let height = |t: Rational| {
        let square = &t * &t;
        (Rational::from_i64(2) * t) / (Rational::one() + square)
    };
    let low = height(Rational::from_ratio(-4, 10_000_000).unwrap());
    let high = height(Rational::from_ratio(-1, 10_000_000).unwrap());
    let mut model = Model::new();
    let mut dec = Decider::default();
    let ring = torus(&mut model, 4.into(), 1.into(), &Pose::identity()).unwrap();
    let pose = Pose::new(
        Rotation::identity(),
        RVector3::new((-6).into(), (-6).into(), low.clone()).into(),
        &mut dec,
    )
    .unwrap();
    let slab = block(&mut model, [12.into(), 12.into(), high - low], &pose).unwrap();
    let body = catallax_boolean::boolean(
        &mut model,
        ring,
        slab,
        catallax_boolean::BoolOp::Intersect,
        &mut dec,
    )
    .unwrap()
    .body;
    assert!(model.check(body).is_empty());
    let mut bytes = Vec::new();
    write_step(
        &model,
        &[body],
        &StepOptions {
            seams,
            ..StepOptions::default()
        },
        &mut dec,
        &mut bytes,
    )
    .unwrap();
    assert!(check_seams(&bytes).is_ok());
    bytes
}

fn fingerprint(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325_u64, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x100_0000_01b3)
    })
}

#[test]
fn legacy_near_equator_torus_keeps_its_pre_offset_chart_bytes() {
    assert_eq!(
        fingerprint(&near_equator_step(Seams::Legacy)),
        0x4be6_d0cd_c872_36c0
    );
}

#[test]
fn required_near_equator_torus_keeps_its_negative_angle_lift() {
    assert_eq!(
        fingerprint(&near_equator_step(Seams::Required)),
        0x91b6_3cbb_78cd_ba55
    );
}
