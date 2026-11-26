//! Legacy compatibility and Required torus lifts at the outer equator.

mod common;

use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_step::{Seams, StepOptions, check_seams, write_step};
use catallax_topology::{Model, block, torus};
use common::{Step, list, number, r, reference};

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

fn torus_angles(seams: Seams) -> Vec<Rational> {
    let step = Step::parse(
        &near_equator_step(seams),
        &StepOptions::default().product_name,
    );
    step.ids("PCURVE")
        .into_iter()
        .filter_map(|id| {
            let pcurve = step.e(id, "PCURVE");
            let surface = &step.entities[&reference(&pcurve.args[1])];
            (surface.kind == "TOROIDAL_SURFACE").then(|| reference(&pcurve.args[2]))
        })
        .flat_map(|id| list(&step.e(id, "DEFINITIONAL_REPRESENTATION").args[1]))
        .flat_map(|id| list(&step.e(reference(id), "B_SPLINE_CURVE_WITH_KNOTS").args[2]))
        .map(|id| number(list(&step.e(reference(id), "CARTESIAN_POINT").args[1])[1]))
        .collect()
}

fn check_lift(seams: Seams, turn: Rational) {
    let angles = torus_angles(seams);
    // The slab cuts both walls of the torus. The inner wall is near pi
    // under either lift; only the outer wall crosses the chart cut.
    let pi = r("3.141592653589793");
    let outer = || angles.iter().filter(|v| (*v - &pi).abs() > Rational::one());
    let low = outer().min().expect("outer torus pcurves");
    let high = outer().max().unwrap();
    let rounding = r("1e-9");
    assert!((low - &(&turn - &r("8e-7"))).abs() < rounding, "{low}");
    assert!((high - &(&turn - &r("2e-7"))).abs() < rounding, "{high}");
}

#[test]
fn legacy_near_equator_torus_keeps_its_positive_angle_lift() {
    check_lift(Seams::Legacy, r("6.283185307179586"));
}

#[test]
fn required_near_equator_torus_keeps_its_negative_angle_lift() {
    check_lift(Seams::Required, Rational::zero());
}
