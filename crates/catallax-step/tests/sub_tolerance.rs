//! Generated straight seams obey the same reader floor as model edges.
use catallax_numbers::{Decider, Rational};
use catallax_shapes::Pose;
use catallax_step::{StepOptions, write_step};
use catallax_topology::{Model, cylinder};

#[test]
fn a_thin_cylinder_refuses_its_generated_line_seam_atomically() {
    let mut model = Model::new();
    let body = cylinder(
        &mut model,
        Rational::one(),
        Rational::parse("3.5e-10").unwrap(),
        &Pose::identity(),
    )
    .unwrap();
    let original = model.counts(body);
    assert!(model.check(body).is_empty());
    let before = b"existing output".to_vec();
    let mut bytes = before.clone();
    let refusal = write_step(
        &model,
        &[body],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap_err();
    assert_eq!(refusal.code, "step.sub-tolerance-edge");
    assert_eq!(bytes, before);
    assert_eq!(model.counts(body), original);
    assert!(model.check(body).is_empty());
}
