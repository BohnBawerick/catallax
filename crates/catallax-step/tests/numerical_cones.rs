//! A tapered bore exercises the reversed cone face and both periodic seam uses.
#[allow(dead_code)]
#[path = "../../catallax-topology/tests/common/cones.rs"]
mod fixture;

use catallax_numbers::Decider;
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;

#[test]
fn inner_cone_seams_export_from_every_cyclic_loop_start() {
    for start in 0..4 {
        let mut d = fixture::positive_bore();
        d.faces.last_mut().unwrap().loops[0].rotate_left(start);
        let body = n::build(d, n::NumericalOptions::default(), &mut Decider::default()).unwrap();
        let file = output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(
            std::str::from_utf8(file.bytes())
                .unwrap()
                .contains("CONICAL_SURFACE(")
        );
        assert!(
            file.report()
                .checks
                .iter()
                .any(|c| c == "file.cone-step-v-axial-height")
        );
        file.check_written(file.bytes(), &mut Decider::default())
            .unwrap();
        if let Ok(directory) = std::env::var("CATALLAX_CONE_EVIDENCE") {
            std::fs::create_dir_all(&directory).unwrap();
            std::fs::write(
                std::path::Path::new(&directory).join(format!("cone-bore-loop-{start}.step")),
                file.bytes(),
            )
            .unwrap();
        }
    }
}

#[test]
fn cone_angle_error_covers_the_original_enclosure() {
    use catallax_numbers::{Ball, Rational};
    use catallax_shapes::Surface;
    let body = n::build(
        fixture::positive_bore(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let Surface::Cone(cone) = &body.definition().faces.last().unwrap().surface else {
        panic!("cone")
    };
    let slope = Decider::default().enclose_at(&cone.slope, 128).unwrap();
    let angle = Ball::atan2(&slope, &Ball::from_rational(&Rational::one(), 128)).unwrap();
    let file = output::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let text = std::str::from_utf8(file.bytes()).unwrap();
    let fields = text
        .split_once("CONICAL_SURFACE(")
        .unwrap()
        .1
        .split_once(')')
        .unwrap()
        .0;
    let written = Rational::parse(
        fields
            .rsplit(',')
            .next()
            .unwrap()
            .trim()
            .trim_end_matches('.'),
    )
    .unwrap();
    let recorded = file
        .report()
        .checks
        .iter()
        .find_map(|s| s.strip_prefix("file.cone-semi-angle-rounding-max-rad="))
        .unwrap();
    let recorded = Rational::parse(recorded).unwrap();
    let expected = (&written - &angle.lo())
        .abs()
        .max((&written - &angle.hi()).abs());
    assert!(
        recorded >= expected,
        "reported angle error omits part of its source enclosure"
    );
}
