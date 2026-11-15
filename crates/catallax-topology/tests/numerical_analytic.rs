//! Retained analytic numerical controls.
#[path = "common/analytic.rs"]
mod fixture;
use catallax_numbers::{Ball, Decider, NumericalCheck};
use catallax_shapes::Surface;
use catallax_topology::{Sense, numerical as n};
use fixture::{q, quarter};

#[test]
fn finite_circular_caps_and_cylinder_keep_numerical_checks() {
    let body = n::build(
        quarter(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let count = body.counts();
    assert_eq!((count.vertices, count.edges, count.faces), (6, 9, 5));
    assert!(body.evidence().checks().iter().any(|check| matches!(check,NumericalCheck::Passed { name,.. } if name=="trim.sampled-simple-oriented-parameter-loop")));
    assert!(
        body.evidence()
            .checks()
            .iter()
            .filter_map(|check| match check {
                NumericalCheck::Samples { samples, .. } => Some(samples),
                _ => None,
            })
            .flatten()
            .all(|sample| sample.distance_upper_mm <= q("0.01"))
    );
    assert!(!body.exact_topology_faults().is_empty());
}

#[test]
fn a_cylinder_residual_over_the_source_tolerance_refuses() {
    for (radius, pass) in [("1.0099", true), ("1.0101", false)] {
        let mut source = quarter();
        let Surface::Cylinder(cylinder) = &mut source.faces[4].surface else {
            unreachable!()
        };
        cylinder.radius_sq = (&q(radius) * &q(radius)).into();
        let result = n::build(
            source,
            n::NumericalOptions::default(),
            &mut Decider::default(),
        );
        assert_eq!(result.is_ok(), pass);
        if !pass {
            assert_eq!(
                result.unwrap_err().code,
                "numbers.numerical-tolerance-exceeded"
            );
        }
    }
}

#[test]
fn reversed_cylinder_winding_and_unresolved_periodic_cover_refuse() {
    let mut source = quarter();
    source.faces[4].sense = Sense::Reversed;
    assert!(
        n::build(
            source,
            n::NumericalOptions::default(),
            &mut Decider::default()
        )
        .is_err()
    );
    let mut source = quarter();
    source.faces[4].loops[0][0].uv[1][0] = Ball::pi(128).hi() * q("2");
    assert!(
        n::build(
            source,
            n::NumericalOptions::default(),
            &mut Decider::default()
        )
        .is_err()
    );
}
