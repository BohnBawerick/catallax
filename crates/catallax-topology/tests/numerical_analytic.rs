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
    // The stored tolerance is 1.1 times the measured gap (rule T4), so a gap
    // of 0.009 mm stores 0.0099 mm and one of 0.0092 mm stores 0.01012 mm.
    for (radius, pass) in [("1.009", true), ("1.0092", false)] {
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

/// A gap whose stored tolerance lies within the binary64 rounding error of
/// the limit is decided at 128 bits (decision Q6, section 4.2 step 7). A gap
/// of 1/110 mm stores exactly the 0.01 mm limit; 1e-16 mm to either side is
/// far below what binary64 resolves there and far above what 128 bits do.
#[test]
fn a_gap_within_rounding_of_the_limit_is_decided_at_higher_precision() {
    let at_limit = &q("1") + &q("1/110");
    for (radius, pass) in [
        (&at_limit - &q("1e-16"), true),
        (&at_limit + &q("1e-16"), false),
    ] {
        let mut source = quarter();
        let Surface::Cylinder(cylinder) = &mut source.faces[4].surface else {
            unreachable!()
        };
        cylinder.radius_sq = (&radius * &radius).into();
        let result = n::build(
            source,
            n::NumericalOptions::default(),
            &mut Decider::default(),
        );
        assert_eq!(result.is_ok(), pass, "{radius}");
        let body = match result {
            Ok(body) => body,
            Err(refusal) => {
                assert_eq!(refusal.code, "numbers.numerical-tolerance-exceeded");
                continue;
            }
        };
        // The largest stored tolerance is the 128-bit value: under the limit
        // by about 1.1e-16 mm, closer than a binary64 measurement could say.
        let largest = body
            .evidence()
            .checks()
            .iter()
            .filter_map(|check| match check {
                NumericalCheck::Samples { name, samples, .. }
                    if name == "trim.edge-surface-distance" =>
                {
                    Some(samples)
                }
                _ => None,
            })
            .flatten()
            .map(|sample| sample.distance_upper_mm.clone())
            .max()
            .expect("samples");
        assert!(largest <= q("0.01"));
        assert!(&q("0.01") - &largest < q("2e-16"), "{largest}");
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
