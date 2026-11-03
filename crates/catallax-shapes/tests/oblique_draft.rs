//! Independent neutral-generator checks for the exact oblique draft carrier.

use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Sign};
use catallax_shapes::{Cylinder, ObliqueDraft, Plane, Point3, Vector3};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn r(s: &str) -> Real {
    q(s).into()
}
fn p(x: &str, y: &str, z: &str) -> Point3 {
    Point3::new(r(x), r(y), r(z))
}
fn v(x: &str, y: &str, z: &str) -> Vector3 {
    Vector3::new(r(x), r(y), r(z))
}
fn source() -> Cylinder {
    Cylinder {
        origin: p("0", "0", "0"),
        axis: v("0", "0", "2"),
        radius_sq: r("64"),
        ref_dir: v("3", "0", "0"),
    }
}
fn neutral() -> Plane {
    Plane {
        origin: p("0", "0", "1"),
        normal: v("1", "2", "8"),
        ref_dir: v("2", "-1", "0"),
    }
}
fn support(tangent: &str) -> ObliqueDraft {
    ObliqueDraft::new(
        &source(),
        v("0", "0", "1"),
        &neutral(),
        r(tangent),
        [r("0"), r("2")],
        &mut Decider::default(),
    )
    .unwrap()
}
fn zero(value: &Real) {
    assert!(
        matches!(Decider::default().sign(value), Outcome::ProvenCoincident(_)),
        "{value:?}"
    );
}

#[test]
fn carrier_equals_neutral_generators_in_every_span_and_at_both_caps() {
    for tangent in ["1", "-1", "1/3"] {
        let draft = support(tangent);
        let tangent = q(tangent);
        for quarter in 0..4 {
            for tau in ["0", "1/7", "1/2", "1"] {
                let tau = q(tau);
                let denominator = Rational::one() + &tau * &tau;
                let c = (Rational::one() - &tau * &tau) / denominator.clone();
                let s = (&tau * &Rational::from_i64(2)) / denominator;
                let [c, s] = match quarter {
                    0 => [c, s],
                    1 => [-s, c],
                    2 => [-c, -s],
                    _ => [s, -c],
                };
                // Independently solve the authored neutral plane on the input
                // radius-8 cylinder, then advance its specified generator.
                let neutral_z = &(&Rational::one() - &c) - &(&Rational::from_i64(2) * &s);
                for z in ["0", "2/3", "2"] {
                    let z = q(z);
                    let axial_displacement = &z - &neutral_z;
                    let radius = Rational::from_i64(8) - &tangent * &axial_displacement;
                    let expected = [&radius * &c, &radius * &s, z.clone()];
                    let u = &Rational::from_i64(quarter) + &tau;
                    let actual = draft
                        .carrier()
                        .evaluate(
                            &Ball::from_rational(&u, 192),
                            &Ball::from_rational(&(&z / &Rational::from_i64(2)), 192),
                            &Decider::default(),
                        )
                        .unwrap();
                    for (got, want) in actual.iter().zip(expected) {
                        assert!(
                            got.lo() <= want && want <= got.hi(),
                            "quarter {quarter}, tau {tau}, z {z}: {got:?}, {want}"
                        );
                        assert!(got.hi() - got.lo() < q("1e-40"));
                    }
                    let at = Point3::new(
                        Real::from(&radius * &c),
                        Real::from(&radius * &s),
                        Real::from(z),
                    );
                    assert_eq!(
                        draft.side_at(&at, &mut Decider::default()).unwrap(),
                        Sign::Zero
                    );
                    let normal = draft.normal_at(&at, &mut Decider::default()).unwrap();
                    let generator = Vector3::new(
                        Real::from(-(&tangent * &c)),
                        Real::from(-(&tangent * &s)),
                        Real::one(),
                    );
                    zero(&normal.dot(&generator));
                }
            }
        }
    }
}

#[test]
fn unsquared_branch_excludes_the_spurious_axis_and_line_roots_beyond_caps() {
    let draft = support("1");
    assert_eq!(
        draft
            .side_at(&p("0", "0", "1"), &mut Decider::default())
            .unwrap(),
        Sign::Negative
    );
    for (point, side) in [
        (p("6", "0", "1"), Sign::Negative),
        (p("7", "0", "1"), Sign::Zero),
        (p("8", "0", "1"), Sign::Positive),
        (p("-9", "0", "1"), Sign::Zero),
    ] {
        assert_eq!(
            draft.side_at(&point, &mut Decider::default()).unwrap(),
            side
        );
    }
    let hits = draft
        .line_hits(
            &p("0", "0", "1"),
            &v("1", "0", "0"),
            &mut Decider::default(),
        )
        .unwrap();
    assert_eq!(hits.len(), 2);
    zero(&(&hits[0].x - &r("-9")));
    zero(&(&hits[1].x - &r("7")));
    assert!(
        draft
            .line_hits(
                &p("0", "0", "3"),
                &v("1", "0", "0"),
                &mut Decider::default()
            )
            .unwrap()
            .is_empty()
    );
    assert!(
        draft
            .line_hits(
                &p("0", "20", "1"),
                &v("1", "0", "0"),
                &mut Decider::default()
            )
            .unwrap()
            .is_empty()
    );
}

#[test]
fn section_integral_contains_the_independent_fourier_volume() {
    // rho(z,u) = 9-z-cos(u)-2sin(u). Integrating rho^2/2 first
    // around u and then from z=0 to 2 gives pi*(401/3).
    let draft = support("1");
    zero(&(&draft.volume_over_pi() - &r("401/3")));
}

#[test]
fn unproved_folds_nontransverse_neutrals_and_invalid_frames_stay_typed() {
    let mut plane = neutral();
    plane.normal = v("1", "0", "0");
    let error = ObliqueDraft::new(
        &source(),
        v("0", "0", "1"),
        &plane,
        r("1"),
        [r("0"), r("2")],
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "draft.unsupported-oblique-arrangement");
    let error = ObliqueDraft::new(
        &source(),
        v("0", "0", "1"),
        &neutral(),
        r("4"),
        [r("0"), r("2")],
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "draft.unsupported-oblique-arrangement");
    let mut cylinder = source();
    cylinder.ref_dir = Vector3::zero();
    let error = ObliqueDraft::new(
        &cylinder,
        v("0", "0", "1"),
        &neutral(),
        r("1"),
        [r("0"), r("2")],
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "shapes.zero-vector");
}
