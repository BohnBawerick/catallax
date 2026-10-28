//! P5 housing junction: between its plane-to-edge transitions the junction
//! ball rests on a countersink mouth edge. Both the authored radius 2 and the
//! edited radius 2.5 are covered.

use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Sign};
use catallax_shapes::{
    GeomStore, Point3, Pose, ROLLING_ENVELOPE_TARGET_NM, RollingEnvelope,
    RollingEnvelopeDefinition, Rotation, Side, Surface, Vector3,
};

fn q(text: &str) -> Real {
    Rational::parse(text).unwrap().into()
}
fn v(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}
fn p(x: Real, y: Real, z: Real) -> Point3 {
    Point3::new(x, y, z)
}
fn definition(
    radius: &str,
    body: &str,
    dec: &mut Decider,
) -> Result<RollingEnvelopeDefinition, catallax_numbers::Refusal> {
    RollingEnvelopeDefinition::new(
        p(q("0"), q("5"), q("0")),
        v(0, 1, 0),
        v(1, 0, 0),
        q(body),
        q("4"),
        q("17"),
        q(radius),
        [q("1"), q("1")],
        dec,
    )
}
fn fit(radius: &str, dec: &mut Decider) -> RollingEnvelope {
    RollingEnvelope::fit(definition(radius, "12", dec).unwrap(), dec).unwrap()
}

#[test]
fn p5_envelopes_meet_the_kernel_target() {
    let target = Rational::from_ratio(ROLLING_ENVELOPE_TARGET_NM, 1_000_000).unwrap();
    for radius in ["2", "2.5"] {
        let mut dec = Decider::default();
        let envelope = fit(radius, &mut dec);
        assert!(envelope.bound().upper <= target);
        assert_eq!(envelope.carrier().degrees(), [8, 2]);
        assert_eq!(
            envelope.carrier().control_count(),
            [8 * envelope.span_count() + 1, 3]
        );
        println!(
            "radius {radius}: {} spans, bound {} nm",
            envelope.span_count(),
            (envelope.bound().upper.clone() * Rational::from_i64(1_000_000)).to_f64_nearest()
        );
    }
}

#[test]
fn boundaries_are_exact_on_their_supports() {
    let mut dec = Decider::default();
    let envelope = fit("2", &mut dec);
    let def = envelope.definition().clone();
    for sign in [Sign::Negative, Sign::Positive] {
        for corner in def.corners(sign) {
            assert!(envelope.on_boundary(&corner, &mut dec).unwrap());
            assert_eq!(envelope.side_at(&corner, &mut dec).unwrap(), Sign::Zero);
        }
    }
    // The mouth foot of the middle station and the rail point represented there.
    let mouth_foot = p(q("13"), q("5"), q("0"));
    assert!(envelope.on_boundary(&mouth_foot, &mut dec).unwrap());
    let middle = envelope.span_count() / 2 * 6;
    let rail = envelope.carrier().controls()[middle * 3].clone();
    assert!(envelope.on_boundary(&rail, &mut dec).unwrap());
    // Off the represented rail height, still on the cylinder: not a boundary.
    let above = p(rail.x.clone(), &rail.y + &q("0.001"), rail.z.clone());
    assert!(!envelope.on_boundary(&above, &mut dec).unwrap());
    // A meridian point halfway round its quarter arc, with its exact normal.
    let [foot, _] = def.corners(Sign::Positive);
    let scale = &q("14") / &q("12");
    let centre = p(&foot.x * &scale, q("7"), &foot.z * &scale);
    let half = &q("2") / &Real::from_i64(2).sqrt();
    let radial = centre.minus(&foot).scaled(&(&Real::one() / &q("2")));
    let point = centre
        .plus(&radial.scaled(&(-&half)))
        .plus(&v(0, -1, 0).scaled(&half));
    assert!(envelope.on_boundary(&point, &mut dec).unwrap());
    let normal = envelope.normal_at(&point, &mut dec).unwrap();
    let residual = normal.cross(&centre.minus(&point)).norm_squared();
    assert!(matches!(dec.sign(&residual), Outcome::ProvenCoincident(_)));
    assert!(envelope.normal_at(&mouth_foot, &mut dec).is_err());
}

#[test]
fn side_follows_the_ball_and_refuses_off_the_chart() {
    for radius in ["2", "2.5"] {
        let mut dec = Decider::default();
        let envelope = fit(radius, &mut dec);
        let r = q(radius);
        let height = (&(&q("2") * &r) - &q("1")).sqrt();
        let centre = p(&q("12") + &r, &q("5") + &height, q("0"));
        assert_eq!(envelope.side_at(&centre, &mut dec).unwrap(), Sign::Positive);
        // Added material near the corner, under the face.
        let corner = p(q("12.2"), q("5.2"), q("0"));
        assert_eq!(envelope.side_at(&corner, &mut dec).unwrap(), Sign::Negative);
        let refusal = envelope
            .side_at(&p(q("0"), q("50"), q("0")), &mut dec)
            .unwrap_err();
        assert_eq!(refusal.code, "shapes.bspline-domain");

        let mut geom = GeomStore::new();
        let id = geom.add_surface(Surface::BSpline(envelope.clone().into()));
        assert!(matches!(
            geom.side_of(id, &centre, &mut dec),
            Outcome::Proven(Side::Front)
        ));
        let cover = envelope
            .parameters_of_box(
                &[0, 1, 2].map(|_| Ball::from_rational(&Rational::zero(), 64)),
                &dec,
            )
            .unwrap();
        assert!(cover.is_none());
    }
}

#[test]
fn patterned_copy_keeps_its_certificate() {
    let mut dec = Decider::default();
    let envelope = fit("2", &mut dec);
    // A quarter turn about the housing axis, as the circular hole pattern applies.
    let rotation = Rotation::from_quaternion(
        Rational::one(),
        Rational::zero(),
        Rational::one(),
        Rational::zero(),
    )
    .unwrap();
    let pose = Pose::new(rotation, Vector3::zero(), &mut dec).unwrap();
    let copy = envelope.transformed(&pose);
    assert_eq!(copy.bound().upper, envelope.bound().upper);
    for sign in [Sign::Negative, Sign::Positive] {
        for corner in copy.definition().corners(sign) {
            assert!(copy.on_boundary(&corner, &mut dec).unwrap());
        }
    }
    let r = q("2");
    let height = (&(&q("2") * &r) - &q("1")).sqrt();
    let centre = pose.apply_point(&p(q("14"), &q("5") + &height, q("0")));
    assert_eq!(copy.side_at(&centre, &mut dec).unwrap(), Sign::Positive);
    let mut geom = GeomStore::new();
    let a = geom.add_surface(Surface::BSpline(envelope.into()));
    let b = geom.add_surface(Surface::BSpline(copy.into()));
    assert!(!geom.surface(a).shares_definition_nodes(geom.surface(b)));
}

#[test]
fn other_contact_regimes_refuse() {
    let mut dec = Decider::default();
    for (radius, body, detail) in [
        // Beyond 2 + sqrt(2) the ball would roll onto the countersink cone.
        ("3.5", "12", "second face"),
        // The ball centre path no longer reaches the mouth's near side.
        ("5", "12", "near side"),
        // A mouth that crosses the cylinder is a different junction.
        ("2", "14", "outside the cylinder"),
    ] {
        let refusal = definition(radius, body, &mut dec).unwrap_err();
        assert_eq!(refusal.code, "shapes.rolling-envelope");
        assert!(refusal.detail.contains(detail), "{}", refusal.detail);
    }
}

// Independent float check, not a certificate: every sampled carrier point is
// within the proved bound of the exact envelope, with rows on their supports.
#[test]
fn sampled_carrier_points_lie_on_the_exact_envelope() {
    for (radius, r) in [("2", 2.0f64), ("2.5", 2.5)] {
        let mut dec = Decider::default();
        let envelope = fit(radius, &mut dec);
        let bound = envelope.bound().upper.to_f64_nearest();
        let b = 12.0 + r;
        let spine = |theta: f64| {
            let rho = (b * b - 34.0 * b * theta.cos() + 289.0).sqrt();
            let h = (r * r - (4.0 - rho).powi(2)).max(0.0).sqrt();
            [b * theta.cos(), 5.0 + h, b * theta.sin()]
        };
        let c0 = (b * b + 273.0) / (34.0 * b);
        let limit = c0.acos();
        let distance = |x: [f64; 3]| {
            let g = |t: f64| {
                let c = spine(t);
                let e = 1e-7;
                let (lo, hi) = (spine(t - e), spine(t + e));
                (0..3).map(|i| (x[i] - c[i]) * (hi[i] - lo[i])).sum::<f64>()
            };
            let (mut lo, mut hi) = (-limit, limit);
            assert!(g(lo) * g(hi) <= 0.0, "no foot station");
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if g(lo) * g(mid) <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let c = spine(0.5 * (lo + hi));
            ((0..3).map(|i| (x[i] - c[i]).powi(2)).sum::<f64>().sqrt() - r).abs()
        };
        let mut worst = 0.0f64;
        for i in 1..40 {
            for j in 0..=10 {
                let u = Rational::from_ratio(i - 20, 20).unwrap();
                let w = Rational::from_ratio(j, 10).unwrap();
                let value = envelope
                    .carrier()
                    .evaluate(
                        &Ball::from_rational(&u, 128),
                        &Ball::from_rational(&w, 128),
                        &dec,
                    )
                    .unwrap()
                    .map(|c| c.mid().to_f64_nearest());
                worst = worst.max(distance(value));
                if j == 0 {
                    let radial = (value[0].powi(2) + value[2].powi(2)).sqrt();
                    assert!((radial - 12.0).abs() < 1e-12);
                }
                if j == 10 {
                    let radial = ((value[0] - 17.0).powi(2) + value[2].powi(2)).sqrt();
                    assert!((radial - 4.0).abs() < 1e-12 && (value[1] - 5.0).abs() < 1e-12);
                }
            }
        }
        assert!(
            worst <= bound + 1e-12,
            "sampled {worst} exceeds bound {bound}"
        );
    }
}

#[test]
fn patterned_holes_reuse_one_fit() {
    let mut dec = Decider::default();
    let envelope = fit("2", &mut dec);
    for toward in [v(0, 0, 1), v(-1, 0, 0), v(0, 0, -1)] {
        let definition = RollingEnvelopeDefinition::new(
            p(q("0"), q("5"), q("0")),
            v(0, 1, 0),
            toward,
            q("12"),
            q("4"),
            q("17"),
            q("2"),
            [q("1"), q("1")],
            &mut dec,
        )
        .unwrap();
        let copy = envelope.reframed(definition.clone(), &mut dec).unwrap();
        for sign in [Sign::Negative, Sign::Positive] {
            for corner in definition.corners(sign) {
                assert!(copy.on_boundary(&corner, &mut dec).unwrap());
            }
        }
    }
    let other = definition("2.5", "12", &mut dec).unwrap();
    let refusal = envelope.reframed(other, &mut dec).unwrap_err();
    assert_eq!(refusal.code, "shapes.rolling-envelope");
}
