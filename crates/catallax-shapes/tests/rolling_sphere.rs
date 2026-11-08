//! The fitted ball surface between a sphere and a cylinder that passes through
//! it off its centre.

use std::sync::OnceLock;

use catallax_numbers::{Ball, BudgetPolicy, Decider, Outcome, Rational, Real, Sign, StrictSign};
use catallax_shapes::{
    BallSide, Cylinder, Point3, Pose, RollingSphere, RollingSphereDefinition, Rotation, Sphere,
    SplineBoundary, Surface, Vector3,
};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn real(value: &str) -> Real {
    Real::from(q(value))
}
fn point(v: [&str; 3]) -> Point3 {
    let [x, y, z] = v.map(real);
    Point3::new(x, y, z)
}
fn vector(v: [&str; 3]) -> Vector3 {
    let [x, y, z] = v.map(real);
    Vector3::new(x, y, z)
}

// The corpus case: a radius-5 sphere at the origin and a radius-2 cylinder on
// the axis x = 1, y = 0, joined, with a radius-1/2 concave round.
fn sphere() -> Sphere {
    Sphere {
        center: point(["0", "0", "0"]),
        radius: real("5"),
        axis: vector(["0", "0", "1"]),
        ref_dir: vector(["1", "0", "0"]),
    }
}
fn cylinder() -> Cylinder {
    Cylinder {
        origin: point(["1", "0", "-10"]),
        axis: vector(["0", "0", "1"]),
        radius_sq: real("4"),
        ref_dir: vector(["1", "0", "0"]),
    }
}
// A point of each ring: x = 3, y = 0 lies on the cylinder, and z^2 = 25 - 9.
fn witness(branch: i64) -> Point3 {
    point(["3", "0", if branch > 0 { "4" } else { "-4" }])
}

fn definition(branch: i64) -> RollingSphereDefinition {
    RollingSphereDefinition::new(
        &sphere(),
        &cylinder(),
        &real("1/2"),
        [BallSide::Outside, BallSide::Outside],
        &witness(branch),
        &mut decider(),
    )
    .unwrap()
}

fn fitted(branch: i64) -> RollingSphere {
    static FITS: OnceLock<[RollingSphere; 2]> = OnceLock::new();
    let fits = FITS.get_or_init(|| {
        [1, -1].map(|branch| RollingSphere::fit(definition(branch), &mut decider()).unwrap())
    });
    fits[usize::from(branch < 0)].clone()
}

fn zero(value: &Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(value), Outcome::ProvenCoincident(_))
}
fn on_sphere(point: &Point3, radius_sq: &Real, dec: &mut Decider) -> bool {
    zero(
        &(&point.minus(&sphere().center).norm_squared() - radius_sq),
        dec,
    )
}
fn on_cylinder(point: &Point3, radius_sq: &Real, dec: &mut Decider) -> bool {
    let c = cylinder();
    let flat = point.minus(&c.origin).cross(&c.axis).norm_squared();
    zero(&(&flat - radius_sq), dec)
}

fn chart_samples() -> Vec<Rational> {
    (0..=40)
        .map(|i| &q("-1") + &(q("1/10") * Rational::from_i64(i)))
        .collect()
}

#[test]
fn the_law_rests_the_ball_on_both_supports_exactly() {
    let mut dec = decider();
    for branch in [1, -1] {
        let definition = definition(branch);
        assert_eq!(definition.branch(), branch);
        for u in chart_samples().iter().step_by(5) {
            let [centre, on_ball_sphere, on_ball_cylinder] = definition.section(u).unwrap();
            assert!(on_sphere(&centre, &real("11/2").square(), &mut dec));
            assert!(on_cylinder(&centre, &real("5/2").square(), &mut dec));
            assert!(on_sphere(&on_ball_sphere, &real("25"), &mut dec));
            assert!(on_cylinder(&on_ball_cylinder, &real("4"), &mut dec));
            let ball = real("1/4");
            for contact in [&on_ball_sphere, &on_ball_cylinder] {
                assert!(zero(
                    &(&contact.minus(&centre).norm_squared() - &ball),
                    &mut dec
                ));
            }
            let height = &centre.z * &Real::from(Rational::from_i64(branch));
            assert!(matches!(
                dec.sign(&height),
                Outcome::Proven(StrictSign::Positive)
            ));
        }
    }
}

#[test]
fn the_rails_lie_exactly_on_their_supports_and_close() {
    let mut dec = decider();
    for branch in [1, -1] {
        let fitted = fitted(branch);
        let [on_ball_sphere, on_ball_cylinder] = fitted.rails();
        assert!(on_ball_sphere.is_boundary_of(fitted.carrier(), SplineBoundary::VMin));
        assert!(on_ball_cylinder.is_boundary_of(fitted.carrier(), SplineBoundary::VMax));
        for u in chart_samples().iter().step_by(3) {
            let u = Real::from(u.clone());
            let at = on_ball_sphere.point_at(&u, &mut dec).unwrap();
            assert!(
                on_sphere(&at, &real("25"), &mut dec),
                "sphere rail at {u:?}"
            );
            let at = on_ball_cylinder.point_at(&u, &mut dec).unwrap();
            assert!(
                on_cylinder(&at, &real("4"), &mut dec),
                "cylinder rail at {u:?}"
            );
        }
        let surface = fitted.carrier();
        let [columns, rows] = surface.control_count();
        for row in 0..rows {
            let last = (columns - 1) * rows + row;
            assert!(zero(
                &(&surface.weights()[row] - &surface.weights()[last]),
                &mut dec
            ));
            let gap = surface.controls()[row].minus(&surface.controls()[last]);
            assert!(zero(&gap.norm_squared(), &mut dec));
        }
    }
}

// Each carrier point sits at the ball's radius from the exact centre at the
// same chart parameter, within the claimed bound, and the bound meets the
// 100 nm target.
#[test]
fn the_carrier_stays_within_its_bound_of_the_rolling_ball() {
    let dec = decider();
    for branch in [1, -1] {
        let fitted = fitted(branch);
        let bound = fitted.bound().upper.to_f64_nearest();
        assert!(bound <= 1e-4, "bound {bound}");
        let definition = definition(branch);
        let ball = |x: &Rational| Ball::from_rational(x, 128);
        for u in chart_samples() {
            let [centre, ..] = definition.section(&u).unwrap();
            let centre = [&centre.x, &centre.y, &centre.z]
                .map(|value| dec.enclose_at(value, 128).unwrap().mid().to_f64_nearest());
            for v in ["0", "1/4", "1/2", "3/4", "1"] {
                let at = fitted
                    .carrier()
                    .evaluate(&ball(&u), &ball(&q(v)), &dec)
                    .unwrap()
                    .map(|x| x.mid().to_f64_nearest());
                let distance = ((at[0] - centre[0]).powi(2)
                    + (at[1] - centre[1]).powi(2)
                    + (at[2] - centre[2]).powi(2))
                .sqrt();
                assert!(
                    (distance - 0.5).abs() <= bound + 1e-9,
                    "u {u:?} v {v}: {distance}"
                );
            }
        }
    }
}

#[test]
fn a_moved_copy_keeps_its_rails_on_the_moved_supports() {
    let fitted = fitted(1);
    let mut dec = decider();
    let rotation = Rotation::from_quaternion(q("1"), q("2"), q("3"), q("4")).unwrap();
    let pose = Pose::new(rotation, vector(["7/3", "-11/7", "13/5"]), &mut dec).unwrap();
    let moved = fitted.transformed(&pose);
    assert_eq!(moved.span_count(), fitted.span_count());
    assert!(moved.shares_definition_nodes(&fitted));
    let Surface::Sphere(round) = Surface::Sphere(sphere()).transformed(&pose) else {
        unreachable!()
    };
    let Surface::Cylinder(tube) = Surface::Cylinder(cylinder()).transformed(&pose) else {
        unreachable!()
    };
    let [on_round, on_tube] = moved.rails();
    for u in ["-1", "1/3", "2"] {
        let at = on_round.point_at(&real(u), &mut dec).unwrap();
        let gap = &at.minus(&round.center).norm_squared() - &real("25");
        assert!(zero(&gap, &mut dec), "{u}");
        let at = on_tube.point_at(&real(u), &mut dec).unwrap();
        let flat = at.minus(&tube.origin).cross(&tube.axis).norm_squared();
        assert!(zero(&(&flat - &real("4")), &mut dec), "{u}");
    }
}

#[test]
fn points_on_the_ball_side_and_on_the_carrier_classify() {
    let fitted = fitted(1);
    let mut dec = decider();
    let definition = definition(1);
    let [centre, ..] = definition.section(&q("1/3")).unwrap();
    // The fitted rails lie exactly on the supports; the exact contacts lie
    // within the bound of them, inside the carrier's enclosure width.
    for rail in fitted.rails() {
        let at = rail.point_at(&real("1/3"), &mut dec).unwrap();
        assert_eq!(fitted.side_at(&at, &mut dec).unwrap(), Sign::Zero);
    }
    assert_eq!(fitted.side_at(&centre, &mut dec).unwrap(), Sign::Positive);
}

#[test]
fn a_cylinder_through_the_centre_or_too_wide_is_refused() {
    let mut dec = decider();
    let mut coaxial = cylinder();
    coaxial.origin = point(["0", "0", "-10"]);
    assert!(
        RollingSphereDefinition::new(
            &sphere(),
            &coaxial,
            &real("1/2"),
            [BallSide::Outside, BallSide::Outside],
            &point(["2", "0", "4"]),
            &mut dec,
        )
        .is_err()
    );
    let mut wide = cylinder();
    wide.radius_sq = real("16");
    assert!(
        RollingSphereDefinition::new(
            &sphere(),
            &wide,
            &real("1/2"),
            [BallSide::Outside, BallSide::Outside],
            &witness(1),
            &mut dec,
        )
        .is_err()
    );
}

// Each table entry is the chart angle and height of the rail's own point at
// that chart parameter, read here with plain trigonometry.
#[test]
fn the_rail_tables_follow_the_rails_on_their_charts() {
    let mut dec = decider();
    let tolerance = q("1/100000000000");
    for branch in [1, -1] {
        let fitted = fitted(branch);
        for rail in 0..2 {
            let table = fitted.parameter_table(rail, &tolerance, &dec).unwrap();
            assert!(table.bound < q("1/1000000000"), "{:?}", table.bound);
            assert!(table.height_bound < q("1/1000000000"));
            let curve = &fitted.rails()[rail];
            let mut last = f64::NEG_INFINITY;
            for (index, knot) in table.knots.iter().enumerate().step_by(7) {
                let at = curve.point_at(&Real::from(knot.clone()), &mut dec).unwrap();
                let [x, y, z] = [&at.x, &at.y, &at.z]
                    .map(|value| dec.enclose_at(value, 96).unwrap().mid().to_f64_nearest());
                // The sphere's centre is the origin and the axis is +z.
                let x = x - if rail == 0 { 0.0 } else { 1.0 };
                let angle = table.angle[index][0].to_f64_nearest();
                let turned = (y.atan2(x) - angle).rem_euclid(std::f64::consts::TAU);
                assert!(
                    turned.min(std::f64::consts::TAU - turned) < 1e-9,
                    "{rail} {index}"
                );
                let height = table.height[index][0].to_f64_nearest();
                let expected = if rail == 0 { (z / 5.0).asin() } else { z };
                assert!(
                    (height - expected).abs() < 1e-9,
                    "{rail} {index}: {height} {expected}"
                );
                assert!(angle > last);
                last = angle;
            }
            let first = table.angle[0][0].to_f64_nearest();
            let end = table.angle.last().unwrap()[0].to_f64_nearest();
            assert!((end - first - std::f64::consts::TAU).abs() < 1e-9);
        }
    }
}

// A vertical line at x = 31/10 in the plane y = 0 crosses the upper band once,
// on the ball about (7/2, 0, sqrt(18)), where its height is sqrt(18) - 3/10.
#[test]
fn a_line_meets_the_band_once_on_the_ball() {
    let mut dec = decider();
    let upper = fitted(1);
    let hits = upper
        .line_hits(
            &point(["31/10", "0", "-20"]),
            &vector(["0", "0", "1"]),
            &mut dec,
        )
        .unwrap();
    assert_eq!(hits.len(), 1);
    let z = dec
        .enclose_at(&hits[0].z, 96)
        .unwrap()
        .mid()
        .to_f64_nearest();
    let bound = upper.bound().upper.to_f64_nearest();
    assert!(
        (z - (18f64.sqrt() - 0.3)).abs() <= 4.0 * bound + 1e-12,
        "{z}"
    );
    // The lower ring's band is this one mirrored in z = 0.
    let lower = fitted(-1)
        .line_hits(
            &point(["31/10", "0", "-20"]),
            &vector(["0", "0", "1"]),
            &mut dec,
        )
        .unwrap();
    assert_eq!(lower.len(), 1);
    let z = dec
        .enclose_at(&lower[0].z, 96)
        .unwrap()
        .mid()
        .to_f64_nearest();
    assert!(
        (z + (18f64.sqrt() - 0.3)).abs() <= 4.0 * bound + 1e-12,
        "{z}"
    );
}
