use super::*;
use catallax_boolean::{BoolOp, boolean};
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Rational, with_operation_budget};
use catallax_shapes::{Pose, Rotation};
use catallax_topology::{block, cone, cylinder, sphere, torus};

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn contours(model: &Model, body: BodyId, d: RVector3) -> Contours {
    contour_generators(model, body, &d, &mut Decider::default())
        .unwrap()
        .value
}

fn is_zero(value: &Real) -> bool {
    sign(value, &mut Decider::default()).unwrap() == 0
}

fn curved_faces(model: &Model, body: BodyId) -> Vec<FaceId> {
    model
        .faces_of_body(body)
        .into_iter()
        .filter(|&face| {
            !matches!(
                model.geom().surface(model.face(face).surface),
                Surface::Plane(_)
            )
        })
        .collect()
}

/// Ends of a line arc, as exact points.
fn ends(c: &Contours, arc: &ContourArc) -> [Point3; 2] {
    let CurveForm::Line { .. } =
        curve_form(c.model.geom(), arc.curve, &mut Decider::default()).unwrap()
    else {
        panic!("expected a line generator");
    };
    [arc.start.unwrap(), arc.end.unwrap()].map(|p| c.model.geom().position(p).clone())
}

/// The surface normal at an exact point of a quadric face, unnormalised.
fn normal_at(surface: &Surface, p: &Point3) -> Vector3 {
    match surface {
        Surface::Cylinder(s) => {
            let r = p.minus(&s.origin);
            r.scaled(&s.axis.norm_squared())
                .minus(&s.axis.scaled(&r.dot(&s.axis)))
        }
        Surface::Cone(s) => {
            let r = p.minus(&s.apex);
            let steep = Real::one() + s.slope.square();
            r.scaled(&s.axis.norm_squared())
                .minus(&s.axis.scaled(&(&steep * &r.dot(&s.axis))))
        }
        Surface::Sphere(s) => p.minus(&s.center),
        other => panic!("no test normal for {:?}", other.kind()),
    }
}

/// Every end of every arc lies on its face's carrier with a normal
/// perpendicular to the view direction.
fn assert_generators(model: &Model, c: &Contours, d: &RVector3) {
    for arc in &c.arcs {
        let surface = model.geom().surface(model.face(arc.face).surface);
        for end in ends(c, arc) {
            assert!(is_zero(&normal_at(surface, &end).dot(&d.to_real())));
        }
    }
}

#[test]
fn cylinder_generators_are_two_axis_lines_and_none_along_the_axis() {
    let mut model = Model::new();
    let body = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let side = curved_faces(&model, body);
    assert_eq!(side.len(), 1);
    let front = contours(&model, body, vector(0, -1, 0));
    assert_eq!(front.arcs.len(), 2);
    let mut xs = Vec::new();
    for arc in &front.arcs {
        assert_eq!(vec![arc.face], side);
        assert_eq!(arc.on_faces, side);
        let [a, b] = ends(&front, arc);
        assert!(is_zero(&a.y) && is_zero(&b.y));
        assert!(is_zero(&(&a.x - &b.x)));
        let mut z = [a.z.as_rational().unwrap(), b.z.as_rational().unwrap()];
        z.sort();
        assert_eq!(z, [Rational::zero(), Rational::from(5)]);
        xs.push(a.x.as_rational().unwrap());
    }
    xs.sort();
    assert_eq!(xs, [Rational::from(-2), Rational::from(2)]);
    assert!(contours(&model, body, vector(0, 0, 1)).arcs.is_empty());
    assert!(contours(&model, body, vector(0, 0, -3)).arcs.is_empty());
    let d = vector(1, 2, 3);
    let oblique = contours(&model, body, d.clone());
    assert_eq!(oblique.arcs.len(), 2);
    assert_generators(&model, &oblique, &d);
}

#[test]
fn cone_generators_meet_the_apex_plane_and_vanish_inside_the_cone_angle() {
    let mut model = Model::new();
    // Radius 3 at z = 0 and 1 at z = 4: slope 1/2, apex at z = 6.
    let body = cone(&mut model, 3.into(), 1.into(), 4.into(), &Pose::identity()).unwrap();
    let front = contours(&model, body, vector(0, -1, 0));
    assert_eq!(front.arcs.len(), 2);
    for arc in &front.arcs {
        let [a, b] = ends(&front, arc);
        let (low, high) = if a.z.as_rational().unwrap().is_zero() {
            (a, b)
        } else {
            (b, a)
        };
        assert!(is_zero(&low.y) && is_zero(&high.y) && is_zero(&low.z));
        assert_eq!(high.z.as_rational().unwrap(), Rational::from(4));
        assert_eq!(low.x.as_rational().unwrap().abs(), Rational::from(3));
        assert_eq!(high.x.as_rational().unwrap().abs(), Rational::one());
    }
    // Along the axis, and inside the half angle (|d_perp| = 1 < 10 / 2), the
    // polar plane meets the nappe only at the apex, which is off the frustum.
    assert!(contours(&model, body, vector(0, 0, 1)).arcs.is_empty());
    assert!(contours(&model, body, vector(0, 1, 10)).arcs.is_empty());
    let d = vector(0, 1, 1);
    let oblique = contours(&model, body, d.clone());
    assert_eq!(oblique.arcs.len(), 2);
    assert_generators(&model, &oblique, &d);
}

#[test]
fn sphere_generator_is_the_whole_great_circle_facing_the_viewer() {
    let mut model = Model::new();
    let body = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let d = vector(1, 2, 2);
    let answer = contours(&model, body, d.clone());
    assert_eq!(answer.arcs.len(), 1);
    let arc = &answer.arcs[0];
    assert_eq!((arc.start, arc.end), (None, None));
    let CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ..
    } = curve_form(answer.model.geom(), arc.curve, &mut Decider::default()).unwrap()
    else {
        panic!("expected a circle");
    };
    assert!(is_zero(&center.minus(&Point3::origin()).norm_squared()));
    assert!(is_zero(&normal.cross(&d.to_real()).norm_squared()));
    assert!(is_zero(&(&radius_sq - &Real::from(25))));
}

#[test]
fn a_generator_on_a_coplanar_face_runs_along_its_edges() {
    // A half cylinder: its flat face lies in the front view's polar plane.
    let mut model = Model::new();
    let round = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let cutter = block(
        &mut model,
        [10.into(), 10.into(), 10.into()],
        &Pose::from_rational(Rotation::identity(), vector(-5, -10, -1)),
    )
    .unwrap();
    let half = boolean(
        &mut model,
        round,
        cutter,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let side = curved_faces(&model, half);
    assert_eq!(side.len(), 1);
    let front = contours(&model, half, vector(0, 1, 0));
    assert_eq!(front.arcs.len(), 2);
    for arc in &front.arcs {
        assert_eq!(vec![arc.face], side);
        assert_eq!(arc.on_faces.len(), 2);
        let [a, b] = ends(&front, arc);
        assert!(is_zero(&a.y) && is_zero(&b.y));
        assert_eq!(a.x.as_rational().unwrap().abs(), Rational::from(2));
    }
}

#[test]
fn planar_bodies_have_no_generators_and_spend_no_work() {
    let mut model = Model::new();
    let body = block(
        &mut model,
        [1.into(), 2.into(), 3.into()],
        &Pose::identity(),
    )
    .unwrap();
    let (answer, work) = with_operation_budget(1_000_000, || {
        contour_generators(&model, body, &vector(1, 1, 1), &mut Decider::default())
    });
    assert!(answer.unwrap().value.arcs.is_empty());
    assert!(!work.exhausted);
}

#[test]
fn zero_directions_and_exhausted_allowances_refuse_with_counts() {
    let mut model = Model::new();
    let round = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let zero =
        contour_generators(&model, round, &vector(0, 0, 0), &mut Decider::default()).unwrap_err();
    assert_eq!(zero.refusal.code, "app.contour-direction");
    let (starved, work) = with_operation_budget(10, || {
        contour_generators(&model, round, &vector(0, 1, 0), &mut Decider::default())
    });
    assert!(starved.is_err());
    assert!(work.exhausted);
}

/// A circle carrier's centre, normal and squared radius.
fn circle_of(c: &Contours, arc: &ContourArc) -> (Point3, Vector3, Real) {
    let CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ..
    } = curve_form(c.model.geom(), arc.curve, &mut Decider::default()).unwrap()
    else {
        panic!("expected a circle");
    };
    (center, normal, radius_sq)
}

fn ring() -> (Model, BodyId) {
    let mut model = Model::new();
    let body = torus(&mut model, 5.into(), 1.into(), &Pose::identity()).unwrap();
    (model, body)
}

/// The ring with everything inside a block removed.
fn ring_minus(corner: RVector3, size: [i64; 3]) -> (Model, BodyId) {
    let (mut model, body) = ring();
    let cutter = block(
        &mut model,
        size.map(Rational::from),
        &Pose::from_rational(Rotation::identity(), corner),
    )
    .unwrap();
    let cut = boolean(
        &mut model,
        body,
        cutter,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    (model, cut)
}

/// The true generator of the ring (centre at the origin, axis z, seam x,
/// R = 5, r = 1) at azimuth `u` on branch `s`, from its normal
/// `cos v e(u) + sin v z` and `N . d = 0`.
fn true_generator(d: [f64; 3], s: f64, u: f64) -> [f64; 3] {
    let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let d = d.map(|x| x / n);
    let e = [u.cos(), u.sin()];
    let (alpha, beta) = (d[2], e[0] * d[0] + e[1] * d[1]);
    let root = (alpha * alpha + beta * beta).sqrt();
    let (cos_v, sin_v) = (s * alpha / root, -s * beta / root);
    [(5.0 + cos_v) * e[0], (5.0 + cos_v) * e[1], sin_v]
}

/// Samples a fitted ring arc exactly at `count` parameters and checks each
/// sample against the true generator at the same azimuth. Returns them.
fn fitted_samples(c: &Contours, arc: &ContourArc, d: [f64; 3], count: usize) -> Vec<[f64; 3]> {
    let fit = arc.fit.as_ref().expect("a fitted arc");
    assert!(fit.deviation_mm <= fit.target_mm);
    assert_eq!(
        fit.tolerance_mm
            .partial_cmp(&(1.1 * fit.deviation_mm).next_up()),
        Some(std::cmp::Ordering::Equal)
    );
    let Curve::BSpline(curve) = c.model.geom().curve(arc.curve) else {
        panic!("a fitted arc is a B-spline");
    };
    let [from, to] = curve.domain().map(Rational::to_f64_nearest);
    let mut dec = Decider::default();
    (0..=count)
        .map(|k| {
            let u = from + (to - from) * k as f64 / count as f64;
            let at = Real::from(Rational::from_f64_exact(u).unwrap());
            let p = curve.point_at(&at, &mut dec).unwrap();
            let p = [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()];
            let gap = [1.0, -1.0]
                .map(|s| {
                    let g = true_generator(d, s, u);
                    ((p[0] - g[0]).powi(2) + (p[1] - g[1]).powi(2) + (p[2] - g[2]).powi(2)).sqrt()
                })
                .into_iter()
                .fold(f64::INFINITY, f64::min);
            assert!(gap <= fit.deviation_mm, "{gap} at {u}");
            p
        })
        .collect()
}

#[test]
fn a_torus_seen_along_its_axis_shows_both_equator_circles_exactly() {
    let (model, body) = ring();
    let answer =
        contour_generators(&model, body, &vector(0, 0, 2), &mut Decider::default()).unwrap();
    assert!(answer.is_exact());
    let c = answer.value;
    assert_eq!(c.arcs.len(), 2);
    let mut radii = Vec::new();
    for arc in &c.arcs {
        assert!(arc.fit.is_none());
        assert_eq!((arc.start, arc.end), (None, None));
        let (center, normal, radius_sq) = circle_of(&c, arc);
        assert!(is_zero(&center.minus(&Point3::origin()).norm_squared()));
        assert!(is_zero(
            &normal.cross(&vector(0, 0, 1).to_real()).norm_squared()
        ));
        radii.push(radius_sq);
    }
    for expected in [36, 16] {
        assert!(radii.iter().any(|r| is_zero(&(r - &Real::from(expected)))));
    }
}

#[test]
fn a_torus_seen_across_its_axis_shows_two_meridians_and_the_top_and_bottom_circles() {
    let (model, body) = ring();
    let answer =
        contour_generators(&model, body, &vector(0, 3, 0), &mut Decider::default()).unwrap();
    assert!(answer.is_exact());
    let c = answer.value;
    assert_eq!(c.arcs.len(), 4);
    let mut meridians = Vec::new();
    let mut parallels = Vec::new();
    for arc in &c.arcs {
        assert!(arc.fit.is_none());
        assert_eq!((arc.start, arc.end), (None, None));
        let (center, normal, radius_sq) = circle_of(&c, arc);
        if is_zero(&(&radius_sq - &Real::one())) {
            // In the plane y = 0 through the axis, centred on the core circle.
            assert!(is_zero(
                &normal.cross(&vector(0, 1, 0).to_real()).norm_squared()
            ));
            assert!(is_zero(&center.y) && is_zero(&center.z));
            meridians.push(center.x.as_rational().unwrap());
        } else {
            assert!(is_zero(&(&radius_sq - &Real::from(25))));
            assert!(is_zero(
                &normal.cross(&vector(0, 0, 1).to_real()).norm_squared()
            ));
            assert!(is_zero(&center.x) && is_zero(&center.y));
            assert!(is_zero(&(&center.z.square() - &Real::one())));
            parallels.push(sign(&center.z, &mut Decider::default()).unwrap());
        }
    }
    meridians.sort();
    assert_eq!(meridians, [Rational::from(-5), Rational::from(5)]);
    parallels.sort();
    assert_eq!(parallels, [-1, 1]);
}

#[test]
fn an_oblique_torus_gives_two_fitted_closed_branches_marked_approximate() {
    let (model, body) = ring();
    let d = vector(0, 1, 1);
    let answer = contour_generators(&model, body, &d, &mut Decider::default()).unwrap();
    assert!(!answer.is_exact());
    let marks: Vec<_> = answer.marks.iter().collect();
    assert_eq!(marks.len(), 1);
    assert_eq!(marks[0].code, CONTOUR_FIT_MARK);
    let c = answer.value;
    assert_eq!(c.arcs.len(), 2);
    let mut worst = 0.0_f64;
    for arc in &c.arcs {
        assert_eq!((arc.start, arc.end), (None, None));
        assert_eq!(arc.on_faces, vec![arc.face]);
        let fit = arc.fit.as_ref().unwrap();
        // Rule T2 for a 6 mm frame and a 1 mm tube: max(1e-6, 10 x 1e-7).
        assert_eq!(
            fit.target_mm.partial_cmp(&1e-6),
            Some(std::cmp::Ordering::Equal)
        );
        assert!(fit.spans > 1);
        worst = worst.max(fit.tolerance_mm);
        let samples = fitted_samples(&c, arc, [0.0, 1.0, 1.0], 48);
        // A whole branch closes on its own first point.
        assert_eq!(
            samples[0].partial_cmp(&samples[48]),
            Some(std::cmp::Ordering::Equal)
        );
    }
    assert_eq!(marks[0].bound, Rational::from_f64_exact(worst).unwrap());
}

#[test]
fn an_oblique_half_torus_is_cut_exactly_on_its_meridians() {
    // Remove x < 0: two meridian circles bound what is left.
    let (model, body) = ring_minus(vector(-10, -10, -10), [10, 20, 20]);
    let d = vector(1, 2, 2);
    let c = contour_generators(&model, body, &d, &mut Decider::default())
        .unwrap()
        .value;
    assert_eq!(c.arcs.len(), 2);
    for arc in &c.arcs {
        for end in [arc.start.unwrap(), arc.end.unwrap()] {
            let p = c.model.geom().position(end);
            assert!(is_zero(&p.x));
            // On a meridian circle, with (p - m).d = 0 for its tube centre m.
            let on_meridian = [5, -5].into_iter().any(|y| {
                let m = Point3::new(Real::zero(), Real::from(y), Real::zero());
                is_zero(&(&p.minus(&m).norm_squared() - &Real::one()))
                    && is_zero(&p.minus(&m).dot(&d.to_real()))
            });
            assert!(on_meridian);
        }
        for p in fitted_samples(&c, arc, [1.0, 2.0, 2.0], 24) {
            assert!(p[0] >= -arc.fit.as_ref().unwrap().deviation_mm);
        }
    }
}

#[test]
fn an_oblique_torus_cap_is_cut_exactly_on_its_parallels() {
    // Remove z < 1/2: the parallels of radius 5 +- sqrt(3)/2 bound the cap.
    let mut model = Model::new();
    let body = torus(&mut model, 5.into(), 1.into(), &Pose::identity()).unwrap();
    let cutter = block(
        &mut model,
        [20.into(), 20.into(), Rational::parse("21/2").unwrap()],
        &Pose::from_rational(Rotation::identity(), vector(-10, -10, -10)),
    )
    .unwrap();
    let body = boolean(
        &mut model,
        body,
        cutter,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let d = vector(0, 1, 1);
    let c = contour_generators(&model, body, &d, &mut Decider::default())
        .unwrap()
        .value;
    // One arc per branch: |sin u| > 1/sqrt(3) on the side the branch rises.
    assert_eq!(c.arcs.len(), 2);
    let half = Real::from(Rational::parse("1/2").unwrap());
    for arc in &c.arcs {
        for end in [arc.start.unwrap(), arc.end.unwrap()] {
            let p = c.model.geom().position(end);
            assert!(is_zero(&(&p.z - &half)));
            // N . d = 0: ((p - c).d) |q| = R (q.d), q the radial part.
            let q = Vector3::new(p.x.clone(), p.y.clone(), Real::zero());
            let dr = d.to_real();
            let lhs = &p.minus(&Point3::origin()).dot(&dr) * &q.norm_squared().sqrt();
            assert!(is_zero(&(&lhs - &(&Real::from(5) * &q.dot(&dr)))));
        }
        let deviation = arc.fit.as_ref().unwrap().deviation_mm;
        for p in fitted_samples(&c, arc, [0.0, 1.0, 1.0], 24) {
            assert!(p[2] >= 0.5 - deviation);
        }
    }
}

#[test]
fn a_generator_along_two_neighbouring_flat_faces_splits_where_its_neighbour_changes() {
    // A half cylinder y >= 0, radius 2, z in [0, 5], joined below z = 2 to a
    // slab y in [-1, 0] wider than it. Seen along y, each generator line
    // x = +-2 runs along an edge with the slab's strip y = 0, |x| >= 2 for
    // z < 2, and along an edge with the half cylinder's own flat face above.
    let mut model = Model::new();
    let round = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let cutter = block(
        &mut model,
        [10.into(), 10.into(), 10.into()],
        &Pose::from_rational(Rotation::identity(), vector(-5, -10, -1)),
    )
    .unwrap();
    let half = boolean(
        &mut model,
        round,
        cutter,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let slab = block(
        &mut model,
        [6.into(), 1.into(), 2.into()],
        &Pose::from_rational(Rotation::identity(), vector(-3, -1, 0)),
    )
    .unwrap();
    let body = boolean(
        &mut model,
        half,
        slab,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let front = contours(&model, body, vector(0, 1, 0));
    assert_eq!(front.arcs.len(), 4);
    let mut dec = Decider::default();
    for x in [2, -2] {
        let mut spans = Vec::new();
        for arc in &front.arcs {
            let [a, b] = ends(&front, arc);
            if a.x.as_rational().unwrap() != Rational::from(x) {
                continue;
            }
            let mut z = [a.z.as_rational().unwrap(), b.z.as_rational().unwrap()];
            z.sort();
            assert_eq!(arc.on_faces.len(), 2);
            let other = arc.on_faces[1];
            // The neighbour holds a point just beside the arc's middle on its
            // own side: the strip beyond |x| = 2 below z = 2, the half
            // cylinder's flat face inside it above.
            let middle = (&z[0] + &z[1]) * Rational::from_ratio(1, 2).unwrap();
            let beside = if z[1] <= Rational::from(2) {
                Rational::from(x) * Rational::from_ratio(5, 4).unwrap()
            } else {
                Rational::from(x) * Rational::from_ratio(3, 4).unwrap()
            };
            let probe = RVector3::new(beside, Rational::zero(), middle).to_real();
            let probe = Point3::origin().plus(&probe);
            assert!(face_contains(&mut model.clone(), other, &probe, &mut dec).unwrap());
            spans.push(z);
        }
        spans.sort();
        assert_eq!(
            spans,
            [
                [Rational::zero(), Rational::from(2)],
                [Rational::from(2), Rational::from(5)]
            ]
        );
    }
}

#[test]
fn tiny_view_directions_fit_the_torus_as_unit_sized_ones_do() {
    // Both directions square to zero in binary64. The exact power-of-two
    // rescale keeps 2^-700 (1, -1, 1) bit for bit the same as (1, -1, 1);
    // 1e-200 is no power of two, so its rounded direction may differ in the
    // last place, and its fit with it.
    let (model, body) = ring();
    let scaled = |factor: Rational| RVector3::new(factor.clone(), -factor.clone(), factor);
    let plain = contours(&model, body, vector(1, -1, 1));
    let binary = contours(&model, body, scaled(Rational::one().mul_pow2(-700)));
    let decimal = contours(&model, body, scaled(Rational::parse("1e-200").unwrap()));
    assert_eq!(plain.arcs.len(), 2);
    for c in [&binary, &decimal] {
        assert_eq!(c.arcs.len(), plain.arcs.len());
    }
    for ((a, b), c) in plain.arcs.iter().zip(&binary.arcs).zip(&decimal.arcs) {
        let fit = a.fit.as_ref().unwrap();
        assert_eq!(a.fit, b.fit);
        let other = c.fit.as_ref().unwrap();
        assert_eq!(other.spans, fit.spans);
        assert!(other.deviation_mm <= other.target_mm);
        assert!((other.deviation_mm - fit.deviation_mm).abs() <= 1e-6 * fit.deviation_mm);
    }
}
