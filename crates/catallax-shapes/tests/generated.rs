#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use catallax_shapes::*;
use common::cases::*;
use common::*;

fn posed(point: &RPoint3, pose: &Pose) -> RPoint3 {
    let point = qp(point);
    let t = qv(&pose.translation);
    rp(std::array::from_fn(|i| {
        eq(&((0..3)
            .map(|j| oq(&pose.rotation.matrix()[i][j]) * &point[j])
            .sum::<Q>()
            + &t[i]))
    }))
}
fn rational_circle(n: i64) -> (Q, Q, i64) {
    // k=n/8 gives c=(64-n²)/(64+n²), s=16n/(64+n²).
    let den = 64 + n * n;
    (q(64 - n * n, den), q(16 * n, den), den)
}
fn point_q(a: [Q; 3]) -> RPoint3 {
    rp(a.map(|x| eq(&x)))
}

blind_test!(cases_18_50_51_generated_32940_math_answered_decisions, {
    let mut rng = Rng(0x25051);
    let mut count = 0;
    let poses: Vec<_> = (0..20).map(|_| rng.pose()).collect();
    for pose in poses {
        for exponent in -30..=30 {
            let scale = r(&format!("1e{exponent}"));
            let e = q(1, 1) / Q::from(rug::Integer::from(10).pow(30));
            let (c, s, _) = rational_circle((rng.next() % 7 + 1) as i64);
            let (c2, s2, _) = rational_circle((rng.next() % 7 + 1) as i64);
            for offset in [-1, 0, 1] {
                let t: Q = q(1, 1) + e.clone() * offset;
                let local = [
                    p("5", "7", &t.to_string()),
                    point_q([q(5, 1) * &c * &t, q(5, 1) * &s * &t, q(7, 1)]),
                    point_q([q(5, 1) * &c * &t, q(5, 1) * &s * &t, q(10, 1)]),
                    point_q([
                        q(3, 1) * &c * &c2 * &t,
                        q(3, 1) * &s * &c2 * &t,
                        q(3, 1) * &s2 * &t,
                    ]),
                    point_q([(q(3, 1) + &c2) * &c, (q(3, 1) + &c2) * &s, s2.clone() * &t]),
                ];
                for (surface, point) in surfaces().into_iter().zip(local) {
                    let expected = sign(&implicit_oracle(&surface, &point));
                    assert_eq!(expected, offset as i8);
                    let point = posed(&scale_point(&point, &scale), &pose);
                    let surface = scale_surface(&surface, &scale).transformed(&pose);
                    let mut store = GeomStore::new();
                    let sid = store.add_surface(surface);
                    let pid = explicit(&mut store, point);
                    side_answer(store.side(sid, pid, &mut Decider::default()), expected);
                    count += 1;
                }
            }
            for (curve, point, on) in curve_cases() {
                let curve = scale_curve(&curve, &scale).transformed(&pose).unwrap();
                let point = posed(&scale_point(&point, &scale), &pose);
                let mut store = GeomStore::new();
                let c = store.add_curve(curve);
                let p = explicit(&mut store, point);
                incidence(store.point_on_curve(c, p, &mut Decider::default()), on);
                count += 1;
            }
        }
    }
    // 20 poses * 61 scales * (5 surfaces * 3 offsets + 12 curve questions).
    assert_eq!(count, 32_940);
});

blind_test!(cases_50_51_named_rational_cases_scaled_and_placed, {
    let mut rng = Rng(0x251);
    let transforms: Vec<_> = (-30..=30)
        .map(|k| (r(&format!("1e{k}")), Pose::identity()))
        .chain((0..20).map(|_| (r("1"), rng.pose())))
        .collect();
    for (scale, pose) in transforms {
        for (surface, point, expected) in side_cases() {
            let mut store = GeomStore::new();
            let s = store.add_surface(scale_surface(&surface, &scale).transformed(&pose));
            let p = explicit(&mut store, posed(&scale_point(&point, &scale), &pose));
            side_answer(store.side(s, p, &mut Decider::default()), expected);
        }
        for (a, b, same) in surface_pairs() {
            let mut store = GeomStore::new();
            let a = store.add_surface(scale_surface(&a, &scale).transformed(&pose));
            let b = store.add_surface(scale_surface(&b, &scale).transformed(&pose));
            incidence(store.same_surface(a, b, &mut Decider::default()), same);
        }
        for (a, b, point, expected) in normal_cases() {
            let mut store = GeomStore::new();
            let a = store.add_surface(scale_surface(&a, &scale).transformed(&pose));
            let b = store.add_surface(scale_surface(&b, &scale).transformed(&pose));
            let p = explicit(&mut store, posed(&scale_point(&point, &scale), &pose));
            let out = store.normals_agree(a, b, p, &mut Decider::default());
            match expected {
                0 => coincident(out, None),
                1 => assert_eq!(out, Outcome::Proven(StrictSign::Positive)),
                -1 => assert_eq!(out, Outcome::Proven(StrictSign::Negative)),
                _ => unreachable!(),
            }
        }
        let mut store = GeomStore::new();
        for (other, same) in [
            (p("1", "2", "3"), true),
            (rp([r("1") + r("1e-30"), r("2"), r("3")]), false),
        ] {
            let a = explicit(
                &mut store,
                posed(&scale_point(&p("1", "2", "3"), &scale), &pose),
            );
            let b = explicit(&mut store, posed(&scale_point(&other, &scale), &pose));
            incidence(store.same_point(a, b, &mut Decider::default()), same);
        }
        let a = store.add_surface(scale_surface(&sphere("5"), &scale).transformed(&pose));
        let b = store.add_surface(
            scale_surface(
                &Surface::Plane(
                    Plane::new(p("0", "0", "3"), v("0", "0", "1"), v("1", "0", "0")).unwrap(),
                ),
                &scale,
            )
            .transformed(&pose),
        );
        let witness = explicit(
            &mut store,
            posed(&scale_point(&p("4", "0", "3"), &scale), &pose),
        );
        let c = store.add_curve(Curve::Intersection(IntersectionCurve {
            surfaces: [a, b],
            witness,
        }));
        for (point, on) in [
            (p("0", "4", "3"), true),
            (rp([r("0"), r("4"), r("3") + r("1e-30")]), false),
        ] {
            let p = explicit(&mut store, posed(&scale_point(&point, &scale), &pose));
            incidence(store.point_on_curve(c, p, &mut Decider::default()), on);
        }
    }
});

blind_test!(case_43_generated_liveness_with_exact_dyadic_boxes, {
    let mut rng = Rng(0x2431);
    for _ in 0..128 {
        let (c, s, d1) = rational_circle((rng.next() % 7 + 1) as i64);
        let (c2, s2, d2) = rational_circle((rng.next() % 7 + 1) as i64);
        let q: [Q; 4] = std::array::from_fn(|_| Q::from((rng.next() % 100 + 1) as i64));
        let norm: Q = q.iter().map(sq).sum();
        let [w, x, y, z] = q.clone().map(|x| eq(&x));
        let pose = Pose::new(
            Rotation::from_quaternion(w, x, y, z).unwrap(),
            v("1000000", "-777777", "314159"),
        );
        // Multiplying by d1*d2*|q|² cancels BOTH rational-circle denominators
        // and the rotation denominator. Every posed coordinate is an integer,
        // hence radius-zero Ball inputs really satisfy the liveness premise.
        let scale = Q::from(d1 * d2) * norm;
        let k = eq(&scale);
        let local = [
            (
                cylinder("5"),
                point_q([qv5(&scale, &c), qv5(&scale, &s), Q::from(7) * &scale]),
                Ball::from_rational(&(&r("7") * &k), 256),
            ),
            (
                cone(),
                point_q([scale.clone() * &c, scale.clone() * &s, Q::from(2) * &scale]),
                Ball::from_rational(&(&r("2") * &k), 256),
            ),
            (
                sphere("3"),
                point_q([
                    Q::from(3) * &scale * &c2 * &c,
                    Q::from(3) * &scale * &c2 * &s,
                    Q::from(3) * &scale * &s2,
                ]),
                Ball::atan2(
                    &Ball::from_rational(&eq(&s2), 256),
                    &Ball::from_rational(&eq(&c2), 256),
                )
                .unwrap(),
            ),
            (
                torus(),
                point_q([
                    scale.clone() * (Q::from(3) + &c2) * &c,
                    scale.clone() * (Q::from(3) + &c2) * &s,
                    scale.clone() * &s2,
                ]),
                Ball::atan2(
                    &Ball::from_rational(&eq(&s2), 256),
                    &Ball::from_rational(&eq(&c2), 256),
                )
                .unwrap(),
            ),
        ];
        let expected_u = Ball::atan2(
            &Ball::from_rational(&eq(&s), 256),
            &Ball::from_rational(&eq(&c), 256),
        )
        .unwrap();
        for (surface, point, expected_v) in local {
            let surface = scale_surface(&surface, &k).transformed(&pose);
            let point = posed(&point, &pose);
            let box_ = balls(&point, 160);
            for b in &box_ {
                assert_eq!(b.rad(), r("0"));
            }
            // c,s and c2,s2 are positive; each sine >=16/65 and cosine >=15/113.
            // Thus distance to every seam/axis exceeds L/100, including sphere poles.
            let uv = surface
                .parameters_of(&box_)
                .expect("liveness: exact box separated from seams and axes");
            overlap(&uv[0], &expected_u);
            overlap(&uv[1], &expected_v);
        }
    }
});
fn qv5(scale: &Q, c: &Q) -> Q {
    Q::from(5) * scale * c
}

fn decision_sequence(
    store: &GeomStore,
    s: SurfaceId,
    c: CurveId,
    a: PointId,
    b: PointId,
) -> (Vec<String>, DecisionStats) {
    let mut d = Decider::default();
    let outcomes = vec![
        format!("{:?}", store.side(s, a, &mut d)),
        format!("{:?}", store.side_of(s, &store.position(b), &mut d)),
        format!("{:?}", store.point_on_curve(c, a, &mut d)),
        format!("{:?}", store.same_point(a, b, &mut d)),
        format!("{:?}", store.same_surface(s, s, &mut d)),
        format!("{:?}", store.normals_agree(s, s, a, &mut d)),
    ];
    (outcomes, d.take_stats())
}
blind_test!(case_52_repeated_and_eight_thread_determinism, {
    let mut store = GeomStore::new();
    let s = store.add_surface(sphere("5"));
    let c = store.add_curve(circle("25"));
    let a = explicit(&mut store, p("3", "4", "0"));
    let b = explicit(&mut store, rp([r("3") + r("1e-30"), r("4"), r("0")]));
    let expected = decision_sequence(&store, s, c, a, b);
    assert_eq!(decision_sequence(&store, s, c, a, b), expected);
    let jobs: Vec<_> = (0..8)
        .map(|_| {
            let store = store.clone();
            std::thread::spawn(move || decision_sequence(&store, s, c, a, b))
        })
        .collect();
    for job in jobs {
        assert_eq!(job.join().unwrap(), expected);
    }
});

use rug::ops::Pow;
