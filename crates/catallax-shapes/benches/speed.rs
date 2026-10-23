//! Reproducible timings without speed thresholds.

use catallax_numbers::{Ball, BoxRoots, Decider, Outcome, PolySystem, Rational, Real, RootBox};
use catallax_shapes::*;
use std::{hint::black_box, sync::Arc, time::Instant};

fn point(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}
fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}
fn ball(x: i64, p: u32) -> Ball {
    Ball::from_rational(&x.into(), p)
}
fn timed(label: &str, n: usize, mut f: impl FnMut()) {
    let start = Instant::now();
    for _ in 0..n {
        f();
    }
    println!(
        "{label}: {:.0} ns/op ({n} iterations)",
        start.elapsed().as_nanos() as f64 / n as f64
    );
}
fn main() {
    let origin = point(0, 0, 0);
    let axis = vector(0, 0, 1);
    let reference = vector(1, 0, 0);
    let surfaces = [
        Surface::Plane(
            Plane::new(
                (origin.clone()).into(),
                (axis.clone()).into(),
                (reference.clone()).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Surface::Cylinder(
            Cylinder::new(
                (origin.clone()).into(),
                (axis.clone()).into(),
                5.into(),
                (reference.clone()).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Surface::Cone(
            Cone::new(
                (origin.clone()).into(),
                (axis.clone()).into(),
                (Rational::from_ratio(1, 2).unwrap()).into(),
                (reference.clone()).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Surface::Sphere(
            Sphere::new(
                (origin.clone()).into(),
                5.into(),
                (axis.clone()).into(),
                (reference.clone()).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Surface::Torus(
            Torus::new(
                (origin.clone()).into(),
                (axis.clone()).into(),
                3.into(),
                1.into(),
                (reference.clone()).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
    ];
    let mut store = GeomStore::new();
    let ids = surfaces.clone().map(|s| store.add_surface(s));
    let samples: Vec<_> = (0..64)
        .map(|i| point(i - 30, 2 * i - 7, 11 - i).to_real())
        .collect();
    let e = Rational::parse("1e-30").unwrap();
    let larger = Rational::one() + e.clone();
    let near = [
        RPoint3::new(5.into(), 7.into(), e),
        RPoint3::new(
            &Rational::from_i64(3) * &larger,
            &Rational::from_i64(4) * &larger,
            7.into(),
        ),
        RPoint3::new(
            3.into(),
            4.into(),
            Rational::from_i64(10) * (Rational::from_i64(2) - larger.clone()),
        ),
        RPoint3::new(
            &Rational::from_i64(3) * &larger,
            &Rational::from_i64(4) * &larger,
            0.into(),
        ),
        RPoint3::new(
            Rational::parse("54/25").unwrap(),
            Rational::parse("72/25").unwrap(),
            Rational::parse("4/5").unwrap() * larger,
        ),
    ]
    .map(|p| p.to_real());
    let sqrt_half = Real::from(Rational::from_ratio(1, 2).unwrap()).sqrt();
    let sqrt_points = [
        (
            ids[1],
            Point3::new(
                &Real::from_i64(5) * &sqrt_half,
                &Real::from_i64(5) * &sqrt_half,
                Real::zero(),
            ),
        ),
        (
            ids[3],
            Point3::new(
                Real::from(Rational::from_ratio(5, 2).unwrap()),
                Real::from(Rational::from_ratio(5, 2).unwrap()),
                Real::from_i64(5) * sqrt_half,
            ),
        ),
    ];
    let pose = Pose::new(
        Rotation::from_quaternion(1.into(), 2.into(), 3.into(), 4.into()).unwrap(),
        (vector(1_000_000, -7, 9)).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    for (index, (surface, id)) in surfaces.iter().zip(ids).enumerate() {
        let mut dec = Decider::default();
        let mut sample = 0;
        timed(&format!("side/{:?}/general", surface.kind()), 256, || {
            let _ = black_box(store.side_of(id, &samples[sample % 64], &mut dec));
            sample += 1;
        });
        println!(
            "  filter fraction: {}/{}",
            dec.stats().filter,
            dec.stats().decisions
        );
        timed(&format!("side/{:?}/near", surface.kind()), 32, || {
            let _ = black_box(store.side_of(id, &near[index], &mut dec));
        });
        println!(
            "  max bits: {}, work: {}",
            dec.stats().max_bits,
            dec.stats().work
        );
        timed(&format!("implicit/{:?}", surface.kind()), 32, || {
            black_box(surface.implicit());
        });
        timed(&format!("transformed/{:?}", surface.kind()), 32, || {
            black_box(surface.transformed(&pose));
        });
        for precision in [64, 1000] {
            let u = ball(1, precision);
            let v = ball(1, precision);
            let p = surface
                .evaluate(&u, &v, &catallax_numbers::Decider::default())
                .unwrap();
            timed(
                &format!("evaluate/{:?}/{precision}", surface.kind()),
                4,
                || {
                    black_box(
                        surface
                            .evaluate(&u, &v, &catallax_numbers::Decider::default())
                            .unwrap(),
                    );
                },
            );
            timed(
                &format!("derivatives/{:?}/{precision}", surface.kind()),
                4,
                || {
                    black_box(
                        surface
                            .derivatives(&u, &v, &catallax_numbers::Decider::default())
                            .unwrap(),
                    );
                },
            );
            timed(
                &format!("parameters/{:?}/{precision}", surface.kind()),
                4,
                || {
                    black_box(
                        surface
                            .parameters_of(&p, &catallax_numbers::Decider::default())
                            .unwrap(),
                    );
                },
            );
        }
    }
    for (id, point) in sqrt_points {
        let mut dec = Decider::default();
        timed(
            &format!("side/{:?}/sqrt", store.surface(id).kind()),
            4,
            || {
                let _ = black_box(store.side_of(id, &point, &mut dec));
            },
        );
        println!(
            "  max bits: {}, work: {}",
            dec.stats().max_bits,
            dec.stats().work
        );
    }
    let p = store.add_point(PointDef::Explicit((point(3, 4, 0)).into()));
    let curves = [
        Curve::Line(
            Line::new(
                (origin.clone()).into(),
                (vector(2, -1, 0)).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Curve::Circle(
            Circle::new(
                (origin.clone()).into(),
                (axis.clone()).into(),
                25.into(),
                (reference.clone()).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Curve::Ellipse(
            Ellipse::new(
                (origin.clone()).into(),
                (axis.clone()).into(),
                (reference.clone()).into(),
                25.into(),
                9.into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Curve::Intersection(IntersectionCurve {
            surfaces: [ids[3], ids[0]],
            witness: p,
        }),
    ];
    let other = store.add_point(PointDef::Explicit((point(3, 4, 1)).into()));
    for curve in curves {
        let kind = curve.kind();
        let id = store.add_curve(curve);
        let mut dec = Decider::default();
        timed(&format!("point_on_curve/{kind:?}"), 64, || {
            let _ = black_box(store.point_on_curve(id, p, &mut dec));
        });
        for precision in [64, 1000] {
            let t = ball(1, precision);
            timed(&format!("curve_evaluate/{kind:?}/{precision}"), 4, || {
                black_box(
                    store
                        .curve(id)
                        .evaluate(&t, &catallax_numbers::Decider::default())
                        .unwrap(),
                );
            });
        }
    }
    let mut dec = Decider::default();
    timed("same_point", 64, || {
        let _ = black_box(store.same_point(p, other, &mut dec));
    });
    for id in ids {
        timed(
            &format!("normals_agree/{:?}", store.surface(id).kind()),
            32,
            || {
                let _ = black_box(store.normals_agree(id, ids[0], p, &mut dec));
            },
        );
    }
    let a = store.add_surface(Surface::Sphere(
        Sphere::new(
            (point(6, 0, 0)).into(),
            5.into(),
            (axis.clone()).into(),
            (reference.clone()).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let c = store.add_surface(Surface::Sphere(
        Sphere::new(
            (point(0, 6, 0)).into(),
            5.into(),
            (axis.clone()).into(),
            (reference.clone()).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let system = Arc::new(
        PolySystem::new(vec![
            store.implicit(ids[3]).unwrap().rational().unwrap().clone(),
            store.implicit(a).unwrap().rational().unwrap().clone(),
            store.implicit(c).unwrap().rational().unwrap().clone(),
        ])
        .unwrap(),
    );
    let candidate = RootBox::new(
        vec![2.into(), 2.into(), 2.into()],
        vec![4.into(), 4.into(), 3.into()],
    )
    .unwrap();
    let root = match dec.certify_root(&system, &candidate) {
        Outcome::Proven(BoxRoots::One(root)) => root,
        other => panic!("root setup: {other:?}"),
    };
    let p = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(vec![ids[3], a, c]),
        position: Position::Isolated(root),
    }));
    timed("side/certified-root", 4, || {
        let _ = black_box(store.side_of(ids[3], &store.position(p), &mut dec));
    });
    timed("side/fact", 256, || {
        let _ = black_box(store.side(ids[3], p, &mut dec));
    });
    let mut mixed = Decider::default();
    timed("side/mixed", 128, || {
        let _ = black_box(store.side(ids[3], p, &mut mixed));
        let _ = black_box(store.side_of(ids[3], &near[3], &mut mixed));
    });
    println!(
        "mixed fact fraction: {}/{}",
        mixed.stats().coincident_fact,
        mixed.stats().decisions
    );
    for bits in [64, 1000] {
        let width = Rational::one().mul_pow2(-bits);
        timed(&format!("enclose_point/2^-{bits}"), 4, || {
            black_box(store.enclose_point(p, &width, &mut dec).unwrap());
        });
    }
    timed("rotation/from_quaternion", 128, || {
        black_box(Rotation::from_quaternion(1.into(), 2.into(), 3.into(), 4.into()).unwrap());
    });
}
