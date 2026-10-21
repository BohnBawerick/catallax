//! Non-gating timings for exact configurations, curve questions and traces.

use catallax_intersect::*;
use catallax_numbers::*;
use catallax_shapes::*;
use std::time::{Duration, Instant};

fn r(x: i64) -> Rational {
    x.into()
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
fn plane(z: i64) -> Surface {
    Surface::Plane(Plane::new(p(0, 0, z), v(0, 0, 1), v(1, 0, 0)).unwrap())
}
fn cylinder(x: i64, y: i64, w: RVector3, radius: i64) -> Surface {
    let reference = if w.cross(&v(1, 0, 0)).is_zero() {
        v(0, 0, 1)
    } else if w.y.is_zero() {
        v(0, 1, 0)
    } else {
        v(0, 0, 1)
    };
    Surface::Cylinder(Cylinder::new(p(x, y, 0), w, r(radius), reference).unwrap())
}
fn sphere(x: i64, z: i64, radius: i64) -> Surface {
    Surface::Sphere(Sphere::new(p(x, 0, z), r(radius), v(0, 0, 1), v(1, 0, 0)).unwrap())
}
fn ns(label: &str, elapsed: Duration, n: u32) {
    println!(
        "{label}: {:.0} ns/call",
        elapsed.as_nanos() as f64 / f64::from(n)
    );
}
fn surfaces() {
    let mut pairs = vec![];
    let mut add = |name, a, b| pairs.push((name, a, b));
    let z = cylinder(0, 0, v(0, 0, 1), 5);
    let xplane = |x| Surface::Plane(Plane::new(p(x, 0, 0), v(1, 0, 0), v(0, 1, 0)).unwrap());
    add("PP1", plane(0), xplane(1));
    add("PP2", plane(0), plane(1));
    add("PP3", plane(0), plane(0));
    add("PC1", plane(0), z.clone());
    add(
        "PC2",
        Surface::Plane(Plane::new(p(0, 0, 0), v(0, 3, 4), v(1, 0, 0)).unwrap()),
        z.clone(),
    );
    for (name, x) in [
        ("PC3 rational", 3),
        ("PC3 irrational", 1),
        ("PC4", 5),
        ("PC5", 6),
    ] {
        add(name, xplane(x), z.clone());
    }
    for (name, distance) in [("PS1", 3), ("PS2", 5), ("PS3", 6)] {
        add(name, plane(distance), sphere(0, 0, 5));
    }
    for (name, x, radius) in [
        ("SS1", 0, 5),
        ("SS2", 0, 4),
        ("SS3", 6, 5),
        ("SS4", 10, 5),
        ("SS5", 2, 3),
        ("SS6", 1, 3),
    ] {
        add(name, sphere(0, 0, 5), sphere(x, 0, radius));
    }
    for (name, x, radius) in [
        ("CP1", 0, 5),
        ("CP2", 0, 4),
        ("CP3 rational", 8, 5),
        ("CP3 irrational", 1, 5),
        ("CP4", 10, 5),
        ("CP5", 2, 3),
        ("CP6", 12, 5),
    ] {
        add(name, z.clone(), cylinder(x, 0, v(0, 0, 1), radius));
    }
    for (name, y) in [("CC1", 4), ("CC2", 3), ("CC3", 2), ("CC4", 1), ("CC5", 0)] {
        add(
            name,
            cylinder(0, 0, v(0, 0, 1), 2),
            cylinder(0, y, v(1, 0, 0), 1),
        );
    }
    add(
        "CC6",
        cylinder(0, 0, v(0, 0, 1), 1),
        cylinder(0, 0, v(3, 0, 4), 1),
    );
    add(
        "CC6 refusal",
        cylinder(0, 0, v(0, 0, 1), 1),
        cylinder(0, 0, v(1, 1, 0), 1),
    );
    for (name, radius) in [
        ("CS1 rational", 13),
        ("CS1 irrational", 6),
        ("CS2", 5),
        ("CS3", 4),
    ] {
        add(name, z.clone(), sphere(0, 1, radius));
    }
    for (name, x, radius) in [
        ("CS4", 4, 2),
        ("CS5", 3, 2),
        ("CS6", 2, 2),
        ("CS7", 1, 2),
        ("CS8", 1, 3),
    ] {
        add(name, cylinder(0, 0, v(0, 0, 1), 1), sphere(x, 0, radius));
    }
    let mut seed = 0x8137ab4du64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % 201) as i64 - 100
    };
    let mut poses = vec![Pose::identity()];
    for _ in 0..8 {
        poses.push(Pose {
            rotation: Rotation::from_quaternion(r(101), r(next()), r(next()), r(next())).unwrap(),
            translation: v(next() * 10000, next() * 10000, next() * 10000),
        });
    }
    for (name, a, b) in pairs {
        for (label, group) in [("identity", &poses[..1]), ("posed", &poses[1..])] {
            let mut elapsed = Duration::ZERO;
            let mut calls = 0;
            for pose in group {
                let mut source = GeomStore::new();
                let a = source.add_surface(a.transformed(pose));
                let b = source.add_surface(b.transformed(pose));
                let mut dec = Decider::new(BudgetPolicy::default());
                for _ in 0..10 {
                    let mut geom = source.clone();
                    let start = Instant::now();
                    let result = intersect_surfaces(&mut geom, a, b, &mut dec);
                    elapsed += start.elapsed();
                    let _ = std::hint::black_box(result);
                    calls += 1;
                }
            }
            ns(&format!("surfaces/{name}/{label}"), elapsed, calls);
        }
    }
}
fn curves() {
    let mut geom = GeomStore::new();
    let mut dec = Decider::new(BudgetPolicy::default());
    let plane = geom.add_surface(plane(1));
    let cylinder = geom.add_surface(cylinder(0, 0, v(0, 0, 1), 1));
    let sphere = geom.add_surface(sphere(1, 0, 2));
    let line = geom.add_curve(Curve::Line(Line::new(p(-5, 0, 0), v(1, 0, 0)).unwrap()));
    let circle = geom.add_curve(Curve::Circle(
        Circle::new(p(0, 0, 0), v(0, 0, 1), r(25), v(1, 0, 0)).unwrap(),
    ));
    let ellipse = geom.add_curve(Curve::Ellipse(
        Ellipse::new(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0), r(25), r(9)).unwrap(),
    ));
    let SurfacePair::Meets(components) = intersect_surfaces(&mut geom, cylinder, sphere, &mut dec)
        .unwrap()
        .value
    else {
        unreachable!()
    };
    let Component::Curve { curve: branch, .. } = components[0] else {
        unreachable!()
    };
    let branch_point = RPoint3::new(
        Rational::from_ratio(-7, 25).unwrap(),
        Rational::from_ratio(24, 25).unwrap(),
        Rational::from_ratio(6, 5).unwrap(),
    );
    for (name, c, positions) in [
        ("line", line, [p(-5, 0, 0), p(0, 0, 0), p(5, 0, 0)]),
        ("circle", circle, [p(5, 0, 0), p(0, 5, 0), p(-5, 0, 0)]),
        ("ellipse", ellipse, [p(5, 0, 0), p(0, 3, 0), p(-5, 0, 0)]),
        ("branch", branch, [p(1, 0, 2), branch_point, p(-1, 0, 0)]),
    ] {
        let points = positions.map(|p| geom.add_point(PointDef::Explicit(p)));
        for (surface_name, s) in [("plane", plane), ("cylinder", cylinder), ("sphere", sphere)] {
            let mut elapsed = Duration::ZERO;
            let _ = dec.take_stats();
            for _ in 0..10 {
                let mut g = geom.clone();
                let start = Instant::now();
                let _ = std::hint::black_box(intersect_curve_surface(&mut g, c, s, &mut dec));
                elapsed += start.elapsed();
            }
            ns(&format!("curve-surface/{name}/{surface_name}"), elapsed, 10);
            let stats = dec.take_stats();
            println!(
                "  decisions={} fact={} filter={} max_bits={} work={}",
                stats.decisions, stats.coincident_fact, stats.filter, stats.max_bits, stats.work
            );
        }
        let start = Instant::now();
        for _ in 0..10 {
            let _ = std::hint::black_box(on_curve(&geom, c, points[0], &mut dec));
        }
        ns(&format!("on-curve/{name}"), start.elapsed(), 10);
        let start = Instant::now();
        for _ in 0..10 {
            let _ = std::hint::black_box(same_curve(&geom, c, c, &mut dec));
        }
        ns(&format!("same-curve/{name}"), start.elapsed(), 10);
        if name == "line" {
            let start = Instant::now();
            for _ in 0..10 {
                let _ =
                    std::hint::black_box(order_on_curve(&geom, c, points[0], points[1], &mut dec));
            }
            ns("order/line", start.elapsed(), 10);
        } else {
            let start = Instant::now();
            for _ in 0..10 {
                let _ = std::hint::black_box(cyclic_order(&geom, c, points, &mut dec));
            }
            ns(&format!("cyclic/{name}"), start.elapsed(), 10);
        }
        let (from, to) = if name == "line" {
            (Some(points[0]), Some(points[2]))
        } else {
            (None, None)
        };
        let start = Instant::now();
        let trace = trace(&geom, c, from, to, &mut dec).unwrap().value;
        ns(&format!("trace/{name}"), start.elapsed(), 1);
        for bits in [64, 1000] {
            for order in [4, 16] {
                let t = Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), bits);
                let start = Instant::now();
                for piece in trace.pieces() {
                    std::hint::black_box(piece.taylor(&t, order, &Decider::default()).unwrap());
                }
                ns(
                    &format!("taylor/{name}/order-{order}/{bits}-bits"),
                    start.elapsed(),
                    trace.pieces().len() as u32,
                );
            }
        }
    }
    let mut max_bits = 0;
    for exponent in [-30, -60] {
        let gap = Rational::parse(&format!("1e{exponent}")).unwrap();
        for sign in [-1, 1] {
            let y = r(5) + gap.clone() * r(sign);
            let line = geom.add_curve(Curve::Line(
                Line::new(RPoint3::new(r(-10), y, r(0)), v(1, 0, 0)).unwrap(),
            ));
            let s = geom.add_surface(Surface::Sphere(
                Sphere::new(p(0, 0, 0), r(5), v(0, 0, 1), v(1, 0, 0)).unwrap(),
            ));
            let mut g = geom.clone();
            let _ = std::hint::black_box(intersect_curve_surface(&mut g, line, s, &mut dec));
            max_bits = max_bits.max(dec.take_stats().max_bits);
        }
    }
    println!("near-miss/max-precision-bits: {max_bits}");
}
fn main() {
    surfaces();
    curves();
}
