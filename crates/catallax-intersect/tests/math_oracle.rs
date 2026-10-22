//! Checks the mathematical fixtures without executing any piece 4 function.
use catallax_numbers::RealPoly;
mod common;
use common::oracle::{Expected, Part, Shape};
use common::*;

// This is a fixture smoke check, not the mathematical oracle. Rational samples
// are substituted exactly; radical samples get a residual enclosure. Exact
// radical derivations live in verify_math.py. The contract tests also compare
// returned points with same_point and default-budget exact decisions.
fn enclose(x: &Real) -> Ball {
    Decider::default().enclose(x, &r("1e-80")).unwrap()
}
fn residual(equation: &RealPoly, p: &Point3) {
    let equation = equation.rational().expect("rational intersection fixture");
    if let Some(p) = rational_point(p) {
        assert_eq!(equation.eval(&[p.x, p.y, p.z]), n(0));
    } else {
        assert!(
            equation
                .eval_ball(&xyz(p).each_ref().map(enclose))
                .contains_zero()
        );
    }
}
fn singular(v: &Vector3) {
    for x in [&v.x, &v.y, &v.z] {
        if let Some(x) = x.as_rational() {
            assert_eq!(x, n(0));
        } else {
            assert!(enclose(x).contains_zero());
        }
    }
}

#[test]
fn surface_fixture_residual_smoke() {
    for case in cases::pairs() {
        let mut g = GeomStore::new();
        let a = g.add_surface(case.a.clone());
        let b = g.add_surface(case.b.clone());
        let Expected::Meets(parts) = oracle::surfaces(&case.a, &case.b) else {
            continue;
        };
        let check = |p: &Point3| {
            for s in [a, b] {
                residual(&g.implicit(s), p);
            }
        };
        for part in parts {
            match part {
                Part::Point(p) => {
                    check(&p);
                    let cross = case.a.gradient(&p).cross(&case.b.gradient(&p));
                    singular(&cross);
                }
                Part::Curve {
                    shape,
                    singular: nodes,
                    ..
                } => {
                    for p in nodes {
                        check(&p);
                        singular(&case.a.gradient(&p).cross(&case.b.gradient(&p)));
                    }
                    match shape {
                        Shape::Analytic(curve) => {
                            let p = curve.parameter_origin(&g);
                            check(&p);
                            match *curve {
                                Curve::Line(c) => check(&p.plus(&c.dir)),
                                Curve::Circle(c) => {
                                    check(&c.center.plus(&p.minus(&c.center).scaled(&rr(-1))))
                                }
                                Curve::Ellipse(c) => {
                                    let minor = c.normal.cross(&c.major_dir);
                                    check(
                                        &c.center.plus(
                                            &minor.scaled(
                                                &(&c.minor_sq / &minor.norm_squared()).sqrt(),
                                            ),
                                        ),
                                    );
                                }
                                _ => unreachable!(),
                            }
                        }
                        Shape::Intersection { witness, .. } => {
                            check(&witness);
                            assert!(
                                enclose(
                                    &case
                                        .a
                                        .gradient(&witness)
                                        .cross(&case.b.gradient(&witness))
                                        .norm_squared()
                                )
                                .lo()
                                    > n(0)
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn curve_hit_fixture_residual_smoke() {
    for case in curves::hits() {
        let mut g = GeomStore::new();
        let c = case.curve.add(&mut g);
        let s = g.add_surface(case.surface);
        if let curves::Hits::Meets(points) = case.expected {
            for (point, _) in points {
                residual(&g.implicit(s), &point);
                for equation in g.curve(c).equations(&g) {
                    residual(&equation, &point);
                }
            }
        }
    }
}
