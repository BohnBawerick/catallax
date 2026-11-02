//! The exact rolling-ball pipe between a plane and an oblique cylinder.
//!
//! The carrier's distance identity has bidegree (8, 4) on each span once the
//! denominators are cleared, so its exact vanishing on a 9 by 5 grid proves it
//! over the whole span, not only at the samples.

use catallax_numbers::{Ball, BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_shapes::{
    BSplineSurface, Cylinder, Ellipse, Plane, Point3, Pose, RollingPipe, Rotation, Vector3,
};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn real(value: &str) -> Real {
    Real::from(q(value))
}
fn vector(v: [&str; 3]) -> Vector3 {
    let [x, y, z] = v.map(real);
    Vector3::new(x, y, z)
}
fn point(v: [&str; 3]) -> Point3 {
    let [x, y, z] = v.map(real);
    Point3::new(x, y, z)
}
fn plane(origin: [&str; 3], normal: [&str; 3]) -> Plane {
    let normal = vector(normal);
    let [a, b, _] = [&normal.x, &normal.y, &normal.z].map(Clone::clone);
    let ref_dir = Vector3::new(b, -a, Real::zero());
    Plane {
        origin: point(origin),
        normal,
        ref_dir,
    }
}
fn cylinder(origin: [&str; 3], axis: [&str; 3], radius: &str) -> Cylinder {
    let axis = vector(axis);
    let [a, b, _] = [&axis.x, &axis.y, &axis.z].map(Clone::clone);
    let ref_dir = Vector3::new(b, -a, Real::zero());
    Cylinder {
        origin: point(origin),
        axis,
        radius_sq: real(radius).square(),
        ref_dir,
    }
}

struct Case {
    name: &'static str,
    plane: Plane,
    cylinder: Cylinder,
    same: [bool; 2],
    radius: &'static str,
    convex: bool,
}
fn tilted() -> Plane {
    plane(["0", "0", "6"], ["3/5", "0", "4/5"])
}
fn upright() -> Cylinder {
    cylinder(["0", "0", "0"], ["0", "0", "1"], "5")
}
fn rational_cases() -> Vec<Case> {
    vec![
        Case {
            name: "rim of a cut cylinder",
            plane: tilted(),
            cylinder: upright(),
            same: [true, true],
            radius: "1",
            convex: true,
        },
        Case {
            name: "boss on an oblique face",
            plane: tilted(),
            cylinder: upright(),
            same: [true, true],
            radius: "1",
            convex: false,
        },
        Case {
            name: "bore through an oblique face",
            plane: tilted(),
            cylinder: upright(),
            same: [true, false],
            radius: "1",
            convex: true,
        },
        Case {
            name: "reversed plane",
            plane: tilted(),
            cylinder: upright(),
            same: [false, true],
            radius: "3/2",
            convex: true,
        },
    ]
}
fn build(case: &Case) -> RollingPipe {
    RollingPipe::new(
        &case.plane,
        &case.cylinder,
        case.same,
        real(case.radius),
        case.convex,
        &mut decider(),
    )
    .unwrap_or_else(|e| panic!("{}: {e:?}", case.name))
}

/// Exact evaluation of the represented carrier, from its controls alone.
fn evaluate(carrier: &BSplineSurface, u: &Rational, v: &Rational) -> Point3 {
    let degrees = carrier.degrees();
    let [count, rows] = carrier.control_count();
    let span = (u.floor().to_i64().unwrap() as usize).min(3);
    let first = span * degrees[0];
    let tau = u.clone() - Rational::from_i64(span as i64);
    let basis = |degree: usize, t: &Rational| -> Vec<Real> {
        (0..=degree)
            .map(|i| {
                let choose = (0..i).fold(1i64, |c, k| c * (degree - k) as i64 / (k as i64 + 1));
                let mut value = Rational::from_i64(choose);
                for _ in 0..i {
                    value = value.clone() * t.clone();
                }
                let rest = Rational::one() - t.clone();
                for _ in i..degree {
                    value = value.clone() * rest.clone();
                }
                Real::from(value)
            })
            .collect()
    };
    let (bu, bv) = (basis(degrees[0], &tau), basis(degrees[1], v));
    let mut sum = [Real::zero(), Real::zero(), Real::zero()];
    let mut total = Real::zero();
    for (i, a) in bu.iter().enumerate() {
        for (j, b) in bv.iter().enumerate() {
            let at = (first + i) * rows + j;
            assert!(first + i < count);
            let weight = &carrier.weights()[at];
            let factor = &(a * b) * weight;
            let c = &carrier.controls()[at];
            sum = [
                &sum[0] + &(&c.x * &factor),
                &sum[1] + &(&c.y * &factor),
                &sum[2] + &(&c.z * &factor),
            ];
            total = &total + &factor;
        }
    }
    let [x, y, z] = sum.map(|value| &value / &total);
    Point3::new(x, y, z)
}
fn zero(value: Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(&value), Outcome::ProvenCoincident(_))
}
fn same_point(a: &Point3, b: &Point3, dec: &mut Decider) -> bool {
    zero(a.minus(b).norm_squared(), dec)
}
fn grid() -> Vec<(Rational, Rational)> {
    let mut samples = Vec::new();
    for span in 0..4 {
        for i in 0..=8 {
            for j in 0..=4 {
                samples.push((
                    Rational::from_i64(span) + Rational::from_ratio(i, 8).unwrap(),
                    Rational::from_ratio(j, 4).unwrap(),
                ));
            }
        }
    }
    samples
}

fn in_ellipse_plane(curve: &Ellipse, p: &Point3, dec: &mut Decider) -> bool {
    zero(p.minus(&curve.center).dot(&curve.normal), dec)
}

#[test]
fn the_rational_carrier_is_the_rolling_ball_definition_on_every_span() {
    for case in rational_cases() {
        let pipe = build(&case);
        assert!(pipe.carrier().has_rational_definition(), "{}", case.name);
        let mut dec = decider();
        for (u, v) in grid() {
            let on_carrier = evaluate(pipe.carrier(), &u, &v);
            assert!(
                same_point(&on_carrier, &pipe.point(&u, &v), &mut dec),
                "{} carrier leaves the law at {u:?} {v:?}",
                case.name
            );
            let gap = on_carrier.minus(&pipe.centre(&u)).norm_squared();
            assert!(
                zero(&gap - &pipe.radius().square(), &mut dec),
                "{}",
                case.name
            );
        }
    }
}

#[test]
fn the_identity_holds_for_irrational_placements() {
    // An axis and a normal with irrational lengths and a tilted cylinder.
    let case = Case {
        name: "skew placement",
        plane: plane(["1", "2", "3"], ["1", "2", "3"]),
        cylinder: cylinder(["1", "1", "1"], ["0", "1", "1"], "4"),
        same: [true, true],
        radius: "3/2",
        convex: true,
    };
    let pipe = build(&case);
    assert!(!pipe.carrier().has_rational_definition());
    let mut dec = decider();
    for (u, v) in [("0", "0"), ("1/3", "1/2"), ("7/4", "1"), ("3", "1/4")] {
        let [u, v] = [q(u), q(v)];
        let on_carrier = evaluate(pipe.carrier(), &u, &v);
        assert!(same_point(&on_carrier, &pipe.point(&u, &v), &mut dec));
    }
}

#[test]
fn both_rails_lie_on_their_supports_and_their_ellipses() {
    for case in rational_cases() {
        let pipe = build(&case);
        let mut dec = decider();
        let (rail, wall, spine) = (
            pipe.plane_rail(&mut dec).unwrap(),
            pipe.wall_rail(&mut dec).unwrap(),
            pipe.spine(&mut dec).unwrap(),
        );
        let radial = |p: &Point3| {
            let d = p.minus(&case.cylinder.origin);
            d.norm_squared()
                - d.dot(&case.cylinder.axis).square() / case.cylinder.axis.norm_squared()
        };
        let height = |p: &Point3| p.minus(&case.plane.origin).dot(&case.plane.normal);
        for (u, _) in grid() {
            let on_plane = evaluate(pipe.carrier(), &u, &Rational::from_i64(0));
            let on_wall = evaluate(pipe.carrier(), &u, &Rational::one());
            assert!(zero(height(&on_plane), &mut dec), "{}", case.name);
            assert!(in_ellipse_plane(&rail, &on_plane, &mut dec));
            assert!(zero(
                radial(&on_wall) - case.cylinder.radius_sq.clone(),
                &mut dec
            ));
            assert!(in_ellipse_plane(&wall, &on_wall, &mut dec));
            let centre = pipe.centre(&u);
            assert!(in_ellipse_plane(&spine, &centre, &mut dec));
            // A plane contact is one radius from its centre along the normal.
            let reach = on_plane.minus(&centre).cross(&case.plane.normal);
            assert!(zero(reach.norm_squared(), &mut dec), "{}", case.name);
        }
    }
}

#[test]
fn the_seam_closes_exactly() {
    for case in rational_cases() {
        let pipe = build(&case);
        let carrier = pipe.carrier();
        let [count, rows] = carrier.control_count();
        let mut dec = decider();
        for j in 0..rows {
            assert!(same_point(
                &carrier.controls()[j],
                &carrier.controls()[(count - 1) * rows + j],
                &mut dec
            ));
        }
    }
}

// The ellipse's own tangent, from four of its points and the chart's angle.
fn spine_tangent(pipe: &RollingPipe, u: &Rational) -> Vector3 {
    let at = |value: i64| pipe.centre(&Rational::from_i64(value));
    let half = Real::from(q("1/2"));
    let major = at(0).minus(&at(2)).scaled(&half);
    let minor = at(1).minus(&at(3)).scaled(&half);
    let span = u.floor().to_i64().unwrap().rem_euclid(4);
    let tau = u.clone() - u.floor().into();
    let wide = Rational::one() + tau.clone() * tau.clone();
    let cos = Real::from((Rational::one() - tau.clone() * tau.clone()) / wide.clone());
    let sin = Real::from(tau * Rational::from_i64(2) / wide);
    let [cos, sin] = match span {
        0 => [cos, sin],
        1 => [-&sin, cos],
        2 => [-&cos, -&sin],
        _ => [sin, -&cos],
    };
    minor.scaled(&cos).minus(&major.scaled(&sin))
}

// With the section in the centre's normal plane and of constant radius, the
// tangent along the path is `c' (1 - r k p.N)`, and the section tangent is normal
// to it with length `r`. `RollingPipe::new` bounds `r k p.N` below one, so this
// identity is what makes the chart regular everywhere.
#[test]
fn every_section_lies_in_the_normal_plane_of_the_centre_path() {
    for case in rational_cases() {
        let pipe = build(&case);
        let mut dec = decider();
        for (u, v) in grid() {
            let ray = pipe.point(&u, &v).minus(&pipe.centre(&u));
            assert!(
                zero(ray.dot(&spine_tangent(&pipe, &u)), &mut dec),
                "{}",
                case.name
            );
        }
    }
}

#[test]
fn a_ball_at_the_tightest_curve_of_its_path_refuses_and_a_smaller_one_builds() {
    // alpha = 3/5 and a wall radius of 2 reach the ellipse's curvature radius at
    // r = 3, because r = (2 + r) 3/5 there.
    let shallow = plane(["0", "0", "1"], ["4/5", "0", "3/5"]);
    let wall = cylinder(["0", "0", "0"], ["0", "0", "1"], "2");
    let at = |radius: &str| {
        RollingPipe::new(
            &shallow,
            &wall,
            [true, true],
            real(radius),
            false,
            &mut decider(),
        )
    };
    assert!(at("3").is_err());
    let pipe = at("29/10").unwrap();
    let dec = decider();
    for k in 0..100 {
        for j in 0..=8 {
            let u = Rational::from_ratio(2 * k + 1, 50).unwrap();
            let v = Rational::from_ratio(j, 8).unwrap();
            let [u, v] = [&u, &v].map(|value| Ball::from_rational(value, 256));
            let [du, dv] = pipe.carrier().derivatives(&u, &v, &dec).unwrap();
            let normal = [
                &(&du[1] * &dv[2]) - &(&du[2] * &dv[1]),
                &(&du[2] * &dv[0]) - &(&du[0] * &dv[2]),
                &(&du[0] * &dv[1]) - &(&du[1] * &dv[0]),
            ];
            let size = &(&normal[0].square() + &normal[1].square()) + &normal[2].square();
            assert!(size.lo() > Rational::from_i64(0), "{k} {j}");
        }
    }
}

fn distance_to_ellipse(p: [f64; 3], pipe: &RollingPipe) -> f64 {
    let mut best = f64::MAX;
    for k in 0..40_000 {
        let u = Rational::from_ratio(k, 10_000).unwrap();
        let c = pipe.centre(&u);
        let d = [
            c.x.to_f64() - p[0],
            c.y.to_f64() - p[1],
            c.z.to_f64() - p[2],
        ];
        best = best.min((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());
    }
    best
}

#[test]
fn a_brute_force_nearest_centre_agrees_in_floating_point() {
    let started = std::time::Instant::now();
    let pipe = build(&rational_cases().remove(0));
    eprintln!("exact carrier in {:?}", started.elapsed());
    for (u, v) in [("1/3", "1/2"), ("5/4", "1/8"), ("7/2", "7/8")] {
        let p = pipe.point(&q(u), &q(v));
        let found = distance_to_ellipse([p.x.to_f64(), p.y.to_f64(), p.z.to_f64()], &pipe);
        assert!((found - 1.0).abs() < 1e-6, "{u} {v}: {found}");
    }
}

#[test]
fn a_rigid_motion_moves_law_and_carrier_together() {
    let pipe = build(&rational_cases().remove(0));
    let mut dec = decider();
    let rotation = Rotation::from_quaternion(q("3/5"), q("4/5"), q("0"), q("0")).unwrap();
    let pose = Pose::new(rotation, vector(["1", "2", "3"]), &mut dec).unwrap();
    let moved = pipe.transformed(&pose);
    let (u, v) = (q("5/4"), q("1/2"));
    assert!(same_point(
        &evaluate(moved.carrier(), &u, &v),
        &moved.point(&u, &v),
        &mut dec
    ));
    assert!(same_point(
        &pose.apply_point(&pipe.point(&u, &v)),
        &moved.point(&u, &v),
        &mut dec
    ));
}

#[test]
fn only_strictly_oblique_regular_contacts_are_accepted() {
    let refused = |plane: &Plane, cylinder: &Cylinder, same, radius: &str, convex| {
        RollingPipe::new(plane, cylinder, same, real(radius), convex, &mut decider()).unwrap_err()
    };
    let cylinder = upright();
    let normal_to_axis = plane(["0", "0", "6"], ["0", "0", "1"]);
    let along_axis = plane(["0", "0", "6"], ["1", "0", "0"]);
    for (plane, radius, convex, code) in [
        (&normal_to_axis, "1", true, "shapes.rolling-pipe"),
        (&along_axis, "1", true, "shapes.rolling-pipe"),
        // The offset reaches the axis.
        (&tilted(), "5", true, "shapes.rolling-pipe"),
        (&tilted(), "0", true, "shapes.non-positive-radius"),
    ] {
        assert_eq!(
            refused(plane, &cylinder, [true, true], radius, convex).code,
            code
        );
    }
    // Nearly along the axis, a ball outside the wall folds at the ellipse tip.
    let shallow = plane(["0", "0", "6"], ["99/101", "0", "20/101"]);
    let thin = cylinder_of_radius("1");
    assert_eq!(
        refused(&shallow, &thin, [true, true], "1", false).code,
        "shapes.rolling-pipe"
    );
}
fn cylinder_of_radius(radius: &str) -> Cylinder {
    cylinder(["0", "0", "0"], ["0", "0", "1"], radius)
}

fn triangle_hit(o: [f64; 3], d: [f64; 3], t: [[f64; 3]; 3]) -> bool {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let (e1, e2) = (sub(t[1], t[0]), sub(t[2], t[0]));
    let p = cross(d, e2);
    let det = dot(e1, p);
    if det.abs() < 1e-14 {
        return false;
    }
    let s = sub(o, t[0]);
    let u = dot(s, p) / det;
    let q = cross(s, e1);
    let v = dot(d, q) / det;
    u >= 0.0 && v >= 0.0 && u + v <= 1.0
}

// Crossings of a whole line with a dense tessellation of the band, from the law.
fn brute_hits(pipe: &RollingPipe, o: [f64; 3], d: [f64; 3]) -> usize {
    const COLUMNS: i64 = 360;
    const LINES: i64 = 24;
    let at = |i: i64, j: i64| {
        let p = pipe.point(
            &Rational::from_ratio(i, COLUMNS / 4).unwrap(),
            &Rational::from_ratio(j, LINES).unwrap(),
        );
        [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()]
    };
    let mut grid = Vec::new();
    for i in 0..=COLUMNS {
        for j in 0..=LINES {
            grid.push(at(i, j));
        }
    }
    let node = |i: i64, j: i64| grid[(i * (LINES + 1) + j) as usize];
    let mut hits = 0;
    for i in 0..COLUMNS {
        for j in 0..LINES {
            let (a, b, c, e) = (
                node(i, j),
                node(i + 1, j),
                node(i, j + 1),
                node(i + 1, j + 1),
            );
            hits += usize::from(triangle_hit(o, d, [a, b, e]));
            hits += usize::from(triangle_hit(o, d, [a, e, c]));
        }
    }
    hits
}

#[test]
fn line_hits_find_every_crossing_of_the_band() {
    let lines: [([&str; 3], [&str; 3]); 8] = [
        (["1/3", "1/7", "4"], ["1", "1/9", "1/7"]),
        (["1/3", "1/7", "4"], ["1/11", "1", "1/7"]),
        (["1/3", "1/7", "4"], ["1/13", "1/9", "1"]),
        (["27/5", "1/17", "-71/20"], ["84", "1/9", "131"]),
        (["-23/5", "2/3", "9"], ["1", "1/3", "1/2"]),
        (["13/3", "1/5", "5"], ["-1", "1/2", "1/4"]),
        (["21/4", "1/19", "-7/2"], ["1", "1/5", "1"]),
        (["-11/2", "3/7", "3"], ["2/3", "1", "-1/3"]),
    ];
    for case in rational_cases() {
        let pipe = build(&case);
        for (origin, direction) in lines {
            let (o, d) = (point(origin), vector(direction));
            let found = pipe.line_hits(&o, &d, &mut decider()).unwrap().unwrap();
            let f = |p: &Point3| [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()];
            let v = [d.x.to_f64(), d.y.to_f64(), d.z.to_f64()];
            assert_eq!(
                found.len(),
                brute_hits(&pipe, f(&o), v),
                "{} line {origin:?} {direction:?}",
                case.name
            );
        }
    }
}

// Each line sits on a special place of the equations, where a grid line of the
// tessellation can sit too, so the count comes from a small shift of the line.
// Transverse crossings keep their count under it.
#[test]
fn line_hits_hold_on_the_special_lines_of_the_law() {
    let lines: [([&str; 3], [&str; 3]); 8] = [
        // Along the axis, in the section plane of t = 0 and of t = pi.
        (["26/5", "0", "0"], ["0", "0", "1"]),
        (["24/5", "0", "0"], ["0", "0", "1"]),
        (["-27/5", "0", "0"], ["0", "0", "1"]),
        (["-9/2", "0", "0"], ["0", "0", "1"]),
        // Along the plane normal, where the foot equation loses the parameter.
        (["23/5", "1/2", "0"], ["3/5", "0", "4/5"]),
        (["1/2", "24/5", "0"], ["3/5", "0", "4/5"]),
        (["-1/2", "-21/4", "0"], ["3/5", "0", "4/5"]),
        // Across the section plane of t = 0.
        (["24/5", "-9", "5/2"], ["0", "1", "1/4"]),
    ];
    let mut reached = [0; 8];
    for case in rational_cases() {
        let pipe = build(&case);
        for (index, (origin, direction)) in lines.into_iter().enumerate() {
            let (o, d) = (point(origin), vector(direction));
            let found = pipe.line_hits(&o, &d, &mut decider()).unwrap().unwrap();
            let f = |p: &Point3| [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()];
            let (a, v) = (f(&o), [d.x.to_f64(), d.y.to_f64(), d.z.to_f64()]);
            let shifted = brute_hits(&pipe, [a[0] + 1e-4, a[1] + 2e-4, a[2]], v);
            assert_eq!(found.len(), shifted, "{} {origin:?}", case.name);
            reached[index] += found.len();
            let radius = pipe.radius().to_f64();
            for hit in &found {
                let dec = decider();
                let near = |value: &Real| dec.enclose_at(value, 64).unwrap().lo().to_f64_nearest();
                let at = [near(&hit.x), near(&hit.y), near(&hit.z)];
                let gap = distance_to_ellipse(at, &pipe) - radius;
                assert!(gap.abs() < 1e-5, "{} {origin:?}: {gap}", case.name);
            }
        }
    }
    assert!(reached.iter().all(|&count| count > 0), "{reached:?}");
}

// A line inside a section plane has no parameter from the foot equation. Lines
// through a known band point fix its crossing exactly.
#[test]
fn a_line_inside_a_section_plane_finds_its_band_point() {
    for case in rational_cases() {
        let pipe = build(&case);
        for u in ["0", "2"] {
            for v in ["1/3", "1/2", "4/5"] {
                let (u, v) = (q(u), q(v));
                let p = pipe.point(&u, &v);
                let ray = p.minus(&pipe.centre(&u));
                let turned = Vector3::new(
                    &ray.x + &(&ray.z * &real("1/5")),
                    Real::zero(),
                    &ray.z - &(&ray.x * &real("1/5")),
                );
                for direction in [ray, turned] {
                    let found = pipe.line_hits(&p, &direction, &mut decider());
                    let found = found.unwrap().unwrap();
                    let mut dec = decider();
                    assert!(
                        found.iter().any(|hit| same_point(hit, &p, &mut dec)),
                        "{} {u:?} {v:?}",
                        case.name
                    );
                    let a = [p.x.to_f64() + 1e-4, p.y.to_f64() + 2e-4, p.z.to_f64()];
                    let d = [
                        direction.x.to_f64(),
                        direction.y.to_f64(),
                        direction.z.to_f64(),
                    ];
                    assert_eq!(found.len(), brute_hits(&pipe, a, d), "{}", case.name);
                }
            }
        }
    }
}

#[test]
fn a_line_tangent_to_the_band_refuses() {
    for case in rational_cases() {
        let pipe = build(&case);
        let half = q("1/2");
        let p = pipe.point(&Rational::zero(), &half);
        let ray = p.minus(&pipe.centre(&Rational::zero()));
        let across = vector(["0", "1", "0"]);
        let around = Vector3::new(ray.z.clone(), Real::zero(), -ray.x.clone());
        for direction in [across, around] {
            let answer = pipe.line_hits(&p, &direction, &mut decider());
            let refused = answer.expect_err(case.name);
            assert_eq!(
                refused.code, "intersect.parametric-surface",
                "{}",
                case.name
            );
        }
    }
}
