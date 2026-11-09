//! The exact rolling-ball pipe between a plane and an oblique cone.
//!
//! On each span the carrier's distance identity has bidegree (12, 4) once its
//! denominators are cleared, so exact vanishing on a 13 by 5 grid proves it
//! over the whole span, not only at the samples.

use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real, Sign};
use catallax_shapes::{BSplineSurface, Cone, Curve, Ellipse, Plane, Point3, RollingPipe, Vector3};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 200_000_000,
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
fn cone(apex: [&str; 3], axis: [&str; 3], slope: &str) -> Cone {
    Cone {
        apex: point(apex),
        axis: vector(axis),
        slope: real(slope),
        ref_dir: vector(["1", "0", "0"]),
    }
}

struct Case {
    name: &'static str,
    plane: Plane,
    cone: Cone,
    same: [bool; 2],
    radius: &'static str,
    convex: bool,
}
// Apex (0, 0, 8) opening downwards with slope 3/4, so radius 6 at the base,
// cut by a plane through (0, 0, 3) with a rational unit normal. Every number
// of the law is rational.
fn cases() -> Vec<Case> {
    let cut = || plane(["0", "0", "3"], ["3/5", "0", "4/5"]);
    let frustum = || cone(["0", "0", "8"], ["0", "0", "-1"], "3/4");
    vec![
        Case {
            name: "rim of a cut cone",
            plane: cut(),
            cone: frustum(),
            same: [true, true],
            radius: "1/2",
            convex: true,
        },
        Case {
            name: "cone standing on an oblique face",
            plane: cut(),
            cone: frustum(),
            same: [false, true],
            radius: "1/2",
            convex: false,
        },
        Case {
            name: "conical bore through an oblique face",
            plane: cut(),
            cone: frustum(),
            same: [true, false],
            radius: "1/3",
            convex: true,
        },
    ]
}
fn build(case: &Case) -> RollingPipe {
    RollingPipe::conical(
        &case.plane,
        &case.cone,
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
        for i in 0..=12 {
            for j in 0..=4 {
                samples.push((
                    Rational::from_i64(span) + Rational::from_ratio(i, 12).unwrap(),
                    Rational::from_ratio(j, 4).unwrap(),
                ));
            }
        }
    }
    samples
}
// Zero on the cone: radial distance squared against slope squared times height
// squared, on the nappe the axis points into.
fn on_cone(cone: &Cone, p: &Point3, dec: &mut Decider) -> bool {
    let d = p.minus(&cone.apex);
    let axial = d.dot(&cone.axis);
    let radial = d.cross(&cone.axis).norm_squared();
    matches!(dec.sign(&axial), Outcome::Proven(_))
        && zero(&radial - &(&cone.slope.square() * &axial.square()), dec)
}

#[test]
fn the_rational_carrier_is_the_rolling_ball_definition_on_every_span() {
    for case in cases() {
        let pipe = build(&case);
        assert!(pipe.is_conical());
        assert!(pipe.carrier().has_rational_definition(), "{}", case.name);
        assert_eq!(pipe.carrier().degrees(), [6, 2]);
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
fn both_rails_lie_on_their_supports_and_centres_on_both_offsets() {
    for case in cases() {
        let pipe = build(&case);
        let mut dec = decider();
        let rail = pipe.plane_rail(&mut dec).unwrap();
        let spine = pipe.spine(&mut dec).unwrap();
        assert!(pipe.wall_rail(&mut dec).is_err());
        let Curve::BSpline(wall) = pipe.wall_curve(&mut dec).unwrap() else {
            panic!("a cone's wall rail is the carrier's last row");
        };
        assert_eq!(wall.degree(), 6);
        let normal = case
            .plane
            .normal
            .scaled(&(&Real::one() / &case.plane.normal.norm_squared().sqrt()));
        let height = |p: &Point3| p.minus(&case.plane.origin).dot(&normal);
        let radius = real(case.radius);
        for (u, _) in grid() {
            let on_plane = evaluate(pipe.carrier(), &u, &Rational::zero());
            let on_wall = evaluate(pipe.carrier(), &u, &Rational::one());
            assert!(zero(height(&on_plane), &mut dec), "{}", case.name);
            assert!(zero(
                on_plane.minus(&rail.center).dot(&rail.normal),
                &mut dec
            ));
            assert!(on_cone(&case.cone, &on_wall, &mut dec), "{}", case.name);
            assert!(
                pipe.on_boundary(&on_plane, &mut dec).unwrap(),
                "{}",
                case.name
            );
            assert!(
                pipe.on_boundary(&on_wall, &mut dec).unwrap(),
                "{}",
                case.name
            );
            let centre = pipe.centre(&u);
            assert!(zero(
                centre.minus(&spine.center).dot(&spine.normal),
                &mut dec
            ));
            // The centre is one radius off the plane and off the cone, so its
            // contacts are the feet of those distances.
            assert!(zero(height(&centre).square() - radius.square(), &mut dec));
            let ray = on_wall.minus(&centre);
            assert!(zero(ray.norm_squared() - radius.square(), &mut dec));
            assert!(zero(
                pipe.normal_at(&on_wall, &mut dec)
                    .unwrap()
                    .cross(&ray)
                    .norm_squared(),
                &mut dec
            ));
        }
        // The plane rail holds its row, and the spine its centres.
        let on = |ellipse: &Ellipse, p: &Point3, dec: &mut Decider| {
            let minor_dir = ellipse.normal.cross(&ellipse.major_dir);
            let d = p.minus(&ellipse.center);
            let x = d.dot(&ellipse.major_dir).square() / ellipse.major_dir.norm_squared();
            let y = d.dot(&minor_dir).square() / minor_dir.norm_squared();
            zero(
                &(&x / &ellipse.major_sq) + &(&y / &ellipse.minor_sq) - Real::one(),
                dec,
            )
        };
        for k in 0..8 {
            let u = Rational::from_ratio(k, 2).unwrap();
            let row = evaluate(pipe.carrier(), &u, &Rational::zero());
            assert!(on(&rail, &row, &mut dec), "{} plane rail at {k}", case.name);
            assert!(
                on(&spine, &pipe.centre(&u), &mut dec),
                "{} spine at {k}",
                case.name
            );
        }
    }
}

#[test]
fn the_seam_closes_exactly() {
    for case in cases() {
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

// The plane's section of the cone, found by hand: its vertices in y = 0 are
// (60/7, 0, -24/7) and (-12/5, 0, 24/5), and its half-width there is sqrt(144/7).
#[test]
fn the_input_edge_is_the_planes_section_of_the_cone() {
    let mut dec = decider();
    let edge = Ellipse::new(
        point(["108/35", "0", "24/35"]),
        vector(["3/5", "0", "4/5"]),
        vector(["4", "0", "-3"]),
        real("2304/49"),
        real("144/7"),
        &mut dec,
    )
    .unwrap();
    for case in cases() {
        let pipe = build(&case);
        assert!(pipe.is_edge(&edge, &mut dec).unwrap(), "{}", case.name);
    }
    let shifted = Ellipse {
        center: point(["108/35", "0", "25/35"]),
        ..edge
    };
    assert!(!build(&cases()[0]).is_edge(&shifted, &mut dec).unwrap());
}

#[test]
fn points_off_the_band_report_the_tube_side() {
    let pipe = build(&cases()[0]);
    let mut dec = decider();
    // The ball centre at u = 0 is inside the tube, and the cone's axis well
    // below the cut is outside it.
    let centre = pipe.centre(&Rational::zero());
    assert_eq!(pipe.side_at(&centre, &mut dec).unwrap(), Sign::Positive);
    assert_eq!(
        pipe.side_at(&point(["0", "0", "-1"]), &mut dec).unwrap(),
        Sign::Negative
    );
    let extrema = pipe.axis_extrema(&mut dec).unwrap();
    assert!(!extrema.is_empty());
    for p in &extrema {
        assert!(matches!(pipe.side_at(p, &mut dec), Ok(Sign::Zero) | Err(_)));
    }
}

#[test]
fn an_open_section_and_a_folding_ball_refuse() {
    let mut dec = decider();
    let frustum = cone(["0", "0", "8"], ["0", "0", "-1"], "3/4");
    // Tilted past the generator: a hyperbola, not a closed ring.
    let steep = plane(["0", "0", "3"], ["4/5", "0", "3/5"]);
    assert!(
        RollingPipe::conical(&steep, &frustum, [true, true], real("1/2"), true, &mut dec).is_err()
    );
    let cut = plane(["0", "0", "3"], ["3/5", "0", "4/5"]);
    assert!(RollingPipe::conical(&cut, &frustum, [true, true], real("4"), true, &mut dec).is_err());
}

// The corpus case: a frustum of radii 10 and 4 over height 24, cut through
// (0, 0, 8) by the plane with normal (336, 0, 527) / 625, rounded by 1/2. Its
// cone has the irrational normal scale sqrt(17) / 4.
#[test]
fn the_corpus_cut_cone_builds_with_an_irrational_scale() {
    let mut dec = decider();
    let pipe = RollingPipe::conical(
        &plane(["0", "0", "8"], ["336/625", "0", "527/625"]),
        &cone(["0", "0", "40"], ["0", "0", "-1"], "1/4"),
        [true, true],
        real("1/2"),
        true,
        &mut dec,
    )
    .unwrap();
    assert!(!pipe.carrier().has_rational_definition());
    for (u, v) in [("0", "0"), ("1/3", "1/2"), ("7/4", "1"), ("3", "1/4")] {
        let [u, v] = [q(u), q(v)];
        let on_carrier = evaluate(pipe.carrier(), &u, &v);
        assert!(same_point(&on_carrier, &pipe.point(&u, &v), &mut dec));
    }
}

// A line through a band point finds it among its crossings, in line order,
// and every crossing lies on the band. The cases' laws are rational but
// conical, and the corpus law irrational, so all take the algebraic path.
#[test]
fn line_hits_find_band_points_on_algebraic_laws() {
    let corpus = RollingPipe::conical(
        &plane(["0", "0", "8"], ["336/625", "0", "527/625"]),
        &cone(["0", "0", "40"], ["0", "0", "-1"], "1/4"),
        [true, true],
        real("1/2"),
        true,
        &mut decider(),
    )
    .unwrap();
    let pipes: Vec<_> = cases().iter().map(build).chain([corpus]).collect();
    let directions = [
        ["1", "1/3", "1/7"],
        ["1/5", "-1", "2/9"],
        ["-1/4", "1/6", "1"],
    ];
    for (index, pipe) in pipes.iter().enumerate() {
        let mut dec = decider();
        for (u, v) in [("1/3", "1/2"), ("9/4", "1/5")] {
            let at = pipe.point(&q(u), &q(v));
            for direction in directions {
                let line = vector(direction);
                let mut dec = decider();
                let hits = pipe
                    .line_hits(&at, &line, &mut dec)
                    .unwrap_or_else(|e| panic!("pipe {index} at ({u}, {v}): {e:?}"))
                    .unwrap();
                assert!(
                    hits.iter().any(|hit| same_point(hit, &at, &mut dec)),
                    "pipe {index} at ({u}, {v}) along {direction:?}"
                );
                for pair in hits.windows(2) {
                    let step = pair[1].minus(&pair[0]).dot(&line);
                    assert_eq!(math_sign(step, &mut dec), Sign::Positive);
                }
                for hit in &hits {
                    // A transverse crossing: the tube's sides differ just before
                    // and just after it.
                    let nudge = line.scaled(&real("1/10000"));
                    let before = pipe.side_at(&hit.plus(&nudge.scaled(&real("-1"))), &mut dec);
                    let after = pipe.side_at(&hit.plus(&nudge), &mut dec);
                    let sides = [before.unwrap(), after.unwrap()];
                    assert!(
                        sides.contains(&Sign::Positive) && sides.contains(&Sign::Negative),
                        "pipe {index} at ({u}, {v}) along {direction:?}"
                    );
                }
            }
        }
        // The cone's axis passes through the ring's hole.
        let axis = pipe
            .line_hits(&point(["0", "0", "0"]), &vector(["0", "0", "1"]), &mut dec)
            .unwrap()
            .unwrap();
        assert!(axis.is_empty(), "pipe {index}: {}", axis.len());
    }
}

fn math_sign(value: Real, dec: &mut Decider) -> Sign {
    match dec.sign(&value) {
        catallax_numbers::Outcome::Proven(catallax_numbers::StrictSign::Positive) => Sign::Positive,
        catallax_numbers::Outcome::Proven(catallax_numbers::StrictSign::Negative) => Sign::Negative,
        other => panic!("{other:?}"),
    }
}

#[test]
fn ball_sweep_bounds_cover_the_conical_band_without_control_enclosures() {
    use catallax_numbers::with_operation_budget;
    let pipe = RollingPipe::conical(
        &plane(["0", "0", "8"], ["336/625", "0", "527/625"]),
        &cone(["0", "0", "40"], ["0", "0", "-1"], "1/4"),
        [true, true],
        real("1/2"),
        true,
        &mut decider(),
    )
    .unwrap();
    let (hull, hull_work) = with_operation_budget(1_000_000, || {
        catallax_numbers::reuse_operation_enclosures();
        pipe.enclosed_controls(&decider())
    });
    if hull_work.exhausted {
        assert_eq!(hull.unwrap_err().code, "numbers.operation-budget");
    } else {
        assert!(hull.is_ok());
    }
    for axis in 0..3 {
        let (range, swept_work) = with_operation_budget(1_000_000, || {
            catallax_numbers::reuse_operation_enclosures();
            pipe.swept_axis_range(axis)
                .map(|value| decider().enclose_at(&value, 256).unwrap())
        });
        assert!(!swept_work.exhausted);
        assert!(swept_work.used < hull_work.used);
        let low = range[0].lo();
        let high = range[1].hi();
        for u in ["0", "1/3", "1", "7/4", "2", "3", "4"] {
            for v in ["0", "1/3", "2/3", "1"] {
                let point = evaluate(pipe.carrier(), &q(u), &q(v));
                let coordinate = [&point.x, &point.y, &point.z][axis];
                let mut dec = decider();
                assert_ne!(
                    math_sign(coordinate - &Real::from(low.clone()), &mut dec),
                    Sign::Negative
                );
                assert_ne!(
                    math_sign(&Real::from(high.clone()) - coordinate, &mut dec),
                    Sign::Negative
                );
            }
        }
    }
}
