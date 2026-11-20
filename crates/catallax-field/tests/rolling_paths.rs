//! 4.2a step 1: the general rolling-ball path routine.
//!
//! The line and circle pieces are derived from their supports alone. Their
//! membership must equal today's compiled field answer at every sampled
//! point, bit for bit, and their distance ranges must contain the exact
//! distance to the fillet face from an independent section-plane formula.
//! The same pieces run through the general panel cover give the same answers.
use catallax_field::rolling_region::elementary::{Flat, Ring, Round, Strip};
use catallax_field::rolling_region::path::{self, Cover, Frame, Path};
use catallax_field::rolling_region::{Region, Relation, query_batch};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{
    BudgetPolicy, Decider, OperationAllowance, Outcome, Rational, Real, Refusal, StrictSign,
    with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_shapes::{Point3, RollingEnvelopeDefinition, Vector3};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn r(s: &str) -> Real {
    q(s).into()
}
fn p3(x: &Rational, y: &Rational, z: &Rational) -> Point3 {
    Point3::new(x.clone().into(), y.clone().into(), z.clone().into())
}
fn v(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(x.into(), y.into(), z.into())
}
fn flat(point: [&str; 3], normal: [i64; 3]) -> Flat {
    Flat {
        point: Point3::new(r(point[0]), r(point[1]), r(point[2])),
        normal: v(normal[0], normal[1], normal[2]),
    }
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}
fn plan(json: &str) -> Plan {
    Plan::compile_checked(&Recipe::parse(json).unwrap(), &mut decider()).unwrap()
}
/// The recipe without its fillet step: the solid the fillet region trims or fills.
fn base(json: &str) -> Plan {
    let mut value: serde_json::Value = serde_json::from_str(json).unwrap();
    let steps = value["steps"].as_array_mut().unwrap();
    let fillet = steps.pop().unwrap();
    assert_eq!(fillet["op"], "fillet");
    value["result"] = fillet["body"].clone();
    plan(&value.to_string())
}
fn grid(values: &[&str]) -> Vec<Rational> {
    values.iter().map(|s| q(s)).collect()
}

/// Whether the face of a convex fillet is subtracted (true) or a concave one added.
fn expected(location: Location, convex: bool) -> Relation {
    match (location, convex) {
        (Location::Boundary, _) => Relation::Boundary,
        (Location::Inside, true) | (Location::Outside, false) => Relation::Exterior,
        (Location::Outside, true) | (Location::Inside, false) => Relation::Interior,
        (other, _) => panic!("unresolved compiled answer {other:?}"),
    }
}

/// Exact distance from a section point to the arc of centre `c`, radius `rad`,
/// between feet `a` and `b` (all in one section plane), plus the squared
/// offset along the path. Independent of the routine under test.
fn arc_distance(
    p: [Rational; 2],
    c: [Rational; 2],
    a: [Rational; 2],
    b: [Rational; 2],
    rad: &Rational,
    axial: Rational,
) -> Real {
    let sub = |x: &[Rational; 2], y: &[Rational; 2]| [&x[0] - &y[0], &x[1] - &y[1]];
    let cross = |x: &[Rational; 2], y: &[Rational; 2]| &x[0] * &y[1] - &x[1] * &y[0];
    let norm = |x: &[Rational; 2]| &x[0] * &x[0] + &x[1] * &x[1];
    let w = sub(&p, &c);
    let (da, db) = (sub(&a, &c), sub(&b, &c));
    let turn = cross(&da, &db);
    let inside = cross(&da, &w) * turn.clone() >= Rational::zero()
        && cross(&w, &db) * turn >= Rational::zero()
        && norm(&w) > Rational::zero();
    let root = |x: Rational| Real::from(x).sqrt();
    let planar = if inside {
        let d = &root(norm(&w)) - &Real::from(rad.clone());
        d.square()
    } else {
        let ea = norm(&sub(&p, &a));
        let eb = norm(&sub(&p, &b));
        Real::from(ea.min(eb))
    };
    (&planar + &Real::from(axial)).sqrt()
}

fn contains(range: &path::Range, exact: &Real, dec: &mut Decider) {
    let above = |x: Real, dec: &mut Decider| {
        matches!(
            dec.sign(&x),
            Outcome::Proven(StrictSign::Positive) | Outcome::ProvenCoincident(_)
        )
    };
    assert!(
        above(exact - &Real::from(range.lower.clone()), dec),
        "lower {:?} above exact {exact:?}",
        range.lower
    );
    let upper = range.upper.clone().expect("a certified witness");
    assert!(
        above(&Real::from(upper) - exact, dec),
        "upper {:?} below exact {exact:?}",
        range.upper
    );
}

/// The same path without its closed forms: every answer comes from the
/// general whole-interval panel cover.
struct Panels<'a, P>(&'a P);
impl<P: Path> Path for Panels<'_, P> {
    type Cover = P::Cover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        self.0.local(point)
    }
    fn walls(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        self.0.walls(local, dec)
    }
    fn walls_are_faces(&self) -> bool {
        self.0.walls_are_faces()
    }
    fn cover(&self, dec: &mut Decider) -> Result<P::Cover, Refusal> {
        self.0.cover(dec)
    }
    fn frames(&self, dec: &mut Decider) -> Result<P::Cover, Refusal> {
        self.0.frames(dec)
    }
}

struct Case<'a> {
    recipe: &'a str,
    convex: bool,
    // World point from section coordinates (s, t) and a station w.
    world: fn(&Rational, &Rational, &Rational) -> Point3,
    s: Vec<Rational>,
    t: Vec<Rational>,
    w: Vec<Rational>,
    centre: [&'a str; 2],
    feet: [[&'a str; 2]; 2],
    radius: &'a str,
    // Run the general panel cover on every `stride`-th grid point.
    stride: usize,
}

fn check<P: Path>(case: &Case<'_>, piece: &P) {
    let plan = plan(case.recipe);
    let base = base(case.recipe);
    let mut workspace = Workspace::default();
    let mut dec = decider();
    let rad = q(case.radius);
    let c = case.centre.map(q);
    let feet = case.feet.map(|f| f.map(q));
    let tiny = q("1/1000000000000");
    let coarse = q("1/1000000");
    let panels = Panels(piece);
    let mut visited = 0usize;
    let mut compared = [0usize; 3];
    for w in &case.w {
        for s in &case.s {
            for t in &case.t {
                let point = (case.world)(s, t, w);
                let rational = point.as_rational().unwrap();
                // Compare where the fillet alone decides: strictly inside the
                // solid it trims, or strictly outside the solid it fills.
                let before = base
                    .proven_checked(&rational, &mut workspace, &mut dec)
                    .unwrap()
                    .location;
                if before
                    != if case.convex {
                        Location::Inside
                    } else {
                        Location::Outside
                    }
                {
                    continue;
                }
                let compiled = plan
                    .proven_checked(&rational, &mut workspace, &mut dec)
                    .unwrap()
                    .location;
                let want = expected(compiled, case.convex);
                let got = path::classify(piece, &point, &mut dec).unwrap();
                assert_eq!(got, want, "at ({s}, {t}, {w}): compiled {compiled:?}");
                compared[match got {
                    Relation::Interior => 0,
                    Relation::Exterior => 1,
                    Relation::Boundary => 2,
                }] += 1;
                let exact = arc_distance(
                    [s.clone(), t.clone()],
                    c.clone(),
                    feet[0].clone(),
                    feet[1].clone(),
                    &rad,
                    Rational::zero(),
                );
                let range = path::distance(piece, &point, &tiny, &mut dec).unwrap();
                assert!(range.converged, "closed form at ({s}, {t}): {range:?}");
                contains(&range, &exact, &mut dec);
                // On every `stride`-th point, the general cover agrees off the
                // boundary and brackets the same exact distance.
                visited += 1;
                if visited.is_multiple_of(case.stride) {
                    if got != Relation::Boundary {
                        assert_eq!(path::classify(&panels, &point, &mut dec).unwrap(), got);
                    }
                    let general = path::distance(&panels, &point, &coarse, &mut dec).unwrap();
                    contains(&general, &exact, &mut dec);
                }
            }
        }
    }
    // Every grid must exercise all three answers.
    assert!(compared.iter().all(|n| *n > 0), "{compared:?}");
}

#[test]
fn box_edge_strip_matches_the_compiled_fillet() {
    let mut dec = decider();
    let piece = Strip::new(
        [
            flat(["0", "0", "10"], [0, 0, -1]),
            flat(["20", "0", "0"], [-1, 0, 0]),
        ],
        r("1"),
        [q("0"), q("20")],
        &mut dec,
    )
    .unwrap();
    check(
        &Case {
            recipe: r#"{"format":"catallax-recipe/4","name":"box-single","result":"round","steps":[
                {"id":"b","op":"block","size":["20","20","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                {"id":"round","op":"fillet","body":"b","radius":"1","edges":[{"between":[{"step":"b","role":"end"},{"step":"b","role":{"side":{"segment":"x_max"}}}]}]}]}"#,
            convex: true,
            world: |s, t, w| p3(s, w, t),
            s: grid(&[
                "35/2", "18", "37/2", "19", "193/10", "39/2", "98/5", "199/10",
            ]),
            t: grid(&[
                "15/2", "8", "17/2", "9", "93/10", "19/2", "48/5", "49/5", "99/10",
            ]),
            w: grid(&["5", "21/2"]),
            centre: ["19", "9"],
            feet: [["19", "10"], ["20", "9"]],
            radius: "1",
            stride: 3,
        },
        &piece,
    );
    // Exactly on the blend: centre (19, 9) plus (3/5, 4/5).
    assert_eq!(
        path::classify(&piece, &p3(&q("98/5"), &q("4"), &q("49/5")), &mut dec).unwrap(),
        Relation::Boundary
    );
}

#[test]
fn oblique_strip_matches_the_compiled_fillet() {
    // Planes x = 10 and 3x + 4y = 50 meet at 126.87 degrees; with r = 1 the
    // centre is (9, 9/2) and the feet (10, 9/2) and (48/5, 53/10).
    let mut dec = decider();
    let piece = Strip::new(
        [
            flat(["10", "0", "0"], [-1, 0, 0]),
            flat(["10", "5", "0"], [-3, -4, 0]),
        ],
        r("1"),
        [q("0"), q("5/2")],
        &mut dec,
    )
    .unwrap();
    check(
        &Case {
            recipe: r#"{"format":"catallax-recipe/4","name":"oblique strip","result":"round","steps":[
                {"id":"section","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
                 "loops":[[{"kind":"line","id":"e0","to":["0","0"]},{"kind":"line","id":"e1","to":["10","0"]},
                           {"kind":"line","id":"e2","to":["10","5"]},{"kind":"line","id":"e3","to":["6","8"]},
                           {"kind":"line","id":"e4","to":["0","8"]}]]},
                {"id":"w","op":"extrude","profile":"section","distance":"10"},
                {"id":"round","op":"fillet","body":"w","radius":"1","edges":[{"between":[{"step":"w","role":{"side":{"segment":"e2"}}},{"step":"w","role":{"side":{"segment":"e3"}}}]}]}]}"#,
            convex: true,
            world: |s, t, w| p3(s, t, w),
            s: grid(&["17/2", "9", "93/10", "19/2", "97/10", "49/5", "99/10"]),
            t: grid(&["4", "9/2", "47/10", "5", "51/10", "21/4"]),
            w: grid(&["3", "15/2"]),
            centre: ["9", "9/2"],
            feet: [["10", "9/2"], ["48/5", "53/10"]],
            radius: "1",
            stride: 3,
        },
        &piece,
    );
}

#[test]
fn concave_strip_matches_the_compiled_fillet() {
    let mut dec = decider();
    let piece = Strip::new(
        [
            flat(["0", "0", "5"], [0, 0, 1]),
            flat(["0", "5", "0"], [0, 1, 0]),
        ],
        r("1"),
        [q("-20"), q("0")],
        &mut dec,
    )
    .unwrap();
    check(
        &Case {
            recipe: r#"{"format":"catallax-recipe/4","name":"l-concave","result":"round","steps":[
                {"id":"base","op":"block","size":["20","20","5"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                {"id":"wall","op":"block","size":["20","5","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                {"id":"j","op":"boolean","kind":"union","target":"base","tool":"wall"},
                {"id":"round","op":"fillet","body":"j","radius":"1","edges":[{"between":[{"step":"base","role":"end"},{"step":"wall","role":{"side":{"segment":"y_max"}}}]}]}]}"#,
            convex: false,
            world: |s, t, w| p3(w, s, t),
            s: grid(&["51/10", "26/5", "27/5", "11/2", "6", "31/5", "13/2", "7"]),
            t: grid(&["51/10", "26/5", "11/2", "6", "31/5", "13/2", "7"]),
            w: grid(&["4", "13"]),
            centre: ["6", "6"],
            feet: [["6", "5"], ["5", "6"]],
            radius: "1",
            stride: 3,
        },
        &piece,
    );
}

fn ring_convex() -> Ring {
    Ring::new(
        flat(["0", "0", "12"], [0, 0, -1]),
        Round {
            point: Point3::origin(),
            axis: v(0, 0, 1),
            radius: r("8"),
            outside: false,
        },
        r("2"),
        &mut decider(),
    )
    .unwrap()
}

// Meridian coordinates (rho, z) at a rational unit direction.
fn meridian(direction: [&'static str; 2]) -> fn(&Rational, &Rational, &Rational) -> Point3 {
    match direction {
        ["1", "0"] => |s, t, _| p3(s, &Rational::zero(), t),
        ["3/5", "4/5"] => |s, t, _| p3(&(s * &q("3/5")), &(s * &q("4/5")), t),
        _ => |s, t, _| p3(&Rational::zero(), &-s, t),
    }
}

#[test]
fn convex_ring_matches_the_compiled_fillet() {
    let piece = ring_convex();
    for direction in [["1", "0"], ["3/5", "4/5"]] {
        check(
            &Case {
                recipe: r#"{"format":"catallax-recipe/3","name":"convex cylinder","steps":[
                    {"id":"base","op":"cylinder","radius":"8","height":"12","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                    {"id":"round","op":"fillet","body":"base","radius":"2","edges":[{"between":[{"step":"base","role":"end"},{"step":"base","role":{"side":{"segment":"wall"}}}]}]}],"result":"round"}"#,
                convex: true,
                world: meridian(direction),
                s: grid(&["5", "23/4", "6", "13/2", "7", "36/5", "15/2", "79/10"]),
                t: grid(&["9", "19/2", "10", "21/2", "11", "58/5", "23/2", "119/10"]),
                w: grid(&["0"]),
                centre: ["6", "10"],
                feet: [["6", "12"], ["8", "10"]],
                radius: "2",
                stride: 6,
            },
            &piece,
        );
    }
}

#[test]
fn concave_boss_ring_matches_the_compiled_fillet() {
    let piece = Ring::new(
        flat(["0", "0", "5"], [0, 0, 1]),
        Round {
            point: Point3::origin(),
            axis: v(0, 0, 1),
            radius: r("5"),
            outside: true,
        },
        r("1"),
        &mut decider(),
    )
    .unwrap();
    for direction in [["1", "0"], ["3/5", "4/5"]] {
        check(
            &Case {
                recipe: r#"{"format":"catallax-recipe/4","name":"boss-rim","result":"round","steps":[
                    {"id":"plate","op":"block","size":["30","30","5"],"pose":{"rotation":["1","0","0","0"],"translation":["-15","-15","0"]}},
                    {"id":"c","op":"cylinder","radius":"5","height":"10","pose":{"rotation":["1","0","0","0"],"translation":["0","0","5"]}},
                    {"id":"j","op":"boolean","kind":"union","target":"plate","tool":"c"},
                    {"id":"round","op":"fillet","body":"j","radius":"1","edges":[{"between":[{"step":"plate","role":"end"},{"step":"c","role":{"side":{"segment":"wall"}}}]}]}]}"#,
                convex: false,
                world: meridian(direction),
                s: grid(&["51/10", "21/4", "27/5", "11/2", "28/5", "6", "13/2", "7"]),
                t: grid(&["51/10", "26/5", "21/4", "27/5", "11/2", "6", "13/2", "7"]),
                w: grid(&["0"]),
                centre: ["6", "6"],
                feet: [["6", "5"], ["5", "6"]],
                radius: "1",
                stride: 6,
            },
            &piece,
        );
    }
}

#[test]
fn derivation_refusals_are_named() {
    let mut dec = decider();
    let parallel = Strip::new(
        [
            flat(["0", "0", "0"], [0, 0, 1]),
            flat(["0", "0", "5"], [0, 0, -1]),
        ],
        r("1"),
        [q("0"), q("1")],
        &mut dec,
    );
    assert_eq!(parallel.unwrap_err().code, "field.rolling-path-reach");
    // A convex ring of radius 2 on a cylinder of radius 4 has C = r: a horn tube.
    let horn = Ring::new(
        flat(["0", "0", "12"], [0, 0, -1]),
        Round {
            point: Point3::origin(),
            axis: v(0, 0, 1),
            radius: r("4"),
            outside: false,
        },
        r("2"),
        &mut dec,
    );
    assert_eq!(horn.unwrap_err().code, "field.rolling-path-reach");
    let tilted = Ring::new(
        flat(["0", "0", "12"], [0, 1, -1]),
        Round {
            point: Point3::origin(),
            axis: v(0, 0, 1),
            radius: r("8"),
            outside: false,
        },
        r("2"),
        &mut dec,
    );
    assert_eq!(tilted.unwrap_err().code, "field.rolling-path-reach");
    // The general cover refuses an exact face point rather than guessing.
    let on_face = p3(&q("36/5"), &Rational::zero(), &q("58/5"));
    assert_eq!(
        path::classify(&Panels(&ring_convex()), &on_face, &mut dec)
            .unwrap_err()
            .code,
        "field.rolling-region-unresolved"
    );
}

#[test]
fn batch_answers_are_identical_at_every_worker_count() {
    let piece = ring_convex();
    let points: Vec<Point3> = [
        ["7", "0", "11"],
        ["38/5", "0", "58/5"],
        ["36/5", "0", "58/5"],
        ["21/5", "28/5", "11"],
        ["5", "0", "9"],
    ]
    .iter()
    .map(|p| p3(&q(p[0]), &q(p[1]), &q(p[2])))
    .collect();
    let width = q("1/1000000");
    let mut first = None;
    for workers in [Some(1), Some(2), Some(4), None] {
        let (batch, _) = with_operation_budget(10_000_000, || {
            query_batch(
                &piece,
                &points,
                &width,
                workers,
                OperationAllowance::reserve(5_000_000).unwrap(),
                &mut decider(),
            )
            .unwrap()
        });
        let records = batch.records;
        assert!(records.iter().all(|r| r.result.is_ok()));
        match &first {
            None => first = Some(records),
            Some(first) => assert_eq!(first, &records),
        }
    }
}

fn rim() -> Region {
    let mut dec = Decider::default();
    Region::new(
        RollingEnvelopeDefinition::new(
            Point3::new(r("0"), r("5"), r("0")),
            v(0, 1, 0),
            v(1, 0, 0),
            r("12"),
            r("4"),
            r("17"),
            r("2"),
            [r("1"), r("1")],
            &mut dec,
        )
        .unwrap(),
        &mut dec,
    )
    .unwrap()
}

/// Dense binary64 samples of the rim face from the M18 chart formulas, written
/// out independently of the routine: centre (C cos, 5 + H, C sin) and the arc
/// from the cylinder foot to the mouth foot in the normal plane.
fn sampled_rim_distance(point: [f64; 3]) -> f64 {
    let (body, mouth, offset, radius) = (12.0f64, 4.0f64, 17.0f64, 2.0f64);
    let c = body + radius;
    let cos0 = (c * c + offset * offset - mouth * mouth) / (2.0 * c * offset);
    let theta0 = cos0.acos();
    let mut best = f64::INFINITY;
    let (steps, arcs) = (1200, 400);
    for i in 0..=steps {
        let theta = -theta0 + 2.0 * theta0 * i as f64 / steps as f64;
        let (x, z) = (c * theta.cos(), c * theta.sin());
        let d = ((x - offset).powi(2) + z * z).sqrt();
        let h = (radius * radius - (mouth - d).powi(2)).sqrt();
        let centre = [x, 5.0 + h, z];
        let foot = [body * theta.cos(), 5.0 + h, body * theta.sin()];
        let mouth_foot = [offset + mouth / d * (x - offset), 5.0, mouth / d * z];
        let a: Vec<f64> = (0..3).map(|k| (foot[k] - centre[k]) / radius).collect();
        let b: Vec<f64> = (0..3)
            .map(|k| (mouth_foot[k] - centre[k]) / radius)
            .collect();
        let angle = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
            .clamp(-1.0, 1.0)
            .acos();
        for j in 0..=arcs {
            let phi = angle * j as f64 / arcs as f64;
            // Slerp from a to b on the characteristic circle.
            let (sa, sb) = ((angle - phi).sin() / angle.sin(), phi.sin() / angle.sin());
            let y: Vec<f64> = (0..3)
                .map(|k| centre[k] + radius * (sa * a[k] + sb * b[k]))
                .collect();
            let dist =
                ((point[0] - y[0]).powi(2) + (point[1] - y[1]).powi(2) + (point[2] - y[2]).powi(2))
                    .sqrt();
            best = best.min(dist);
        }
    }
    best
}

#[test]
fn rim_distance_brackets_dense_samples() {
    let region = rim();
    let mut dec = decider();
    let width = q("1/1000000");
    for point in [
        ["12.2", "5.2", "0"],
        ["12.5", "5.4", "0"],
        ["12.5", "5.42", "0"],
        ["13", "6", "1"],
        ["12.3", "5.1", "-1.5"],
        ["14", "6", "0"],
    ] {
        let world = Point3::new(r(point[0]), r(point[1]), r(point[2]));
        let range = region.distance(&world, &width, &mut dec).unwrap();
        let sampled = sampled_rim_distance(point.map(|s| s.parse::<f64>().unwrap()));
        let lower = range.lower.to_f64_nearest();
        let upper = range.upper.clone().unwrap().to_f64_nearest();
        // The true distance is at most the sampled minimum and within the
        // sampling resolution of it (chord 0.03 mm, arc 0.006 mm).
        assert!(
            lower <= sampled + 1e-12,
            "{point:?}: {range:?} vs {sampled}"
        );
        assert!(upper >= sampled - 0.05, "{point:?}: {range:?} vs {sampled}");
        assert!(range.converged, "{point:?}: {range:?}");
        println!("rim {point:?}: [{lower}, {upper}] sampled {sampled}");
    }
}

#[test]
fn frame_bounds_are_used_through_the_public_cover() {
    // The rim cover exposes its frames for distance queries only.
    let region = rim();
    let mut dec = decider();
    let cover = path::Path::frames(&region, &mut dec).unwrap();
    let frame: Frame = cover
        .frame(&catallax_numbers::Ball::from_rational(&q("0"), 96))
        .unwrap();
    // At u = 0 the cylinder foot is (12, 5 + sqrt 3, 0) in local coordinates
    // (axis origin at height 5), so local height sqrt 3.
    let h = frame.feet[0][1].clone();
    assert!(h.lo() < q("1.7321") && h.hi() > q("1.7320"));
    assert!(frame.feet[1][1].contains(&Rational::zero()));
}

/// A unit-radius straight centre line `c(t) = (t, 0, 0)` whose feet turn: at
/// station `t` the face arc spans 45 degrees either side of the angle
/// `2 pi t - pi / 2` in the normal plane. It offers the naive nearest hook
/// (the projection onto the line) but is not uniform.
struct Turning;
#[derive(Clone)]
struct TurningCover(catallax_numbers::Ball);
impl Cover for TurningCover {
    fn radius(&self) -> &catallax_numbers::Ball {
        &self.0
    }
    fn centre(&self, u: &catallax_numbers::Ball) -> Result<[catallax_numbers::Ball; 3], Refusal> {
        let zero = catallax_numbers::Ball::from_rational(&Rational::zero(), 96);
        Ok([u.clone(), zero.clone(), zero])
    }
    fn frame(&self, u: &catallax_numbers::Ball) -> Result<Frame, Refusal> {
        use catallax_numbers::Ball;
        let pi = Ball::pi(96);
        let ball = |s: &str| Ball::from_rational(&q(s), 96);
        let middle = &(&(&ball("2") * &pi) * u) - &(&pi * &ball("1/2"));
        let quarter = &pi * &ball("1/4");
        let foot = |angle: Ball| [u.clone(), angle.cos(), angle.sin()];
        Ok(Frame {
            centre: self.centre(u)?,
            tangent: [ball("1"), ball("0"), ball("0")],
            feet: [foot(&middle - &quarter), foot(&middle + &quarter)],
        })
    }
}
impl Path for Turning {
    type Cover = TurningCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        [point.x.clone(), point.y.clone(), point.z.clone()]
    }
    fn walls(&self, _: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        Ok(vec![])
    }
    fn cover(&self, _: &mut Decider) -> Result<TurningCover, Refusal> {
        Ok(TurningCover(catallax_numbers::Ball::from_rational(
            &Rational::one(),
            96,
        )))
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        let x = dec.enclose_budgeted_at(&local[0], 96)?;
        self.cover(dec)?.frame(&x).map(Some)
    }
}

/// Post-merge review of step 1: the arc of the nearest centre bounds the
/// distance to the whole face only on uniform paths. At `(0, 0, 1)` the
/// projection's arc lies within 45 degrees of `-z`, at least
/// `sqrt(2 + sqrt 2) > 1.8` away, while the arc at `t = 3/8` reaches `+z` at
/// `(3/8, 0, 1)`. The true distance, about 0.3703, lies a little nearer still.
#[test]
fn a_nearest_centre_bounds_the_face_only_on_uniform_paths() {
    let point = Point3::new(r("0"), r("0"), r("1"));
    let mut dec = decider();
    let range = path::distance(&Turning, &point, &q("1/1000000"), &mut dec).unwrap();
    assert!(range.converged, "{range:?}");
    assert!(range.upper.as_ref().unwrap() <= &q("3/8"), "{range:?}");
    assert!(range.lower >= q("36/100"), "{range:?}");
    // A dense independent sampling of the arcs agrees.
    let mut best = f64::MAX;
    for i in 0..=200_000 {
        let t = -1.0 + 2.0 * f64::from(i) / 200_000.0;
        let middle = 2.0 * std::f64::consts::PI * t - std::f64::consts::FRAC_PI_2;
        let lo = middle - std::f64::consts::FRAC_PI_4;
        let hi = middle + std::f64::consts::FRAC_PI_4;
        let reach = if lo <= std::f64::consts::FRAC_PI_2 && std::f64::consts::FRAC_PI_2 <= hi {
            1.0
        } else {
            lo.sin().max(hi.sin())
        };
        best = best.min((t * t + 2.0 - 2.0 * reach).sqrt());
    }
    assert!(
        range.lower.to_f64_nearest() <= best + 1e-12,
        "{range:?} {best}"
    );
    assert!(
        range.upper.as_ref().unwrap().to_f64_nearest() >= best - 1e-6,
        "{range:?} {best}"
    );
}
