//! 4.2a step 7.2: the patch-shell leaf on the M27 loft.
//!
//! The loft is that of `catallax-features/tests/m11_spline_fillet.rs`: ruled,
//! from the square of half-width 1 at height 0 to the quadrilateral
//! `(1,0), (0,2), (-2,0), (0,-1)` at height 2. The field's shell is built
//! here from those numbers alone, four bilinear sides and two polygon caps,
//! and every answer is checked against the B-rep loft.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{Built, evaluate};
use catallax_field::Location;
use catallax_field::patch_shell::{Face, PatchFace, Polygon, Shell};
use catallax_numbers::{BudgetPolicy, Decider, Interval, Outcome, Rational, Real};
use catallax_recipe::Recipe;
use catallax_shapes::{BSplineSurface, Isometry, Point3, RPoint3, RVector3, RationalPatch};
use serde_json::{Value, json};

type Q = Rational;

fn q(n: i64) -> Q {
    Q::from_i64(n)
}

fn r(text: &str) -> Q {
    Q::parse(text).unwrap()
}

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}

const BOTTOM: [(i64, i64); 4] = [(-1, -1), (1, -1), (1, 1), (-1, 1)];
const TOP: [(i64, i64); 4] = [(0, -1), (1, 0), (0, 2), (-2, 0)];

fn polygon(prefix: &str, points: &[(&str, &str)]) -> Value {
    json!(
        points
            .iter()
            .enumerate()
            .map(|(i, (x, y))| json!({"kind":"line","id":format!("{prefix}_{i}"),"to":[x,y]}))
            .collect::<Vec<_>>()
    )
}

/// The M11 loft's recipe: alone, with the second section starting at
/// `top_start`, or with the M11 fillet of radius 1/10 after it.
fn recipe(top_start: &str, fillet: bool) -> Recipe {
    let bottom = polygon("a", &[("1", "-1"), ("1", "1"), ("-1", "1"), ("-1", "-1")]);
    let top = polygon("b", &[("1", "0"), ("0", "2"), ("-2", "0"), ("0", "-1")]);
    let section = |i: usize, height: &str, loops: Value| {
        json!({"id":format!("section{i}"),"op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0",height]},
            "loops":[loops]})
    };
    let side = |segment: &str| json!({"step":"loft","role":{"side":{"segment":segment}}});
    let mut steps = vec![
        section(0, "0", bottom),
        section(1, "2", top),
        json!({"id":"loft","op":"loft","shape":"ruled","sections":[
            {"profile":"section0","start":["a_0"]},
            {"profile":"section1","start":[top_start]}]}),
    ];
    if fillet {
        steps.push(
            json!({"id":"round","op":"fillet","body":"loft","radius":"1/10",
            "edges":[{"between":[side("a_0"), side("a_1")]}]}),
        );
    }
    let value = json!({
        "format":"catallax-recipe/5","name":"m11 loft",
        "steps":steps,
        "result": if fillet { "round" } else { "loft" }
    });
    Recipe::parse(&value.to_string()).unwrap()
}

/// The M11 loft alone, built by the B-rep.
fn brep() -> Built {
    let bottom = polygon("a", &[("1", "-1"), ("1", "1"), ("-1", "1"), ("-1", "-1")]);
    let top = polygon("b", &[("1", "0"), ("0", "2"), ("-2", "0"), ("0", "-1")]);
    let section = |i: usize, height: &str, loops: Value| {
        json!({"id":format!("section{i}"),"op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0",height]},
            "loops":[loops]})
    };
    let value = json!({
        "format":"catallax-recipe/5","name":"m11 loft",
        "steps":[
            section(0, "0", bottom),
            section(1, "2", top),
            {"id":"loft","op":"loft","shape":"ruled","sections":[
                {"profile":"section0","start":["a_0"]},
                {"profile":"section1","start":["b_0"]}]}
        ],
        "result":"loft"
    });
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    evaluate(&recipe, &mut decider()).unwrap()
}

fn at(corner: (i64, i64), height: i64) -> Point3 {
    Point3::new(
        Real::from_i64(corner.0),
        Real::from_i64(corner.1),
        Real::from_i64(height),
    )
}

/// The bilinear side between the `k`-th bottom and top segments: `u` along
/// the segments, `v` from bottom to top, so its natural normal points out of
/// the counterclockwise sections.
fn side(k: usize, dec: &mut Decider) -> RationalPatch {
    let next = (k + 1) % 4;
    let controls = vec![
        at(BOTTOM[k], 0),
        at(TOP[k], 2),
        at(BOTTOM[next], 0),
        at(TOP[next], 2),
    ];
    let knots = vec![q(0), q(0), q(1), q(1)];
    let carrier = BSplineSurface::new(
        [1, 1],
        [knots.clone(), knots],
        [2, 2],
        controls,
        vec![Real::one(); 4],
        dec,
    )
    .unwrap();
    RationalPatch::new(carrier, dec).unwrap()
}

fn cap(corners: &[(i64, i64)], height: i64) -> Polygon {
    Polygon::new(
        corners
            .iter()
            .map(|&(x, y)| RPoint3::new(q(x), q(y), q(height)))
            .collect(),
    )
    .unwrap()
}

/// The field's shell of the loft.
fn shell() -> Shell {
    let mut dec = decider();
    let mut faces: Vec<Face> = (0..4)
        .map(|k| {
            Face::Patch(Box::new(
                PatchFace::new(&side(k, &mut dec), true, &mut dec).unwrap(),
            ))
        })
        .collect();
    let mut bottom = BOTTOM.to_vec();
    bottom.reverse();
    faces.push(Face::Polygon(Box::new(cap(&bottom, 0))));
    faces.push(Face::Polygon(Box::new(cap(&TOP, 2))));
    Shell::new(faces).unwrap()
}

/// The B-rep's answer, or `None` where it refuses (it does at a few grid
/// points by a loft corner, where its ray is tangent to a side).
fn brep_location(built: &mut Built, p: &RPoint3) -> Option<Location> {
    match classify_point(&mut built.model, built.result, p, &mut decider()) {
        Outcome::Proven(Containment::Inside) => Some(Location::Inside),
        Outcome::Proven(Containment::Outside) => Some(Location::Outside),
        Outcome::ProvenCoincident(_) => Some(Location::Boundary),
        Outcome::Refused(refusal) if refusal.code == "shapes.rational-patch-unproven" => None,
        other => panic!("the B-rep does not settle {p:?}: {other:?}"),
    }
}

/// A grid of points round the loft and past it, 0.3 apart across and
/// 0.325 apart in height, none on a cap's plane.
fn grid() -> Vec<RPoint3> {
    let mut out = Vec::new();
    for i in 0..13 {
        for j in 0..13 {
            for k in 0..9 {
                let x = &r("-23/10") + &(&q(i) * &r("3/10"));
                let y = &r("-13/10") + &(&q(j) * &r("3/10"));
                let z = &r("-3/10") + &(&q(k) * &r("13/40"));
                out.push(RPoint3::new(x, y, z));
            }
        }
    }
    out
}

/// The point at parameters `(u, v)` of side `k`, from the corners.
fn on_side(k: usize, u: &Q, v: &Q) -> RPoint3 {
    let next = (k + 1) % 4;
    let lerp = |a: i64, b: i64, t: &Q| &q(a) + &(t * &q(b - a));
    let bottom = [
        lerp(BOTTOM[k].0, BOTTOM[next].0, u),
        lerp(BOTTOM[k].1, BOTTOM[next].1, u),
    ];
    let top = [
        lerp(TOP[k].0, TOP[next].0, u),
        lerp(TOP[k].1, TOP[next].1, u),
    ];
    let mix = |a: &Q, b: &Q| a + &(v * &(b - a));
    RPoint3::new(
        mix(&bottom[0], &top[0]),
        mix(&bottom[1], &top[1]),
        v * &q(2),
    )
}

/// The third check: the smallest binary64 squared distance from `p` to the
/// four sides (a 32 by 32 grid, then zoomed 5 by 5 grids) and the two caps
/// (as polygons in the plane, or their edges).
fn sampled(p: [f64; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for k in 0..4 {
        let next = (k + 1) % 4;
        let point = |u: f64, v: f64| -> [f64; 3] {
            let b = [
                BOTTOM[k].0 as f64 + u * (BOTTOM[next].0 - BOTTOM[k].0) as f64,
                BOTTOM[k].1 as f64 + u * (BOTTOM[next].1 - BOTTOM[k].1) as f64,
            ];
            let t = [
                TOP[k].0 as f64 + u * (TOP[next].0 - TOP[k].0) as f64,
                TOP[k].1 as f64 + u * (TOP[next].1 - TOP[k].1) as f64,
            ];
            [b[0] + v * (t[0] - b[0]), b[1] + v * (t[1] - b[1]), 2.0 * v]
        };
        let gap = |u: f64, v: f64| {
            let s = point(u, v);
            (0..3).map(|i| (s[i] - p[i]).powi(2)).sum::<f64>()
        };
        let mut here = (f64::INFINITY, 0.0, 0.0);
        for i in 0..=32 {
            for j in 0..=32 {
                let (u, v) = (f64::from(i) / 32.0, f64::from(j) / 32.0);
                if gap(u, v) < here.0 {
                    here = (gap(u, v), u, v);
                }
            }
        }
        let mut step = 1.0 / 64.0;
        for _ in 0..40 {
            let (_, cu, cv) = here;
            for i in -2..=2 {
                for j in -2..=2 {
                    let u = (cu + f64::from(i) * step).clamp(0.0, 1.0);
                    let v = (cv + f64::from(j) * step).clamp(0.0, 1.0);
                    if gap(u, v) < here.0 {
                        here = (gap(u, v), u, v);
                    }
                }
            }
            step /= 2.0;
        }
        best = best.min(here.0);
    }
    for (corners, height) in [(&BOTTOM, 0.0), (&TOP, 2.0)] {
        let inside = (0..4).all(|i| {
            let (a, b) = (corners[i], corners[(i + 1) % 4]);
            let (ax, ay, bx, by) = (a.0 as f64, a.1 as f64, b.0 as f64, b.1 as f64);
            (bx - ax) * (p[1] - ay) - (by - ay) * (p[0] - ax) >= 0.0
        });
        if inside {
            best = best.min((p[2] - height).powi(2));
        }
    }
    best
}

#[test]
fn the_shell_classifies_as_the_brep_loft() {
    let mut built = brep();
    let shell = shell();
    let mut counts = [0usize; 4];
    let mut refused = Vec::new();
    for p in grid() {
        let (location, distance) = shell.classify(&p, &mut decider()).unwrap();
        let Some(expected) = brep_location(&mut built, &p) else {
            refused.push((p.clone(), location));
            continue;
        };
        if location == Location::Indeterminate {
            // Only a point on a side may stay open.
            assert_eq!(distance.lo(), 0.0, "{p:?}");
            assert_eq!(expected, Location::Boundary, "{p:?}");
        } else {
            assert_eq!(location, expected, "{p:?}");
        }
        counts[match location {
            Location::Inside => 0,
            Location::Outside => 1,
            Location::Boundary => 2,
            _ => 3,
        }] += 1;
        // The third check on the distance.
        let fast = [&p.x, &p.y, &p.z].map(|x| x.to_f64_nearest());
        let found = sampled(fast);
        let squared = distance.square();
        assert!(
            found >= squared.lo() - 1e-12 && found <= squared.hi() + 1e-12,
            "{p:?}: sampled {found}, shell {distance:?}"
        );
        assert!(distance.hi() - distance.lo() <= catallax_field::patch_distance::WIDTH);
    }
    println!("grid: {counts:?} (inside, outside, boundary, indeterminate)");
    println!("the B-rep refuses {} points: {refused:?}", refused.len());
    assert!(refused.len() <= 10, "{refused:?}");
    assert!(counts[0] > 100 && counts[1] > 100, "{counts:?}");
}

/// Points of the caps are on the boundary exactly; points of the sides at
/// rational parameters have the distance bounds `[0, w]` and stay
/// `Indeterminate`, never a wrong side. Points a millimetre off a side's
/// centre along its normal are inside and outside as in the B-rep.
#[test]
fn boundary_points_are_settled_or_left_open() {
    let mut built = brep();
    let shell = shell();
    for p in [
        RPoint3::new(r("1/3"), r("1/5"), q(0)),
        RPoint3::new(r("-1/2"), r("1/2"), q(2)),
    ] {
        assert_eq!(
            shell.classify(&p, &mut decider()).unwrap().0,
            Location::Boundary
        );
        assert_eq!(brep_location(&mut built, &p), Some(Location::Boundary));
    }
    for k in 0..4 {
        for (u, v) in [("1/2", "1/2"), ("1/7", "5/6"), ("3/4", "1/9")] {
            let p = on_side(k, &r(u), &r(v));
            let (location, distance) = shell.classify(&p, &mut decider()).unwrap();
            assert_eq!(location, Location::Indeterminate, "side {k} at {u}, {v}");
            assert_eq!(distance.lo(), 0.0);
            assert_eq!(brep_location(&mut built, &p), Some(Location::Boundary));
        }
    }
}

/// A rotated and reflected copy answers as the original at the moved point.
#[test]
fn a_moved_shell_answers_as_the_original() {
    let shell = shell();
    let mut moved = shell.clone();
    let mut dec = decider();
    let turn = Isometry::around_axis(
        &RPoint3::new(q(1), q(2), q(3)),
        &RVector3::new(q(0), q(0), q(1)),
        &q(90),
        &mut dec,
    )
    .unwrap();
    let mirror = Isometry::reflection(
        &RPoint3::new(q(5), q(0), q(0)),
        &RVector3::new(q(1), q(0), q(0)),
    )
    .unwrap();
    let map = turn.then(&mirror);
    moved.transform(&map).unwrap();
    for p in grid().into_iter().step_by(7) {
        let image = map.apply_rational_point(&p).unwrap();
        let a = shell.classify(&p, &mut decider()).unwrap();
        let b = moved.classify(&image, &mut decider()).unwrap();
        assert_eq!(a, b, "{p:?}");
    }
}

/// Answers and work are the same on 1, 2 and 4 threads.
#[test]
fn shell_answers_are_the_same_at_every_worker_count() {
    let shell = shell();
    let points: Vec<RPoint3> = grid().into_iter().step_by(5).collect();
    let run = |workers: usize| -> Vec<(Location, Interval, u64)> {
        let mut out = vec![None; points.len()];
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers)
                .map(|worker| {
                    let (shell, points) = (&shell, &points);
                    scope.spawn(move || {
                        let mut mine = Vec::new();
                        for (i, p) in points.iter().enumerate().skip(worker).step_by(workers) {
                            let mut dec = decider();
                            let (l, d) = shell.classify(p, &mut dec).unwrap();
                            mine.push((i, (l, d, dec.stats().work)));
                        }
                        mine
                    })
                })
                .collect();
            for handle in handles {
                for (i, answer) in handle.join().unwrap() {
                    out[i] = Some(answer);
                }
            }
        });
        out.into_iter().map(Option::unwrap).collect()
    };
    let one = run(1);
    assert_eq!(one, run(2));
    assert_eq!(one, run(4));
}

/// The field compiles the M11 loft through its ruled-loft route as a patch
/// shell, and its plan answers as the shell and the B-rep do.
#[test]
fn the_field_builds_the_m11_loft_as_a_patch_shell() {
    use catallax_field::{Plan, Value, Workspace};
    let plan = Plan::compile(&recipe("b_0", false)).unwrap();
    let shell = shell();
    let mut built = brep();
    let mut agreed = 0;
    for p in grid().into_iter().step_by(3) {
        let query = plan
            .proven_checked(&p, &mut Workspace::default(), &mut decider())
            .unwrap();
        let (location, distance) = shell.classify(&p, &mut decider()).unwrap();
        assert_eq!(query.location, location, "{p:?}");
        let Value::TrueDistance(value) = query.value else {
            panic!("{p:?}: {:?}", query.value);
        };
        let signed = if location == Location::Inside {
            -distance.hi()..=-distance.lo()
        } else {
            distance.lo()..=distance.hi()
        };
        assert!(
            value.lo() <= *signed.start() && *signed.end() <= value.hi(),
            "{p:?}: {value:?} against {distance:?}"
        );
        if let Some(expected) = brep_location(&mut built, &p)
            && location != Location::Indeterminate
        {
            assert_eq!(location, expected, "{p:?}");
            agreed += 1;
        }
    }
    println!("plan agrees with the B-rep at {agreed} points");
    assert!(agreed > 400);
}

/// A loft whose correspondence turns the top half a turn collapses its
/// middle slice, so no slice proof holds and the old refusal stands.
#[test]
fn a_loft_whose_slices_are_not_convex_still_refuses() {
    use catallax_field::Plan;
    let refusal = Plan::compile(&recipe("b_2", false)).unwrap_err();
    assert_eq!(refusal.code, "field.loft-correspondence-unsupported");
}

#[test]
fn the_m11_fillet_replaces_its_named_loft_edge() {
    use catallax_field::{Plan, Workspace};
    use catallax_numbers::with_operation_budget;
    use catallax_recipe::Role;
    let (plan, work) = with_operation_budget(20_000_000, || Plan::compile(&recipe("b_0", true)));
    let plan = plan.unwrap();
    let origin = plan.spline_blend_origin().expect("rounded loft face");
    assert!(matches!(origin.role, Role::Blend { .. }));
    assert_eq!(origin.step.as_str(), "round");
    println!("M27 rounded loft: {} work", work.used);
    for ([x, y, z], location) in [
        (["1", "-1/2", "1"], Location::Outside),
        (["1", "0", "2"], Location::Outside),
        (["0", "0", "0"], Location::Boundary),
        (["0", "0", "1"], Location::Inside),
    ] {
        let p = RPoint3::new(r(x), r(y), r(z));
        let result = plan
            .proven_checked(&p, &mut Workspace::default(), &mut decider())
            .unwrap();
        assert_eq!(result.location, location, "{p:?}: {result:?}");
    }
}

/// The print classifier settles a box of the loft only where the box meets
/// no face's box, and then as the shell does at its corners and centre; a
/// box across a side stays unknown.
#[test]
fn print_boxes_of_the_loft_agree_with_the_shell() {
    use catallax_field::Plan;
    use catallax_field::printing::{BoxState, classify_box};
    let plan = Plan::compile(&recipe("b_0", false)).unwrap();
    let shell = shell();
    let mut settled = [0usize; 3];
    let size = r("1/4");
    for i in 0..13 {
        for j in 0..13 {
            for k in 0..9 {
                let lo = [
                    &r("-9/4") + &(&q(i) * &size),
                    &r("-5/4") + &(&q(j) * &size),
                    &r("-1/4") + &(&q(k) * &size),
                ];
                let hi: [Q; 3] = std::array::from_fn(|a| &lo[a] + &size);
                let state = classify_box(&plan, &lo, &hi, &mut decider()).unwrap();
                let expected = match state {
                    BoxState::Inside => Location::Inside,
                    BoxState::Outside => Location::Outside,
                    _ => {
                        settled[2] += 1;
                        continue;
                    }
                };
                settled[usize::from(expected == Location::Outside)] += 1;
                let half: [Q; 3] = std::array::from_fn(|a| (&lo[a] + &hi[a]).mul_pow2(-1));
                for corner in 0..9 {
                    let p: [Q; 3] = std::array::from_fn(|a| match corner {
                        8 => half[a].clone(),
                        _ if corner >> a & 1 == 1 => hi[a].clone(),
                        _ => lo[a].clone(),
                    });
                    let p = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
                    assert_eq!(
                        shell.classify(&p, &mut decider()).unwrap().0,
                        expected,
                        "box {lo:?} at {p:?}"
                    );
                }
            }
        }
    }
    println!("print boxes: {settled:?} (inside, outside, unknown)");
    assert!(settled[0] > 0 && settled[1] > 0);
    let centre = on_side(0, &r("1/2"), &r("1/2"));
    let c = [centre.x, centre.y, centre.z];
    let lo: [Q; 3] = std::array::from_fn(|a| &c[a] - &r("1/100"));
    let hi: [Q; 3] = std::array::from_fn(|a| &c[a] + &r("1/100"));
    assert!(!matches!(
        classify_box(&plan, &lo, &hi, &mut decider()).unwrap(),
        BoxState::Inside | BoxState::Outside
    ));
}
