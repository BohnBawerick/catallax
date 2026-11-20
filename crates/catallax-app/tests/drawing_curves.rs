//! Exact sections and drawing curves: plane sections of solids, including
//! planes through faces, and the binary64 lines, arcs, circles and NURBS
//! that edges, face boundaries, sections and hidden-line strokes become.

use catallax_app::DEFAULT_BUILD_WORK;
use catallax_app::queries::{SectionCurves, SectionPlane, section_curves};
use catallax_app::views::{
    DRAWING_FIT_MARK, DrawingCurve, DrawnCurve, drawing_curve, edge_drawing, face_drawing_loops,
    hidden_line_view, stroke_drawing,
};
use catallax_boolean::{BoolOp, boolean};
use catallax_features::evaluate;
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Point3, Pose, RPoint3, RVector3, Rotation};
use catallax_topology::{BodyId, Model, block, cone, cylinder, sphere};

const TARGET: u64 = 20_000_000;

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn at(x: i64, y: i64, z: i64) -> Pose {
    Pose::from_rational(Rotation::identity(), vector(x, y, z))
}

fn part(source: &str) -> (Model, BodyId) {
    let recipe = Recipe::parse(source).unwrap();
    let (built, _) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let built = built.unwrap();
    (built.model, built.result)
}

fn section(
    model: &Model,
    body: BodyId,
    origin: [i64; 3],
    normal: [i64; 3],
) -> (SectionCurves, u64) {
    let plane = SectionPlane {
        origin: RPoint3::new(origin[0].into(), origin[1].into(), origin[2].into()),
        normal: vector(normal[0], normal[1], normal[2]),
    };
    let (answer, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        section_curves(model, body, &plane, &mut Decider::default())
    });
    let answer = answer.unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    (answer.value, work.used)
}

/// A binary64 image of an exact point.
fn close(p: &Point3) -> [f64; 3] {
    let floor = Rational::from_f64_exact(1e-12).unwrap();
    let mut dec = Decider::default();
    [&p.x, &p.y, &p.z].map(|x| dec.enclose(x, &floor).unwrap().mid().to_f64_nearest())
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// A point of a drawing curve at share `s` of its parameter range, by its
/// own definition: de Boor's algorithm for a NURBS (Piegl and Tiller,
/// algorithm A4.1), independent of the code under test.
fn eval(curve: &DrawingCurve, s: f64) -> [f64; 3] {
    match curve {
        DrawingCurve::Line { a, b } => std::array::from_fn(|k| a[k] + s * (b[k] - a[k])),
        DrawingCurve::Arc {
            center,
            normal,
            x_axis,
            radius,
            start_deg,
            end_deg,
        } => {
            let t = (start_deg + s * (end_deg - start_deg)).to_radians();
            let y = cross(*normal, *x_axis);
            std::array::from_fn(|k| center[k] + radius * (t.cos() * x_axis[k] + t.sin() * y[k]))
        }
        DrawingCurve::Circle {
            center,
            normal,
            radius,
        } => {
            let helper = if normal[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            let x = cross(*normal, helper);
            let x = x.map(|c| c / length(x));
            let y = cross(*normal, x);
            let t = s * std::f64::consts::TAU;
            std::array::from_fn(|k| center[k] + radius * (t.cos() * x[k] + t.sin() * y[k]))
        }
        DrawingCurve::Nurbs {
            degree,
            knots,
            points,
            weights,
        } => {
            let p = *degree;
            let (lo, hi) = (knots[p], knots[knots.len() - p - 1]);
            let u = lo + s * (hi - lo);
            let n = points.len();
            let k = (p..n).rfind(|&k| knots[k] <= u).unwrap_or(p).min(n - 1);
            let w = |i: usize| if weights.is_empty() { 1.0 } else { weights[i] };
            let mut d: Vec<[f64; 4]> = (0..=p)
                .map(|j| {
                    let i = j + k - p;
                    let q = points[i];
                    [w(i) * q[0], w(i) * q[1], w(i) * q[2], w(i)]
                })
                .collect();
            for r in 1..=p {
                for j in (r..=p).rev() {
                    let left = knots[j + k - p];
                    let right = knots[j + 1 + k - r];
                    let alpha = if right > left {
                        (u - left) / (right - left)
                    } else {
                        0.0
                    };
                    for c in 0..4 {
                        d[j][c] = (1.0 - alpha) * d[j - 1][c] + alpha * d[j][c];
                    }
                }
            }
            [d[p][0] / d[p][3], d[p][1] / d[p][3], d[p][2] / d[p][3]]
        }
    }
}

fn first(curve: &DrawingCurve) -> [f64; 3] {
    eval(curve, 0.0)
}

fn last(curve: &DrawingCurve) -> [f64; 3] {
    eval(curve, 1.0)
}

fn curve_length(curve: &DrawingCurve) -> f64 {
    let n = 2000;
    (0..n)
        .map(|i| {
            length(sub(
                eval(curve, (i + 1) as f64 / n as f64),
                eval(curve, i as f64 / n as f64),
            ))
        })
        .sum()
}

/// Every section piece drawn exactly, in section order.
fn drawn(section: &SectionCurves) -> Vec<DrawnCurve> {
    section
        .curves
        .iter()
        .map(|arc| {
            let answer = drawing_curve(
                &section.model,
                arc.curve,
                arc.start,
                arc.end,
                &mut Decider::default(),
            )
            .unwrap();
            assert!(answer.marks.is_empty());
            assert!(answer.value.fit.is_none());
            answer.value
        })
        .collect()
}

#[test]
fn p1_sections_through_the_bore_and_along_the_plate_top_are_exact() {
    let (model, body) = part(include_str!("../../../parts/P1-bracket.json"));
    // Through the boss and bore axis: plate, flange with its slot, boss and
    // bore, all cut across their faces.
    let (across, work) = section(&model, body, [30, 0, 0], [1, 0, 0]);
    eprintln!(
        "P1 section x = 30: {} curves, {work} work",
        across.curves.len()
    );
    assert!(work < TARGET, "{work}");
    assert!(across.points.is_empty());
    assert!(across.curves.iter().all(|arc| arc.edge.is_none()));
    let curves = drawn(&across);
    assert!(
        curves
            .iter()
            .all(|c| matches!(c.curve, DrawingCurve::Line { .. }))
    );
    // The pieces close up: every end meets another piece's end.
    let ends: Vec<[f64; 3]> = curves
        .iter()
        .flat_map(|c| [first(&c.curve), last(&c.curve)])
        .collect();
    for (i, e) in ends.iter().enumerate() {
        let partners = ends
            .iter()
            .enumerate()
            .filter(|(j, f)| *j != i && length(sub(*e, **f)) < 1e-9)
            .count();
        assert_eq!(partners % 2, 1, "{e:?}");
    }
    // Through the plate's top face, which holds the bolt holes' rims and
    // meets the boss and the flange: every edge in the plane once, and the
    // bore's wall cut across.
    let (top, work) = section(&model, body, [0, 0, 5], [0, 0, 1]);
    eprintln!("P1 section z = 5: {work} work");
    assert!(work < TARGET, "{work}");
    let edges = top.curves.iter().filter(|arc| arc.edge.is_some()).count();
    let mut seen = std::collections::BTreeSet::new();
    assert!(
        top.curves
            .iter()
            .filter_map(|arc| arc.edge)
            .all(|e| seen.insert(e))
    );
    assert!(edges > 0);
    let circles: Vec<f64> = drawn(&top)
        .iter()
        .filter_map(|c| match c.curve {
            DrawingCurve::Circle { radius, .. } => Some(radius),
            _ => None,
        })
        .collect();
    // Four bolt holes and the boss stand on the plane; the bore crosses it.
    let count = |r: f64| circles.iter().filter(|&&x| (x - r).abs() < 1e-12).count();
    assert_eq!(
        (count(2.0), count(8.0), count(4.0)),
        (4, 1, 1),
        "{circles:?}"
    );
}

#[test]
fn a_step_in_the_plane_gives_its_edges_and_cuts_once() {
    let mut model = Model::new();
    let base = block(
        &mut model,
        [4.into(), 4.into(), 2.into()],
        &Pose::identity(),
    )
    .unwrap();
    let upper = block(&mut model, [2.into(), 2.into(), 2.into()], &at(0, 0, 2)).unwrap();
    let body = boolean(
        &mut model,
        base,
        upper,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let (step, _) = section(&model, body, [0, 0, 2], [0, 0, 1]);
    assert!(step.points.is_empty());
    // The outer square and the two inner edges of the step: 16 plus 4.
    let total: f64 = drawn(&step).iter().map(|c| curve_length(&c.curve)).sum();
    assert!((total - 20.0).abs() < 1e-9, "{total}");
}

#[test]
fn conic_sections_draw_as_exact_conics() {
    let mut model = Model::new();
    let frustum = cone(&mut model, 4.into(), 1.into(), 6.into(), &Pose::identity()).unwrap();
    // `r(z) = 4 - z / 2` on the side.
    let on_cone = |p: [f64; 3]| ((p[0] * p[0] + p[1] * p[1]).sqrt() - (4.0 - p[2] / 2.0)).abs();
    for (normal, offset, kind) in [
        ([1, 0, 0], [2, 0, 0], "hyperbola"),
        ([2, 0, 1], [3, 0, 0], "parabola"),
    ] {
        let (cut, _) = section(&model, frustum, offset, normal);
        let mut found = 0;
        for (arc, curve) in cut.curves.iter().zip(drawn(&cut)) {
            let form = curve_form(cut.model.geom(), arc.curve, &mut Decider::default()).unwrap();
            let conic = match form {
                CurveForm::Hyperbola { .. } => "hyperbola",
                CurveForm::Parabola { .. } => "parabola",
                _ => continue,
            };
            assert_eq!(conic, kind);
            found += 1;
            let DrawingCurve::Nurbs {
                degree, weights, ..
            } = &curve.curve
            else {
                panic!("{conic}: {:?}", curve.curve);
            };
            assert_eq!(*degree, 2);
            assert_eq!(weights.is_empty(), kind == "parabola", "{weights:?}");
            for i in 0..=100 {
                let p = eval(&curve.curve, i as f64 / 100.0);
                assert!(on_cone(p) < 1e-12, "{kind} {p:?}");
                let n = normal.map(|c| c as f64);
                let o = offset.map(|c| c as f64);
                let off = n.iter().zip(sub(p, o)).map(|(a, b)| a * b).sum::<f64>();
                assert!(off.abs() < 1e-12, "{kind} {p:?}");
            }
        }
        assert!(found > 0, "{kind}");
    }
    // An oblique plane through a cylinder: a whole ellipse.
    let tube = cylinder(&mut model, 2.into(), 10.into(), &Pose::identity()).unwrap();
    let (oblique, _) = section(&model, tube, [0, 0, 5], [1, 0, 2]);
    let ellipses: Vec<DrawnCurve> = drawn(&oblique)
        .into_iter()
        .filter(|c| matches!(c.curve, DrawingCurve::Nurbs { .. }))
        .collect();
    assert_eq!(ellipses.len(), 1);
    let curve = &ellipses[0].curve;
    assert!(length(sub(first(curve), last(curve))) < 1e-12);
    for i in 0..=200 {
        let p = eval(curve, i as f64 / 200.0);
        assert!(
            ((p[0] * p[0] + p[1] * p[1]).sqrt() - 2.0).abs() < 1e-12,
            "{p:?}"
        );
        assert!((p[0] + 2.0 * p[2] - 10.0).abs() < 1e-12, "{p:?}");
    }
}

#[test]
fn face_loops_chain_and_strokes_keep_their_exact_ends_on_the_p1_bracket() {
    let (model, body) = part(include_str!("../../../parts/P1-bracket.json"));
    let mut arcs = 0;
    for face in model.faces_of_body(body) {
        let loops = face_drawing_loops(&model, body, face, &mut Decider::default()).unwrap();
        assert!(loops.marks.is_empty());
        for chain in &loops.value {
            for (i, e) in chain.iter().enumerate() {
                assert!(e.drawn.fit.is_none());
                if matches!(e.drawn.curve, DrawingCurve::Arc { .. }) {
                    arcs += 1;
                }
                let next = &chain[(i + 1) % chain.len()];
                let gap = length(sub(last(&e.drawn.curve), first(&next.drawn.curve)));
                assert!(gap < 1e-9, "{face:?} {i} {gap}");
            }
        }
    }
    assert!(arcs > 0);
    // Each edge drawn from its start vertex to its end vertex.
    for edge in model.edges_of_body(body) {
        let drawn = edge_drawing(&model, body, edge, &mut Decider::default()).unwrap();
        let stored = model.edge(edge);
        if let (Some(a), Some(b)) = (stored.start, stored.end) {
            let a = close(&model.geom().position(model.vertex(a).point));
            let b = close(&model.geom().position(model.vertex(b).point));
            assert!(length(sub(first(&drawn.value.curve), a)) < 1e-9);
            assert!(length(sub(last(&drawn.value.curve), b)) < 1e-9);
        }
    }
    let view = hidden_line_view(&model, body, &vector(1, -1, 1), &mut Decider::default()).unwrap();
    for stroke in &view.value.strokes {
        let drawn = stroke_drawing(&view.value, stroke, &mut Decider::default()).unwrap();
        if let (Some(a), Some(b)) = (stroke.start, stroke.end) {
            let a = close(&view.value.model.geom().position(a));
            let b = close(&view.value.model.geom().position(b));
            assert!(length(sub(first(&drawn.value.curve), a)) < 1e-9);
            assert!(length(sub(last(&drawn.value.curve), b)) < 1e-9);
        }
    }
}

#[test]
fn intersection_branches_are_fitted_and_marked() {
    let mut model = Model::new();
    let ball = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let hole = cylinder(&mut model, 2.into(), 20.into(), &at(1, 0, -10)).unwrap();
    let body = boolean(
        &mut model,
        ball,
        hole,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let mut branches = 0;
    for edge in model.edges_of_body(body) {
        let form = curve_form(
            model.geom(),
            model.edge(edge).curve,
            &mut Decider::default(),
        )
        .unwrap();
        if !matches!(form, CurveForm::Branch) {
            continue;
        }
        branches += 1;
        let drawn = edge_drawing(&model, body, edge, &mut Decider::default()).unwrap();
        assert!(drawn.marks.iter().any(|m| m.code == DRAWING_FIT_MARK));
        let fit = drawn.value.fit.clone().expect("a fitted branch");
        assert!(fit.deviation_mm <= fit.target_mm, "{fit:?}");
        // Independently: the fitted curve stays within its tolerance of
        // both the sphere and the cylinder.
        for i in 0..=1000 {
            let p = eval(&drawn.value.curve, i as f64 / 1000.0);
            let ball = (length(p) - 5.0).abs();
            let wall = (((p[0] - 1.0).powi(2) + p[1] * p[1]).sqrt() - 2.0).abs();
            assert!(
                ball <= fit.tolerance_mm && wall <= fit.tolerance_mm,
                "{p:?} {fit:?}"
            );
        }
    }
    assert_eq!(branches, 2);
}
