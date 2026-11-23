//! Hidden-line views of the P1 bracket and the P3 filleted boss, of
//! cylinders, cones, spheres and tori, of spline faces, and of edges on
//! intersection branches.

use std::num::NonZeroUsize;

use catallax_app::DEFAULT_BUILD_WORK;
use catallax_app::views::{
    CONTOUR_FIT_MARK, DRAWING_FIT_MARK, DrawingCurve, HiddenLineView, StrokeSource, VIEW_RAY_MARK,
    VIEW_SHORT_PIECE_MARK, VIEW_SPLIT_MARK, ViewStroke, contour_generators, hidden_line_view,
    stroke_drawing,
};
use catallax_boolean::{BoolOp, boolean};
use catallax_features::evaluate;
use catallax_numbers::{Answer, Decider, Rational, Real, set_threads, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, Point3, PointId, Pose, RVector3, Rotation, Surface};
use catallax_topology::{BodyId, EdgeId, FaceId, Model, block, cone, cylinder, sphere, torus};

#[allow(dead_code)]
mod rational_patch_cases;
mod view_records;
use rational_patch_cases::{Roof, solid};

const TARGET: u64 = 20_000_000;

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn part(source: &str) -> (Model, BodyId) {
    let recipe = Recipe::parse(source).unwrap();
    let (built, _) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let built = built.unwrap();
    (built.model, built.result)
}

/// A binary64 image of an exact point, good to the last few bits.
fn close(p: &Point3) -> [f64; 3] {
    let floor = Rational::from_f64_exact(1e-12).unwrap();
    let mut dec = Decider::default();
    [&p.x, &p.y, &p.z].map(|x| dec.enclose(x, &floor).unwrap().mid().to_f64_nearest())
}

fn has(answer: &Answer<HiddenLineView>, code: &str) -> bool {
    answer.marks.iter().any(|m| m.code == code)
}

/// What a test compares across thread counts: each stroke's source,
/// visibility, smoothness and ends.
type Record = Vec<(StrokeSource, bool, bool, Option<[f64; 3]>, Option<[f64; 3]>)>;

fn record(view: &HiddenLineView) -> Record {
    view.strokes
        .iter()
        .map(|s| {
            let end = |p: Option<PointId>| p.map(|p| close(&view.model.geom().position(p)));
            (s.source, s.visible, s.smooth, end(s.start), end(s.end))
        })
        .collect()
}

/// The view at 1, 2 and 4 threads: the same strokes and work each time,
/// inside the target.
fn steady(model: &Model, body: BodyId, d: &RVector3) -> Answer<HiddenLineView> {
    let mut first: Option<(Answer<HiddenLineView>, Record, u64)> = None;
    for threads in [1, 2, 4] {
        set_threads(NonZeroUsize::new(threads));
        let (answer, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
            hidden_line_view(model, body, d, &mut Decider::default())
        });
        assert!(!work.exhausted);
        let answer = answer.unwrap();
        assert!(work.used < TARGET, "{} work units", work.used);
        let r = record(&answer.value);
        match &first {
            None => first = Some((answer, r, work.used)),
            Some((_, r0, w0)) => {
                assert_eq!(*w0, work.used, "work at {threads} threads");
                assert_eq!(*r0, r, "strokes at {threads} threads");
            }
        }
    }
    set_threads(None);
    let (answer, _, work) = first.unwrap();
    eprintln!("{d:?}: {:?}, {work} work", counts(&answer.value));
    view_records::emit(&answer, d);
    answer
}

/// (visible, hidden, contour strokes, smooth strokes)
fn counts(view: &HiddenLineView) -> (usize, usize, usize, usize) {
    let visible = view.strokes.iter().filter(|s| s.visible).count();
    (
        visible,
        view.strokes.len() - visible,
        view.strokes
            .iter()
            .filter(|s| matches!(s.source, StrokeSource::Contour(_)))
            .count(),
        view.strokes.iter().filter(|s| s.smooth).count(),
    )
}

fn radius_sq(model: &Model, face: FaceId) -> Option<f64> {
    match model.geom().surface(model.face(face).surface) {
        Surface::Cylinder(c) => c.radius_sq.as_rational().map(|r| r.to_f64_nearest()),
        _ => None,
    }
}

/// The edge between two faces.
fn edge_between(model: &Model, body: BodyId, f: FaceId, g: FaceId) -> EdgeId {
    model
        .edges_of_body(body)
        .into_iter()
        .find(|&e| {
            let faces: Vec<FaceId> = model
                .coedges_of_edge(e)
                .into_iter()
                .filter_map(|c| model.owner_of_coedge(c).map(|o| o.1))
                .collect();
            faces.contains(&f) && faces.contains(&g)
        })
        .expect("adjacent faces share an edge")
}

#[test]
fn p1_bracket_views_are_exact_inside_the_target() {
    let (model, body) = part(include_str!("../../../parts/P1-bracket.json"));
    // (visible, hidden, contour strokes, smooth strokes)
    for (name, d, expected) in [
        ("front", vector(0, -1, 0), (26, 21, 12, 0)),
        ("top", vector(0, 0, 1), (32, 6, 2, 4)),
        ("right", vector(1, 0, 0), (31, 17, 12, 0)),
    ] {
        let answer = steady(&model, body, &d);
        assert_eq!(counts(&answer.value), expected, "{name}");
        assert!(
            !has(&answer, VIEW_SPLIT_MARK) && !has(&answer, VIEW_RAY_MARK),
            "{name}"
        );
        // The boss stands behind the flange from the front and in the open
        // from the right; the bore is inside the boss either way.
        for s in &answer.value.strokes {
            if let StrokeSource::Contour(face) = s.source {
                match radius_sq(&model, face) {
                    Some(64.0) => assert_eq!(s.visible, name == "right", "boss {name}"),
                    Some(16.0) => assert!(!s.visible, "bore {name}"),
                    _ => {}
                }
            }
        }
    }
    // The final-output control: in the right view four generators of the
    // slot ends run along their tangent edges with the slot walls. Each is
    // drawn as that edge, under the edge's id, as outline rather than a
    // smooth edge, hidden inside the flange; from the top the same edges are
    // smooth.
    let right = vector(1, 0, 0);
    let generators = contour_generators(&model, body, &right, &mut Decider::default())
        .unwrap()
        .value;
    let along: Vec<EdgeId> = generators
        .arcs
        .iter()
        .filter(|a| a.on_faces.len() > 1)
        .map(|a| edge_between(&model, body, a.on_faces[0], a.on_faces[1]))
        .collect();
    assert_eq!(along.len(), 4);
    let side = steady(&model, body, &right).value;
    let top = steady(&model, body, &vector(0, 0, 1)).value;
    for edge in &along {
        let pieces: Vec<&ViewStroke> = side
            .strokes
            .iter()
            .filter(|s| s.source == StrokeSource::Edge(*edge))
            .collect();
        assert!(!pieces.is_empty());
        assert!(pieces.iter().all(|s| !s.smooth && !s.visible), "{edge:?}");
        let from_top: Vec<&ViewStroke> = top
            .strokes
            .iter()
            .filter(|s| s.source == StrokeSource::Edge(*edge))
            .collect();
        assert!(
            !from_top.is_empty() && from_top.iter().all(|s| s.smooth),
            "{edge:?}"
        );
    }
    // Twelve drawn generators: sixteen less the four drawn as edges.
    assert_eq!(generators.arcs.len(), 16);
    let iso = steady(&model, body, &vector(1, -1, 1));
    assert_eq!(counts(&iso.value), (28, 36, 14, 5));
    // Circles in perpendicular planes cross on the iso sheet: binary64 cuts.
    assert!(has(&iso, VIEW_SPLIT_MARK));
}

#[test]
fn a_cylinder_seen_from_above_hides_the_back_of_its_bottom_circle_exactly() {
    let mut model = Model::new();
    let body = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let answer = steady(&model, body, &vector(1, 0, 1));
    assert!(!has(&answer, VIEW_SPLIT_MARK) && !has(&answer, VIEW_RAY_MARK));
    let view = &answer.value;
    // Top circle whole and visible, two visible generators x = 0, and the
    // bottom circle cut exactly where they end: the front half visible, the
    // back half (through x = -2) hidden.
    assert_eq!(counts(view), (4, 1, 2, 0));
    let geom = view.model.geom();
    let zero = |x: &Real| x.as_rational().is_some_and(|q| q.is_zero());
    for s in &view.strokes {
        let ends = [s.start, s.end].map(|p| p.map(|p| geom.position(p)));
        match (s.source, &ends) {
            (StrokeSource::Edge(_), [Some(a), Some(b)]) => {
                assert!(zero(&a.x) && zero(&b.x) && zero(&a.z) && zero(&b.z));
                // Counter-clockwise from y = 2 to y = -2 passes x = -2.
                assert_eq!(s.visible, a.y.to_f64() < 0.0);
            }
            (StrokeSource::Edge(_), [None, None]) => assert!(s.visible),
            (StrokeSource::Contour(_), _) => {
                assert!(s.visible);
                assert!(ends.iter().flatten().all(|p| zero(&p.x)));
            }
            other => panic!("unexpected stroke {other:?}"),
        }
    }
}

#[test]
fn cones_and_spheres_show_their_outline_and_hide_what_lies_behind() {
    let mut model = Model::new();
    let k = cone(&mut model, 3.into(), 1.into(), 4.into(), &Pose::identity()).unwrap();
    let s = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    for d in [vector(1, 0, 1), vector(1, -2, 2)] {
        let answer = steady(&model, k, &d);
        assert!(!has(&answer, VIEW_SPLIT_MARK) && !has(&answer, VIEW_RAY_MARK));
        // Top circle and both generators visible; the base circle cut where
        // the generators end, its back arc hidden.
        assert_eq!(counts(&answer.value), (4, 1, 2, 0));
        let answer = steady(&model, s, &d);
        let view = &answer.value;
        assert_eq!(counts(view), (1, 0, 1, 0));
        assert_eq!((view.strokes[0].start, view.strokes[0].end), (None, None));
    }
}

/// The torus `(|x|^2 + R^2 - r^2)^2 = 4 R^2 (x^2 + y^2)` about z.
fn torus_value(p: [f64; 3], major: f64, minor: f64) -> f64 {
    let s = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + major * major - minor * minor;
    s * s - 4.0 * major * major * (p[0] * p[0] + p[1] * p[1])
}

fn torus_gradient(p: [f64; 3], major: f64, minor: f64) -> [f64; 3] {
    let s = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + major * major - minor * minor;
    [
        4.0 * s * p[0] - 8.0 * major * major * p[0],
        4.0 * s * p[1] - 8.0 * major * major * p[1],
        4.0 * s * p[2],
    ]
}

/// An independent judge: from the true contour point nearest `p` (Newton's
/// method on `F = 0`, `grad F . d = 0`), does the ray toward the viewer pass
/// inside the torus?
fn torus_hides(p: [f64; 3], d: [f64; 3], major: f64, minor: f64) -> bool {
    let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let d = d.map(|x| x / n);
    let along = |q: [f64; 3]| {
        let g = torus_gradient(q, major, minor);
        g[0] * d[0] + g[1] * d[1] + g[2] * d[2]
    };
    let mut x = p;
    for _ in 0..20 {
        let g = torus_gradient(x, major, minor);
        let (f, c) = (torus_value(x, major, minor), along(x));
        // The gradient of c by central differences.
        let h = 1e-6;
        let gc: [f64; 3] = std::array::from_fn(|k| {
            let (mut a, mut b) = (x, x);
            a[k] += h;
            b[k] -= h;
            (along(a) - along(b)) / (2.0 * h)
        });
        // The minimum-norm Newton step on the two equations.
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let (gg, gq, qq) = (dot(g, g), dot(g, gc), dot(gc, gc));
        let det = gg * qq - gq * gq;
        let (a, b) = ((qq * f - gq * c) / det, (gg * c - gq * f) / det);
        x = std::array::from_fn(|k| x[k] - a * g[k] - b * gc[k]);
    }
    let mut t = 1e-4;
    while t < 4.0 * (major + minor) {
        let q = std::array::from_fn(|k| x[k] + t * d[k]);
        if torus_value(q, major, minor) < -1e-9 {
            return true;
        }
        t *= 1.001;
    }
    false
}

/// The binary64 middle of a fitted stroke on its carrier.
fn middle(view: &HiddenLineView, s: &ViewStroke) -> [f64; 3] {
    let Curve::BSpline(curve) = view.model.geom().curve(s.curve) else {
        panic!("a fitted stroke lies on a B-spline");
    };
    let [lo, hi] = curve.domain().map(Rational::to_f64_nearest);
    let mut dec = Decider::default();
    let mut at = |t: f64| {
        let t = Rational::from_f64_exact(t.clamp(lo, hi)).unwrap();
        close(&curve.point_at(&Real::from(t), &mut dec).unwrap())
    };
    let count = 400;
    let samples: Vec<(f64, [f64; 3])> = (0..=count)
        .map(|k| {
            let t = lo + (hi - lo) * k as f64 / count as f64;
            (t, at(t))
        })
        .collect();
    let nearest = |p: [f64; 3]| {
        let gap = |q: &[f64; 3]| (0..3).map(|k| (q[k] - p[k]).powi(2)).sum::<f64>();
        samples
            .iter()
            .min_by(|a, b| gap(&a.1).total_cmp(&gap(&b.1)))
            .unwrap()
            .0
    };
    let (Some(a), Some(b)) = (s.start, s.end) else {
        return at((lo + hi) / 2.0);
    };
    let ta = nearest(close(&view.model.geom().position(a)));
    let mut tb = nearest(close(&view.model.geom().position(b)));
    if tb <= ta {
        tb += hi - lo;
    }
    let mut t = (ta + tb) / 2.0;
    if t > hi {
        t -= hi - lo;
    }
    at(t)
}

#[test]
fn torus_outlines_agree_with_an_independent_ray_march() {
    let mut model = Model::new();
    let body = torus(&mut model, 3.into(), 1.into(), &Pose::identity()).unwrap();
    // Steep and shallow views; from low down the near tube hides part of
    // the inner outline, which turns at four cusps.
    for (d, expected) in [
        (vector(1, 0, 1), (2, 0, 2, 0)),
        (vector(1, -2, 2), (2, 0, 2, 0)),
        (vector(3, 0, 1), (2, 1, 3, 0)),
        (vector(2, -3, 1), (2, 1, 3, 0)),
    ] {
        let answer = steady(&model, body, &d);
        assert!(has(&answer, CONTOUR_FIT_MARK));
        let view = &answer.value;
        assert_eq!(counts(view), expected, "{d:?}");
        let d64 = [
            d.x.to_f64_nearest(),
            d.y.to_f64_nearest(),
            d.z.to_f64_nearest(),
        ];
        for s in &view.strokes {
            let p = middle(view, s);
            assert_eq!(!torus_hides(p, d64, 3.0, 1.0), s.visible, "{d:?} at {p:?}");
        }
    }
}

#[test]
fn p3_fillet_views_inside_the_target() {
    let (model, body) = part(include_str!("../../../parts/P3-filleted-boss.json"));
    for (name, d, expected) in [
        ("front", vector(0, -1, 0), (21, 4, 4, 4)),
        ("top", vector(0, 0, 1), (19, 0, 0, 4)),
        ("iso", vector(1, -1, 1), (21, 13, 7, 12)),
    ] {
        let answer = steady(&model, body, &d);
        assert_eq!(counts(&answer.value), expected, "{name}");
    }
}

#[test]
fn spline_faces_draw_their_rulings_and_traced_outlines() {
    // The extruded Bezier arch: an exact ruling, Bezier edges cut by
    // binary64 crossings.
    let arch = catallax_app::build_soliddraft(
        &serde_json::json!({"format":"soliddraft-recipe/1","name":"spline-profile","steps":[
            {"id":"profile","op":"profile","plane":"XY","loops":{"outer":[
                {"type":"bezier","from":[0,0],"controls":[[0,4],[4,4]],"to":[4,0]},
                {"type":"line","from":[4,0],"to":[0,0]}
            ]}},
            {"id":"body","op":"extrude","profile":"profile","distance":2}
        ],"result":"body"})
        .to_string(),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    for (d, expected) in [
        (vector(1, 1, 1), (5, 3, 1, 0)),
        (vector(-1, 2, 1), (5, 3, 1, 0)),
    ] {
        let answer = steady(&arch.built.model, arch.built.result, &d);
        assert_eq!(counts(&answer.value), expected, "arch {d:?}");
        assert!(has(&answer, VIEW_SPLIT_MARK));
    }
    // Rational patch roofs: traced, fitted generators; rays against the
    // patch are binary64 meetings, and the view says so.
    for (kind, d) in [
        (Roof::Bilinear, vector(2, 2, 1)),
        (Roof::JoinedNonlinear, vector(4, 4, 3)),
    ] {
        let fixture = solid(kind);
        let answer = steady(&fixture.model, fixture.body, &d);
        eprintln!("roof {kind:?}: {:?}", counts(&answer.value));
        assert!(has(&answer, VIEW_RAY_MARK) && has(&answer, CONTOUR_FIT_MARK));
        for s in &answer.value.strokes {
            if let StrokeSource::Contour(face) = s.source {
                assert!(fixture.roofs.contains(&face));
                assert!(s.fit.is_some());
            }
        }
    }
}

/// The ball of radius 5 about the origin less the cylinder of radius 2
/// along z through (1, 0): its two edges are closed intersection branches.
fn holed_ball() -> (Model, BodyId) {
    let mut model = Model::new();
    let ball = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let hole = cylinder(
        &mut model,
        2.into(),
        20.into(),
        &Pose::from_rational(Rotation::identity(), vector(1, 0, -10)),
    )
    .unwrap();
    let body = boolean(
        &mut model,
        ball,
        hole,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    (model, body)
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// An independent judge for the holed ball: strictly inside the ball and
/// outside the hole.
fn in_holed_ball(q: [f64; 3]) -> bool {
    dot3(q, q) < 25.0 - 1e-9 && (q[0] - 1.0).powi(2) + q[1] * q[1] > 4.0 + 1e-9
}

/// The point of the true edge `|x|^2 = 25`, `(x - 1)^2 + y^2 = 4` nearest
/// `p`, by minimum-norm Newton steps on the two equations.
fn on_branch(p: [f64; 3]) -> [f64; 3] {
    let mut x = p;
    for _ in 0..20 {
        let (f, g) = (dot3(x, x) - 25.0, (x[0] - 1.0).powi(2) + x[1] * x[1] - 4.0);
        let (gf, gg) = (
            [2.0 * x[0], 2.0 * x[1], 2.0 * x[2]],
            [2.0 * (x[0] - 1.0), 2.0 * x[1], 0.0],
        );
        let (ff, fg, qq) = (dot3(gf, gf), dot3(gf, gg), dot3(gg, gg));
        let det = ff * qq - fg * fg;
        let (a, b) = ((qq * f - fg * g) / det, (ff * g - fg * f) / det);
        x = std::array::from_fn(|k| x[k] - a * gf[k] - b * gg[k]);
    }
    x
}

/// Whether the ray from `p` toward the viewer passes inside the holed
/// ball, marched in binary64.
fn holed_ball_hides(p: [f64; 3], d: [f64; 3]) -> bool {
    let n = dot3(d, d).sqrt();
    let d = d.map(|x| x / n);
    let mut t = 1e-4;
    while t < 20.0 {
        if in_holed_ball(std::array::from_fn(|k| p[k] + t * d[k])) {
            return true;
        }
        t *= 1.001;
    }
    false
}

/// A binary64 point in the middle of a stroke on a line or circle carrier.
fn exact_middle(view: &HiddenLineView, s: &ViewStroke) -> [f64; 3] {
    let geom = view.model.geom();
    let ends = [s.start, s.end].map(|p| p.map(|p| close(&geom.position(p))));
    match geom.curve(s.curve) {
        Curve::Line(_) => {
            let [Some(a), Some(b)] = ends else {
                panic!("a line stroke has ends");
            };
            std::array::from_fn(|k| (a[k] + b[k]) / 2.0)
        }
        Curve::Circle(c) => {
            let centre = close(&c.center);
            let n = close(&Point3::new(
                c.normal.x.clone(),
                c.normal.y.clone(),
                c.normal.z.clone(),
            ));
            let n = n.map(|x| x / dot3(n, n).sqrt());
            let radius = c.radius_sq.to_f64().sqrt();
            let r = close(&Point3::new(
                c.ref_dir.x.clone(),
                c.ref_dir.y.clone(),
                c.ref_dir.z.clone(),
            ));
            let (u, sweep) = match ends {
                [Some(a), Some(b)] => {
                    let u: [f64; 3] = std::array::from_fn(|k| a[k] - centre[k]);
                    let v = [
                        n[1] * u[2] - n[2] * u[1],
                        n[2] * u[0] - n[0] * u[2],
                        n[0] * u[1] - n[1] * u[0],
                    ];
                    let w: [f64; 3] = std::array::from_fn(|k| b[k] - centre[k]);
                    let mut sweep = dot3(w, v).atan2(dot3(w, u));
                    if sweep <= 1e-12 {
                        sweep += std::f64::consts::TAU;
                    }
                    (u, sweep)
                }
                _ => {
                    let scale = radius / dot3(r, r).sqrt();
                    (r.map(|x| x * scale), 0.0)
                }
            };
            let v = [
                n[1] * u[2] - n[2] * u[1],
                n[2] * u[0] - n[0] * u[2],
                n[0] * u[1] - n[1] * u[0],
            ];
            let (sin, cos) = (sweep / 2.0).sin_cos();
            std::array::from_fn(|k| centre[k] + cos * u[k] + sin * v[k])
        }
        other => panic!("unexpected carrier {other:?}"),
    }
}

#[test]
fn intersection_branch_edges_agree_with_an_independent_ray_march() {
    let (model, body) = holed_ball();
    // The last view runs along the ball's plane of symmetry y = 0.
    for (d, symmetric) in [
        (vector(1, 0, 1), false),
        (vector(1, -2, 2), false),
        (vector(3, 0, 1), false),
        (vector(-1, 1, -2), false),
        (vector(0, 0, 1), false),
        (vector(0, -1, 0), true),
    ] {
        let answer = steady(&model, body, &d);
        assert!(has(&answer, DRAWING_FIT_MARK), "{d:?}");
        let view = &answer.value;
        let d64 = [
            d.x.to_f64_nearest(),
            d.y.to_f64_nearest(),
            d.z.to_f64_nearest(),
        ];
        let mut fitted = 0;
        for s in &view.strokes {
            let p = match (s.source, &s.fit) {
                (StrokeSource::Edge(_), Some(fit)) => {
                    fitted += 1;
                    assert!(fit.deviation_mm <= fit.target_mm, "{fit:?}");
                    assert!(!s.smooth);
                    on_branch(middle(view, s))
                }
                _ => exact_middle(view, s),
            };
            // Seen along its plane of symmetry, each branch's far half lies
            // on the sheet exactly over its near half, and the far half's
            // ray only touches the near edge, which leaves it visible. The
            // fitted edge's ray starts beside it and meets the near hole
            // wall, as OpenCASCADE's hider counts a point on a face's
            // outline as covered by the face; the far half is drawn under
            // the visible near half either way.
            if symmetric && s.fit.is_some() && p[1] > 0.0 {
                assert!(!s.visible, "{d:?} at {p:?}");
                continue;
            }
            assert_eq!(!holed_ball_hides(p, d64), s.visible, "{d:?} at {p:?}");
        }
        assert!(fitted >= 2, "{d:?}");
    }
}

/// The holed ball with a square post of half-width `1 / n` about the line
/// x = -6/5, y = 0, from z = 4 inside the ball to z = 8 above it, just
/// outside the top branch edge.
fn holed_ball_with_post(n: i64) -> (Model, BodyId, f64) {
    let (mut model, holed) = holed_ball();
    let h = Rational::from(1) / Rational::from(n);
    let w = &h + &h;
    let x = Rational::from(-6) / Rational::from(5);
    let post = block(
        &mut model,
        [w.clone(), w, 4.into()],
        &Pose::from_rational(
            Rotation::identity(),
            RVector3::new(&x - &h, -h.clone(), 4.into()),
        ),
    )
    .unwrap();
    let body = boolean(
        &mut model,
        holed,
        post,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    (model, body, h.to_f64_nearest())
}

/// Whether the open ray from `p` along `d`, which has no zero component,
/// passes inside the post of half-width `h`: the slab method on the post
/// shrunk by 1e-12, so that a ray leaving one of its faces stays out.
fn post_hides(p: [f64; 3], d: [f64; 3], h: f64) -> bool {
    let lo = [-1.2 - h, -h, 4.0];
    let hi = [-1.2 + h, h, 8.0];
    let (mut enter, mut leave) = (0.0_f64, f64::INFINITY);
    for (((lo, hi), p), d) in lo.iter().zip(hi).zip(p).zip(d) {
        let a = (lo + 1e-12 - p) / d;
        let b = (hi - 1e-12 - p) / d;
        enter = enter.max(a.min(b));
        leave = leave.min(a.max(b));
    }
    leave > enter
}

/// The point at `t` of a clamped NURBS, by de Boor's algorithm on
/// homogeneous points (Piegl and Tiller, "The NURBS Book", algorithm A4.1).
fn nurbs_point(
    degree: usize,
    knots: &[f64],
    points: &[[f64; 3]],
    weights: &[f64],
    t: f64,
) -> [f64; 3] {
    let span = (degree..points.len())
        .rev()
        .find(|&i| knots[i] <= t)
        .unwrap_or(degree);
    let weight = |i: usize| if weights.is_empty() { 1.0 } else { weights[i] };
    let mut h: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let i = span - degree + j;
            let w = weight(i);
            [points[i][0] * w, points[i][1] * w, points[i][2] * w, w]
        })
        .collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = span - degree + j;
            let alpha = (t - knots[i]) / (knots[i + degree + 1 - r] - knots[i]);
            h[j] = std::array::from_fn(|k| (1.0 - alpha) * h[j - 1][k] + alpha * h[j][k]);
        }
    }
    let q = h[degree];
    [q[0] / q[3], q[1] / q[3], q[2] / q[3]]
}

/// A binary64 point in the middle of a stroke, from its drawing curve; for a
/// fitted stroke, the true branch point nearest it.
fn drawn_middle(view: &HiddenLineView, s: &ViewStroke) -> [f64; 3] {
    let drawn = stroke_drawing(view, s, &mut Decider::default())
        .unwrap()
        .value;
    match drawn.curve {
        DrawingCurve::Line { a, b } => std::array::from_fn(|k| (a[k] + b[k]) / 2.0),
        DrawingCurve::Arc {
            center,
            normal,
            x_axis,
            radius,
            start_deg,
            end_deg,
        } => {
            let y_axis = [
                normal[1] * x_axis[2] - normal[2] * x_axis[1],
                normal[2] * x_axis[0] - normal[0] * x_axis[2],
                normal[0] * x_axis[1] - normal[1] * x_axis[0],
            ];
            let (sin, cos) = ((start_deg + end_deg) / 2.0).to_radians().sin_cos();
            std::array::from_fn(|k| center[k] + radius * (cos * x_axis[k] + sin * y_axis[k]))
        }
        DrawingCurve::Circle {
            center,
            normal,
            radius,
        } => {
            // Any point: a whole circle has one visibility.
            let x_axis = if normal[0].abs() < 0.9 {
                [0.0, normal[2], -normal[1]]
            } else {
                [-normal[2], 0.0, normal[0]]
            };
            let n = dot3(x_axis, x_axis).sqrt();
            std::array::from_fn(|k| center[k] + radius * x_axis[k] / n)
        }
        DrawingCurve::Nurbs {
            degree,
            knots,
            points,
            weights,
        } => {
            let t = (knots[degree] + knots[points.len()]) / 2.0;
            let p = nurbs_point(degree, &knots, &points, &weights, t);
            if s.fit.is_some() { on_branch(p) } else { p }
        }
    }
}

/// Every stroke of a view of the holed ball with a post agrees with the
/// independent judge at its middle.
fn agrees_with_the_judge(view: &HiddenLineView, d64: [f64; 3], h: f64) {
    for s in &view.strokes {
        let p = drawn_middle(view, s);
        let hidden = holed_ball_hides(p, d64) || post_hides(p, d64, h);
        assert_eq!(!hidden, s.visible, "{:?} at {p:?}", s.source);
    }
}

#[test]
fn a_thin_post_hides_a_thin_piece_of_a_fitted_edge() {
    // The post stands between the viewer and the top branch edge and hides
    // a piece of it a few times its half-width long. The piece's middle lies
    // a few first offsets of a ray's start (twice the fit's tolerance) from
    // its ends; the start stays within the piece, so the piece stays hidden
    // between its visible neighbours.
    let (model, body, h) = holed_ball_with_post(200_000);
    let d = vector(-4, 1, 4);
    let d64 = [-4.0, 1.0, 4.0];
    let answer = steady(&model, body, &d);
    let fit = answer
        .marks
        .iter()
        .find(|m| m.code == DRAWING_FIT_MARK)
        .expect("fitted edges")
        .bound
        .to_f64_nearest();
    assert!(!has(&answer, VIEW_SHORT_PIECE_MARK));
    let view = &answer.value;
    agrees_with_the_judge(view, d64, h);
    let geom = view.model.geom();
    let hidden_under_post: Vec<&ViewStroke> = view
        .strokes
        .iter()
        .filter(|s| s.fit.is_some() && !s.visible)
        .filter(|s| post_hides(drawn_middle(view, s), d64, h))
        .collect();
    assert_eq!(hidden_under_post.len(), 1);
    let s = hidden_under_post[0];
    let [a, b] = [s.start, s.end].map(|p| close(&geom.position(p.unwrap())));
    let length = (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt();
    assert!(length > 2.0 * h && length < 8.0 * h, "{length}");
    assert!(
        length / 2.0 > 2.0 * fit && length / 2.0 < 8.0 * fit,
        "{length} {fit}"
    );
}

#[test]
fn a_piece_too_short_for_the_fit_takes_its_neighbours_visibility() {
    // A post ten times thinner hides a piece of the top branch edge shorter
    // than the first offset that clears the fit: no exact face point beside
    // it lies within its ends. It takes the visibility of the edge around
    // it, and the view says so with a bound covering the piece.
    let (model, body, h) = holed_ball_with_post(2_000_000);
    let d = vector(-4, 1, 4);
    let d64 = [-4.0, 1.0, 4.0];
    let answer = steady(&model, body, &d);
    let fit = answer
        .marks
        .iter()
        .find(|m| m.code == DRAWING_FIT_MARK)
        .expect("fitted edges")
        .bound
        .to_f64_nearest();
    let short = answer
        .marks
        .iter()
        .find(|m| m.code == VIEW_SHORT_PIECE_MARK)
        .expect("a short piece")
        .bound
        .to_f64_nearest();
    assert!(short > 2.0 * h && short < 8.0 * fit, "{short} {fit}");
    let view = &answer.value;
    agrees_with_the_judge(view, d64, h);
    // The top branch edge is one visible stroke.
    let top: Vec<&ViewStroke> = view
        .strokes
        .iter()
        .filter(|s| s.fit.is_some() && drawn_middle(view, s)[2] > 0.0)
        .collect();
    assert_eq!(top.len(), 1);
    assert!(top[0].visible && top[0].start.is_none());
}
