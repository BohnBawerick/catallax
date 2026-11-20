//! Hidden-line views of the P1 bracket and the P3 filleted boss, of
//! cylinders, cones, spheres and tori, and of spline faces.

use std::num::NonZeroUsize;

use catallax_app::DEFAULT_BUILD_WORK;
use catallax_app::views::{
    CONTOUR_FIT_MARK, HiddenLineView, StrokeSource, VIEW_RAY_MARK, VIEW_SPLIT_MARK, ViewStroke,
    contour_generators, hidden_line_view,
};
use catallax_features::evaluate;
use catallax_numbers::{Answer, Decider, Rational, Real, set_threads, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, Point3, PointId, Pose, RVector3, Surface};
use catallax_topology::{BodyId, EdgeId, FaceId, Model, cone, cylinder, sphere, torus};

#[allow(dead_code)]
mod rational_patch_cases;
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
