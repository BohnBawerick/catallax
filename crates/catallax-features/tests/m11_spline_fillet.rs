//! M11: a fillet on an edge between two spline faces of a loft. The M27
//! loft's first two side faces are bilinear patches; the ball rolls inside
//! the material between them, the fitted face is proved against the
//! rolling-ball definition, and the loft's caps cut it.

use std::num::NonZeroUsize;

use catallax_features::{Built, evaluate};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Curve, SplineSupport, Surface};
use catallax_topology::{EdgeId, FaceId, Loop};
use serde_json::{Value, json};

/// The M10 application allowance.
const ALLOWANCE: u64 = 20_000_000;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn polygon(prefix: &str, points: &[(&str, &str)]) -> Value {
    json!(
        points
            .iter()
            .enumerate()
            .map(|(i, (x, y))| json!({"kind":"line","id":format!("{prefix}_{i}"),"to":[x,y]}))
            .collect::<Vec<_>>()
    )
}

/// The ruled loft from the square of half-width 1 at height 0 to the
/// quadrilateral `(1,0), (0,2), (-2,0), (0,-1)` at height 2, rounded along
/// the edge between its first two sides, `(1, -1 + h, 2h)`.
fn recipe(radius: &str) -> Recipe {
    let bottom = polygon("a", &[("1", "-1"), ("1", "1"), ("-1", "1"), ("-1", "-1")]);
    let top = polygon("b", &[("1", "0"), ("0", "2"), ("-2", "0"), ("0", "-1")]);
    let section = |i: usize, height: &str, loops: Value| {
        json!({"id":format!("section{i}"),"op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0",height]},
            "loops":[loops]})
    };
    let side = |segment: &str| json!({"step":"loft","role":{"side":{"segment":segment}}});
    let value = json!({
        "format":"catallax-recipe/5","name":"m11 loft fillet",
        "steps":[
            section(0, "0", bottom),
            section(1, "2", top),
            {"id":"loft","op":"loft","shape":"ruled","sections":[
                {"profile":"section0","start":["a_0"]},
                {"profile":"section1","start":["b_0"]}]},
            {"id":"round","op":"fillet","body":"loft","radius":radius,
             "edges":[{"between":[side("a_0"), side("a_1")]}]}
        ],
        "result":"round"
    });
    Recipe::parse(&value.to_string()).unwrap()
}

fn build(radius: &str) -> (Built, u64) {
    let (built, work) = with_operation_budget(ALLOWANCE, || {
        evaluate(&recipe(radius), &mut Decider::default())
    });
    assert!(!work.exhausted, "{} units", work.used);
    (built.unwrap(), work.used)
}

/// Every face's origins and every edge's name, sorted.
fn names(built: &Built) -> (Vec<String>, Vec<String>) {
    let model = &built.model;
    let mut faces: Vec<String> = model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| model.face_origins(face).iter().map(|o| format!("{o:?}")))
        .collect();
    faces.sort();
    let mut edges: Vec<String> = model
        .edges_of_body(built.result)
        .into_iter()
        .map(|edge| format!("{:?}", model.edge_name(built.result, edge)))
        .collect();
    edges.sort();
    (faces, edges)
}

/// A bound in millimetres, rounded up to whole nanometres.
fn nanometres(upper: &Rational) -> i64 {
    (upper * &q("1000000")).ceil().to_string().parse().unwrap()
}

/// The faces using an edge.
fn users(built: &Built, edge: EdgeId) -> Vec<FaceId> {
    let model = &built.model;
    model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|face| {
            model.face(*face).loops.iter().any(|lp| {
                matches!(model.loop_(*lp), Loop::Coedges(cs)
                    if cs.iter().any(|c| model.coedge(*c).edge == edge))
            })
        })
        .collect()
}

/// The edge is replaced by a proved fillet face. Its rails trim the two
/// side faces, and each cap gains the cut that joins the rails at its
/// height.
#[test]
fn m27_loft_edge_is_rounded_between_its_caps() {
    let (built, used) = build("1/10");
    println!("M11 loft fillet work: {used}");
    let model = &built.model;
    assert_eq!(model.check(built.result), vec![]);
    // Six loft faces and the blend; one edge became two rails and two cuts,
    // and its two vertices became four rail points.
    let faces = model.faces_of_body(built.result);
    let edges = model.edges_of_body(built.result);
    assert_eq!((faces.len(), edges.len()), (7, 15));
    let blends: Vec<FaceId> = faces
        .iter()
        .copied()
        .filter(|face| {
            model
                .face_origins(*face)
                .iter()
                .any(|origin| matches!(origin.role, Role::Blend { .. }))
        })
        .collect();
    let [blend] = blends[..] else {
        panic!("one blend face, not {}", blends.len());
    };
    let Surface::BSpline(SplineSupport::SplineFillet(fillet)) =
        model.geom().surface(model.face(blend).surface)
    else {
        panic!("the blend is a spline fillet face");
    };
    // The loft's charts run `u` over a quarter of the profile each, so the
    // fit differs from the shapes tests' unit charts: 18 spans and a bound
    // of 59 nm, rounded up.
    assert_eq!(
        (fillet.spans(), nanometres(&fillet.bound().upper)),
        (18, 59)
    );
    let Loop::Coedges(boundary) = model.loop_(model.face(blend).loops[0]) else {
        panic!("the blend has an edge loop");
    };
    assert_eq!(model.face(blend).loops.len(), 1);
    let (mut rails, mut cuts) = (0, Vec::new());
    for coedge in boundary {
        let edge = model.coedge(*coedge).edge;
        let others: Vec<FaceId> = users(&built, edge)
            .into_iter()
            .filter(|face| *face != blend)
            .collect();
        let [other] = others[..] else {
            panic!("each blend edge has one other face");
        };
        match (
            model.geom().curve(model.edge(edge).curve),
            model.geom().surface(model.face(other).surface),
        ) {
            (Curve::BSpline(_), Surface::BSpline(SplineSupport::RationalAtlas(_))) => rails += 1,
            (Curve::Intersection(_), Surface::Plane(_)) => {
                // Both ends of a cut lie exactly on its cap.
                let e = model.edge(edge);
                let heights: Vec<Rational> = [e.start, e.end]
                    .into_iter()
                    .map(|v| {
                        let p = model.geom().position(model.vertex(v.unwrap()).point);
                        p.z.as_rational().expect("an exact rail point")
                    })
                    .collect();
                assert_eq!(heights[0], heights[1]);
                cuts.push(heights[0].clone());
            }
            other => panic!("unexpected blend boundary {other:?}"),
        }
    }
    cuts.sort();
    assert_eq!((rails, cuts), (2, vec![q("0"), q("2")]));
}

/// The same body at 1, 2 and 4 workers and the default, and the names
/// unchanged under a radius edit.
#[test]
fn m27_loft_fillet_matrix_and_radius_edit() {
    catallax_numbers::set_run_ahead_wait(true);
    let mut first: Option<(Vec<String>, Vec<String>, u64)> = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let (built, used) = build("1/10");
        let (faces, edges) = names(&built);
        match &first {
            Some(seen) => assert_eq!(seen, &(faces, edges, used), "{threads:?} workers"),
            None => first = Some((faces, edges, used)),
        }
    }
    let (built, used) = build("3/32");
    let model = &built.model;
    let fillet = model
        .faces_of_body(built.result)
        .into_iter()
        .find_map(
            |face| match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::SplineFillet(fillet)) => Some(fillet.clone()),
                _ => None,
            },
        )
        .expect("a fillet face");
    let edit = (fillet.spans(), nanometres(&fillet.bound().upper));
    println!("M11 loft fillet at radius 3/32: {used} units, {edit:?}");
    assert!(edit.1 <= 100);
    let (faces, edges, _) = first.unwrap();
    assert_eq!(names(&built), (faces, edges));
}
