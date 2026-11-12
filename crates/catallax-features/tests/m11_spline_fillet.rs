//! M11: a fillet on an edge between two spline faces of a loft. The M27
//! loft's first two side faces are bilinear patches; the ball rolls inside
//! the material between them, the fitted face is proved against the
//! rolling-ball definition, and the loft's caps cut it.

use std::num::NonZeroUsize;

use catallax_features::{BuildError, Built, ImportFiles, evaluate, evaluate_with_files};
use catallax_import::sha256_hex;
use catallax_numbers::{Decider, Rational, Refusal, with_operation_budget};
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

/// The loft rounded at the first edge, then at the next, between its second
/// and third sides.
fn twice(radius: &str) -> Recipe {
    let mut value = serde_json::to_value(recipe(radius)).unwrap();
    let side = |segment: &str| json!({"step":"loft","role":{"side":{"segment":segment}}});
    value["steps"].as_array_mut().unwrap().push(json!(
        {"id":"again","op":"fillet","body":"round","radius":radius,
         "edges":[{"between":[side("a_1"), side("a_2")]}]}));
    value["result"] = json!("again");
    Recipe::parse(&value.to_string()).unwrap()
}

fn build(radius: &str) -> (Built, u64) {
    let (built, work) = with_operation_budget(ALLOWANCE, || {
        evaluate(&recipe(radius), &mut Decider::default())
    });
    assert!(!work.exhausted, "{} units", work.used);
    (built.unwrap(), work.used)
}

/// The body exactly: each face's origins, surface kind and loop of edge names
/// and senses, each edge's curve and exact end points, and the fillet
/// carrier's knots, controls and weights.
fn signature(built: &Built) -> Vec<String> {
    let model = &built.model;
    let body = built.result;
    let mut out = Vec::new();
    for face in model.faces_of_body(body) {
        let f = model.face(face);
        let surface = model.geom().surface(f.surface);
        let mut line = format!("face {:?} {:?}", model.face_origins(face), f.sense);
        if let Surface::BSpline(SplineSupport::SplineFillet(fillet)) = surface {
            let carrier = fillet.carrier();
            line += &format!(
                " fillet {:?} {:?} {:?}",
                carrier.knots(),
                carrier.controls(),
                carrier.weights()
            );
        }
        for lp in &f.loops {
            if let Loop::Coedges(cs) = model.loop_(*lp) {
                for c in cs {
                    let coedge = model.coedge(*c);
                    line += &format!(
                        " {:?}{:?}",
                        model.edge_name(body, coedge.edge),
                        coedge.sense
                    );
                }
            }
        }
        out.push(line);
    }
    for edge in model.edges_of_body(body) {
        let e = model.edge(edge);
        let ends: Vec<String> = [e.start, e.end]
            .into_iter()
            .flatten()
            .map(|v| format!("{:?}", model.geom().position(model.vertex(v).point)))
            .collect();
        let curve = match model.geom().curve(e.curve) {
            Curve::BSpline(curve) => format!(
                "bspline {:?} {:?} {:?}",
                curve.knots(),
                curve.controls().collect::<Vec<_>>(),
                curve.weights().collect::<Vec<_>>()
            ),
            Curve::Intersection(_) => "intersection".into(),
            other => format!("{other:?}"),
        };
        out.push(format!(
            "edge {:?} {:?} {curve} {ends:?}",
            model.edge_name(body, edge),
            e.sense
        ));
    }
    out
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
/// height. Returns the fit's spans and its bound in nanometres, rounded up.
fn rounded(built: &Built) -> (usize, i64) {
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
    let Loop::Coedges(boundary) = model.loop_(model.face(blend).loops[0]) else {
        panic!("the blend has an edge loop");
    };
    assert_eq!(model.face(blend).loops.len(), 1);
    let (mut rails, mut cuts) = (0, Vec::new());
    for coedge in boundary {
        let edge = model.coedge(*coedge).edge;
        let others: Vec<FaceId> = users(built, edge)
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
            (
                Curve::BSpline(_),
                Surface::BSpline(SplineSupport::RationalAtlas(_) | SplineSupport::RationalPatch(_)),
            ) => rails += 1,
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
    (fillet.spans(), nanometres(&fillet.bound().upper))
}

#[test]
fn m27_loft_edge_is_rounded_between_its_caps() {
    let (built, used) = build("1/10");
    println!("M11 loft fillet work: {used}");
    // The loft's charts run `u` over a quarter of the profile each, so the
    // fit differs from the shapes tests' unit charts: 18 spans and a bound
    // of 59 nm, rounded up.
    assert_eq!(rounded(&built), (18, 59));
}

/// The second round of [`twice`]: its failing step, its whole refusal record
/// and the work used.
fn second_round(radius: &str) -> (String, Refusal, u64) {
    let (result, work) = with_operation_budget(ALLOWANCE, || {
        evaluate(&twice(radius), &mut Decider::default())
    });
    let Err(BuildError::Geometry { step, refusal }) = result else {
        panic!("the second round must refuse, not {:?}", result.map(|_| ()));
    };
    (step.as_str().to_owned(), refusal, work.used)
}

/// The same body at 1, 2 and 4 workers and the default: its exact
/// signature, names and work, and the second round's whole refusal record and
/// work. The names are unchanged under a radius edit.
#[test]
fn m27_loft_fillet_matrix_and_radius_edit() {
    catallax_numbers::set_run_ahead_wait(true);
    let mut first: Option<(Vec<String>, Vec<String>, u64)> = None;
    let mut body: Option<Vec<String>> = None;
    let mut refused: Option<(String, Refusal, u64)> = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let (built, used) = build("1/10");
        let (faces, edges) = names(&built);
        let exact = signature(&built);
        let again = second_round("1/10");
        match (&first, &body, &refused) {
            (Some(seen), Some(exactly), Some(refusal)) => {
                assert_eq!(seen, &(faces, edges, used), "{threads:?} workers");
                assert!(exactly == &exact, "{threads:?} workers: the body differs");
                assert_eq!(refusal, &again, "{threads:?} workers: the refusal differs");
            }
            _ => {
                first = Some((faces, edges, used));
                body = Some(exact);
                refused = Some(again);
            }
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

/// A second round on a face the first round trimmed is refused: the face is
/// no longer its whole patch, so a rail inside the patch's rectangle is not
/// proved inside the face.
#[test]
fn a_face_trimmed_by_an_earlier_round_is_refused() {
    let (step, refusal, _) = second_round("1/10");
    assert_eq!(
        (step.as_str(), refusal.code, refusal.detail.as_str()),
        (
            "again",
            "fillet.unsupported-spline",
            "a spline fillet's side face must be its whole patch, with no other trims; \
             maximum unknown"
        )
    );
}

/// A rounded loft written to STEP: the fillet carrier, the trimmed patch
/// faces with their rails' exact chart curves, and both cuts written as the
/// carrier composed with a chart curve, within the requested tolerance and
/// under the seam contract. The file goes to `out` when that variable is set.
fn written(built: &Built, label: &str, out: &str) {
    let mut bytes = Vec::new();
    let (report, work) = with_operation_budget(ALLOWANCE, || {
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &catallax_step::StepOptions::default(),
            &mut Decider::default(),
            &mut bytes,
        )
    });
    let report = report.expect("a STEP file").value;
    println!(
        "{label}: {} bytes, {} faces, {} edges, {} approximated curves, deviation {:.4} nm, {} units",
        bytes.len(),
        report.faces,
        report.edges,
        report.approximated_curves,
        (&report.max_deviation * &q("1000000")).to_f64_nearest(),
        work.used
    );
    if let Some(path) = std::env::var_os(out) {
        std::fs::write(path, &bytes).unwrap();
    }
    assert_eq!(
        (report.faces, report.edges, report.approximated_curves),
        (7, 15, 2)
    );
    assert!(report.max_deviation <= catallax_step::StepOptions::default().curve_tolerance);
    // Every edge of every curved face carries its chart curve there.
    catallax_step::check_seams(&bytes).expect("the seam contract");
}

#[test]
fn m27_loft_fillet_writes_step() {
    let (built, _) = build("1/10");
    written(&built, "M11 loft fillet STEP", "CATALLAX_M11_STEP_OUT");
}

/// The M27 loft as OpenCASCADE writes it: four bilinear B-spline side faces
/// between two planes. The importer's tests share the file.
const OCCT_LOFT: &[u8] = include_bytes!("../../catallax-import/tests/fixtures/m27-loft-occt.step");

/// The OpenCASCADE loft imported and rounded along the edge `(1, -1 + h, 2h)`,
/// between the side faces the file names `#277` and `#17`.
fn imported(radius: &str) -> (Built, u64) {
    let hash = sha256_hex(OCCT_LOFT);
    let face = |key: &str| json!({"step":"part","role":{"imported":{"key":key}}});
    let value = json!({
        "format":"catallax-recipe/5","name":"m11 imported fillet",
        "steps":[
            {"id":"part","op":"import_step","file":{"name":"m27-loft-occt.step","sha256":hash},
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"part","radius":radius,
             "edges":[{"between":[face("#277"), face("#17")]}]}
        ],
        "result":"round"
    });
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let files = ImportFiles::from([(hash, OCCT_LOFT)]);
    let (built, work) = with_operation_budget(ALLOWANCE, || {
        evaluate_with_files(&recipe, &files, &mut Decider::default())
    });
    assert!(!work.exhausted, "{} units", work.used);
    (built.unwrap(), work.used)
}

/// The same edge on the imported part: the importer proves that its spline
/// faces meet only along their edges, and the fillet rounds the edge between
/// two of them as on the loft body.
#[test]
fn occt_loft_edge_is_rounded_between_its_caps() {
    let (built, used) = imported("1/10");
    println!("M11 imported fillet work: {used}");
    // The file's charts are its own, so the fit differs from the loft
    // body's: 18 spans and a bound of 25 nm, rounded up.
    assert_eq!(rounded(&built), (18, 25));
}

/// The imported part's rounded body is the same at 1, 2 and 4 workers and
/// the default: its exact signature, names and work. A radius edit keeps
/// every name.
#[test]
fn occt_loft_fillet_matrix_and_radius_edit() {
    catallax_numbers::set_run_ahead_wait(true);
    type Run = (Vec<String>, (Vec<String>, Vec<String>), u64);
    let mut first: Option<Run> = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let (built, used) = imported("1/10");
        let run = (signature(&built), names(&built), used);
        match &first {
            Some(seen) => assert!(seen == &run, "{threads:?} workers: the body differs"),
            None => first = Some(run),
        }
    }
    let (built, used) = imported("3/32");
    let edit = rounded(&built);
    println!("M11 imported fillet at radius 3/32: {used} units, {edit:?}");
    assert!(edit.1 <= 100);
    assert_eq!(names(&built), first.unwrap().1);
}

#[test]
fn occt_loft_fillet_writes_step() {
    let (built, _) = imported("1/10");
    written(
        &built,
        "M11 imported fillet STEP",
        "CATALLAX_M11_IMPORTED_STEP_OUT",
    );
}

#[test]
fn m11_point_query_probe() {
    use catallax_shapes::Point3;
    for (label, (mut built, _)) in [("loft", build("1/10")), ("import", imported("1/10"))] {
        for coordinates in [
            ["0", "0", "1"],
            ["4", "0", "1"],
            ["1", "-1/2", "1"],
            ["0", "0", "0"],
        ] {
            let [x, y, z] = coordinates.map(|v| q(v).into());
            let point = Point3::new(x, y, z);
            let (answer, work) = with_operation_budget(ALLOWANCE, || {
                catallax_boolean::classify_point(
                    &mut built.model,
                    built.result,
                    &point,
                    &mut Decider::default(),
                )
            });
            println!("{label} {coordinates:?}: {answer:?}, {} units", work.used);
        }
        let faces = built.model.faces_of_body(built.result);
        let blend = faces
            .iter()
            .find_map(
                |face| match built.model.geom().surface(built.model.face(*face).surface) {
                    Surface::BSpline(SplineSupport::SplineFillet(fillet)) => {
                        Some((*face, fillet.clone()))
                    }
                    _ => None,
                },
            )
            .unwrap();
        let tau = blend.1.nodes().find(|tau| **tau > q("1/2")).unwrap();
        let rail = blend.1.rail_point(0, tau).unwrap();
        for face in faces {
            let (answer, work) = with_operation_budget(ALLOWANCE, || {
                catallax_boolean::face_contains(
                    &mut built.model,
                    face,
                    &rail,
                    &mut Decider::default(),
                )
            });
            println!("{label} rail on {face:?}: {answer:?}, {} units", work.used);
        }
    }
}
