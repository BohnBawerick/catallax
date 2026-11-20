//! Curvature-continuous (G2) fillets (roadmap 3.4, RECIPE.md section 4.5.3)
//! through the build the application calls, `evaluate_incremental_outcome`:
//! build, save, reopen, edit, display and STEP.
//!
//! Set `CATALLAX_G2_PROOF` to a directory to keep the STEP files for the raw
//! OpenCASCADE check in `tools/step-regressions/test_g2_fillets.py`.
use std::collections::BTreeMap;

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, Built, EvaluationOutcome, evaluate_incremental_outcome};
use catallax_import::sha256_hex;
use catallax_measure::Accuracy;
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{
    Decider, OperationWork, Outcome, Rational, Real, StrictSign, with_operation_budget,
};
use catallax_recipe::{Recipe, RecipeError, Role};
use catallax_shapes::{Point3, RPoint3, SplineSupport, Surface};
use catallax_topology::RationalPatchFace;
use serde_json::{Value, json};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn side(segment: &str) -> Value {
    json!({"step": "b", "role": {"side": {"segment": segment}}})
}

fn top() -> Value {
    json!({"step": "b", "role": "end"})
}

/// A 20 by 20 by 10 block on the origin and a G2 fillet on the given edges.
fn rounded(edges: &[&str], radius: &str, fullness: Option<&str>) -> Value {
    let edges: Vec<_> = edges
        .iter()
        .map(|segment| json!({"between": [top(), side(segment)]}))
        .collect();
    let mut round = json!({"id": "round", "op": "fillet_g2", "body": "b",
        "edges": edges, "radius": radius});
    if let Some(fullness) = fullness {
        round["fullness"] = json!(fullness);
    }
    json!({"format": "catallax-recipe/5", "name": "g2 fillet", "result": "round", "steps": [
        {"id": "b", "op": "block", "size": ["20", "20", "10"],
         "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}},
        round]})
}

fn parse(doc: &Value) -> Result<Recipe, RecipeError> {
    Recipe::parse(&doc.to_string())
}

/// The application's 100M allowance, as the ordinary no-display build gives.
fn build(recipe: &Recipe) -> (Result<EvaluationOutcome, BuildError>, OperationWork) {
    with_operation_budget(100_000_000, || {
        evaluate_incremental_outcome(recipe, None, &Default::default(), &mut Decider::default())
    })
}

/// The exact body, its STEP bytes and the work of the build.
fn built(recipe: &Recipe) -> (Built, Vec<u8>, OperationWork) {
    let (outcome, work) = build(recipe);
    assert!(!work.exhausted, "{work:?}");
    let EvaluationOutcome::Exact(checkpoint) = outcome.unwrap() else {
        panic!("a G2 fillet between planes builds exactly");
    };
    let built = checkpoint.into_built();
    assert!(built.model.check(built.result).is_empty());
    let mut step = Vec::new();
    let (report, export) = with_operation_budget(100_000_000, || {
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut step,
        )
    });
    assert!(!export.exhausted, "{export:?}");
    let report = report.unwrap();
    assert!(report.marks.is_empty(), "{:?}", report.marks);
    println!(
        "{}: build work {}, STEP export work {}, STEP {} bytes, sha256 {}",
        recipe.name,
        work.used,
        export.used,
        step.len(),
        sha256_hex(&step)
    );
    (built, step, work)
}

fn keep(name: &str, bytes: &[u8]) {
    if let Some(directory) = std::env::var_os("CATALLAX_G2_PROOF") {
        std::fs::write(std::path::Path::new(&directory).join(name), bytes).unwrap();
    }
}

/// SMLib's tangent length rule (`SmBSplineCurve::FairBlendCurveDerivatives`)
/// solved as SMLib solves it, from its quadratic, for contacts `s` from the
/// edge at a material angle `angle`. The construction uses the closed form;
/// this is the independent check that the two agree.
fn smlib_tangent(s: f64, angle: f64) -> f64 {
    // The edge station at the origin, the contacts along unit rays at `angle`.
    let p0 = [s, 0.0];
    let p1 = [s * angle.cos(), s * angle.sin()];
    let v0 = [-1.0, 0.0];
    let v1 = [-angle.cos(), -angle.sin()];
    let a = [p1[0] - p0[0], p1[1] - p0[1]];
    let b = [v1[0] - v0[0], v1[1] - v0[1]];
    let qa = b[0] * b[0] + b[1] * b[1] - 1.0;
    let qb = 2.0 * (a[0] * b[0] + a[1] * b[1]);
    let qc = a[0] * a[0] + a[1] * a[1];
    let roots = if qa.abs() < 1e-14 {
        vec![-qc / qb]
    } else {
        let d = (qb * qb - 4.0 * qa * qc).sqrt();
        vec![(-qb - d) / (2.0 * qa), (-qb + d) / (2.0 * qa)]
    };
    let t = roots
        .into_iter()
        .filter(|t| *t > 0.0)
        .fold(f64::INFINITY, f64::min);
    3.0 * t
}

/// The section's six controls in the plane across a right-angled edge, with
/// the edge station at the origin and the contacts at (-s, 0) and (0, -s).
fn section(s: f64, fullness: f64) -> [[f64; 2]; 6] {
    let step = smlib_tangent(s, std::f64::consts::FRAC_PI_2) * fullness / 5.0 / s;
    [
        [-s, 0.0],
        [-s * (1.0 - step), 0.0],
        [-s * (1.0 - 2.0 * step), 0.0],
        [0.0, -s * (1.0 - 2.0 * step)],
        [0.0, -s * (1.0 - step)],
        [0.0, -s],
    ]
}

/// The area between the section and its edge station. The two straight sides
/// run through the station, so only the curve adds to the area integral. Its
/// integrand has degree 9, which five-point Gauss-Legendre integrates exactly.
fn removed_area(s: f64, fullness: f64) -> f64 {
    let c = section(s, fullness);
    let binomial = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];
    let point = |t: f64| {
        let mut p = [0.0; 2];
        for (j, control) in c.iter().enumerate() {
            let basis = binomial[j] * t.powi(j as i32) * (1.0 - t).powi(5 - j as i32);
            p[0] += basis * control[0];
            p[1] += basis * control[1];
        }
        p
    };
    let tangent = |t: f64| {
        let mut d = [0.0; 2];
        for j in 0..5 {
            let basis =
                [1.0, 4.0, 6.0, 4.0, 1.0][j] * t.powi(j as i32) * (1.0 - t).powi(4 - j as i32);
            d[0] += 5.0 * basis * (c[j + 1][0] - c[j][0]);
            d[1] += 5.0 * basis * (c[j + 1][1] - c[j][1]);
        }
        d
    };
    let nodes = [
        (0.0, 128.0 / 225.0),
        (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (0.906_179_845_938_664, 0.236_926_885_056_189_08),
        (-0.906_179_845_938_664, 0.236_926_885_056_189_08),
    ];
    let mut twice = 0.0;
    for (x, w) in nodes {
        let t = 0.5 * (x + 1.0);
        let (p, d) = (point(t), tangent(t));
        twice += 0.5 * w * (p[0] * d[1] - p[1] * d[0]);
    }
    (twice / 2.0).abs()
}

fn volume(built: &Built) -> (f64, f64) {
    let enclosure = catallax_measure::volume(
        &built.model,
        built.result,
        &Accuracy {
            relative: Rational::zero(),
            absolute: q("1e-7"),
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(enclosure.marks.is_empty(), "{:?}", enclosure.marks);
    (
        enclosure.value.lo.to_f64_nearest(),
        enclosure.value.hi.to_f64_nearest(),
    )
}

fn assert_volume(built: &Built, expected: f64) {
    let (lo, hi) = volume(built);
    println!("volume [{lo}, {hi}] mm^3, expected {expected}");
    assert!(
        lo - 1e-9 <= expected && expected <= hi + 1e-9,
        "[{lo}, {hi}] misses {expected}"
    );
}

/// Exact point membership, at the application's 100M allowance.
fn classify(built: &Built, p: [&str; 3]) -> Containment {
    let mut model = built.model.clone();
    let point = RPoint3::new(q(p[0]), q(p[1]), q(p[2]));
    let (answer, work) = with_operation_budget(100_000_000, || {
        classify_point(&mut model, built.result, &point, &mut Decider::default())
    });
    println!("{p:?}: work {}", work.used);
    match answer {
        Outcome::Proven(found) => found,
        other => panic!("{p:?}: {other:?}"),
    }
}

/// Whether a point lies exactly on one of the block's planes.
type OnPlane<'a> = &'a dyn Fn(&Point3) -> bool;

fn exactly_zero(value: &Real) -> bool {
    matches!(Decider::default().sign(value), Outcome::ProvenCoincident(_))
}

fn positive(value: &Real) -> bool {
    matches!(
        Decider::default().sign(value),
        Outcome::Proven(StrictSign::Positive)
    )
}

#[test]
fn smlib_fairing_matches_the_closed_form_at_every_angle() {
    // The construction's step 6 w sin(a/2) / (5 (1 + 2 sin(a/2))) of the
    // contact distance is SMLib's quadratic root rule, a fifth of `3 t`.
    for degrees in [20.0f64, 45.0, 60.0, 90.0, 120.0, 150.0, 175.0] {
        let angle = degrees.to_radians();
        let sine = (angle / 2.0).sin();
        let closed = 6.0 * sine / (1.0 + 2.0 * sine);
        let smlib = smlib_tangent(1.0, angle);
        assert!(
            (closed - smlib).abs() <= 1e-12,
            "{degrees}: {closed} {smlib}"
        );
    }
}

#[test]
fn a_box_edge_builds_exactly_with_curvature_continuous_contact() {
    let recipe = parse(&rounded(&["x_max"], "2", None)).unwrap();
    let (built, step, work) = built(&recipe);
    assert!(work.used <= 20_000_000, "{work:?}");
    keep("g2-box.step", &step);
    let faces = built.model.faces_of_body(built.result);
    assert_eq!(faces.len(), 7);
    let blend = faces
        .iter()
        .copied()
        .find(|f| {
            matches!(
                built.model.geom().surface(built.model.face(*f).surface),
                Surface::BSpline(SplineSupport::RationalPatch(_))
            )
        })
        .expect("one blend face");
    assert!(
        built
            .model
            .face_origins(blend)
            .iter()
            .all(|o| o.step.as_str() == "round" && matches!(o.role, Role::Blend { .. }))
    );
    let checked = RationalPatchFace::of_face(&built.model, blend, &mut Decider::default())
        .unwrap()
        .expect("a checked four-row patch face");
    // G1 and G2 at both rails, on the exact controls. Each section starts on
    // one face and ends on the other. Its first three controls lie in the
    // first face's plane, so the first and second derivatives across that
    // rail lie in it: the blend is tangent to the plane there and its normal
    // curvature is zero, as the plane's is. The last three lie in the other
    // plane. The patch is an extrusion along y, so no other second
    // derivative remains.
    let carrier = checked.support().carrier();
    assert_eq!(carrier.degrees(), [5, 1]);
    let controls = carrier.controls();
    let on_top = |p: &Point3| exactly_zero(&(&p.z - &Real::from_i64(10)));
    let on_side = |p: &Point3| exactly_zero(&(&p.x - &Real::from_i64(20)));
    let top_first = on_top(&controls[0]);
    for station in 0..2 {
        let row = |j: usize| &controls[2 * j + station];
        let (first, last): (OnPlane, OnPlane) = if top_first {
            (&on_top, &on_side)
        } else {
            (&on_side, &on_top)
        };
        for j in 0..3 {
            assert!(first(row(j)), "control {j} leaves its rail's plane");
            assert!(
                last(row(5 - j)),
                "control {} leaves its rail's plane",
                5 - j
            );
        }
        // The section turns at the middle: its third control is short of the
        // edge on the first plane, its fourth on the second. No control
        // leaves the block.
        assert!(!last(row(2)) && !first(row(3)));
        for j in 0..6 {
            assert!(!positive(&(&row(j).z - &Real::from_i64(10))));
            assert!(!positive(&(&row(j).x - &Real::from_i64(20))));
        }
    }
    // The material between the section and the edge is gone along all 20 mm.
    assert_volume(&built, 4000.0 - 20.0 * removed_area(2.0, 1.0));
    // The section runs nearer the edge than the circular fillet of the same
    // radius, whose middle is 2 - sqrt(2) = 0.586 from the edge on each axis:
    // this point, 0.52 from the edge, is material only under the G2 section.
    // The removed side is probed in the OpenCASCADE check; an exact probe
    // there projects onto the band, which this first build leaves for later.
    assert_eq!(
        classify(&built, ["19.48", "10", "9.48"]),
        Containment::Inside
    );
    assert_eq!(classify(&built, ["19.6", "10", "11"]), Containment::Outside);
    assert_eq!(classify(&built, ["10", "10", "5"]), Containment::Inside);
    // Display: every face meshed, closed and consistently oriented.
    let options = MeshOptions {
        chord_tolerance: q("0.01"),
        ..Default::default()
    };
    let mesh = mesh_body(
        &built.model,
        built.result,
        &options,
        &mut Decider::default(),
    )
    .unwrap();
    assert!(mesh.marks.is_empty(), "{:?}", mesh.marks);
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    for triangle in &mesh.value.triangles {
        let [a, b, c] = triangle.vertices;
        for (a, b) in [(a, b), (b, c), (c, a)] {
            let uses = incidence.entry((a.min(b), a.max(b))).or_default();
            uses.0 += 1;
            uses.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(
        incidence.values().all(|v| *v == (2, 0)),
        "open or misoriented display"
    );
    assert_eq!(mesh.value.faces.len(), 7);
    assert!(
        mesh.value
            .faces
            .iter()
            .all(|f| f.max_chord_error <= options.chord_tolerance)
    );
}

#[test]
fn a_g2_fillet_saves_reopens_and_edits() {
    let recipe = parse(&rounded(&["x_max"], "2", None)).unwrap();
    let (_, step, work) = built(&recipe);
    let reopened = Recipe::parse(&recipe.to_text()).unwrap();
    assert_eq!(reopened, recipe);
    let (_, again, again_work) = built(&reopened);
    assert_eq!(
        again, step,
        "a reopened recipe must write the same STEP bytes"
    );
    assert_eq!(again_work.used, work.used);
    // Edits: a larger ball, then a thinner section.
    let (larger, step, _) = built(&parse(&rounded(&["x_max"], "3", None)).unwrap());
    keep("g2-box-radius-3.step", &step);
    assert_volume(&larger, 4000.0 - 20.0 * removed_area(3.0, 1.0));
    let (thinner, step, _) = built(&parse(&rounded(&["x_max"], "3", Some("1/2"))).unwrap());
    keep("g2-box-thin.step", &step);
    assert_volume(&thinner, 4000.0 - 20.0 * removed_area(3.0, 0.5));
}

#[test]
fn separate_edges_each_get_their_band() {
    let (built, step, _) = built(&parse(&rounded(&["x_max", "x_min"], "2", Some("5/4"))).unwrap());
    keep("g2-two-edges.step", &step);
    assert_eq!(built.model.faces_of_body(built.result).len(), 8);
    assert_volume(&built, 4000.0 - 40.0 * removed_area(2.0, 1.25));
}

#[test]
fn unsupported_and_invalid_g2_fillets_refuse_by_name() {
    let refusal = |doc: &Value| {
        let (outcome, work) = build(&parse(doc).unwrap());
        assert!(!work.exhausted, "{work:?}");
        let BuildError::Geometry { refusal, .. } = outcome.unwrap_err() else {
            panic!("expected a geometric refusal");
        };
        refusal.code
    };
    // Two bands meeting at a corner need a junction this first build lacks.
    assert!(
        refusal(&rounded(&["x_max", "y_min"], "2", None)).starts_with("fillet.unsupported-corner")
    );
    // A band wider than the side face reaches the bottom face.
    assert!(refusal(&rounded(&["x_max"], "12", None)).starts_with("fillet.unsupported-overflow"));
    // Curved faces are coverage follow-on item 3, not this first build.
    let rod = json!({"format": "catallax-recipe/5", "name": "rod", "result": "round", "steps": [
        {"id": "b", "op": "cylinder", "radius": "5", "height": "10",
         "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}},
        {"id": "round", "op": "fillet_g2", "body": "b", "edges": "all", "radius": "1"}]});
    assert!(refusal(&rod).starts_with("fillet.unsupported-g2-contact"));
    // Structural rules: format 5, a positive radius, fullness in (0, 5/4].
    let mut old = rounded(&["x_max"], "2", None);
    old["format"] = json!("catallax-recipe/4");
    assert!(
        parse(&old)
            .unwrap_err()
            .to_string()
            .contains("catallax-recipe/5")
    );
    assert!(parse(&rounded(&["x_max"], "0", None)).is_err());
    assert!(parse(&rounded(&["x_max"], "2", Some("3/2"))).is_err());
    assert!(parse(&rounded(&["x_max"], "2", Some("0"))).is_err());
}
