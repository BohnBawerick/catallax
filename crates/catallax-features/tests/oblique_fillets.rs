//! Rolling-ball fillets on the ellipse where a plane cuts a cylinder obliquely.

use catallax_boolean::Containment;
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_field::Plan;
use catallax_field::measurement::{CentrePath, FaceDefinition, MeasurementOptions, measure_face};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, with_operation_log};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Curve, RPoint3, SplineSupport, Surface};
use serde_json::{Value, json};

/// A half turn about y of `2 atan(7/24)`: its z axis is `(336, 0, 527) / 625`.
const TILT: [&str; 4] = ["24/25", "0", "7/25", "0"];
const LEVEL: [&str; 4] = ["1", "0", "0", "0"];

fn cylinder(id: &str, height: &str, z: &str) -> Value {
    json!({"id":id,"op":"cylinder","radius":"5","height":height,
        "pose":{"rotation":LEVEL,"translation":["0","0",z]}})
}
/// A slab 60 wide in its own tilted frame, starting at height `at` on the axis.
fn slab(id: &str, at: &str, distance: &str) -> [Value; 2] {
    let section = format!("{id}-section");
    [
        json!({"id":section,"op":"sketch",
            "plane":{"rotation":TILT,"translation":["0","0",at]},
            "loops":[[{"kind":"line","id":"s","to":["30","-30"]},
                {"kind":"line","id":"e","to":["30","30"]},
                {"kind":"line","id":"n","to":["-30","30"]},
                {"kind":"line","id":"w","to":["-30","-30"]}]]}),
        json!({"id":id,"op":"extrude","profile":section,"distance":distance}),
    ]
}
fn edge(role: &str) -> Value {
    json!([{"between":[{"step":"wall","role":{"side":{"segment":"wall"}}},
        {"step":"cut","role":role}]}])
}
fn recipe(mut steps: Vec<Value>, radius: Option<(&str, Value)>) -> Recipe {
    let mut result = "joined";
    if let Some((radius, edges)) = radius {
        steps.push(json!({"id":"round","op":"fillet","body":"joined",
            "radius":radius,"edges":edges}));
        result = "round";
    }
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"oblique round",
            "steps":steps,"result":result})
        .to_string(),
    )
    .unwrap()
}
/// A cylinder cut off by a tilted plane through (0, 0, 8).
fn rim() -> Vec<Value> {
    let mut steps = vec![cylinder("wall", "16", "0")];
    steps.extend(slab("cut", "8", "30"));
    steps
        .push(json!({"id":"joined","op":"boolean","kind":"subtract","target":"wall","tool":"cut"}));
    steps
}
/// A boss standing on a tilted slab whose top passes through the origin.
fn boss() -> Vec<Value> {
    let mut steps = vec![cylinder("wall", "22", "-10")];
    steps.extend(slab("cut", "-25", "527/25"));
    steps.push(json!({"id":"joined","op":"boolean","kind":"union","target":"wall","tool":"cut"}));
    steps
}
/// A bore through a tilted slab whose top passes through the origin.
fn bore() -> Vec<Value> {
    let mut steps = vec![cylinder("wall", "100", "-60")];
    steps.extend(slab("cut", "-25", "527/25"));
    steps
        .push(json!({"id":"joined","op":"boolean","kind":"subtract","target":"cut","tool":"wall"}));
    steps
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn classify(built: &mut Built, points: &[[&str; 3]], want: Containment) {
    for p in points {
        let [x, y, z] = p.map(q);
        let actual = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &RPoint3::new(x, y, z),
            &mut decider(),
        );
        assert!(
            matches!(actual, Outcome::Proven(found) if found == want),
            "at {p:?}, wanted {want:?}: {actual:?}"
        );
    }
}

// The divergence theorem over the closed triangle mesh.
fn volume(mesh: &Mesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|triangle| {
            let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}
fn fine(built: &Built) -> Mesh {
    let options = MeshOptions {
        chord_tolerance: q("1/500"),
        max_vertices: 4_000_000,
        ..Default::default()
    };
    mesh_body(&built.model, built.result, &options, &mut decider())
        .unwrap()
        .value
}

fn pipes(built: &Built) -> usize {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            matches!(
                built.model.geom().surface(built.model.face(*f).surface),
                Surface::BSpline(SplineSupport::RollingPipe(_))
            )
        })
        .count()
}

fn built(steps: Vec<Value>, radius: &str, role: &str) -> Built {
    let built = evaluate(&recipe(steps, Some((radius, edge(role)))), &mut decider()).unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    assert_eq!(built.model.body(built.result).lumps.len(), 1);
    assert_eq!(built.fillets.len(), 1);
    assert_eq!(pipes(&built), 1);
    built
}
// The band's volume, from an independent section integral over the centre
// ellipse, against the change in a fine closed mesh. A mesh is only near the
// solid, so this bound is wide. The exact witnesses decide which side is which.
fn assert_change(steps: Vec<Value>, rounded: &Built, want: f64) {
    let input = evaluate(&recipe(steps, None), &mut decider()).unwrap();
    let change = volume(&fine(rounded)) - volume(&fine(&input));
    assert!((change - want).abs() < 0.03, "{change} against {want}");
}

#[test]
fn a_cut_cylinder_rim_loses_exactly_its_band() {
    for (radius, band) in [("1", 9.8242), ("2", 37.6296)] {
        let mut rounded = built(rim(), radius, "start");
        assert!(rounded.fillets.iter().all(|f| f.convex));
        assert_change(rim(), &rounded, -band);
        classify(
            &mut rounded,
            &[["0", "0", "4"], ["0", "4", "6"], ["-3", "-3", "7"]],
            Containment::Inside,
        );
        classify(&mut rounded, &[["0", "0", "9"]], Containment::Outside);
    }
    let mut rounded = built(rim(), "1", "start");
    classify(
        &mut rounded,
        &[["24/5", "0", "47/10"], ["-43/10", "0", "991/100"]],
        Containment::Inside,
    );
    classify(
        &mut rounded,
        &[["99/20", "0", "239/50"], ["-24/5", "0", "541/50"]],
        Containment::Outside,
    );
}

#[test]
fn a_boss_on_an_oblique_face_gains_exactly_its_band() {
    let mut rounded = built(boss(), "1", "end");
    assert!(rounded.fillets.iter().all(|f| !f.convex));
    assert_change(boss(), &rounded, 10.6578);
    classify(
        &mut rounded,
        &[["101/20", "0", "-79/25"], ["-26/5", "0", "71/20"]],
        Containment::Inside,
    );
    classify(
        &mut rounded,
        &[["26/5", "0", "-77/25"], ["-57/10", "0", "223/50"]],
        Containment::Outside,
    );
}

#[test]
fn a_bore_through_an_oblique_face_loses_exactly_its_lip() {
    let mut rounded = built(bore(), "1", "end");
    assert!(rounded.fillets.iter().all(|f| f.convex));
    assert_eq!(
        rounded.fillets[0]
            .same_sense
            .iter()
            .filter(|s| !**s)
            .count(),
        1
    );
    assert_change(bore(), &rounded, -10.6578);
    classify(
        &mut rounded,
        &[["26/5", "0", "-71/20"], ["-101/20", "0", "79/25"]],
        Containment::Outside,
    );
    classify(
        &mut rounded,
        &[["57/10", "0", "-223/50"], ["-26/5", "0", "77/25"]],
        Containment::Inside,
    );
}

#[test]
fn a_cached_rim_records_its_phases() {
    let prefix = evaluate_incremental(&recipe(rim(), None), None, &mut decider()).unwrap();
    let rounded = recipe(rim(), Some(("1", edge("start"))));
    let mut dec = decider();
    let (result, phases) =
        with_operation_log(|| evaluate_incremental(&rounded, Some(prefix), &mut dec));
    let mut built = result.unwrap().into_built();
    for name in ["fillet.contact", "fillet.sections", "fillet.trimming"] {
        let measured: Vec<_> = phases.iter().filter(|p| p.operation == name).collect();
        assert!(!measured.is_empty(), "{name} was not recorded");
        let total: f64 = measured.iter().map(|p| p.wall_ms).sum();
        eprintln!("{name}={total:.3}ms over {} spans", measured.len());
    }
    let start = std::time::Instant::now();
    assert_eq!(built.model.check(built.result), vec![]);
    classify(&mut built, &[["24/5", "0", "47/10"]], Containment::Inside);
    classify(
        &mut built,
        &[["99/20", "0", "239/50"]],
        Containment::Outside,
    );
    eprintln!("certification={:?}", start.elapsed());
}

#[test]
fn a_radius_edit_keeps_the_blend_name() {
    let names = |radius: &str| {
        let rounded = built(rim(), radius, "start");
        let mut found: Vec<_> = rounded
            .model
            .faces_of_body(rounded.result)
            .into_iter()
            .flat_map(|f| rounded.model.face_origins(f).to_vec())
            .filter(|o| matches!(o.role, Role::Blend { .. }))
            .collect();
        found.sort();
        found
    };
    let first = names("1");
    assert_eq!(first.len(), 1);
    assert_eq!(first, names("2"));
}

#[test]
fn a_mirrored_copy_keeps_its_band() {
    let mut steps = rim();
    steps.push(
        json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":edge("start")}),
    );
    steps.push(json!({"id":"flipped","op":"mirror","body":"round",
        "plane":{"point":["20","0","0"],"normal":["3","4","0"]},"keep_original":false}));
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"mirrored round",
            "steps":steps,"result":"flipped"})
        .to_string(),
    )
    .unwrap();
    let built = evaluate(&recipe, &mut decider()).unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    assert_eq!(pipes(&built), 1);
    let mesh = fine(&built);
    let want = 25.0 * std::f64::consts::PI * 8.0 - 9.8242;
    assert!(
        (volume(&mesh).abs() - want).abs() < 0.05,
        "{}",
        volume(&mesh)
    );
}

fn refusal(recipe: &Recipe) -> (String, String) {
    let BuildError::Geometry { refusal, .. } = evaluate(recipe, &mut decider()).unwrap_err() else {
        panic!("expected a geometric refusal");
    };
    (refusal.code.to_string(), refusal.detail.to_string())
}

#[test]
fn unsupported_oblique_contacts_stay_typed_refusals() {
    // The offset reaches the axis.
    let (code, _) = refusal(&recipe(rim(), Some(("5", edge("start")))));
    assert_eq!(code, "fillet.unsupported-contact");
    // The ball folds at the tightest curve of its centre path.
    let (code, detail) = refusal(&recipe(boss(), Some(("27", edge("end")))));
    assert_eq!(code, "fillet.unsupported-contact", "{detail}");
    // A short cylinder: the wall contact would pass the bottom cap.
    let mut short = vec![cylinder("wall", "16", "0")];
    short.extend(slab("cut", "4", "30"));
    short
        .push(json!({"id":"joined","op":"boolean","kind":"subtract","target":"wall","tool":"cut"}));
    let (code, _) = refusal(&recipe(short, Some(("2", edge("start")))));
    assert_eq!(code, "fillet.unsupported-overflow");
    // A bore inside the wall is not a plane, and nothing clears the band of it.
    let mut hollow = rim();
    hollow.push(
        json!({"id":"hole","op":"cylinder","radius":"1","height":"40",
        "pose":{"rotation":LEVEL,"translation":["0","0","-10"]}}),
    );
    hollow.push(
        json!({"id":"open","op":"boolean","kind":"subtract","target":"joined","tool":"hole"}),
    );
    hollow
        .push(json!({"id":"round","op":"fillet","body":"open","radius":"1","edges":edge("start")}));
    let hollow = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"hollow round",
            "steps":hollow,"result":"round"})
        .to_string(),
    )
    .unwrap();
    assert_eq!(refusal(&hollow).0, "fillet.unsupported-overflow");
    // Every edge includes the bottom rim, which is a separate circular ring.
    assert!(
        refusal(&recipe(rim(), Some(("20", json!("all")))))
            .0
            .starts_with("fillet.")
    );
}

#[test]
fn consumers_without_a_band_path_refuse_with_a_code() {
    let rounded = built(rim(), "1", "start");
    let volume = catallax_measure::volume(
        &rounded.model,
        rounded.result,
        &Default::default(),
        &mut decider(),
    );
    assert_eq!(volume.unwrap_err().code, "measure.unsupported-surface");
    let mut bytes = Vec::new();
    assert!(
        catallax_step::write_step(
            &rounded.model,
            &[rounded.result],
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .is_err()
    );
    assert!(bytes.is_empty());
    let mut rounded = rounded;
    let (band, pipe) = rounded
        .model
        .faces_of_body(rounded.result)
        .into_iter()
        .find_map(|f| {
            let surface = rounded.model.geom().surface(rounded.model.face(f).surface);
            match surface {
                Surface::BSpline(SplineSupport::RollingPipe(pipe)) => Some((f, pipe.clone())),
                _ => None,
            }
        })
        .unwrap();
    let spine = pipe.spine(&mut decider()).unwrap();
    let curve = rounded.model.geom_mut().add_curve(Curve::Ellipse(spine));
    let definition = FaceDefinition::Fillet {
        path: CentrePath::Arc {
            curve,
            from: None,
            to: None,
        },
        radius: pipe.radius().clone(),
    };
    let measured = measure_face(
        &mut rounded.model,
        band,
        Some(&definition),
        &MeasurementOptions::default(),
        &mut decider(),
    );
    assert_eq!(measured.unwrap_err().code, "field.face-measurement");
    let rounded = recipe(rim(), Some(("1", edge("start"))));
    let field = Plan::compile_checked(&rounded, &mut decider()).unwrap_err();
    assert_eq!(field.code, "fillet.lost-edge");
}

// A through cut splits the ellipse, so the edge has vertices and two open ends.
#[test]
fn a_partial_ellipse_stays_a_typed_refusal() {
    let mut steps = rim();
    steps.push(json!({"id":"notch-section","op":"sketch",
        "plane":{"rotation":LEVEL,"translation":["0","0","-5"]},
        "loops":[[{"kind":"line","id":"s","to":["9","-9"]},
            {"kind":"line","id":"e","to":["9","9"]},
            {"kind":"line","id":"n","to":["3","9"]},
            {"kind":"line","id":"w","to":["3","-9"]}]]}));
    steps.push(json!({"id":"notch","op":"extrude","profile":"notch-section","distance":"40"}));
    steps.push(
        json!({"id":"open","op":"boolean","kind":"subtract","target":"joined","tool":"notch"}),
    );
    steps
        .push(json!({"id":"round","op":"fillet","body":"open","radius":"1","edges":edge("start")}));
    let notched = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"notched round",
            "steps":steps,"result":"round"})
        .to_string(),
    )
    .unwrap();
    assert_eq!(refusal(&notched).0, "fillet.unsupported-topology");
}
