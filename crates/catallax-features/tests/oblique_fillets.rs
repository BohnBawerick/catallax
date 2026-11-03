//! Rolling-ball fillets on the ellipse where a plane cuts a cylinder obliquely.

use catallax_boolean::Containment;
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_field::Plan;
use catallax_field::measurement::{CentrePath, FaceDefinition, MeasurementOptions, measure_face};
use catallax_measure::{
    Accuracy, Enclosure, bounding_box, face_area, face_bounding_box, face_centroid,
    mass_properties, volume,
};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real, with_operation_log};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Curve, RPoint3, SplineSupport, Surface};
use catallax_step::write_step;
use catallax_topology::FaceId;
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
        max_work: 1_000_000_000,
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
fn mesh_volume(mesh: &Mesh) -> f64 {
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
// Centroid x and z, then inertia xx, zz and xz about it, by tetrahedra.
fn mesh_moments(mesh: &Mesh) -> [f64; 5] {
    let mut size = 0.0;
    let mut first = [0.0; 3];
    let mut second = [[0.0; 3]; 3];
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
        let d = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
        size += d;
        for i in 0..3 {
            first[i] += d * (a[i] + b[i] + c[i]) / 4.0;
            for j in 0..3 {
                let sum = (a[i] + b[i] + c[i]) * (a[j] + b[j] + c[j])
                    + a[i] * a[j]
                    + b[i] * b[j]
                    + c[i] * c[j];
                second[i][j] += d * sum / 20.0;
            }
        }
    }
    let centre = first.map(|m| m / size);
    let central = |i: usize, j: usize| second[i][j] - size * centre[i] * centre[j];
    [
        centre[0],
        centre[2],
        central(1, 1) + central(2, 2),
        central(0, 0) + central(1, 1),
        -central(0, 2),
    ]
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
    let change = mesh_volume(&fine(rounded)) - mesh_volume(&fine(&input));
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
        (mesh_volume(&mesh).abs() - want).abs() < 0.05,
        "{}",
        mesh_volume(&mesh)
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
    assert_eq!(code, "fillet.radius-too-large");
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

// A chart quadrature of the rolling-ball law, sharing no formula with the
// kernel: band area, then its centroid.
const BANDS: [(&str, &str, [f64; 4]); 4] = [
    ("1", "rim", [49.51037253, -28.788419, 0.0, 389.102937]),
    ("2", "rim", [89.97352207, -31.109023, 0.0, 642.068250]),
    ("1", "boss", [58.55759553, -67.329841, 0.0, 70.154621]),
    ("1", "bore", [58.55759553, 67.329841, 0.0, -70.154621]),
];

fn fine_accuracy() -> Accuracy {
    Accuracy {
        relative: q("1/100000000"),
        absolute: q("1/100000000"),
    }
}
fn mid(value: &Enclosure) -> f64 {
    ((&value.lo + &value.hi) * q("1/2")).to_f64_nearest()
}
fn band_of(built: &Built) -> FaceId {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|f| {
            matches!(
                built.model.geom().surface(built.model.face(*f).surface),
                Surface::BSpline(SplineSupport::RollingPipe(_))
            )
        })
        .unwrap()
}
fn recipe_of(kind: &str) -> (Vec<Value>, &'static str) {
    match kind {
        "rim" => (rim(), "start"),
        "boss" => (boss(), "end"),
        _ => (bore(), "end"),
    }
}

#[test]
fn a_band_measures_as_the_independent_quadrature() {
    for (radius, kind, [area, x, y, z]) in BANDS {
        let (steps, role) = recipe_of(kind);
        let built = built(steps, radius, role);
        let band = band_of(&built);
        let found = face_area(&built.model, band, &fine_accuracy(), &mut decider()).unwrap();
        assert!(found.marks.is_empty(), "{kind}");
        assert!(found.value.width() <= q("1/10000000"), "{kind}");
        assert!((mid(&found.value) - area).abs() < 1e-7, "{kind}");
        let centre = face_centroid(&built.model, band, &fine_accuracy(), &mut decider()).unwrap();
        for (axis, want) in centre.value.iter().zip([x, y, z]) {
            assert!((mid(axis) - want / area).abs() < 1e-5, "{kind}");
        }
    }
}

// The volume change against the section integral, and the mass properties
// against the divergence theorem on a fine closed mesh. A mesh is only near the
// solid, so its bounds are wide.
#[test]
fn a_rounded_solid_has_exact_mass_properties() {
    for (radius, kind, change) in [
        ("1", "rim", -9.8242),
        ("2", "rim", -37.6296),
        ("1", "boss", 10.6578),
        ("1", "bore", -10.6578),
    ] {
        let (steps, role) = recipe_of(kind);
        let input = evaluate(&recipe(steps.clone(), None), &mut decider()).unwrap();
        let built = built(steps, radius, role);
        let whole = volume(&built.model, built.result, &fine_accuracy(), &mut decider());
        let before = volume(&input.model, input.result, &fine_accuracy(), &mut decider());
        let (whole, before) = (whole.unwrap(), before.unwrap());
        assert!(whole.marks.is_empty(), "{kind}");
        let found = mid(&whole.value) - mid(&before.value);
        assert!((found - change).abs() < 1e-3, "{kind}: {found} {change}");
        let mesh = fine(&built);
        let properties =
            mass_properties(&built.model, built.result, &fine_accuracy(), &mut decider()).unwrap();
        let properties = properties.value;
        assert!(
            (mid(&properties.volume) - mesh_volume(&mesh)).abs() < 0.05,
            "{kind}"
        );
        let centre = properties.centroid.unwrap();
        let tensor = properties.inertia.unwrap();
        let [cx, cz, xx, zz, xz] = mesh_moments(&mesh);
        assert!((mid(&centre[0]) - cx).abs() < 1e-2, "{kind}");
        assert!((mid(&centre[2]) - cz).abs() < 1e-2, "{kind}");
        for (found, want) in [
            (&tensor[0][0], xx),
            (&tensor[2][2], zz),
            (&tensor[0][2], xz),
        ] {
            let found = mid(found);
            assert!(
                (found - want).abs() < 2e-3 * xx.abs(),
                "{kind}: {found} {want}"
            );
        }
    }
}

// The box's enclosures have no slack to speak of, so a fine mesh sits inside it
// and reaches within the mesh's own chord of each side.
#[test]
fn a_rounded_solid_has_an_exact_box() {
    for (radius, kind, _) in BANDS {
        let (steps, role) = recipe_of(kind);
        let built = built(steps, radius, role);
        let found = bounding_box(&built.model, built.result, &fine_accuracy(), &mut decider());
        let found = found.unwrap();
        assert!(found.marks.is_empty(), "{kind}");
        let band = band_of(&built);
        let face = face_bounding_box(&built.model, band, &fine_accuracy(), &mut decider());
        let face = face.unwrap().value;
        let whole = found.value.unwrap();
        let mesh = fine(&built);
        for axis in 0..3 {
            let (low, high) = mesh
                .positions
                .iter()
                .map(|p| p[axis])
                .fold((f64::MAX, f64::MIN), |(a, b), v| (a.min(v), b.max(v)));
            let (box_low, box_high) = (
                whole.min[axis].to_f64_nearest(),
                whole.max[axis].to_f64_nearest(),
            );
            assert!(
                (low - box_low).abs() < 0.01,
                "{kind} {axis}: {low} {box_low}"
            );
            assert!(
                (high - box_high).abs() < 0.01,
                "{kind} {axis}: {high} {box_high}"
            );
            assert!(face.min[axis] >= whole.min[axis], "{kind}");
            assert!(face.max[axis] <= whole.max[axis], "{kind}");
        }
    }
}

#[test]
fn measuring_a_band_against_its_definition_is_exact() {
    let mut rounded = built(rim(), "1", "start");
    let band = band_of(&rounded);
    let Surface::BSpline(SplineSupport::RollingPipe(pipe)) = rounded
        .model
        .geom()
        .surface(rounded.model.face(band).surface)
        .clone()
    else {
        unreachable!()
    };
    let spine = pipe.spine(&mut decider()).unwrap();
    let curve = rounded.model.geom_mut().add_curve(Curve::Ellipse(spine));
    let mut measured = |radius: Real| {
        let definition = FaceDefinition::Fillet {
            path: CentrePath::Arc {
                curve,
                from: None,
                to: None,
            },
            radius,
        };
        measure_face(
            &mut rounded.model,
            band,
            Some(&definition),
            &MeasurementOptions::default(),
            &mut decider(),
        )
    };
    let exact = measured(pipe.radius().clone()).unwrap();
    assert_eq!(exact.upper, Rational::zero());
    assert!(exact.converged);
    let wrong = measured(Real::from(q("5/4"))).unwrap();
    assert_eq!(wrong.lower, q("1/4"));
    assert_eq!(wrong.upper, q("1/4"));
}

#[test]
fn the_field_still_refuses_an_oblique_band() {
    let rounded = recipe(rim(), Some(("1", edge("start"))));
    let field = Plan::compile_checked(&rounded, &mut decider()).unwrap_err();
    assert_eq!(field.code, "fillet.lost-edge");
}

#[test]
fn a_rounded_solid_writes_step() {
    for (radius, kind, _) in BANDS {
        let (steps, role) = recipe_of(kind);
        let rounded = built(steps, radius, role);
        let mut bytes = Vec::new();
        let answer = write_step(
            &rounded.model,
            &[rounded.result],
            &Default::default(),
            &mut decider(),
            &mut bytes,
        )
        .unwrap();
        assert!(answer.marks.is_empty(), "{kind}");
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("RATIONAL_B_SPLINE_SURFACE"), "{kind}");
        assert!(text.contains("SEAM_CURVE"), "{kind}");
    }
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
    // The notch cuts the band open, and its ends are not built yet.
    assert_eq!(refusal(&notched).0, "fillet.unsupported-end");
}
