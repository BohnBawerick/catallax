//! Rolling-ball fillets on the ellipse where a plane cuts a cone obliquely.

use catallax_boolean::Containment;
use catallax_features::{Built, evaluate};
use catallax_measure::{Accuracy, Enclosure, volume, volume_area_centroid};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, StrictSign};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, RPoint3, SplineSupport, Surface};
use catallax_step::write_step;
use serde_json::{Value, json};

/// A half turn about y of `2 atan(7/24)`: its z axis is `(336, 0, 527) / 625`.
const TILT: [&str; 4] = ["24/25", "0", "7/25", "0"];
const LEVEL: [&str; 4] = ["1", "0", "0", "0"];

// The corpus's frustum of radii 10 and 4 over height 24, cut off above the
// tilted plane through (0, 0, 8).
fn cut_cone(radius: Option<&str>) -> Recipe {
    recipe(radius, false)
}
// The same, mirrored in x = 0 after the round, keeping only the copy.
fn mirrored_cut_cone() -> Recipe {
    recipe(Some("1/2"), true)
}
fn recipe(radius: Option<&str>, mirror: bool) -> Recipe {
    let mut steps = vec![
        json!({"id":"wall","op":"cone","base_radius":"10","top_radius":"4","height":"24",
            "pose":{"rotation":LEVEL,"translation":["0","0","0"]}}),
        json!({"id":"section","op":"sketch",
            "plane":{"rotation":TILT,"translation":["0","0","8"]},
            "loops":[[{"kind":"line","id":"s","to":["30","-30"]},
                {"kind":"line","id":"e","to":["30","30"]},
                {"kind":"line","id":"n","to":["-30","30"]},
                {"kind":"line","id":"w","to":["-30","-30"]}]]}),
        json!({"id":"cut","op":"extrude","profile":"section","distance":"30"}),
        json!({"id":"joined","op":"boolean","kind":"subtract","target":"wall","tool":"cut"}),
    ];
    let mut result = "joined";
    if let Some(radius) = radius {
        let edges: Value = json!([{"between":[{"step":"wall","role":{"side":{"segment":"wall"}}},
            {"step":"cut","role":"start"}]}]);
        steps.push(json!({"id":"round","op":"fillet","body":"joined",
            "radius":radius,"edges":edges}));
        result = "round";
    }
    if mirror {
        steps.push(json!({"id":"reflected","op":"mirror","body":result,
            "plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}));
        result = "reflected";
    }
    Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"oblique cone round",
            "steps":steps,"result":result})
        .to_string(),
    )
    .unwrap()
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
fn bands(built: &Built) -> Vec<catallax_shapes::RollingPipe> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter_map(
            |f| match built.model.geom().surface(built.model.face(f).surface) {
                Surface::BSpline(SplineSupport::RollingPipe(pipe)) => Some(pipe.clone()),
                _ => None,
            },
        )
        .collect()
}
fn measured(built: &Built) -> Enclosure {
    let accuracy = Accuracy {
        relative: q("1/100000000"),
        absolute: q("1/100000000"),
    };
    let answer = volume(&built.model, built.result, &accuracy, &mut decider()).unwrap();
    assert!(answer.marks.is_empty());
    answer.value
}

// The volume a fine display mesh encloses, after checking that every mesh edge
// is used once in each direction, so the mesh is closed and consistently wound.
fn closed_mesh_volume(built: &Built) -> f64 {
    let mesh = mesh_body(
        &built.model,
        built.result,
        &MeshOptions {
            chord_tolerance: q("1/500"),
            max_vertices: 4_000_000,
            ..Default::default()
        },
        &mut decider(),
    )
    .unwrap()
    .value;
    let mut uses = std::collections::BTreeMap::new();
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices;
        for (from, to) in [(a, b), (b, c), (c, a)] {
            *uses.entry((from, to)).or_insert(0) += 1;
        }
    }
    for (&(from, to), &count) in &uses {
        assert_eq!(count, 1, "edge {from}-{to} repeats its direction");
        assert_eq!(uses.get(&(to, from)), Some(&1), "edge {from}-{to} is open");
    }
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
// The sign of the plane rail's ellipse normal against the band's chart normal,
// which fixes the direction mesh and STEP give that rail.
fn rail_turn(built: &Built) -> StrictSign {
    let pipes = bands(built);
    assert_eq!(pipes.len(), 1);
    let turn = pipes[0].chart_normal();
    let mut found = None;
    for edge in built.model.edges_of_body(built.result) {
        if let Curve::Ellipse(rail) = built.model.geom().curve(built.model.edge(edge).curve) {
            let Outcome::Proven(sign) = decider().sign(&rail.normal.dot(&turn)) else {
                panic!("the plane rail's turn is not separated from zero");
            };
            assert!(found.replace(sign).is_none(), "one plane rail");
        }
    }
    found.expect("a plane rail")
}

// Steiner's formula on the cut frustum's opening, with the bottom rim's band
// added back, gives 1957.4896872686486130 for a radius of 1/2. The same
// derivation reproduces Pappus's value for a square cut.
#[test]
fn a_cut_cone_rim_loses_exactly_its_band() {
    let mut rounded = evaluate(&cut_cone(Some("1/2")), &mut decider()).unwrap();
    assert_eq!(rounded.model.check(rounded.result), vec![]);
    assert_eq!(rounded.model.body(rounded.result).lumps.len(), 1);
    assert_eq!(rounded.fillets.len(), 1);
    assert!(rounded.fillets[0].convex);
    let pipes = bands(&rounded);
    assert_eq!(pipes.len(), 1);
    assert!(pipes[0].is_conical());
    // The plane keeps an exact ellipse; the cone's contact is the band's row.
    let contacts = &rounded.fillets[0].contacts;
    assert_eq!(
        contacts
            .iter()
            .filter(|c| matches!(c, Curve::Ellipse(_)))
            .count(),
        1
    );
    assert_eq!(
        contacts
            .iter()
            .filter(|c| matches!(c, Curve::BSpline(_)))
            .count(),
        1
    );
    // The exact value to twenty places, far inside the enclosure's width.
    let exact = q("1957.48968726864861303");
    let Enclosure { lo, hi } = measured(&rounded);
    assert!(lo <= exact && exact <= hi, "{lo:?} {hi:?}");
    assert!(&hi - &lo <= q("1/10000"));
    let properties = volume_area_centroid(
        &rounded.model,
        rounded.result,
        &Accuracy {
            relative: q("1/1000000"),
            absolute: q("1/1000000"),
        },
        &mut decider(),
    )
    .unwrap()
    .value;
    assert!(properties.volume.lo <= exact && exact <= properties.volume.hi);
    // The frustum is symmetric across y = 0, and the cut tilts towards +x
    // where it sits lowest, which moves the centroid to -x.
    let centre = properties.centroid.unwrap();
    assert!(centre[1].lo <= Rational::zero() && Rational::zero() <= centre[1].hi);
    assert!(centre[0].hi < Rational::zero());
    // A fine display mesh closes, and each of its points lies within the chord
    // tolerance of the surface, so its volume misses by at most about the
    // area times that tolerance.
    let shown = closed_mesh_volume(&rounded);
    let bound = properties.area.hi.to_f64_nearest() / 500.0;
    assert!((shown - exact.to_f64_nearest()).abs() < bound, "{shown}");
    classify(
        &mut rounded,
        &[
            ["0", "0", "4"],
            ["0", "5", "5"],
            ["-5", "-3", "6"],
            // Inside the ball as it rolls past the rim's lowest and highest
            // points.
            ["47/5", "0", "37/20"],
            ["-673/100", "0", "1197/100"],
        ],
        Containment::Inside,
    );
    classify(
        &mut rounded,
        &[
            ["0", "0", "9"],
            ["10", "0", "1"],
            ["0", "0", "-1/10"],
            // Inside the sharp rim but outside the ball: the blend cut these.
            ["4749/500", "0", "1919/1000"],
            ["-6821/1000", "0", "12201/1000"],
        ],
        Containment::Outside,
    );
}

#[test]
fn a_cut_cone_rim_writes_its_band_to_step() {
    let rounded = evaluate(&cut_cone(Some("1/2")), &mut decider()).unwrap();
    let mut out = Vec::new();
    write_step(
        &rounded.model,
        &[rounded.result],
        &Default::default(),
        &mut decider(),
        &mut out,
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.contains("CONICAL_SURFACE"));
    assert!(text.contains("B_SPLINE_SURFACE"));
}

#[test]
fn a_radius_edit_keeps_the_blend_name() {
    let names = |radius: &str| {
        let rounded = evaluate(&cut_cone(Some(radius)), &mut decider()).unwrap();
        let mut found: Vec<_> = rounded
            .model
            .faces_of_body(rounded.result)
            .into_iter()
            .flat_map(|f| rounded.model.face_origins(f).to_vec())
            .filter(|o| matches!(o.role, catallax_recipe::Role::Blend { .. }))
            .collect();
        found.sort();
        found
    };
    let first = names("1/2");
    assert_eq!(first.len(), 1);
    assert_eq!(first, names("1"));
}

// A mirror after the round reverses the band's frame. The chart must turn with
// the copied carrier, so the plane rail keeps its sense against the chart and
// the mesh still closes with one winding.
#[test]
fn a_mirrored_cut_cone_rim_keeps_its_rail_sense() {
    let original = evaluate(&cut_cone(Some("1/2")), &mut decider()).unwrap();
    let mut mirrored = evaluate(&mirrored_cut_cone(), &mut decider()).unwrap();
    assert_eq!(mirrored.model.check(mirrored.result), vec![]);
    assert_eq!(rail_turn(&mirrored), rail_turn(&original));
    let exact = q("1957.48968726864861303");
    let Enclosure { lo, hi } = measured(&mirrored);
    assert!(lo <= exact && exact <= hi, "{lo:?} {hi:?}");
    let shown = closed_mesh_volume(&mirrored);
    assert!((shown - exact.to_f64_nearest()).abs() < 1.0, "{shown}");
    let mut out = Vec::new();
    write_step(
        &mirrored.model,
        &[mirrored.result],
        &Default::default(),
        &mut decider(),
        &mut out,
    )
    .unwrap();
    // The original's witnesses, mirrored in x = 0.
    classify(
        &mut mirrored,
        &[
            ["0", "0", "4"],
            ["-47/5", "0", "37/20"],
            ["673/100", "0", "1197/100"],
        ],
        Containment::Inside,
    );
    classify(
        &mut mirrored,
        &[
            ["0", "0", "9"],
            ["-4749/500", "0", "1919/1000"],
            ["6821/1000", "0", "12201/1000"],
        ],
        Containment::Outside,
    );
}
