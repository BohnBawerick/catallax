//! Equal-cylinder fillets retain exact circular sections and shared pole fans.

use catallax_boolean::Containment;
use catallax_features::{
    BuildError, certify_junction, evaluate,
    fillet::{UNSUPPORTED_CORNER, radius_limit},
};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, SaddleRing, SplineSupport, Surface};
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn recipe(kind: &str, radius: &str, second_radius: &str, gap: &str) -> Value {
    json!({"format":"catallax-recipe/3","name":"equal-radius crossing", "steps":[
        {"id":"a","op":"cylinder","radius":"8","height":"32",
         "pose":{"rotation":["1","0","1","0"],"translation":["-16","0","0"]}},
        {"id":"b","op":"cylinder","radius":second_radius,"height":"32",
         "pose":{"rotation":["1","0","0","0"],"translation":["0",gap,"-16"]}},
        {"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"b"},
        {"id":"round","op":"fillet","body":"joined","radius":radius,"edges":[
            {"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
                        {"step":"b","role":{"side":{"segment":"wall"}}}]}]}
    ],"result":"round"})
}
fn parsed(value: &Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn refusal(value: &Value) -> catallax_numbers::Refusal {
    match evaluate(&parsed(value), &mut decider()) {
        Err(BuildError::Geometry { step, refusal }) => {
            assert_eq!(step.as_str(), "round");
            refusal
        }
        other => panic!("expected a fillet refusal: {other:?}"),
    }
}

#[test]
fn equal_crossing_bands_build_and_certify_both_material_sides() {
    for kind in ["union", "intersect"] {
        let mut built = evaluate(&parsed(&recipe(kind, "1", "8", "0")), &mut decider()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let mut poles = BTreeSet::new();
        let mut bands = 0;
        for face in built.model.faces_of_body(built.result) {
            let Surface::BSpline(SplineSupport::RollingSaddle(band)) =
                built.model.geom().surface(built.model.face(face).surface)
            else {
                continue;
            };
            assert!(matches!(
                band.definition().law().ring,
                SaddleRing::Crossing { .. }
            ));
            assert_eq!(band.bound().upper, Rational::zero());
            for (rail, _) in catallax_topology::saddle_rows(&built.model, face).unwrap() {
                let edge = built.model.edge(rail);
                poles.extend([edge.start.unwrap(), edge.end.unwrap()]);
            }
            bands += 1;
        }
        assert_eq!(bands, 4);
        assert_eq!(poles.len(), 2);
        let points = if kind == "union" {
            [
                (["81/10", "0", "81/10"], Containment::Inside),
                (["17/2", "0", "17/2"], Containment::Outside),
            ]
        } else {
            [
                (["15/2", "0", "15/2"], Containment::Inside),
                (["79/10", "0", "79/10"], Containment::Outside),
            ]
        };
        for (xyz, expected) in points {
            let [x, y, z] = xyz.map(|value| Rational::parse(value).unwrap());
            let outcome = catallax_boolean::classify_point(
                &mut built.model,
                built.result,
                &RPoint3::new(x, y, z),
                &mut decider(),
            );
            assert!(
                matches!(outcome, Outcome::Proven(found) if found == expected),
                "{kind}: {outcome:?}"
            );
        }
        let certified = certify_junction(built, &mut decider()).unwrap();
        certified.verify(&mut decider()).unwrap();
        let certificate = certified.junction().unwrap();
        assert_eq!(
            certificate.record().uniform_bound.as_ref().unwrap().upper,
            Rational::zero()
        );
        assert!(
            certificate
                .record()
                .regions
                .values()
                .all(|region| region.exception.is_none())
        );
    }
}

#[test]
fn oblique_equal_crossing_keeps_complete_pole_fans() {
    for (rotation, translation, spatial) in [
        (["3", "0", "1", "0"], ["-96/5", "0", "-128/5"], false),
        (["2", "1", "1", "0"], ["-64/3", "64/3", "-32/3"], true),
    ] {
        let mut value = recipe("union", "3/2", "8", "0");
        value["steps"][0]["height"] = json!("64");
        value["steps"][0]["pose"]["translation"] = json!(["-32", "0", "0"]);
        value["steps"][1]["height"] = json!("64");
        value["steps"][1]["pose"] = json!({"rotation": rotation, "translation": translation});
        let built = evaluate(&parsed(&value), &mut decider()).unwrap();
        let certified = certify_junction(built, &mut decider()).unwrap();
        certified.verify(&mut decider()).unwrap();
        let certificate = certified.junction().unwrap();
        assert_eq!(
            certificate.record().uniform_bound.as_ref().unwrap().upper,
            Rational::zero()
        );
        if spatial {
            // Outside both original cylinders, inside the added fillet pocket.
            // The common perpendicular has length sqrt(5)/3, not a rational
            // frame that a line solver could obtain by snapping the axes.
            let [x, y, z] =
                ["2401/300", "1/150", "481/60"].map(|value| Rational::parse(value).unwrap());
            let point = RPoint3::new(x, y, z);
            let mut model = certified.built().model.clone();
            for (body, expected) in [
                (certificate.base(), Containment::Outside),
                (certificate.body(), Containment::Inside),
            ] {
                let outcome =
                    catallax_boolean::classify_point(&mut model, body, &point, &mut decider());
                assert!(
                    matches!(outcome, Outcome::Proven(found) if found == expected),
                    "{outcome:?}"
                );
            }
        }
    }
}

#[test]
fn a_crossing_clipped_at_its_input_ends_keeps_a_corner_refusal() {
    let mut value = recipe("union", "3/2", "8", "0");
    value["steps"][1]["pose"] = json!({"rotation":["3","0","1","0"],
        "translation":["-48/5","0","-64/5"]});
    let error = refusal(&value);
    assert_eq!(error.code, UNSUPPORTED_CORNER);
    assert_eq!(radius_limit(&error), None);
}

#[test]
fn different_ball_sides_and_clipped_arcs_keep_their_contact_refusals() {
    let mut clipped = recipe("union", "1", "8", "0");
    let steps = clipped["steps"].as_array_mut().unwrap();
    let mut fillet = steps.pop().unwrap();
    steps.extend([
        json!({"id":"window","op":"block","size":["40","40","2"],
            "pose":{"rotation":["1","0","0","0"],"translation":["-20","-20","2"]}}),
        json!({"id":"trimmed","op":"boolean","kind":"intersect","target":"joined","tool":"window"}),
    ]);
    fillet["body"] = json!("trimmed");
    steps.push(fillet);
    for value in [recipe("subtract", "1", "8", "0"), clipped] {
        assert_eq!(refusal(&value).code, "fillet.unsupported-contact");
    }
}

#[test]
fn equal_radius_offset_axes_and_unequal_crossing_axes_still_build() {
    for value in [
        recipe("union", "1", "8", "1"),
        recipe("union", "1", "4", "0"),
    ] {
        let mut built = evaluate(&parsed(&value), &mut decider()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        for (point, expected) in [
            (["12", "0", "0"], Containment::Inside),
            (["0", "1", "12"], Containment::Inside),
            (["12", "12", "12"], Containment::Outside),
        ] {
            let [x, y, z] = point.map(|v| Rational::parse(v).unwrap());
            let p = RPoint3::new(x, y, z);
            let outcome = catallax_boolean::classify_point(
                &mut built.model,
                built.result,
                &p,
                &mut decider(),
            );
            assert!(
                matches!(outcome, Outcome::Proven(found) if found == expected),
                "{p:?}: {outcome:?}"
            );
        }
    }
}

#[test]
fn reflected_crossing_preserves_pole_fans_and_contact_queries() {
    let mut value = recipe("union", "1", "8", "0");
    value["steps"].as_array_mut().unwrap().push(json!({
        "id":"reflected", "op":"mirror", "body":"round",
        "plane":{"point":["5","0","0"],"normal":["1","0","0"]},
        "keep_original":false
    }));
    value["result"] = json!("reflected");
    let mut built = evaluate(&parsed(&value), &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let mut bands = 0;
    for face in built.model.faces_of_body(built.result) {
        if matches!(
            built.model.geom().surface(built.model.face(face).surface),
            Surface::BSpline(SplineSupport::RollingSaddle(_))
        ) {
            assert!(catallax_topology::saddle_rows(&built.model, face).is_some());
            bands += 1;
        }
    }
    assert_eq!(bands, 4);
    for (xyz, expected) in [
        (["197/100", "1", "803/100"], Containment::Inside),
        (["3/2", "1", "17/2"], Containment::Outside),
    ] {
        let [x, y, z] = xyz.map(|value| Rational::parse(value).unwrap());
        let outcome = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &RPoint3::new(x, y, z),
            &mut decider(),
        );
        assert!(
            matches!(outcome, Outcome::Proven(found) if found == expected),
            "{outcome:?}"
        );
    }
}

#[test]
fn crossing_display_is_an_outward_manifold_through_both_poles() {
    let built = evaluate(&parsed(&recipe("union", "1", "8", "0")), &mut decider()).unwrap();
    let options = catallax_mesh::MeshOptions {
        chord_tolerance: Rational::from_ratio(6, 125).unwrap(),
        angular_deviation: Some(Rational::from_ratio(314159265, 1200000000).unwrap()),
        ..Default::default()
    };
    let mesh = catallax_mesh::mesh_body(&built.model, built.result, &options, &mut decider())
        .unwrap()
        .value;
    assert!(
        mesh.faces
            .iter()
            .all(|face| face.max_chord_error <= options.chord_tolerance)
    );
    let mut edges = Vec::with_capacity(mesh.triangles.len() * 3);
    let mut links = vec![Vec::new(); mesh.positions.len()];
    let mut volume = 0.0;
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices;
        assert!(a != b && b != c && c != a, "collapsed pole triangle");
        let [p, q, r] = [a, b, c].map(|index| mesh.positions[index]);
        let ab: [f64; 3] = std::array::from_fn(|axis| q[axis] - p[axis]);
        let ac: [f64; 3] = std::array::from_fn(|axis| r[axis] - p[axis]);
        let normal = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        assert!(normal.iter().map(|v| v * v).sum::<f64>() > 0.0);
        volume += (p[0] * (q[1] * r[2] - q[2] * r[1]) - p[1] * (q[0] * r[2] - q[2] * r[0])
            + p[2] * (q[0] * r[1] - q[1] * r[0]))
            / 6.0;
        for (a, b) in [(a, b), (b, c), (c, a)] {
            edges.push((a.min(b), a.max(b), if a < b { 1i8 } else { -1i8 }));
        }
        for (at, a, b) in [(a, b, c), (b, c, a), (c, a, b)] {
            links[at].push((a, b));
        }
    }
    edges.sort_unstable();
    for group in edges.chunk_by(|a, b| a.0 == b.0 && a.1 == b.1) {
        assert_eq!(group.len(), 2, "open or multiply used mesh edge");
        assert_eq!(group[0].2 + group[1].2, 0, "inconsistent shell winding");
    }
    for pairs in links.iter_mut().filter(|pairs| !pairs.is_empty()) {
        pairs.sort_unstable();
        assert!(pairs.windows(2).all(|pair| pair[0].0 != pair[1].0));
        let first = pairs[0].0;
        let mut at = first;
        for visited in 1..=pairs.len() {
            let index = pairs
                .binary_search_by_key(&at, |pair| pair.0)
                .expect("open vertex link");
            at = pairs[index].1;
            if at == first {
                assert_eq!(visited, pairs.len(), "disconnected pole fan");
                break;
            }
        }
        assert_eq!(at, first, "vertex link is not a cycle");
    }
    assert!(volume > 0.0, "inward shell");
}

#[test]
fn crossing_poles_have_enclosed_volume_area_and_centroid() {
    let built = evaluate(&parsed(&recipe("union", "1", "8", "0")), &mut decider()).unwrap();
    let q = |text: &str| Rational::parse(text).unwrap();
    let accuracy = catallax_measure::Accuracy {
        relative: q("1/1000000"),
        absolute: q("1/1000"),
    };
    let properties = catallax_measure::volume_area_centroid(
        &built.model,
        built.result,
        &accuracy,
        &mut decider(),
    )
    .unwrap()
    .value;
    // Independent OCCT 8 import gives 10145.886972282246 mm^3 and
    // 2981.2028800349353 mm^2. The entire native enclosure must agree.
    assert!(properties.volume.lo > q("10145.88697"));
    assert!(properties.volume.hi < q("10145.88698"));
    assert!(properties.area.lo > q("2981.20287"));
    assert!(properties.area.hi < q("2981.20289"));
    for coordinate in properties.centroid.unwrap() {
        assert!(coordinate.lo <= Rational::zero());
        assert!(coordinate.hi >= Rational::zero());
        assert!(coordinate.hi - coordinate.lo < q("1/1000000"));
    }
}
