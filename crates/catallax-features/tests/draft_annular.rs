//! Independent section, volume, finite-domain and display controls for a hole.
use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, Surface};
use catallax_topology::{BodyKind, Loop, Sense};
use serde_json::{Value as Json, json};
use std::collections::BTreeMap;

fn document() -> Json {
    serde_json::from_str(include_str!("fixtures/draft-annular/both.json")).unwrap()
}
fn recipe(source: &Json) -> Recipe {
    Recipe::parse(&source.to_string()).unwrap()
}
fn r(x: &str) -> Rational {
    Rational::parse(x).unwrap()
}
fn point(xyz: [&str; 3]) -> RPoint3 {
    let [x, y, z] = xyz.map(r);
    RPoint3::new(x, y, z)
}
fn both(source: &Json, cases: &[([&str; 3], Location)]) {
    let source = recipe(source);
    let mut built = evaluate(&source, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let plan = Plan::compile(&source).unwrap();
    for (xyz, expected) in cases {
        let q = plan.proven(&point(*xyz), &mut Workspace::default());
        assert_eq!(q.location, *expected, "field {xyz:?}: {q:?}");
        let q = classify_point(
            &mut built.model,
            built.result,
            &point(*xyz),
            &mut Decider::default(),
        );
        assert!(
            match expected {
                Location::Inside => matches!(q, Outcome::Proven(Containment::Inside)),
                Location::Outside => matches!(q, Outcome::Proven(Containment::Outside)),
                Location::Boundary => matches!(q, Outcome::ProvenCoincident(_)),
                _ => false,
            },
            "B-rep {xyz:?}: {q:?}"
        );
    }
}

#[test]
fn either_or_both_walls_keep_exact_sections_origins_and_volume() {
    // These are independently integrated squared linear radius laws, in pi mm3.
    for (segments, numerator) in [
        (vec!["outer"], 272),
        (vec!["inner"], 328),
        (vec!["outer", "inner"], 240),
    ] {
        let mut source = document();
        source["steps"][2]["faces"] = json!(
            segments
                .iter()
                .map(|s| json!({"face":{"step":"base","role":{"side":{"segment":s}}}}))
                .collect::<Vec<_>>()
        );
        let built = evaluate(&recipe(&source), &mut Decider::default()).unwrap();
        assert_eq!(built.model.counts(built.result).faces, 4);
        assert_eq!(built.model.counts(built.result).edges, 4);
        assert!(built.consumed.is_empty());
        for f in built.model.faces_of_body(built.result) {
            let face = built.model.face(f);
            assert_eq!(face.loops.len(), 2);
            for lp in &face.loops {
                assert!(matches!(built.model.loop_(*lp), Loop::Coedges(cs) if cs.len() == 1));
            }
            assert!(
                built
                    .model
                    .face_origins(f)
                    .iter()
                    .all(|o| o.step.as_str() == "base")
            );
            let inner = built.model.face_origins(f).iter().any(|o| matches!(&o.role,catallax_recipe::Role::Side{segment,..} if segment.as_str()=="inner"));
            if inner {
                assert_eq!(face.sense, Sense::Reversed);
            }
            match built.model.geom().surface(face.surface) {
                Surface::Plane(_) => (),
                Surface::Cone(_) => {
                    assert!(segments.contains(&if inner { "inner" } else { "outer" }))
                }
                Surface::Cylinder(_) => {
                    assert!(!segments.contains(&if inner { "inner" } else { "outer" }))
                }
                surface => panic!("unexpected surface: {surface:?}"),
            }
        }
        let volume = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let oracle =
            Ball::pi(192) * Ball::from_rational(&Rational::from_ratio(numerator, 3).unwrap(), 192);
        assert!(
            volume.value.lo <= oracle.lo() && volume.value.hi >= oracle.hi(),
            "{volume:?} {oracle:?}"
        );
        assert!(volume.value.width() <= r("1e-9"));
        let outer = if segments.contains(&"outer") {
            "7"
        } else {
            "8"
        };
        let inner = if segments.contains(&"inner") {
            "3"
        } else {
            "2"
        };
        both(
            &source,
            &[
                (["0", "0", "1"], Location::Outside),
                (["5", "0", "1"], Location::Inside),
                ([outer, "0", "1"], Location::Boundary),
                ([inner, "0", "1"], Location::Boundary),
                (["5", "0", "0"], Location::Boundary),
                (["5", "0", "2"], Location::Boundary),
                (["0", "0", "0"], Location::Outside),
                (["5", "0", "3"], Location::Outside),
            ],
        );
    }
}

#[test]
fn signed_angles_and_pull_reversal_follow_material_normals() {
    for (angle, pull, top_inner, top_outer) in [
        ("45", "1", "3", "7"),
        ("-45", "1", "1", "9"),
        ("45", "-1", "1", "9"),
    ] {
        let mut source = document();
        source["steps"][2]["neutral"]["point"][2] = json!("1");
        source["steps"][2]["angle"] = json!(angle);
        source["steps"][2]["pull"][2] = json!(pull);
        both(
            &source,
            &[
                ([top_inner, "0", "2"], Location::Boundary),
                ([top_outer, "0", "2"], Location::Boundary),
                (["5", "0", "1"], Location::Inside),
                (["0", "0", "1"], Location::Outside),
            ],
        );
        let built = evaluate(&recipe(&source), &mut Decider::default()).unwrap();
        let v = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        let oracle = Ball::pi(192) * Ball::from_rational(&r("120"), 192);
        assert!(v.lo <= oracle.lo() && v.hi >= oracle.hi());
    }
}

#[test]
fn retention_requires_material_interior_and_positive_neutral_section_area() {
    for (neutral, keep, empty) in [
        ("10", None, true),
        ("10", Some(json!([])), true),
        ("10", Some(json!([["5", "0", "1"], ["5", "0", "1"]])), false),
        ("0", None, false),
        ("2", None, false),
    ] {
        let mut source = document();
        source["steps"][2]["angle"] = json!("0");
        source["steps"][2]["neutral"]["point"][2] = json!(neutral);
        if let Some(keep) = keep {
            source["steps"][2]["keep_regions"] = keep;
        }
        let built = evaluate(&recipe(&source), &mut Decider::default()).unwrap();
        assert_eq!(
            built.model.body(built.result).kind == BodyKind::Empty,
            empty
        );
        assert_eq!(
            Plan::compile(&recipe(&source)).unwrap().is_proven_empty(),
            empty
        );
    }
    // With no moved faces, a plane tangent along a generator has no area.
    for (offset, empty) in [("8", true), ("7", false), ("0", false)] {
        let mut source = document();
        source["steps"][2]["faces"] = json!([]);
        source["steps"][2]["neutral"] = json!({"point":[offset,"0","0"],"normal":["1","0","0"]});
        let built = evaluate(&recipe(&source), &mut Decider::default()).unwrap();
        assert_eq!(
            built.model.body(built.result).kind == BodyKind::Empty,
            empty
        );
        assert_eq!(
            Plan::compile(&recipe(&source)).unwrap().is_proven_empty(),
            empty
        );
    }
    for keep in [
        json!([["0", "0", "1"]]),
        json!([["3", "0", "1"]]),
        json!([["5", "0", "0"]]),
    ] {
        let mut source = document();
        source["steps"][2]["keep_regions"] = keep;
        refuses(&source, "draft.invalid-keep-region");
    }
}

fn refuses(source: &Json, code: &str) {
    assert!(
        matches!(evaluate(&recipe(source),&mut Decider::default()),Err(BuildError::Geometry{refusal,..}) if refusal.code == code),
        "B-rep expected {code}"
    );
    assert_eq!(Plan::compile(&recipe(source)).unwrap_err().code, code);
}

#[test]
fn crossings_folds_noncoaxial_and_oblique_curves_refuse_honestly() {
    let mut source = document();
    source["steps"][1]["distance"] = json!("4");
    refuses(&source, "draft.unsupported-arrangement");
    let mut source = document();
    source["steps"][2]["angle"] = json!("-45");
    refuses(&source, "draft.unsupported-arrangement");
    let mut source = document();
    source["steps"][0]["loops"][1][0]["centre"][0] = json!("1");
    refuses(&source, "draft.unsupported-arrangement");
    let mut source = document();
    source["steps"][2]["neutral"]["normal"] = json!(["1", "0", "1"]);
    refuses(&source, "draft.unsupported-oblique-neutral");
    let mut source = document();
    source["steps"][2]["faces"] = json!([{"face":{"step":"base","role":"end"}}]);
    refuses(&source, "draft.parallel-face-normal");
}

#[test]
fn saved_scalar_edits_rebuild_radii_and_do_not_retarget_keep_points() {
    let mut source = document();
    source["steps"][0]["loops"][1][0]["radius"] = json!("3");
    both(
        &source,
        &[
            (["4", "0", "1"], Location::Boundary),
            (["3", "0", "1"], Location::Outside),
            (["5", "0", "1"], Location::Inside),
        ],
    );
    source["steps"][2]["keep_regions"] = json!([["3", "0", "1"]]);
    refuses(&source, "draft.invalid-keep-region");
}

#[test]
fn full_requested_properties_contain_independent_annular_integrals() {
    let built = evaluate(&recipe(&document()), &mut Decider::default()).unwrap();
    let properties = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(properties.marks.is_empty());
    let pi = Ball::pi(192);
    let b = |n: &str| Ball::from_rational(&r(n), 192);
    let check = |e: &catallax_measure::Enclosure, oracle: Ball| {
        assert!(
            e.lo <= oracle.lo() && e.hi >= oracle.hi(),
            "{e:?} vs {oracle:?}"
        );
        assert!(e.width() <= r("1e-7"));
    };
    // R(z)=8-z, r(z)=2+z, z in [0,2]. Integrate disks with
    // the inner disk subtracted, including the parallel-axis correction.
    check(&properties.value.volume, pi.clone() * b("80"));
    check(
        &properties.value.area,
        pi.clone() * (b("80") + b("40") * b("2").sqrt().unwrap()),
    );
    let centroid = properties.value.centroid.unwrap();
    assert!(centroid[0].contains(&r("0")) && centroid[1].contains(&r("0")));
    assert!(centroid[2].contains(&r("5/6")));
    let inertia = properties.value.inertia.unwrap();
    check(&inertia[0][0], pi.clone() * b("11020/9"));
    check(&inertia[1][1], pi.clone() * b("11020/9"));
    check(&inertia[2][2], pi * b("2400"));
    for (i, row) in inertia.iter().enumerate() {
        for (j, e) in row.iter().enumerate() {
            if i != j {
                assert!(e.contains(&r("0")));
            }
        }
    }
}

#[test]
fn algebraic_signed_angles_and_placed_scaled_cells_keep_exact_definitions() {
    for angle in ["30", "-30"] {
        let mut source = document();
        source["steps"][1]["distance"] = json!("1");
        source["steps"][2]["angle"] = json!(angle);
        both(
            &source,
            &[
                (["0", "0", "1/2"], Location::Outside),
                (["5", "0", "1/2"], Location::Inside),
            ],
        );
        let built = evaluate(&recipe(&source), &mut Decider::default()).unwrap();
        let volume = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        let b = |n: &str| Ball::from_rational(&r(n), 192);
        let correction = b("10").checked_div(&b("3").sqrt().unwrap()).unwrap();
        let oracle = (if angle == "30" {
            b("60") - correction
        } else {
            b("60") + correction
        }) * Ball::pi(192);
        assert!(volume.lo <= oracle.lo() && volume.hi >= oracle.hi());
    }
    let mut source = document();
    source["steps"][0]["plane"] =
        json!({"rotation":["1","1","0","0"],"translation":["11","-7","3"]});
    source["steps"][2]["pull"] = json!(["0", "-1", "0"]);
    source["steps"][2]["neutral"] = json!({"point":["11","-7","3"],"normal":["0","-1","0"]});
    both(
        &source,
        &[
            (["11", "-8", "3"], Location::Outside),
            (["16", "-8", "3"], Location::Inside),
            (["14", "-8", "3"], Location::Boundary),
            (["18", "-8", "3"], Location::Boundary),
        ],
    );
    let mut source = document();
    source["steps"][0]["loops"][0][0]["radius"] = json!("16");
    source["steps"][0]["loops"][1][0]["radius"] = json!("4");
    source["steps"][1]["distance"] = json!("4");
    both(
        &source,
        &[
            (["10", "0", "2"], Location::Inside),
            (["6", "0", "2"], Location::Boundary),
            (["14", "0", "2"], Location::Boundary),
        ],
    );
}

#[test]
fn brep_display_is_closed_outward_and_has_finite_annular_caps() {
    let built = evaluate(&recipe(&document()), &mut Decider::default()).unwrap();
    let options = catallax_mesh::MeshOptions {
        chord_tolerance: r("0.005"),
        ..Default::default()
    };
    let mesh = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &options,
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(mesh.faces.len(), 4);
    let mut edges = BTreeMap::new();
    for triangle in &mesh.triangles {
        let face = built.model.face(triangle.face);
        let p = triangle.vertices.map(|i| mesh.positions[i]);
        let centre: [f64; 3] = std::array::from_fn(|i| (p[0][i] + p[1][i] + p[2][i]) / 3.0);
        let a: [f64; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
        let b: [f64; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
        let n = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        if let Surface::Plane(plane) = built.model.geom().surface(face.surface) {
            let inner = if centre[2] < 1.0 { 2.0 } else { 4.0 };
            assert!(centre[0].hypot(centre[1]) >= inner - 0.005);
            assert!(n[2] * plane.normal.z.to_f64() > 0.0);
        } else {
            let radial = n[0] * centre[0] + n[1] * centre[1];
            assert!(if face.sense == Sense::Same {
                radial > 0.0
            } else {
                radial < 0.0
            });
        }
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let counts = edges.entry((a.min(b), a.max(b))).or_insert([0, 0]);
            counts[usize::from(a > b)] += 1;
        }
    }
    assert!(edges.values().all(|counts| *counts == [1, 1]));
}
