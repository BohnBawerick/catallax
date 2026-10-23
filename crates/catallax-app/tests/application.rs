//! Application recipes, editing, picking provenance and refusal boundaries.
use catallax_app::{
    Error, build_soliddraft,
    soliddraft::{ConversionError, convert},
};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Operation, Role};
use serde_json::{Value, json};

const BRACKET: &str = include_str!("fixtures/bracket.json");

fn recipe(steps: Vec<Value>, result: &str) -> String {
    json!({"format":"soliddraft-recipe/1","name":"test","steps":steps,"result":result}).to_string()
}
fn box_steps() -> Vec<Value> {
    vec![
        json!({"id":"Sketch","op":"profile","plane":"XY","rectangle":{"center":[0,0],"width":4,"height":2}}),
        json!({"id":"Body","op":"extrude","profile":"Sketch","distance":3}),
    ]
}

#[test]
fn bracket_mesh_pick_handles_keep_merged_origins_and_hole_names_after_edit() {
    let mut input: Value = serde_json::from_str(BRACKET).unwrap();
    for height in [60, 72] {
        input["steps"][3]["distance"] = height.into();
        let result = build_soliddraft(
            &input.to_string(),
            &MeshOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let body = &result.body;
        let mesh = &body.display.value;
        let model = &body.built.model;
        let picked_origins: Vec<_> = mesh
            .triangles
            .iter()
            .flat_map(|t| &body.names.faces[&t.face])
            .collect();
        for step in [
            "plate",
            "flange",
            "boss",
            "bore_tool",
            "hole1_tool",
            "hole2_tool",
            "hole3_tool",
            "hole4_tool",
            "slot_tool",
        ] {
            assert!(
                picked_origins.iter().any(|o| o.step.as_str() == step),
                "missing picked origin {step}"
            );
        }
        assert!(
            body.names
                .faces
                .values()
                .any(|os| os.iter().any(|o| o.step.as_str() == "plate")
                    && os.iter().any(|o| o.step.as_str() == "flange")),
            "coplanar merge lost one maker"
        );
        assert_eq!(
            mesh.faces.len(),
            model.faces_of_body(body.built.result).len()
        );
        assert_eq!(
            mesh.edges.len(),
            model.edges_of_body(body.built.result).len()
        );
        for edge in &mesh.edges {
            let name = body.names.edges[&edge.edge]
                .as_ref()
                .expect("bracket edge needs origin pair");
            assert!(name.between().iter().all(|o| picked_origins.contains(&o)));
            assert!(edge.vertices.windows(2).all(|v| {
                let a = mesh.positions[v[0]];
                let b = mesh.positions[v[1]];
                (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>() > 0.0
            }));
            assert!(edge.max_chord_error <= MeshOptions::default().chord_tolerance);
        }
        let flange_top = mesh
            .triangles
            .iter()
            .find(|t| {
                body.names.faces[&t.face]
                    .iter()
                    .any(|o| o.step.as_str() == "flange" && o.role == Role::End)
            })
            .unwrap();
        assert!(
            flange_top
                .vertices
                .iter()
                .all(|v| (mesh.positions[*v][2] - f64::from(height)).abs() < 1e-10)
        );
    }
}

#[test]
fn unsupported_operations_are_typed_and_never_return_the_input_body() {
    for (op, params) in [
        ("fillet", json!({"edges":"all","radius":0.5})),
        (
            "chamfer",
            json!({"edges":{"nearest":[0,0,0]},"distance":0.5}),
        ),
        (
            "shell",
            json!({"open_faces":{"created_by":"Body"},"thickness":0.5}),
        ),
    ] {
        let mut steps = box_steps();
        let mut step = json!({"id":"finish","op":op,"body":"Body"});
        step.as_object_mut()
            .unwrap()
            .extend(params.as_object().unwrap().clone());
        steps.push(step);
        assert!(
            matches!(build_soliddraft(&recipe(steps,"finish"), &MeshOptions::default(), &mut Decider::default()),
            Err(Error::Conversion(ConversionError::NotYet { step, code:"features.unsupported-operation", .. })) if step=="finish")
        );
    }
}

#[test]
fn ids_are_collision_free_and_numbers_are_decimal_intent_not_binary_expansions() {
    let mut steps = box_steps();
    steps.push(json!({"id":"body","op":"extrude","profile":"Sketch","distance":0.1}));
    let converted = convert(&recipe(steps, "body")).unwrap();
    assert_ne!(converted.steps["Body"], converted.steps["body"]);
    assert_eq!(converted.recipe.result, converted.steps["body"]);
    match &converted.recipe.steps[2].op {
        Operation::Extrude { distance, .. } => {
            assert_eq!(distance.value(), &Rational::parse("1/10").unwrap())
        }
        _ => panic!("wrong converted operation"),
    }
    assert_eq!(
        catallax_recipe::Recipe::parse(&converted.recipe.to_text()).unwrap(),
        converted.recipe
    );
}

#[test]
fn malformed_recipes_are_not_silently_changed() {
    let mut forward = box_steps();
    forward.swap(0, 1);
    let mut unknown = box_steps();
    unknown[0]["ignored_geometry"] = 1.into();
    let mut duplicate = box_steps();
    duplicate.push(duplicate[0].clone());
    let mut two_shapes = box_steps();
    two_shapes[0]["circle"] = json!({"center":[0,0],"radius":1});
    for steps in [forward, unknown, duplicate, two_shapes] {
        assert!(matches!(
            convert(&recipe(steps, "Body")),
            Err(ConversionError::Invalid { .. })
        ));
    }
}

#[test]
fn joined_lines_and_arc_endpoints_are_made_exact_before_the_kernel() {
    let steps = vec![
        json!({"id":"p","op":"profile","plane":"XY","loops":{"outer":[
            {"type":"arc","center":[0,0],"from":[1,0],"to":[0,1.000000001]},
            {"type":"line","from":[0,1],"to":[0,0]},
            {"type":"line","from":[0,0],"to":[1,0]}
        ]}}),
        json!({"id":"b","op":"extrude","profile":"p","distance":2}),
    ];
    let converted = convert(&recipe(steps, "b")).unwrap();
    let Operation::Sketch { loops, .. } = &converted.recipe.steps[0].op else {
        panic!()
    };
    let catallax_recipe::Segment::Arc { to, centre, .. } = &loops[0].0[0] else {
        panic!()
    };
    let x = to.0[0].value() - centre.0[0].value();
    let y = to.0[1].value() - centre.0[1].value();
    assert_eq!(&x * &x + &y * &y, Rational::one());
    let built = catallax_app::build(
        &converted.recipe,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(
        built
            .names
            .faces
            .values()
            .flatten()
            .any(|o| matches!(&o.role,Role::Side { segment,.. } if segment.as_str()=="l0s0"))
    );
}

#[test]
fn full_plane_frames_keep_world_axes_and_offset() {
    for (normal, x_dir, expected) in [
        ([0, 0, 1], [1, 0, 0], [1, 0, 0, 0]),
        ([1, 0, 0], [0, 1, 0], [1, 1, 1, 1]),
        ([0, 1, 0], [0, 0, 1], [1, -1, -1, -1]),
    ] {
        let mut steps = box_steps();
        steps[0]["plane"] = json!({"origin":[1,2,3],"normal":normal,"x_dir":x_dir});
        let c = convert(&recipe(steps, "Body")).unwrap();
        let Operation::Sketch { plane, .. } = &c.recipe.steps[0].op else {
            panic!()
        };
        let q = plane.rotation.each_ref().map(|v| v.value());
        for i in 0..4 {
            assert_eq!(
                q[i] * &Rational::from_i64(expected[0]),
                q[0] * &Rational::from_i64(expected[i])
            );
        }
        assert_eq!(
            plane.translation.0.each_ref().map(|v| v.canonical()),
            ["1", "2", "3"]
        );
    }
}

#[test]
fn empty_boolean_is_a_valid_empty_display() {
    let mut steps = box_steps();
    steps
        .push(json!({"id":"empty","op":"boolean","kind":"subtract","target":"Body","tool":"Body"}));
    let result = build_soliddraft(
        &recipe(steps, "empty"),
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(result.body.display.value, catallax_mesh::Mesh::default());
    assert_eq!(result.body.names, catallax_app::Names::default());
}

#[test]
fn repeated_operands_preserve_set_identity_and_names_through_alias_steps() {
    let baseline = build_soliddraft(
        &recipe(box_steps(), "Body"),
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    for kind in ["union", "intersect"] {
        let mut steps = box_steps();
        steps.push(json!({"id":"alias","op":"boolean","kind":kind,"target":"Body","tool":"Body"}));
        steps.push(json!({"id":"again","op":"boolean","kind":kind,"target":"Body","tool":"alias"}));
        let result = build_soliddraft(
            &recipe(steps, "again"),
            &MeshOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(result.body.names, baseline.body.names);
        assert_eq!(result.body.display.value, baseline.body.display.value);
    }
}

#[test]
fn display_errors_are_separate_from_recipe_and_geometry_errors() {
    let options = MeshOptions {
        chord_tolerance: Rational::zero(),
        ..MeshOptions::default()
    };
    let error = build_soliddraft(
        &recipe(box_steps(), "Body"),
        &options,
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code(), "mesh.invalid-options");
    assert!(matches!(error, Error::Display(_)));
    let mut steps = box_steps();
    steps[1]["direction"] = json!([0, 0, 0]);
    let error = build_soliddraft(
        &recipe(steps, "Body"),
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code(), "extrude.zero-direction");
    assert!(matches!(
        error,
        Error::Build(catallax_features::BuildError::Geometry { .. })
    ));
}

#[test]
fn zx_revolution_keeps_axial_height_radial_bounds_and_segment_names() {
    let input = recipe(
        vec![
            json!({"id":"section","op":"profile","plane":"ZX","rectangle":{"center":[3,25],"width":6,"height":10}}),
            json!({"id":"collar","op":"revolve","profile":"section","axis_origin":[0,0,0],"axis_direction":[0,0,1],"angle_deg":360}),
        ],
        "collar",
    );
    let result =
        build_soliddraft(&input, &MeshOptions::default(), &mut Decider::default()).unwrap();
    let mesh = &result.body.display.value;
    assert_eq!(result.body.names.faces.len(), 4);
    for origins in result.body.names.faces.values() {
        assert_eq!(origins.len(), 1);
        assert_eq!(origins[0].step.as_str(), "collar");
        assert!(
            matches!(&origins[0].role,Role::Side { segment,.. } if ["l0s0","l0s1","l0s2","l0s3"].contains(&segment.as_str()))
        );
    }
    for p in &mesh.positions {
        assert!(p[2] >= -1e-10 && p[2] <= 6.0 + 1e-10);
        let radius_sq = p[0] * p[0] + p[1] * p[1];
        assert!((400.0 - 1e-8..=900.0 + 1e-8).contains(&radius_sq));
    }
    assert!(
        mesh.edges
            .iter()
            .all(|e| e.vertices.first() == e.vertices.last())
    );
}

#[test]
fn p3_fillet_faces_are_pickable_watertight_and_named_after_radius_edit() {
    use catallax_recipe::Recipe;
    use catallax_shapes::Surface;
    use std::collections::BTreeMap;
    let mut recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let mut previous = None;
    for edited_radius in ["1", "3/2"] {
        let Operation::Fillet { radius, .. } = &mut recipe.steps[4].op else {
            panic!("P3 must retain its selected top-edge fillets");
        };
        *radius = serde_json::from_value(json!(edited_radius)).unwrap();
        let app =
            catallax_app::build(&recipe, &MeshOptions::default(), &mut Decider::default()).unwrap();
        let model = &app.built.model;
        let mesh = &app.display.value;
        assert_eq!(
            mesh.faces.len(),
            model.faces_of_body(app.built.result).len()
        );
        let mut origins = Vec::new();
        let (mut cylinders, mut tori) = (0, 0);
        for face in &mesh.faces {
            assert!(!face.triangles.is_empty());
            assert!(face.max_chord_error <= MeshOptions::default().chord_tolerance);
            for origin in &app.names.faces[&face.face] {
                if matches!(origin.role, Role::Blend { .. }) {
                    origins.push(origin.clone());
                    match model.geom().surface(model.face(face.face).surface) {
                        Surface::Cylinder(_) => {
                            assert_eq!(origin.step.as_str(), "rounded");
                            cylinders += 1;
                        }
                        Surface::Torus(_) => {
                            assert_eq!(origin.step.as_str(), "junction");
                            tori += 1;
                        }
                        other => panic!("unexpected blend carrier {other:?}"),
                    }
                }
            }
        }
        assert_eq!((cylinders, tori), (4, 1));
        origins.sort();
        if let Some(previous) = previous {
            assert_eq!(origins, previous);
        }
        previous = Some(origins);
        let mut uses = BTreeMap::<(usize, usize), (usize, i32)>::new();
        for triangle in &mesh.triangles {
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let entry = uses.entry((a.min(b), a.max(b))).or_default();
                entry.0 += 1;
                entry.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(uses.values().all(|&count| count == (2, 0)));
        for edge in &mesh.edges {
            assert!(edge.max_chord_error <= MeshOptions::default().chord_tolerance);
            for pair in edge.vertices.windows(2) {
                assert_eq!(uses[&(pair[0].min(pair[1]), pair[0].max(pair[1]))], (2, 0));
            }
        }
    }
}
