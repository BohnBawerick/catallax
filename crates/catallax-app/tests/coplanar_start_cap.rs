//! Authored sketch joints must survive arc lifting and opposite sheets must not publish.
use catallax_app::soliddraft::convert;
use catallax_boolean::{Containment, check_body, classify_point};
use catallax_measure::{Accuracy, face_area};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::Role;
use catallax_shapes::RPoint3;
use catallax_topology::{BodyId, Model};
use serde_json::json;

#[test]
fn released_application_sketch_opens_the_coplanar_top_without_a_lid() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let mut application = catallax_app::build_soliddraft(
            include_str!("fixtures/coplanar-start-cap.json"),
            &MeshOptions::default(),
            &mut dec,
        )
        .unwrap();
        let built = &mut application.body.built;
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 13);
        let flat: Vec<_> = faces.iter().copied().filter(|&face| {
            built.model.face_origins(face).iter().any(|origin| {
                origin.step == application.steps["op_3"]
                    && matches!(&origin.role, Role::Side { segment, .. } if segment.as_str() == "l0s3")
            })
        }).collect();
        assert_eq!(flat.len(), 1);
        let accuracy = Accuracy {
            relative: Rational::zero(),
            absolute: Rational::parse("0.000001").unwrap(),
        };
        let area = face_area(&built.model, flat[0], &accuracy, &mut dec)
            .unwrap()
            .value;
        assert!(
            area.contains(&Rational::parse("1894.61").unwrap()),
            "{area:?}"
        );
        assert!(area.width() <= accuracy.absolute);
        assert!(!faces.iter().any(|&face| {
            built.model.face_origins(face).iter().any(|origin| {
                origin.step == application.steps["op_45"] && origin.role == Role::Start
            })
        }));
        let mouth = RPoint3::new(0.into(), (-38).into(), 60.into());
        assert!(matches!(
            classify_point(&mut built.model, built.result, &mouth, &mut dec),
            Outcome::Proven(Containment::Outside)
        ));
        assert!(
            check_body(&built.model, built.result, &mut dec)
                .unwrap()
                .is_empty()
        );
    });
    assert!(!work.exhausted, "{work:?}");
}

fn two_boxes(offset: [i32; 3]) -> (Model, BodyId) {
    let source = json!({
        "format":"soliddraft-recipe/1", "name":"opposite sheets", "result":"a",
        "steps":[
            {"id":"p","op":"profile","plane":"XY","rectangle":{"center":[0,0],"width":4,"height":4}},
            {"id":"a","op":"extrude","profile":"p","distance":1},
            {"id":"q","op":"profile","plane":{"origin":offset,"normal":[0,0,1],"x_dir":[1,0,0]},"rectangle":{"center":[0,0],"width":4,"height":4}},
            {"id":"b","op":"extrude","profile":"q","distance":1}
        ]
    });
    let converted = convert(&source.to_string()).unwrap();
    let built = catallax_features::evaluate(&converted.recipe, &mut Decider::default()).unwrap();
    let a = built.bodies[&converted.steps["a"]];
    let b = built.bodies[&converted.steps["b"]];
    let mut model = built.model;
    let mut combined = model.body(a).clone();
    combined.lumps.extend_from_slice(&model.body(b).lumps);
    let body = model.add_body(combined);
    assert!(model.check(body).is_empty());
    (model, body)
}

#[test]
fn body_check_rejects_full_and_partial_opposite_interior_overlap() {
    for offset in [[0, 0, 1], [3, 1, 1]] {
        let (model, body) = two_boxes(offset);
        let (faults, work) = with_operation_budget(100_000_000, || {
            check_body(&model, body, &mut Decider::default()).unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        assert!(
            faults
                .iter()
                .any(|fault| fault.code == "topology.coincident-opposite-faces"),
            "{offset:?}: {faults:?}"
        );
    }
}

#[test]
fn body_check_permits_boundary_contact_and_disjoint_coplanar_trims() {
    for offset in [[4, 0, 1], [5, 0, 1], [0, 0, 2]] {
        let (model, body) = two_boxes(offset);
        let (faults, work) = with_operation_budget(100_000_000, || {
            check_body(&model, body, &mut Decider::default()).unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        assert!(faults.is_empty(), "{offset:?}: {faults:?}");
    }
}

#[test]
fn body_check_handles_curved_carrier_trims_and_their_boundary_contacts() {
    for (height, overlaps) in [(0, true), (1, false)] {
        let source = json!({
            "format":"soliddraft-recipe/1", "name":"cylindrical sheets", "result":"a",
            "steps":[
                {"id":"p","op":"profile","plane":"XY","circle":{"center":[0,0],"radius":2}},
                {"id":"a","op":"extrude","profile":"p","distance":1},
                {"id":"q","op":"profile","plane":{"origin":[0,0,height],"normal":[0,0,1],"x_dir":[1,0,0]},"circle":{"center":[0,0],"radius":3}},
                {"id":"outer","op":"extrude","profile":"q","distance":1},
                {"id":"r","op":"profile","plane":{"origin":[0,0,height],"normal":[0,0,1],"x_dir":[1,0,0]},"circle":{"center":[0,0],"radius":2}},
                {"id":"inner","op":"extrude","profile":"r","distance":1},
                {"id":"b","op":"boolean","kind":"subtract","target":"outer","tool":"inner"}
            ]
        });
        let ((), work) = with_operation_budget(100_000_000, || {
            let converted = convert(&source.to_string()).unwrap();
            let built =
                catallax_features::evaluate(&converted.recipe, &mut Decider::default()).unwrap();
            let a = built.bodies[&converted.steps["a"]];
            let b = built.bodies[&converted.steps["b"]];
            let mut model = built.model;
            let mut combined = model.body(a).clone();
            combined.lumps.extend_from_slice(&model.body(b).lumps);
            let body = model.add_body(combined);
            assert!(model.check(body).is_empty());
            let faults = check_body(&model, body, &mut Decider::default()).unwrap();
            assert_eq!(
                faults
                    .iter()
                    .any(|fault| fault.code == "topology.coincident-opposite-faces"),
                overlaps,
                "{height}: {faults:?}"
            );
        });
        assert!(!work.exhausted, "{work:?}");
    }
}
