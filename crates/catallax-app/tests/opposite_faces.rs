//! Explicit geometric validation distinguishes opposite sheets from boundary contacts.
use catallax_app::soliddraft::convert;
use catallax_boolean::check_body;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model, block};
use serde_json::json;

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

#[test]
fn thousands_of_separated_lumps_validate_within_the_default_fixed_allowance() {
    const LUMPS: usize = 2_358;
    let mut model = Model::new();
    let first = block(
        &mut model,
        [4.into(), 4.into(), 1.into()],
        &Pose::identity(),
    )
    .unwrap();
    let mut combined = model.body(first).clone();
    combined.lumps.reserve(LUMPS - 1);
    let mut dec = Decider::default();
    for index in 1..LUMPS {
        let pose = Pose::new(
            Rotation::identity(),
            RVector3::new(
                Rational::zero(),
                Rational::zero(),
                Rational::from_i64(10 * index as i64),
            )
            .into(),
            &mut dec,
        )
        .unwrap();
        let body = block(&mut model, [4.into(), 4.into(), 1.into()], &pose).unwrap();
        combined.lumps.extend_from_slice(&model.body(body).lumps);
    }
    let body = model.add_body(combined);
    assert_eq!(model.faces_of_body(body).len(), 6 * LUMPS);
    let (faults, work) = with_operation_budget(100_000_000, || {
        check_body(&model, body, &mut Decider::default()).unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    assert!(faults.is_empty(), "{faults:?}");
    eprintln!(
        "validated {LUMPS} separated lumps / {} faces using {} work",
        6 * LUMPS,
        work.used
    );
}
