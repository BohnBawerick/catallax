//! The field's support-sign query is independent of height evaluation and B-rep.

use catallax_features::end_condition::curved::{AuthoredCap, Location};
use catallax_field::CurvedField;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::RPoint3;

fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn compile() -> (AuthoredCap, CurvedField) {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/curved-limits/sphere-source.json"
    ))
    .unwrap();
    let cap = AuthoredCap::compile(
        &recipe,
        &id("target"),
        &id("profile"),
        None,
        None,
        &mut Decider::default(),
    )
    .unwrap();
    let field = CurvedField::new(cap.clone(), id("capped"));
    (cap, field)
}

#[test]
fn support_signs_match_the_closed_first_hit_set_and_authored_names() {
    let (cap, field) = compile();
    let mut dec = Decider::default();
    for (p, expected) in [
        (["0", "0", "2"], Location::Inside),
        (["0", "0", "4"], Location::Outside),
        (["0", "0", "3"], Location::Boundary),
        (["1/2", "0", "2"], Location::Boundary),
        (["3/4", "0", "2"], Location::Outside),
        (["0", "0", "0"], Location::Boundary),
        (["0", "0", "7"], Location::Outside),
    ] {
        let [x, y, z] = p.map(q);
        let point = RPoint3::new(x.clone(), y.clone(), z.clone());
        assert_eq!(field.location(&point, &mut dec).unwrap(), expected);
        assert_eq!(
            cap.definition()
                .location(&[x.into(), y.into()], &z.into(), &mut dec)
                .unwrap(),
            expected
        );
    }
    let start = field
        .origins(&RPoint3::new(q("0"), q("0"), q("0")), &mut dec)
        .unwrap();
    assert_eq!(start.len(), 1);
    assert_eq!(start[0].role, Role::Start);
    let end = field
        .origins(&RPoint3::new(q("0"), q("0"), q("3")), &mut dec)
        .unwrap();
    assert_eq!(end.len(), 1);
    assert_eq!(end[0].role, Role::End);
    let side = field
        .origins(&RPoint3::new(q("1/2"), q("0"), q("2")), &mut dec)
        .unwrap();
    assert_eq!(side.len(), 1);
    assert!(matches!(&side[0].role,Role::Side { segment,.. } if segment.as_str()=="right"));
    assert!(
        start
            .iter()
            .chain(&end)
            .chain(&side)
            .all(|o| o.step == id("capped"))
    );
}

#[test]
fn query_cannot_reuse_a_definition_to_escape_request_exhaustion() {
    let (_, field) = compile();
    let mut dec = Decider::default();
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        field.location(&RPoint3::new(q("0"), q("0"), q("2")), &mut dec)
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}

#[test]
fn inside_and_backward_queries_keep_the_first_crossing_and_world_origins() {
    for (start, direction, inside, outside, end) in [
        ("5", None, "6", "8", "7"),
        ("10", Some(["0", "0", "-2"]), "8", "6", "7"),
    ] {
        let mut value: serde_json::Value = serde_json::from_str(include_str!(
            "../../catallax-features/tests/fixtures/curved-limits/sphere-source.json"
        ))
        .unwrap();
        value["steps"][0]["plane"]["translation"][2] = start.into();
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let direction = direction.map(|v| {
            serde_json::from_value::<catallax_recipe::Vector3>(serde_json::json!(v)).unwrap()
        });
        let mut dec = Decider::default();
        let cap = AuthoredCap::compile(
            &recipe,
            &id("target"),
            &id("profile"),
            direction.as_ref(),
            None,
            &mut dec,
        )
        .unwrap();
        let field = CurvedField::new(cap, id("capped"));
        for (z, want) in [
            (inside, Location::Inside),
            (outside, Location::Outside),
            (end, Location::Boundary),
        ] {
            assert_eq!(
                field
                    .location(&RPoint3::new(q("0"), q("0"), q(z)), &mut dec)
                    .unwrap(),
                want
            );
        }
        let origins = field
            .origins(&RPoint3::new(q("0"), q("0"), q(end)), &mut dec)
            .unwrap();
        assert_eq!(origins.len(), 1);
        assert_eq!(origins[0].role, Role::End);
        assert_eq!(origins[0].step, id("capped"));
    }
}

#[test]
fn native_plan_dispatch_classifies_curved_entry_exit_and_selected_end_faces() {
    use catallax_field::{FaceState, Location as PlanLocation, Plan, Workspace};
    for (text, inside, outside, end) in [
        (
            include_str!(
                "../../catallax-features/tests/fixtures/curved-limits/sphere-entry-to.json"
            ),
            ["0", "0", "2"],
            ["0", "0", "4"],
            ["0", "0", "3"],
        ),
        (
            include_str!(
                "../../catallax-features/tests/fixtures/curved-limits/sphere-exit-next.json"
            ),
            ["0", "0", "6"],
            ["0", "0", "8"],
            ["0", "0", "7"],
        ),
        (
            include_str!(
                "../../catallax-features/tests/fixtures/curved-limits/sphere-backward-next.json"
            ),
            ["0", "0", "8"],
            ["0", "0", "6"],
            ["0", "0", "7"],
        ),
        (
            include_str!(
                "../../catallax-features/tests/fixtures/curved-limits/cylinder-entry-to.json"
            ),
            ["1", "0", "5"],
            ["4", "0", "5"],
            ["3", "0", "5"],
        ),
        (
            include_str!(
                "../../catallax-features/tests/fixtures/curved-limits/cone-entry-next.json"
            ),
            ["1", "0", "5"],
            ["4", "0", "5"],
            ["2", "0", "5"],
        ),
    ] {
        let recipe = Recipe::parse(text).unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        let point = |p: [&str; 3]| {
            let [x, y, z] = p.map(q);
            RPoint3::new(x, y, z)
        };
        for (p, want) in [
            (inside, PlanLocation::Inside),
            (outside, PlanLocation::Outside),
            (end, PlanLocation::Boundary),
        ] {
            assert_eq!(
                plan.proven(&point(p), &mut workspace).location,
                want,
                "{} at {p:?}",
                recipe.name
            );
        }
        let end_origin = catallax_recipe::Origin {
            step: id("cap"),
            role: Role::End,
        };
        // Several construction supports carry End. Native handles retain those
        // supports; face_query deliberately resolves their named union.
        let faces: Vec<_> = plan
            .faces()
            .iter()
            .enumerate()
            .filter(|(_, face)| face.origin == end_origin)
            .filter_map(|(index, _)| plan.result_face(index))
            .collect();
        assert!(!faces.is_empty());
        let endpoint = point(end).to_real();
        let mut incident_supports = 0;
        let mut dec = Decider::default();
        for face in faces {
            assert_eq!(
                plan.face_query(face, &point(end), &mut workspace)
                    .unwrap()
                    .state,
                FaceState::OnFace
            );
            let residual = plan.faces()[face.index()]
                .surface
                .implicit()
                .expect("analytic cap support")
                .eval_real(&[endpoint.x.clone(), endpoint.y.clone(), endpoint.z.clone()]);
            match dec.sign(&residual) {
                catallax_numbers::Outcome::ProvenCoincident(_) => incident_supports += 1,
                catallax_numbers::Outcome::Proven(_) => (),
                other => panic!("uncertified end-support incidence: {other:?}"),
            }
            assert_eq!(
                plan.face_query(face, &point(inside), &mut workspace)
                    .unwrap()
                    .state,
                FaceState::NotOnFace
            );
        }
        assert_eq!(
            incident_supports, 1,
            "{} must retain the unique physical end",
            recipe.name
        );
    }
}

#[test]
fn native_field_display_keeps_its_default_certificate() {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/curved-limits/sphere-entry-next.json"
    ))
    .unwrap();
    let plan = catallax_field::Plan::compile(&recipe).unwrap();
    let options = catallax_field::display::Options::default();
    let (result, work) = catallax_numbers::with_operation_budget(100_000_000, || {
        catallax_field::display::mesh_plan(&plan, &options, &mut Decider::default())
    });
    let mesh = result.unwrap();
    assert!(!work.exhausted);
    assert!(!mesh.triangles.is_empty());
    assert!(mesh.max_deviation <= options.tolerance);
    let mut roles = Vec::new();
    for triangle in &mesh.triangles {
        let origin = &plan.faces()[triangle.face.index()].origin;
        assert_eq!(origin.step, id("cap"));
        if !roles.contains(&origin.role) {
            roles.push(origin.role.clone());
        }
    }
    for role in [Role::Start, Role::End] {
        assert!(roles.contains(&role), "missing retained cap {role:?}");
    }
    for side in ["left", "right", "bottom", "top"] {
        assert!(
            roles
                .iter()
                .any(|role| matches!(role, Role::Side {segment,..} if segment.as_str()==side)),
            "missing crossing-boundary wall {side}"
        );
    }
    // The cells cross at the trim edges. The display optimization must not
    // classify the retained interior or closed side boundary as outside.
    let mut workspace = catallax_field::Workspace::default();
    for (p, want) in [
        (["0", "0", "2"], catallax_field::Location::Inside),
        (["0", "0", "3"], catallax_field::Location::Boundary),
        (["1/2", "0", "2"], catallax_field::Location::Boundary),
        (["1/2", "1/2", "0"], catallax_field::Location::Boundary),
        (["0", "0", "4"], catallax_field::Location::Outside),
    ] {
        let [x, y, z] = p.map(q);
        assert_eq!(
            plan.proven(&RPoint3::new(x, y, z), &mut workspace).location,
            want,
            "{p:?}"
        );
    }
}
