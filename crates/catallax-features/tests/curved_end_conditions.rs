//! Independent whole-domain certificates and real Boolean boundary assembly.
//! Tests consume the production feature module on the merged dependency base.

use catallax_features::end_condition::curved;

use catallax_numbers::{Decider, Outcome, Rational, Real};
use catallax_recipe::{ElementId, Recipe, Role, StepId};
use catallax_shapes::{Point3, RPoint3, Vector3};
use catallax_topology::Model;
use curved::{AuthoredCap, Branch, CurvedCapped, Frame, Location, Scene};

fn q(x: &str) -> Rational {
    Rational::parse(x).unwrap()
}
fn x(x: &str) -> Real {
    q(x).into()
}
fn vector(a: [&str; 3]) -> Vector3 {
    let [x, y, z] = a.map(x);
    Vector3::new(x, y, z)
}
fn frame(start: &str, direction: &str) -> Frame {
    Frame {
        origin: Point3::new(x("0"), x("0"), x(start)),
        u: vector(["1", "0", "0"]),
        v: vector(["0", "1", "0"]),
        direction: vector(["0", "0", direction]),
    }
}
fn scene(dec: &mut Decider) -> Scene {
    let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"sphere source","steps":[
        {"id":"target","op":"sphere","radius":"2","pose":{"rotation":["1","0","0","0"],"translation":["0","0","5"]}}
    ],"result":"target"}"#).unwrap();
    Scene::compile(&recipe, &StepId::new("target").unwrap(), dec).unwrap()
}
fn cap(start: &str, direction: &str, dec: &mut Decider) -> CurvedCapped {
    CurvedCapped::rectangle(
        scene(dec),
        frame(start, direction),
        [[q("-1/2"), q("1/2")], [q("-1/2"), q("1/2")]],
        dec,
    )
    .unwrap()
}
fn ids() -> [ElementId; 4] {
    ["left", "right", "bottom", "top"].map(|s| ElementId::new(s).unwrap())
}
fn equal(a: &Real, b: &Real, dec: &mut Decider) {
    assert!(matches!(dec.sign(&(a - b)), Outcome::ProvenCoincident(_)));
}

#[test]
fn first_entries_exits_and_raw_direction_scaling_keep_exact_spherical_roots() {
    let mut dec = Decider::default();
    for (start, w, branch, center) in [
        ("0", "1", Branch::Entry, "3"),
        ("5", "1", Branch::Exit, "2"),
        ("10", "-2", Branch::Entry, "3/2"),
    ] {
        let c = cap(start, w, &mut dec);
        assert_eq!(c.branch(), branch);
        equal(
            &c.height(&[x("0"), x("0")], &mut dec).unwrap(),
            &x(center),
            &mut dec,
        );
        let t = c.height(&[x("1/2"), x("1/2")], &mut dec).unwrap();
        let point = c
            .frame()
            .point(&[x("1/2"), x("1/2")])
            .plus(&c.frame().direction.scaled(&t));
        let residual = c
            .scene()
            .surface()
            .implicit()
            .expect("stock quadric support")
            .eval_real(&[point.x, point.y, point.z]);
        equal(&residual, &Real::zero(), &mut dec);
    }
}

#[test]
fn field_coordinates_classify_actual_boundary_and_keep_closed_edges() {
    let mut dec = Decider::default();
    let c = cap("0", "1", &mut dec);
    for (uv, t, want) in [
        (["0", "0"], "2", Location::Inside),
        (["0", "0"], "4", Location::Outside),
        (["0", "0"], "3", Location::Boundary),
        (["1/2", "0"], "2", Location::Boundary),
        (["3/4", "0"], "2", Location::Outside),
        (["0", "0"], "-1", Location::Outside),
    ] {
        assert_eq!(c.location(&uv.map(x), &x(t), &mut dec).unwrap(), want);
    }
}

#[test]
fn silhouette_and_start_crossings_refuse_as_coverage_without_claiming_impossibility() {
    let mut dec = Decider::default();
    let error = CurvedCapped::rectangle(
        scene(&mut dec),
        frame("0", "1"),
        [[q("-3"), q("3")], [q("-1/2"), q("1/2")]],
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(error.code, "features.unsupported-curved-limit");
    let error = CurvedCapped::rectangle(
        scene(&mut dec),
        frame("3", "1"),
        [[q("-1/2"), q("1/2")], [q("-1/2"), q("1/2")]],
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(error.code, "features.unsupported-curved-limit");
}

#[test]
fn analytic_boundary_build_keeps_only_extrusion_origins_and_matches_point_set() {
    for (start, inside_z, outside_z) in [("0", "2", "4"), ("5", "6", "8")] {
        let mut dec = Decider::default();
        let c = cap(start, "1", &mut dec);
        let mut model = Model::new();
        let step = StepId::new("capped").unwrap();
        let body = c.build(&mut model, &step, ids(), &mut dec).unwrap();
        assert!(model.check(body).is_empty());
        let mut ends = 0;
        for face in model.faces_of_body(body) {
            let names = model.face_origins(face);
            assert!(!names.is_empty());
            assert!(names.iter().all(|o| o.step == step));
            ends += usize::from(names.iter().any(|o| o.role == Role::End));
        }
        assert_eq!(ends, 1);
        for (z, want) in [
            (inside_z, catallax_boolean::Containment::Inside),
            (outside_z, catallax_boolean::Containment::Outside),
        ] {
            let p = RPoint3::new(q("0"), q("0"), q(z));
            assert!(
                matches!(catallax_boolean::classify_point(&mut model,body,&p,&mut dec),Outcome::Proven(found) if found==want)
            );
        }
    }
}

#[test]
fn irrational_branch_is_above_the_axis_hit_without_rounding() {
    let mut dec = Decider::default();
    let c = cap("0", "1", &mut dec);
    let t = c.height(&[x("1/2"), x("0")], &mut dec).unwrap();
    assert!(matches!(
        dec.sign(&(t - x("3"))),
        Outcome::Proven(catallax_numbers::StrictSign::Positive)
    ));
}

#[test]
fn complete_curved_build_and_step_bytes_keep_work_at_one_two_and_four_workers() {
    fn artifact(case: usize) -> (Vec<u8>, catallax_topology::EntityCounts, u64) {
        let (result, work) = catallax_numbers::with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let (start, w) = [("0", "1"), ("5", "1"), ("10", "-2")][case];
            let c = cap(start, w, &mut dec);
            let mut model = Model::new();
            let body = c
                .build(&mut model, &StepId::new("capped").unwrap(), ids(), &mut dec)
                .unwrap();
            let mut bytes = Vec::new();
            let report = catallax_step::write_step(
                &model,
                &[body],
                &Default::default(),
                &mut dec,
                &mut bytes,
            )
            .unwrap();
            assert!(report.marks.is_empty());
            assert_eq!(report.value.solids, 1);
            (bytes, model.counts(body))
        });
        assert!(!work.exhausted);
        (result.0, result.1, work.used)
    }
    let mut baseline = None;
    for workers in [1, 2, 4] {
        let handles: Vec<_> = (0..workers)
            .map(|worker| {
                std::thread::spawn(move || {
                    (worker..3)
                        .step_by(workers)
                        .map(|case| {
                            // A fresh numeric thread isolates request-local caches just as
                            // the verification CLI does; workers only schedule requests.
                            let result = std::thread::spawn(move || artifact(case)).join().unwrap();
                            (case, result)
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        let mut actual: Vec<_> = handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect();
        actual.sort_by_key(|(case, _)| *case);
        if let Some(expected) = &baseline {
            assert_eq!(&actual, expected);
        } else {
            baseline = Some(actual);
        }
    }
}

#[test]
fn transverse_cylinder_and_cone_caps_use_the_closed_finite_wall() {
    for shape in ["cylinder", "cone"] {
        let dimensions = if shape == "cylinder" {
            r#""radius":"2""#
        } else {
            r#""base_radius":"4","top_radius":"2""#
        };
        let json = format!(
            r#"{{"format":"catallax-recipe/3","name":"transverse source","steps":[{{"id":"target","op":"{shape}",{dimensions},"height":"10","pose":{{"rotation":["1","0","0","0"],"translation":["5","0","0"]}}}}],"result":"target"}}"#
        );
        let recipe = Recipe::parse(&json).unwrap();
        let mut dec = Decider::default();
        let target = Scene::compile(&recipe, &StepId::new("target").unwrap(), &mut dec).unwrap();
        let f = Frame {
            origin: Point3::new(x("0"), x("0"), x("5")),
            u: vector(["0", "1", "0"]),
            v: vector(["0", "0", "1"]),
            direction: vector(["1", "0", "0"]),
        };
        let c = CurvedCapped::rectangle(
            target,
            f,
            [[q("-1/2"), q("1/2")], [q("-1/2"), q("1/2")]],
            &mut dec,
        )
        .unwrap();
        let want = if shape == "cylinder" { "3" } else { "2" };
        equal(
            &c.height(&[x("0"), x("0")], &mut dec).unwrap(),
            &x(want),
            &mut dec,
        );
        let mut model = Model::new();
        let body = c
            .build(&mut model, &StepId::new("capped").unwrap(), ids(), &mut dec)
            .unwrap();
        assert!(model.check(body).is_empty());
        let sample = RPoint3::new(q("1"), q("0"), q("5"));
        assert!(matches!(
            catallax_boolean::classify_point(&mut model, body, &sample, &mut dec),
            Outcome::Proven(catallax_boolean::Containment::Inside)
        ));
    }
}

#[test]
fn zero_request_allowance_cannot_publish_a_curved_definition_or_body() {
    let mut dec = Decider::default();
    let s = scene(&mut dec);
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        CurvedCapped::rectangle(
            s,
            frame("0", "1"),
            [[q("-1/2"), q("1/2")], [q("-1/2"), q("1/2")]],
            &mut dec,
        )
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}

#[test]
fn clockwise_authored_rectangle_keeps_physical_segment_ids() {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/curved-limits/sphere-source.json")).unwrap();
    value["steps"][0]["loops"][0] = serde_json::json!([
        {"kind":"line","id":"bottom","to":["-1/2","-1/2"]},
        {"kind":"line","id":"left","to":["-1/2","1/2"]},
        {"kind":"line","id":"top","to":["1/2","1/2"]},
        {"kind":"line","id":"right","to":["1/2","-1/2"]}
    ]);
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let cap = AuthoredCap::compile(
        &recipe,
        &StepId::new("target").unwrap(),
        &StepId::new("profile").unwrap(),
        None,
        None,
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(cap.side_ids(), &ids());
}

#[test]
fn proper_placement_transports_the_exact_branch_and_boundary() {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/curved-limits/sphere-source.json")).unwrap();
    value["steps"][0]["plane"] =
        serde_json::json!({"rotation":["1","1","0","0"],"translation":["7","-11","3"]});
    value["steps"][1]["pose"] =
        serde_json::json!({"rotation":["1","1","0","0"],"translation":["7","-16","3"]});
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let mut dec = Decider::default();
    let cap = AuthoredCap::compile(
        &recipe,
        &StepId::new("target").unwrap(),
        &StepId::new("profile").unwrap(),
        None,
        None,
        &mut dec,
    )
    .unwrap();
    equal(
        &cap.definition()
            .height(&[x("0"), x("0")], &mut dec)
            .unwrap(),
        &x("3"),
        &mut dec,
    );
    let mut model = Model::new();
    let body = cap
        .build(&mut model, &StepId::new("capped").unwrap(), &mut dec)
        .unwrap();
    assert!(model.check(body).is_empty());
    assert!(model.body(body).marks.is_empty());
    for (y, want) in [
        ("-13", catallax_boolean::Containment::Inside),
        ("-15", catallax_boolean::Containment::Outside),
    ] {
        let p = RPoint3::new(q("7"), q(y), q("3"));
        assert!(
            matches!(catallax_boolean::classify_point(&mut model,body,&p,&mut dec),Outcome::Proven(found) if found==want)
        );
    }
}

#[test]
fn unproved_trim_coverage_does_not_become_a_missing_cap_claim() {
    let recipe=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"trim event","steps":[
        {"id":"target","op":"cylinder","radius":"2","height":"10","pose":{"rotation":["1","0","0","0"],"translation":["5","0","0"]}}
    ],"result":"target"}"#).unwrap();
    let mut dec = Decider::default();
    let scene = Scene::compile(&recipe, &StepId::new("target").unwrap(), &mut dec).unwrap();
    let frame = Frame {
        origin: Point3::new(x("0"), x("0"), x("1/4")),
        u: vector(["0", "1", "0"]),
        v: vector(["0", "0", "1"]),
        direction: vector(["1", "0", "0"]),
    };
    let error = CurvedCapped::rectangle(
        scene,
        frame,
        [[q("-1/2"), q("1/2")], [q("-1/2"), q("1/2")]],
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(error.code, "features.unsupported-curved-limit");
}

#[test]
fn recipe_dispatch_builds_entry_exit_and_selected_stock_walls() {
    for text in [
        include_str!("fixtures/curved-limits/sphere-entry-next.json"),
        include_str!("fixtures/curved-limits/sphere-entry-to.json"),
        include_str!("fixtures/curved-limits/sphere-exit-next.json"),
        include_str!("fixtures/curved-limits/sphere-backward-next.json"),
        include_str!("fixtures/curved-limits/cylinder-entry-next.json"),
        include_str!("fixtures/curved-limits/cylinder-entry-to.json"),
        include_str!("fixtures/curved-limits/cone-entry-next.json"),
        include_str!("fixtures/curved-limits/cone-entry-to.json"),
    ] {
        let recipe = Recipe::parse(text).unwrap();
        assert_eq!(Recipe::parse(&recipe.to_text()).unwrap(), recipe);
        let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        assert!(built.model.body(built.result).marks.is_empty());
        let names: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f))
            .collect();
        assert!(names.iter().all(|o| o.step.as_str() == "cap"));
        for role in [Role::Start, Role::End] {
            assert!(names.iter().any(|o| o.role == role));
        }
        for side in ids() {
            assert!(
                names
                    .iter()
                    .any(|o| matches!(&o.role, Role::Side { segment, .. } if *segment == side))
            );
        }
    }
}
