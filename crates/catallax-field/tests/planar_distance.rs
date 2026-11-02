//! Exact distances use independently reconstructed physical planar trim chains.

use catallax_field::{FaceState, Plan, Workspace};
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, distance::Options};
use serde_json::{Value, json};

fn pose() -> Value {
    json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]})
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}
fn build(steps: Vec<Value>) -> Plan {
    let result = steps.last().unwrap()["id"].clone();
    let source = Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"trimmed plane","steps":steps,"result":result})
            .to_string(),
    )
    .unwrap();
    let (plan, work) = with_operation_budget(100_000_000, || Plan::compile(&source).unwrap());
    assert!(!work.exhausted, "{work:?}");
    plan
}
fn check(plan: &Plan, at: RPoint3, squared: i64) {
    let cap = plan
        .faces()
        .iter()
        .find(|face| face.origin.step.as_str() == "stock" && face.origin.role == Role::End)
        .unwrap();
    let id = plan.select(&cap.origin).unwrap();
    let mut dec = Decider::default();
    let controls = Options {
        width: Rational::parse("1e-8").unwrap(),
        max_cells: 1,
    };
    let (answer, work) = with_operation_budget(100_000_000, || {
        plan.prepare_distances(&mut dec)
            .unwrap()
            .point_to_face(id, &at.to_real(), &controls, &mut dec)
            .unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    let expected = dec
        .enclose_at(&Real::from_i64(squared).sqrt(), 256)
        .unwrap();
    let upper = answer
        .upper
        .as_ref()
        .expect("complete physical trim witness");
    assert!(
        answer.lower <= expected.lo() && upper >= &expected.hi(),
        "{answer:?}"
    );
    assert!(
        answer.converged && upper - &answer.lower <= controls.width,
        "{answer:?}"
    );
    let witness = answer.witnesses.unwrap()[1]
        .as_rational()
        .expect("these independent oracles have rational closest points");
    assert_eq!(
        plan.closed_face_query_checked(id, &witness, &mut Workspace::default(), &mut dec)
            .unwrap()
            .state,
        FaceState::OnFace
    );
}

#[test]
fn reentrant_and_holed_profiles_have_complete_planar_distances() {
    let loops = |points: &[[i64; 2]], prefix: &str| {
        json!(points.iter().enumerate().map(|(i,_)| json!({"kind":"line","id":format!("{prefix}-{i}"),"to":points[(i+1)%points.len()].map(|n| n.to_string())})).collect::<Vec<_>>())
    };
    let l = build(vec![
        json!({"id":"profile","op":"sketch","plane":pose(),"loops":[loops(&[[0,0],[4,0],[4,1],[1,1],[1,4],[0,4]],"l")]}),
        json!({"id":"stock","op":"extrude","profile":"profile","distance":"1"}),
    ]);
    check(&l, p(3, 3, 3), 8);
    let hole = build(vec![
        json!({"id":"profile","op":"sketch","plane":pose(),"loops":[loops(&[[0,0],[6,0],[6,6],[0,6]],"outer"),loops(&[[1,1],[1,3],[3,3],[3,1]],"hole")]}),
        json!({"id":"stock","op":"extrude","profile":"profile","distance":"1"}),
    ]);
    check(&hole, p(2, 2, 2), 2);
}

#[test]
fn circular_arcs_and_boolean_cut_rims_bound_the_actual_named_cap() {
    let mut steps = vec![
        json!({"id":"profile","op":"sketch","plane":pose(),"loops":[[
            {"kind":"arc","id":"arc","to":["-2","0"],"centre":["0","0"],"turn":"ccw"},
            {"kind":"line","id":"chord","to":["2","0"]}
        ]]}),
        json!({"id":"stock","op":"extrude","profile":"profile","distance":"2"}),
    ];
    let cap = build(steps.clone());
    check(&cap, p(0, -1, 5), 10);
    steps.extend([
        json!({"id":"tool","op":"block","size":["3","6","4"],"pose":{"rotation":["1","0","0","0"],"translation":["-3","-3","-1"]}}),
        json!({"id":"cut","op":"boolean","kind":"subtract","target":"stock","tool":"tool"}),
    ]);
    let cut = build(steps);
    check(&cut, p(-1, 1, 3), 2);
    check(&cut, p(-3, 0, 2), 9);
}

#[test]
fn named_coplanar_face_pairs_measure_the_finite_gap_between_disjoint_bodies() {
    use catallax_shapes::Surface;
    let plan = build(vec![
        json!({"id":"a","op":"block","size":["1","1","1"],"pose":pose()}),
        json!({"id":"b","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["3","3","0"]}}),
        json!({"id":"union","op":"boolean","kind":"union","target":"a","tool":"b"}),
    ]);
    let top = |name: &str| {
        let face = plan
            .faces()
            .iter()
            .find(|face| {
                face.origin.step.as_str() == name
                    && matches!(&face.surface,Surface::Plane(plane)
            if plane.origin.z.rational_ref()==Some(&Rational::one())
            && plane.normal.x.rational_ref()==Some(&Rational::zero())
            && plane.normal.y.rational_ref()==Some(&Rational::zero()))
            })
            .unwrap();
        plan.select(&face.origin).unwrap()
    };
    let ids = [top("a"), top("b")];
    let mut dec = Decider::default();
    let controls = Options {
        width: Rational::parse("1e-8").unwrap(),
        max_cells: 1,
    };
    let (answer, work) = with_operation_budget(100_000_000, || {
        plan.prepare_distances(&mut dec)
            .unwrap()
            .face_to_face(ids[0], ids[1], &controls, &mut dec)
            .unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    let oracle = dec.enclose_at(&Real::from_i64(8).sqrt(), 256).unwrap();
    let upper = answer.upper.as_ref().unwrap();
    assert!(
        answer.lower <= oracle.lo() && upper >= &oracle.hi() && answer.converged,
        "{answer:?}"
    );
    assert!(upper - &answer.lower <= controls.width, "{answer:?}");
    for (id, witness) in ids.into_iter().zip(answer.witnesses.unwrap()) {
        assert_eq!(
            plan.closed_face_query_checked(
                id,
                &witness.as_rational().unwrap(),
                &mut Workspace::default(),
                &mut dec
            )
            .unwrap()
            .state,
            FaceState::OnFace
        );
    }
}
