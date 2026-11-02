//! Finite kept-face tubes, meeting cavities and names on nonconvex planar bodies.

use catallax_features::{Built, evaluate};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

const WORK: u64 = 100_000_000;
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn point(xyz: [&str; 3]) -> RPoint3 {
    let [x, y, z] = xyz.map(q);
    RPoint3::new(x, y, z)
}
fn document() -> Value {
    serde_json::from_str(include_str!("fixtures/nonconvex-shell.json")).unwrap()
}
fn parse(source: &Value) -> Recipe {
    Recipe::parse(&source.to_string()).unwrap()
}
fn build(source: &Value) -> (Built, Plan) {
    let recipe = parse(source);
    let (built, work) =
        with_operation_budget(WORK, || evaluate(&recipe, &mut Decider::default()).unwrap());
    assert!(!work.exhausted, "{work:?}");
    assert!(built.model.check(built.result).is_empty());
    let (plan, work) = with_operation_budget(WORK, || Plan::compile(&recipe).unwrap());
    assert!(!work.exhausted, "{work:?}");
    (built, plan)
}
fn material(built: &mut Built, plan: &Plan, cases: &[([&str; 3], bool)]) {
    for (xyz, inside) in cases {
        let p = point(*xyz);
        let expected = if *inside {
            catallax_boolean::Containment::Inside
        } else {
            catallax_boolean::Containment::Outside
        };
        let actual = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &p,
            &mut Decider::default(),
        );
        assert!(
            matches!(actual, Outcome::Proven(value) if value == expected),
            "B-rep {xyz:?}: {actual:?}"
        );
        let query = plan.proven(&p, &mut Workspace::default());
        assert_eq!(
            query.location,
            if *inside {
                Location::Inside
            } else {
                Location::Outside
            },
            "field {xyz:?}: {query:?}"
        );
    }
}
fn volume(built: &Built, constant: &str, pi: &str) {
    let (answer, work) = with_operation_budget(WORK, || {
        catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    let oracle =
        Ball::from_rational(&q(constant), 256) + Ball::pi(256) * Ball::from_rational(&q(pi), 256);
    assert!(
        answer.value.volume.lo <= oracle.lo() && answer.value.volume.hi >= oracle.hi(),
        "{:?}",
        answer.value.volume
    );
    assert!(answer.value.volume.width() < q("1e-6"));
}
fn loops(points: &[[i64; 2]], prefix: &str) -> Value {
    json!(points.iter().enumerate().map(|(i,_)| json!({"kind":"line","id":format!("{prefix}-{i}"),"to":points[(i+1)%points.len()].map(|n| n.to_string())})).collect::<Vec<_>>())
}
fn keep_only(source: &mut Value, sides: &[&str], kept: &str) {
    let mut open = vec![
        json!({"face":{"step":"stock","role":"start"}}),
        json!({"face":{"step":"stock","role":"end"}}),
    ];
    open.extend(
        sides
            .iter()
            .filter(|side| **side != kept)
            .map(|side| json!({"face":{"step":"stock","role":{"side":{"segment":side}}}})),
    );
    source["steps"][2]["open"] = json!(open);
}

#[test]
fn reentrant_shell_uses_corner_distance_not_untrimmed_plane_clearances() {
    let (mut built, plan) = build(&document());
    material(
        &mut built,
        &plan,
        &[
            (["21/10", "21/10", "2"], false),
            (["5/2", "5/2", "2"], true),
            (["4", "3/2", "2"], false),
            (["4", "1/2", "2"], true),
            (["4", "3/2", "1/2"], true),
            (["4", "4", "2"], false),
        ],
    );
    volume(&built, "120", "3/4");
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .any(|face| built
                .model
                .face_origins(face)
                .iter()
                .any(|origin| matches!(origin.role, Role::Round { .. })))
    );
}

#[test]
fn a_kept_finite_face_rounds_its_open_rims_without_extending_the_plane() {
    let mut source = document();
    keep_only(
        &mut source,
        &[
            "bottom",
            "right",
            "notch-bottom",
            "notch-left",
            "top",
            "left",
        ],
        "notch-bottom",
    );
    let (mut built, plan) = build(&source);
    material(
        &mut built,
        &plan,
        &[
            (["5/2", "7/2", "2"], true),
            (["2", "7/2", "2"], false),
            (["4", "5/2", "2"], true),
            (["1", "5/2", "2"], false),
        ],
    );
    volume(&built, "20", "2");
}

#[test]
fn kept_face_normal_images_reach_material_on_either_side_of_the_source_plane() {
    let mut source = document();
    source["steps"][0]["loops"] = json!([loops(
        &[
            [0, 0],
            [8, 0],
            [8, 8],
            [5, 8],
            [5, 2],
            [3, 2],
            [3, 8],
            [0, 8]
        ],
        "u"
    )]);
    source["steps"][2]["thickness"] = json!("3");
    keep_only(
        &mut source,
        &["u-0", "u-1", "u-2", "u-3", "u-4", "u-5", "u-6", "u-7"],
        "u-5",
    );
    let (mut built, plan) = build(&source);
    material(
        &mut built,
        &plan,
        &[
            (["11/2", "5", "2"], true),
            (["13/2", "5", "2"], false),
            (["1", "5", "2"], true),
            (["4", "5", "2"], false),
        ],
    );
}

#[test]
fn hole_envelopes_meet_split_the_cavity_and_then_consume_it() {
    let mut source = document();
    source["steps"][0]["loops"] = json!([
        loops(&[[0, 0], [10, 0], [10, 10], [0, 10]], "outer"),
        loops(&[[3, 3], [3, 7], [7, 7], [7, 3]], "hole")
    ]);
    source["steps"][2]["open"] = json!([
        {"face":{"step":"stock","role":"start"}}, {"face":{"step":"stock","role":"end"}}
    ]);
    source["steps"][2]["thickness"] = json!("3/2");
    let (mut built, plan) = build(&source);
    material(
        &mut built,
        &plan,
        &[
            (["8/5", "8/5", "2"], false),
            (["7/5", "5", "2"], true),
            (["5", "5", "2"], false),
        ],
    );
    volume(&built, "300", "9");
    source["steps"][2]["thickness"] = json!("2");
    let (mut built, plan) = build(&source);
    material(
        &mut built,
        &plan,
        &[(["8/5", "8/5", "2"], true), (["5", "5", "2"], false)],
    );
    volume(&built, "336", "0");
    assert!(
        !built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|face| built.model.face_origins(face))
            .any(|origin| matches!(
                origin.role,
                Role::Inner { .. } | Role::Round { .. } | Role::RoundCorner { .. }
            ))
    );
}
