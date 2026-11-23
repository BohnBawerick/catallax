//! Graded wall geometry, face ownership, edits and verifier replay.

use catallax_app::{DEFAULT_BUILD_WORK, TARGET_BUILD_WORK};
use catallax_boolean::{Containment, classify_point};
use catallax_cli::verification::{ExportedQuestion, Request, answer};
use catallax_features::{BuildError, Built};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use catallax_topology::BodyKind;
use serde_json::{Value, json};
use std::fs;

fn pose(translation: [&str; 3]) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":translation})
}

fn stock() -> Value {
    json!({"id":"stock","op":"block","size":["10","8","6"],"pose":pose(["0","0","0"])})
}

fn grade(from: &str, to: &str, end: &str) -> Value {
    json!({"linear":{"from":[from,"0","0"],"to":[to,"0","0"],"start":"1","end":end}})
}

fn thicken(side: &str, thickness: Value) -> Value {
    json!({"id":"wall","op":"thicken_graded","body":"stock","side":side,
        "faces":[{"face":{"step":"stock","role":"end"}}],"thickness":thickness})
}

fn shell(open: Value, thickness: Value) -> Value {
    json!({"id":"wall","op":"shell_graded","body":"stock","open":open,"thickness":thickness})
}

fn recipe(steps: Vec<Value>, result: &str) -> Recipe {
    Recipe::parse(
        &json!({
            "format":"catallax-recipe/5","name":"graded wall","steps":steps,"result":result
        })
        .to_string(),
    )
    .unwrap()
}

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn volume(built: &Built, expected: &str) {
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert!(
        mass.volume.contains(&q(expected)),
        "volume {expected} is outside {:?}",
        mass.volume
    );
    assert!(mass.volume.width() <= q("1e-12"), "{:?}", mass.volume);
}

fn build(recipe: &Recipe) -> Built {
    let (built, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        catallax_features::evaluate(recipe, &mut Decider::default()).unwrap()
    });
    assert!(
        !work.exhausted && work.used <= TARGET_BUILD_WORK,
        "{work:?}"
    );
    assert!(built.model.check(built.result).is_empty());
    built
}

fn probes(recipe: &Recipe, built: &mut Built, points: &[([&str; 3], Location)]) {
    let plan = Plan::compile(recipe).unwrap();
    let mut dec = Decider::default();
    let mut workspace = Workspace::default();
    for &(point, expected) in points {
        let [x, y, z] = point.map(q);
        let query = RPoint3::new(x, y, z);
        let field = plan
            .proven_checked(&query, &mut workspace, &mut dec)
            .unwrap();
        assert_eq!(field.location, expected, "field point {point:?}");
        let actual = classify_point(&mut built.model, built.result, &query, &mut dec);
        let matches = match expected {
            Location::Inside => matches!(actual, Outcome::Proven(Containment::Inside)),
            Location::Outside => matches!(actual, Outcome::Proven(Containment::Outside)),
            Location::Boundary => matches!(actual, Outcome::ProvenCoincident(_)),
            _ => panic!("expected a settled location"),
        };
        assert!(
            matches,
            "B-rep point {point:?} should be {expected:?}: {actual:?}"
        );
    }
}

#[test]
fn thicken_preserves_exact_end_thicknesses_and_editable_face_names() {
    let front = recipe(
        vec![stock(), thicken("front", grade("0", "10", "3"))],
        "wall",
    );
    let (application, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        catallax_app::build_incremental_without_display(&front, None, &mut Decider::default())
            .unwrap()
    });
    assert!(
        !work.exhausted && work.used <= TARGET_BUILD_WORK,
        "{work:?}"
    );
    volume(application.checkpoint.built(), "160");
    assert!(
        application
            .checkpoint
            .built()
            .model
            .faces_of_body(application.checkpoint.built().result)
            .into_iter()
            .flat_map(|f| application.checkpoint.built().model.face_origins(f))
            .any(|o| o.step.as_str() == "wall"
                && matches!(&o.role, Role::Offset {of} if of.role == Role::End))
    );

    let mut front_body = build(&front);
    probes(
        &front,
        &mut front_body,
        &[
            (["0", "4", "7"], Location::Boundary),
            (["0", "4", "71/10"], Location::Outside),
            (["10", "4", "9"], Location::Boundary),
            (["10", "4", "91/10"], Location::Outside),
            (["1", "4", "8"], Location::Outside),
            (["9", "4", "8"], Location::Inside),
        ],
    );

    let edited = recipe(
        vec![stock(), thicken("front", grade("0", "10", "5"))],
        "wall",
    );
    let edited = Recipe::parse(&edited.to_text()).unwrap();
    let incremental = catallax_app::build_incremental_without_display(
        &edited,
        Some(application.checkpoint),
        &mut Decider::default(),
    )
    .unwrap();
    volume(incremental.checkpoint.built(), "240");
    let mut rebuilt = build(&edited);
    probes(
        &edited,
        &mut rebuilt,
        &[
            (["0", "4", "7"], Location::Boundary),
            (["0", "4", "71/10"], Location::Outside),
            (["10", "4", "11"], Location::Boundary),
            (["10", "4", "111/10"], Location::Outside),
        ],
    );
}

#[test]
fn a_later_extrude_can_attach_to_the_graded_offset_face() {
    let supported = recipe(
        vec![
            stock(),
            thicken("front", grade("0", "10", "3")),
            json!({"id":"section","op":"sketch","plane":pose(["0","0","0"]),"loops":[[
                {"kind":"line","id":"s0","to":["4","3"]},
                {"kind":"line","id":"s1","to":["4","4"]},
                {"kind":"line","id":"s2","to":["3","4"]},
                {"kind":"line","id":"s3","to":["3","3"]}
            ]]}),
            json!({"id":"boss","op":"extrude_on_face","profile":"section","body":"wall",
            "face":{"step":"wall","role":{"offset":{"of":{"step":"stock","role":"end"}}}},
            "normal_direction":["-1","0","5"],"x_direction":["5","0","1"],
            "offset":"0","distance":"1"}),
        ],
        "boss",
    );
    let reopened = Recipe::parse(&supported.to_text()).unwrap();
    let built = build(&reopened);
    volume(&built, "1");
    Plan::compile(&reopened).unwrap();
}

#[test]
fn reversed_anchors_and_internal_clamps_work_on_the_back() {
    let thickness =
        json!({"linear":{"from":["8","0","0"],"to":["2","0","0"],"start":"3","end":"1"}});
    let back = recipe(vec![stock(), thicken("back", thickness)], "wall");
    let mut built = build(&back);
    volume(&built, "160");
    probes(
        &back,
        &mut built,
        &[
            (["1", "4", "5"], Location::Boundary),
            (["2", "4", "11/2"], Location::Inside),
            (["5", "4", "4"], Location::Boundary),
            (["8", "4", "11/2"], Location::Inside),
            (["9", "4", "3"], Location::Boundary),
            (["9", "4", "13/2"], Location::Outside),
        ],
    );
}

#[test]
fn open_and_closed_shell_have_exact_cavities_and_keep_outer_faces() {
    for (open, expected) in [
        (json!([{"near":["5","4","7"]}]), "9192/25"),
        (json!([]), "10344/25"),
    ] {
        let tray = recipe(vec![stock(), shell(open, grade("0", "10", "3"))], "wall");
        let mut built = build(&tray);
        volume(&built, expected);
        let roles: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f))
            .collect();
        assert!(
            roles
                .iter()
                .any(|o| o.step.as_str() == "stock" && o.role == Role::End)
        );
        assert!(roles.iter().any(|o| o.step.as_str() == "wall"
            && matches!(&o.role, Role::Inner {of} if of.role == Role::Start)));
        probes(
            &tray,
            &mut built,
            &[
                (["2", "9/5", "3"], Location::Outside),
                (["6", "9/5", "3"], Location::Inside),
                (["6", "11/5", "3"], Location::Boundary),
                (["6", "5/2", "3"], Location::Outside),
                (
                    ["5", "4", "6"],
                    if expected == "9192/25" {
                        Location::Outside
                    } else {
                        Location::Boundary
                    },
                ),
            ],
        );
    }
    let empty = recipe(
        vec![stock(), shell(json!([{}]), grade("0", "10", "3"))],
        "wall",
    );
    let built = build(&empty);
    assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
}

#[test]
fn the_shell_cavity_stays_open_across_both_internal_clamps() {
    let tray = recipe(
        vec![
            stock(),
            shell(
                json!([{"face":{"step":"stock","role":"end"}}]),
                grade("3", "6", "3"),
            ),
        ],
        "wall",
    );
    let mut built = build(&tray);
    volume(&built, "364");
    probes(
        &tray,
        &mut built,
        &[
            (["3", "4", "5"], Location::Outside),
            (["6", "4", "5"], Location::Outside),
            (["2", "3/2", "5"], Location::Outside),
            (["13/2", "3/2", "5"], Location::Inside),
        ],
    );
}

#[test]
fn the_saved_graded_wall_accepts_a_later_feature_and_replays_exact_step() {
    let cut = recipe(
        vec![
            stock(),
            shell(
                json!([{"face":{"step":"stock","role":"end"}}]),
                grade("0", "10", "3"),
            ),
            json!({"id":"tool","op":"block","size":["1","4","1"],"pose":pose(["3","-1","3"])}),
            json!({"id":"cut","op":"boolean","kind":"subtract","target":"wall","tool":"tool"}),
        ],
        "cut",
    );
    let cut = Recipe::parse(&cut.to_text()).unwrap();
    let mut built = build(&cut);
    volume(&built, "18299/50");
    probes(
        &cut,
        &mut built,
        &[
            (["7/2", "1", "7/2"], Location::Outside),
            (["7/2", "1", "5"], Location::Inside),
        ],
    );
    let out = std::env::temp_dir().join(format!("catallax-graded-wall-{}", std::process::id()));
    let request = Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem: serde_json::to_value(&cut).unwrap(),
        out: out.join("original"),
        max_work: DEFAULT_BUILD_WORK,
    };
    let original = answer(&request).unwrap();
    assert_eq!(original["manifest"]["status"], "answered", "{original}");
    assert_eq!(original["manifest"]["exact"], true);
    assert_eq!(original["work"]["over_target"], false, "{original}");
    let exported: ExportedQuestion =
        serde_json::from_slice(&fs::read(request.out.join("question.json")).unwrap()).unwrap();
    let request_again = exported.request(out.join("replay")).unwrap();
    let again = answer(&request_again).unwrap();
    assert_eq!(again["manifest"]["status"], "answered", "{again}");
    assert_eq!(original["work"], again["work"]);
    let step = fs::read(request.out.join("answer.step")).unwrap();
    assert_eq!(
        step,
        fs::read(request_again.out.join("answer.step")).unwrap()
    );
    fs::remove_dir_all(out).unwrap();
}

#[test]
fn uncertified_selection_grade_and_source_have_named_refusals() {
    let mut multi = thicken("front", grade("0", "10", "3"));
    multi["faces"] = json!([{}]);
    let mut diagonal = grade("0", "10", "3");
    diagonal["linear"]["to"][1] = json!("4");
    let cylinder =
        json!({"id":"stock","op":"cylinder","radius":"4","height":"6","pose":pose(["0","0","0"])});
    for (source, op, code) in [
        (stock(), multi, "distance.grade-selection-unsupported"),
        (
            stock(),
            thicken("front", diagonal),
            "distance.grade-function-unsupported",
        ),
        (
            stock(),
            thicken("front", grade("0", "0", "3")),
            "distance.grade-coincident-anchors",
        ),
        (
            stock(),
            thicken("front", grade("0", "10", "0")),
            "distance.non-positive-thickness",
        ),
        (
            stock(),
            shell(json!([]), grade("0", "10", "5")),
            "distance.grade-clearance",
        ),
        (
            cylinder,
            shell(json!([]), grade("0", "10", "3")),
            "distance.grade-source-unsupported",
        ),
    ] {
        let recipe = recipe(vec![source, op], "wall");
        let error = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap_err();
        assert!(
            matches!(error, BuildError::Geometry { ref refusal, .. } if refusal.code == code),
            "{error}"
        );
        assert_eq!(Plan::compile(&recipe).unwrap_err().code, code);
    }
}
