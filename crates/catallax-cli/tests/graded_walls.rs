//! Graded wall geometry, face ownership, edits and verifier replay.

use catallax_app::{DEFAULT_BUILD_WORK, TARGET_BUILD_WORK};
use catallax_boolean::{Containment, classify_point};
use catallax_cli::verification::{ExportedQuestion, Request, answer};
use catallax_features::{BuildError, Built};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{
    Decider, Outcome, Rational, TestDeadline, TestDeadlinePoint, with_operation_budget,
};
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

fn enclosing() -> Value {
    json!({"id":"enclosing","op":"block","size":["14","12","12"],"pose":pose(["-2","-2","-2"])})
}

fn intersection(id: &str, target: &str, tool: &str) -> Value {
    json!({"id":id,"op":"boolean","kind":"intersect","target":target,"tool":tool})
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

fn open_sides() -> Vec<Value> {
    ["x_min", "x_max", "y_min", "y_max"]
        .map(|segment| json!({"face":{"step":"stock","role":{"side":{"segment":segment}}}}))
        .into()
}

fn point_grade() -> Value {
    json!({"point":{"centre":["5","4","-2"],"r0":"7","r1":"8","start":"1","end":"2"}})
}

fn radial_grade() -> Value {
    json!({"radial":{"origin":["0","4","6"],"direction":["3","0","0"],
        "r0":"0","r1":"8","start":"1","end":"2"}})
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

fn face_mount(
    mut steps: Vec<Value>,
    body: &str,
    role: Value,
    normal: [&str; 3],
    x: [&str; 3],
) -> Recipe {
    steps.extend([
        json!({"id":"section","op":"sketch","plane":pose(["0","0","0"]),"loops":[[
            {"kind":"line","id":"s0","to":["4","3"]},
            {"kind":"line","id":"s1","to":["4","4"]},
            {"kind":"line","id":"s2","to":["3","4"]},
            {"kind":"line","id":"s3","to":["3","3"]}
        ]]}),
        json!({"id":"boss","op":"extrude_on_face","profile":"section","body":body,
            "face":{"step":"wall","role":role},"normal_direction":normal,"x_direction":x,
            "offset":"0","distance":"1"}),
    ]);
    let mounted = recipe(steps, "boss");
    Recipe::parse(&mounted.to_text()).unwrap()
}

fn offset_role() -> Value {
    json!({"offset":{"of":{"step":"stock","role":"end"}}})
}

#[test]
fn a_later_extrude_can_attach_to_the_graded_offset_face() {
    for body in ["wall", "trimmed"] {
        let supported = face_mount(
            vec![
                stock(),
                thicken("front", grade("0", "10", "3")),
                enclosing(),
                intersection("trimmed", "wall", "enclosing"),
            ],
            body,
            offset_role(),
            ["-1", "0", "5"],
            ["5", "0", "1"],
        );
        let built = build(&supported);
        volume(&built, "1");
        Plan::compile(&supported).unwrap();
    }
}

#[test]
fn a_later_extrude_can_attach_to_a_graded_rim_or_inner_wall() {
    let rim = json!({"rim":{"edge":{"between":[
        {"step":"stock","role":"end"},
        {"step":"stock","role":{"side":{"segment":"x_min"}}}
    ]}}});
    let inner = json!({"inner":{"of":{"step":"stock","role":"start"}}});
    for (wall, role, normal, x) in [
        (
            thicken("front", grade("2", "8", "3")),
            rim,
            ["-1", "0", "0"],
            ["0", "1", "0"],
        ),
        (
            shell(
                json!([{"face":{"step":"stock","role":"end"}}]),
                grade("0", "10", "3"),
            ),
            inner,
            ["-1", "0", "5"],
            ["5", "0", "1"],
        ),
    ] {
        let mounted = face_mount(vec![stock(), wall], "wall", role, normal, x);
        let built = build(&mounted);
        volume(&built, "1");
        Plan::compile(&mounted).unwrap();
    }
}

#[test]
fn a_coincident_replacement_does_not_restore_a_removed_graded_offset_name() {
    let wall = thicken("front", grade("0", "10", "3"));
    let mut replacement = wall.clone();
    replacement["id"] = json!("replacement");
    let mounted = face_mount(
        vec![
            stock(),
            wall,
            json!({"id":"low","op":"block","size":["14","12","17/2"],"pose":pose(["-2","-2","-2"])}),
            intersection("trimmed", "wall", "low"),
            replacement,
            json!({"id":"joined","op":"boolean","kind":"union","target":"trimmed","tool":"replacement"}),
        ],
        "joined",
        offset_role(),
        ["-1", "0", "5"],
        ["5", "0", "1"],
    );
    let error = catallax_features::evaluate(&mounted, &mut Decider::default()).unwrap_err();
    assert!(
        matches!(error, BuildError::Geometry { ref refusal, .. } if refusal.code == "reference.lost-face"),
        "{error}",
    );
    assert_eq!(
        Plan::compile(&mounted).unwrap_err().code,
        "reference.lost-face"
    );
}

#[test]
fn enclosing_intersection_preserves_the_six_graded_shell_walls() {
    let mut wall = shell(json!([]), grade("0", "10", "3"));
    wall["body"] = json!("trimmed");
    let tray = recipe(
        vec![
            stock(),
            enclosing(),
            intersection("trimmed", "stock", "enclosing"),
            wall,
        ],
        "wall",
    );
    let mut built = build(&tray);
    volume(&built, "10344/25");
    probes(
        &tray,
        &mut built,
        &[
            (["6", "9/5", "3"], Location::Inside),
            (["6", "5/2", "3"], Location::Outside),
        ],
    );
}

#[test]
fn graded_walls_compile_a_shared_intersection_dag_within_linear_visit_budget() {
    let mut second = stock();
    second["id"] = json!("second");
    second["size"] = json!(["11", "8", "6"]);
    let mut steps = vec![stock(), second];
    let (mut a, mut b) = ("stock".to_string(), "second".to_string());
    for level in 0..18 {
        let next_a = format!("a{level}");
        let next_b = format!("b{level}");
        steps.push(intersection(&next_a, &a, &b));
        steps.push(intersection(&next_b, &b, &a));
        (a, b) = (next_a, next_b);
    }
    for mut wall in [
        shell(json!([]), grade("0", "10", "3")),
        thicken("front", grade("0", "10", "3")),
    ] {
        wall["body"] = json!(a);
        let mut steps = steps.clone();
        steps.push(wall);
        let recipe = recipe(steps, "wall");
        let deadline = TestDeadline::after(TestDeadlinePoint::Enclosure, 100_000);
        let mut dec = Decider::default();
        dec.set_deadline(Some(deadline.instant()));
        let (plan, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
            Plan::compile_checked(&recipe, &mut dec)
        });
        assert!(!plan.unwrap().is_proven_empty());
        assert!(
            !work.exhausted && work.used <= TARGET_BUILD_WORK,
            "{work:?}"
        );
    }
}

#[test]
fn irrational_radial_wall_planes_keep_comparison_coverage_available() {
    let grade = json!({"radial":{"origin":["0","0","6"],"direction":["1","1","0"],
        "r0":"0","r1":"8","start":"1","end":"2"}});
    let mut open = open_sides();
    open.push(json!({"face":{"step":"stock","role":"start"}}));
    for wall in [thicken("front", grade.clone()), shell(json!(open), grade)] {
        let recipe = recipe(vec![stock(), wall], "wall");
        let controller = catallax_compare::Controller::default();
        assert_eq!(
            controller.begin().check(&recipe, &Default::default()),
            catallax_compare::Completion::Published
        );
        let report = controller.latest().unwrap();
        assert!(report.refusal.is_none(), "{report:?}");
        assert!(report.reverse_samples > 0, "{report:?}");
        assert_ne!(
            report.verdict,
            catallax_compare::Verdict::Disagree,
            "{report:?}"
        );
    }
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
fn point_grade_certifies_each_face_in_its_own_clamp() {
    let tray = recipe(
        vec![stock(), shell(json!(open_sides()), point_grade())],
        "wall",
    );
    let mut built = build(&tray);
    volume(&built, "240");
    probes(
        &tray,
        &mut built,
        &[
            (["5", "4", "1"], Location::Boundary),
            (["5", "4", "2"], Location::Outside),
            (["5", "4", "4"], Location::Boundary),
            (["5", "4", "5"], Location::Inside),
        ],
    );
    let sweep = recipe(vec![stock(), thicken("front", point_grade())], "wall");
    let mut built = build(&sweep);
    volume(&built, "160");
    probes(
        &sweep,
        &mut built,
        &[
            (["5", "4", "8"], Location::Boundary),
            (["1", "1", "15/2"], Location::Inside),
            (["5", "4", "81/10"], Location::Outside),
        ],
    );
}

#[test]
fn radial_grade_varies_on_both_sides_of_an_axis_in_the_face() {
    let mut open = open_sides();
    open.push(json!({"face":{"step":"stock","role":"start"}}));
    let tray = recipe(vec![stock(), shell(json!(open), radial_grade())], "wall");
    let mut built = build(&tray);
    volume(&built, "100");
    probes(
        &tray,
        &mut built,
        &[
            (["5", "4", "5"], Location::Boundary),
            (["5", "2", "19/4"], Location::Boundary),
            (["5", "6", "19/4"], Location::Boundary),
            (["5", "4", "19/4"], Location::Outside),
        ],
    );
    let sweep = recipe(vec![stock(), thicken("back", radial_grade())], "wall");
    let mut built = build(&sweep);
    volume(&built, "100");
    probes(
        &sweep,
        &mut built,
        &[
            (["5", "4", "19/4"], Location::Outside),
            (["5", "2", "19/4"], Location::Boundary),
            (["5", "6", "19/4"], Location::Boundary),
            (["5", "1", "19/4"], Location::Inside),
        ],
    );
}

#[test]
fn oblique_shell_projects_onto_each_wall() {
    let grade = json!({"linear":{"from":["0","0","0"],"to":["10","10","0"],
        "start":"1","end":"3"}});
    let tray = recipe(vec![stock(), shell(json!([]), grade)], "wall");
    let mut built = build(&tray);
    probes(
        &tray,
        &mut built,
        &[
            (["5", "3", "41/10"], Location::Outside),
            (["5", "3", "21/5"], Location::Boundary),
            (["5", "3", "43/10"], Location::Inside),
            (["5", "3/2", "3"], Location::Boundary),
            (["7/5", "4", "3"], Location::Boundary),
        ],
    );
}

#[test]
fn oblique_thicken_certifies_both_clamps_and_accepts_a_later_cut() {
    let grade = json!({"linear":{"from":["2","2","0"],"to":["6","6","0"],
        "start":"1","end":"2"}});
    let sweep = thicken("front", grade);
    let thickened = recipe(vec![stock(), sweep.clone()], "wall");
    let mut built = build(&thickened);
    probes(
        &thickened,
        &mut built,
        &[
            (["1", "1", "7"], Location::Boundary),
            (["4", "4", "15/2"], Location::Boundary),
            (["8", "6", "8"], Location::Boundary),
            (["4", "4", "151/20"], Location::Outside),
        ],
    );
    let cut = recipe(
        vec![
            stock(),
            sweep,
            json!({"id":"tool","op":"block","size":["1","1","3"],"pose":pose(["4","4","6"])}),
            json!({"id":"cut","op":"boolean","kind":"subtract","target":"wall","tool":"tool"}),
        ],
        "cut",
    );
    let cut = Recipe::parse(&cut.to_text()).unwrap();
    let mut built = build(&cut);
    probes(
        &cut,
        &mut built,
        &[
            (["9/2", "9/2", "7"], Location::Outside),
            (["9/2", "11/2", "7"], Location::Inside),
        ],
    );
}

#[test]
fn uncertified_selection_grade_and_source_have_named_refusals() {
    let mut multi = thicken("front", grade("0", "10", "3"));
    multi["faces"] = json!([{}]);
    let mut curved = point_grade();
    curved["point"]["r1"] = json!("10");
    let cylinder =
        json!({"id":"stock","op":"cylinder","radius":"4","height":"6","pose":pose(["0","0","0"])});
    for (source, op, code) in [
        (stock(), multi, "distance.grade-selection-unsupported"),
        (
            stock(),
            thicken("front", curved),
            "distance.grade-curved-wall-unsupported",
        ),
        (
            stock(),
            shell(json!([]), radial_grade()),
            "distance.grade-curved-wall-unsupported",
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
