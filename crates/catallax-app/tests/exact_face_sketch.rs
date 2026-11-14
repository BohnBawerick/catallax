//! Face-relative tools retain their supporting plane in both kernel interpretations.
use catallax_boolean::{Containment, classify_point};
use catallax_numbers::{Decider, Outcome, Rational, Real};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::{Point3, Surface};
use serde_json::{Value, json};

fn lines(points: &[[i32; 2]]) -> Value {
    json!([(0..points.len())
        .map(|i| json!({
            "kind":"line", "id":format!("s{i}"),
            "to":points[(i+1)%points.len()].map(|v| v.to_string())
        }))
        .collect::<Vec<_>>()])
}

fn recipe(angled: bool, sides: usize, offset: i32) -> Recipe {
    let points = if angled {
        vec![
            [-55, -10],
            [-65, -15],
            [-65, -25],
            [-55, -30],
            [-45, -25],
            [-45, -15],
        ]
    } else {
        vec![[45, 45], [35, 40], [35, 30], [45, 25], [55, 30], [55, 40]]
    };
    let face = if angled {
        json!({"side":{"segment":"s1"}})
    } else {
        json!("end")
    };
    Recipe::parse(&json!({"format":"catallax-recipe/5","name":"face cut","steps":[
        {"id":"base","op":"sketch","plane":{"rotation":["1","1","0","0"],"translation":["0","0","0"]},
            "loops":[[
                {"kind":"line","id":"s0","to":["150.08645629882812","0"]},
                {"kind":"line","id":"s1","to":["0","120"]},
                {"kind":"line","id":"s2","to":["0","0"]}
            ]]},
        {"id":"boss","op":"extrude","profile":"base","distance":"105.5"},
        {"id":"section","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
            "loops":lines(&points[..sides])},
        {"id":"tool","op":"extrude_on_face","profile":"section","body":"boss",
            "face":{"step":"boss","role":face},
            "normal_direction":if angled { ["1","0","1"] } else { ["0","-1","0"] },
            "x_direction":if angled { ["0","1","0"] } else { ["1","0","0"] },
            "offset":offset.to_string(),"distance":"-150"},
        {"id":"cut","op":"boolean","kind":"subtract","target":"boss","tool":"tool"}
    ],"result":"cut"}).to_string()).unwrap()
}

fn zero(value: &Real, dec: &mut Decider) {
    assert!(
        matches!(dec.sign(value), Outcome::ProvenCoincident(_)),
        "not exactly zero: {value:?}"
    );
}

fn opens(angled: bool, sides: usize) {
    let recipe = recipe(angled, sides, 0);
    let mut dec = Decider::default();
    let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
    let tool_step = StepId::new("tool").unwrap();
    let boss = built.bodies[&StepId::new("boss").unwrap()];
    let tool = built.bodies[&tool_step];
    let cap = built
        .model
        .faces_of_body(tool)
        .into_iter()
        .find(|id| {
            built
                .model
                .face_origins(*id)
                .iter()
                .any(|o| o.step == tool_step && o.role == Role::Start)
        })
        .unwrap();
    let Surface::Plane(start) = built.model.geom().surface(built.model.face(cap).surface) else {
        panic!()
    };
    let start = start.clone();
    let source = built
        .model
        .faces_of_body(boss)
        .into_iter()
        .find(|id| {
            built.model.face_origins(*id).iter().any(|o| {
                if angled {
                    matches!(&o.role, Role::Side { segment, .. } if segment.as_str()=="s1")
                } else {
                    o.role == Role::End
                }
            })
        })
        .unwrap();
    let Surface::Plane(support) = built.model.geom().surface(built.model.face(source).surface)
    else {
        panic!()
    };
    zero(
        &support.normal.cross(&start.normal).norm_squared(),
        &mut dec,
    );
    zero(
        &support.normal.dot(&start.origin.minus(&support.origin)),
        &mut dec,
    );
    assert!(
        !built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .any(|id| {
                built
                    .model
                    .face_origins(id)
                    .iter()
                    .any(|o| o.step == tool_step && o.role == Role::Start)
            }),
        "the cut kept its start cap"
    );

    let x = if angled { [0, 1, 0] } else { [1, 0, 0] };
    let map = catallax_features::face_placement(
        std::slice::from_ref(support),
        &if angled {
            catallax_shapes::RVector3::new(1.into(), 0.into(), 1.into())
        } else {
            catallax_shapes::RVector3::new(0.into(), (-1).into(), 0.into())
        },
        &catallax_shapes::RVector3::new(x[0].into(), x[1].into(), x[2].into()),
        &Rational::zero(),
        &catallax_recipe::Pose::identity(),
        &mut dec,
    )
    .unwrap();
    let (u, v) = if angled { (-55, -20) } else { (45, 35) };
    for depth in [0, -1] {
        let point = map.apply_point(&Point3::new(u.into(), v.into(), depth.into()));
        if depth == -1 {
            assert!(matches!(
                classify_point(&mut built.model, boss, &point, &mut dec),
                Outcome::Proven(Containment::Inside)
            ));
        }
        assert!(
            matches!(
                classify_point(&mut built.model, built.result, &point, &mut dec),
                Outcome::Proven(Containment::Outside)
            ),
            "material remains across the sketch face at depth {depth}"
        );
    }
    let field = catallax_field::Plan::compile_checked(&recipe, &mut Decider::default()).unwrap();
    let point = if angled {
        [74, -55, 60]
    } else {
        [45, -105, 35]
    };
    let point = catallax_shapes::RPoint3::new(point[0].into(), point[1].into(), point[2].into());
    let query = field
        .proven_checked(
            &point,
            &mut catallax_field::Workspace::default(),
            &mut Decider::default(),
        )
        .unwrap();
    assert_eq!(query.location, catallax_field::Location::Outside);
}

#[test]
fn six_sided_cut_opens_an_irrationally_oriented_face() {
    opens(true, 6);
}

#[test]
fn five_sided_cut_opens_an_irrationally_oriented_face() {
    opens(true, 5);
}

#[test]
fn flat_face_control_stays_open() {
    opens(false, 6);
}

#[test]
fn a_pocket_floor_keeps_the_selected_facing_in_both_kernels() {
    let mut input = serde_json::to_value(recipe(false, 6, 0)).unwrap();
    input["steps"][3]["distance"] = json!("-10");
    let steps = input["steps"].as_array_mut().unwrap();
    steps.push(json!({"id":"floor_section","op":"sketch",
        "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
        "loops":lines(&[[43,33],[47,33],[47,37],[43,37]])}));
    steps.push(
        json!({"id":"floor_tool","op":"extrude_on_face","profile":"floor_section",
        "body":"cut","face":{"step":"tool","role":"end"},
        "normal_direction":["0","-1","0"],"x_direction":["1","0","0"],
        "offset":"0","distance":"-5"}),
    );
    input["result"] = json!("floor_tool");
    let recipe = Recipe::parse(&input.to_string()).unwrap();
    let mut dec = Decider::default();
    let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
    let point = catallax_shapes::RPoint3::new(45.into(), (-93).into(), 35.into());
    assert!(matches!(
        classify_point(&mut built.model, built.result, &point, &mut dec),
        Outcome::Proven(Containment::Inside)
    ));
    let field = catallax_field::Plan::compile_checked(&recipe, &mut dec).unwrap();
    let query = field
        .proven_checked(&point, &mut catallax_field::Workspace::default(), &mut dec)
        .unwrap();
    assert_eq!(query.location, catallax_field::Location::Inside);
}

#[test]
fn a_real_inward_offset_keeps_the_outer_face() {
    let mut dec = Decider::default();
    let built = catallax_features::evaluate(&recipe(true, 6, -2), &mut dec).unwrap();
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .any(|id| {
                built
                    .model
                    .face_origins(id)
                    .iter()
                    .any(|o| o.step.as_str() == "tool" && o.role == Role::Start)
            })
    );
}
