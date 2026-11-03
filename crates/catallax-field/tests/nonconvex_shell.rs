//! Definition-only shell material, finite faces and surviving names.

use catallax_field::{Direction, FaceId, FaceState, Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_recipe::{Origin, Recipe};
use catallax_shapes::{RPoint3, distance::Options};
use serde_json::{Value, json};

const WORK: u64 = 100_000_000;

fn document() -> Value {
    serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/nonconvex-shell.json"
    )).unwrap()
}
fn compile(source: &Value) -> Plan {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let (answer, work) = with_operation_budget(WORK, || Plan::compile(&recipe));
    assert!(!work.exhausted, "{work:?}");
    answer.unwrap()
}
fn point(xyz: [&str; 3]) -> RPoint3 {
    let [x, y, z] = xyz.map(|s| Rational::parse(s).unwrap());
    RPoint3::new(x, y, z)
}
fn material(plan: &Plan, cases: &[([&str; 3], Location)]) {
    for &(xyz, expected) in cases {
        let (answer, work) = with_operation_budget(WORK, || {
            plan.proven_checked(&point(xyz), &mut Workspace::default(), &mut Decider::default())
        });
        assert!(!work.exhausted, "{xyz:?}: {work:?}");
        let answer = answer.unwrap();
        assert_eq!(answer.location, expected, "{xyz:?}: {answer:?}");
    }
}
fn origin(value: Value) -> Origin {
    serde_json::from_value(value).unwrap()
}
fn round(a: &str, b: &str) -> Origin {
    origin(json!({"step":"shell","role":{"round":{"edge":{"between":[
        {"step":"stock","role":{"side":{"segment":a}}},
        {"step":"stock","role":{"side":{"segment":b}}}
    ]}}}}))
}
fn named_on(plan: &Plan, name: &Origin, xyz: [&str; 3]) -> (FaceId, Direction) {
    for (index, face) in plan.faces().iter().enumerate() {
        if &face.origin != name {
            continue;
        }
        let Some(handle) = plan.result_face(index) else { continue };
        let query = plan.face_query(handle, &point(xyz), &mut Workspace::default()).unwrap();
        if query.state == FaceState::OnFace {
            return (handle, query.direction);
        }
    }
    panic!("{name:?} has no surviving face at {xyz:?}");
}
fn open_except(source: &mut Value, segment: &str) {
    let mut open = vec![json!({"face":{"step":"stock","role":"start"}}),
                        json!({"face":{"step":"stock","role":"end"}})];
    for edge in source["steps"][0]["loops"].as_array().unwrap().iter()
        .flat_map(|cycle| cycle.as_array().unwrap())
    {
        if edge["id"] != segment {
            open.push(json!({"face":{"step":"stock","role":{"side":{"segment":edge["id"]}}}}));
        }
    }
    source["steps"][2]["open"] = json!(open);
}
fn chain(points: &[[i64; 2]], prefix: &str) -> Value {
    json!(points.iter().enumerate().map(|(i, _)| json!({
        "kind":"line", "id":format!("{prefix}-{i}"),
        "to":points[(i + 1) % points.len()].map(|n| n.to_string())
    })).collect::<Vec<_>>())
}

#[test]
fn reentrant_tray_compiles_independently_with_finite_corner_material() {
    let plan = compile(&document());
    material(&plan, &[
        (["2.1", "2.1", "2"], Location::Outside),
        (["2.5", "2.5", "2"], Location::Inside),
        (["4", "1.5", "2"], Location::Outside),
        (["4", "0.5", "2"], Location::Inside),
        (["4", "1.5", "0.5"], Location::Inside),
        (["4", "4", "2"], Location::Outside),
        (["12/5", "11/5", "2"], Location::Boundary),
    ]);
    let (_, direction) = named_on(&plan, &round("notch-bottom", "notch-left"), ["12/5", "11/5", "2"]);
    let Direction::Unique(normal) = direction else { panic!("{direction:?}") };
    assert!(normal[0].lo() > -0.601 && normal[0].hi() < -0.599);
    assert!(normal[1].lo() > -0.801 && normal[1].hi() < -0.799);
    assert!(normal[2].contains(0.0));
}

#[test]
fn kept_planar_face_includes_material_on_both_normal_sides() {
    let mut source = document();
    open_except(&mut source, "notch-bottom");
    let plan = compile(&source);
    material(&plan, &[
        (["2.5", "3.5", "2"], Location::Inside),
        (["2", "3.5", "2"], Location::Outside),
        (["4", "2.5", "2"], Location::Inside),
        (["1", "2.5", "2"], Location::Outside),
    ]);
    named_on(&plan, &round("notch-bottom", "notch-left"), ["12/5", "19/5", "2"]);
}

#[test]
fn material_across_a_gap_keeps_the_positive_parallel_face_and_finite_trim() {
    let mut source = document();
    source["steps"][0]["loops"] = json!([chain(
        &[[0,0],[8,0],[8,8],[5,8],[5,2],[3,2],[3,8],[0,8]], "u"
    )]);
    open_except(&mut source, "u-5");
    source["steps"][2]["thickness"] = json!("3");
    let plan = compile(&source);
    material(&plan, &[
        (["5.5", "5", "2"], Location::Inside),
        (["6.5", "5", "2"], Location::Outside),
        (["1", "5", "2"], Location::Inside),
        (["4", "5", "2"], Location::Outside),
        (["6", "5", "2"], Location::Boundary),
    ]);
    let name = origin(json!({"step":"shell","role":{"inner":{"of":{
        "step":"stock","role":{"side":{"segment":"u-5"}}
    }}}}));
    let (handle, direction) = named_on(&plan, &name, ["6", "5", "2"]);
    let Direction::Unique(normal) = direction else { panic!("{direction:?}") };
    assert!(normal[0].lo() > 0.999 && normal[0].hi() < 1.001);
    assert!(normal[1].contains(0.0) && normal[2].contains(0.0));
    let options = Options { width: Rational::parse("1/100").unwrap(), max_cells: 16_384 };
    let mut dec = Decider::default();
    let prepared = plan.prepare_distances(&mut dec).unwrap();
    let (answer, work) = with_operation_budget(WORK, || {
        prepared.point_to_face(handle, &point(["7", "9", "2"]).to_real(), &options, &mut dec)
    });
    assert!(!work.exhausted, "{work:?}");
    let answer = answer.unwrap();
    let expected = dec.enclose_at(&Real::from_i64(2).sqrt(), 128).unwrap();
    let upper = answer.upper.as_ref().expect("finite physical witness");
    assert!(answer.lower <= expected.lo() && upper >= &expected.hi(), "{answer:?}");
    assert!(answer.converged && upper - &answer.lower <= options.width, "{answer:?}");
}

#[test]
fn holed_caps_and_consumed_shells_keep_the_original_void() {
    let mut source = document();
    source["steps"][0]["loops"] = json!([
        chain(&[[0,0],[10,0],[10,10],[0,10]], "outer"),
        chain(&[[3,3],[3,7],[7,7],[7,3]], "hole")
    ]);
    source["steps"][2]["open"] = json!([
        {"face":{"step":"stock","role":"start"}},
        {"face":{"step":"stock","role":"end"}}
    ]);
    source["steps"][2]["thickness"] = json!("3/2");
    material(&compile(&source), &[
        (["1.6", "1.6", "2"], Location::Outside),
        (["1.4", "5", "2"], Location::Inside),
        (["5", "5", "2"], Location::Outside),
    ]);
    source["steps"][2]["thickness"] = json!("2");
    let plan = compile(&source);
    material(&plan, &[
        (["1.6", "1.6", "2"], Location::Inside),
        (["5", "5", "2"], Location::Outside),
        (["0", "5", "2"], Location::Boundary),
    ]);
    let inner = origin(json!({"step":"shell","role":{"inner":{"of":{
        "step":"stock","role":{"side":{"segment":"outer-3"}}
    }}}}));
    assert_eq!(plan.face_query(plan.select(&inner).unwrap(), &point(["2", "5", "2"]),
        &mut Workspace::default()).unwrap().state, FaceState::NotOnFace);
    let curve = round("hole-0", "hole-3");
    assert_eq!(plan.face_query(plan.select(&curve).unwrap(), &point(["1", "3", "2"]),
        &mut Workspace::default()).unwrap().state, FaceState::NotOnFace);
}
