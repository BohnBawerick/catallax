//! Distances must use surviving trimmed faces, not merely carrier residuals.
use catallax_field::Plan;
use catallax_numbers::{Decider, Rational, Real};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, distance::Options};
use serde_json::{Value, json};
fn pose(x: &str) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":[x,"0","0"]})
}
fn recipe(steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"trimmed-distance","steps":steps,"result":result}).to_string()).unwrap()
}
fn point(x: i64, y: i64, z: i64) -> catallax_shapes::Point3 {
    RPoint3::new(x.into(), y.into(), z.into()).into()
}
fn controls() -> Options {
    Options {
        width: Rational::from_ratio(1, 1000).unwrap(),
        max_cells: 4096,
    }
}
#[test]
fn cylindrical_wall_distance_respects_height_trims_and_cap_pair_distance() {
    let recipe = recipe(vec![
        json!({"id":"stock","op":"cylinder","radius":"1","height":"2","pose":pose("0")}),
    ]);
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let prepared = plan.prepare_distances(&mut dec).unwrap();
    let face = |role: Role| {
        plan.faces()
            .iter()
            .find(|f| f.origin.role == role)
            .map(|f| plan.select(&f.origin).unwrap())
            .unwrap()
    };
    let wall = plan
        .faces()
        .iter()
        .find(|f| matches!(f.origin.role, Role::Side { .. }))
        .map(|f| plan.select(&f.origin).unwrap())
        .unwrap();
    let answer = prepared
        .point_to_face(wall, &point(1, 0, 3), &controls(), &mut dec)
        .unwrap();
    assert!(answer.lower <= Rational::one(), "{answer:?}");
    assert!(
        answer.lower >= Rational::parse("0.999").unwrap(),
        "{answer:?}"
    );
    assert!(answer.upper.as_ref().unwrap() >= &Rational::one());
    assert!(
        answer.upper.as_ref().unwrap() <= &Rational::parse("1.001").unwrap(),
        "{answer:?}"
    );
    let pair = prepared
        .face_to_face(face(Role::Start), face(Role::End), &controls(), &mut dec)
        .unwrap();
    assert_eq!(pair.lower, Rational::from_i64(2));
    assert_eq!(pair.upper, Some(Rational::from_i64(2)));
}
#[test]
fn cap_hole_rejects_the_plane_projection() {
    let recipe = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0"),"loops":[[{"kind":"circle","id":"outer","centre":["0","0"],"radius":"2"}],[{"kind":"circle","id":"inner","centre":["0","0"],"radius":"1"}]]}),
        json!({"id":"feature","op":"extrude","profile":"profile","distance":"2"}),
    ]);
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let cap = plan
        .faces()
        .iter()
        .find(|f| f.origin.role == Role::Start)
        .map(|f| plan.select(&f.origin).unwrap())
        .unwrap();
    let answer = plan
        .prepare_distances(&mut dec)
        .unwrap()
        .point_to_face(cap, &point(0, 0, -1), &controls(), &mut dec)
        .unwrap();
    let expected = dec.enclose_at(&Real::from_i64(2).sqrt(), 128).unwrap();
    assert!(answer.lower <= expected.lo(), "{answer:?}");
    assert!(
        answer.upper.as_ref().unwrap() >= &expected.hi(),
        "{answer:?}"
    );
    assert!(answer.lower > Rational::parse("1.4").unwrap(), "{answer:?}");
}
#[test]
fn a_buried_tangent_source_cannot_supply_a_face_distance_witness() {
    let recipe = recipe(vec![
        json!({"id":"outer","op":"sphere","radius":"2","pose":pose("0")}),
        json!({"id":"inner","op":"sphere","radius":"1","pose":pose("1")}),
        json!({"id":"result","op":"boolean","kind":"union","target":"outer","tool":"inner"}),
    ]);
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let inner = plan
        .faces()
        .iter()
        .find(|f| f.origin.step.as_str() == "inner")
        .map(|f| plan.select(&f.origin).unwrap())
        .unwrap();
    let answer = plan
        .prepare_distances(&mut dec)
        .unwrap()
        .point_to_face(
            inner,
            &point(2, 0, 0),
            &Options {
                max_cells: 64,
                ..controls()
            },
            &mut dec,
        )
        .unwrap();
    assert!(answer.upper.is_none(), "{answer:?}");
    assert!(!answer.converged);
}

#[test]
fn discarded_copy_alias_cannot_supply_a_distance_witness() {
    let recipe = recipe(vec![
        json!({"id":"seed","op":"sphere","radius":"1","pose":pose("0")}),
        json!({"id":"reflected","op":"mirror","body":"seed","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}),
        json!({"id":"image","op":"mirror","body":"reflected","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}),
    ]);
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let named = |step: &str| {
        plan.faces()
            .iter()
            .find(|f| f.origin.step.as_str() == step)
            .map(|f| plan.select(&f.origin).unwrap())
            .unwrap()
    };
    let prepared = plan.prepare_distances(&mut dec).unwrap();
    let error = prepared
        .point_to_face(named("seed"), &point(2, 0, 0), &controls(), &mut dec)
        .unwrap_err();
    assert_eq!(error.code, "field.face-id");
    let image = prepared
        .point_to_face(named("image"), &point(2, 0, 0), &controls(), &mut dec)
        .unwrap();
    assert_eq!(image.lower, Rational::one());
    assert_eq!(image.upper, Some(Rational::one()));
}

#[test]
fn canonical_cone_cap_distance_uses_its_generated_trim_role() {
    let recipe = recipe(vec![
        json!({"id":"stock","op":"cone","base_radius":"1","top_radius":"1/2","height":"1","pose":pose("0")}),
    ]);
    let step = &recipe.steps[0];
    let mut dec = Decider::default();
    let plan = Plan::compile_canonical_stock(&step.op, &step.id, &mut dec)
        .unwrap()
        .unwrap();
    let cap = plan
        .faces()
        .iter()
        .find(|face| face.origin.role == Role::End)
        .map(|face| plan.select(&face.origin).unwrap())
        .unwrap();
    let target = RPoint3::new(
        Rational::zero(),
        Rational::from_ratio(1, 4).unwrap(),
        Rational::from_i64(2),
    )
    .to_real();
    let answer = plan
        .prepare_distances(&mut dec)
        .unwrap()
        .point_to_face(cap, &target, &controls(), &mut dec)
        .unwrap();
    assert_eq!(answer.lower, Rational::one());
    assert_eq!(answer.upper, Some(Rational::one()));
}

#[test]
fn drilled_wall_distance_admits_the_surviving_algebraic_radial_witness() {
    let recipe = recipe(vec![
        json!({"id":"stock","op":"block","size":["10","10","10"],"pose":pose("-5")}),
        json!({"id":"drill","op":"hole","target":"stock","at":[{"id":"centre","point":["0","5","0"],"direction":["0","0","1"]}],"diameter":"2","depth":"through","tip":"flat"}),
    ]);
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let wall = plan
        .faces()
        .iter()
        .find(|f| {
            matches!(
                f.origin.role,
                Role::Hole {
                    part: catallax_recipe::HolePart::Wall,
                    ..
                }
            )
        })
        .map(|f| plan.select(&f.origin).unwrap())
        .unwrap();
    let options = controls();
    let answer = plan
        .prepare_distances(&mut dec)
        .unwrap()
        .point_to_face(wall, &point(2, 7, 5), &options, &mut dec)
        .unwrap();
    let expected = dec
        .enclose_at(&(Real::from_i64(8).sqrt() - Real::one()), 128)
        .unwrap();
    let upper = answer
        .upper
        .as_ref()
        .expect("surviving drilled-wall witness");
    assert!(answer.lower <= expected.hi(), "{answer:?}");
    assert!(upper >= &expected.lo(), "{answer:?}");
    assert!(upper - &answer.lower <= options.width, "{answer:?}");
}
