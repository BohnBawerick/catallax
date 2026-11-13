//! Exact moments of coaxial meridians, including crossings inside an axial slab.
use catallax_field::{
    Plan,
    assembly::{self, PlacedPart},
    mass,
};
use catallax_measure::Enclosure;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};
use serde_json::{Value, json};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn revolve(id: &str, points: [[&str; 2]; 4]) -> Vec<Value> {
    let sketch = format!("{id}-section");
    vec![
        json!({"id":sketch,"op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[points.iter().enumerate().map(|(i,p)|json!({"kind":"line","id":format!("{id}-e{i}"),"to":p})).collect::<Vec<_>>()]}),
        json!({"id":id,"op":"revolve","profile":sketch,"axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":"360"}),
    ]
}
fn plan(steps: Vec<Value>, result: &str) -> Plan {
    Plan::compile(
        &Recipe::parse(
            &json!({"format":"catallax-recipe/3","name":"meridians","steps":steps,"result":result})
                .to_string(),
        )
        .unwrap(),
    )
    .unwrap()
}
fn ring() -> Plan {
    plan(
        revolve("ring", [["1", "0"], ["2", "0"], ["2", "3"], ["1", "3"]]),
        "ring",
    )
}
fn pi_multiple(e: &Enclosure, factor: &str) {
    let factor = q(factor);
    let lo = &q(
        "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679",
    ) * &factor;
    let hi = &q(
        "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170680",
    ) * &factor;
    assert!(
        e.lo <= lo.clone().min(hi.clone()) && e.hi >= lo.max(hi),
        "{e:?} versus {factor}*pi"
    );
    assert!(e.width() <= q("1e-20"));
}

#[test]
fn annulus_volume_centroid_and_full_inertia_are_independent_analytic_values() {
    let p = mass::revolved_properties(
        &ring(),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    pi_multiple(&p.volume, "9");
    let centre = p.centroid.as_ref().unwrap();
    assert_eq!(
        centre.clone().map(|v| {
            assert!(v.is_exact());
            v.lo
        }),
        [q("0"), q("1.5"), q("0")]
    );
    let tensor = p.inertia.as_ref().unwrap();
    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                pi_multiple(&tensor[i][j], if i == 1 { "22.5" } else { "18" });
            } else {
                assert!(tensor[i][j].is_exact() && tensor[i][j].lo.is_zero());
            }
        }
    }
    assert!(p.contains(p.estimate.as_ref().unwrap()));
}

#[test]
fn boolean_meridians_split_at_the_internal_radial_order_crossing() {
    for (kind, volume, centre) in [
        ("intersect", "14/3", "1"),
        ("union", "38/3", "1"),
        ("subtract", "4", "5/3"),
    ] {
        let mut steps = revolve("a", [["0", "0"], ["1", "0"], ["3", "2"], ["0", "2"]]);
        steps.extend(revolve(
            "b",
            [["0", "0"], ["3", "0"], ["1", "2"], ["0", "2"]],
        ));
        steps.push(json!({"id":"result","op":"boolean","kind":kind,"target":"a","tool":"b"}));
        let p = mass::revolved_properties(
            &plan(steps, "result"),
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default(),
        )
        .unwrap();
        pi_multiple(&p.volume, volume);
        assert_eq!(p.centroid.unwrap()[1].lo, q(centre));
    }
}

#[test]
fn nested_meridian_loop_preserves_an_enclosed_ring_cavity() {
    let mut steps = revolve("ring", [["1", "0"], ["4", "0"], ["4", "3"], ["1", "3"]]);
    let hole = revolve("hole", [["2", "1"], ["2", "2"], ["3", "2"], ["3", "1"]]);
    steps[0]["loops"]
        .as_array_mut()
        .unwrap()
        .push(hole[0]["loops"][0].clone());
    let p = mass::revolved_properties(
        &plan(steps, "ring"),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    pi_multiple(&p.volume, "40");
    assert_eq!(p.centroid.unwrap()[1].lo, q("1.5"));
}

#[test]
fn rigid_placement_transforms_the_centroid_and_tensor_without_changing_mass() {
    let pose = Pose::from_rational(
        Rotation::from_quaternion(q("1"), q("0"), q("0"), q("1")).unwrap(),
        RVector3::new(q("10"), q("20"), q("30")),
    );
    let placed = PlacedPart::new(&ring(), &pose, &mut Decider::default()).unwrap();
    let p = mass::revolved_properties(
        placed.plan(),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    pi_multiple(&p.volume, "9");
    assert_eq!(
        p.centroid.unwrap().map(|v| v.lo),
        [q("8.5"), q("20"), q("30")]
    );
    let inertia = p.inertia.unwrap();
    for (i, value) in ["22.5", "18", "18"].iter().enumerate() {
        pi_multiple(&inertia[i][i], value);
    }
    let m = assembly::revolved_uniform_mass(
        &placed,
        &q("2"),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    pi_multiple(&m.grams, "18");
}

#[test]
fn oblique_rational_placement_preserves_all_tensor_components() {
    let pose = Pose::from_rational(
        Rotation::from_quaternion(q("2"), q("1"), q("2"), q("0")).unwrap(),
        RVector3::new(q("-3"), q("7"), q("2")),
    );
    let placed = PlacedPart::new(&ring(), &pose, &mut Decider::default()).unwrap();
    let p = mass::revolved_properties(
        placed.plan(),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let expected = pose
        .apply_rational_point(&RPoint3::new(q("0"), q("1.5"), q("0")))
        .unwrap();
    assert_eq!(
        p.centroid.unwrap().map(|v| v.lo),
        [expected.x, expected.y, expected.z]
    );
    let axis = pose
        .apply_rational_vector(&RVector3::new(q("0"), q("1"), q("0")))
        .unwrap();
    let axis = [axis.x, axis.y, axis.z];
    let tensor = p.inertia.unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let coefficient =
                &q(if i == j { "18" } else { "0" }) + &(&q("4.5") * &(&axis[i] * &axis[j]));
            pi_multiple(&tensor[i][j], &coefficient.to_string());
        }
    }
}

#[test]
fn exhausted_work_small_cell_limits_and_invalid_density_refuse() {
    let plan = ring();
    let result = mass::revolved_properties(
        &plan,
        &mass::Accuracy::default(),
        &mass::Options { max_cells: 1 },
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(result.code, "field.mass-budget");
    let (result, work) = with_operation_budget(0, || {
        mass::revolved_properties(
            &plan,
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default(),
        )
    });
    assert!(result.is_err() && work.exhausted);
    let placed = PlacedPart::new(&plan, &Pose::identity(), &mut Decider::default()).unwrap();
    assert_eq!(
        assembly::revolved_uniform_mass(
            &placed,
            &q("0"),
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "assembly.density"
    );
}

#[test]
fn empty_material_is_zero_and_partial_or_noncoaxial_revolutions_refuse() {
    let mut steps = revolve("a", [["1", "0"], ["2", "0"], ["2", "3"], ["1", "3"]]);
    steps.push(json!({"id":"empty","op":"boolean","kind":"subtract","target":"a","tool":"a"}));
    let p = mass::revolved_properties(
        &plan(steps, "empty"),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(p.volume.is_exact() && p.volume.lo.is_zero());
    assert!(p.centroid.is_none() && p.inertia.is_none());
    let mut steps = revolve("a", [["1", "0"], ["2", "0"], ["2", "3"], ["1", "3"]]);
    steps[1]["angle"] = json!("180");
    assert_eq!(
        mass::revolved_properties(
            &plan(steps, "a"),
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "field.mass-revolution-coverage"
    );
    let mut steps = revolve("a", [["1", "0"], ["2", "0"], ["2", "3"], ["1", "3"]]);
    let mut b = revolve("b", [["1", "0"], ["2", "0"], ["2", "3"], ["1", "3"]]);
    b[0]["plane"]["translation"] = json!(["1", "0", "0"]);
    b[1]["axis"]["point"] = json!(["1", "0", "0"]);
    steps.extend(b);
    steps.push(json!({"id":"both","op":"boolean","kind":"union","target":"a","tool":"b"}));
    assert_eq!(
        mass::revolved_properties(
            &plan(steps, "both"),
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "field.mass-revolution-coverage"
    );
}
