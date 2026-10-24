//! Exact periodic beam membership regressions.
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn plan(translation: &str) -> Plan {
    let recipe = format!(
        r#"{{"format":"catallax-recipe/3","name":"lattice classification","steps":[{{"id":"region","op":"block","size":["4","4","4"],"pose":{{"rotation":["1","0","0","0"],"translation":["{translation}","0","0"]}}}},{{"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2","frame":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},"radius":{{"constant":"0.25"}}}}],"result":"beams"}}"#
    );
    Plan::compile(&Recipe::parse(&recipe).unwrap()).unwrap()
}

#[test]
fn periodic_beams_are_clipped_to_the_region() {
    let plan = plan("0");
    let mut storage = Workspace::default();
    let mut dec = Decider::default();
    for (xyz, expected) in [
        (["2", "1", "2"], Location::Inside),
        (["1", "1", "1"], Location::Outside),
        (["0.25", "1", "2"], Location::Boundary),
        (["-2", "1", "2"], Location::Outside),
    ] {
        let [x, y, z] = xyz.map(q);
        let answer = plan
            .proven_checked(&RPoint3::new(x, y, z), &mut storage, &mut dec)
            .unwrap();
        assert_eq!(answer.location, expected, "{xyz:?}");
    }
}

#[test]
fn periodic_cell_indices_are_not_limited_to_machine_integers() {
    let offset = q("100000000000000000000");
    let plan = plan("100000000000000000000");
    let mut storage = Workspace::default();
    let point = RPoint3::new(&offset + &q("2"), q("1"), q("2"));
    assert_eq!(plan.proven(&point, &mut storage).location, Location::Inside);
}

#[test]
fn linear_radius_is_read_at_the_query_point_and_clamped() {
    let recipe = Recipe::parse(r#"{
        "format":"catallax-recipe/3","name":"graded beams",
        "steps":[
            {"id":"region","op":"block","size":["4","4","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
             "radius":{"linear":{"from":["0","0","0"],"to":["0","1","0"],"start":"0.2","end":"0.6"}}}
        ],"result":"beams"
    }"#).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let mut storage = Workspace::default();
    for (y, expected) in [
        ("0.25", Location::Inside),
        ("1/3", Location::Boundary),
        ("3", Location::Outside),
    ] {
        let point = RPoint3::new(q("1"), q(y), q("2"));
        assert_eq!(
            plan.proven(&point, &mut storage).location,
            expected,
            "y={y}"
        );
    }
}

#[test]
fn cubic_grid_line_reduction_matches_periodic_segments_across_negative_cells() {
    let mut recipe = serde_json::json!({
        "format":"catallax-recipe/3","name":"cubic equivalence",
        "steps":[
            {"id":"region","op":"block","size":["8","8","8"],
             "pose":{"rotation":["1","0","0","0"],"translation":["-4","-4","-4"]}},
            {"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2",
             "frame":{"rotation":["3","0","0","1"],"translation":["0.25","-0.5","0.75"]},
             "radius":{"constant":"0.25"}}],"result":"beams"
    });
    let cubic = Plan::compile(&Recipe::parse(&recipe.to_string()).unwrap()).unwrap();
    recipe["steps"][1]["cell"] = serde_json::json!({"custom":{
        "nodes":[["0","0","0"],["1","0","0"],["0","1","0"],["0","0","1"]],
        "beams":[[0,1],[0,2],[0,3]]
    }});
    let generic = Plan::compile(&Recipe::parse(&recipe.to_string()).unwrap()).unwrap();
    let values = [
        "-2.25", "-2", "-1.75", "-1", "-0.25", "0", "0.25", "1", "1.75", "2", "2.25",
    ]
    .map(q);
    let mut a = Workspace::default();
    let mut b = Workspace::default();
    for index in 0..48 {
        let x = &values[index % 8];
        let y = &values[index / 8];
        let z = &values[index * 7 % values.len()];
        let point = RPoint3::new(
            &(&(&q("0.8") * x) - &(&q("0.6") * y)) + &q("0.25"),
            &(&(&q("0.6") * x) + &(&q("0.8") * y)) - &q("0.5"),
            z + &q("0.75"),
        );
        assert_eq!(
            cubic.proven(&point, &mut a).location,
            generic.proven(&point, &mut b).location,
            "local {x},{y},{z}"
        );
    }
}
