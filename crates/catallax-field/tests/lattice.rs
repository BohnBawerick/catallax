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

fn mirrored_pair(cell: serde_json::Value, radius: serde_json::Value, twice: bool) -> (Plan, Plan) {
    let lattice = |extra: Vec<serde_json::Value>, result: &str| {
        let mut steps = vec![
            serde_json::json!({"id":"region","op":"block","size":["4","4","4"],
                "pose":{"rotation":["1","0","0","0"],"translation":["0.5","0","0"]}}),
            serde_json::json!({"id":"beams","op":"lattice","region":"region","cell":cell,"size":"2",
                "frame":{"rotation":["3","0","0","1"],"translation":["0.25","-0.5","0.75"]},
                "radius":radius}),
        ];
        steps.extend(extra);
        let recipe = serde_json::json!({"format":"catallax-recipe/3","name":"lattice mirror",
            "steps":steps,"result":result});
        Plan::compile(&Recipe::parse(&recipe.to_string()).unwrap()).unwrap()
    };
    let mirror = |id: &str, body: &str| {
        serde_json::json!({"id":id,"op":"mirror","body":body,
            "plane":{"point":["1","0","0"],"normal":["1","2","0"]},"keep_original":false})
    };
    let mut extra = vec![mirror("image", "beams")];
    let mut result = "image";
    if twice {
        extra.push(mirror("back", "image"));
        result = "back";
    }
    (lattice(vec![], "beams"), lattice(extra, result))
}

fn reflect(p: &[Rational; 3]) -> [Rational; 3] {
    // x - 2 ((x - a) . m / (m . m)) m with a = (1,0,0), m = (1,2,0)
    let k = &(&(&(&p[0] - &q("1")) + &(&q("2") * &p[1])) * &q("2")) / &q("5");
    [&p[0] - &k, &p[1] - &(&k * &q("2")), p[2].clone()]
}

fn assert_mirrored(cell: serde_json::Value, radius: serde_json::Value) {
    let (original, image) = mirrored_pair(cell, radius, false);
    let mut a = Workspace::default();
    let mut b = Workspace::default();
    // Local frame coordinates, turned and shifted into the world by the frame.
    let values = ["0", "0.1", "0.25", "1/2", "1", "1.1", "2", "3"].map(q);
    let mut inside = 0;
    let mut outside = 0;
    for index in 0..512 {
        let l = [
            values[index % 8].clone(),
            values[index / 8 % 8].clone(),
            values[index / 64].clone(),
        ];
        let p = [
            &(&(&q("0.8") * &l[0]) - &(&q("0.6") * &l[1])) + &q("0.25"),
            &(&(&q("0.6") * &l[0]) + &(&q("0.8") * &l[1])) - &q("0.5"),
            &l[2] + &q("0.75"),
        ];
        let expected = original.proven(
            &RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()),
            &mut a,
        );
        let r = reflect(&p);
        let got = image.proven(
            &RPoint3::new(r[0].clone(), r[1].clone(), r[2].clone()),
            &mut b,
        );
        assert_eq!(got.location, expected.location, "{p:?}");
        inside += (expected.location == Location::Inside) as usize;
        outside += (expected.location == Location::Outside) as usize;
    }
    assert!(inside > 0 && outside > 0);
}

#[test]
fn mirrored_cubic_lattice_is_the_reflected_point_set() {
    assert_mirrored(
        serde_json::json!("cubic"),
        serde_json::json!({"constant":"0.25"}),
    );
}

#[test]
fn mirrored_asymmetric_custom_beams_reflect_with_their_linear_grade() {
    assert_mirrored(
        serde_json::json!({"custom":{
            "nodes":[["0","0","0"],["1","0","0"],["0","1","1/2"],["1/2","1/2","1"]],
            "beams":[[0,1],[0,2],[1,3]]}}),
        serde_json::json!({"linear":{"from":["0","0","0"],"to":["1","2","1"],"start":"0.1","end":"0.4"}}),
    );
}

#[test]
fn mirroring_a_lattice_twice_restores_the_original_answers() {
    let (original, twice) = mirrored_pair(
        serde_json::json!("body_centred"),
        serde_json::json!({"linear":{"from":["0","0","0"],"to":["0","1","0"],"start":"0.2","end":"0.4"}}),
        true,
    );
    let mut a = Workspace::default();
    let mut b = Workspace::default();
    let values = ["-1", "0", "0.25", "1", "1.5", "2.25"].map(q);
    for index in 0..72 {
        let p = RPoint3::new(
            values[index % 6].clone(),
            values[index / 6 % 6].clone(),
            values[index * 5 % 6].clone(),
        );
        assert_eq!(
            original.proven(&p, &mut a).location,
            twice.proven(&p, &mut b).location,
            "{p:?}"
        );
    }
}

#[test]
fn mirrored_lattice_keeps_boundary_contact_exact() {
    // Beam axis x-grid line through the frame origin at radius 1/4: the point
    // 1/4 off the line is on the boundary, and so is its image.
    let (original, image) = mirrored_pair(
        serde_json::json!("cubic"),
        serde_json::json!({"constant":"0.25"}),
        false,
    );
    let mut a = Workspace::default();
    let mut b = Workspace::default();
    // Local (1, 1/4, 0) on the frame's x beam, turned by (0.8, 0.6) and shifted.
    let p = [q("0.9"), q("0.3"), q("0.75")];
    let o = original.proven(
        &RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()),
        &mut a,
    );
    let r = reflect(&p);
    let m = image.proven(
        &RPoint3::new(r[0].clone(), r[1].clone(), r[2].clone()),
        &mut b,
    );
    assert_eq!(o.location, Location::Boundary);
    assert_eq!(m.location, Location::Boundary);
}
