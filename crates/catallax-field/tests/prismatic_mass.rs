//! Independent analytic controls for line/arc CSG moments.
use catallax_field::{Plan, mass};
use catallax_measure::Enclosure;
use catallax_numbers::{Ball, Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn b(s: &str) -> Ball {
    Ball::from_rational(&q(s), 512)
}
fn cylinder(id: &str, xyz: [&str; 3], radius: &str, height: &str) -> Value {
    json!({"id":id,"op":"cylinder","radius":radius,"height":height,"pose":{"rotation":["1","0","0","0"],"translation":xyz}})
}
fn block(id: &str, xyz: [&str; 3], size: [&str; 3]) -> Value {
    json!({"id":id,"op":"block","size":size,"pose":{"rotation":["1","0","0","0"],"translation":xyz}})
}
fn boolean(id: &str, kind: &str, a: &str, c: &str) -> Value {
    json!({"id":id,"op":"boolean","kind":kind,"target":a,"tool":c})
}
fn plan(steps: Vec<Value>, result: &str) -> Plan {
    Plan::compile(&Recipe::parse(&json!({"format":"catallax-recipe/3","name":"analytic CSG","steps":steps,"result":result}).to_string()).unwrap()).unwrap()
}
fn properties(steps: Vec<Value>, result: &str) -> mass::Properties {
    let plan = plan(steps, result);
    let (p, work) = with_operation_budget(20_000_000, || {
        mass::prismatic_properties(
            &plan,
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    p.unwrap()
}
fn check(e: &Enclosure, expected: &Ball) {
    assert!(
        e.lo <= expected.hi() && expected.lo() <= e.hi,
        "{e:?} versus {expected:?}"
    );
    assert!(e.width() <= q("1e-20"));
}
fn multiple(pi: &Ball, coefficient: &str) -> Ball {
    pi * &b(coefficient)
}

#[test]
fn translated_cylinder_has_analytic_volume_centroid_and_all_inertia_entries() {
    let p = properties(vec![cylinder("c", ["5", "7", "11"], "2", "3")], "c");
    let pi = Ball::pi(512);
    check(&p.volume, &multiple(&pi, "12"));
    let centre = p.centroid.as_ref().unwrap();
    for (v, expected) in centre.iter().zip(["5", "7", "12.5"]) {
        assert!(v.contains(&q(expected)));
    }
    let tensor = p.inertia.as_ref().unwrap();
    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                check(
                    &tensor[i][j],
                    &multiple(&pi, if i == 2 { "24" } else { "21" }),
                );
            } else {
                assert!(tensor[i][j].contains(&q("0")));
            }
        }
    }
    assert!(p.contains(p.estimate.as_ref().unwrap()));
}

#[test]
fn clipped_quarter_disk_matches_centroid_and_mixed_inertia_formula() {
    let p = properties(
        vec![
            cylinder("c", ["0", "0", "0"], "1", "1"),
            block("b", ["0", "0", "0"], ["1", "1", "1"]),
            boolean("cut", "intersect", "c", "b"),
        ],
        "cut",
    );
    let pi = Ball::pi(512);
    check(&p.volume, &multiple(&pi, "1/4"));
    let mean = b("4").checked_div(&(&b("3") * &pi)).unwrap();
    let centre = p.centroid.unwrap();
    check(&centre[0], &mean);
    check(&centre[1], &mean);
    assert!(centre[2].contains(&q("0.5")));
    let inverse = b("1").checked_div(&pi).unwrap();
    let tensor = p.inertia.unwrap();
    let perpendicular = &multiple(&pi, "1/12") - &(&b("4/9") * &inverse);
    check(&tensor[0][0], &perpendicular);
    check(&tensor[1][1], &perpendicular);
    check(
        &tensor[2][2],
        &(&multiple(&pi, "1/8") - &(&b("8/9") * &inverse)),
    );
    let mixed = &(&b("4/9") * &inverse) - &b("1/8");
    check(&tensor[0][1], &mixed);
    check(&tensor[1][0], &mixed);
    assert!(tensor[0][2].contains(&q("0")) && tensor[1][2].contains(&q("0")));
}

#[test]
fn lens_union_intersection_and_difference_have_correct_signed_moments() {
    let pi = Ball::pi(512);
    let root = b("3").sqrt().unwrap();
    let intersection = &multiple(&pi, "4/3") - &root;
    for kind in ["intersect", "union", "subtract"] {
        let p = properties(
            vec![
                cylinder("a", ["0", "0", "0"], "1", "2"),
                cylinder("b", ["1", "0", "0"], "1", "2"),
                boolean("result", kind, "a", "b"),
            ],
            "result",
        );
        let volume = match kind {
            "intersect" => intersection.clone(),
            "union" => &multiple(&pi, "8/3") + &root,
            _ => &multiple(&pi, "2/3") + &root,
        };
        check(&p.volume, &volume);
        let centre = p.centroid.unwrap();
        if kind == "subtract" {
            check(
                &centre[0],
                &(-(&intersection * &b("0.5"))).checked_div(&volume).unwrap(),
            );
        } else {
            assert!(centre[0].contains(&q("0.5")));
        }
        assert!(centre[1].contains(&q("0")) && centre[2].contains(&q("1")));
    }
}

#[test]
fn redundant_collinear_drills_and_bridge_make_the_exact_capsule() {
    let mut steps = Vec::new();
    let mut root = String::new();
    for (i, x) in ["-2", "-1", "0", "1", "2"].iter().enumerate() {
        let id = format!("c{i}");
        steps.push(cylinder(&id, [x, "0", "0"], "1", "1"));
        if i == 0 {
            root = id;
        } else {
            let next = format!("union{i}");
            steps.push(boolean(&next, "union", &root, &id));
            root = next;
        }
    }
    steps.push(block("bridge", ["-2", "-1", "0"], ["4", "2", "1"]));
    steps.push(boolean("capsule", "union", &root, "bridge"));
    let p = properties(steps, "capsule");
    let pi = Ball::pi(512);
    check(&p.volume, &(&b("8") + &pi));
    let centre = p.centroid.unwrap();
    assert!(
        centre[0].contains(&q("0")) && centre[1].contains(&q("0")) && centre[2].contains(&q("0.5"))
    );
    let tensor = p.inertia.unwrap();
    for (i, (constant, coefficient)) in [("10/3", "1/3"), ("50/3", "13/3"), ("56/3", "9/2")]
        .iter()
        .enumerate()
    {
        check(&tensor[i][i], &(&b(constant) + &multiple(&pi, coefficient)));
    }
}

#[test]
fn empty_material_and_resource_limits_keep_their_own_outcomes() {
    let p = properties(
        vec![
            cylinder("a", ["0", "0", "0"], "1", "1"),
            cylinder("b", ["2", "0", "0"], "1", "1"),
            boolean("touch", "intersect", "a", "b"),
        ],
        "touch",
    );
    assert!(p.volume.is_exact() && p.volume.lo.is_zero() && p.centroid.is_none());
    let plan = plan(vec![cylinder("a", ["0", "0", "0"], "1", "1")], "a");
    assert_eq!(
        mass::prismatic_properties(
            &plan,
            &mass::Accuracy::default(),
            &mass::Options { max_cells: 1 },
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "field.mass-budget"
    );
    let (p, work) = with_operation_budget(0, || {
        mass::prismatic_properties(
            &plan,
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default(),
        )
    });
    assert!(p.is_err() && work.exhausted);
}

#[test]
fn disjoint_orthogonal_cylinders_have_correct_world_tensor_and_overlap_refuses() {
    let a = cylinder("a", ["0", "0", "0"], "1", "2");
    let mut other = cylinder("b", ["10", "0", "0"], "1", "2");
    other["pose"]["rotation"] = json!(["1", "1", "0", "0"]);
    let p = properties(
        vec![
            a.clone(),
            other.clone(),
            boolean("union", "union", "a", "b"),
        ],
        "union",
    );
    let pi = Ball::pi(512);
    check(&p.volume, &multiple(&pi, "4"));
    let centre = p.centroid.unwrap();
    for (v, want) in centre.iter().zip(["5", "-0.5", "0.5"]) {
        assert!(v.contains(&q(want)));
    }
    let tensor = p.inertia.unwrap();
    for (i, coefficient) in ["13/3", "619/6", "619/6"].iter().enumerate() {
        check(&tensor[i][i], &multiple(&pi, coefficient));
    }
    check(&tensor[0][1], &multiple(&pi, "10"));
    check(&tensor[0][2], &multiple(&pi, "10"));
    check(&tensor[1][2], &multiple(&pi, "-1"));
    other["pose"]["translation"] = json!(["0", "0", "0"]);
    let plan = plan(
        vec![a, other, boolean("overlap", "intersect", "a", "b")],
        "overlap",
    );
    assert_eq!(
        mass::prismatic_properties(
            &plan,
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "field.mass-prismatic-coverage"
    );
}

#[test]
fn unequal_disks_cover_internal_tangency_and_a_partial_lens() {
    let pi = Ball::pi(512);
    let p = properties(
        vec![
            cylinder("a", ["0", "0", "0"], "2", "1"),
            cylinder("b", ["1", "0", "0"], "1", "1"),
            boolean("inside", "intersect", "a", "b"),
        ],
        "inside",
    );
    check(&p.volume, &pi);
    assert!(p.centroid.unwrap()[0].contains(&q("1")));
    let p = properties(
        vec![
            cylinder("a", ["0", "0", "0"], "2", "1"),
            cylinder("b", ["2", "0", "0"], "1", "1"),
            boolean("lens", "intersect", "a", "b"),
        ],
        "lens",
    );
    let expected = &(&b("4") * &b("7/8").acos().unwrap()) + &b("1/4").acos().unwrap()
        - &b("15").sqrt().unwrap() * &b("1/2");
    check(&p.volume, &expected);
}

#[test]
fn rational_box_boolean_moments_remain_exact_at_zero_requested_width() {
    let plan = plan(
        vec![
            block("a", ["1/3", "2/3", "1/5"], ["1", "1", "1"]),
            block("b", ["1/3", "2/3", "1/5"], ["1/2", "1/2", "1"]),
            boolean("cut", "subtract", "a", "b"),
        ],
        "cut",
    );
    let exact = catallax_measure::Accuracy {
        relative: q("0"),
        absolute: q("0"),
    };
    let accuracy = mass::Accuracy {
        volume: exact.clone(),
        centroid: exact.clone(),
        inertia: exact,
    };
    let p = mass::prismatic_properties(
        &plan,
        &accuracy,
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(p.volume.is_exact());
    assert_eq!(p.volume.lo, q("3/4"));
    assert_eq!(
        p.centroid.unwrap().map(|v| {
            assert!(v.is_exact());
            v.lo
        }),
        [q("11/12"), q("5/4"), q("7/10")]
    );
    assert!(p.inertia.unwrap().iter().flatten().all(Enclosure::is_exact));
}
