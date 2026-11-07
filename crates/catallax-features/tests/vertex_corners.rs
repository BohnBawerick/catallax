//! Fully selected planar vertices with more than three incident faces.

use catallax_boolean::Containment;
use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Curve, RPoint3, Surface};
use catallax_topology::Loop;
use serde_json::{Value, json};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

fn assert_zero(value: &Real, dec: &mut Decider) {
    assert!(
        matches!(dec.sign(value), Outcome::ProvenCoincident(_)),
        "{value:?}"
    );
}

fn pyramid(radius: &str, shift: [i64; 3], height: i64) -> Recipe {
    let row: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == "vertex-4-all")
        .unwrap();
    let mut recipe = row["recipe"].clone();
    for step in recipe["steps"].as_array_mut().unwrap() {
        if step["op"] == "fillet" {
            step["radius"] = json!(radius);
        }
        if step["op"] == "sketch" {
            step["loops"][0][2]["to"][1] = json!(height.to_string());
            for (value, delta) in step["plane"]["translation"]
                .as_array_mut()
                .unwrap()
                .iter_mut()
                .zip(shift)
            {
                let old = value.as_str().unwrap().parse::<i64>().unwrap();
                *value = json!((old + delta).to_string());
            }
        }
    }
    Recipe::parse(&recipe.to_string()).unwrap()
}

#[test]
fn four_contact_corner_is_one_exact_shared_spherical_polygon() {
    let mut previous = None;
    for (radius, shift, height) in [
        ("1/2", [0, 0, 0], 10),
        ("1/4", [3, 5, 7], 10),
        ("1/3", [0, 0, 0], 5),
    ] {
        let mut dec = decider();
        let built = evaluate(&pyramid(radius, shift, height), &mut dec).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let faces = built.model.faces_of_body(built.result);
        let apex = faces
            .iter()
            .copied()
            .find(|face| {
                built.model.face_origins(*face).iter().any(|origin| {
                matches!(&origin.role, Role::Corner { vertex } if vertex.faces().len() == 4)
            })
            })
            .expect("the four-face input vertex has a corner origin");
        let origin = built.model.face_origins(apex).to_vec();
        if let Some(previous) = previous {
            assert_eq!(origin, previous);
        }
        previous = Some(origin);
        let Surface::Sphere(sphere) = built.model.geom().surface(built.model.face(apex).surface)
        else {
            panic!("a common rolling ball must produce an exact sphere")
        };
        let r = Real::from(Rational::parse(radius).unwrap());
        assert_zero(&(&sphere.radius - &r), &mut dec);
        assert_zero(&(&sphere.center.x - &Real::from_i64(shift[0])), &mut dec);
        assert_zero(&(&sphere.center.y - &Real::from_i64(shift[1])), &mut dec);
        let slope = Real::from(Rational::from_i64(height) / Rational::from_i64(10));
        let expected_z =
            Real::from_i64(height + shift[2]) - &r * &(Real::one() + slope.square()).sqrt();
        assert_zero(&(&sphere.center.z - &expected_z), &mut dec);
        let [lp] = built.model.face(apex).loops.as_slice() else {
            panic!("one polygon loop")
        };
        let Loop::Coedges(cs) = built.model.loop_(*lp) else {
            panic!("corner has boundary arcs")
        };
        assert_eq!(cs.len(), 4);
        for c in cs {
            let edge = built.model.coedge(*c).edge;
            let Curve::Circle(circle) = built.model.geom().curve(built.model.edge(edge).curve)
            else {
                panic!("spherical hull boundaries are great-circle arcs")
            };
            assert_zero(
                &circle.center.minus(&sphere.center).norm_squared(),
                &mut dec,
            );
            assert_zero(&(&circle.radius_sq - &r.square()), &mut dec);
            let neighbours: Vec<_> = faces
                .iter()
                .copied()
                .filter(|face| {
                    built.model.face(*face).loops.iter().any(|lp| {
                        matches!(built.model.loop_(*lp), Loop::Coedges(cs)
                        if cs.iter().any(|c| built.model.coedge(*c).edge == edge))
                    })
                })
                .collect();
            assert_eq!(neighbours.len(), 2);
            let strip = neighbours.into_iter().find(|face| *face != apex).unwrap();
            assert!(matches!(
                built.model.geom().surface(built.model.face(strip).surface),
                Surface::Cylinder(_)
            ));
        }
    }
}

#[test]
fn four_contact_corner_removes_the_apex_without_cutting_the_core() {
    let mut dec = decider();
    let mut built = evaluate(&pyramid("1/2", [0, 0, 0], 10), &mut dec).unwrap();
    for (coordinates, expected) in [
        (["0", "0", "97/10"], Containment::Inside),
        (["0", "0", "99/10"], Containment::Outside),
        (["5", "5", "5"], Containment::Outside),
        (["0", "0", "5"], Containment::Inside),
        (["0", "0", "-1/10"], Containment::Outside),
    ] {
        let [x, y, z] = coordinates.map(|v| Rational::parse(v).unwrap());
        let p = RPoint3::new(x, y, z);
        assert!(
            matches!(catallax_boolean::classify_point(
            &mut built.model, built.result, &p, &mut dec,
        ), Outcome::Proven(actual) if actual == expected),
            "at {p:?}"
        );
    }
}
