//! Planar multi-face vertices and partial selected-strip ends.

use catallax_boolean::Containment;
use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Curve, RPoint3, Surface};
use catallax_topology::Loop;
use serde_json::{Value, json};
use std::collections::BTreeMap;

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

fn pyramid(case: &str, radius: &str, shift: [i64; 3], height: i64, wall: &str) -> Recipe {
    let row: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == case)
        .unwrap();
    let mut recipe = row["recipe"].clone();
    for step in recipe["steps"].as_array_mut().unwrap() {
        if step["op"] == "fillet" {
            step["radius"] = json!(radius);
            if let Some(edges) = step["edges"].as_array_mut() {
                for edge in edges {
                    for face in edge["between"].as_array_mut().unwrap() {
                        if face["step"] == "b" {
                            face["role"]["side"]["segment"] = json!(wall);
                        }
                    }
                }
            }
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
        let built = evaluate(
            &pyramid("vertex-4-all", radius, shift, height, "right"),
            &mut dec,
        )
        .unwrap();
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
    let mut built = evaluate(
        &pyramid("vertex-4-all", "1/2", [0, 0, 0], 10, "right"),
        &mut dec,
    )
    .unwrap();
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

#[test]
fn one_edge_at_a_four_face_vertex_has_two_exact_end_wall_arcs() {
    for (radius, shift, height) in [
        ("1/2", [0, 0, 0], 10),
        ("1/4", [3, 5, 7], 10),
        ("1/3", [0, 0, 0], 5),
    ] {
        let mut dec = decider();
        let mut built = evaluate(
            &pyramid("vertex-4-one", radius, shift, height, "right"),
            &mut dec,
        )
        .unwrap();
        assert!(built.model.check(built.result).is_empty());
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 6);
        let strips: Vec<_> = faces
            .iter()
            .copied()
            .filter(|face| {
                built
                    .model
                    .face_origins(*face)
                    .iter()
                    .any(|origin| matches!(origin.role, Role::Blend { .. }))
            })
            .collect();
        let [strip] = strips.as_slice() else {
            panic!("only the selected edge has a blend face")
        };
        assert!(matches!(
            built.model.geom().surface(built.model.face(*strip).surface),
            Surface::Cylinder(_)
        ));
        assert!(faces.iter().all(|face| {
            !built
                .model
                .face_origins(*face)
                .iter()
                .any(|origin| matches!(origin.role, Role::Corner { .. }))
        }));
        let [lp] = built.model.face(*strip).loops.as_slice() else {
            panic!("one strip loop")
        };
        let Loop::Coedges(cs) = built.model.loop_(*lp) else {
            panic!("strip loop has edges")
        };
        assert_eq!(cs.len(), 5);
        let mut ellipses = 0;
        for coedge in cs {
            let edge = built.model.coedge(*coedge).edge;
            match built.model.geom().curve(built.model.edge(edge).curve) {
                Curve::Ellipse(_) => ellipses += 1,
                Curve::Line(_) => {}
                curve => panic!("unexpected end curve: {curve:?}"),
            }
            let uses = faces
                .iter()
                .filter(|face| {
                    built.model.face(**face).loops.iter().any(|lp| {
                        matches!(built.model.loop_(*lp), Loop::Coedges(cs)
                    if cs.iter().any(|c| built.model.coedge(*c).edge == edge))
                    })
                })
                .count();
            assert_eq!(uses, 2);
        }
        assert_eq!(ellipses, 3);
        if radius == "1/2" {
            for (coordinates, expected) in [
                (["0", "0", "99/10"], Containment::Inside),
                (["0", "0", "199/20"], Containment::Outside),
                (["5", "5", "5"], Containment::Outside),
                (["49/10", "49/10", "5"], Containment::Inside),
                (["499/100", "-499/100", "5"], Containment::Inside),
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
    }
}

#[test]
fn partial_four_face_corner_keeps_a_sphere_a_cliff_and_both_lower_strips() {
    for (radius, shift, height, wall) in [
        ("1/2", [0, 0, 0], 10, "right"),
        ("1/2", [0, 0, 0], 10, "left"),
        ("1/4", [3, 5, 7], 10, "right"),
        ("1/3", [0, 0, 0], 5, "right"),
    ] {
        let mut dec = decider();
        let mut built = evaluate(
            &pyramid("vertex-4-two", radius, shift, height, wall),
            &mut dec,
        )
        .unwrap();
        assert!(built.model.check(built.result).is_empty());
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 11);
        let strips: Vec<_> = faces
            .iter()
            .copied()
            .filter(|face| {
                built
                    .model
                    .face_origins(*face)
                    .iter()
                    .any(|origin| matches!(origin.role, Role::Blend { .. }))
            })
            .collect();
        assert_eq!(strips.len(), 4);
        let mut by_origin = BTreeMap::new();
        for face in &strips {
            *by_origin
                .entry(built.model.face_origins(*face).to_vec())
                .or_insert(0) += 1;
        }
        assert_eq!(by_origin.len(), 2);
        assert!(by_origin.values().all(|count| *count == 2));
        let edges: Vec<_> = strips
            .iter()
            .copied()
            .map(|face| {
                assert!(matches!(
                    built.model.geom().surface(built.model.face(face).surface),
                    Surface::Cylinder(_)
                ));
                let [lp] = built.model.face(face).loops.as_slice() else {
                    panic!("one strip loop")
                };
                let Loop::Coedges(cs) = built.model.loop_(*lp) else {
                    panic!("strip has an edged boundary")
                };
                assert_eq!(cs.len(), 4);
                cs.iter()
                    .map(|c| built.model.coedge(*c).edge)
                    .collect::<Vec<_>>()
            })
            .collect();
        let mut counts = BTreeMap::new();
        for edge in edges.iter().flatten() {
            *counts.entry(*edge).or_insert(0) += 1;
        }
        let shared: Vec<_> = counts
            .into_iter()
            .filter_map(|(edge, count)| (count == 2).then_some(edge))
            .collect();
        let [shared] = shared.as_slice() else {
            panic!("one exact strip-strip junction")
        };
        assert!(matches!(
            built.model.geom().curve(built.model.edge(*shared).curve),
            Curve::Ellipse(_)
        ));
        for edge in edges.iter().flatten() {
            let uses = faces
                .iter()
                .filter(|face| {
                    built.model.face(**face).loops.iter().any(|lp| {
                        matches!(built.model.loop_(*lp), Loop::Coedges(cs)
                    if cs.iter().any(|c| built.model.coedge(*c).edge == *edge))
                    })
                })
                .count();
            assert_eq!(uses, 2);
        }
        let corner_faces: Vec<_> = faces
            .iter()
            .copied()
            .filter(|face| {
                built.model.face_origins(*face).iter().any(|origin| {
                matches!(&origin.role, Role::Corner { vertex } if vertex.faces().len() == 4)
            })
            })
            .collect();
        assert_eq!(corner_faces.len(), 2);
        assert_eq!(
            built.model.face_origins(corner_faces[0]),
            built.model.face_origins(corner_faces[1])
        );
        let apex = *corner_faces
            .iter()
            .find(|face| {
                matches!(
                    built.model.geom().surface(built.model.face(**face).surface),
                    Surface::Sphere(_)
                )
            })
            .expect("the selected-normal hull is a sphere");
        let cliff = *corner_faces
            .iter()
            .find(|face| {
                matches!(
                    built.model.geom().surface(built.model.face(**face).surface),
                    Surface::Plane(_)
                )
            })
            .expect("the unselected cone boundary is a plane");
        let Surface::Sphere(sphere) = built.model.geom().surface(built.model.face(apex).surface)
        else {
            unreachable!()
        };
        let r = Real::from(Rational::parse(radius).unwrap());
        assert_zero(&(&sphere.radius - &r), &mut dec);
        assert_zero(&(&sphere.center.x - &Real::from_i64(shift[0])), &mut dec);
        assert_zero(&(&sphere.center.y - &Real::from_i64(shift[1])), &mut dec);
        let slope = Real::from(Rational::from_i64(height) / Rational::from_i64(10));
        let expected_z =
            Real::from_i64(height + shift[2]) - &r * &(Real::one() + slope.square()).sqrt();
        assert_zero(&(&sphere.center.z - &expected_z), &mut dec);
        let Surface::Plane(plane) = built.model.geom().surface(built.model.face(cliff).surface)
        else {
            unreachable!()
        };
        assert_zero(&plane.normal.x, &mut dec);
        assert_zero(&plane.normal.z, &mut dec);
        assert_zero(
            &plane.origin.minus(&sphere.center).dot(&plane.normal),
            &mut dec,
        );
        let [lp] = built.model.face(apex).loops.as_slice() else {
            panic!("one spherical triangle")
        };
        let Loop::Coedges(cs) = built.model.loop_(*lp) else {
            panic!("the sphere has great-circle boundaries")
        };
        assert_eq!(cs.len(), 3);
        for c in cs {
            let edge = built.model.coedge(*c).edge;
            let Curve::Circle(circle) = built.model.geom().curve(built.model.edge(edge).curve)
            else {
                panic!("a spherical hull boundary is a great circle")
            };
            assert_zero(
                &circle.center.minus(&sphere.center).norm_squared(),
                &mut dec,
            );
            assert_zero(&(&circle.radius_sq - &r.square()), &mut dec);
        }
        if radius == "1/2" {
            let end = built.model.edge(*shared).end.unwrap();
            let hub = built.model.geom().position(built.model.vertex(end).point);
            let r = Real::from(Rational::parse(radius).unwrap());
            let y = &r
                * &(Real::one() / Real::from_i64(3).sqrt()
                    - Real::one() / Real::from_i64(2).sqrt());
            let sign = if wall == "right" { 1 } else { -1 };
            let expected_y = y.clone() * Real::from_i64(sign);
            assert_zero(&hub.x, &mut dec);
            assert_zero(&(&hub.y - &expected_y), &mut dec);
            assert_zero(&(&hub.z - &(Real::from_i64(10) + y)), &mut dec);
            for (coordinates, expected) in [
                (["0", "0", "199/20"], Containment::Outside),
                (["0", "1/10", "49/5"], Containment::Outside),
                (["0", "-1/10", "49/5"], Containment::Inside),
                (["0", "1/10", "97/10"], Containment::Inside),
                (["5", "5", "5"], Containment::Outside),
                (["-5", "5", "5"], Containment::Outside),
                (["49/10", "49/10", "5"], Containment::Inside),
                (["-49/10", "49/10", "5"], Containment::Inside),
                (["499/100", "-499/100", "5"], Containment::Inside),
                (["-499/100", "-499/100", "5"], Containment::Inside),
                (["0", "0", "5"], Containment::Inside),
                (["0", "0", "-1/10"], Containment::Outside),
            ] {
                let [x, y, z] = coordinates.map(|v| Rational::parse(v).unwrap());
                let p = RPoint3::new(x, y * Rational::from_i64(sign), z);
                let mut point_dec = decider();
                assert!(
                    matches!(catallax_boolean::classify_point(
                    &mut built.model, built.result, &p, &mut point_dec,
                ), Outcome::Proven(actual) if actual == expected),
                    "at {p:?}"
                );
            }
            let p = RPoint3::new(
                Rational::zero(),
                Rational::zero(),
                Rational::parse("49/5").unwrap(),
            );
            assert!(matches!(
                catallax_boolean::classify_point(
                    &mut built.model,
                    built.result,
                    &p,
                    &mut decider(),
                ),
                Outcome::ProvenCoincident(_)
            ));
        }
    }
}

#[test]
fn smooth_partial_corner_keeps_full_radius_and_three_tangent_spherical_joins() {
    for (radius, shift, height, wall, moving, extra) in [
        ("1/2", [0, 0, 0], 10, "right", false, false),
        ("1/2", [0, 0, 0], 10, "left", false, false),
        ("1/4", [3, 5, 7], 10, "right", false, false),
        ("1/3", [0, 0, 0], 5, "right", false, false),
        ("1/2", [0, 0, 0], 10, "right", true, false),
        ("1/2", [0, 0, 0], 10, "right", false, true),
    ] {
        let mut value =
            serde_json::to_value(pyramid("vertex-4-two", radius, shift, height, wall)).unwrap();
        value["steps"].as_array_mut().unwrap().last_mut().unwrap()["corner"] = json!("smooth");
        if moving {
            for step in value["steps"].as_array_mut().unwrap() {
                if step["id"] == "a" {
                    step["distance"] = json!("24");
                }
                if step["id"] == "bs" {
                    step["loops"][0][0]["to"] = json!(["-14", "0"]);
                    step["loops"][0][2]["to"] = json!(["-2", "10"]);
                }
            }
        }
        if extra {
            value["steps"].as_array_mut().unwrap().last_mut().unwrap()["edges"]
                .as_array_mut()
                .unwrap()
                .push(json!({"between": [
                    {"step": "b", "role": {"side": {"segment": "left"}}},
                    {"step": "a", "role": {"side": {"segment": "base"}}}
                ]}));
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let mut dec = decider();
        let mut built = evaluate(&recipe, &mut dec).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), if extra { 10 } else { 9 });
        let cylinders: Vec<_> = faces
            .iter()
            .copied()
            .filter(|face| {
                matches!(
                    built.model.geom().surface(built.model.face(*face).surface),
                    Surface::Cylinder(_)
                )
            })
            .collect();
        assert_eq!(cylinders.len(), if extra { 4 } else { 3 });
        let apex = *faces
            .iter()
            .find(|face| {
                matches!(
                    built.model.geom().surface(built.model.face(**face).surface),
                    Surface::Sphere(_)
                )
            })
            .unwrap();
        let origin = built.model.face_origins(apex).to_vec();
        assert!(matches!(&origin[0].role, Role::Corner { vertex } if vertex.faces().len() == 4));
        assert!(
            cylinders
                .iter()
                .any(|face| built.model.face_origins(*face) == origin.as_slice())
        );
        let Surface::Sphere(sphere) = built.model.geom().surface(built.model.face(apex).surface)
        else {
            unreachable!()
        };
        let r = Real::from(Rational::parse(radius).unwrap());
        assert_zero(&(&sphere.radius - &r), &mut dec);
        if moving {
            let expected = Real::from_i64(-2)
                + &r * &((Real::from_i64(6) * Real::from_i64(2).sqrt()
                    - Real::from_i64(61).sqrt())
                    / Real::from_i64(5));
            assert_zero(&(&sphere.center.y - &expected), &mut dec);
        }
        let [lp] = built.model.face(apex).loops.as_slice() else {
            panic!("one spherical triangle")
        };
        let Loop::Coedges(coedges) = built.model.loop_(*lp) else {
            panic!("three great-circle joins")
        };
        assert_eq!(coedges.len(), 3);
        for coedge in coedges {
            let coedge = built.model.coedge(*coedge);
            let Curve::Circle(circle) = built
                .model
                .geom()
                .curve(built.model.edge(coedge.edge).curve)
            else {
                panic!("the whole sphere/cylinder join is a great circle")
            };
            assert_zero(
                &circle.center.minus(&sphere.center).norm_squared(),
                &mut dec,
            );
            assert_zero(&(&circle.radius_sq - &r.square()), &mut dec);
            let other = cylinders
                .iter()
                .find_map(|face| {
                    built.model.face(*face).loops.iter().find_map(|lp| {
                        let Loop::Coedges(edges) = built.model.loop_(*lp) else {
                            return None;
                        };
                        edges
                            .iter()
                            .find(|edge| built.model.coedge(**edge).edge == coedge.edge)
                            .map(|edge| (*face, built.model.coedge(*edge).sense))
                    })
                })
                .expect("each sphere edge is shared with a cylinder");
            assert_ne!(other.1, coedge.sense);
            let Surface::Cylinder(cylinder) = built
                .model
                .geom()
                .surface(built.model.face(other.0).surface)
            else {
                unreachable!()
            };
            assert_zero(&(&cylinder.radius() - &r), &mut dec);
            assert_zero(
                &sphere
                    .center
                    .minus(&cylinder.origin)
                    .cross(&cylinder.axis)
                    .norm_squared(),
                &mut dec,
            );
            assert_zero(
                &circle.normal.cross(&cylinder.axis).norm_squared(),
                &mut dec,
            );
        }
        if radius == "1/2" {
            let mirror = Rational::from_i64(if wall == "right" { 1 } else { -1 });
            for (coordinates, expected) in [
                (["0", "1/100", "99/10"], Containment::Outside),
                (["0", "-1/100", "99/10"], Containment::Outside),
                (["0", "-1/5", "1959/200"], Containment::Outside),
                (["0", "-1/5", "979/100"], Containment::Inside),
                (["1/3", "-1/3", "29/3"], Containment::Outside),
                (["0", "0", "5"], Containment::Inside),
            ] {
                let [x, y, z] = coordinates.map(|v| Rational::parse(v).unwrap());
                let p = RPoint3::new(
                    x,
                    &y * &mirror + Rational::from_i64(if moving { -2 } else { 0 }),
                    z,
                );
                assert!(
                    matches!(catallax_boolean::classify_point(
                    &mut built.model, built.result, &p, &mut decider()
                ), Outcome::Proven(actual) if actual == expected),
                    "at {p:?}"
                );
            }
            if extra {
                let p = RPoint3::new(Rational::zero(), Rational::from_i64(-10), Rational::zero());
                assert!(matches!(
                    catallax_boolean::classify_point(
                        &mut built.model,
                        built.result,
                        &p,
                        &mut decider()
                    ),
                    Outcome::Proven(Containment::Outside)
                ));
                let p = RPoint3::new(Rational::zero(), Rational::from_i64(-9), Rational::one());
                assert!(matches!(
                    catallax_boolean::classify_point(
                        &mut built.model,
                        built.result,
                        &p,
                        &mut decider()
                    ),
                    Outcome::ProvenCoincident(_)
                ));
            }
        }
    }
}

fn network_recipe(radius: &str, lower: &str, rotated: bool, all: bool, shift: [i64; 3]) -> Recipe {
    let row: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == "vertex-4-all-unequal-slopes")
        .unwrap();
    let mut value = row["recipe"].clone();
    for step in value["steps"].as_array_mut().unwrap() {
        if step["id"] == "bs" {
            step["loops"][0][0]["to"][0] = json!(lower);
        }
        if step["id"] == "a" && lower == "-15" {
            step["distance"] = json!("40");
        }
        if step["op"] == "sketch" {
            if rotated {
                if step["id"] == "as" {
                    step["plane"]["rotation"] = json!(["1", "1", "1", "-1"]);
                } else {
                    step["plane"]["rotation"] = json!(["0", "1", "1", "0"]);
                    step["plane"]["translation"] = json!(["0", "0", "10"]);
                }
            }
            for (coordinate, delta) in step["plane"]["translation"]
                .as_array_mut()
                .unwrap()
                .iter_mut()
                .zip(shift)
            {
                *coordinate = json!(
                    (coordinate.as_str().unwrap().parse::<i64>().unwrap() + delta).to_string()
                );
            }
        }
        if step["op"] == "fillet" {
            step["radius"] = json!(radius);
            if all {
                step["edges"] = json!("all");
            }
        }
    }
    Recipe::parse(&value.to_string()).unwrap()
}

#[test]
fn nonconcurrent_four_face_offsets_form_two_balls_and_one_internal_sheet() {
    let mut previous_origin = None;
    for (radius, lower, rotated, all, shift) in [
        ("1/2", "-5", false, false, [0, 0, 0]),
        ("1/4", "-5", true, false, [3, 5, 7]),
        ("1/2", "-15", false, false, [0, 0, 0]),
        ("1/2", "-5", false, true, [0, 0, 0]),
    ] {
        let mut dec = decider();
        let mut built = evaluate(
            &network_recipe(radius, lower, rotated, all, shift),
            &mut dec,
        )
        .unwrap();
        assert!(built.model.check(built.result).is_empty());
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), if all { 20 } else { 12 });
        let junction: Vec<_> = faces
            .iter()
            .copied()
            .filter(|face| {
                built.model.face_origins(*face).iter().any(|origin| {
                matches!(&origin.role, Role::Corner { vertex } if vertex.faces().len() == 4)
            })
            })
            .collect();
        assert_eq!(junction.len(), 3);
        let origin = built.model.face_origins(junction[0]).to_vec();
        if let Some(previous) = previous_origin {
            assert_eq!(origin, previous);
        }
        for face in &junction {
            assert_eq!(built.model.face_origins(*face), origin.as_slice());
        }
        previous_origin = Some(origin);
        let spheres: Vec<_> = junction
            .iter()
            .copied()
            .filter(|face| {
                matches!(
                    built.model.geom().surface(built.model.face(*face).surface),
                    Surface::Sphere(_)
                )
            })
            .collect();
        assert_eq!(spheres.len(), 2);
        let r = Real::from(Rational::parse(radius).unwrap());
        let sqrt2 = Real::from_i64(2).sqrt();
        let mut expected = if lower == "-5" {
            let sqrt5 = Real::from_i64(5).sqrt();
            let t = r.clone() * (sqrt5.clone() - sqrt2.clone()) / Real::from_i64(3);
            let z = Real::from_i64(10)
                - r.clone() * (Real::from_i64(2) * sqrt2 + sqrt5) / Real::from_i64(3);
            [
                catallax_shapes::Point3::new(t.clone(), t.clone(), z.clone()),
                catallax_shapes::Point3::new(-t.clone(), t, z),
            ]
        } else {
            let y = r.clone() * (Real::from_i64(13).sqrt() - Real::from_i64(3) * sqrt2.clone())
                / Real::from_i64(2);
            let z = Real::from_i64(10) - r.clone() * sqrt2;
            [
                catallax_shapes::Point3::new(Real::zero(), Real::zero(), z.clone()),
                catallax_shapes::Point3::new(Real::zero(), y, z),
            ]
        };
        for point in &mut expected {
            let [x, y, z] = if rotated {
                [point.z.clone(), point.y.clone(), -point.x.clone()]
            } else {
                [point.x.clone(), point.y.clone(), point.z.clone()]
            };
            *point = catallax_shapes::Point3::new(
                x + Real::from_i64(shift[0]),
                y + Real::from_i64(shift[1]),
                z + Real::from_i64(shift[2]),
            );
        }
        let mut matched = [false; 2];
        for face in spheres {
            let Surface::Sphere(sphere) =
                built.model.geom().surface(built.model.face(face).surface)
            else {
                unreachable!()
            };
            assert_zero(&(&sphere.radius - &r), &mut dec);
            let index = expected
                .iter()
                .position(|point| {
                    matches!(
                        dec.sign(&sphere.center.minus(point).norm_squared()),
                        Outcome::ProvenCoincident(_)
                    )
                })
                .expect("each ball has its independently solved triple-contact centre");
            matched[index] = true;
            let [boundary] = built.model.face(face).loops.as_slice() else {
                panic!("one triangular spherical sector")
            };
            let Loop::Coedges(coedges) = built.model.loop_(*boundary) else {
                panic!("three shared great-circle arcs")
            };
            assert_eq!(coedges.len(), 3);
            for coedge in coedges {
                let coedge = built.model.coedge(*coedge);
                let Curve::Circle(circle) = built
                    .model
                    .geom()
                    .curve(built.model.edge(coedge.edge).curve)
                else {
                    panic!("a triple-contact sphere must meet its strips on exact circles")
                };
                assert_zero(
                    &circle.center.minus(&sphere.center).norm_squared(),
                    &mut dec,
                );
                assert_zero(&(&circle.radius_sq - &r.square()), &mut dec);
                let neighbour = faces
                    .iter()
                    .find_map(|other| {
                        if *other == face {
                            return None;
                        }
                        built.model.face(*other).loops.iter().find_map(|boundary| {
                            let Loop::Coedges(edges) = built.model.loop_(*boundary) else {
                                return None;
                            };
                            edges
                                .iter()
                                .find(|edge| built.model.coedge(**edge).edge == coedge.edge)
                                .map(|edge| (*other, built.model.coedge(*edge).sense))
                        })
                    })
                    .expect("every spherical cap has one cylindrical neighbour");
                assert_ne!(neighbour.1, coedge.sense);
                let Surface::Cylinder(cylinder) = built
                    .model
                    .geom()
                    .surface(built.model.face(neighbour.0).surface)
                else {
                    panic!("sphere edges join physical strips or the internal contact sheet")
                };
                assert_zero(&(&cylinder.radius() - &r), &mut dec);
                assert_zero(
                    &sphere
                        .center
                        .minus(&cylinder.origin)
                        .cross(&cylinder.axis)
                        .norm_squared(),
                    &mut dec,
                );
                assert_zero(
                    &circle.normal.cross(&cylinder.axis).norm_squared(),
                    &mut dec,
                );
            }
        }
        assert_eq!(matched, [true, true]);
        let tube = junction
            .iter()
            .find_map(|face| {
                if let Surface::Cylinder(cylinder) =
                    built.model.geom().surface(built.model.face(*face).surface)
                {
                    Some((*face, cylinder))
                } else {
                    None
                }
            })
            .expect("the two contact balls share one finite internal cylinder");
        let [boundary] = built.model.face(tube.0).loops.as_slice() else {
            panic!("one sheet loop")
        };
        let Loop::Coedges(edges) = built.model.loop_(*boundary) else {
            panic!("two rails and two shared spherical caps")
        };
        assert_eq!(edges.len(), 4);
        let mut rails = 0;
        for edge in edges {
            let edge = built.model.coedge(*edge).edge;
            let Curve::Line(line) = built.model.geom().curve(built.model.edge(edge).curve) else {
                continue;
            };
            rails += 1;
            let plane = faces
                .iter()
                .find_map(|face| {
                    let Surface::Plane(plane) =
                        built.model.geom().surface(built.model.face(*face).surface)
                    else {
                        return None;
                    };
                    built
                        .model
                        .face(*face)
                        .loops
                        .iter()
                        .any(|boundary| {
                            matches!(built.model.loop_(*boundary), Loop::Coedges(edges)
                        if edges.iter().any(|own| built.model.coedge(*own).edge == edge))
                        })
                        .then_some(plane)
                })
                .expect("each internal-sheet rail retains its tangent input face");
            assert_zero(&line.dir.cross(&tube.1.axis).norm_squared(), &mut dec);
            assert_zero(
                &line.origin.minus(&plane.origin).dot(&plane.normal),
                &mut dec,
            );
            let depth = tube.1.origin.minus(&plane.origin).dot(&plane.normal);
            assert_zero(
                &(depth.square() - r.square() * plane.normal.norm_squared()),
                &mut dec,
            );
        }
        assert_eq!(rails, 2);
        if radius == "1/2" && lower == "-5" {
            for (p, expected) in [
                (["0", "0.14", "9.65"], Containment::Inside),
                (["0", "0.14", "9.66"], Containment::Outside),
                (["0.31", "0.09", "9.63"], Containment::Outside),
                (["-0.31", "0.09", "9.63"], Containment::Outside),
                (["0", "0", "5"], Containment::Inside),
            ] {
                let [x, y, z] = p.map(|value| Rational::parse(value).unwrap());
                let p = RPoint3::new(x, y, z);
                assert!(
                    matches!(
                        catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut decider()),
                        Outcome::Proven(actual) if actual == expected
                    ),
                    "at {p:?}"
                );
            }
        }
    }
}
