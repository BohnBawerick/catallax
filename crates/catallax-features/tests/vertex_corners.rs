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
