//! Independent analytic controls for finite signed algebraic extrusion offsets.
use catallax_features::evaluate;
use catallax_numbers::{Ball, Decider, Outcome, Rational};
use catallax_recipe::{Operation, Recipe};
use catallax_shapes::RPoint3;

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&r(n), 192)
}
fn recipes() -> [Recipe; 3] {
    [
        include_str!("fixtures/algebraic-extrusions/rectangle.json"),
        include_str!("fixtures/algebraic-extrusions/circle.json"),
        include_str!("fixtures/algebraic-extrusions/slot.json"),
    ]
    .map(|s| Recipe::parse(s).unwrap())
}
fn signed(mut recipe: Recipe, distance: i64) -> Recipe {
    let Operation::Extrude { distance: d, .. } = &mut recipe.steps[1].op else {
        unreachable!()
    };
    *d = distance.into();
    recipe
}
fn encloses(value: &catallax_measure::Enclosure, expected: Ball) {
    assert!(
        value.lo <= expected.lo() && value.hi >= expected.hi(),
        "{value:?} versus {expected:?}"
    );
    assert!(value.width() <= Rational::parse("1e-9").unwrap());
}
#[test]
fn signed_line_circle_and_arc_sweeps_retain_exact_volumes_topology_and_step() {
    for (recipe, area) in
        recipes()
            .into_iter()
            .zip([b(8), b(4) * Ball::pi(192), b(8) + Ball::pi(192)])
    {
        for distance in [5, -5] {
            let recipe = signed(recipe.clone(), distance);
            let built = evaluate(&recipe, &mut Decider::default()).unwrap();
            assert!(built.model.check(built.result).is_empty());
            let mass = catallax_measure::mass_properties(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
            assert!(mass.marks.is_empty());
            encloses(
                &mass.value.volume,
                &area * &(b(5) * b(2).sqrt().unwrap().recip().unwrap()),
            );
            let mut bytes = Vec::new();
            let report = catallax_step::write_step(
                &built.model,
                &[built.result],
                &Default::default(),
                &mut Decider::default(),
                &mut bytes,
            )
            .unwrap();
            assert!(report.marks.is_empty());
            assert!(bytes.starts_with(b"ISO-10303-21;"));
            let mesh = catallax_mesh::mesh_body(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
            assert!(mesh.marks.is_empty());
            assert!(!mesh.value.triangles.is_empty());
        }
    }
}
#[test]
fn signed_rectangle_contains_exactly_the_analytic_prism() {
    for distance in [5, -5] {
        let mut built = evaluate(
            &signed(recipes()[0].clone(), distance),
            &mut Decider::default(),
        )
        .unwrap();
        for x in [-1, 1, 3, 5, 8] {
            for y in [-1, 1, 3] {
                for z in [-5, -3, -1, 1, 3, 5] {
                    // At height z, the rectangle has translated by (z,0).
                    // For distance +/-5 its cap is +/-5/sqrt(2), between 3 and 4.
                    let inside = y == 1
                        && x > z
                        && x < z + 4
                        && if distance > 0 {
                            z > 0 && z < 4
                        } else {
                            z < 0 && z > -4
                        };
                    let boundary = y == 1
                        && (x == z || x == z + 4)
                        && if distance > 0 {
                            z > 0 && z < 4
                        } else {
                            z < 0 && z > -4
                        };
                    let actual = catallax_boolean::classify_point(
                        &mut built.model,
                        built.result,
                        &RPoint3::new(r(x), r(y), r(z)),
                        &mut Decider::default(),
                    );
                    if boundary {
                        assert!(
                            matches!(actual, Outcome::ProvenCoincident(_)),
                            "{distance} ({x},{y},{z}) {actual:?}"
                        );
                    } else {
                        let expected = if inside {
                            catallax_boolean::Containment::Inside
                        } else {
                            catallax_boolean::Containment::Outside
                        };
                        assert!(
                            matches!(actual, Outcome::Proven(value) if value == expected),
                            "{distance} ({x},{y},{z}) {actual:?} expected {expected:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn direction_scale_sign_and_a_rigid_profile_pose_preserve_analytic_volume() {
    use catallax_shapes::{RVector3, Rotation};
    let pose = catallax_shapes::Pose::from_rational(
        Rotation::from_quaternion(r(1), r(2), r(3), r(4)).unwrap(),
        RVector3::new(r(7), r(-3), r(11)),
    );
    for components in [[2, 0, 2], [1, 1, 1], [1, 0, -1], [-1, 0, -1]] {
        let direction = RVector3::new(r(components[0]), r(components[1]), r(components[2]));
        let mapped = pose.apply_rational_vector(&direction).unwrap();
        for distance in [3, -3] {
            let mut source = signed(recipes()[0].clone(), distance);
            let Operation::Sketch { plane, .. } = &mut source.steps[0].op else {
                unreachable!()
            };
            plane.rotation = [1, 2, 3, 4].map(Into::into);
            plane.translation = catallax_recipe::Vector3([7, -3, 11].map(Into::into));
            let Operation::Extrude {
                direction: value, ..
            } = &mut source.steps[1].op
            else {
                unreachable!()
            };
            *value = Some(catallax_recipe::Vector3(
                [mapped.x.clone(), mapped.y.clone(), mapped.z.clone()].map(Into::into),
            ));
            let built = evaluate(&source, &mut Decider::default()).unwrap();
            assert!(built.model.check(built.result).is_empty());
            let volume = catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
            let expected = b(24 * components[2].abs())
                * Ball::from_rational(&direction.norm_squared(), 192)
                    .sqrt()
                    .unwrap()
                    .recip()
                    .unwrap();
            encloses(&volume.value, expected);
        }
    }
}

#[test]
fn zero_direction_and_zero_thickness_keep_their_distinct_refusals() {
    for (direction, code) in [
        ([0, 0, 0], "extrude.zero-direction"),
        ([1, 1, 0], "extrude.zero-thickness"),
    ] {
        let mut recipe = recipes()[0].clone();
        let Operation::Extrude {
            direction: value, ..
        } = &mut recipe.steps[1].op
        else {
            unreachable!()
        };
        *value = Some(catallax_recipe::Vector3(direction.map(Into::into)));
        let error = evaluate(&recipe, &mut Decider::default()).unwrap_err();
        let catallax_features::BuildError::Geometry { refusal, .. } = error else {
            panic!("{error:?}")
        };
        assert_eq!(refusal.code, code);
    }
}

#[test]
fn translated_and_rotated_brep_copies_retain_volume_and_exact_membership() {
    use catallax_shapes::{Isometry, RVector3, Rotation};
    let mut built = evaluate(&signed(recipes()[0].clone(), -5), &mut Decider::default()).unwrap();
    for rotation in [[1, 0, 0, 0], [1, 2, 3, 4]] {
        let [w, x, y, z] = rotation.map(r);
        let pose = catallax_shapes::Pose::from_rational(
            Rotation::from_quaternion(w, x, y, z).unwrap(),
            RVector3::new(r(7), r(-3), r(11)),
        );
        let copy = catallax_topology::transformed_body(
            &mut built.model,
            built.result,
            &Isometry::from(&pose),
            Clone::clone,
        );
        assert!(built.model.check(copy).is_empty());
        let volume = catallax_measure::volume(
            &built.model,
            copy,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        encloses(&volume.value, b(40) * b(2).sqrt().unwrap().recip().unwrap());
        for (p, expected) in [
            ([0, 1, -1], catallax_boolean::Containment::Inside),
            ([0, 1, -4], catallax_boolean::Containment::Outside),
            ([4, 1, -1], catallax_boolean::Containment::Outside),
        ] {
            let [x, y, z] = p.map(r);
            let mapped = pose.apply_rational_point(&RPoint3::new(x, y, z)).unwrap();
            let actual = catallax_boolean::classify_point(
                &mut built.model,
                copy,
                &mapped,
                &mut Decider::default(),
            );
            assert!(
                matches!(actual, Outcome::Proven(value) if value == expected),
                "{rotation:?} {p:?}: {actual:?}"
            );
        }
    }
}

#[test]
fn nested_circle_loops_and_irrational_arc_radii_keep_their_exact_material() {
    use serde_json::json;
    let annulus = json!([[{"kind":"circle","id":"outer","centre":["0","0"],"radius":"2"}],
                        [{"kind":"circle","id":"hole","centre":["0","0"],"radius":"1"}]]);
    // The short arc has radius sqrt(2) and subtends pi/2; subtract its
    // right triangle to obtain profile area pi/2 - 1 independently.
    let segment = json!([[{"kind":"arc","id":"round","centre":["0","0"],"to":["1","1"],"turn":"ccw"},
                         {"kind":"line","id":"chord","to":["1","-1"]}]]);
    for (loops, area) in [
        (annulus, b(3) * Ball::pi(192)),
        (segment, Ball::pi(192) * b(2).recip().unwrap() - b(1)),
    ] {
        for distance in [5, -5] {
            let mut source: serde_json::Value =
                serde_json::from_str(include_str!("fixtures/algebraic-extrusions/rectangle.json"))
                    .unwrap();
            source["steps"][0]["loops"] = loops.clone();
            source["steps"][1]["distance"] = distance.to_string().into();
            let recipe = Recipe::parse(&source.to_string()).unwrap();
            let built = evaluate(&recipe, &mut Decider::default()).unwrap();
            assert!(built.model.check(built.result).is_empty());
            let volume = catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
            encloses(
                &volume.value,
                &area * &(b(5) * b(2).sqrt().unwrap().recip().unwrap()),
            );
        }
    }
}
