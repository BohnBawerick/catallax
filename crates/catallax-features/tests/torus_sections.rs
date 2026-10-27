//! Finite torus trims must survive sections of the unchanged rounded boss.

use catallax_boolean::{BoolOp, boolean};
use catallax_features::evaluate;
use catallax_measure::{Accuracy, face_area};
use catallax_numbers::{Ball, BudgetPolicy, Decider, Rational};
use catallax_recipe::{Operation, Recipe};
use catallax_shapes::{Pose, RVector3, Rotation, Surface};
use catallax_topology::block;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

#[test]
fn p3_original_and_radius_edit_have_the_analytic_center_sections() {
    let mut recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    for radius in ["2", "5/2"] {
        let Operation::Fillet { radius: value, .. } = &mut recipe.steps[3].op else {
            panic!("P3 junction must remain a fillet");
        };
        *value = serde_json::from_value(serde_json::json!(radius)).unwrap();
        for axis in [0, 1] {
            let mut built = evaluate(&recipe, &mut decider()).unwrap();
            let coordinate = if axis == 0 { q("20") } else { q("15") };
            let mut offset = [q("-10"), q("-10"), q("-10")];
            offset[axis] = coordinate.clone();
            let [x, y, z] = offset;
            let pose = Pose::new(
                Rotation::identity(),
                RVector3::new(x, y, z).into(),
                &mut decider(),
            )
            .unwrap();
            let clip = block(&mut built.model, [q("60"), q("60"), q("50")], &pose).unwrap();
            let cut = boolean(
                &mut built.model,
                built.result,
                clip,
                BoolOp::Intersect,
                &mut decider(),
            )
            .unwrap_or_else(|error| panic!("radius={radius}, axis={axis}: {error:?}"));
            assert!(built.model.check(cut.body).is_empty());
            let section: Vec<_> = built
                .model
                .faces_of_body(cut.body)
                .into_iter()
                .filter(|face| {
                    let Surface::Plane(plane) =
                        built.model.geom().surface(built.model.face(*face).surface)
                    else {
                        return false;
                    };
                    let origin = if axis == 0 {
                        &plane.origin.x
                    } else {
                        &plane.origin.y
                    };
                    let normal = if axis == 0 {
                        &plane.normal.x
                    } else {
                        &plane.normal.y
                    };
                    origin.as_rational() == Some(coordinate.clone())
                        && normal.as_rational().is_some_and(|v| v != q("0"))
                })
                .collect();
            assert_eq!(
                section.len(),
                1,
                "the section is one connected material region"
            );
            let accuracy = Accuracy {
                relative: q("0"),
                absolute: q("1e-8"),
            };
            let area = face_area(&built.model, section[0], &accuracy, &mut decider()).unwrap();
            assert!(area.marks.is_empty());
            // Width*10 + boss width*12, less the two convex corner deficits,
            // plus the two concave quarter-circle complements.
            let b = |value: Rational| Ball::from_rational(&value, 256);
            let radius2 = q(radius) * q(radius);
            let width = if axis == 0 { q("30") } else { q("40") };
            let expected = b(width * q("10") + q("190") + &q("2") * &radius2)
                + b((q("1") - radius2) / q("2")) * Ball::pi(256);
            assert!(
                area.value.lo <= expected.lo() && expected.hi() <= area.value.hi,
                "radius={radius}, axis={axis}: area={:?}, expected={expected:?}",
                area.value
            );
            assert!(area.value.width() <= q("1e-8"));
            let options = catallax_mesh::MeshOptions {
                chord_tolerance: q("0.005"),
                angular_deviation: Some(Ball::pi(256).lo() / q("90")),
                ..Default::default()
            };
            let display =
                catallax_mesh::mesh_body(&built.model, cut.body, &options, &mut decider())
                    .unwrap_or_else(|error| panic!("radius={radius}, axis={axis}: {error}"));
            assert!(display.marks.is_empty());
            assert!(display.value.unavailable_controls.is_empty());
            for face in &display.value.faces {
                assert!(face.max_chord_error <= options.chord_tolerance);
                assert!(
                    face.max_angular_deviation.as_ref().unwrap()
                        <= options.angular_deviation.as_ref().unwrap()
                );
            }
        }
    }
}

#[test]
fn finite_torus_arcs_do_not_turn_support_extensions_into_boundaries() {
    let recipe = Recipe::parse(
        r#"{
        "format":"catallax-recipe/3","name":"trimmed torus",
        "steps":[
            {"id":"ring","op":"torus","major":"3","minor":"1",
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"clip","op":"block","size":["5","10","5"],
             "pose":{"rotation":["1","0","0","0"],"translation":["0","-5","0"]}},
            {"id":"cut","op":"boolean","kind":"intersect","target":"ring","tool":"clip"}
        ],"result":"cut"
    }"#,
    )
    .unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let faces: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|face| {
            matches!(
                built.model.geom().surface(built.model.face(*face).surface),
                Surface::Torus(_)
            )
        })
        .collect();
    assert_eq!(
        faces.len(),
        1,
        "the azimuth seam must not split the half-ring"
    );
    for (coordinates, expected) in [
        (["18/5", "0", "4/5"], true),
        (["12/5", "0", "4/5"], true),
        (["54/25", "72/25", "4/5"], true),
        (["54/25", "-72/25", "4/5"], true),
        (["-18/5", "0", "4/5"], false),
        (["18/5", "0", "-4/5"], false),
        (["0", "18/5", "-4/5"], false),
        (["0", "-18/5", "-4/5"], false),
        (["-4", "0", "0"], false),
        (["0", "18/5", "4/5"], true),
        (["4", "0", "0"], true),
    ] {
        let [x, y, z] = coordinates.map(q);
        let point = catallax_shapes::RPoint3::new(x, y, z);
        assert_eq!(
            catallax_boolean::face_contains(&mut built.model, faces[0], &point, &mut decider())
                .unwrap(),
            expected,
            "point={coordinates:?}",
        );
    }
}
