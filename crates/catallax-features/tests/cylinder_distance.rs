//! Circular cylinder joins retain true distance, physical names and material sense.

use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{Recipe, Role};
use serde_json::json;

fn recipe(distance: &str, rotation: [&str; 4], translation: [&str; 3]) -> Recipe {
    Recipe::parse(
        &json!({
            "format":"catallax-recipe/3", "name":"outward cylinder",
            "steps":[
                {"id":"stock","op":"cylinder","radius":"2","height":"4",
                 "pose":{"rotation":rotation,"translation":translation}},
                {"id":"result","op":"offset","body":"stock","distance":distance}
            ], "result":"result"
        })
        .to_string(),
    )
    .unwrap()
}

#[test]
fn outward_cylinder_has_tangent_torus_sectors_not_an_enlarged_sharp_cylinder() {
    let recipe = recipe("1/5", ["1", "0", "0", "0"], ["0", "0", "0"]);
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (point, want) in [
        (["0", "0", "-1/5"], Location::Boundary),
        (["11/5", "0", "2"], Location::Boundary),
        (["53/25", "0", "-4/25"], Location::Boundary),
        (["53/25", "0", "104/25"], Location::Boundary),
        (["211/100", "0", "-4/25"], Location::Inside),
        (["213/100", "0", "-4/25"], Location::Outside),
        (["0", "0", "2"], Location::Inside),
        (["11/5", "0", "21/5"], Location::Outside),
    ] {
        let [x, y, z] = point.map(|v| Rational::parse(v).unwrap());
        let point = catallax_shapes::RPoint3::new(x, y, z);
        assert_eq!(
            plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                .unwrap()
                .location,
            want
        );
        if want != Location::Boundary {
            let actual = classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default(),
            );
            assert_eq!(
                actual.proven(),
                Some(&if want == Location::Inside {
                    Containment::Inside
                } else {
                    Containment::Outside
                }),
                "{actual:?}"
            );
        }
    }
    let rounds: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| matches!(o.role, Role::Round { .. }))
        })
        .collect();
    assert_eq!(rounds.len(), 2);
    for face in rounds {
        assert!(matches!(
            built.model.geom().surface(built.model.face(face).surface),
            catallax_shapes::Surface::Torus(_)
        ));
    }
}

#[test]
fn posed_cylinder_offset_mass_contains_independent_parallel_body_formula() {
    let recipe = recipe("1/5", ["1", "2", "3", "4"], ["3", "-2", "5"]);
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let answer = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    let b = |n, d| Ball::from_rational(&Rational::from_ratio(n, d).unwrap(), 256);
    let pi = Ball::pi(256);
    // Core, planar slabs, cylindrical wall slab, two quarter-torus joins.
    let volume = &pi * &b(7864, 375) + pi.square() * b(2, 25);
    let area = &pi * &b(644, 25) + pi.square() * b(4, 5);
    assert!(answer.value.volume.lo <= volume.lo() && answer.value.volume.hi >= volume.hi());
    assert!(answer.value.area.lo <= area.lo() && answer.value.area.hi >= area.hi());
}

#[test]
fn horn_and_spindle_edge_tubes_keep_a_typed_refusal() {
    for distance in ["2", "3"] {
        let recipe = recipe(distance, ["1", "0", "0", "0"], ["0", "0", "0"]);
        let catallax_features::BuildError::Geometry { refusal, .. } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected geometric refusal")
        };
        assert_eq!(refusal.code, "offset.spindle-carrier");
        assert_eq!(
            Plan::compile(&recipe).unwrap_err().code,
            "offset.spindle-carrier"
        );
    }
}

#[test]
fn sphere_cylinder_meetings_use_spherical_caps_and_convex_circular_rounds() {
    let recipe = Recipe::parse(include_str!(
        "fixtures/cylinder-distance/sphere-cylinder.json"
    ))
    .unwrap();
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (p, want) in [
        (["0", "0", "27/10"], Location::Boundary),
        (["54/25", "0", "81/50"], Location::Boundary),
        (["11/5", "0", "0"], Location::Boundary),
        (["211/100", "0", "83/50"], Location::Inside),
        (["11/5", "0", "17/10"], Location::Outside),
        (["0", "0", "269/100"], Location::Inside),
        (["0", "0", "271/100"], Location::Outside),
    ] {
        let [x, y, z] = p.map(|v| Rational::parse(v).unwrap());
        let point = catallax_shapes::RPoint3::new(x, y, z);
        assert_eq!(
            plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                .unwrap()
                .location,
            want,
            "at {p:?}"
        );
        if want != Location::Boundary {
            let actual = classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default(),
            );
            assert_eq!(
                actual.proven(),
                Some(&if want == Location::Inside {
                    Containment::Inside
                } else {
                    Containment::Outside
                }),
                "{actual:?}"
            );
        }
    }
}

#[test]
fn coaxial_intersection_offsets_keep_only_the_active_wall_and_caps() {
    let recipe = Recipe::parse(include_str!(
        "fixtures/cylinder-distance/coaxial-cylinders.json"
    ))
    .unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (p, want) in [
        (["0", "0", "11/5"], Location::Boundary),
        (["11/5", "0", "0"], Location::Boundary),
        (["5/2", "0", "0"], Location::Outside),
        (["0", "0", "3"], Location::Outside),
    ] {
        let [x, y, z] = p.map(|v| Rational::parse(v).unwrap());
        assert_eq!(
            plan.proven_checked(
                &catallax_shapes::RPoint3::new(x, y, z),
                &mut workspace,
                &mut Decider::default()
            )
            .unwrap()
            .location,
            want
        );
    }
}

#[test]
fn circular_shell_and_single_face_sweeps_retain_their_finite_normal_domains() {
    for (text, inside, outside, rim) in [
        (
            include_str!("fixtures/cylinder-distance/shell-closed.json"),
            ["19/10", "0", "0"],
            ["0", "0", "0"],
            ["9/5", "0", "0"],
        ),
        (
            include_str!("fixtures/cylinder-distance/wall-front.json"),
            ["21/10", "0", "0"],
            ["21/10", "0", "8/5"],
            ["21/10", "0", "3/2"],
        ),
        (
            include_str!("fixtures/cylinder-distance/wall-back.json"),
            ["19/10", "0", "0"],
            ["19/10", "0", "8/5"],
            ["19/10", "0", "3/2"],
        ),
        (
            include_str!("fixtures/cylinder-distance/sphere-front.json"),
            ["0", "0", "13/5"],
            ["21/10", "0", "3/2"],
            ["52/25", "0", "39/25"],
        ),
        (
            include_str!("fixtures/cylinder-distance/sphere-back.json"),
            ["0", "0", "12/5"],
            ["39/20", "0", "7/5"],
            ["48/25", "0", "36/25"],
        ),
    ] {
        let recipe = Recipe::parse(text).unwrap();
        let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert!(
            built.model.check(built.result).is_empty(),
            "{}",
            recipe.name
        );
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        for (p, want) in [
            (inside, Location::Inside),
            (outside, Location::Outside),
            (rim, Location::Boundary),
        ] {
            let [x, y, z] = p.map(|v| Rational::parse(v).unwrap());
            let point = catallax_shapes::RPoint3::new(x, y, z);
            assert_eq!(
                plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                    .unwrap()
                    .location,
                want,
                "{} at {p:?}",
                recipe.name
            );
            if want != Location::Boundary {
                let actual = classify_point(
                    &mut built.model,
                    built.result,
                    &point,
                    &mut Decider::default(),
                );
                assert_eq!(
                    actual.proven(),
                    Some(&if want == Location::Inside {
                        Containment::Inside
                    } else {
                        Containment::Outside
                    }),
                    "{}: {actual:?}",
                    recipe.name
                );
            }
        }
    }
}
