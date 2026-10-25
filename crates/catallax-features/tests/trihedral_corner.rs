//! Exact spherical corners where three perpendicular convex strips meet.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::Containment;
use catallax_features::{Built, evaluate};
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Operation, Origin, Recipe, Role};
use catallax_shapes::{Curve, RPoint3, Surface};
use catallax_topology::Loop;
use serde_json::{Value, json};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    })
}
fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn corner(radius: &str) -> Recipe {
    let mut recipe = Recipe::parse(include_str!("fixtures/trihedral-corner.json")).unwrap();
    if let Operation::Fillet { radius: value, .. } = &mut recipe.steps[1].op {
        *value = serde_json::from_value(json!(radius)).unwrap();
    }
    recipe
}
fn rounded_block(radius: &str) -> Recipe {
    Recipe::parse(
        &json!({
            "format": "catallax-recipe/3", "name": "rounded block",
            "steps": [
                {"id":"block","op":"block","size":["10","10","10"],
                 "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                {"id":"round","op":"fillet","body":"block","radius":radius,"edges":"all"}
            ],
            "result": "round"
        })
        .to_string(),
    )
    .unwrap()
}
fn build(recipe: &Recipe) -> Built {
    let built = evaluate(recipe, &mut decider()).unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    built
}
// The corner offset `c + radius * (1/3, 2/3, 2/3)` is a rational sphere point.
fn sphere_point(centre: &str, radius: &str) -> RPoint3 {
    let c = r(centre);
    let radius = r(radius);
    RPoint3::new(
        &c + &(&radius * &r("1/3")),
        &c + &(&radius * &r("2/3")),
        &c + &(&radius * &r("2/3")),
    )
}
fn diagonal(t: &str) -> RPoint3 {
    RPoint3::new(r(t), r(t), r(t))
}
// Points near the 10,10,10 corner of a radius-2 rounding, with the Steinmetz
// point that lies in all three strip cylinders but outside the ball.
fn controls() -> Vec<(RPoint3, Containment)> {
    vec![
        (diagonal("5"), Containment::Inside),
        (diagonal("91/10"), Containment::Inside),
        (diagonal("92/10"), Containment::Outside),
        (diagonal("93/10"), Containment::Outside),
        (diagonal("99/10"), Containment::Outside),
        (
            RPoint3::new(r("93/10"), r("93/10"), r("5")),
            Containment::Inside,
        ),
        (
            RPoint3::new(r("95/10"), r("95/10"), r("5")),
            Containment::Outside,
        ),
    ]
}
fn names(built: &Built) -> Vec<Origin> {
    let mut names: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| built.model.face_origins(face).to_vec())
        .filter(|o| matches!(o.role, Role::Blend { .. } | Role::Corner { .. }))
        .collect();
    names.sort();
    names
}
// 1000 minus the edge strips outside the corner cubes and each cube less its
// ball octant: rational and pi coefficients of the exact volume.
fn volume(edges: i64, length: &Rational, corners: i64, radius: &Rational) -> (Rational, Rational) {
    let r2 = radius * radius;
    let r3 = &r2 * radius;
    let strips = &(&Rational::from_i64(edges) * length) * &r2;
    let cubes = &Rational::from_i64(corners) * &r3;
    (
        &(&Rational::from_i64(1000) - &strips) - &cubes,
        &(&strips / &Rational::from_i64(4)) + &(&cubes / &Rational::from_i64(6)),
    )
}
fn contains_volume(lo: &Rational, hi: &Rational, exact: &(Rational, Rational)) {
    let below = &exact.0 + &(&exact.1 * &r("3.14159265358979"));
    let above = &exact.0 + &(&exact.1 * &r("3.14159265358980"));
    assert!(lo <= &above && hi >= &below, "[{lo}, {hi}]");
    assert!(hi - lo <= r("1e-6"), "[{lo}, {hi}]");
}

#[test]
fn block_corner_builds_a_shared_exact_sphere_sector() {
    for (radius, centre) in [("2", "8"), ("3/2", "17/2")] {
        let built = build(&corner(radius));
        let model = &built.model;
        let faces = model.faces_of_body(built.result);
        assert_eq!(faces.len(), 10);
        let corners: Vec<_> = faces
            .iter()
            .filter(|f| {
                model
                    .face_origins(**f)
                    .iter()
                    .any(|o| matches!(o.role, Role::Corner { .. }))
            })
            .copied()
            .collect();
        let [corner] = corners.as_slice() else {
            panic!("one corner face expected: {corners:?}");
        };
        let origins = model.face_origins(*corner);
        assert_eq!(origins.len(), 1);
        let Role::Corner { vertex } = &origins[0].role else {
            unreachable!()
        };
        let roles: Vec<_> = vertex.faces().iter().map(|o| o.role.clone()).collect();
        assert_eq!(roles.len(), 3);
        assert!(roles.contains(&Role::End));
        let Surface::Sphere(sphere) = model.geom().surface(model.face(*corner).surface) else {
            panic!("corner is not a sphere");
        };
        for value in [&sphere.center.x, &sphere.center.y, &sphere.center.z] {
            assert_eq!(value.as_rational(), Some(r(centre)));
        }
        assert_eq!(sphere.radius.as_rational(), Some(r(radius)));
        // Every corner boundary is an exact quarter circle shared with one strip.
        let [lp] = model.face(*corner).loops.as_slice() else {
            panic!("one corner loop expected");
        };
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            panic!("corner loop has edges");
        };
        assert_eq!(coedges.len(), 3);
        let mut neighbours = BTreeSet::new();
        for c in coedges {
            let edge = model.coedge(*c).edge;
            let Curve::Circle(circle) = model.geom().curve(model.edge(edge).curve) else {
                panic!("corner boundary is not a circle");
            };
            assert_eq!(
                circle.radius_sq.as_rational(),
                Some(&r(radius) * &r(radius))
            );
            let sharing: Vec<_> = faces
                .iter()
                .filter(|f| {
                    model.face(**f).loops.iter().any(|lp| {
                        matches!(model.loop_(*lp), Loop::Coedges(cs)
                            if cs.iter().any(|c| model.coedge(*c).edge == edge))
                    })
                })
                .copied()
                .collect();
            assert_eq!(sharing.len(), 2);
            let strip = *sharing.iter().find(|f| *f != corner).unwrap();
            assert!(matches!(
                model.geom().surface(model.face(strip).surface),
                Surface::Cylinder(_)
            ));
            assert!(
                model
                    .face_origins(strip)
                    .iter()
                    .any(|o| matches!(o.role, Role::Blend { .. }))
            );
            neighbours.insert(strip);
        }
        assert_eq!(neighbours.len(), 3);
    }
}

#[test]
fn corner_names_survive_a_radius_edit_in_both_evaluators() {
    let first = names(&build(&corner("2")));
    assert_eq!(first.len(), 4);
    assert_eq!(first, names(&build(&corner("3/2"))));
    let field: Vec<_> = Plan::compile_checked(&corner("2"), &mut decider())
        .unwrap()
        .faces()
        .iter()
        .map(|f| f.origin.clone())
        .filter(|o| matches!(o.role, Role::Blend { .. } | Role::Corner { .. }))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    assert_eq!(field, first);
}

#[test]
fn corner_material_agrees_with_independent_field_and_analytic_controls() {
    let mut built = build(&corner("2"));
    let plan = Plan::compile_checked(&corner("2"), &mut decider()).unwrap();
    let mut workspace = Workspace::default();
    let mut dec = decider();
    for (p, want) in controls() {
        assert!(
            matches!(
                catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut dec),
                Outcome::Proven(actual) if actual == want
            ),
            "B-rep at {p:?}"
        );
        let field = plan.proven_checked(&p, &mut workspace, &mut dec).unwrap();
        let want = match want {
            Containment::Inside => Location::Inside,
            Containment::Outside => Location::Outside,
        };
        assert_eq!(field.location, want, "field at {p:?}: {field:?}");
    }
    let on = sphere_point("8", "2");
    let corner_face = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| matches!(o.role, Role::Corner { .. }))
        })
        .unwrap();
    assert!(catallax_boolean::face_contains(&mut built.model, corner_face, &on, &mut dec).unwrap());
    let origin = built.model.face_origins(corner_face)[0].clone();
    let field = plan.proven_checked(&on, &mut workspace, &mut dec).unwrap();
    assert_eq!(field.location, Location::Boundary);
    let face = plan.select(&origin).unwrap();
    let query = plan
        .face_query_checked(face, &on, &mut workspace, &mut dec)
        .unwrap();
    assert_eq!(query.state, FaceState::OnFace);
    // A sphere point outside the sector is removed by the strips and the cube.
    let off = RPoint3::new(r("8"), r("8"), r("6"));
    assert_eq!(
        plan.proven_checked(&off, &mut workspace, &mut dec)
            .unwrap()
            .location,
        Location::Inside
    );
}

#[test]
fn corner_volume_display_and_step_are_exact_closed_and_complete() {
    for (radius, edges, corners) in [("2", 3, 1), ("3/2", 3, 1)] {
        let built = build(&corner(radius));
        let radius = r(radius);
        let mass = catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut decider(),
        )
        .unwrap();
        assert!(mass.marks.is_empty(), "{:?}", mass.marks);
        let exact = volume(edges, &(&r("10") - &radius), corners, &radius);
        contains_volume(&mass.value.volume.lo, &mass.value.volume.hi, &exact);
        check_display(&built);
        check_step(&built, 10);
    }
}

#[test]
fn fully_rounded_block_closes_eight_corners_in_both_evaluators() {
    let built = build(&rounded_block("1"));
    let faces = built.model.faces_of_body(built.result);
    assert_eq!(faces.len(), 26);
    let spheres = faces
        .iter()
        .filter(|f| {
            matches!(
                built.model.geom().surface(built.model.face(**f).surface),
                Surface::Sphere(_)
            )
        })
        .count();
    assert_eq!(spheres, 8);
    let corner_names: BTreeSet<_> = names(&built)
        .into_iter()
        .filter(|o| matches!(o.role, Role::Corner { .. }))
        .collect();
    assert_eq!(corner_names.len(), 8);
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap();
    let exact = volume(12, &r("8"), 8, &r("1"));
    contains_volume(&mass.value.volume.lo, &mass.value.volume.hi, &exact);
    check_display(&built);
    check_step(&built, 26);
    let plan = Plan::compile_checked(&rounded_block("1"), &mut decider()).unwrap();
    let field_corners: BTreeSet<_> = plan
        .faces()
        .iter()
        .filter(|f| matches!(f.origin.role, Role::Corner { .. }))
        .map(|f| f.origin.clone())
        .collect();
    assert_eq!(field_corners, corner_names);
    let mut workspace = Workspace::default();
    let mut dec = decider();
    let mut built = built;
    for (p, want) in [
        (diagonal("2/5"), Location::Outside),
        (diagonal("1/2"), Location::Inside),
        (diagonal("48/5"), Location::Outside),
        (diagonal("19/2"), Location::Inside),
        (sphere_point("1", "-1"), Location::Boundary),
        (sphere_point("9", "1"), Location::Boundary),
    ] {
        let query = plan.proven_checked(&p, &mut workspace, &mut dec).unwrap();
        assert_eq!(query.location, want, "field at {p:?}");
        let brep = catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut dec);
        match want {
            Location::Inside => assert!(
                matches!(brep, Outcome::Proven(Containment::Inside)),
                "B-rep at {p:?}"
            ),
            Location::Outside => assert!(
                matches!(brep, Outcome::Proven(Containment::Outside)),
                "B-rep at {p:?}"
            ),
            _ => assert!(!matches!(brep, Outcome::Proven(_)), "B-rep at {p:?}"),
        }
    }
}

fn check_display(built: &Built) {
    let options = catallax_mesh::MeshOptions::default();
    let mesh = catallax_mesh::mesh_body(&built.model, built.result, &options, &mut decider())
        .unwrap()
        .value;
    let expected: BTreeSet<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .collect();
    let actual: BTreeSet<_> = mesh.triangles.iter().map(|t| t.face).collect();
    assert_eq!(actual, expected);
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let uses = incidence.entry((a.min(b), a.max(b))).or_default();
            uses.0 += 1;
            uses.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(incidence.values().all(|uses| *uses == (2, 0)));
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
    }
}

fn check_step(built: &Built, faces: usize) {
    let mut bytes = Vec::new();
    let answer = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    let text = String::from_utf8(bytes).unwrap();
    assert_eq!(text.matches("ADVANCED_FACE(").count(), faces);
    assert!(text.contains("SPHERICAL_SURFACE("));
}

#[test]
fn out_of_domain_corners_are_refused_in_both_evaluators() {
    // A sheared prism meets its end and floor at an oblique side with a
    // rational unit normal, so only the corner itself is out of domain.
    let sheared: Value = json!({
        "format": "catallax-recipe/3", "name": "sheared corner",
        "steps": [
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
             "loops":[[
                {"kind":"line","id":"bottom","to":["10","0"]},
                {"kind":"line","id":"right","to":["13","4"]},
                {"kind":"line","id":"top","to":["3","4"]},
                {"kind":"line","id":"left","to":["0","0"]}
             ]]},
            {"id":"prism","op":"extrude","profile":"profile","distance":"10"},
            {"id":"round","op":"fillet","body":"prism","radius":"1","edges":[
                {"between":[{"step":"prism","role":"end"},{"step":"prism","role":{"side":{"segment":"bottom"}}}]},
                {"between":[{"step":"prism","role":"end"},{"step":"prism","role":{"side":{"segment":"right"}}}]},
                {"between":[{"step":"prism","role":{"side":{"segment":"bottom"}}},{"step":"prism","role":{"side":{"segment":"right"}}}]}
            ]}
        ],
        "result": "round"
    });
    let recipe = Recipe::parse(&sheared.to_string()).unwrap();
    let Err(catallax_features::BuildError::Geometry { refusal, .. }) =
        evaluate(&recipe, &mut decider())
    else {
        panic!("sheared corner built");
    };
    assert_eq!(refusal.code, "fillet.unsupported-corner", "{refusal:?}");
    let refusal = Plan::compile_checked(&recipe, &mut decider()).unwrap_err();
    assert_eq!(refusal.code, "fillet.unsupported-corner", "{refusal:?}");
}
