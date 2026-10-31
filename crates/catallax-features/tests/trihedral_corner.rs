//! Exact spherical corners where three convex strips meet.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::Containment;
use catallax_features::{
    Built, evaluate, fillet::FilletDefinition, fillet_corner::CornerDefinition,
};
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_recipe::{Operation, Origin, Recipe, Role, StepId};
use catallax_shapes::{Circle, Curve, Cylinder, Line, Plane, Point3, RPoint3, Surface, Vector3};
use catallax_topology::{FaceId, Loop};
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

// The all-edges triangular prism `(0,0) (8,0) (0,6)` extruded 10, with one case
// per placement. `map` carries a point of the unit-scale unplaced prism, and
// lengths of volume and area scale by `scale`.
struct Skew {
    recipe: Recipe,
    radius: Rational,
    scale: Rational,
    map: fn(&RPoint3) -> RPoint3,
}
fn triangle_prism(radius: &str, scale: &Rational, rotation: [&str; 4], shift: [&str; 3]) -> Recipe {
    let at = |x: i64| (&Rational::from_i64(x) * scale).to_string();
    Recipe::parse(
        &json!({
            "format": "catallax-recipe/4", "name": "skew prism",
            "steps": [
                {"id":"section","op":"sketch",
                 "plane":{"rotation":rotation,"translation":shift},
                 "loops":[[
                    {"kind":"line","id":"base","to":[at(8),at(0)]},
                    {"kind":"line","id":"slope","to":[at(0),at(6)]},
                    {"kind":"line","id":"upright","to":[at(0),at(0)]}
                 ]]},
                {"id":"base","op":"extrude","profile":"section","distance":at(10)},
                {"id":"rounded","op":"fillet","body":"base","radius":radius,"edges":"all"}
            ],
            "result": "rounded"
        })
        .to_string(),
    )
    .unwrap()
}
fn local(p: &RPoint3) -> RPoint3 {
    p.clone()
}
fn halved(p: &RPoint3) -> RPoint3 {
    let half = r("1/2");
    RPoint3::new(&p.x * &half, &p.y * &half, &p.z * &half)
}
// The report's exact rotation for quaternion (1,2,3,4), with its translation.
fn placed(p: &RPoint3) -> RPoint3 {
    let row = |a: &str, b: &str, c: &str, t: &str| {
        &(&(&r(a) * &p.x) + &(&r(b) * &p.y)) + &(&(&r(c) * &p.z) + &r(t))
    };
    RPoint3::new(
        row("-2/3", "2/15", "11/15", "7"),
        row("2/3", "-1/3", "2/3", "11"),
        row("1/3", "14/15", "2/15", "-13"),
    )
}
fn skew_cases() -> [Skew; 4] {
    let plain = (["1", "0", "0", "0"], ["0", "0", "0"]);
    let turned = (["1", "2", "3", "4"], ["7", "11", "-13"]);
    let case = |radius: &str, unit: &str, scale: &str, (q, t), map| Skew {
        recipe: triangle_prism(radius, &r(scale), q, t),
        radius: r(unit),
        scale: r(scale),
        map,
    };
    [
        case("1", "1", "1", plain, local),
        case("1/2", "1/2", "1", plain, local),
        case("1", "1", "1", turned, placed),
        case("1/2", "1", "1/2", plain, halved),
    ]
}
// Inner triangle vertices of the rolling-centre prism, with the two outward
// side normals active at each, then its two heights with the end normal.
fn skew_centres(radius: &Rational) -> Vec<(RPoint3, [[Rational; 3]; 3])> {
    let (zero, one, minus) = (r("0"), r("1"), r("-1"));
    let base = [zero.clone(), minus.clone(), zero.clone()];
    let slope = [r("3/5"), r("4/5"), zero.clone()];
    let upright = [minus.clone(), zero.clone(), zero.clone()];
    let mut centres = Vec::new();
    for (x, y, a, b) in [
        (radius.clone(), radius.clone(), &base, &upright),
        (&r("8") - &(radius * &r("3")), radius.clone(), &base, &slope),
        (
            radius.clone(),
            &r("6") - &(radius * &r("2")),
            &slope,
            &upright,
        ),
    ] {
        for (z, n) in [(radius.clone(), &minus), (&r("10") - radius, &one)] {
            centres.push((
                RPoint3::new(x.clone(), y.clone(), z),
                [
                    a.clone(),
                    b.clone(),
                    [zero.clone(), zero.clone(), n.clone()],
                ],
            ));
        }
    }
    centres
}
fn key(p: &RPoint3) -> [Rational; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
fn exact(p: &Point3) -> RPoint3 {
    p.as_rational().expect("rational point")
}
fn skew_volume(radius: &Rational) -> (Rational, Rational) {
    let r2 = radius * radius;
    let r3 = &r2 * radius;
    (
        &(&r("240") - &(&r("108") * &r2)) + &(&r("24") * &r3),
        &(&r("22") * &r2) - &(&(&r("20") * &r3) / &r("3")),
    )
}
fn skew_area(radius: &Rational) -> (Rational, Rational) {
    let r2 = radius * radius;
    (
        &(&r("288") - &(&r("216") * radius)) + &(&r("36") * &r2),
        &(&r("44") * radius) - &(&r("12") * &r2),
    )
}
fn pi_bounds() -> (Rational, Rational) {
    let lo = r("3.1415926535897932384626433832795028841971693993751");
    let hi = &lo + &r("1e-49");
    (lo, hi)
}
// The measured enclosure overlaps the reference enclosure and is narrow enough.
fn overlaps(
    enclosure: &catallax_measure::Enclosure,
    exact: &(Rational, Rational),
    scale: &Rational,
    width: &str,
) {
    let (lo, hi) = pi_bounds();
    let below = &(&exact.0 + &(&exact.1 * &lo)) * scale;
    let above = &(&exact.0 + &(&exact.1 * &hi)) * scale;
    assert!(
        enclosure.lo <= above && enclosure.hi >= below,
        "[{}, {}] misses [{below}, {above}]",
        enclosure.lo,
        enclosure.hi
    );
    assert!(&enclosure.hi - &enclosure.lo <= r(width), "{enclosure:?}");
}
fn kinds(built: &Built) -> [usize; 5] {
    let mut count = [0; 5];
    for face in built.model.faces_of_body(built.result) {
        let index = match built.model.geom().surface(built.model.face(face).surface) {
            Surface::Plane(_) => 0,
            Surface::Cylinder(_) => 1,
            Surface::Sphere(_) => 2,
            Surface::Torus(_) => 3,
            _ => 4,
        };
        count[index] += 1;
    }
    count
}
fn roles(built: &Built) -> BTreeSet<Origin> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| built.model.face_origins(face).to_vec())
        .collect()
}
fn sphere_at(built: &Built, centre: &RPoint3) -> FaceId {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|f| {
            matches!(built.model.geom().surface(built.model.face(*f).surface),
                Surface::Sphere(s) if key(&exact(&s.center)) == key(centre))
        })
        .expect("sphere about the expected centre")
}
// A sphere sector's boundary: exact great-circle arcs shared with three
// different cylinder strips, through the three expected contact feet.
fn check_sector(built: &Built, face: FaceId, radius: &Rational, feet: &BTreeSet<[Rational; 3]>) {
    let model = &built.model;
    let [lp] = model.face(face).loops.as_slice() else {
        panic!("one sector loop expected");
    };
    let Loop::Coedges(coedges) = model.loop_(*lp) else {
        panic!("sector loop has edges");
    };
    assert_eq!(coedges.len(), 3);
    let mut neighbours = BTreeSet::new();
    let mut corners = BTreeSet::new();
    for c in coedges {
        let edge = model.coedge(*c).edge;
        let Curve::Circle(circle) = model.geom().curve(model.edge(edge).curve) else {
            panic!("sector boundary is not a circle");
        };
        assert_eq!(circle.radius_sq.as_rational(), Some(radius * radius));
        let uses = model.coedges_of_edge(edge);
        assert_eq!(uses.len(), 2);
        let strip = uses
            .iter()
            .map(|c| model.owner_of_coedge(*c).unwrap().1)
            .find(|f| *f != face)
            .unwrap();
        assert!(matches!(
            model.geom().surface(model.face(strip).surface),
            Surface::Cylinder(_)
        ));
        neighbours.insert(strip);
        let edge = model.edge(edge);
        for v in [edge.start.unwrap(), edge.end.unwrap()] {
            corners.insert(key(&exact(&model.geom().position(model.vertex(v).point))));
        }
    }
    assert_eq!(neighbours.len(), 3);
    assert_eq!(&corners, feet);
}

#[test]
fn skew_prism_builds_exact_sectors_in_every_placement() {
    let mut first_names = None;
    for case in skew_cases() {
        let built = build(&case.recipe);
        let model = &built.model;
        assert_eq!(kinds(&built), [5, 9, 6, 0, 0]);
        let body = model.body(built.result);
        let [lump] = body.lumps.as_slice() else {
            panic!("one lump expected");
        };
        assert_eq!(model.lump(*lump).shells.len(), 1);
        let edges = model.edges_of_body(built.result);
        assert_eq!(edges.len(), 36);
        assert_eq!(model.vertices_of_body(built.result).len(), 18);
        assert_eq!(18 - 36 + model.faces_of_body(built.result).len() as i64, 2);
        let lines = edges
            .iter()
            .filter(|e| matches!(model.geom().curve(model.edge(**e).curve), Curve::Line(_)))
            .count();
        assert_eq!(lines, 18);
        let actual = &case.radius * &case.scale;
        let mut vertices = BTreeSet::new();
        for (centre, normals) in skew_centres(&case.radius) {
            let face = sphere_at(&built, &(case.map)(&centre));
            let Surface::Sphere(sphere) = model.geom().surface(model.face(face).surface) else {
                unreachable!()
            };
            assert_eq!(sphere.radius.as_rational(), Some(actual.clone()));
            let feet: BTreeSet<_> = normals
                .iter()
                .map(|n| {
                    key(&(case.map)(&RPoint3::new(
                        &centre.x + &(&case.radius * &n[0]),
                        &centre.y + &(&case.radius * &n[1]),
                        &centre.z + &(&case.radius * &n[2]),
                    )))
                })
                .collect();
            check_sector(&built, face, &actual, &feet);
            vertices.extend(feet);
        }
        let found: BTreeSet<_> = model
            .vertices_of_body(built.result)
            .into_iter()
            .map(|v| key(&exact(&model.geom().position(model.vertex(v).point))))
            .collect();
        assert_eq!(found, vertices);
        for e in &edges {
            if let Curve::Line(rail) = model.geom().curve(model.edge(*e).curve) {
                let on_strip = model.coedges_of_edge(*e).iter().any(|c| {
                    let face = model.owner_of_coedge(*c).unwrap().1;
                    matches!(
                        model.geom().surface(model.face(face).surface),
                        Surface::Cylinder(c)
                            if rail.dir.cross(&c.axis).as_rational()
                                .is_some_and(|v| [v.x, v.y, v.z].iter().all(Rational::is_zero))
                    )
                });
                assert!(on_strip, "rail without a parallel strip");
            }
        }
        let origins = roles(&built);
        let count = |f: fn(&Role) -> bool| origins.iter().filter(|o| f(&o.role)).count();
        assert_eq!(count(|o| matches!(o, Role::Blend { .. })), 9);
        assert_eq!(count(|o| matches!(o, Role::Corner { .. })), 6);
        assert_eq!(origins.len(), 20);
        assert!(
            model
                .faces_of_body(built.result)
                .iter()
                .all(|f| model.face_origins(*f).len() == 1)
        );
        let names = first_names.get_or_insert_with(|| origins.clone());
        assert_eq!(*names, origins);
        let mass = catallax_measure::mass_properties(
            model,
            built.result,
            &Default::default(),
            &mut decider(),
        )
        .unwrap();
        assert!(mass.marks.is_empty(), "{:?}", mass.marks);
        let cubic = &(&case.scale * &case.scale) * &case.scale;
        overlaps(
            &mass.value.volume,
            &skew_volume(&case.radius),
            &cubic,
            "1e-4",
        );
        overlaps(
            &mass.value.area,
            &skew_area(&case.radius),
            &(&case.scale * &case.scale),
            "1e-5",
        );
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        let field: BTreeSet<_> = plan
            .faces()
            .iter()
            .filter(|f| matches!(f.origin.role, Role::Blend { .. } | Role::Corner { .. }))
            .map(|f| f.origin.clone())
            .collect();
        let brep: BTreeSet<_> = origins
            .iter()
            .filter(|o| matches!(o.role, Role::Blend { .. } | Role::Corner { .. }))
            .cloned()
            .collect();
        assert_eq!(field, brep);
    }
}

#[test]
fn skew_prism_body_box_is_the_rounded_extent() {
    let cases = skew_cases();
    for (index, low, high) in [
        (0, ["0", "0", "0"], ["6", "5", "10"]),
        (1, ["0", "0", "0"], ["7", "11/2", "10"]),
        (2, ["53/15", "10", "-63/5"], ["217/15", "21", "-101/15"]),
    ] {
        let built = build(&cases[index].recipe);
        let bounds = catallax_measure::bounding_box(
            &built.model,
            built.result,
            &Default::default(),
            &mut decider(),
        )
        .unwrap()
        .value
        .unwrap();
        let margin = r("1e-6");
        for (found, want) in bounds.min.iter().zip(low) {
            assert!(found <= &r(want) && found + &margin >= r(want), "{found}");
        }
        for (found, want) in bounds.max.iter().zip(high) {
            assert!(
                found >= &r(want) && found <= &(&r(want) + &margin),
                "{found}"
            );
        }
    }
}

#[test]
fn skew_material_agrees_with_analytic_controls_in_both_evaluators() {
    use Containment::{Inside, Outside};
    let rows = [
        (("2", "2", "5"), Inside, Inside),
        (("8", "0", "0"), Outside, Outside),
        (("11/2", "1/2", "1/2"), Inside, Inside),
        (("29/5", "2/5", "1/5"), Outside, Inside),
        (("13/3", "1/3", "2/3"), Inside, Inside),
    ];
    let mut workspace = Workspace::default();
    let mut dec = decider();
    for case in skew_cases() {
        let mut built = build(&case.recipe);
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        for ((x, y, z), at_one, at_half) in rows {
            let want = if case.radius == r("1") {
                at_one
            } else {
                at_half
            };
            let p = (case.map)(&RPoint3::new(r(x), r(y), r(z)));
            assert!(
                matches!(
                    catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut dec),
                    Outcome::Proven(actual) if actual == want
                ),
                "B-rep at {p:?}"
            );
            let location = plan
                .proven_checked(&p, &mut workspace, &mut dec)
                .unwrap()
                .location;
            let want = match want {
                Inside => Location::Inside,
                Outside => Location::Outside,
            };
            assert_eq!(location, want, "field at {p:?}");
        }
        if case.radius != r("1") {
            continue;
        }
        // The first point is in the normal cone at (5,1,1) although a dot
        // product with the slope's normal is negative; the second is on the
        // same sphere's support outside the sector.
        let on = (case.map)(&RPoint3::new(r("17/3"), r("1/3"), r("2/3")));
        let off = (case.map)(&RPoint3::new(r("13/3"), r("1/3"), r("2/3")));
        let centre = (case.map)(&RPoint3::new(r("5"), r("1"), r("1")));
        let face = sphere_at(&built, &centre);
        assert!(matches!(
            catallax_boolean::classify_point(&mut built.model, built.result, &on, &mut dec),
            Outcome::ProvenCoincident(_)
        ));
        assert!(catallax_boolean::face_contains(&mut built.model, face, &on, &mut dec).unwrap());
        assert!(!catallax_boolean::face_contains(&mut built.model, face, &off, &mut dec).unwrap());
        assert_eq!(
            plan.proven_checked(&on, &mut workspace, &mut dec)
                .unwrap()
                .location,
            Location::Boundary
        );
        let origin = built.model.face_origins(face)[0].clone();
        let field_face = plan.select(&origin).unwrap();
        let query = |p: &RPoint3, workspace: &mut Workspace, dec: &mut Decider| {
            plan.face_query_checked(field_face, p, workspace, dec)
                .unwrap()
                .state
        };
        assert_eq!(query(&on, &mut workspace, &mut dec), FaceState::OnFace);
        assert_eq!(query(&off, &mut workspace, &mut dec), FaceState::NotOnFace);
    }
}

#[test]
fn skew_prism_display_is_closed_over_every_region() {
    let cases = skew_cases();
    let built = build(&cases[0].recipe);
    assert_eq!(built.model.faces_of_body(built.result).len(), 20);
    check_display(&built);
}

#[test]
fn zero_radius_skew_prism_keeps_its_invalid_radius_refusal() {
    let recipe = triangle_prism("0", &r("1"), ["1", "0", "0", "0"], ["0", "0", "0"]);
    let Err(catallax_features::BuildError::Geometry { refusal, .. }) =
        evaluate(&recipe, &mut decider())
    else {
        panic!("zero radius built");
    };
    assert_eq!(refusal.code, "fillet.non-positive-radius", "{refusal:?}");
    let refusal = Plan::compile_checked(&recipe, &mut decider()).unwrap_err();
    assert_eq!(refusal.code, "fillet.non-positive-radius", "{refusal:?}");
}

#[test]
fn sheared_corner_is_a_sphere_sector_in_both_evaluators() {
    // A sheared prism meets its end and floor at an oblique side with a
    // rational unit normal. Inward offsets give the centre (19/2,1,9).
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
    let mut built = build(&recipe);
    assert_eq!(kinds(&built), [6, 3, 1, 0, 0]);
    let centre = RPoint3::new(r("19/2"), r("1"), r("9"));
    let face = sphere_at(&built, &centre);
    let feet = [
        RPoint3::new(r("19/2"), r("0"), r("9")),
        RPoint3::new(r("103/10"), r("2/5"), r("9")),
        RPoint3::new(r("19/2"), r("1"), r("10")),
    ]
    .iter()
    .map(key)
    .collect();
    check_sector(&built, face, &r("1"), &feet);
    let on = RPoint3::new(r("61/6"), r("1/3"), r("28/3"));
    let mut dec = decider();
    assert!(matches!(
        catallax_boolean::classify_point(&mut built.model, built.result, &on, &mut dec),
        Outcome::ProvenCoincident(_)
    ));
    assert!(catallax_boolean::face_contains(&mut built.model, face, &on, &mut dec).unwrap());
    let plan = Plan::compile_checked(&recipe, &mut decider()).unwrap();
    let mut workspace = Workspace::default();
    assert_eq!(
        plan.proven_checked(&on, &mut workspace, &mut dec)
            .unwrap()
            .location,
        Location::Boundary
    );
    let origin = built.model.face_origins(face)[0].clone();
    let query = plan
        .face_query_checked(plan.select(&origin).unwrap(), &on, &mut workspace, &mut dec)
        .unwrap();
    assert_eq!(query.state, FaceState::OnFace);
}

// Strips of one radius between planes through the origin of the listed
// outward normals, with the listed edge directions. No B-rep is involved.
fn strip(
    supports: [&Surface; 2],
    dir: [i64; 3],
    radius: &str,
    convex: bool,
    dec: &mut Decider,
) -> FilletDefinition {
    let origin = Origin {
        step: StepId::new("round").unwrap(),
        role: Role::End,
    };
    FilletDefinition::new(
        origin,
        supports.map(Surface::clone),
        [true, true],
        Curve::Line(Line {
            origin: Point3::origin(),
            dir: vector(dir[0], dir[1], dir[2]),
        }),
        r(radius).into(),
        convex,
        dec,
    )
    .unwrap()
}
fn vector(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}
fn through_origin(normal: [&str; 3], reference: [i64; 3]) -> Surface {
    let [x, y, z] = normal.map(|v| Real::from(r(v)));
    Surface::Plane(Plane {
        origin: Point3::origin(),
        normal: Vector3::new(x, y, z),
        ref_dir: vector(reference[0], reference[1], reference[2]),
    })
}

#[test]
fn skew_sector_is_the_positive_cone_of_the_contact_normals() {
    // The acute prism corner moved to the origin: its ball centre is (-3,1,1),
    // and both points are on the ball with displacements of length one.
    let mut dec = decider();
    let [floor, slope, end] = [
        through_origin(["0", "-1", "0"], [1, 0, 0]),
        through_origin(["3/5", "4/5", "0"], [0, 0, 1]),
        through_origin(["0", "0", "-1"], [1, 0, 0]),
    ];
    let strips = [
        strip([&floor, &slope], [0, 0, 1], "1", true, &mut dec),
        strip([&slope, &end], [-4, 3, 0], "1", true, &mut dec),
        strip([&end, &floor], [1, 0, 0], "1", true, &mut dec),
    ];
    let origin = Origin {
        step: StepId::new("round").unwrap(),
        role: Role::End,
    };
    let corner =
        CornerDefinition::new(origin, [&strips[0], &strips[1], &strips[2]], &mut dec).unwrap();
    assert_eq!(
        key(&exact(corner.centre())),
        key(&RPoint3::new(r("-3"), r("1"), r("1")))
    );
    assert_eq!(key(&exact(&corner.vertex)), [r("0"), r("0"), r("0")]);
    let inside = RPoint3::new(r("-7/3"), r("1/3"), r("2/3")).to_real();
    assert!(corner.in_sector(&inside, &mut dec).unwrap());
    let outside = RPoint3::new(r("-11/3"), r("1/3"), r("2/3")).to_real();
    assert!(!corner.in_sector(&outside, &mut dec).unwrap());
}

#[test]
fn corner_definition_keeps_its_refusals() {
    let mut dec = decider();
    let [x, y, z] = [
        through_origin(["1", "0", "0"], [0, 1, 0]),
        through_origin(["0", "1", "0"], [1, 0, 0]),
        through_origin(["0", "0", "1"], [1, 0, 0]),
    ];
    let name = || Origin {
        step: StepId::new("round").unwrap(),
        role: Role::End,
    };
    let good = [
        strip([&x, &y], [0, 0, 1], "1", true, &mut dec),
        strip([&y, &z], [1, 0, 0], "1", true, &mut dec),
        strip([&z, &x], [0, 1, 0], "1", true, &mut dec),
    ];
    let reject = |strips: [&FilletDefinition; 3], detail: &str, dec: &mut Decider| {
        let refusal = CornerDefinition::new(name(), strips, dec).unwrap_err();
        assert_eq!(refusal.code, "fillet.unsupported-corner", "{refusal:?}");
        assert!(refusal.detail.contains(detail), "{refusal:?}");
    };
    assert!(CornerDefinition::new(name(), [&good[0], &good[1], &good[2]], &mut dec).is_ok());
    let concave = strip([&x, &y], [0, 0, 1], "1", false, &mut dec);
    reject([&concave, &good[1], &good[2]], "concave", &mut dec);
    let wider = strip([&x, &y], [0, 0, 1], "2", true, &mut dec);
    reject([&good[0], &wider, &good[2]], "unequal", &mut dec);
    let wall = Surface::Cylinder(Cylinder {
        origin: Point3::origin(),
        axis: vector(0, 0, 1),
        radius_sq: Real::from_i64(25),
        ref_dir: vector(1, 0, 0),
    });
    let curved = FilletDefinition::new(
        name(),
        [z.clone(), wall],
        [true, true],
        Curve::Circle(Circle {
            center: Point3::origin(),
            normal: vector(0, 0, 1),
            radius_sq: Real::from_i64(25),
            ref_dir: vector(1, 0, 0),
        }),
        Real::from_i64(1),
        true,
        &mut dec,
    )
    .unwrap();
    reject([&curved, &good[1], &good[2]], "planar", &mut dec);
    // Three planes through one line, in the planes of a triangular prism.
    let slope = through_origin(["3/5", "4/5", "0"], [0, 0, 1]);
    let dependent = [
        strip([&x, &y], [0, 0, 1], "1", true, &mut dec),
        strip([&y, &slope], [0, 0, 1], "1", true, &mut dec),
        strip([&slope, &x], [0, 0, 1], "1", true, &mut dec),
    ];
    reject(
        [&dependent[0], &dependent[1], &dependent[2]],
        "independent",
        &mut dec,
    );
    // The ball touching these three planes lies past the vertex along the
    // edge of x and y, so that strip has no finite run to the corner.
    let lean = through_origin(["2/3", "2/3", "1/3"], [1, -1, 0]);
    let past = [
        strip([&x, &y], [0, 0, 1], "1", true, &mut dec),
        strip([&y, &lean], [1, 0, -2], "1", true, &mut dec),
        strip([&lean, &x], [0, 1, -2], "1", true, &mut dec),
    ];
    reject(
        [&past[0], &past[1], &past[2]],
        "beyond the vertex",
        &mut dec,
    );
}
