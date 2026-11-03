//! Exact rolling strips along straight generators: between parallel cylinders,
//! and between a cylinder and a plane parallel to its axis. Expected volumes come
//! from an independent arc-length integration of the section, and the convex
//! cases also agree with the Steiner area of the offset section grown by the
//! radius.

use std::time::Instant;

use catallax_boolean::Containment;
use catallax_features::fillet::FilletDefinition;
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_field::measurement::{CentrePath, FaceDefinition, MeasurementOptions, measure_face};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real, with_operation_log};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::{Curve, Cylinder, Line, Plane, Point3, RPoint3, Surface, Vector3};
use serde_json::{Value, json};

const IDENTITY: [&str; 4] = ["1", "0", "0", "0"];
/// A half turn of `2 atan(4/3)` about x, with a rational rotation matrix.
const TILT: [&str; 4] = ["3/5", "4/5", "0", "0"];

fn cylinder(id: &str, radius: &str, height: &str, rotation: [&str; 4], at: [&str; 3]) -> Value {
    json!({"id":id,"op":"cylinder","radius":radius,"height":height,
        "pose":{"rotation":rotation,"translation":at}})
}
fn walls() -> Value {
    json!({"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
        {"step":"b","role":{"side":{"segment":"wall"}}}]})
}
fn near(y: &str) -> Value {
    let mut edge = walls();
    edge["near"] = json!(["0", y, "6"]);
    edge
}
fn rounded(steps: Vec<Value>, radius: &str, edges: Value) -> Recipe {
    let mut steps = steps;
    let body = steps.last().unwrap()["id"].clone();
    steps.push(json!({"id":"round","op":"fillet","body":body,"radius":radius,"edges":edges}));
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"parallel cylinder round",
            "steps":steps,"result":"round"})
        .to_string(),
    )
    .unwrap()
}
fn joined(kind: &str, a: Value, b: Value) -> Vec<Value> {
    vec![
        a,
        b,
        json!({"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"b"}),
    ]
}
fn pair(kind: &str) -> Vec<Value> {
    joined(
        kind,
        cylinder("a", "5", "12", IDENTITY, ["-3", "0", "0"]),
        cylinder("b", "5", "12", IDENTITY, ["3", "0", "0"]),
    )
}
/// A block through which two overlapping bores leave two cusps of material.
fn bores() -> Vec<Value> {
    vec![
        json!({"id":"section","op":"sketch",
            "plane":{"rotation":IDENTITY,"translation":["0","0","0"]},
            "loops":[[{"kind":"line","id":"s","to":["12","-10"]},
                {"kind":"line","id":"e","to":["12","10"]},
                {"kind":"line","id":"n","to":["-12","10"]},
                {"kind":"line","id":"w","to":["-12","-10"]}]]}),
        json!({"id":"block","op":"extrude","profile":"section","distance":"10"}),
        cylinder("a", "5", "12", IDENTITY, ["-3", "0", "-1"]),
        cylinder("b", "5", "12", IDENTITY, ["3", "0", "-1"]),
        json!({"id":"one","op":"boolean","kind":"subtract","target":"block","tool":"a"}),
        json!({"id":"two","op":"boolean","kind":"subtract","target":"one","tool":"b"}),
    ]
}
/// A rectangular block over z = 0 to 12, whose north side is the flat.
fn block(id: &str, [west, east, south, north]: [&str; 4], rotation: [&str; 4]) -> [Value; 2] {
    let section = format!("{id}-section");
    [
        json!({"id":section,"op":"sketch",
            "plane":{"rotation":rotation,"translation":["0","0","0"]},
            "loops":[[{"kind":"line","id":"s","to":[east,south]},
                {"kind":"line","id":"e","to":[east,north]},
                {"kind":"line","id":"n","to":[west,north]},
                {"kind":"line","id":"w","to":[west,south]}]]}),
        json!({"id":id,"op":"extrude","profile":section,"distance":"12"}),
    ]
}
fn flat() -> Value {
    json!({"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
        {"step":"block","role":{"side":{"segment":"n"}}}]})
}
/// A cylinder of radius 5 about the z axis, joined to a block.
fn with_block(kind: &str, bounds: [&str; 4], rotation: [&str; 4]) -> Vec<Value> {
    let mut steps = vec![cylinder("a", "5", "12", IDENTITY, ["0", "0", "0"])];
    steps.extend(block("block", bounds, rotation));
    steps.push(json!({"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"block"}));
    steps
}
const FLAT: [&str; 4] = ["-12", "12", "-10", "3"];
/// The cylinder with a flat at y = 3, whose two generators are (+-4, 3).
fn d_shaft() -> Vec<Value> {
    with_block("intersect", FLAT, IDENTITY)
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn point(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(q);
    RPoint3::new(x, y, z)
}
fn classify(built: &mut Built, points: &[[&str; 3]], want: Containment) {
    for p in points {
        let actual = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &point(*p),
            &mut decider(),
        );
        assert!(
            matches!(actual, Outcome::Proven(found) if found == want),
            "at {p:?}, wanted {want:?}: {actual:?}"
        );
    }
}
fn assert_volume(built: &Built, want: &str) {
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value;
    let want = q(want);
    let tolerance = q("1e-5");
    assert!(&volume.hi - &volume.lo <= tolerance, "{volume:?}");
    assert!(
        &want - &tolerance <= volume.lo && volume.hi <= &want + &tolerance,
        "{} not within 1e-5 of {}",
        volume.lo.to_f64_nearest(),
        want.to_f64_nearest()
    );
}

/// The two-sided proof of every strip: its contacts lie on its supports, its
/// centre line is the signed offset from both axes, and the whole trimmed face
/// is exactly one radius from that line.
fn prove_strips(built: &mut Built, dec: &mut Decider) {
    let started = Instant::now();
    let faces = built.model.faces_of_body(built.result);
    let mut proved = Vec::new();
    for definition in &built.fillets {
        let Surface::Cylinder(strip) = &definition.surface else {
            panic!("a generator strip must be a cylinder");
        };
        let shift = if definition.convex { 1 } else { -1 };
        for (i, (support, contact)) in definition
            .supports
            .iter()
            .zip(&definition.contacts)
            .enumerate()
        {
            let Curve::Line(contact) = contact else {
                panic!("a generator strip needs straight contacts");
            };
            let side = if definition.same_sense[i] {
                shift
            } else {
                -shift
            };
            let inset = definition.radius.clone() * from(side);
            let checks = match support {
                Surface::Cylinder(support) => {
                    let radial = |p: &Point3| {
                        let d = p.minus(&support.origin);
                        d.norm_squared()
                            - d.dot(&support.axis).square() / support.axis.norm_squared()
                    };
                    let offset = support.radius() - inset;
                    [
                        (radial(&contact.origin), support.radius_sq.clone()),
                        (radial(&strip.origin), offset.square()),
                    ]
                }
                Surface::Plane(plane) => {
                    let height = |p: &Point3| p.minus(&plane.origin).dot(&plane.normal);
                    [
                        (height(&contact.origin), Real::zero()),
                        (
                            height(&strip.origin).square() / plane.normal.norm_squared(),
                            inset.square(),
                        ),
                    ]
                }
                other => panic!("no generator strip on {other:?}"),
            };
            let reach = (
                contact.origin.minus(&strip.origin).norm_squared(),
                definition.radius.square(),
            );
            for (value, want) in checks.into_iter().chain([reach]) {
                assert!(dec.sign(&(value - want)).is_coincident());
            }
        }
        let face = faces
            .iter()
            .copied()
            .find(|f| {
                !proved.contains(f)
                    && built.model.face_origins(*f).contains(&definition.origin)
                    && matches!(built.model.geom().surface(built.model.face(*f).surface),
                        Surface::Cylinder(c)
                            if dec.sign(&c.origin.minus(&strip.origin).cross(&strip.axis)
                                .norm_squared()).is_coincident())
            })
            .expect("each strip definition has its own face");
        proved.push(face);
        let measurement = measure_face(
            &mut built.model,
            face,
            Some(&FaceDefinition::Fillet {
                path: CentrePath::Line(Line {
                    origin: strip.origin.clone(),
                    dir: strip.axis.clone(),
                }),
                radius: definition.radius.clone(),
            }),
            &MeasurementOptions::default(),
            dec,
        )
        .unwrap();
        assert!(measurement.converged);
        assert_eq!((measurement.lower, measurement.upper), (q("0"), q("0")));
    }
    eprintln!(
        "certification_ms={:.3}; fitting=not-applicable",
        started.elapsed().as_secs_f64() * 1000.0
    );
}
fn from(value: i64) -> Real {
    Real::from_i64(value)
}

fn check(recipe: &Recipe, strips: usize, volume: &str) -> Built {
    let mut built = evaluate(recipe, &mut decider()).unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    assert_eq!(built.model.body(built.result).lumps.len(), 1);
    assert_eq!(built.fillets.len(), strips);
    assert_volume(&built, volume);
    prove_strips(&mut built, &mut decider());
    built
}

#[test]
fn intersection_removes_exactly_the_two_tips_on_a_cached_input() {
    let recipe = rounded(pair("intersect"), "5/4", json!([walls()]));
    let mut input = recipe.clone();
    input.steps.pop();
    input.result = StepId::new("joined").unwrap();
    let prefix = evaluate_incremental(&input, None, &mut decider()).unwrap();
    let mut dec = decider();
    let (result, phases) =
        with_operation_log(|| evaluate_incremental(&recipe, Some(prefix), &mut dec));
    let mut built = result.unwrap().into_built();
    assert_eq!(built.model.check(built.result), vec![]);
    assert_eq!(built.fillets.len(), 2);
    assert_volume(&built, "258.87423595");
    classify(
        &mut built,
        &[["0", "17/5", "6"], ["0", "-17/5", "6"]],
        Containment::Inside,
    );
    classify(
        &mut built,
        &[["0", "18/5", "6"], ["0", "-18/5", "6"]],
        Containment::Outside,
    );
    prove_strips(&mut built, &mut dec);
    for name in ["fillet.contact", "fillet.sections", "fillet.trimming"] {
        let measured: Vec<_> = phases.iter().filter(|p| p.operation == name).collect();
        assert!(!measured.is_empty(), "{name} was not recorded");
        let total: f64 = measured.iter().map(|p| p.wall_ms).sum();
        eprintln!("{name}={total:.3}ms over {} spans", measured.len());
    }
}

#[test]
fn union_adds_exactly_the_two_notches() {
    let mut built = check(
        &rounded(pair("union"), "11/8", json!([walls()])),
        2,
        "1619.22037815",
    );
    classify(
        &mut built,
        &[["0", "41/10", "6"], ["0", "-41/10", "6"]],
        Containment::Inside,
    );
    classify(
        &mut built,
        &[["0", "22/5", "6"], ["0", "-22/5", "6"]],
        Containment::Outside,
    );
}

#[test]
fn overlapping_bores_lose_their_cusps_on_reversed_walls() {
    let mut built = check(
        &rounded(bores(), "11/8", json!([walls()])),
        2,
        "3450.64968487",
    );
    assert!(built.fillets.iter().all(|f| f.same_sense == [false, false]));
    assert!(built.fillets.iter().all(|f| f.convex));
    classify(&mut built, &[["0", "41/10", "5"]], Containment::Outside);
    classify(&mut built, &[["0", "22/5", "5"]], Containment::Inside);
}

#[test]
fn a_crescent_mixes_wall_senses_and_an_irrational_centre() {
    let mut built = check(
        &rounded(pair("subtract"), "1", json!([walls()])),
        2,
        "663.74442801",
    );
    assert!(
        built
            .fillets
            .iter()
            .all(|f| f.same_sense[0] != f.same_sense[1])
    );
    // The sharp tip is at (0, +-4). Both witnesses lie on the line from the ball
    // centre through it, 1.514 and 0.9 from a centre of radius one.
    classify(
        &mut built,
        &[["-1/6", "3977/1000", "6"]],
        Containment::Outside,
    );
    classify(&mut built, &[["-1", "4", "6"]], Containment::Inside);
}

#[test]
fn unequal_radii_keep_their_own_contacts() {
    check(
        &rounded(
            joined(
                "intersect",
                cylinder("a", "5", "12", IDENTITY, ["0", "0", "0"]),
                cylinder("b", "29/5", "12", IDENTITY, ["36/5", "0", "0"]),
            ),
            "1",
            json!([walls()]),
        ),
        2,
        "232.1495117",
    );
}

#[test]
fn a_tilted_pose_keeps_the_same_solid() {
    check(
        &rounded(
            joined(
                "intersect",
                cylinder("a", "5", "12", TILT, ["-3", "0", "0"]),
                cylinder("b", "5", "12", TILT, ["3", "0", "0"]),
            ),
            "5/4",
            json!([walls()]),
        ),
        2,
        "258.87423595",
    );
}

#[test]
fn a_near_point_rounds_one_generator_and_refuses_a_tie() {
    let recipe = rounded(pair("intersect"), "5/4", json!([near("3")]));
    let mut built = check(&recipe, 1, "263.62568338");
    classify(&mut built, &[["0", "18/5", "6"]], Containment::Outside);
    classify(&mut built, &[["0", "-39/10", "6"]], Containment::Inside);
    let BuildError::Geometry { refusal, .. } = evaluate(
        &rounded(pair("intersect"), "5/4", json!([near("0")])),
        &mut decider(),
    )
    .unwrap_err() else {
        panic!("a tie must refuse");
    };
    assert_eq!(refusal.code, "recipe.ambiguous-reference");
}

#[test]
fn radius_edits_keep_blend_names_and_step_output_stays_exact() {
    for (kind, radii) in [("intersect", ["5/4", "1"]), ("union", ["11/8", "1"])] {
        let mut previous = None;
        for radius in radii {
            let built = evaluate(
                &rounded(pair(kind), radius, json!([walls()])),
                &mut decider(),
            )
            .unwrap();
            let mut names: Vec<_> = built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .flat_map(|f| built.model.face_origins(f).to_vec())
                .filter(|o| matches!(o.role, Role::Blend { .. }))
                .collect();
            names.sort();
            assert_eq!(names.len(), 2);
            if let Some(old) = previous.replace(names.clone()) {
                assert_eq!(names, old);
            }
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
            assert_eq!(text.matches("ADVANCED_FACE(").count(), 6);
            assert_eq!(text.matches("CYLINDRICAL_SURFACE(").count(), 4);
        }
    }
}

fn refusal(recipe: &Recipe) -> String {
    // A built model is far too large for a panic message.
    let error = evaluate(recipe, &mut decider()).err();
    let Some(BuildError::Geometry { refusal, .. }) = error else {
        panic!("expected a geometric refusal");
    };
    refusal.code.to_string()
}

#[test]
fn unsupported_cases_stay_typed_refusals() {
    // The offset circles of radius 3 touch at the origin; beyond that they
    // miss, and a radius that passes a wall's axis has no contact either.
    assert_eq!(
        refusal(&rounded(pair("intersect"), "2", json!([walls()]))),
        "fillet.radius-too-large"
    );
    for radius in ["3", "6"] {
        assert_eq!(
            refusal(&rounded(pair("intersect"), radius, json!([walls()]))),
            "fillet.radius-too-large"
        );
    }
    // A shorter cylinder ends the strip, but the longer one's caps meet the
    // carrier beyond it. Bounded strip clearance is not proved yet.
    let shorter = joined(
        "union",
        cylinder("a", "5", "12", IDENTITY, ["-3", "0", "0"]),
        cylinder("b", "5", "8", IDENTITY, ["3", "0", "2"]),
    );
    assert_eq!(
        refusal(&rounded(shorter, "11/8", json!([walls()]))),
        "fillet.unsupported-overflow"
    );
    // Every edge of the lens includes the plane rims, whose corners need a patch.
    assert!(refusal(&rounded(pair("intersect"), "1", json!("all"))).starts_with("fillet."));
}

// Crossed cylinders meet in rings that are not generators, and each gets a band.
#[test]
fn crossed_cylinders_round_both_rings() {
    let crossed = joined(
        "union",
        cylinder("a", "5", "12", IDENTITY, ["0", "0", "0"]),
        cylinder(
            "b",
            "3",
            "20",
            ["1/2", "1/2", "1/2", "1/2"],
            ["-10", "0", "6"],
        ),
    );
    let built = evaluate(&rounded(crossed, "1", json!([walls()])), &mut decider()).unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    assert_eq!(built.fillets.len(), 2);
}

#[test]
fn the_independent_field_refuses_cylinder_pairs_rather_than_guessing_a_section() {
    let recipe = rounded(pair("intersect"), "5/4", json!([walls()]));
    let refusal = Plan::compile_checked(&recipe, &mut decider()).unwrap_err();
    assert_eq!(refusal.code, "field.fillet-contact", "{refusal:?}");
}

#[test]
fn a_flat_on_a_cylinder_loses_exactly_its_two_edges_on_a_cached_input() {
    let recipe = rounded(d_shaft(), "7/4", json!([flat()]));
    let mut input = recipe.clone();
    input.steps.pop();
    input.result = StepId::new("joined").unwrap();
    let prefix = evaluate_incremental(&input, None, &mut decider()).unwrap();
    let mut dec = decider();
    let (result, phases) =
        with_operation_log(|| evaluate_incremental(&recipe, Some(prefix), &mut dec));
    let mut built = result.unwrap().into_built();
    assert_eq!(built.model.check(built.result), vec![]);
    assert_eq!(built.fillets.len(), 2);
    assert_volume(&built, "803.89442531");
    classify(
        &mut built,
        &[["18/5", "13/5", "6"], ["-18/5", "13/5", "6"]],
        Containment::Inside,
    );
    classify(
        &mut built,
        &[["39/10", "29/10", "6"], ["-39/10", "29/10", "6"]],
        Containment::Outside,
    );
    prove_strips(&mut built, &mut dec);
    for name in ["fillet.contact", "fillet.sections", "fillet.trimming"] {
        let measured: Vec<_> = phases.iter().filter(|p| p.operation == name).collect();
        assert!(!measured.is_empty(), "{name} was not recorded");
        let total: f64 = measured.iter().map(|p| p.wall_ms).sum();
        eprintln!("{name}={total:.3}ms over {} spans", measured.len());
    }
}

#[test]
fn a_flat_with_an_irrational_centre_matches_its_steiner_volume() {
    check(
        &rounded(d_shaft(), "3/4", json!([flat()])),
        2,
        "807.69782011",
    );
}

#[test]
fn an_extruded_d_profile_rounds_like_the_boolean_one() {
    let section = json!({"id":"section","op":"sketch",
        "plane":{"rotation":IDENTITY,"translation":["0","0","0"]},
        "loops":[[{"kind":"line","id":"flat","to":["4","3"]},
            {"kind":"arc","id":"round","centre":["0","0"],"to":["-4","3"],"turn":"cw"}]]});
    let body = json!({"id":"body","op":"extrude","profile":"section","distance":"12"});
    let edge = json!({"between":[{"step":"body","role":{"side":{"segment":"flat"}}},
        {"step":"body","role":{"side":{"segment":"round"}}}]});
    check(
        &rounded(vec![section, body], "7/4", json!([edge])),
        2,
        "803.89442531",
    );
}

#[test]
fn a_boss_through_a_flat_gains_exactly_its_two_notches() {
    let mut built = check(
        &rounded(with_block("union", FLAT, IDENTITY), "9/4", json!([flat()])),
        2,
        "3881.07230145",
    );
    assert!(built.fillets.iter().all(|f| !f.convex));
    classify(&mut built, &[["81/20", "61/20", "6"]], Containment::Inside);
    classify(&mut built, &[["9/2", "9/2", "6"]], Containment::Outside);
}

#[test]
fn a_bore_through_a_flat_loses_its_two_lips_on_a_reversed_wall() {
    let mut built = check(
        &rounded(with_block_cut(), "9/16", json!([flat()])),
        2,
        "2930.61816418",
    );
    assert!(built.fillets.iter().all(|f| f.convex));
    classify(&mut built, &[["41/10", "59/20", "6"]], Containment::Outside);
    classify(&mut built, &[["5", "13/5", "6"]], Containment::Inside);
}
/// A block with a bore that breaks through its north side.
fn with_block_cut() -> Vec<Value> {
    let mut steps = block("block", FLAT, IDENTITY).to_vec();
    steps.push(cylinder("a", "5", "12", IDENTITY, ["0", "0", "0"]));
    steps.push(json!({"id":"joined","op":"boolean","kind":"subtract","target":"block","tool":"a"}));
    steps
}

#[test]
fn an_oblique_plane_normal_keeps_the_same_solid() {
    check(
        &rounded(
            with_block("intersect", FLAT, ["3/5", "0", "0", "4/5"]),
            "7/4",
            json!([flat()]),
        ),
        2,
        "803.89442531",
    );
}

#[test]
fn a_near_point_rounds_one_edge_of_a_flat() {
    let mut edge = flat();
    edge["near"] = json!(["4", "3", "6"]);
    let mut built = check(&rounded(d_shaft(), "7/4", json!([edge])), 1, "806.09182799");
    classify(&mut built, &[["39/10", "29/10", "6"]], Containment::Outside);
    classify(&mut built, &[["-39/10", "29/10", "6"]], Containment::Inside);
}

#[test]
fn contacts_follow_the_support_order() {
    let wall = Surface::Cylinder(Cylinder {
        origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
        axis: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        radius_sq: from(25),
        ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
    });
    let plane = Surface::Plane(Plane {
        origin: Point3::new(Real::zero(), from(3), Real::zero()),
        normal: Vector3::new(Real::zero(), Real::one(), Real::zero()),
        ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
    });
    let edge = Curve::Line(Line {
        origin: Point3::new(from(4), from(3), Real::zero()),
        dir: Vector3::new(Real::zero(), Real::zero(), Real::one()),
    });
    let origin = Origin {
        step: StepId::new("round").unwrap(),
        role: Role::Start,
    };
    let definition = |supports: [Surface; 2]| {
        FilletDefinition::new(
            origin.clone(),
            supports,
            [true, true],
            edge.clone(),
            Real::from(q("7/4")),
            true,
            &mut decider(),
        )
        .unwrap()
    };
    let forward = definition([plane.clone(), wall.clone()]);
    let reverse = definition([wall, plane]);
    let origin_of = |curve: &Curve| match curve {
        Curve::Line(line) => line.origin.clone(),
        other => panic!("{other:?}"),
    };
    let mut dec = decider();
    for (a, b) in [(0, 1), (1, 0)] {
        let gap = origin_of(&forward.contacts[a]).minus(&origin_of(&reverse.contacts[b]));
        assert!(dec.sign(&gap.norm_squared()).is_coincident());
    }
    let rail = origin_of(&forward.contacts[0]);
    assert!(dec.sign(&(rail.x - from(3))).is_coincident());
    assert!(dec.sign(&(rail.y - from(3))).is_coincident());
}

#[test]
fn a_flat_stays_exact_in_step() {
    let built = evaluate(&rounded(d_shaft(), "7/4", json!([flat()])), &mut decider()).unwrap();
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
    assert_eq!(text.matches("ADVANCED_FACE(").count(), 6);
    assert_eq!(text.matches("CYLINDRICAL_SURFACE(").count(), 3);
}

#[test]
fn flat_refusals_stay_typed() {
    // The offset line and circle touch at r = 4, and miss beyond it.
    assert_eq!(
        refusal(&rounded(d_shaft(), "4", json!([flat()]))),
        "fillet.radius-too-large"
    );
    for radius in ["5", "6"] {
        assert_eq!(
            refusal(&rounded(d_shaft(), radius, json!([flat()]))),
            "fillet.radius-too-large"
        );
    }
    // A plane that is not parallel to the axis meets the wall in an ellipse.
    let tilted = with_block("intersect", FLAT, TILT);
    assert!(refusal(&rounded(tilted, "1", json!([flat()]))).starts_with("fillet."));
    // Every edge includes the rim arcs, whose corners need a patch.
    assert!(refusal(&rounded(d_shaft(), "1", json!("all"))).starts_with("fillet."));
}

/// A corner of the cylinder between x = 0 and y = 3, where the wall and the flat
/// share one generator, at (-4, 3).
fn corner() -> Recipe {
    rounded(
        with_block("intersect", ["-12", "0", "-10", "3"], IDENTITY),
        "7/4",
        json!([flat()]),
    )
}

#[test]
fn the_field_agrees_on_a_single_generator_between_a_cylinder_and_a_flat() {
    let recipe = corner();
    let mut built = check(&recipe, 1, "401.94721265");
    let mut dec = decider();
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    for (p, want, containment) in [
        (
            ["-18/5", "13/5", "6"],
            Location::Inside,
            Containment::Inside,
        ),
        (
            ["-39/10", "29/10", "6"],
            Location::Outside,
            Containment::Outside,
        ),
        (["-3", "1", "6"], Location::Inside, Containment::Inside),
        (["-1", "1", "6"], Location::Inside, Containment::Inside),
        (["1", "1", "6"], Location::Outside, Containment::Outside),
    ] {
        let query = plan
            .proven_checked(&point(p), &mut workspace, &mut dec)
            .unwrap();
        assert_eq!(query.location, want, "field at {p:?}: {query:?}");
        classify(&mut built, &[p], containment);
    }
    // The strip itself: the ball centre (-3, 5/4) plus 7/4 of (-3/5, 4/5).
    let on_strip = ["-81/20", "53/20", "6"];
    let query = plan
        .proven_checked(&point(on_strip), &mut workspace, &mut dec)
        .unwrap();
    assert_eq!(query.location, Location::Boundary, "{query:?}");
}

#[test]
fn the_field_refuses_what_it_cannot_section_exactly() {
    // Both generators of a flat share one blend name.
    let both = rounded(d_shaft(), "7/4", json!([flat()]));
    let refusal = Plan::compile_checked(&both, &mut decider()).unwrap_err();
    assert_eq!(refusal.code, "field.fillet-contact", "{refusal:?}");
    // An irrational ball centre has no rational section.
    let irrational = rounded(
        with_block("intersect", ["-12", "0", "-10", "3"], IDENTITY),
        "3/4",
        json!([flat()]),
    );
    assert!(Plan::compile_checked(&irrational, &mut decider()).is_err());
}
