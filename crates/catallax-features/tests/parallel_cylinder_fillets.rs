//! Exact rolling strips between parallel cylinders. Expected volumes come from
//! an independent arc-length integration of the section, and the convex cases
//! also agree with the Steiner area of the offset lens grown by the radius.

use std::time::Instant;

use catallax_boolean::Containment;
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_field::Plan;
use catallax_field::measurement::{CentrePath, FaceDefinition, MeasurementOptions, measure_face};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, with_operation_log};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::{Curve, Line, RPoint3, Surface};
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
            panic!("a parallel strip must be a cylinder");
        };
        let shift = if definition.convex { 1 } else { -1 };
        for (i, (support, contact)) in definition
            .supports
            .iter()
            .zip(&definition.contacts)
            .enumerate()
        {
            let (Surface::Cylinder(support), Curve::Line(contact)) = (support, contact) else {
                panic!("a parallel strip needs straight contacts on cylinders");
            };
            let radial = |p: &catallax_shapes::Point3| {
                let d = p.minus(&support.origin);
                d.norm_squared() - d.dot(&support.axis).square() / support.axis.norm_squared()
            };
            let side = if definition.same_sense[i] {
                shift
            } else {
                -shift
            };
            let offset = support.radius() - definition.radius.clone() * from(side);
            for (value, want) in [
                (radial(&contact.origin), support.radius_sq.clone()),
                (radial(&strip.origin), offset.square()),
                (
                    contact.origin.minus(&strip.origin).norm_squared(),
                    definition.radius.square(),
                ),
            ] {
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
fn from(value: i64) -> catallax_numbers::Real {
    catallax_numbers::Real::from_i64(value)
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
    let BuildError::Geometry { refusal, .. } = evaluate(recipe, &mut decider()).unwrap_err() else {
        panic!("expected a geometric refusal");
    };
    refusal.code.to_string()
}

#[test]
fn unsupported_cases_stay_typed_refusals() {
    // The offset circles of radius 3 touch at the origin; beyond that they miss.
    for radius in ["2", "3"] {
        assert_eq!(
            refusal(&rounded(pair("intersect"), radius, json!([walls()]))),
            "fillet.unsupported-contact"
        );
    }
    // A radius that passes a wall's axis has no regular contact.
    assert_eq!(
        refusal(&rounded(pair("intersect"), "6", json!([walls()]))),
        "fillet.unsupported-contact"
    );
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
    // Crossed cylinders meet in curves that are not generators.
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
    assert!(refusal(&rounded(crossed, "1", json!([walls()]))).starts_with("fillet."));
    // Every edge of the lens includes the plane rims, whose corners need a patch.
    assert!(refusal(&rounded(pair("intersect"), "1", json!("all"))).starts_with("fillet."));
}

#[test]
fn the_independent_field_refuses_rather_than_guessing_a_straight_section() {
    let recipe = rounded(pair("intersect"), "5/4", json!([walls()]));
    let refusal = Plan::compile_checked(&recipe, &mut decider()).unwrap_err();
    assert_eq!(refusal.code, "field.fillet-contact", "{refusal:?}");
}
