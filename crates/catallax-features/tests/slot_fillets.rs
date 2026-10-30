//! Extruded slot top boundary: a straight strip runs past its rail ends and is
//! cut by the two cylindrical end walls, with exact volumes, terminal curves,
//! radius-independent names and an immutable input after refusal.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real, StrictSign};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::{Curve, Point3, RPoint3, Surface};
use catallax_topology::{Loop, Sense};
use serde_json::json;

fn recipe(segments: &[&str], radius: &str, version: u8) -> Recipe {
    let mut steps = vec![
        json!({"id":"profile","op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
            "loops":[[
                {"kind":"line","id":"bottom","to":["20","-5"]},
                {"kind":"arc","id":"right","centre":["20","0"],"to":["20","5"],"turn":"ccw"},
                {"kind":"line","id":"top","to":["-20","5"]},
                {"kind":"arc","id":"left","centre":["-20","0"],"to":["-20","-5"],"turn":"ccw"}]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"10"}),
    ];
    let mut result = "body";
    if !segments.is_empty() {
        let edges: Vec<_> = segments
            .iter()
            .map(|segment| {
                json!({"between":[{"step":"body","role":"end"},
                    {"step":"body","role":{"side":{"segment":segment}}}]})
            })
            .collect();
        steps.push(json!({"id":"round","op":"fillet","body":"body","edges":edges,"radius":radius}));
        result = "round";
    }
    Recipe::parse(
        &json!({"format":format!("catallax-recipe/{version}"),"name":"slot tangent ends",
            "steps":steps,"result":result})
        .to_string(),
    )
    .unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn zero(value: Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(&value), Outcome::ProvenCoincident(_))
}
fn positive(value: &Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(value), Outcome::Proven(StrictSign::Positive))
}
fn int(v: i64) -> Real {
    Real::from_i64(v)
}
fn export(built: &Built) -> Vec<u8> {
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    bytes
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
    let want = Rational::parse(want).unwrap();
    let tolerance = Rational::parse("1e-5").unwrap();
    assert!(&volume.hi - &volume.lo <= tolerance, "{volume:?}");
    assert!(&want - &tolerance <= volume.lo && volume.hi <= &want + &tolerance);
}
fn classify(built: &mut Built, points: &[([&str; 3], Containment)]) {
    for (p, want) in points {
        let [x, y, z] = p.map(|v| Rational::parse(v).unwrap());
        let point = RPoint3::new(x, y, z);
        assert!(
            matches!(classify_point(&mut built.model, built.result, &point, &mut decider()),
                Outcome::Proven(actual) if actual == *want),
            "at {p:?}"
        );
    }
}
/// The sorted origins of every face, the blend's included.
fn origins(built: &Built) -> Vec<Vec<catallax_recipe::Origin>> {
    let mut all: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .map(|f| built.model.face_origins(f).to_vec())
        .collect();
    all.sort();
    all
}

#[test]
fn a_straight_strip_ends_on_both_cylindrical_end_walls() {
    for (radius, volume) in [
        (1, "4776.244436115815"),
        (2, "4747.891984181232"),
        (4, "4630.739296215008"),
    ] {
        let mut built = evaluate(&recipe(&["top"], &radius.to_string(), 3), &mut decider())
            .unwrap_or_else(|e| panic!("r={radius}: {e}"));
        let model = &built.model;
        assert_eq!(model.check(built.result), vec![]);
        let lumps = &model.body(built.result).lumps;
        assert_eq!(lumps.len(), 1);
        assert_eq!(model.lump(lumps[0]).shells.len(), 1);
        let faces = model.faces_of_body(built.result);
        let mut kinds = BTreeMap::new();
        let mut blends = Vec::new();
        for face in &faces {
            let surface = model.geom().surface(model.face(*face).surface);
            *kinds.entry(format!("{:?}", surface.kind())).or_insert(0) += 1;
            if model
                .face_origins(*face)
                .iter()
                .any(|o| matches!(o.role, Role::Blend { .. }))
            {
                blends.push(*face);
            }
        }
        assert_eq!(
            kinds,
            BTreeMap::from([("Cylinder".into(), 3), ("Plane".into(), 4)])
        );
        let [blend] = blends[..] else {
            panic!("one blend face, not {blends:?}");
        };
        let mut dec = decider();
        let r = int(radius);
        let Surface::Cylinder(c) = model.geom().surface(model.face(blend).surface) else {
            panic!("the strip is a cylinder");
        };
        assert!(zero(&c.radius_sq - &r.square(), &mut dec));
        assert!(zero(c.axis.y.clone(), &mut dec) && zero(c.axis.z.clone(), &mut dec));
        assert!(zero(&c.origin.y - &(int(5) - r.clone()), &mut dec));
        assert!(zero(&c.origin.z - &(int(10) - r.clone()), &mut dec));

        // Every vertex lies on every face it bounds; every edge has two
        // opposite uses.
        let mut uses: BTreeMap<_, Vec<Sense>> = BTreeMap::new();
        for face in &faces {
            let surface = model.face(*face).surface;
            for lp in &model.face(*face).loops {
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    continue;
                };
                for c in cs {
                    let c = model.coedge(*c);
                    uses.entry(c.edge).or_default().push(c.sense);
                    let e = model.edge(c.edge);
                    for v in [e.start, e.end].into_iter().flatten() {
                        let p = model.geom().position(model.vertex(v).point);
                        assert!(matches!(
                            model.geom().side_of(surface, &p, &mut dec),
                            Outcome::ProvenCoincident(_)
                        ));
                    }
                }
            }
        }
        for senses in uses.values() {
            assert_eq!(senses.len(), 2);
            assert_ne!(senses[0], senses[1]);
        }

        // Both terminal curves are blend/end-wall intersections. They run
        // from the wall contact (+-20, 5, 10-r) to the top contact
        // (+-(20 + sqrt(10r - r^2)), 5-r, 10), beyond the straight segment.
        let mut terminals = 0;
        let mut ends = BTreeSet::new();
        let Loop::Coedges(cs) = model.loop_(model.face(blend).loops[0]) else {
            panic!("the strip has a coedge boundary");
        };
        for c in cs {
            let e = model.edge(model.coedge(*c).edge);
            if matches!(model.geom().curve(e.curve), Curve::Intersection(_)) {
                terminals += 1;
            }
            for v in [e.start, e.end].into_iter().flatten() {
                let p: Point3 = model.geom().position(model.vertex(v).point);
                let side = if positive(&p.x, &mut dec) { 1 } else { -1 };
                let offset = &(&p.x * &int(side)) - &int(20);
                let wall = zero(&p.y - &int(5), &mut dec)
                    && zero(&p.z - &(int(10) - r.clone()), &mut dec)
                    && zero(offset.clone(), &mut dec);
                let top = zero(&p.y - &(int(5) - r.clone()), &mut dec)
                    && zero(&p.z - &int(10), &mut dec)
                    && positive(&offset, &mut dec)
                    && zero(
                        &offset.square() - &(&int(10 * radius) - &r.square()),
                        &mut dec,
                    );
                assert!(wall ^ top, "r={radius}: {p:?}");
                ends.insert((side, top));
            }
        }
        assert_eq!(terminals, 2);
        assert_eq!(ends.len(), 4);
        assert_volume(&built, volume);
        // The strip extends past x=+-20, while the opposite edge stays sharp.
        classify(
            &mut built,
            &[
                (["0", "49/10", "99/10"], Containment::Outside),
                (["22", "9/2", "99/10"], Containment::Outside),
                (["0", "-49/10", "99/10"], Containment::Inside),
                (["249/10", "0", "99/10"], Containment::Inside),
                (["0", "0", "1"], Containment::Inside),
            ],
        );
    }
}

#[test]
fn strip_names_survive_radius_edits_and_the_rounded_edge_is_gone() {
    let mut previous = None;
    for radius in ["1", "4", "1"] {
        let built = evaluate(&recipe(&["top"], radius, 3), &mut decider()).unwrap();
        let names = origins(&built);
        if let Some(previous) = &previous {
            assert_eq!(previous, &names);
        }
        previous = Some(names);
    }
    // The rounded edge no longer exists, although both its faces survive.
    let mut again = recipe(&["top"], "1", 3);
    let mut second = again.steps[2].clone();
    second.id = StepId::new("again").unwrap();
    let catallax_recipe::Operation::Fillet { body, .. } = &mut second.op else {
        unreachable!();
    };
    *body = StepId::new("round").unwrap();
    again.steps.push(second);
    again.result = StepId::new("again").unwrap();
    match evaluate(&again, &mut decider()) {
        Err(BuildError::Geometry { step, refusal }) => {
            assert_eq!(step.as_str(), "again");
            assert_eq!(refusal.code, "fillet.lost-edge");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn oversized_slot_rounds_refuse_without_touching_the_retained_input() {
    let base = evaluate_incremental(&recipe(&[], "0", 3), None, &mut decider()).unwrap();
    let bytes = export(base.built());
    let names = origins(base.built());
    for (segments, code) in [
        (
            &["bottom", "right", "top", "left"][..],
            "fillet.unsupported-contact",
        ),
        (&["top"][..], "fillet.unsupported-overflow"),
    ] {
        match evaluate_incremental(
            &recipe(segments, "6", 3),
            Some(base.clone()),
            &mut decider(),
        ) {
            Err(BuildError::Geometry { step, refusal }) => {
                assert_eq!(step.as_str(), "round");
                assert_eq!(refusal.code, code, "{segments:?}");
            }
            other => panic!("{segments:?}: {other:?}"),
        }
        let retained = base.built();
        assert_eq!(retained.model.check(retained.result), vec![]);
        assert_eq!(export(retained), bytes);
        assert_eq!(origins(retained), names);
    }
    assert_volume(base.built(), "4785.398163397448");
}
