//! P5 flanged housing: the junction fillet continues over the four countersink
//! mouths as torus pieces between the holes and one checked rolling-envelope
//! face per hole, for the authored radius 2 and its edit to 5/2.

use std::collections::BTreeSet;

use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Operation, Origin, Recipe, Role};
use catallax_shapes::{ROLLING_ENVELOPE_TARGET_NM, SplineSupport, Surface};
use catallax_topology::Loop;

fn recipe(radius: &str) -> Recipe {
    let mut recipe = Recipe::parse(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    let last = recipe.steps.len() - 1;
    let Operation::Fillet { radius: r, .. } = &mut recipe.steps[last].op else {
        panic!("P5 ends with its junction fillet");
    };
    *r = serde_json::from_value(serde_json::json!(radius)).unwrap();
    recipe
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}

#[test]
fn p5_junction_continues_over_countersink_mouths() {
    let target = Rational::from_ratio(ROLLING_ENVELOPE_TARGET_NM, 1_000_000).unwrap();
    let mut names: Option<BTreeSet<Origin>> = None;
    for radius in ["2", "5/2"] {
        let mut dec = decider();
        let built = evaluate(&recipe(radius), &mut dec).unwrap();
        let model = &built.model;
        assert_eq!(model.check(built.result), vec![]);
        let faces = model.faces_of_body(built.result);
        assert_eq!(faces.len(), 22);
        let (mut tori, mut envelopes) = (0, 0);
        let mut blends = BTreeSet::new();
        for face in &faces {
            let surface = model.geom().surface(model.face(*face).surface);
            let blend = model
                .face_origins(*face)
                .iter()
                .any(|o| matches!(o.role, Role::Blend { .. }));
            match surface {
                Surface::Torus(_) => {
                    assert!(blend);
                    tori += 1;
                }
                Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                    assert!(blend);
                    assert!(envelope.bound().upper <= target);
                    envelopes += 1;
                }
                _ => assert!(!blend, "only the junction strip carries the blend origin"),
            }
            if blend {
                blends.extend(model.face_origins(*face).iter().cloned());
            }
        }
        assert_eq!((tori, envelopes), (4, 4));
        // Every piece keeps the one blend name, unchanged by the radius edit.
        assert_eq!(blends.len(), 1);
        if let Some(previous) = &names {
            assert_eq!(previous, &blends);
        }
        names = Some(blends);

        // Every boundary vertex lies exactly on every face it bounds, and the
        // shell is the genus-five surface of a bored housing with four holes.
        let mut vertices = BTreeSet::new();
        let mut edges = BTreeSet::new();
        let mut rings = BTreeSet::new();
        let mut loops = 0i64;
        let geom = model.geom();
        for face in &faces {
            let surface = model.face(*face).surface;
            for lp in &model.face(*face).loops {
                loops += 1;
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    continue;
                };
                for c in cs {
                    let edge = model.coedge(*c).edge;
                    edges.insert(edge);
                    let e = model.edge(edge);
                    if e.start.is_none() {
                        rings.insert(edge);
                    }
                    for v in [e.start, e.end].into_iter().flatten() {
                        vertices.insert(v);
                        let point = geom.position(model.vertex(v).point);
                        assert!(matches!(
                            geom.side_of(surface, &point, &mut dec),
                            Outcome::ProvenCoincident(_)
                        ));
                    }
                }
            }
        }
        let euler = vertices.len() as i64 + rings.len() as i64 - edges.len() as i64
            + 2 * faces.len() as i64
            - loops;
        assert_eq!(euler, 2 * (1 - 5));
    }
}

#[test]
fn a_ball_that_would_roll_onto_the_countersink_cone_refuses() {
    // Beyond 2 + sqrt(2) the contact leaves the mouth edge for the cone face,
    // a different continuation; it must refuse rather than build this one.
    let refusal = evaluate(&recipe("7/2"), &mut decider()).unwrap_err();
    let text = format!("{refusal:?}");
    assert!(text.contains("shapes.rolling-envelope"), "{text}");
}
