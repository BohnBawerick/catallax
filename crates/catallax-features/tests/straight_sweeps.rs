//! Straight-path sweeps preserve the rigid-motion definition and physical names.

use catallax_features::{BuildError, evaluate, evaluate_incremental};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Operation, PathSegment, Recipe, Role, SweepOrientation};
use catallax_topology::BodyKind;

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn recipe() -> Recipe {
    Recipe::parse(include_str!("fixtures/straight-sweep.json")).unwrap()
}
fn set_end(recipe: &mut Recipe, end: [i64; 3]) {
    let Operation::Path { segments, .. } = &mut recipe.steps[1].op else {
        unreachable!()
    };
    let PathSegment::Line { to, .. } = &mut segments[0] else {
        unreachable!()
    };
    to.0 = end.map(Into::into);
}

#[test]
fn oblique_signed_sweep_uses_displacement_not_path_position_or_unit_direction() {
    for negative in [false, true] {
        for orientation in [
            SweepOrientation::Fixed,
            SweepOrientation::RotationMinimising,
        ] {
            let mut recipe = recipe();
            if negative {
                set_end(&mut recipe, [98, 199, 297]);
            }
            let Operation::Sweep {
                orientation: mode, ..
            } = &mut recipe.steps[2].op
            else {
                unreachable!()
            };
            *mode = orientation;
            let built = evaluate(&recipe, &mut Decider::default()).unwrap();
            assert!(built.model.check(built.result).is_empty());
            assert!(built.model.body(built.result).marks.is_empty());
            let mass = catallax_measure::mass_properties(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
            assert!(mass.marks.is_empty());
            // Area(R) * |displacement dot normal| = 8 * 3, not 8 * sqrt(14).
            assert_eq!(mass.value.volume.lo, r(24));
            assert_eq!(mass.value.volume.hi, r(24));
            let sign = if negative { -1 } else { 1 };
            for (value, expected) in mass.value.centroid.unwrap().iter().zip([
                r(9 + sign),
                Rational::from_ratio(-4 + sign, 2).unwrap(),
                Rational::from_ratio(22 + 3 * sign, 2).unwrap(),
            ]) {
                assert!(value.contains(&expected), "{value:?} excludes {expected}");
                assert!(value.width() <= Rational::parse("1e-9").unwrap());
            }
            let bounds = catallax_measure::bounding_box(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
            .value
            .unwrap();
            let (low, high) = if negative {
                ([5, -4, 8], [11, -1, 11])
            } else {
                ([7, -3, 11], [13, 0, 14])
            };
            assert_eq!(bounds.min, low.map(r));
            assert_eq!(bounds.max, high.map(r));
            let mut roles = Vec::new();
            for face in built.model.faces_of_body(built.result) {
                for origin in built.model.face_origins(face) {
                    assert_eq!(origin.step.as_str(), "body");
                    if let Role::Side { segment, path } = &origin.role {
                        assert!(["a", "b", "c", "d"].contains(&segment.as_str()));
                        assert_eq!(path.as_ref().unwrap().as_str(), "travel");
                    }
                    roles.push(origin.role.clone());
                }
            }
            assert!(roles.contains(&Role::Start));
            assert!(roles.contains(&Role::End));
            assert_eq!(roles.len(), 6);
            assert!(
                built
                    .model
                    .edges_of_body(built.result)
                    .iter()
                    .all(|e| built.model.edge_name(built.result, *e).is_some())
            );
        }
    }
}

#[test]
fn empty_path_and_in_plane_motion_regularise_to_empty_but_zero_length_line_is_impossible() {
    let mut empty = recipe();
    let Operation::Path { segments, .. } = &mut empty.steps[1].op else {
        unreachable!()
    };
    segments.clear();
    for recipe in [empty, {
        let mut in_plane = recipe();
        set_end(&mut in_plane, [102, 201, 300]);
        in_plane
    }] {
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
        assert!(built.model.check(built.result).is_empty());
    }
    let mut invalid = recipe();
    set_end(&mut invalid, [100, 200, 300]);
    let BuildError::Geometry { step, refusal } =
        evaluate(&invalid, &mut Decider::default()).unwrap_err()
    else {
        panic!("expected path failure")
    };
    assert_eq!(step.as_str(), "rail");
    assert_eq!(refusal.code, "path.zero-length");
}

#[test]
fn reused_prefix_retains_path_and_edit_rebuilds_its_displacement() {
    let source = recipe();
    let saved = evaluate_incremental(&source, None, &mut Decider::default()).unwrap();
    let mut suffix = source.clone();
    suffix.steps.push(catallax_recipe::Step {
        id: catallax_recipe::StepId::new("second").unwrap(),
        op: source.steps[2].op.clone(),
    });
    suffix.result = suffix.steps[3].id.clone();
    let continued =
        evaluate_incremental(&suffix, Some(saved.clone()), &mut Decider::default()).unwrap();
    assert_eq!(
        continued
            .built()
            .model
            .counts(continued.built().result)
            .faces,
        6
    );
    let mut edited = source;
    set_end(&mut edited, [100, 200, 306]);
    let changed =
        evaluate_incremental(&edited, Some(saved.clone()), &mut Decider::default()).unwrap();
    for (built, expected) in [
        (saved.built(), 24),
        (continued.built(), 24),
        (changed.built(), 48),
    ] {
        let mass = catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(mass.value.volume.lo, r(expected));
        assert_eq!(mass.value.volume.hi, r(expected));
    }
}

#[test]
fn curved_and_multi_segment_paths_keep_explicit_coverage_refusals() {
    for segments in [
        serde_json::json!([{"kind":"arc","id":"turn","centre":["100","200","303"],"to":["103","200","303"],"axis":["0","1","0"]}]),
        serde_json::json!([{"kind":"line","id":"one","to":["100","200","301"]},{"kind":"line","id":"two","to":["100","200","302"]}]),
    ] {
        let mut value = serde_json::to_value(recipe()).unwrap();
        value["steps"][1]["segments"] = segments;
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let BuildError::Geometry { step, refusal } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected coverage refusal")
        };
        assert_eq!(step.as_str(), "rail");
        assert_eq!(refusal.code, "features.unsupported-path");
    }
}
