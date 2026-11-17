//! Exact coplanar pockets spanning a tangent planar and cylindrical boundary.
use catallax_boolean::{Containment, classify_point};
use catallax_measure::{Accuracy, face_area};
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;

#[test]
fn tangent_top_pocket_opens_the_flat_face_without_a_downward_start_cap() {
    let recipe = Recipe::parse(include_str!("fixtures/coplanar-start-cap.json")).unwrap();
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 13);
        let flat: Vec<_> = faces.iter().copied().filter(|&face| {
            built.model.face_origins(face).iter().any(|origin| {
                origin.step.as_str() == "base"
                    && matches!(&origin.role, Role::Side { segment, .. } if segment.as_str() == "flat_top")
            })
        }).collect();
        assert_eq!(flat.len(), 1);
        let accuracy = Accuracy {
            relative: Rational::zero(),
            absolute: Rational::parse("0.000001").unwrap(),
        };
        let area = face_area(&built.model, flat[0], &accuracy, &mut dec)
            .unwrap()
            .value;
        assert!(
            area.contains(&Rational::parse("1894.61").unwrap()),
            "{area:?}"
        );
        assert!(area.width() <= accuracy.absolute);
        assert!(!faces.iter().any(|&face| {
            built
                .model
                .face_origins(face)
                .iter()
                .any(|origin| origin.step.as_str() == "pocket" && origin.role == Role::Start)
        }));
        let mouth = RPoint3::new(0.into(), (-38).into(), 60.into());
        assert!(matches!(
            classify_point(&mut built.model, built.result, &mouth, &mut dec),
            Outcome::Proven(Containment::Outside)
        ));
    });
    assert!(!work.exhausted, "{work:?}");
}
