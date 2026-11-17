//! Authored sketch joints must survive arc lifting without creating a cut lid.
use catallax_boolean::{Containment, classify_point};
use catallax_measure::{Accuracy, face_area};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::Role;
use catallax_shapes::RPoint3;

#[test]
fn released_application_sketch_opens_the_coplanar_top_without_a_lid() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let mut application = catallax_app::build_soliddraft(
            include_str!("fixtures/coplanar-start-cap.json"),
            &MeshOptions::default(),
            &mut dec,
        )
        .unwrap();
        let built = &mut application.body.built;
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 13);
        let flat: Vec<_> = faces.iter().copied().filter(|&face| {
            built.model.face_origins(face).iter().any(|origin| {
                origin.step == application.steps["op_3"]
                    && matches!(&origin.role, Role::Side { segment, .. } if segment.as_str() == "l0s3")
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
            built.model.face_origins(face).iter().any(|origin| {
                origin.step == application.steps["op_45"] && origin.role == Role::Start
            })
        }));
        let mouth = RPoint3::new(0.into(), (-38).into(), 60.into());
        assert!(matches!(
            classify_point(&mut built.model, built.result, &mouth, &mut dec),
            Outcome::Proven(Containment::Outside)
        ));
    });
    assert!(!work.exhausted, "{work:?}");
}
