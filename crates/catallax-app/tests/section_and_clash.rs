//! Exact application queries on the release collar and shelled tray.
use catallax_app::queries::{
    SectionPlane, contact_faces, overlap_volume, section_curves, silhouette,
};
use catallax_features::evaluate;
use catallax_intersect::{CurveForm, curve_form};
use catallax_measure::Accuracy;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{RPoint3, RVector3};

#[test]
fn collar_section_overlap_and_shared_faces_use_retained_definitions() {
    let recipe = Recipe::parse(include_str!("../../../parts/P2-collar.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let ring = built.bodies[&StepId::new("ring").unwrap()];
    let section = section_curves(
        &built.model,
        built.result,
        &SectionPlane {
            origin: RPoint3::new(0.into(), 6.into(), 0.into()),
            normal: RVector3::new(0.into(), 1.into(), 0.into()),
        },
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(section.value.curves.len(), 2);
    let mut radii = section
        .value
        .curves
        .iter()
        .map(|edge| {
            assert_eq!((edge.start, edge.end), (None, None));
            let CurveForm::Circle { radius_sq, .. } = curve_form(
                section.value.model.geom(),
                edge.curve,
                &mut Decider::default(),
            )
            .unwrap() else {
                panic!("circle")
            };
            radius_sq.rational_ref().unwrap().clone()
        })
        .collect::<Vec<_>>();
    radii.sort();
    assert_eq!(radii, vec![Rational::from(36), Rational::from(144)]);
    let overlap = overlap_volume(
        &built.model,
        built.result,
        ring,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(overlap.value.lo > Rational::from(3998));
    assert!(overlap.value.hi < Rational::from(3999));
    let contacts =
        contact_faces(&built.model, built.result, ring, &mut Decider::default()).unwrap();
    assert_eq!(contacts.value.first.len(), 4);
    assert_eq!(contacts.value.second.len(), 4);
    assert!(
        contacts
            .value
            .first
            .iter()
            .chain(&contacts.value.second)
            .all(|p| {
                p.orientation == catallax_boolean::FaceClass::OnSame
                    && !contacts.value.model.face_origins(p.face).is_empty()
            })
    );
    let refusal = silhouette(
        &built.model,
        built.result,
        &RVector3::new(0.into(), 1.into(), 0.into()),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.refusal.code, "app.silhouette-surface");
    assert_eq!(refusal.operands[0].faces, 6);
    assert_eq!(refusal.operands[0].curved_faces, 4);
}

#[test]
fn tray_section_retains_outer_and_inner_trimmed_branches() {
    let recipe = Recipe::parse(include_str!("../../../parts/P4-shelled-tray.json")).unwrap();
    let (built, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let built = built.unwrap();
    assert!(!work.exhausted);
    let answer = section_curves(
        &built.model,
        built.result,
        &SectionPlane {
            origin: RPoint3::new(0.into(), 0.into(), 6.into()),
            normal: RVector3::new(0.into(), 0.into(), 1.into()),
        },
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(answer.value.curves.len(), 12);
    assert!(answer.value.points.is_empty());
    assert!(
        answer
            .value
            .curves
            .iter()
            .all(|e| e.start.is_some() && e.end.is_some())
    );
    let mut lines = 0;
    let mut circles = 0;
    for edge in &answer.value.curves {
        match curve_form(
            answer.value.model.geom(),
            edge.curve,
            &mut Decider::default(),
        )
        .unwrap()
        {
            CurveForm::Line { .. } => lines += 1,
            CurveForm::Circle { .. } => circles += 1,
            form => panic!("unexpected section: {form:?}"),
        }
    }
    assert_eq!((lines, circles), (8, 4));
}
