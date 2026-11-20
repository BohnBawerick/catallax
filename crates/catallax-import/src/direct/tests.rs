use super::*;
use crate::tolerant::{healing::HealingChange, opened};
use catallax_recipe::StepId;

use crate::incidence_features::common::{BoxFixture, box_step, r};

fn max_x(solid: &DirectSolid) -> String {
    let BodyView::Numerical(body) = solid.body() else {
        panic!("edited material must be numerical");
    };
    let mut dec = Decider::default();
    let written = catallax_step::numerical::build(body, &Default::default(), &mut dec).unwrap();
    let file = crate::read_step_with_decider(written.bytes(), &dec).unwrap();
    file.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start.point[0], &e.end.point[0]])
        .max()
        .unwrap()
        .to_string()
}

#[test]
fn nonzero_sew_is_canonical_for_a_tighter_second_edit() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let text = box_step(BoxFixture::new(["2", "3", "4"]));
        let opened = opened::open_step(
            text.as_bytes(),
            &crate::sha256_hex(text.as_bytes()),
            &StepId::new("chained-sew-source").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let original =
            DirectSolid::from_opened(Arc::new(opened), 0, &mut Decider::default()).unwrap();
        let source_faces = original.source_faces().to_vec();
        let left = original
            .faces()
            .iter()
            .find(|face| {
                matches!(&face.surface,
            Surface::Plane(p) if p.origin[0].is_zero() && !p.axis[0].is_zero())
            })
            .unwrap()
            .name
            .clone();
        let right = original
            .faces()
            .iter()
            .find(|face| {
                matches!(&face.surface,
            Surface::Plane(p) if p.origin[0] == r("2") && !p.axis[0].is_zero())
            })
            .unwrap()
            .name
            .clone();
        let translation = [r("1"), r("0"), r("0")];
        let mut proposal = original.faces().to_vec();
        let construction =
            support::move_face(&mut proposal, &right, &translation, &mut Decider::default())
                .unwrap();
        // Initial import already heals a gapped STEP fixture. Instead give the
        // real publication boundary a bounded, disconnected first-edit proposal;
        // admission and sewing are production code, not substituted evidence.
        let imperfect = proposal.iter_mut().find(|face| face.name == left).unwrap();
        for boundary in &mut imperfect.loops {
            for edge in &mut boundary.edges {
                edge.edge = format!("proposal:{}", edge.edge);
                for vertex in [&mut edge.start, &mut edge.end] {
                    vertex.name = format!("proposal:{}", vertex.name);
                    vertex.point[2] = &vertex.point[2] + &r("0.001");
                }
                let Curve::Line(origin, _) = &mut edge.curve else {
                    panic!("box line");
                };
                origin[2] = &origin[2] + &r("0.001");
            }
        }
        let first = original
            .admit_edit(
                Edit::MoveFace {
                    face: right.clone(),
                    translation_mm: translation.clone(),
                },
                proposal,
                construction,
                Default::default(),
                Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
        assert!(first.admission().unwrap().healing().report().changes.iter().any(|change| matches!(change,
            HealingChange::VerticesJoined { distance_upper_mm, .. } if r(distance_upper_mm) > r("0"))));
        assert!(first.history()[0].modified_faces.contains(&left));
        let retained_left = first
            .faces()
            .iter()
            .find(|face| face.name == left)
            .unwrap()
            .clone();
        let second = first
            .apply(
                Edit::MoveFace {
                    face: right,
                    translation_mm: translation,
                },
                ToleranceOptions {
                    initial_mm: r("0.000000001"),
                    maximum_mm: r("0.000001"),
                },
                n::NumericalOptions {
                    tolerance_mm: r("0.000001"),
                    vertex_tolerance_mm: None,
                },
                &mut Decider::default(),
            )
            .unwrap();
        assert_eq!(
            second
                .faces()
                .iter()
                .find(|face| face.name == left)
                .unwrap(),
            &retained_left
        );
        assert!(!second.history()[1].modified_faces.contains(&left));
        assert_eq!(max_x(&first), "3");
        assert_eq!(max_x(&second), "4");
        assert_eq!(original.source_faces(), source_faces);
        assert_eq!(first.source_faces(), second.source_faces());
    });
    assert!(!work.exhausted);
}
