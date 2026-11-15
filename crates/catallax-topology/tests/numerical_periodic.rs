//! Complete circular rims need an explicit paired seam before a cylinder chart closes.
#[path = "common/periodic.rs"]
mod fixture;

use catallax_numbers::Decider;
use catallax_topology::numerical::{self as n, NumericalOptions};

#[test]
fn complete_rims_and_one_paired_seam_admit_a_closed_bore() {
    let body = n::build(
        fixture::bore(),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (7, 15, 10));
    assert_eq!((counts.shells, counts.lumps), (1, 1));
    assert_eq!(body.evidence().tolerance_mm(), &fixture::q("1/100"));
}

#[test]
fn closed_partial_arcs_extra_turns_and_implicit_closure_refuse() {
    for kind in 0..3 {
        let mut definition = fixture::bore();
        let rim = definition
            .edges
            .iter()
            .position(|e| e.key == "rim-0")
            .unwrap();
        match kind {
            0 => definition.edges[rim].parameters[1] = fixture::q("1"),
            1 => {
                definition.edges[rim].parameters[1] =
                    &definition.edges[rim].parameters[1] * &fixture::q("2")
            }
            _ => {
                let mut coincident = definition.vertices[8].clone();
                coincident.key = "unsewn-coincident-anchor".into();
                definition.edges[rim].vertices[1] = definition.vertices.len();
                definition.vertices.push(coincident);
            }
        }
        assert!(
            n::build(
                definition,
                NumericalOptions::default(),
                &mut Decider::default()
            )
            .is_err()
        );
    }
}

#[test]
fn separate_seam_edges_wrong_winding_and_double_period_cells_refuse() {
    for kind in 0..3 {
        let mut definition = fixture::bore();
        let wall = definition.faces.len() - 1;
        match kind {
            0 => {
                let mut duplicate = definition.edges.last().unwrap().clone();
                duplicate.key = "unshared-seam".into();
                definition.faces[wall].loops[0][3].edge = definition.edges.len();
                definition.edges.push(duplicate);
            }
            1 => definition.faces[wall].loops[0][3].sense = catallax_topology::Sense::Same,
            _ => {
                for trim in &mut definition.faces[wall].loops[0] {
                    for uv in &mut trim.uv {
                        uv[0] = &uv[0] * &fixture::q("2");
                    }
                }
            }
        }
        assert!(
            n::build(
                definition,
                NumericalOptions::default(),
                &mut Decider::default()
            )
            .is_err()
        );
    }
}

#[test]
fn translated_seam_branches_may_slope_within_the_numerical_tolerance() {
    let mut definition = fixture::bore();
    for trim in definition
        .faces
        .last_mut()
        .unwrap()
        .loops
        .iter_mut()
        .flatten()
    {
        for uv in &mut trim.uv {
            uv[0] = &uv[0] + &(&uv[1] * &fixture::q("1/100000000"));
        }
    }
    n::build(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
}

#[test]
fn a_shared_seam_cannot_fall_back_to_a_near_turn_chart() {
    let mut definition = fixture::bore();
    for trim in definition
        .faces
        .last_mut()
        .unwrap()
        .loops
        .iter_mut()
        .flatten()
    {
        for uv in &mut trim.uv {
            if uv[0] > fixture::q("1") {
                uv[0] = &uv[0] - &fixture::q("1/1000");
            }
        }
    }
    let error = n::build(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "numerical.invalid-trim");
    assert!(error.detail.contains("paired seam"));
}

#[test]
fn continuous_circle_clearance_respects_the_requested_tolerance() {
    let definition = fixture::bore_with_radius("399/200");
    assert!(
        n::build(
            definition.clone(),
            NumericalOptions::default(),
            &mut Decider::default()
        )
        .is_err()
    );
    let body = n::build(
        definition,
        NumericalOptions {
            tolerance_mm: fixture::q("1/1000"),
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(body.evidence().checks().iter().any(|check| matches!(check,
        catallax_numbers::NumericalCheck::Passed { name,.. } if name=="trim.continuous-circle-hole-clearance")));
}
