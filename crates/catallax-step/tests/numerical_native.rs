//! Analytic export must not discard retained native parameter curves.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;
#[path = "../../catallax-topology/tests/common/periodic.rs"]
mod periodic;

use catallax_numbers::{Decider, Real, with_operation_budget};
use catallax_shapes::{BSplineCurve, Point3};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;
use fixture::{q, quarter};

#[test]
fn an_analytic_body_with_native_uv_uses_the_retaining_writer() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let definition = quarter();
        let face = &definition.faces[4];
        let trim = &face.loops[0][0];
        let edge = &definition.edges[trim.edge];
        let mut domain = edge.parameters.clone();
        domain.sort();
        let curve = BSplineCurve::new(
            1,
            vec![
                domain[0].clone(),
                domain[0].clone(),
                domain[1].clone(),
                domain[1].clone(),
            ],
            domain
                .iter()
                .map(|t| Point3::new(t.clone().into(), 0.into(), 0.into()))
                .collect(),
            vec![Real::one(); 2],
            &mut Decider::default(),
        )
        .unwrap();
        let binding = n::NumericalNativeTrim {
            face: face.key.clone(),
            boundary: 0,
            use_index: 0,
            edge: edge.key.clone(),
            curve,
        };
        let body = n::build_with_native_trims(
            definition,
            vec![binding],
            n::NumericalOptions {
                tolerance_mm: q("1e-7"),
            },
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(body.native_trims().len(), 1);
        // The retaining writer currently requires rational patch/plane faces.
        // Refuse this cylinder instead of routing it through an empty-sidecar path.
        let refusal = match output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        ) {
            Ok(file) => panic!(
                "analytic writer discarded native trims: source={}, emitted={}",
                file.report().source_native_trims.len(),
                file.report().emitted_native_trims.len()
            ),
            Err(refusal) => refusal,
        };
        assert_eq!(refusal.code, "step.numerical-unsupported");
    });
    assert!(!work.exhausted, "{work:?}");
}

#[test]
fn native_outer_edges_keep_continuous_circle_clearance() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let definition = periodic::bore_with_radius("399/200");
        let face = &definition.faces[0];
        let trim = &face.loops[0][0];
        let edge = &definition.edges[trim.edge];
        let mut endpoints = trim.uv.clone();
        if trim.sense == catallax_topology::Sense::Reversed {
            endpoints.reverse();
        }
        let mut domain = edge.parameters.clone();
        if domain[0] > domain[1] {
            domain.reverse();
            endpoints.reverse();
        }
        let curve = BSplineCurve::new(
            1,
            vec![
                domain[0].clone(),
                domain[0].clone(),
                domain[1].clone(),
                domain[1].clone(),
            ],
            endpoints
                .into_iter()
                .map(|[u, v]| Point3::new(u.into(), v.into(), 0.into()))
                .collect(),
            vec![Real::one(); 2],
            &mut Decider::default(),
        )
        .unwrap();
        let binding = n::NumericalNativeTrim {
            face: face.key.clone(),
            boundary: 0,
            use_index: 0,
            edge: edge.key.clone(),
            curve,
        };
        let bad = n::build_with_native_trims(
            definition.clone(),
            vec![binding.clone()],
            n::NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(bad.code, "numerical.invalid-trim");
        assert!(bad.detail.contains("continuous clearance"));
        let body = n::build_with_native_trims(
            definition,
            vec![binding],
            n::NumericalOptions {
                tolerance_mm: q("1/1000"),
            },
            &mut Decider::default(),
        )
        .unwrap();
        assert!(body.evidence().checks().iter().any(|check| matches!(check,
            catallax_numbers::NumericalCheck::Passed { name,.. }
            if name=="trim.continuous-circle-hole-clearance")));
    });
    assert!(!work.exhausted, "{work:?}");
}
