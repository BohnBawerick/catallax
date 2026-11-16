//! Cone chart constraints through the public native-trim build entry.
#[allow(dead_code)]
#[path = "common/cones.rs"]
mod fixture;

use catallax_numbers::{Decider, Real};
use catallax_shapes::{BSplineCurve, Curve, Point3, Surface};
use catallax_topology::{Sense, numerical as n};

fn native_seam(d: &n::NumericalDefinition) -> n::NumericalNativeTrim {
    let face = d
        .faces
        .iter()
        .find(|f| matches!(f.surface, Surface::Cone(_)))
        .unwrap();
    let (use_index, trim) = face.loops[0]
        .iter()
        .enumerate()
        .find(|(_, t)| matches!(d.edges[t.edge].curve, Curve::Line(_)))
        .unwrap();
    let edge = &d.edges[trim.edge];
    let mut ends = trim.uv.clone();
    if trim.sense == Sense::Reversed {
        ends.reverse();
    }
    let mut interval = edge.parameters.clone();
    if interval[0] > interval[1] {
        interval.reverse();
        ends.reverse();
    }
    let curve = BSplineCurve::new(
        1,
        vec![
            interval[0].clone(),
            interval[0].clone(),
            interval[1].clone(),
            interval[1].clone(),
        ],
        ends.map(|[u, v]| Point3::new(u.into(), v.into(), Real::zero()))
            .to_vec(),
        vec![Real::one(); 2],
        &mut Decider::default(),
    )
    .unwrap();
    n::NumericalNativeTrim {
        face: face.key.clone(),
        boundary: 0,
        use_index,
        edge: edge.key.clone(),
        curve,
    }
}

#[test]
fn reproduce_native_cone_height_guard_bypass() {
    let positive = fixture::positive_bore();
    n::build(
        positive.clone(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let negative = fixture::negative_height_bore();
    let error = n::build(
        negative.clone(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "numerical.invalid-trim");
    assert!(error.detail.contains("above its apex"));
    for d in [positive, negative] {
        let binding = native_seam(&d);
        let body = n::build_with_native_trims(
            d,
            vec![binding],
            n::NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(body.native_trims().len(), 1);
        assert!(
            body.definition()
                .faces
                .iter()
                .any(|f| matches!(f.surface, Surface::Cone(_)))
        );
    }
}
