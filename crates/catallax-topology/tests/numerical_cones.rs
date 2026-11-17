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
fn native_cone_trims_keep_complete_positive_height_and_periodic_domains() {
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
    let body = n::build_with_native_trims(
        positive.clone(),
        vec![native_seam(&positive)],
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (7, 15, 10));
    let error = n::build_with_native_trims(
        negative.clone(), vec![native_seam(&negative)],
        n::NumericalOptions::default(), &mut Decider::default(),
    ).unwrap_err();
    assert_eq!(error.code, "numerical.native-cone-domain-unverified");

    // Positive endpoints do not authorize an interior apex crossing or a
    // second angular turn hidden in the native parameter-curve hull.
    for apex in [true, false] {
        let mut binding = native_seam(&positive);
        let controls: Vec<_> = binding.curve.controls().cloned().collect();
        let mut middle = controls[0].clone();
        if apex {
            middle.y = Real::from_i64(-8);
        } else {
            middle.x = (catallax_numbers::Ball::pi(128).mid().mul_pow2(2)).into();
        }
        let [lo, hi] = binding.curve.domain();
        binding.curve = BSplineCurve::new(
            2, vec![lo.clone(), lo.clone(), lo.clone(), hi.clone(), hi.clone(), hi.clone()],
            vec![controls[0].clone(), middle, controls[1].clone()],
            vec![Real::one(); 3], &mut Decider::default(),
        ).unwrap();
        let error = n::build_with_native_trims(
            positive.clone(), vec![binding],
            n::NumericalOptions::default(), &mut Decider::default(),
        ).unwrap_err();
        assert_eq!(error.code, "numerical.native-cone-domain-unverified");
    }
}
