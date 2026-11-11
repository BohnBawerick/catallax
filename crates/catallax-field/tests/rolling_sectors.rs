//! Angular coverage is an explicit geometric prerequisite, not P5 admission.
use catallax_field::rolling_region::{Region, Sector, Sectors, Transition};
use catallax_numbers::{Decider, Outcome, Real, with_operation_budget};
use catallax_shapes::{Point3, RollingEnvelopeDefinition, Vector3};

fn v(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(x.into(), y.into(), z.into())
}
fn region(toward: Vector3, moved: bool) -> Region {
    let mut dec = Decider::default();
    Region::new(
        RollingEnvelopeDefinition::new(
            Point3::new(if moved { 7 } else { 0 }.into(), 5.into(), 0.into()),
            v(0, 1, 0),
            toward,
            12.into(),
            4.into(),
            17.into(),
            2.into(),
            [1.into(), 1.into()],
            &mut dec,
        )
        .unwrap(),
        &mut dec,
    )
    .unwrap()
}
fn four(moved: bool) -> Sectors {
    Sectors::new(
        [v(1, 0, 0), v(0, 0, 1), v(-1, 0, 0), v(0, 0, -1)]
            .into_iter()
            .map(|t| region(t, moved))
            .collect(),
        &mut Decider::default(),
    )
    .unwrap()
}
fn zero(value: Real, dec: &mut Decider) {
    assert!(matches!(dec.sign(&value), Outcome::ProvenCoincident(_)));
}

#[test]
fn p5_four_sectors_and_exact_endpoint_contacts() {
    // Bind this control to the authored radius, not the admitted 1/2 edit.
    let source: serde_json::Value =
        serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    assert_eq!(
        source["steps"]
            .as_array()
            .unwrap()
            .iter()
            .find(|s| s["id"] == "filleted")
            .unwrap()["radius"],
        "2"
    );
    for moved in [false, true] {
        let sectors = four(moved);
        let mut dec = Decider::default();
        for (index, region) in sectors.regions().iter().enumerate() {
            let d = region.definition();
            assert_eq!(
                sectors.locate(d.toward(), &mut dec).unwrap(),
                Sector::Mouth(index)
            );
            for end in [Transition::Start, Transition::End] {
                let s = sectors.transition(index, end, &mut dec).unwrap();
                let delta = s.centre.minus(d.origin());
                let height = delta.dot(d.axis());
                let radial = delta.minus(&d.axis().scaled(&height));
                assert_eq!(
                    sectors.locate(&radial, &mut dec).unwrap(),
                    Sector::Transition(index, end)
                );
                // Exact directions on each side of the transition. The small
                // rational displacement is a test input, not a predicate tolerance.
                let epsilon: Real = catallax_numbers::Rational::parse("1/1000000")
                    .unwrap()
                    .into();
                let outward = d.across().scaled(&match end {
                    Transition::Start => -epsilon,
                    Transition::End => epsilon,
                });
                assert_eq!(
                    sectors.locate(&radial.plus(&outward), &mut dec).unwrap(),
                    Sector::Torus
                );
                assert_eq!(
                    sectors.locate(&radial.minus(&outward), &mut dec).unwrap(),
                    Sector::Mouth(index)
                );
                // Ordinary torus centre circle, plane and cylinder contacts.
                zero(height - d.radius().clone(), &mut dec);
                zero(
                    radial.norm_squared() - (d.body() + d.radius()).square(),
                    &mut dec,
                );
                zero(
                    s.mouth
                        .minus(&s.centre.plus(&d.axis().scaled(&(-d.radius().clone()))))
                        .norm_squared(),
                    &mut dec,
                );
                let torus_cylinder = s
                    .centre
                    .plus(&radial.scaled(&(-d.radius().clone() / (d.body() + d.radius()))));
                zero(s.cylinder.minus(&torus_cylinder).norm_squared(), &mut dec);
            }
        }
        for direction in [v(1, 0, 1), v(-1, 0, 1), v(-1, 0, -1), v(1, 0, -1)] {
            assert_eq!(sectors.locate(&direction, &mut dec).unwrap(), Sector::Torus);
            assert_eq!(
                sectors
                    .locate(&direction.scaled(&7.into()), &mut dec)
                    .unwrap(),
                Sector::Torus
            );
        }
    }
}

#[test]
fn closed_sector_touch_overlap_and_incompatible_ring_refuse() {
    let first = region(v(1, 0, 0), false);
    let mut dec = Decider::default();
    let station = first
        .station(&catallax_numbers::Rational::one(), &mut dec)
        .unwrap();
    let delta = station.centre.minus(first.definition().origin());
    let c = delta.x / Real::from_i64(14);
    let s = delta.z / Real::from_i64(14);
    // Rotate by precisely twice the half-angle: the two closed sectors touch.
    let touching = Vector3::new(
        &c.square() - &s.square(),
        Real::zero(),
        Real::from_i64(2) * c * s,
    );
    for regions in [
        vec![first.clone(), region(touching, false)],
        vec![first.clone(), first.clone()],
        vec![first.clone(), region(v(0, 0, 1), true)],
        vec![],
    ] {
        assert_eq!(
            Sectors::new(regions, &mut dec).unwrap_err().code,
            "field.rolling-sectors-unresolved"
        );
    }
    let sectors = four(false);
    for bad in [v(0, 0, 0), v(1, 1, 0)] {
        assert_eq!(
            sectors.locate(&bad, &mut dec).unwrap_err().code,
            "field.rolling-sectors-unresolved"
        );
    }
    assert_eq!(
        sectors
            .transition(4, Transition::Start, &mut dec)
            .unwrap_err()
            .code,
        "field.rolling-sectors-unresolved"
    );
}

#[test]
fn cancellation_and_exhaustion_preserve_refusals() {
    let sectors = four(false);
    let mut dec = Decider::default();
    dec.set_deadline(Some(std::time::Instant::now()));
    let expected = dec.check_deadline().unwrap_err();
    assert_eq!(sectors.locate(&v(1, 0, 0), &mut dec).unwrap_err(), expected);
    assert_eq!(
        sectors
            .transition(0, Transition::Start, &mut dec)
            .unwrap_err(),
        expected
    );
    assert_eq!(
        Sectors::new(sectors.regions().to_vec(), &mut dec).unwrap_err(),
        expected
    );
    let (result, work) =
        with_operation_budget(0, || sectors.locate(&v(1, 0, 0), &mut Decider::default()));
    assert!(result.is_err());
    assert!(work.exhausted);
}
