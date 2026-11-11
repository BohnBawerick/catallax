//! Exact finite rolling-region controls, distinct from recipe acceptance.
use catallax_field::rolling_region::{Region, Relation};
use catallax_numbers::{
    Decider, OperationAllowance, Outcome, Rational, Real, Sign, with_operation_budget,
};
use catallax_shapes::{Point3, RollingEnvelopeDefinition, Vector3};
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn r(s: &str) -> Real {
    q(s).into()
}
fn p(x: &str, y: &str, z: &str) -> Point3 {
    Point3::new(r(x), r(y), r(z))
}
fn v(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(x.into(), y.into(), z.into())
}
fn region(radius: &str, placement: usize) -> Region {
    let (origin, axis, toward) = match placement {
        0 => (p("0", "5", "0"), v(0, 1, 0), v(1, 0, 0)),
        1 => (p("7", "16", "-13"), v(0, 1, 0), v(1, 0, 0)),
        2 => (p("-5", "0", "0"), v(-1, 0, 0), v(0, 1, 0)),
        _ => unreachable!(),
    };
    let mut dec = Decider::default();
    Region::new(
        RollingEnvelopeDefinition::new(
            origin,
            axis,
            toward,
            r("12"),
            r("4"),
            r("17"),
            r(radius),
            [r("1"), r("1")],
            &mut dec,
        )
        .unwrap(),
        &mut dec,
    )
    .unwrap()
}
fn place(p: Point3, i: usize) -> Point3 {
    match i {
        0 => p,
        1 => p.plus(&v(7, 11, -13)),
        2 => Point3::new(-p.y, p.x, p.z),
        _ => unreachable!(),
    }
}
fn zero(value: Real, dec: &mut Decider) {
    assert!(
        matches!(dec.sign(&value), Outcome::ProvenCoincident(_)),
        "{value:?}"
    );
}
#[test]
fn exact_contacts_and_both_transitions() {
    for radius in ["2", "5/2"] {
        for placement in 0..3 {
            let region = region(radius, placement);
            let d = region.definition();
            let mut dec = Decider::default();
            for u in ["-1", "-1/2", "0", "1/2", "1"] {
                let station = region.station(&q(u), &mut dec).unwrap();
                for point in [&station.cylinder, &station.mouth] {
                    zero(
                        point.minus(&station.centre).norm_squared() - d.radius().square(),
                        &mut dec,
                    );
                }
                let a = station.cylinder.minus(d.origin());
                let axial = a.dot(d.axis());
                zero(
                    a.norm_squared() - axial.square() - d.body().square(),
                    &mut dec,
                );
                let mouth = d.origin().plus(&d.toward().scaled(d.offset()));
                zero(
                    station.mouth.minus(&mouth).norm_squared() - d.mouth().square(),
                    &mut dec,
                );
                zero(station.mouth.minus(d.origin()).dot(d.axis()), &mut dec);
            }
            for (u, sign) in [("-1", Sign::Negative), ("1", Sign::Positive)] {
                let station = region.station(&q(u), &mut dec).unwrap();
                let [cylinder, mouth] = d.corners(sign);
                zero(station.cylinder.minus(&cylinder).norm_squared(), &mut dec);
                zero(station.mouth.minus(&mouth).norm_squared(), &mut dec);
                zero(
                    station.centre.minus(d.origin()).dot(d.axis()) - d.radius().clone(),
                    &mut dec,
                );
            }
            assert!(region.station(&q("1001/1000"), &mut dec).is_err());
        }
    }
}
#[test]
fn finite_material_side_is_not_the_full_torus() {
    for placement in 0..3 {
        let region = region("2", placement);
        let mut dec = Decider::default();
        for (point, want) in [
            (p("12.2", "5.2", "0"), Relation::Interior),
            (p("12.5", "5.5", "0"), Relation::Exterior),
            (p("14", "6", "0"), Relation::Exterior),
            (p("12", "4", "0"), Relation::Exterior),
            (p("12", "8", "0"), Relation::Exterior),
            (p("0", "6", "12"), Relation::Exterior),
            (p("-12.2", "5.2", "0"), Relation::Exterior),
        ] {
            assert_eq!(
                region.classify(&place(point, placement), &mut dec).unwrap(),
                want
            );
        }
    }
    // At (12.5,5.5,0), the full torus adds material: distance^2=4.5>r^2=4.
    // The native middle ball has centre (14,5+sqrt(3),0), and the point is
    // strictly inside it: 2.25+(sqrt(3)-.5)^2=5.5-sqrt(3)<4.
    let mut dec = Decider::default();
    assert!(matches!(
        dec.sign(&(r("3").sqrt() - r("3/2"))),
        Outcome::Proven(_)
    ));
}
#[test]
fn boundary_and_normal_cone_remain_explicit() {
    let region = region("2", 0);
    let mut dec = Decider::default();
    // Clipping-wall incidence is not promoted into a boundary certificate.
    assert_eq!(
        region
            .classify(&p("12.2", "5", "0"), &mut dec)
            .unwrap_err()
            .code,
        "field.rolling-region-unresolved"
    );
    let def = RollingEnvelopeDefinition::new(
        p("0", "5", "0"),
        v(0, 1, 0),
        v(1, 0, 0),
        r("12"),
        r("4"),
        r("17"),
        r("7/2"),
        [r("1"), r("1")],
        &mut dec,
    );
    assert_eq!(def.unwrap_err().code, "shapes.rolling-envelope");
}
#[test]
fn production_batch_full_semantics_and_allowances() {
    let region = region("2", 0);
    let points = [
        p("12.5", "5.4", "0"),
        p("12.5", "5.5", "0"),
        p("14", "6", "0"),
        p("0", "6", "12"),
        p("12.2", "5", "0"),
        p("12", "8", "0"),
        p("-12.2", "5.2", "0"),
    ];
    for cap in [700_000, 7, 0] {
        let mut expected = None;
        for workers in [Some(1), Some(2), Some(3), Some(4), Some(9), None] {
            let (result, work) = with_operation_budget(1_000_000, || {
                region
                    .classify_batch(
                        &points,
                        workers,
                        OperationAllowance::reserve(cap).unwrap(),
                        &mut Decider::default(),
                    )
                    .unwrap()
            });
            assert_eq!(
                result.workers,
                workers
                    .unwrap_or_else(|| std::thread::available_parallelism().unwrap().get())
                    .min(points.len())
            );
            assert_eq!(result.records.len(), points.len());
            assert_eq!(work.used, cap + 2);
            assert!(!work.exhausted);
            if cap == 0 {
                assert!(
                    result
                        .records
                        .iter()
                        .all(|r| r.result.is_err() && r.exhausted)
                );
            }
            if let Some(expected) = &expected {
                assert_eq!(&result.records, expected);
            } else {
                expected = Some(result.records.clone());
            }
            println!(
                "cap={cap} requested={workers:?} spawned={} work={:?} exhausted={:?} outcomes={:?}",
                result.workers,
                result.records.iter().map(|r| r.work).collect::<Vec<_>>(),
                result
                    .records
                    .iter()
                    .map(|r| r.exhausted)
                    .collect::<Vec<_>>(),
                result.records.iter().map(|r| &r.result).collect::<Vec<_>>()
            );
        }
    }
    let mut dec = Decider::default();
    assert!(
        region
            .classify_batch(
                &[],
                Some(0),
                OperationAllowance::reserve(0).unwrap(),
                &mut dec
            )
            .is_err()
    );
    assert_eq!(
        region
            .classify_batch(&[], None, OperationAllowance::reserve(0).unwrap(), &mut dec)
            .unwrap()
            .workers,
        0
    );
    let (result, work) = with_operation_budget(10, || OperationAllowance::reserve(11));
    assert!(result.is_err());
    assert!(work.exhausted);
}

#[test]
fn caller_deadline_is_not_dropped() {
    let region = region("2", 0);
    let mut dec = Decider::default();
    let deadline = std::time::Instant::now();
    dec.set_deadline(Some(deadline));
    assert_eq!(dec.deadline(), Some(deadline));
    let scalar = region
        .classify(&p("12.2", "5.2", "0"), &mut dec)
        .unwrap_err();
    let batch = region
        .classify_batch(
            &[p("12.2", "5.2", "0")],
            None,
            OperationAllowance::reserve(1000).unwrap(),
            &mut dec,
        )
        .unwrap_err();
    assert_eq!(scalar, batch);
}

#[test]
fn analytic_middle_envelope_and_precision_refusal() {
    let region = region("2", 0);
    let mut dec = Decider::default();
    // On the symmetry meridian, both centre height and planar distance from
    // a positive-x query increase with |u|. Thus u=0 gives the lowest graph.
    assert_eq!(
        region.classify(&p("12.5", "5.4", "0"), &mut dec).unwrap(),
        Relation::Interior
    );
    assert_eq!(
        region.classify(&p("12.5", "5.42", "0"), &mut dec).unwrap(),
        Relation::Exterior
    );
    let boundary = Point3::new(
        r("12.5"),
        r("5") + r("3").sqrt() - r("7").sqrt() / r("2"),
        r("0"),
    );
    assert_eq!(
        region.classify(&boundary, &mut dec).unwrap_err().code,
        "field.rolling-region-unresolved"
    );
    let mut limited = Decider::new(catallax_numbers::BudgetPolicy {
        max_precision_bits: 64,
        max_work: 1_000_000,
    });
    let batch = region
        .classify_batch(
            &[p("12.5", "5.4", "0")],
            Some(1),
            OperationAllowance::reserve(100_000).unwrap(),
            &mut limited,
        )
        .unwrap();
    assert_eq!(
        batch.records[0].result.as_ref().unwrap_err().code,
        "numbers.precision-limit"
    );
}
