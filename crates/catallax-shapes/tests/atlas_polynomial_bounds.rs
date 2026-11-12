use std::num::NonZeroUsize;

use catallax_numbers::{Decider, Rational, Real, Refusal, with_operation_budget};
use catallax_shapes::{BSplineSurface, Point3, RationalAtlas, SplineBoundary};

fn q(n: i64, d: i64) -> Real {
    Rational::from_ratio(n, d).unwrap().into()
}

fn cylinder(gap: bool, dec: &mut Decider) -> Result<BSplineSurface, Refusal> {
    let mut scale = Real::one();
    let mut controls = Vec::with_capacity(18);
    let mut weights = Vec::with_capacity(18);
    for (i, (x, y)) in [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
        (1, 0),
    ]
    .into_iter()
    .enumerate()
    {
        let weight = if i % 2 == 0 {
            scale.clone()
        } else {
            &scale * &q(7, 10)
        };
        let x = if gap && i == 8 {
            &Real::from_i64(x) + &q(1, 1_000_000)
        } else {
            x.into()
        };
        let y = Real::from_i64(y);
        for axial in [-17, 17] {
            controls.push(Point3::new(
                &q(1, 2) * &x,
                &q(1, 2) * &y,
                Real::from_i64(axial),
            ));
            weights.push(weight.clone());
        }
        if i % 2 == 1 {
            scale = &scale * &q(49, 50);
        }
    }
    BSplineSurface::new(
        [2, 1],
        [
            [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
                .map(Rational::from_i64)
                .to_vec(),
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        ],
        [9, 2],
        controls,
        weights,
        dec,
    )
}

#[test]
fn legacy_cylinder_admission_keeps_its_exact_work_boundary() {
    let (result, work) = with_operation_budget(71_656, || {
        let mut dec = Decider::default();
        let source = cylinder(false, &mut dec)?;
        RationalAtlas::partition(source, [true, false], &mut dec)
    });
    let atlases = result.unwrap();
    assert_eq!(atlases.len(), 1);
    assert!(atlases[0].periodic_u() && !atlases[0].periodic_v());
    assert!(!work.exhausted);
    assert_eq!(work.used, 71_656);
}

fn topology(atlas: &RationalAtlas) -> String {
    let mut record = String::new();
    for chart in atlas.charts() {
        let domain = chart.patch().domain();
        let neighbors = [
            SplineBoundary::UMin,
            SplineBoundary::UMax,
            SplineBoundary::VMin,
            SplineBoundary::VMax,
        ]
        .map(|side| chart.neighbor(side));
        assert!(neighbors[0].is_some() && neighbors[1].is_some());
        assert_eq!(neighbors[2].is_none(), domain[1][0] == Rational::zero());
        assert_eq!(neighbors[3].is_none(), domain[1][1] == Rational::one());
        if domain[0][0] == Rational::zero() {
            let across = atlas.charts()[neighbors[0].unwrap()].patch().domain();
            assert_eq!(across[0][1], Rational::from_i64(4));
            assert_eq!(across[1], domain[1]);
        }
        if domain[0][1] == Rational::from_i64(4) {
            let across = atlas.charts()[neighbors[1].unwrap()].patch().domain();
            assert_eq!(across[0][0], Rational::zero());
            assert_eq!(across[1], domain[1]);
        }
        record.push_str(&format!("{domain:?}:{neighbors:?};"));
    }
    record
}

#[test]
fn opt_in_periodic_topology_and_open_seam_refusals_are_thread_invariant() {
    let mut baseline = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let records: Vec<_> = [false, true]
            .into_iter()
            .map(|gap| {
                let (result, work) = with_operation_budget(20_000_000, || {
                    let mut dec = Decider::default();
                    let source = cylinder(gap, &mut dec)?;
                    let atlases = RationalAtlas::partition_with_polynomial_bounds(
                        source,
                        [true, false],
                        &mut dec,
                    )?;
                    assert_eq!(atlases.len(), 1);
                    Ok::<_, Refusal>(topology(&atlases[0]))
                });
                assert!(!work.exhausted);
                if gap {
                    assert_eq!(
                        result.as_ref().unwrap_err().code,
                        "shapes.rational-atlas-open"
                    );
                } else {
                    result.as_ref().unwrap();
                }
                (format!("{result:?}"), work.used)
            })
            .collect();
        println!(
            "polynomial admission threads={threads:?}, work={:?}",
            records.iter().map(|(_, work)| *work).collect::<Vec<_>>()
        );
        if let Some(baseline) = &baseline {
            assert_eq!(baseline, &records);
        } else {
            baseline = Some(records);
        }
    }
    catallax_numbers::set_threads(None);
}

#[test]
fn nonrational_coefficients_refuse_only_the_opt_in_constructor() {
    let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
        let mut dec = Decider::default();
        let height = Real::from_i64(2).sqrt();
        let source = BSplineSurface::new(
            [1, 1],
            std::array::from_fn(|_| [0, 0, 1, 1].map(Rational::from_i64).to_vec()),
            [2, 2],
            [[0, 0], [0, 1], [1, 0], [1, 1]]
                .map(|[x, y]| Point3::new(x.into(), y.into(), height.clone()))
                .to_vec(),
            vec![Real::one(); 4],
            &mut dec,
        )?;
        let error =
            RationalAtlas::partition_with_polynomial_bounds(source.clone(), [false; 2], &mut dec)
                .unwrap_err();
        assert_eq!(error.code, "shapes.rational-atlas-polynomial-coefficients");
        let legacy = RationalAtlas::partition(source, [false; 2], &mut dec)?;
        assert_eq!(legacy.len(), 1);
        Ok(())
    });
    assert!(!work.exhausted);
    result.unwrap();
}
