//! Exact chart and periodic ownership through the existing cell-pair engine.

use catallax_numbers::{Decider, Rational, Real, Refusal, with_operation_budget};
use catallax_shapes::{BSplineSurface, Point3, RationalAtlas};
use std::num::NonZeroUsize;

fn q(n: i64, d: i64) -> Real {
    Real::from_rational(Rational::from_ratio(n, d).unwrap())
}

fn cylinder(
    periodic: bool,
    transpose: bool,
    gap: bool,
    dec: &mut Decider,
) -> Result<RationalAtlas, Refusal> {
    let circle = [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
        (1, 0),
    ];
    let mut rows = Vec::new();
    let mut scale = Real::one();
    for (i, (x, y)) in circle.into_iter().enumerate() {
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
        rows.push((x, Real::from_i64(y), weight));
        if i % 2 == 1 {
            scale = &scale * &q(49, 50);
        }
    }
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    for i in 0..18 {
        let (circle, axial) = if transpose {
            (i % 9, i / 9)
        } else {
            (i / 2, i % 2)
        };
        let (x, y, w) = &rows[circle];
        controls.push(Point3::new(
            x.clone(),
            y.clone(),
            Real::from_i64(if axial == 0 { -1 } else { 1 }),
        ));
        weights.push(w.clone());
    }
    let circle_knots = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
        .map(Rational::from_i64)
        .to_vec();
    let linear_knots = [0, 0, 1, 1].map(Rational::from_i64).to_vec();
    let (degree, count, knots, periodic) = if transpose {
        (
            [1, 2],
            [2, 9],
            [linear_knots, circle_knots],
            [false, periodic],
        )
    } else {
        (
            [2, 1],
            [9, 2],
            [circle_knots, linear_knots],
            [periodic, false],
        )
    };
    let source = BSplineSurface::new(degree, knots, count, controls, weights, dec)?;
    let mut atlases = RationalAtlas::partition(source, periodic, dec)?;
    assert_eq!(atlases.len(), 1);
    Ok(atlases.pop().unwrap())
}

fn plane(dec: &mut Decider) -> RationalAtlas {
    let source = BSplineSurface::new(
        [1, 1],
        std::array::from_fn(|_| [0, 0, 1, 1].map(Rational::from_i64).to_vec()),
        [2, 2],
        [[-2, -2], [-2, 2], [2, -2], [2, 2]]
            .map(|p| Point3::new(p[0].into(), p[1].into(), q(1, 5)))
            .to_vec(),
        vec![Real::one(); 4],
        dec,
    )
    .unwrap();
    RationalAtlas::partition(source, [false, false], dec)
        .unwrap()
        .pop()
        .unwrap()
}

fn run(periodic: bool, transpose: bool, reversed: bool) -> (String, u64) {
    let (result, work) = with_operation_budget(20_000_000, || -> Result<String, Refusal> {
        let mut dec = Decider::default();
        let cylinder = cylinder(periodic, transpose, false, &mut dec)?;
        let plane = plane(&mut dec);
        let found = if reversed {
            plane.patch_meeting(&cylinder, &mut dec)?
        } else {
            cylinder.patch_meeting(&plane, &mut dec)?
        };
        assert_eq!(
            found.branches().len(),
            1,
            "chart cuts cannot fragment the circle"
        );
        let branch = &found.branches()[0];
        assert_eq!(
            branch.is_closed(),
            periodic,
            "only a proved periodic identification closes the parameter graph"
        );
        let mut summary = format!("closed={};", branch.is_closed());
        for arc in branch.arcs() {
            summary.push_str(&format!(
                "{:?};",
                (
                    arc.first_chart,
                    arc.first_span,
                    &arc.first_cell,
                    arc.second_chart,
                    arc.second_span,
                    &arc.second_cell,
                    arc.certificate
                )
            ));
            for end in &arc.ends {
                summary.push_str(&format!("{:?};", end.enclosure()));
                let point = found.point(end);
                let on_cylinder = &point.x.square() + &point.y.square() - Real::one();
                let on_plane = &point.z - &q(1, 5);
                for equation in [on_cylinder, on_plane] {
                    let bound = dec.enclose_budgeted_at(&equation, 128)?;
                    assert!(
                        bound.contains_zero(),
                        "retained endpoint leaves an original support: {bound:?}"
                    );
                    assert!(bound.hi() - bound.lo() < Rational::one().mul_pow2(-80));
                }
            }
        }
        Ok(summary)
    });
    assert!(!work.exhausted);
    (result.unwrap(), work.used)
}

#[test]
fn periodic_rows_join_exactly_in_either_parameter_and_either_pair_order() {
    let mut baseline = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let records: Vec<_> = [false, true]
            .into_iter()
            .flat_map(|transpose| {
                [false, true]
                    .into_iter()
                    .map(move |reversed| run(true, transpose, reversed))
            })
            .collect();
        println!(
            "atlas pair threads={threads:?} work={:?}",
            records.iter().map(|(_, work)| work).collect::<Vec<_>>()
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
fn unjoined_rows_remain_physical_even_when_their_world_points_agree() {
    for transpose in [false, true] {
        run(false, transpose, false);
    }
}

#[test]
fn a_nearly_closing_row_is_not_a_periodic_identification() {
    let (result, work) = with_operation_budget(20_000_000, || {
        cylinder(true, false, true, &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert_eq!(result.unwrap_err().code, "shapes.rational-atlas-open");
}
