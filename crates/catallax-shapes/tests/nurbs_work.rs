//! Construction must stop at the caller's deterministic work allowance.

use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::spline_tools::{RationalPoint, RationalSplineSurface};

#[test]
fn dense_grid_construction_refuses_when_its_rational_work_exhausts_the_allowance() {
    let nodes: Vec<_> = (0..4)
        .map(|i| Rational::from_ratio(i, 3).unwrap())
        .collect();
    let points: Vec<RationalPoint> = nodes
        .iter()
        .flat_map(|u| nodes.iter().map(move |v| [u.clone(), v.clone(), u * v]))
        .collect();
    // The smaller limit cannot admit even the dense matrix storage; the larger
    // one admits that storage but must stop the rational solve or tensor work.
    for allowance in [16, 512] {
        let (result, work) = with_operation_budget(allowance, || {
            RationalSplineSurface::through_grid(
                [&nodes, &nodes],
                [3, 3],
                &points,
                &mut Decider::default(),
            )
        });
        assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
        assert!(work.exhausted);
        assert!(work.used <= allowance);
    }
    let (surface, work) = with_operation_budget(100_000, || {
        RationalSplineSurface::through_grid(
            [&nodes, &nodes],
            [3, 3],
            &points,
            &mut Decider::default(),
        )
    });
    let surface = surface.unwrap();
    assert!(!work.exhausted);
    let middle = Rational::from_ratio(1, 2).unwrap();
    assert_eq!(
        surface
            .point_at(&middle, &middle, &Decider::default())
            .unwrap(),
        [middle.clone(), middle, Rational::from_ratio(1, 4).unwrap()]
    );
}

#[test]
fn prescribed_grid_basis_refuses_nonclamped_zero_width_endpoint_spans() {
    let q = |text: &str| Rational::parse(text).unwrap();
    let u = ["0", "1/3", "2/3", "1"].map(q);
    let v = ["0", "1"].map(q);
    let points: Vec<_> = u
        .iter()
        .flat_map(|u| v.iter().map(move |v| [u.clone(), v.clone(), u * v]))
        .collect();
    let result = RationalSplineSurface::through_grid_on_knots(
        [&u, &v],
        [2, 1],
        [
            ["0", "0", "0", "1", "1", "2", "2"].map(q).to_vec(),
            ["0", "0", "1", "1"].map(q).to_vec(),
        ],
        &points,
        &mut Decider::default(),
    );
    assert_eq!(result.unwrap_err().code, "shapes.spline-exchange-data");
    let surface = RationalSplineSurface::through_grid_on_knots(
        [&u, &v],
        [2, 1],
        [
            ["0", "0", "0", "1/2", "1", "1", "1"].map(q).to_vec(),
            ["0", "0", "1", "1"].map(q).to_vec(),
        ],
        &points,
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        surface
            .point_at(&q("1/2"), &q("1/2"), &Decider::default())
            .unwrap(),
        [q("1/2"), q("1/2"), q("1/4")]
    );
}
