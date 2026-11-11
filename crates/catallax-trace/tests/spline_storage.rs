//! Checks the owned-storage boundary that replaces C++ pointers and contexts.
use catallax_trace::{
    Fault,
    spline::{BSplineCurve, BSplineSurface, ControlPoint, NO_WEIGHT},
    vector::Vector3,
};

fn point(x: f64, y: f64, z: f64) -> ControlPoint {
    ControlPoint {
        coordinates: [x, y, z, NO_WEIGHT],
    }
}

#[test]
fn malformed_storage_is_refused_before_indexing() {
    let points = vec![point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0)];
    for knots in [
        vec![],
        vec![0.0, 1.0],
        vec![0.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 1.0],
        vec![0.0, 0.0, f64::NAN, 1.0],
    ] {
        assert_eq!(
            BSplineCurve::new(1, knots, points.clone()),
            Err(Fault::InvalidInput)
        );
    }
    assert_eq!(
        BSplineSurface::new(
            [1, 1],
            std::array::from_fn(|_| vec![0.0, 0.0, 1.0, 1.0]),
            [usize::MAX, 2],
            points
        ),
        Err(Fault::InvalidInput)
    );
    let bad = ControlPoint {
        coordinates: [1.0, 2.0, 3.0, 0.0],
    };
    assert_eq!(
        BSplineCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![bad, bad]),
        Err(Fault::InvalidInput)
    );
}

#[test]
fn unsupported_derivative_orders_and_nonfinite_parameters_are_typed_refusals() {
    let curve = BSplineCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0)],
    )
    .unwrap();
    assert_eq!(
        curve.evaluate(0.5, 5, true, false),
        Err(Fault::InvalidInput)
    );
    assert_eq!(curve.evaluate_point(f64::NAN), Err(Fault::InvalidInput));
    assert_eq!(
        curve.control_point(2, catallax_trace::spline::PointForm::Euclidean),
        Err(Fault::InvalidInput)
    );
    assert_eq!(
        curve.evaluate_point(-2.0).unwrap(),
        Vector3::new(0.0, 0.0, 0.0)
    );
    assert_eq!(
        curve.evaluate_point(3.0).unwrap(),
        Vector3::new(1.0, 0.0, 0.0)
    );
}

#[test]
fn a_zero_rational_denominator_is_refused() {
    let curve = BSplineCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![
            ControlPoint {
                coordinates: [0.0, 0.0, 0.0, 1.0],
            },
            ControlPoint {
                coordinates: [-1.0, 0.0, 0.0, -1.0],
            },
        ],
    )
    .unwrap();
    assert_eq!(curve.evaluate_point(0.5), Err(Fault::Degenerate));
    assert_eq!(curve.evaluate(0.5, 1, true, false), Err(Fault::Degenerate));
}

#[test]
fn overflowing_endpoint_projection_is_refused_even_when_clamped() {
    let overflowing = ControlPoint {
        coordinates: [f64::MAX, 0.0, 0.0, 0.5],
    };
    let curve =
        BSplineCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![overflowing, overflowing]).unwrap();
    assert_eq!(curve.evaluate_point(0.5), Err(Fault::Degenerate));
    for parameter in [-1.0, 0.0, 1.0, 2.0] {
        assert_eq!(curve.evaluate_point(parameter), Err(Fault::Degenerate));
        assert_eq!(
            curve.evaluate(parameter, 0, true, false),
            Err(Fault::Degenerate)
        );
    }
}

#[test]
fn a_failed_point_evaluation_cannot_be_used_to_make_a_box() {
    let width = 5e-25;
    let curve = BSplineCurve::new(
        2,
        vec![0.0, 0.0, 0.0, width, width, width],
        vec![
            point(0.0, 0.0, 0.0),
            point(1.0, 0.0, 0.0),
            point(2.0, 0.0, 0.0),
        ],
    )
    .unwrap();
    assert_eq!(curve.evaluate_point(width / 2.0), Err(Fault::Degenerate));
    assert_eq!(
        curve.bounding_box(catallax_trace::extent::Extent1::new(
            width * 0.2,
            width * 0.8
        )),
        Err(Fault::Degenerate)
    );
}
