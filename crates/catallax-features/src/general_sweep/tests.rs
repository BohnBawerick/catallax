use super::*;

fn helix(turns: f64, rise: f64) -> Path {
    Path::Helix {
        radius: 4.0,
        pitch: rise / turns,
        turns,
        origin: [0.0; 3],
        axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    }
}
fn frame(path: &Path, t: f64, twist: f64, initial: Option<[V; 2]>) -> Frame {
    Frame::at(
        path,
        None,
        [1.0, 0.0, 0.0],
        initial,
        0.0,
        twist,
        1.0,
        t,
        catallax_recipe::SweepOrientation::RotationMinimising,
        &mut Decider::default(),
    )
    .unwrap()
}
fn close(a: V, b: V, tolerance: f64) {
    let delta = sub(a, b);
    assert!(dot(delta, delta).sqrt() <= tolerance, "{a:?} != {b:?}");
}

#[test]
fn helix_transport_has_only_the_authored_twist_about_the_tangent() {
    for turns in [-1.0, 1.0] {
        for rise in [-2.0, 2.0] {
            for twist in [0.0, std::f64::consts::PI] {
                let path = helix(turns, rise);
                for t in [0.2, 0.5, 0.8] {
                    let h = 1e-5;
                    let derivative = mul(
                        sub(
                            frame(&path, t + h, twist, None).x,
                            frame(&path, t - h, twist, None).x,
                        ),
                        0.5 / h,
                    );
                    let spin = dot(derivative, frame(&path, t, twist, None).y);
                    assert!((spin - twist).abs() < 1e-6, "spin {spin}, authored {twist}");
                }
            }
        }
    }
}

#[test]
fn oblique_authored_section_is_not_realigned_to_the_spine() {
    let path = helix(1.0, 2.0);
    let axes = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let start = frame(&path, 0.0, 0.0, Some(axes));
    close(start.x, axes[0], 1e-12);
    close(start.y, axes[1], 1e-12);
    let mut dec = Decider::default();
    let initial_angle = dot(
        cross(start.x, start.y),
        unit(path.at(0.0, &mut dec).unwrap().1).unwrap(),
    );
    let end = frame(&path, 0.7, 0.0, Some(axes));
    let final_angle = dot(
        cross(end.x, end.y),
        unit(path.at(0.7, &mut dec).unwrap().1).unwrap(),
    );
    assert!((final_angle - initial_angle).abs() < 1e-12);
    assert!(initial_angle < 0.0);
}

#[test]
fn guide_rotation_and_authored_twist_both_move_the_section() {
    let path = Path::Line([0.0; 3], [0.0, 0.0, 4.0]);
    let guide = Path::Line([1.0, 0.0, 0.0], [-1.0, 1.0, 4.0]);
    let frame = Frame::at(
        &path,
        Some(&guide),
        [1.0, 0.0, 0.0],
        None,
        0.0,
        std::f64::consts::FRAC_PI_2,
        1.0,
        0.5,
        catallax_recipe::SweepOrientation::RotationMinimising,
        &mut Decider::default(),
    )
    .unwrap();
    close(frame.apply([2.0, 0.0, 0.0]), [0.0, 2.0, 2.0], 1e-12);
}

#[test]
fn a_projected_reference_is_not_mislabelled_as_spline_parallel_transport() {
    let mut dec = Decider::default();
    let path = Path::Spline(
        BSplineCurve::new(
            3,
            [0, 0, 0, 0, 1, 1, 1, 1].map(q).to_vec(),
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 2.0],
                [1.0, 1.0, 3.0],
            ]
            .map(|p| point(p).unwrap())
            .to_vec(),
            vec![Real::from(1); 4],
            &mut dec,
        )
        .unwrap(),
    );
    let refusal = Frame::at(
        &path,
        None,
        [1.0, 0.0, 0.0],
        None,
        0.0,
        0.0,
        1.0,
        0.5,
        catallax_recipe::SweepOrientation::RotationMinimising,
        &mut dec,
    )
    .err()
    .expect("a general spline needs an integrated frame, not a projected reference");
    assert_eq!(refusal.code, "sweep.spline-transport");
}

#[test]
fn profile_chart_subdivision_preserves_rational_geometry_and_parameter() {
    let mut dec = Decider::default();
    let chart = crate::loft::charts::Chart {
        bounds: [q(0), q(1)],
        curve: BSplineCurve::new(
            2,
            [0, 0, 0, 1, 1, 1].map(q).to_vec(),
            [[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]]
                .map(|p| point(p).unwrap())
                .to_vec(),
            [q(1), q(1).mul_pow2(-1), q(1)].map(Real::from).to_vec(),
            &mut dec,
        )
        .unwrap(),
    };
    let pieces = profile_charts(chart, &mut dec).unwrap();
    let ratio = |a, b| Rational::from_ratio(a, b).unwrap();
    for (piece, parameter, expected) in [
        (0, ratio(1, 4), [ratio(5, 13), ratio(6, 13)]),
        (1, ratio(3, 4), [ratio(21, 13), ratio(6, 13)]),
        (0, ratio(1, 2), [q(1), ratio(2, 3)]),
        (1, ratio(1, 2), [q(1), ratio(2, 3)]),
    ] {
        let point = pieces[piece]
            .as_ref()
            .unwrap()
            .curve
            .evaluate_box([&parameter, &parameter], 128, &dec)
            .unwrap();
        for axis in 0..2 {
            assert!(point[axis].lo() <= expected[axis]);
            assert!(point[axis].hi() >= expected[axis]);
            assert!(point[axis].hi() - point[axis].lo() < ratio(1, 1_000_000_000_000_000_000));
        }
    }
}
