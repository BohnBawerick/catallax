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
