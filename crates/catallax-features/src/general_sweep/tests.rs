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
    FrameLaw {
        path,
        guide: None,
        reference: [1.0, 0.0, 0.0],
        initial,
        normal_offset: 0.0,
        twist,
        end_scale: 1.0,
        orientation: catallax_recipe::SweepOrientation::RotationMinimising,
    }
    .at(t, &mut Decider::default())
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
    let frame = FrameLaw {
        path: &path,
        guide: Some(&guide),
        reference: [1.0, 0.0, 0.0],
        initial: None,
        normal_offset: 0.0,
        twist: std::f64::consts::FRAC_PI_2,
        end_scale: 1.0,
        orientation: catallax_recipe::SweepOrientation::RotationMinimising,
    }
    .at(0.5, &mut Decider::default())
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
    let refusal = FrameLaw {
        path: &path,
        guide: None,
        reference: [1.0, 0.0, 0.0],
        initial: None,
        normal_offset: 0.0,
        twist: 0.0,
        end_scale: 1.0,
        orientation: catallax_recipe::SweepOrientation::RotationMinimising,
    }
    .at(0.5, &mut Decider::default())
    .err()
    .expect("a general spline needs an integrated frame, not a projected reference");
    assert_eq!(refusal.code, "sweep.spline-transport");
}

#[test]
fn standalone_helix_heading_places_an_asymmetric_section() {
    let path = helix(1.0, 2.0);
    let start = FrameLaw {
        path: &path,
        guide: None,
        reference: [0.0, 0.0, 1.0],
        initial: None,
        normal_offset: 0.0,
        twist: 0.0,
        end_scale: 1.0,
        orientation: catallax_recipe::SweepOrientation::RotationMinimising,
    }
    .at(0.0, &mut Decider::default())
    .unwrap();
    let speed = libm::sqrt((4.0 * std::f64::consts::TAU).powi(2) + 4.0);
    let horizontal = 4.0 * std::f64::consts::TAU / speed;
    let vertical = 2.0 / speed;
    for [x, y] in [[2.0, 0.0], [0.0, 1.0], [-1.0, 0.0]] {
        close(
            start.apply([x, y, 0.0]),
            [4.0 + y, -vertical * x, horizontal * x],
            1e-12,
        );
    }
}

fn line_definition() -> SweepRecipe {
    serde_json::from_value(serde_json::json!({
        "name": "checked profile",
        "loops": [[
            {"kind": "line", "id": "a", "to": ["4", "0"]},
            {"kind": "line", "id": "b", "to": ["4", "2"]},
            {"kind": "line", "id": "c", "to": ["0", "2"]},
            {"kind": "line", "id": "d", "to": ["0", "0"]}
        ]],
        "path": {"kind": "line", "start": ["0", "0", "0"], "end": ["0", "0", "4"]},
        "reference": ["1", "0", "0"],
        "twist_degrees": "0",
        "end_scale": "1",
        "tolerance_mm": "0.001",
        "max_sections": 1
    }))
    .unwrap()
}

#[test]
fn translated_profile_rounding_cannot_hide_a_sample_displacement() {
    for axis in [0, 2] {
        let mut recipe = line_definition();
        let mut plane = RecipePose::identity();
        plane.translation.0[axis] = Rational::from_ratio(100_000_000_000_000_001, 10)
            .unwrap()
            .into();
        recipe.plane = Some(plane);
        let SweepPath::Line { start, end } = &mut recipe.path else {
            unreachable!();
        };
        start.0[axis] = q(10_000_000_000_000_000).into();
        end.0[axis] = q(if axis == 2 {
            10_000_000_000_000_004
        } else {
            10_000_000_000_000_000
        })
        .into();
        let error = build(&recipe, &mut Decider::default())
            .expect_err("the emitted section lost 0.1 mm at a 0.001 mm tolerance");
        assert_eq!(error.code, "sweep.numerical-definition");
    }
    let offset = Rational::from_ratio(10_000_000_001, 10).unwrap();
    let mut recipe = line_definition();
    let mut plane = RecipePose::identity();
    plane.translation.0[0] = offset.clone().into();
    recipe.plane = Some(plane);
    let SweepPath::Line { start, end } = &mut recipe.path else {
        unreachable!();
    };
    start.0[0] = q(1_000_000_000).into();
    end.0[0] = q(1_000_000_000).into();
    let swept = build(&recipe, &mut Decider::default()).unwrap();
    let conversion = (rational(offset.to_f64_nearest()).unwrap() - offset)
        .abs()
        .to_f64_nearest();
    assert!(conversion > 0.0);
    assert!(swept.checks().max_profile_sample_error_mm >= conversion);
    assert!(swept.checks().max_profile_sample_error_mm <= 0.001 / 4.0);
}

#[test]
fn standalone_build_and_edit_preserve_the_callers_work_allowance() {
    let recipe = line_definition();
    let (result, work) =
        catallax_numbers::with_operation_budget(1, || build(&recipe, &mut Decider::default()));
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
    assert!(work.used <= 1);
    let swept = build(&recipe, &mut Decider::default()).unwrap();
    let (result, work) =
        catallax_numbers::with_operation_budget(1, || swept.edit(recipe, &mut Decider::default()));
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
    assert!(work.used <= 1);
}

#[test]
fn authored_spine_endpoint_rounding_is_checked_before_admission() {
    let mut lost = line_definition();
    let SweepPath::Line { end, .. } = &mut lost.path else {
        unreachable!()
    };
    end.0[2] = Rational::from_ratio(100_000_000_000_000_001, 10)
        .unwrap()
        .into();
    let error = build(&lost, &mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "sweep.motion-sample-tolerance");

    let mut lost_guide = line_definition();
    let SweepPath::Line { start, end } = &mut lost_guide.path else {
        unreachable!()
    };
    start.0[0] = q(10_000_000_000_000_000).into();
    end.0[0] = q(10_000_000_000_000_000).into();
    lost_guide.guide = Some(
        serde_json::from_value(serde_json::json!({
            "kind": "line",
            "start": ["10000000000000000.1", "0.0000001", "0"],
            "end": ["10000000000000000.1", "0.0000001", "4"]
        }))
        .unwrap(),
    );
    assert_eq!(
        build(&lost_guide, &mut Decider::default())
            .unwrap_err()
            .code,
        "sweep.motion-sample-tolerance",
    );

    let mut retained = line_definition();
    let endpoint = Rational::from_ratio(41, 10).unwrap();
    let SweepPath::Line { end, .. } = &mut retained.path else {
        unreachable!()
    };
    end.0[2] = endpoint.clone().into();
    let swept = build(&retained, &mut Decider::default()).unwrap();
    let conversion = (rational(endpoint.to_f64_nearest()).unwrap() - endpoint)
        .abs()
        .to_f64_nearest();
    assert!(conversion > 0.0);
    let error = swept.checks().max_motion_endpoint_error_mm.unwrap();
    assert!(error >= conversion);
    assert!(error <= 0.001 / 4.0);
    assert_eq!(
        swept.edit(lost, &mut Decider::default()).unwrap_err().code,
        "sweep.motion-sample-tolerance",
    );
}

#[test]
fn authored_line_difference_preserves_the_section_heading() {
    let mut recipe = line_definition();
    recipe.path = serde_json::from_value(serde_json::json!({
        "kind": "line",
        "start": ["10000000000000000.1", "0", "0"],
        "end": ["10000000000000000.2", "0", "4"]
    }))
    .unwrap();
    let mut dec = Decider::default();
    let path = Path::new(&recipe.path, &mut dec).unwrap();
    let frame = FrameLaw {
        path: &path,
        guide: None,
        reference: [1.0, 0.0, 0.0],
        initial: None,
        normal_offset: 0.0,
        twist: 0.0,
        end_scale: 1.0,
        orientation: recipe.orientation,
    }
    .at(0.5, &mut dec)
    .unwrap();
    let normal = cross(frame.x, frame.y);
    assert!(normal[0] > 0.02 && normal[0] < 0.03);
    assert!(normal[2] > 0.99);
}
