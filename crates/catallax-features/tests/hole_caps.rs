//! Native cap/wall curves from the compiled E07 Hole frames.

use catallax_features::hole::{CylinderSpan, HoleTool, cap_curve};
use catallax_numbers::{Ball, Decider, NumericalCheck, Rational, Real, with_operation_budget};
use catallax_recipe::{HoleDepth, HolePlacement, HoleTip};
use catallax_shapes::{BSplineSurface, Plane, Point3, Vector3};
use catallax_topology::numerical::NumericalOptions;
use serde_json::json;

#[path = "hole_caps/shell.rs"]
mod shell;
#[path = "hole_caps/skew.rs"]
mod skew;

fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn vector(moved: bool, values: [i64; 3]) -> Vector3 {
    let [x, y, z] = values;
    let values = if moved { [x, -z, y] } else { values };
    let [x, y, z] = values.map(Real::from_i64);
    Vector3::new(x, y, z)
}

fn point(moved: bool, [x, y, z]: [&str; 3]) -> Point3 {
    let [x, y, z] = [q(x), q(y), q(z)];
    let values = if moved {
        [x + q("7"), -z + q("11"), y - q("13")]
    } else {
        [x, y, z]
    };
    let [x, y, z] = values.map(Real::from);
    Point3::new(x, y, z)
}

fn span(moved: bool, skew: bool, dec: &mut Decider) -> CylinderSpan {
    span_with_diameter(moved, skew, &q("2"), dec)
}

fn span_with_diameter(
    moved: bool,
    skew: bool,
    diameter: &Rational,
    dec: &mut Decider,
) -> CylinderSpan {
    let point = point(moved, ["21/2", "12", "0"]);
    let direction = vector(moved, if skew { [1, -1, 1] } else { [0, -1, 0] });
    let at: HolePlacement = serde_json::from_value(json!({
        "id":"native-cap", "point":[point.x.as_rational().unwrap().to_string(),
            point.y.as_rational().unwrap().to_string(), point.z.as_rational().unwrap().to_string()],
        "direction":[direction.x.as_rational().unwrap().to_string(),
            direction.y.as_rational().unwrap().to_string(), direction.z.as_rational().unwrap().to_string()]
    }))
    .unwrap();
    HoleTool::new(
        &at,
        diameter,
        &HoleDepth::Through,
        &HoleTip::Flat,
        None,
        &q("100"),
        dec,
    )
    .unwrap()
    .cylinder_span()
    .unwrap()
}

fn quarter(span: &CylinderSpan, index: i64, dec: &mut Decider) -> BSplineSurface {
    let axis = span
        .cylinder
        .axis
        .scaled(&(&Real::one() / &span.cylinder.axis.norm_squared().sqrt()));
    let x = span.cylinder.radial_reference();
    let y = axis.cross(&x);
    let minus = Real::from_i64(-1);
    let (x, y) = match index {
        0 => (x, y),
        1 => (y, x.scaled(&minus)),
        2 => (x.scaled(&minus), y.scaled(&minus)),
        3 => (y.scaled(&minus), x),
        _ => unreachable!(),
    };
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    for (radial, weight) in [
        (x.clone(), Real::one()),
        (x.plus(&y), Real::from(q("7/10"))),
        (y, Real::from(q("49/50"))),
    ] {
        for height in &span.interval_mm {
            let p = span
                .cylinder
                .origin
                .plus(&axis.scaled(&Real::from(height.clone())))
                .plus(&radial);
            let [x, y, z] = [&p.x, &p.y, &p.z]
                .map(|value| Real::from(dec.enclose_budgeted_at(value, 128).unwrap().mid()));
            controls.push(Point3::new(x, y, z));
            weights.push(weight.clone());
        }
    }
    BSplineSurface::new(
        [2, 1],
        [
            [index, index, index, index + 1, index + 1, index + 1]
                .map(Rational::from_i64)
                .to_vec(),
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        ],
        [3, 2],
        controls,
        weights,
        dec,
    )
    .unwrap()
}

fn cap(moved: bool, top: bool) -> Plane {
    Plane {
        origin: point(moved, ["0", if top { "12" } else { "0" }, "0"]),
        normal: vector(moved, [0, 1, 0]),
        ref_dir: vector(moved, [1, 0, 0]),
    }
}

fn floats(values: [Ball; 3]) -> [f64; 3] {
    values.map(|value| value.mid().to_f64_nearest())
}

#[test]
fn native_cap_curves_match_independent_surface_section_evaluations() {
    let options = NumericalOptions {
        tolerance_mm: q("1e-7"),
    };
    for moved in [false, true] {
        for skew in [false, true] {
            let mut dec = Decider::default();
            let span = span(moved, skew, &mut dec);
            for top in [false, true] {
                let plane = cap(moved, top);
                let n = [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(Real::to_f64);
                let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z].map(Real::to_f64);
                for index in 0..4 {
                    let source = quarter(&span, index, &mut dec);
                    let (output, work) = with_operation_budget(20_000_000, || {
                        cap_curve("native-cap", &span, &source, &plane, &options, &mut dec)
                    });
                    let output = output.unwrap().unwrap();
                    assert!(!work.exhausted);
                    println!(
                        "cap moved={moved} skew={skew} top={top} quarter={index} work={}",
                        work.used
                    );
                    assert_eq!(output.curve().domain(), source.domain()[0]);
                    assert_eq!(output.plane_uv().domain(), source.domain()[0]);
                    assert_eq!(output.cylinder_uv().domain(), source.domain()[0]);
                    assert_eq!(output.cylinder_uv().degree(), 3);
                    assert_eq!(output.evidence().checks().len(), 6);
                    for check in output.evidence().checks() {
                        if let NumericalCheck::Passed { name, .. } = check {
                            assert_eq!(name, "hole.cap.finite-v-positive-weight-hulls");
                            continue;
                        }
                        let NumericalCheck::Samples { name, samples, .. } = check else {
                            panic!("expected actual composition samples")
                        };
                        assert_eq!(
                            samples.len(),
                            if name.ends_with("source-axial-coordinates") {
                                6
                            } else {
                                17
                            }
                        );
                        assert!(samples.iter().all(|s| s.distance_upper_mm <= q("1e-7")));
                    }
                    for share in ["0", "1/7", "1/2", "4/5", "1"] {
                        let t = Rational::from_i64(index) + q(share);
                        let at = Ball::from_rational(&t, 128);
                        let low = floats(
                            source
                                .evaluate(&at, &Ball::from_rational(&q("0"), 128), &dec)
                                .unwrap(),
                        );
                        let high = floats(
                            source
                                .evaluate(&at, &Ball::from_rational(&q("1"), 128), &dec)
                                .unwrap(),
                        );
                        let numerator: f64 = (0..3).map(|j| n[j] * (origin[j] - low[j])).sum();
                        let denominator: f64 = (0..3).map(|j| n[j] * (high[j] - low[j])).sum();
                        let fraction = numerator / denominator;
                        let expected: [f64; 3] =
                            std::array::from_fn(|j| low[j] + fraction * (high[j] - low[j]));
                        let actual = floats(output.curve().evaluate(&at, &dec).unwrap());
                        for (a, b) in actual.into_iter().zip(expected) {
                            assert!(
                                (a - b).abs() < 1e-8,
                                "{moved} {skew} {top} {index}: {a} {b}"
                            );
                        }
                        let uv = floats(output.cylinder_uv().evaluate(&at, &dec).unwrap());
                        assert!((uv[0] - t.to_f64_nearest()).abs() < 1e-12);
                        assert!((uv[1] - fraction).abs() < 1e-12);
                    }
                }
            }
        }
    }
}

#[test]
fn wrong_carrier_or_interval_refuses_instead_of_publishing_a_cap_curve() {
    let mut dec = Decider::default();
    let span = span(false, true, &mut dec);
    let source = quarter(&span, 0, &mut dec);
    let options = NumericalOptions {
        tolerance_mm: q("1e-7"),
    };
    let plane = cap(false, true);
    let mut wrong = span.clone();
    wrong.cylinder.radius_sq = Real::from_i64(4);
    assert_eq!(
        cap_curve("wrong-radius", &wrong, &source, &plane, &options, &mut dec)
            .unwrap_err()
            .code,
        "numbers.numerical-tolerance-exceeded"
    );
    wrong = span.clone();
    wrong.interval_mm = [q("-1"), q("1")];
    assert!(
        cap_curve(
            "wrong-interval",
            &wrong,
            &source,
            &plane,
            &options,
            &mut dec
        )
        .is_err()
    );
    wrong.interval_mm = [
        &span.interval_mm[0] + &q("1"),
        &span.interval_mm[1] + &q("1"),
    ];
    assert!(
        cap_curve(
            "shifted-interval",
            &wrong,
            &source,
            &plane,
            &options,
            &mut dec
        )
        .is_err()
    );
    let (result, work) = with_operation_budget(0, || {
        cap_curve("cancelled", &span, &source, &plane, &options, &mut dec)
    });
    assert!(work.exhausted && result.is_err());
}

#[test]
fn a_plane_parallel_to_the_axis_has_no_ellipse_template() {
    let mut dec = Decider::default();
    let span = span(false, false, &mut dec);
    let source = quarter(&span, 0, &mut dec);
    let plane = Plane {
        origin: point(false, ["0", "0", "0"]),
        normal: vector(false, [1, 0, 0]),
        ref_dir: vector(false, [0, 1, 0]),
    };
    assert!(
        cap_curve(
            "parallel",
            &span,
            &source,
            &plane,
            &NumericalOptions {
                tolerance_mm: q("1e-7")
            },
            &mut dec
        )
        .unwrap()
        .is_none()
    );
}

fn finite_v_fixture(extent: &str, dec: &mut Decider) -> (CylinderSpan, BSplineSurface, Plane) {
    let placement = serde_json::from_value(json!({
        "id":"finite-v", "point":["0","0","0"], "direction":["0","0","1"]
    }))
    .unwrap();
    let span = HoleTool::new(
        &placement,
        &q("2"),
        &HoleDepth::Through,
        &HoleTip::Flat,
        None,
        &(&q(extent) - &q("1")),
        dec,
    )
    .unwrap()
    .cylinder_span()
    .unwrap();
    assert_eq!(span.interval_mm, [-q(extent), q(extent)]);
    let axis = span.cylinder.axis.clone();
    let x = span.cylinder.radial_reference();
    let y = axis.cross(&x);
    let plane = Plane::new(
        span.cylinder.origin.clone(),
        axis.minus(&x.scaled(&Real::from_i64(4)))
            .minus(&y.scaled(&Real::from_i64(3))),
        x.scaled(&Real::from_i64(3))
            .minus(&y.scaled(&Real::from_i64(4))),
        dec,
    )
    .unwrap();
    let source = quarter(&span, 0, dec);
    (span, source, plane)
}

fn projected_v(source: &BSplineSurface, plane: &Plane, t: Rational, dec: &Decider) -> f64 {
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(Real::to_f64);
    let at = Ball::from_rational(&t, 128);
    let a = floats(
        source
            .evaluate(&at, &Ball::from_rational(&q("0"), 128), dec)
            .unwrap(),
    );
    let b = floats(
        source
            .evaluate(&at, &Ball::from_rational(&q("1"), 128), dec)
            .unwrap(),
    );
    let numerator: f64 = (0..3).map(|j| -normal[j] * a[j]).sum();
    let denominator: f64 = (0..3).map(|j| normal[j] * (b[j] - a[j])).sum();
    numerator / denominator
}

#[test]
fn finite_v_excursion_between_seventeen_stations_is_rejected() {
    let mut dec = Decider::default();
    let (span, source, plane) = finite_v_fixture("4.999", &mut dec);
    assert!((0..=16).all(|i| (0.0..=1.0).contains(&projected_v(
        &source,
        &plane,
        Rational::from_ratio(i, 16).unwrap(),
        &dec
    ))));
    assert!(projected_v(&source, &plane, q("5/12"), &dec) > 1.0);
    let result = cap_curve(
        "between-stations",
        &span,
        &source,
        &plane,
        &NumericalOptions {
            tolerance_mm: q("1e-7"),
        },
        &mut dec,
    );
    assert!(
        result.is_err(),
        "a cap leaving finite V between stations must refuse"
    );
    assert_eq!(result.unwrap_err().code, "hole.cap-finite-v");
}

#[test]
fn finite_v_subdivision_accepts_a_safe_curve_with_an_inconclusive_initial_hull() {
    let mut dec = Decider::default();
    let (span, source, plane) = finite_v_fixture("5.001", &mut dec);
    assert!(projected_v(&source, &plane, q("5/12"), &dec) < 1.0);
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(Real::to_f64);
    let a = &source.controls()[2];
    let b = &source.controls()[3];
    let a = [&a.x, &a.y, &a.z].map(Real::to_f64);
    let b = [&b.x, &b.y, &b.z].map(Real::to_f64);
    let numerator: f64 = (0..3).map(|j| -normal[j] * a[j]).sum();
    let denominator: f64 = (0..3).map(|j| normal[j] * (b[j] - a[j])).sum();
    assert!(numerator / denominator > 1.0);
    let output = cap_curve(
        "subdivided-hulls",
        &span,
        &source,
        &plane,
        &NumericalOptions {
            tolerance_mm: q("1e-7"),
        },
        &mut dec,
    )
    .unwrap()
    .unwrap();
    assert!(output.evidence().checks().iter().any(|check|
        matches!(check, NumericalCheck::Passed { name, .. } if name == "hole.cap.finite-v-positive-weight-hulls")));
}
