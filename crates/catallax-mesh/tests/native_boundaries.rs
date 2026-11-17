//! Native previews keep repaired boundaries in physical millimetres.
#[path = "../../catallax-import/tests/common/mod.rs"]
mod step_fixtures;

use catallax_import::{sha256_hex, tolerant};
use catallax_mesh::{PreviewOptions, preview_tolerant_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};

#[test]
fn reconstructed_pole_edges_preserve_the_half_sphere_display() {
    let source =
        include_bytes!("../../catallax-import/tests/fixtures/half-sphere-missing-poles.step");
    let (import, import_work) = with_operation_budget(20_000_000, || {
        tolerant::admission::import_step(
            source,
            &sha256_hex(source),
            Default::default(),
            &mut Decider::default(),
        )
    });
    let import = import.unwrap();
    assert!(!import_work.exhausted);
    let options = PreviewOptions {
        chord_tolerance: Rational::parse("0.02").unwrap(),
        angular_deviation: None,
        max_vertices: 20_000,
    };
    // Preview has no implicit work allowance. This checks geometry, not latency.
    let preview = preview_tolerant_body(&import, 0, &options, &mut Decider::default()).unwrap();
    let mut sphere_area = 0.0;
    let mut plane_area = 0.0;
    let mut highest_y: f64 = 0.0;
    for triangle in &preview.triangles {
        let points = triangle.vertices.map(|index| preview.positions[index]);
        let a: [f64; 3] = std::array::from_fn(|axis| points[1][axis] - points[0][axis]);
        let b: [f64; 3] = std::array::from_fn(|axis| points[2][axis] - points[0][axis]);
        let cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let area = cross[0].hypot(cross[1]).hypot(cross[2]) * 0.5;
        assert!(area.is_finite());
        if triangle.face == "#17" {
            sphere_area += area;
            for point in points {
                assert!((point[0].hypot(point[1]).hypot(point[2]) - 1.0).abs() < 1e-6);
                assert!(point[1] >= -1e-6);
                highest_y = highest_y.max(point[1]);
            }
        } else {
            assert!(points.iter().all(|point| point[1].abs() < 1e-6));
            plane_area += area;
        }
    }
    assert!(highest_y > 0.95);
    assert!(
        (sphere_area - std::f64::consts::TAU).abs() < 0.3,
        "sphere area {sphere_area}"
    );
    assert!(
        (plane_area - std::f64::consts::PI).abs() < 0.15,
        "cut-plane area {plane_area}"
    );
    if let Ok(folder) = std::env::var("CATALLAX_HEALING_EVIDENCE") {
        let mut obj = String::new();
        for point in &preview.positions {
            obj.push_str(&format!("v {} {} {}\n", point[0], point[1], point[2]));
        }
        for triangle in &preview.triangles {
            obj.push_str(&format!(
                "f {} {} {}\n",
                triangle.vertices[0] + 1,
                triangle.vertices[1] + 1,
                triangle.vertices[2] + 1
            ));
        }
        std::fs::write(std::path::Path::new(&folder).join("native-poles.obj"), obj).unwrap();
    }
}

#[test]
fn repaired_line_previews_ignore_vector_scale_and_source_length_units() {
    use step_fixtures::{BoxFixture, Unit, box_step};
    for (unit, size) in [
        (Unit::Millimetre, ["10", "6", "4"]),
        (Unit::Metre, ["0.01", "0.006", "0.004"]),
    ] {
        let mut fixture = BoxFixture::new(size);
        fixture.unit = unit;
        let source = box_step(fixture)
            .lines()
            .map(|line| {
                if line.contains("=VECTOR(") {
                    let (prefix, _) = line.rsplit_once(',').unwrap();
                    format!("{prefix},2.);")
                } else {
                    line.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let (import, work) = with_operation_budget(20_000_000, || {
            tolerant::admission::import_step(
                source.as_bytes(),
                &sha256_hex(source.as_bytes()),
                Default::default(),
                &mut Decider::default(),
            )
        });
        let import = import.unwrap();
        assert!(!work.exhausted);
        let preview = preview_tolerant_body(
            &import,
            0,
            &PreviewOptions {
                chord_tolerance: Rational::parse("0.02").unwrap(),
                angular_deviation: None,
                max_vertices: 2_000,
            },
            &mut Decider::default(),
        )
        .unwrap();
        for (axis, extent) in [10.0, 6.0, 4.0].into_iter().enumerate() {
            let low = preview
                .positions
                .iter()
                .map(|point| point[axis])
                .fold(f64::INFINITY, f64::min);
            let high = preview
                .positions
                .iter()
                .map(|point| point[axis])
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(low.abs() < 1e-9 && (high - extent).abs() < 1e-9);
        }
        let area = preview
            .triangles
            .iter()
            .map(|triangle| {
                let points = triangle.vertices.map(|index| preview.positions[index]);
                let a: [f64; 3] = std::array::from_fn(|axis| points[1][axis] - points[0][axis]);
                let b: [f64; 3] = std::array::from_fn(|axis| points[2][axis] - points[0][axis]);
                let cross = [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ];
                cross[0].hypot(cross[1]).hypot(cross[2]) * 0.5
            })
            .sum::<f64>();
        assert!((area - 248.0).abs() < 1e-6, "physical box area {area}");
    }
}
