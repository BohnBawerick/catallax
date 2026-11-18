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

fn volume(preview: &catallax_mesh::PreviewMesh<String, String>) -> f64 {
    preview
        .triangles
        .iter()
        .map(|triangle| {
            let [a, b, c] = triangle.vertices.map(|index| preview.positions[index]);
            (a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

// A vertex of a tolerant body is where its edges were made to meet, here
// 0.1 mm off the plane of one of its faces. The triangles at it miss that
// plane by ten times the chord target however often they are split, and each
// carrier sample used to be enclosed at 96 bits: the preview spent the whole
// request allowance on one face and refused.
#[test]
fn a_vertex_admitted_off_a_face_previews_within_the_request_allowance() {
    use step_fixtures::{BoxFixture, box_step};
    let source = box_step(BoxFixture::new(["10", "8", "6"]));
    let vertex = source
        .lines()
        .find_map(|line| line.split_once("=VERTEX_POINT('',")?.1.strip_suffix(");"))
        .unwrap();
    let mut moved = None;
    let source = source
        .lines()
        .map(|line| {
            let Some(rest) = line.strip_prefix(&format!("{vertex}=CARTESIAN_POINT(")) else {
                return line.to_owned();
            };
            let (name, coordinates) = rest.split_once(",(").unwrap();
            let mut point: Vec<f64> = coordinates
                .trim_end_matches(");")
                .trim_end_matches(')')
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            // Along the box's X edge, so the vertex leaves only the X face.
            point[0] += if point[0] == 0.0 { 0.1 } else { -0.1 };
            moved = Some([point[0], point[1], point[2]]);
            format!(
                "{vertex}=CARTESIAN_POINT({name},({:?},{:?},{:?}));",
                point[0], point[1], point[2]
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    let moved = moved.unwrap();
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
    for angular_deviation in [None, Some(Rational::parse("0.2617993875").unwrap())] {
        let options = PreviewOptions {
            chord_tolerance: Rational::parse("0.01").unwrap(),
            angular_deviation,
            max_vertices: 2_000,
        };
        let (preview, work) = with_operation_budget(100_000_000, || {
            preview_tolerant_body(&import, 0, &options, &mut Decider::default())
        });
        let preview = preview.unwrap();
        assert!(work.used < 1_000_000, "preview work {}", work.used);
        assert_eq!(preview.faces.len(), 6);
        // The vertex is drawn where it was admitted, not pulled onto a face.
        assert!(preview.positions.contains(&moved), "{moved:?}");
        // The box, less the sliver the moved corner cuts off.
        let volume = volume(&preview);
        assert!((479.0..480.0).contains(&volume), "volume {volume}");
    }
}

// Curved native carriers: every sample of the sphere is a sine and cosine
// evaluation, and the face takes thousands of them.
#[test]
fn a_native_sphere_previews_within_a_fraction_of_the_request_allowance() {
    let source =
        include_bytes!("../../catallax-import/tests/fixtures/half-sphere-missing-poles.step");
    let (import, work) = with_operation_budget(20_000_000, || {
        tolerant::admission::import_step(
            source,
            &sha256_hex(source),
            Default::default(),
            &mut Decider::default(),
        )
    });
    let import = import.unwrap();
    assert!(!work.exhausted);
    let options = PreviewOptions {
        chord_tolerance: Rational::parse("0.02").unwrap(),
        angular_deviation: None,
        max_vertices: 20_000,
    };
    let (preview, work) = with_operation_budget(100_000_000, || {
        preview_tolerant_body(&import, 0, &options, &mut Decider::default())
    });
    let preview = preview.unwrap();
    assert!(work.used < 5_000_000, "preview work {}", work.used);
    // A sphere's point is its unit normal times the radius: every drawn
    // vertex of the curved face lies on the unit sphere it samples.
    let sphere = preview
        .faces
        .iter()
        .map(|face| &preview.triangles[face.triangles.clone()])
        .max_by_key(|triangles| triangles.len())
        .unwrap();
    for triangle in sphere {
        for index in triangle.vertices {
            let radius = preview.positions[index]
                .iter()
                .fold(0.0f64, |radius, value| radius.hypot(*value));
            assert!((radius - 1.0).abs() < 1e-6, "radius {radius}");
        }
    }
}

// SolidDraft's STEP export of a pocket rim rounded at R2 by the rolling ball:
// 19 faces, among them torus strips and fitted spline cuts. Opened by default
// it is a native tolerant body, and its preview used to stop at the 100M
// request allowance. The source preview of the same file takes 4.6M.
#[test]
#[ignore = "opening the rim costs 1.09G work, 47 s in release; run on demand"]
fn a_reopened_rolling_ball_rim_previews_within_the_request_allowance() {
    use catallax_import::tolerant::opened::{BodyView, open_step};
    let source = include_bytes!("fixtures/reopened-rolling-ball-rim-r2.step");
    // The default opening installs no allowance for native admission.
    let opened = open_step(
        source,
        &sha256_hex(source),
        &catallax_recipe::StepId::new("import").unwrap(),
        Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let Some(BodyView::Tolerant { import, body }) = opened.body(0) else {
        panic!("the rim opens as a native tolerant body");
    };
    // SolidDraft's display request: the default chord and 15 degrees.
    let options = PreviewOptions::from(&catallax_mesh::MeshOptions {
        angular_deviation: Some(Rational::parse("0.2617993875").unwrap()),
        ..Default::default()
    });
    let (preview, work) = with_operation_budget(100_000_000, || {
        preview_tolerant_body(import, body, &options, &mut Decider::default())
    });
    let preview = preview.unwrap();
    assert!(work.used < 20_000_000, "preview work {}", work.used);
    assert_eq!(preview.faces.len(), 19);
    // The rim's volume is 317,998.676 mm3; the triangles keep it to a chord.
    let volume = volume(&preview);
    assert!((volume - 317_998.676).abs() < 0.5, "volume {volume}");
}
