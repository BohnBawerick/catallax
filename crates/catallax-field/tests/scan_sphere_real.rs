//! Phase 4.4 acceptance on measured points from a GOM ATOS calibration scan.

use catallax_field::reconstruction::scan::sphere::fit_sphere_scan;
use catallax_field::{Location, Workspace};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{FileRef, Operation, Pose};
use catallax_shapes::RPoint3;

const POINTS: &[u8] = include_bytes!("fixtures/gom-atos-target-11.obj");
const WORK: u64 = 20_000_000;
const MAX_DEVIATION_MM: f64 = 0.06;

fn step() -> Operation {
    Operation::ScanImport {
        file: FileRef {
            name: "gom-atos-target-11.obj".into(),
            sha256: catallax_import::sha256_hex(POINTS),
        },
        pose: Pose::identity(),
        scale: 1.into(),
    }
}

#[test]
fn measured_sphere_rebuilds_in_both_kernels_and_has_a_display_mesh() {
    assert_eq!(
        catallax_import::sha256_hex(POINTS),
        "c81d85fe2f45e70923d9ca6d2b4dee6ef2491ea23e6925ea92799dc4a3e94ce4"
    );
    let mut dec = Decider::default();
    let (fit, work) = with_operation_budget(WORK, || {
        fit_sphere_scan(&step(), POINTS, MAX_DEVIATION_MM, &mut dec)
    });
    let fit = fit.unwrap();
    assert!(!work.exhausted);
    let Operation::Sphere { radius, pose } = &fit.recipe.steps[1].op else {
        panic!("the clean recipe must contain an analytic sphere");
    };
    assert!((radius.value().to_f64_nearest() - 18.928).abs() <= 0.005);
    let [x, y, z] = pose.translation.0.each_ref().map(|n| n.value().clone());
    let point = RPoint3::new(x, y, z);
    assert_eq!(
        fit.field.proven(&point, &mut Workspace::default()).location,
        Location::Inside
    );
    assert!(fit.brep.model.check(fit.brep.result).is_empty());
    assert_eq!(fit.scan_to_part.samples, 13_154);
    assert!(fit.scan_to_part.max_mm <= MAX_DEVIATION_MM);
    assert!(fit.scan_to_part.rms_mm < 0.018);
    assert_eq!(fit.part_to_scan.samples, 512);
    // The rear of the mounted target has no samples.
    assert!((15.0..18.0).contains(&fit.part_to_scan.max_mm));
    println!(
        "fit work {}. scan to part {:?}. part to scan {:?}",
        work.used, fit.scan_to_part, fit.part_to_scan
    );

    let chord = Rational::from_ratio(1, 1000).unwrap();
    let (display, work) = with_operation_budget(WORK, || {
        mesh_body(
            &fit.brep.model,
            fit.brep.result,
            &MeshOptions {
                chord_tolerance: chord.clone(),
                ..Default::default()
            },
            &mut Decider::default(),
        )
    });
    let mesh = display.expect("certified display mesh").value;
    assert!(!work.exhausted);
    assert!(!mesh.triangles.is_empty());
    assert_eq!(mesh.faces.len(), 1);
    assert!(mesh.faces[0].max_chord_error <= chord);
    println!(
        "display work {}. {} triangles. chord {} mm",
        work.used,
        mesh.triangles.len(),
        mesh.faces[0].max_chord_error.to_f64_nearest()
    );
    if let Some(path) = std::env::var_os("CATALLAX_SCAN_DISPLAY_OBJ") {
        mesh.write_obj(&mut std::fs::File::create(path).unwrap())
            .unwrap();
    }
}

#[test]
fn measured_sphere_refuses_when_the_points_exceed_the_requested_deviation() {
    let (fit, work) = with_operation_budget(WORK, || {
        fit_sphere_scan(&step(), POINTS, 0.01, &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert_eq!(fit.err().unwrap().code, "reconstruction.sphere-deviation");
}

#[test]
fn translated_sphere_refuses_a_point_hidden_by_binary64_rounding() {
    use std::collections::BTreeSet;
    use std::fmt::Write;

    let mut points = BTreeSet::new();
    for [a, b, c] in [[0, 0, 9], [1, 4, 8], [3, 6, 6], [4, 4, 7]] {
        for [x, y, z] in [
            [a, b, c],
            [a, c, b],
            [b, a, c],
            [b, c, a],
            [c, a, b],
            [c, b, a],
        ] {
            for sx in [-1, 1] {
                for sy in [-1, 1] {
                    for sz in [-1, 1] {
                        points.insert([sx * x, sy * y, sz * z]);
                    }
                }
            }
        }
    }

    let mut cloud = String::new();
    for [x, y, z] in points {
        let suffix = if [x, y, z] == [9, 0, 0] { ".003" } else { "" };
        writeln!(cloud, "v {}{suffix} {y} {z}", (1_i64 << 45) + x).unwrap();
    }
    let bytes = cloud.as_bytes();
    let scan = Operation::ScanImport {
        file: FileRef {
            name: "translated-sphere.obj".into(),
            sha256: catallax_import::sha256_hex(bytes),
        },
        pose: Pose::identity(),
        scale: 1.into(),
    };
    let (fit, work) = with_operation_budget(WORK, || {
        fit_sphere_scan(&scan, bytes, 0.001, &mut Decider::default())
    });
    assert!(!work.exhausted);
    let error = fit.err().expect("the 0.003 mm outlier must be refused");
    assert_eq!(error.code, "reconstruction.sphere-deviation", "{error:?}");
}
