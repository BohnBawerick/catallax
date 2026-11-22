//! Application measurement and numerical STEP of a multi-turn round-wire sweep.

use catallax_app::{NumericalOptions, numerical_sweeps};
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_step::numerical as step;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::ffi::OsStr;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const APP_WORK: u64 = 100_000_000;
const TARGET_WORK: u64 = 20_000_000;
static SPRING_DUMP_SEQUENCE: AtomicUsize = AtomicUsize::new(0);

fn spring(turns: u32, pitch_mm: u32) -> Recipe {
    serde_json::from_value(serde_json::json!({
        "format": "catallax-recipe/5", "name": "spring", "steps": [
            {"id": "profile", "op": "sketch", "plane": {
                "rotation": ["1", "1", "0", "0"], "translation": ["0", "0", "0"]
            }, "loops": [[{"kind": "circle", "id": "wire", "centre": ["17", "2.5"], "radius": "1"}]]},
            {"id": "path", "op": "path", "start": ["17", "0", "2.5"], "segments": [
                {"kind": "helix", "id": "coil", "centre": ["0", "0", "2.5"],
                 "axis": ["0", "0", "1"], "angle": (360 * turns).to_string(),
                 "rise": (pitch_mm * turns).to_string()}
            ]},
            {"id": "spring", "op": "sweep", "profile": "profile", "path": "path",
             "orientation": "rotation_minimising"}
        ], "result": "spring"
    }))
    .unwrap()
}

fn near(enclosure: &Enclosure, expected: f64, difference: f64) {
    assert!(
        enclosure.lo.to_f64_nearest() >= expected - difference
            && enclosure.hi.to_f64_nearest() <= expected + difference,
        "{enclosure:?} should lie within {difference} of {expected}"
    );
}

fn write_step_dump(requested: &OsStr, turns: u32, pitch_mm: u32, bytes: &[u8]) -> PathBuf {
    let requested = Path::new(requested);
    let stem = requested
        .file_stem()
        .unwrap_or_else(|| requested.file_name().unwrap_or_default())
        .to_string_lossy();
    let extension = requested
        .extension()
        .and_then(OsStr::to_str)
        .unwrap_or("step");
    let (output, mut file) = loop {
        let sequence = SPRING_DUMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let output = requested.with_file_name(format!(
            "{stem}-{turns}turn-{pitch_mm}mm-{}-{sequence}.{extension}",
            std::process::id()
        ));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&output)
        {
            Ok(file) => break (output, file),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => panic!("create {}: {error}", output.display()),
        }
    };
    file.write_all(bytes)
        .unwrap_or_else(|error| panic!("write {}: {error}", output.display()));
    output
}

fn step_vertex_z_bounds(text: &str) -> (f64, f64) {
    let mut points = HashMap::new();
    let mut vertex_points = Vec::new();
    for entity in text.split(';') {
        let Some((id, definition)) = entity.trim().split_once('=') else {
            continue;
        };
        let id = id.trim();
        if definition.starts_with("CARTESIAN_POINT(") {
            let Some((_, coordinates)) = definition.split_once(",(") else {
                continue;
            };
            let Some(z) = coordinates
                .split(')')
                .next()
                .and_then(|values| values.split(',').nth(2))
                .and_then(|value| value.trim().parse().ok())
            else {
                continue;
            };
            points.insert(id.to_owned(), z);
        } else if definition.starts_with("VERTEX_POINT(")
            && let Some((_, point)) = definition.split_once(",#")
        {
            vertex_points.push(format!("#{}", point.trim_end_matches(')')));
        }
    }
    let mut z_values = vertex_points
        .iter()
        .filter_map(|point| points.get(point).copied());
    let first = z_values.next().expect("STEP has vertex points");
    z_values.fold((first, first), |(min, max), z| (min.min(z), max.max(z)))
}

#[test]
fn three_turn_spring_measures_and_exports() {
    measure_and_export(3, 3, TARGET_WORK);
}

fn measure_and_export(turns: u32, pitch_mm: u32, work_target: u64) {
    let recipe = spring(turns, pitch_mm);
    let (body, build_work) = with_operation_budget(APP_WORK, || {
        numerical_sweeps::build(
            &recipe,
            NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    let body = body.unwrap();
    assert!(!build_work.exhausted && build_work.used <= work_target);

    let (properties, property_work) = with_operation_budget(APP_WORK, || {
        body.properties(&Accuracy::default(), &mut Decider::default())
    });
    let properties = properties.unwrap();
    assert!(
        !property_work.exhausted && property_work.used <= work_target,
        "{property_work:?}"
    );
    assert!(properties.marks.is_empty(), "{:?}", properties.marks);

    let rise = f64::from(turns * pitch_mm);
    let length = f64::from(turns) * (2. * PI * 17.).hypot(f64::from(pitch_mm));
    let ideal_volume = PI * length;
    let ideal_area = 2. * PI * (length + 1.);
    near(&properties.fitted.volume, ideal_volume, ideal_volume / 100.);
    near(&properties.fitted.area, ideal_area, ideal_area / 100.);
    let centroid = properties.fitted.centroid.as_ref().unwrap();
    for (value, expected) in centroid.iter().zip([0., 0., 2.5 + rise / 2.]) {
        near(value, expected, 0.01);
    }

    let (file, export_work) = with_operation_budget(APP_WORK, || {
        let dec = &mut Decider::default();
        let file = step::build(body.body(), &Default::default(), dec)?;
        if let Some(path) = std::env::var_os("CATALLAX_SPRING_STEP") {
            let path = write_step_dump(&path, turns, pitch_mm, file.bytes());
            file.check_written(&fs::read(path).unwrap(), dec)?;
        } else {
            file.check_written(file.bytes(), dec)?;
        }
        Ok::<_, catallax_numbers::Refusal>(file)
    });
    let file = file.unwrap();
    assert!(!export_work.exhausted && export_work.used <= work_target);
    let text = std::str::from_utf8(file.bytes()).unwrap();
    let provenance = &file.report().source_evidence;
    let provenance_stamp = format!(
        "catallax numerical approximate;source={};source_tolerance_mm={};source_checks={};source_repairs={};",
        file.report().source_name,
        provenance.tolerance_mm(),
        provenance.checks().len(),
        provenance.repairs().len(),
    );
    assert_eq!(file.report().source_name, "spring");
    assert!(
        text.contains(&provenance_stamp),
        "missing provenance: {provenance_stamp}"
    );
    assert_eq!(text.matches("MANIFOLD_SOLID_BREP(").count(), 1);
    let (min_z, max_z) = step_vertex_z_bounds(text);
    let mut source_z = file
        .report()
        .source_definition
        .vertices
        .iter()
        .map(|vertex| vertex.position[2].to_f64_nearest());
    let first = source_z.next().expect("spring source has vertices");
    let (source_min_z, source_max_z) =
        source_z.fold((first, first), |(min, max), z| (min.min(z), max.max(z)));
    assert!(
        (min_z - source_min_z).abs() <= 0.01 && (max_z - source_max_z).abs() <= 0.01,
        "STEP vertex z bounds {min_z}..{max_z} differ from source bounds {source_min_z}..{source_max_z}"
    );
    eprintln!("sweep properties {:?}", properties.fitted);
    eprintln!(
        "spring work: build {}, properties {}, STEP and reread {}",
        build_work.used, property_work.used, export_work.used
    );
}

#[test]
fn spring_step_dump_preserves_existing_file_and_uses_unique_paths() {
    let sequence = SPRING_DUMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let requested = std::env::temp_dir().join(format!(
        "catallax-spring-step-{}-{sequence}.step",
        std::process::id()
    ));
    fs::write(&requested, b"existing artifact").unwrap();
    let first = write_step_dump(requested.as_os_str(), 3, 3, b"first STEP");
    let second = write_step_dump(requested.as_os_str(), 3, 3, b"second STEP");

    assert_eq!(fs::read(&requested).unwrap(), b"existing artifact");
    assert_ne!(first, second);
    assert_eq!(fs::read(&first).unwrap(), b"first STEP");
    assert_eq!(fs::read(&second).unwrap(), b"second STEP");

    fs::remove_file(requested).unwrap();
    fs::remove_file(first).unwrap();
    fs::remove_file(second).unwrap();
}

#[test]
#[ignore = "near-allowance spring qualification"]
fn near_limit_spring_measures_and_excess_turn_refuses() {
    measure_and_export(20, 3, APP_WORK);

    let (result, work) = with_operation_budget(APP_WORK, || {
        numerical_sweeps::build(
            &spring(21, 3),
            NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    let Err(catallax_app::Error::Build(catallax_features::BuildError::Geometry {
        refusal, ..
    })) = result
    else {
        panic!("an over-allowance spring must return a geometry refusal");
    };
    assert_eq!(refusal.code, "numbers.operation-budget");
    assert!(work.exhausted);
}
