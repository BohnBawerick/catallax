//! The file-aware stage keeps the body source through preview and print.

mod import_support;

use catallax_features::ImportFiles;
use catallax_field::lattice_stage::{
    GyroidSpec, LatticeSpec, Proof, Settings, Stage, TerminalField, default_preview_options,
};
use catallax_numbers::{Decider, Rational, current_operation_work, with_operation_budget};
use catallax_recipe::{
    FileRef, GyroidKind, LatticeCell, Number, Operation, Pose, Recipe, Step, StepId,
};

const TARGET: u64 = 20_000_000;
const ALLOWANCE: u64 = 100_000_000;
const CUBE: &[u8] = b"v -1 -1 -1\nv 1 -1 -1\nv 1 1 -1\nv -1 1 -1\n\
v -1 -1 1\nv 1 -1 1\nv 1 1 1\nv -1 1 1\n\
f 1 4 3\nf 1 3 2\nf 5 6 7\nf 5 7 8\nf 1 2 6\nf 1 6 5\n\
f 2 3 7\nf 2 7 6\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\n";

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn settings(accuracy: &str) -> Settings {
    Settings {
        accuracy: q(accuracy),
        minimum_detail: q("1/2"),
        max_cells: 2_000_000,
        max_vertices: 1_000_000,
    }
}

fn pose(translation: [&str; 3]) -> Pose {
    serde_json::from_value(serde_json::json!({
        "rotation": ["1", "0", "0", "0"],
        "translation": translation,
    }))
    .unwrap()
}

fn lattice() -> LatticeSpec {
    LatticeSpec {
        cell: LatticeCell::Cubic,
        size: Number::parse("4").unwrap(),
        frame: pose(["0", "0", "1"]),
        radius: serde_json::from_str(r#"{"constant":"0.35"}"#).unwrap(),
    }
}

fn gyroid(level: &str) -> GyroidSpec {
    GyroidSpec {
        period: Number::parse("8").unwrap(),
        level: serde_json::from_value(serde_json::json!({"constant": level})).unwrap(),
        kind: GyroidKind::Network,
        frame: Pose::identity(),
    }
}

fn source(bytes: &[u8], step: bool) -> (Recipe, ImportFiles<'_>) {
    let sha256 = catallax_import::sha256_hex(bytes);
    let op = if step {
        Operation::ImportStep {
            file: FileRef {
                name: "trimmed.step".into(),
                sha256: sha256.clone(),
            },
            pose: pose(["1/8", "0", "0"]),
            source_solid: Some("#15".into()),
            decisions: Default::default(),
        }
    } else {
        Operation::ImportMesh {
            file: FileRef {
                name: "closed.obj".into(),
                sha256: sha256.clone(),
            },
            pose: pose(["0", "0", "1/2"]),
            scale: Number::parse("1").unwrap(),
        }
    };
    let id = StepId::new("source").unwrap();
    (
        Recipe::new("Imported solid", vec![Step { id: id.clone(), op }], id),
        ImportFiles::from([(sha256, bytes)]),
    )
}

fn stage_preview<S: TerminalField>(
    body: &Recipe,
    files: &ImportFiles<'_>,
    spec: S,
    settings: Settings,
    dec: &mut Decider,
) -> Stage<S> {
    let stage = Stage::new_with_files(body, files, spec, settings, dec).unwrap();
    assert_eq!(&stage.recipe().steps[..body.steps.len()], &body.steps);
    assert_eq!(stage.recipe().result, *stage.lattice_step());
    let preview = stage
        .preview(&default_preview_options(stage.plan(), dec).unwrap(), dec)
        .unwrap();
    assert_eq!(preview.proof, Proof::Unproved);
    assert!(!preview.mesh.triangles.is_empty());
    if matches!(body.steps[0].op, Operation::ImportStep { .. }) {
        // All placed loft controls satisfy y - x <= 15/8. The preview's
        // 1/16 mm pitch allows at most 1/8 beyond that line in a mixed cell.
        assert!(preview.mesh.positions.iter().all(|p| p[1] - p[0] <= 2.0));
    }
    assert!(current_operation_work().used <= TARGET);
    stage
}

fn stage_print<S: TerminalField>(
    body: &Recipe,
    files: &ImportFiles<'_>,
    spec: S,
    settings: Settings,
) -> (Stage<S>, Vec<u8>) {
    let ((stage, bytes), work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        let stage = stage_preview(body, files, spec, settings.clone(), &mut dec);
        let (bytes, output) = stage.export(&mut dec).unwrap().labelled_stl().unwrap();
        let verdict = output.verdict();
        assert!(verdict.closed && verdict.connected && verdict.non_self_intersecting);
        assert_eq!(verdict.components, 1);
        assert!(verdict.max_deviation <= settings.accuracy);
        assert_eq!(bytes.len(), 84 + 50 * output.mesh().triangles.len());
        assert_eq!(&bytes[..12], b"CXPRINT1\x02\x01\x01\x01");
        (stage, bytes)
    });
    assert!(!work.exhausted);
    println!(
        "{} at {}: {} work",
        stage.lattice_step(),
        stage.settings().accuracy,
        work.used
    );
    (stage, bytes)
}

type Triangle = [[f64; 3]; 3];

fn stl_triangles(bytes: &[u8]) -> Vec<Triangle> {
    let count = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
    assert!(count > 0);
    assert_eq!(bytes.len(), 84 + 50 * count);
    bytes[84..]
        .as_chunks::<50>()
        .0
        .iter()
        .map(|record| {
            std::array::from_fn(|vertex| {
                std::array::from_fn(|axis| {
                    let offset = 12 + 12 * vertex + 4 * axis;
                    f64::from(f32::from_le_bytes(
                        record[offset..offset + 4].try_into().unwrap(),
                    ))
                })
            })
        })
        .collect()
}

fn obeys_step_clip(triangles: &[Triangle], accuracy: f64) -> bool {
    // The placed trim controls reach y - x = 15/8 and z = 0 and 2.
    // Divide by the oblique plane's normal length to get millimetres.
    triangles.iter().flatten().all(|&[x, y, z]| {
        x.is_finite()
            && y.is_finite()
            && (y - x - 15.0 / 8.0) / 2.0_f64.sqrt() <= accuracy
            && (-accuracy..=2.0 + accuracy).contains(&z)
    })
}

fn obeys_cube_clip(triangles: &[Triangle], accuracy: f64) -> bool {
    triangles.iter().flatten().all(|&[x, y, z]| {
        x.abs() <= 1.0 + accuracy && y.abs() <= 1.0 + accuracy && (z - 0.5).abs() <= 1.0 + accuracy
    })
}

fn crossing_error(triangles: &[Triangle], p: [f64; 3], axis: usize) -> f64 {
    let [u, v] = [(axis + 1) % 3, (axis + 2) % 3];
    triangles
        .iter()
        .filter_map(|&[a, b, c]| {
            let det = |x: [f64; 3], y: [f64; 3]| {
                (x[u] - a[u]) * (y[v] - a[v]) - (x[v] - a[v]) * (y[u] - a[u])
            };
            let area = det(b, c);
            if area == 0.0 {
                return None;
            }
            let s = det(p, c) / area;
            let t = det(b, p) / area;
            (s >= 0.0 && t >= 0.0 && s + t <= 1.0).then(|| {
                (a[axis] + s * (b[axis] - a[axis]) + t * (c[axis] - a[axis]) - p[axis]).abs()
            })
        })
        .fold(f64::INFINITY, f64::min)
}

fn cube_lattice_crossings(triangles: &[Triangle]) -> [f64; 3] {
    // Along the vertical beam's axis, both crossings come from the cube.
    // Across a horizontal beam, the cylinder radius supplies its own crossing.
    [
        crossing_error(triangles, [0.0, 0.0, -0.5], 2),
        crossing_error(triangles, [0.0, 0.0, 1.5], 2),
        crossing_error(triangles, [0.75, 0.35, 1.0], 1),
    ]
}

fn matches_cube_lattice(errors: [f64; 3], accuracy: f64) -> bool {
    // The cap windows end halfway to the corresponding unplaced cube plane.
    let cap_window = accuracy.min(0.25);
    let [bottom, top, beam] = errors;
    bottom < cap_window && top < cap_window && beam < accuracy
}

fn cube_gyroid_crossings(triangles: &[Triangle], accuracy: f64) -> bool {
    // Along the z axis the level-zero gyroid is sin(pi*z/4) <= 0.
    // The translated cube clips that material to [-1/2, 0].
    crossing_error(triangles, [0.0, 0.0, -0.5], 2) < accuracy.min(0.25)
        && crossing_error(triangles, [0.0, 0.0, 0.0], 2) < accuracy.min(0.25)
}

fn translate_stl(bytes: &mut [u8], axis: usize, distance: f32) {
    for record in bytes[84..].as_chunks_mut::<50>().0 {
        for vertex in 0..3 {
            let offset = 12 + 12 * vertex + 4 * axis;
            let coordinate = &mut record[offset..offset + 4];
            let moved = f32::from_le_bytes(coordinate.try_into().unwrap()) + distance;
            coordinate.copy_from_slice(&moved.to_le_bytes());
        }
    }
}

#[test]
fn rectangular_step_clips_lattice_to_a_certified_stl() {
    let step = import_support::trimmed_loft();
    let (body, files) = source(step.as_bytes(), true);
    let (_, mut bytes) = stage_print(&body, &files, lattice(), settings("1/2"));
    assert!(obeys_step_clip(&stl_triangles(&bytes), 0.5));
    translate_stl(&mut bytes, 1, 4.0);
    assert!(!obeys_step_clip(&stl_triangles(&bytes), 0.5));
}

#[test]
fn closed_mesh_clips_lattice_and_restages_at_a_new_accuracy() {
    let (body, files) = source(CUBE, false);
    let (first, original) = stage_print(&body, &files, lattice(), settings("1/2"));
    let recipe = first.recipe().clone();
    let (returned, spec, mut options) = first.return_before_lattice();
    assert_eq!(body, returned);
    options.accuracy = q("3/8");
    let (restaged, replacement) = stage_print(&returned, &files, spec, options);
    assert_eq!(*restaged.recipe(), recipe);
    assert_ne!(original, replacement);
    for (mut bytes, accuracy) in [(original, 0.5), (replacement, 0.375)] {
        let triangles = stl_triangles(&bytes);
        assert!(obeys_cube_clip(&triangles, accuracy));
        let errors = cube_lattice_crossings(&triangles);
        assert!(
            matches_cube_lattice(errors, accuracy),
            "{errors:?} at {accuracy}"
        );
        translate_stl(&mut bytes, 2, -0.5);
        assert!(!matches_cube_lattice(
            cube_lattice_crossings(&stl_triangles(&bytes)),
            accuracy,
        ));
    }
}

#[test]
fn rectangular_step_gyroid_preserves_the_disconnected_refusal() {
    let step = import_support::trimmed_loft();
    let (body, files) = source(step.as_bytes(), true);
    let (refusal, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        let stage = stage_preview(&body, &files, gyroid("0"), settings("1/2"), &mut dec);
        stage.export(&mut dec).unwrap_err()
    });
    assert!(!work.exhausted);
    assert_eq!(refusal.code, "field.print-disconnected", "{refusal:?}");
    println!("{} work before refusal: {refusal:?}", work.used);
}

#[test]
fn rectangular_step_clips_a_connected_network_gyroid_to_a_certified_stl() {
    let step = import_support::trimmed_loft();
    let (body, files) = source(step.as_bytes(), true);
    let (_, bytes) = stage_print(&body, &files, gyroid("1/2"), settings("1/2"));
    assert!(obeys_step_clip(&stl_triangles(&bytes), 0.5));
}

#[test]
fn closed_mesh_clips_a_network_gyroid_to_a_certified_stl() {
    let (body, files) = source(CUBE, false);
    let (_, mut bytes) = stage_print(&body, &files, gyroid("0"), settings("1/2"));
    let triangles = stl_triangles(&bytes);
    assert!(obeys_cube_clip(&triangles, 0.5));
    assert!(cube_gyroid_crossings(&triangles, 0.5));
    translate_stl(&mut bytes, 2, -0.5);
    assert!(!cube_gyroid_crossings(&stl_triangles(&bytes), 0.5));
}

#[test]
fn a_refused_trim_keeps_the_plan_refusal() {
    let step = import_support::trimmed_loft().replace("(1.,-2.,-2.)", "(1.01,-2.,-2.)");
    let (body, files) = source(step.as_bytes(), true);
    let mut dec = Decider::default();
    let refusal =
        Stage::new_with_files(&body, &files, lattice(), settings("1/2"), &mut dec).unwrap_err();
    let expected =
        catallax_field::Plan::compile_with_files_checked(&body, &files, &mut dec).unwrap_err();
    assert_eq!(expected.code, "field.import-spline");
    assert_eq!(
        (refusal.code, refusal.detail),
        (expected.code, expected.detail)
    );
}

#[test]
fn compilation_and_preview_share_the_callers_work() {
    let (body, files) = source(CUBE, false);
    let compile =
        |dec: &mut Decider| Stage::new_with_files(&body, &files, lattice(), settings("1/2"), dec);
    let (stage, work) = with_operation_budget(TARGET, || compile(&mut Decider::default()).unwrap());
    assert!(!work.exhausted);
    let controls = default_preview_options(stage.plan(), &mut Decider::default()).unwrap();
    let (result, work) = with_operation_budget(work.used, || {
        let mut dec = Decider::default();
        let stage = compile(&mut dec).unwrap();
        stage.preview(&controls, &mut dec)
    });
    let refusal = result.unwrap_err();
    assert!(work.exhausted, "{refusal:?}");
    assert_eq!(refusal.code, "numbers.operation-budget");
}
