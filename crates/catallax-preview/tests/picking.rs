//! General display picking from the field's own named chart triangles.
use catallax_field::{
    Plan,
    display::{PreviewMesh, PreviewOptions, Triangle, preview_plan},
    distance::{NearestOptions, Ownership},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_preview::{gpu::Gpu, picking::Picker};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;

fn fixture() -> PreviewMesh {
    let plan = Plan::compile(
        &Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap(),
    )
    .unwrap();
    PreviewMesh {
        positions: vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 4.0],
            [2.0, 0.0, 4.0],
            [0.0, 2.0, 4.0],
        ],
        triangles: vec![
            Triangle {
                vertices: [0, 1, 2],
                face: plan.result_face(0).unwrap(),
            },
            Triangle {
                vertices: [3, 4, 5],
                face: plan.result_face(1).unwrap(),
            },
            Triangle {
                vertices: [0, 0, 0],
                face: plan.result_face(0).unwrap(),
            },
        ],
        cells: 7,
        unresolved_cells: 3,
        skipped_checks: 0,
    }
}

#[test]
fn nearest_triangle_controls_keep_names_gaps_and_worker_bits() {
    let mesh = fixture();
    let picker = Picker::new(&mesh).unwrap();
    let points = [
        [0.5, 0.5, -1.0],
        [0.5, 0.5, 5.0],
        [3.0, 0.0, 0.0],
        [0.5, 0.5, 2.0],
    ];
    let selected = picker.points(&points, None, 1, 100).unwrap();
    assert_eq!(selected.unresolved_cells, 3);
    assert_eq!(selected.omitted_triangles, 1);
    assert_eq!(
        selected.points.iter().map(|p| p.face).collect::<Vec<_>>(),
        vec![
            mesh.triangles[0].face,
            mesh.triangles[1].face,
            mesh.triangles[0].face,
            mesh.triangles[0].face
        ]
    );
    assert_eq!(
        selected
            .points
            .iter()
            .map(|p| p.distance)
            .collect::<Vec<_>>(),
        vec![1.0, 1.0, 1.0, 2.0]
    );
    for workers in [2, 4] {
        assert_eq!(
            picker.points(&points, None, workers, 100).unwrap().points,
            selected.points
        );
    }
    assert!(picker.points(&[], None, 4, 0).unwrap().points.is_empty());
    assert_eq!(
        picker.points(&points, None, 0, 100).unwrap_err().code,
        "preview.ownership-workers"
    );
    assert_eq!(
        picker.points(&points, None, 1, 7).unwrap_err().code,
        "preview.ownership-size"
    );
    assert_eq!(
        picker
            .points(&[[f64::NAN, 0.0, 0.0]], None, 1, 100)
            .unwrap_err()
            .code,
        "preview.ownership-range"
    );
    let mut broken = mesh.clone();
    broken.triangles[0].vertices[0] = 100;
    assert_eq!(
        Picker::new(&broken).err().unwrap().code,
        "preview.ownership-mesh"
    );
    let mut empty = mesh;
    empty.triangles.clear();
    assert_eq!(
        Picker::new(&empty).err().unwrap().code,
        "preview.ownership-empty"
    );
}

#[test]
fn picking_kernel_validates_without_gpu() {
    let source = include_str!("../src/picking.wgsl");
    let module = naga::front::wgsl::parse_str(source)
        .unwrap_or_else(|e| panic!("{}", e.emit_to_string(source)));
    naga::valid::Validator::new(
        naga::valid::ValidationFlags::all(),
        naga::valid::Capabilities::empty(),
    )
    .validate(&module)
    .unwrap();
}

#[test]
fn target_field_meshes_pick_named_faces_with_exact_f8_checks() {
    let gpu = match Gpu::new(false) {
        Ok(gpu) => Some(gpu),
        Err(reason) => {
            assert!(
                std::env::var_os("CATALLAX_REQUIRE_GPU").is_none(),
                "{reason}"
            );
            eprintln!("GPU picking skipped: {reason}");
            None
        }
    };
    for (name, text) in [
        ("P1", include_str!("../../../parts/P1-bracket.json")),
        ("P2", include_str!("../../../parts/P2-collar.json")),
        ("P3", include_str!("../../../parts/P3-filleted-boss.json")),
        ("P4", include_str!("../../../parts/P4-shelled-tray.json")),
    ] {
        let plan = Plan::compile(&Recipe::parse(text).unwrap()).unwrap();
        let (mesh, work) = with_operation_budget(20_000_000, || {
            preview_plan(
                &plan,
                &PreviewOptions {
                    tolerance: Rational::from_ratio(1, 2).unwrap(),
                    ..PreviewOptions::default()
                },
                &mut Decider::default(),
            )
        });
        let mesh = mesh.unwrap_or_else(|e| panic!("{name}: {e:?}; work={}", work.used));
        let picker = Picker::new(&mesh).unwrap();
        let samples: Vec<_> = mesh
            .triangles
            .iter()
            .step_by(mesh.triangles.len().div_ceil(16))
            .take(16)
            .map(|t| {
                std::array::from_fn(|axis| {
                    t.vertices
                        .iter()
                        .map(|&i| mesh.positions[i][axis])
                        .sum::<f64>()
                        / 3.0
                })
            })
            .collect();
        let reference = picker.points(&samples, None, 1, 100_000_000).unwrap();
        for workers in [2, 4] {
            assert_eq!(
                picker
                    .points(&samples, None, workers, 100_000_000)
                    .unwrap()
                    .points,
                reference.points,
                "{name} workers={workers}"
            );
        }
        if let Some(gpu) = &gpu {
            let actual = picker.points(&samples, Some(gpu), 1, 100_000_000).unwrap();
            let repeat = picker.points(&samples, Some(gpu), 4, 100_000_000).unwrap();
            assert_eq!(actual.points, repeat.points, "{name} GPU repeat");
            for (a, b) in actual.points.iter().zip(&reference.points) {
                assert_eq!(a.face, b.face, "{name}");
            }
        }
        let prepared = plan.prepare_distances(&mut Decider::default()).unwrap();
        prepared.prepare_nearest(&mut Decider::default()).unwrap();
        let (mut proved, mut refused, mut mismatches) = (0, 0, 0);
        for (point, pick) in samples.iter().zip(&reference.points) {
            let [x, y, z] = point.map(|v| Rational::from_f64_exact(f64::from(v as f32)).unwrap());
            let (answer, work) = with_operation_budget(100_000_000, || {
                let point = RPoint3::new(x, y, z).to_real();
                let owner = prepared.nearest_face(
                    &point,
                    &NearestOptions::default(),
                    &mut Decider::default(),
                )?;
                let selected = prepared.point_to_face(
                    pick.face,
                    &point,
                    &catallax_shapes::distance::Options::default(),
                    &mut Decider::default(),
                )?;
                Ok::<_, catallax_numbers::Refusal>((owner, selected))
            });
            let (owner, selected) = answer.unwrap();
            assert!(!work.exhausted);
            // A chart triangle can cross an approximate trim border. Its name
            // must still have a certified witness within this visual target.
            assert!(
                selected
                    .upper
                    .as_ref()
                    .is_some_and(|u| u <= &Rational::from_ratio(1, 2).unwrap()),
                "{name} selected face {:?} at {point:?}: {selected:?}",
                pick.face
            );
            match owner {
                Ownership::Proved(owner) => {
                    if plan.faces()[pick.face.index()].origin != owner.origin {
                        mismatches += 1;
                        eprintln!(
                            "{name} unproved picture mismatch at {point:?}: selected={:?}; exact={:?}",
                            plan.faces()[pick.face.index()].origin,
                            owner
                        );
                    }
                    proved += 1;
                }
                Ownership::Refused(ambiguity) => {
                    eprintln!(
                        "{name} F8 refused at {point:?}: {ambiguity:?}; work={}",
                        work.used
                    );
                    refused += 1;
                }
            }
        }
        assert!(
            proved > mismatches,
            "{name}: no picture sample agrees with a proved owner"
        );
        eprintln!(
            "{name}: {} triangles; {} missing cells; {} omitted; chart work={}; exact F8 {proved} proved, {refused} refused, {mismatches} picture mismatches; 1/2/4 processor picks agree",
            mesh.triangles.len(),
            mesh.unresolved_cells,
            reference.omitted_triangles,
            work.used
        );
    }
}
