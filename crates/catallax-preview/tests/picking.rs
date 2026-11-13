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
        let (resolved, proof_work) = with_operation_budget(100_000_000, || {
            reference.clone().resolve(
                &prepared,
                &samples,
                &NearestOptions::default(),
                &mut Decider::default(),
            )
        });
        let resolved = resolved.unwrap();
        assert!(!proof_work.exhausted);
        let (mut proved, mut refused, mut raw_mismatches, mut accepted_mismatches) = (0, 0, 0, 0);
        for (index, (point, pick)) in samples.iter().zip(&reference.points).enumerate() {
            let [x, y, z] = point.map(|v| Rational::from_f64_exact(f64::from(v as f32)).unwrap());
            let (answer, work) = with_operation_budget(100_000_000, || {
                prepared.nearest_face(
                    &RPoint3::new(x, y, z).to_real(),
                    &NearestOptions::default(),
                    &mut Decider::default(),
                )
            });
            assert!(!work.exhausted);
            match answer.unwrap() {
                Ownership::Proved(owner) => {
                    let actual = resolved
                        .face(index)
                        .expect("every proved reference sample must have an accepted owner");
                    if plan.faces()[pick.face.index()].origin != owner.origin {
                        raw_mismatches += 1;
                    }
                    if plan.faces()[actual.index()].origin != owner.origin {
                        accepted_mismatches += 1;
                    }
                    assert_eq!(
                        plan.faces()[actual.index()].origin,
                        owner.origin,
                        "{name} at {point:?}"
                    );
                    let Ownership::Proved(actual) = &resolved.ownership[index] else {
                        panic!("proved reference sample refused");
                    };
                    assert_eq!(actual.lower, owner.lower);
                    assert_eq!(actual.upper, owner.upper);
                    proved += 1;
                }
                Ownership::Refused(expected) => {
                    assert!(
                        resolved.face(index).is_none(),
                        "{name}: an ambiguous sample must have no accepted name"
                    );
                    let Ownership::Refused(actual) = &resolved.ownership[index] else {
                        panic!("ambiguous reference sample accepted");
                    };
                    assert_eq!(*actual, expected);
                    assert_eq!(actual.refusal().code, "field.ownership-ambiguous");
                    assert!(!actual.contenders.is_empty());
                    eprintln!(
                        "{name} refused at {point:?}: {actual:?}; work={}",
                        work.used
                    );
                    refused += 1;
                }
            }
        }
        assert_eq!(proved + refused, samples.len());
        assert_eq!(accepted_mismatches, 0);
        eprintln!(
            "{name}: {} triangles; {} missing cells; {} omitted; chart work={}; samples={}; proved={proved}; refused={refused}; raw triangle disagreements={raw_mismatches}; accepted mismatches={accepted_mismatches}; resolution work={}; 1/2/4 processor candidates agree",
            mesh.triangles.len(),
            mesh.unresolved_cells,
            reference.omitted_triangles,
            work.used,
            samples.len(),
            proof_work.used
        );
    }
}

fn ownership_fingerprint(answer: &Ownership) -> String {
    match answer {
        Ownership::Proved(o) => format!(
            "{:?} face={} [{:?},{:?}] witness={:?} direction={:?} margin={:?} cells={}",
            o.origin,
            o.face.index(),
            o.lower,
            o.upper,
            o.witness.as_rational(),
            o.direction,
            o.margin,
            o.cells
        ),
        Ownership::Refused(a) => format!("refused {a:?}"),
    }
}

#[test]
fn exact_resolution_corrects_candidates_and_refuses_ties_at_every_worker_count() {
    let plan = Plan::compile(
        &Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap(),
    )
    .unwrap();
    let prepared = plan.prepare_distances(&mut Decider::default()).unwrap();
    prepared.prepare_nearest(&mut Decider::default()).unwrap();
    let picker = Picker::new(&fixture()).unwrap();
    let points = [[5.0, 6.0, -3.0], [0.0, 0.0, 0.0]];
    let mut expected = None;
    for workers in [1, 2, 4] {
        let candidates = picker.points(&points, None, workers, 100).unwrap();
        let (resolved, work) = with_operation_budget(100_000_000, || {
            candidates.resolve(
                &prepared,
                &points,
                &NearestOptions::default(),
                &mut Decider::default(),
            )
        });
        let resolved = resolved.unwrap();
        assert!(!work.exhausted);
        let face = resolved.face(0).unwrap();
        assert_eq!(
            plan.faces()[face.index()].origin.role,
            catallax_recipe::Role::Start
        );
        assert!(resolved.face(1).is_none());
        let Ownership::Refused(ambiguity) = &resolved.ownership[1] else {
            panic!("a corner must refuse");
        };
        assert_eq!(ambiguity.refusal().code, "field.ownership-ambiguous");
        assert!(ambiguity.contenders.len() >= 2);
        assert!(
            ambiguity
                .contenders
                .iter()
                .all(|c| c.lower == Rational::zero()
                    && c.upper.as_ref().is_some_and(|u| u == &Rational::zero()))
        );
        let fingerprint = (
            resolved
                .ownership
                .iter()
                .map(ownership_fingerprint)
                .collect::<Vec<_>>(),
            work.used,
        );
        match &expected {
            None => expected = Some(fingerprint),
            Some(expected) => assert_eq!(&fingerprint, expected),
        }
    }
    let candidates = picker.points(&points, None, 1, 100).unwrap();
    let (answer, work) = with_operation_budget(1000, || {
        candidates.resolve(
            &prepared,
            &points,
            &NearestOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted);
    let refusal = answer.unwrap_err();
    assert_eq!(refusal.code, "numbers.operation-budget");
    assert!(refusal.detail.contains("point 0 of 2"));
    let candidates = picker.points(&points, None, 1, 100).unwrap();
    assert_eq!(
        candidates
            .resolve(
                &prepared,
                &points[..1],
                &NearestOptions::default(),
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
        "preview.ownership-points"
    );
}
