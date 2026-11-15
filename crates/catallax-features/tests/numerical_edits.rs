//! N04 source -> sew -> dimension edit -> cut -> saved-definition replay.
#[path = "../../catallax-step/tests/common/n03_tube.rs"]
mod fixture;
use catallax_features::{
    evaluate_numerical,
    numerical::{NumericalEdit, NumericalOperation, NumericalRecipe},
};
use catallax_numbers::{BudgetPolicy, Decider, NumericalCheck, with_operation_budget};
use catallax_shapes::{Curve, Surface};
use catallax_topology::numerical::NumericalOptions;
use fixture::{q, tube};

fn recipe() -> NumericalRecipe {
    let mut source = tube();
    let keeper = source.edges[0].vertices[0];
    let removed = source.vertices.len();
    let mut vertex = source.vertices[keeper].clone();
    vertex.key = "split-endpoint".into();
    vertex.position[0] = &vertex.position[0] + &q("0.005");
    source.vertices.push(vertex);
    source.edges[0].vertices[0] = removed;
    NumericalRecipe {
        source,
        sew: vec![[keeper, removed]],
        options: NumericalOptions::default(),
        edits: vec![],
    }
}
fn edit(id: &str, operation: NumericalOperation) -> NumericalEdit {
    NumericalEdit {
        id: id.into(),
        options: NumericalOptions {
            tolerance_mm: q("0.002"),
        },
        operation,
    }
}
fn resize(axis: usize) -> NumericalEdit {
    edit(
        "height",
        NumericalOperation::Resize {
            axis,
            length_mm: q("2"),
        },
    )
}
fn cut(axis: usize, at: &str, keep_lower: bool) -> NumericalEdit {
    edit(
        "trim",
        NumericalOperation::Cut {
            axis,
            at_mm: q(at),
            keep_lower,
        },
    )
}
fn step(body: &catallax_topology::numerical::NumericalBody) -> Vec<u8> {
    catallax_step::numerical::build(
        body,
        &catallax_step::numerical::NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .bytes()
    .to_vec()
}
fn volume(d: &catallax_topology::numerical::NumericalDefinition) -> catallax_numbers::Rational {
    // Independent oriented polygon tetrahedra, including the negative hole loops.
    let mut six = q("0");
    for face in &d.faces {
        for boundary in &face.loops {
            let points: Vec<_> = boundary
                .iter()
                .map(|t| {
                    let e = &d.edges[t.edge];
                    &d.vertices
                        [e.vertices[usize::from(t.sense == catallax_topology::Sense::Reversed)]]
                    .position
                })
                .collect();
            let a = points[0];
            for i in 1..points.len() - 1 {
                let b = points[i];
                let c = points[i + 1];
                six = six
                    + &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
                    + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
                    + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0]);
            }
        }
    }
    six / q("6")
}

#[test]
fn dimension_then_cut_replays_origins_repairs_current_checks_and_step() {
    let original = evaluate_numerical(&recipe(), &mut Decider::default()).unwrap();
    assert_eq!(volume(original.body().definition()), q("12"));
    let resized = original.edit(resize(2), &mut Decider::default()).unwrap();
    assert_eq!(volume(resized.body().definition()), q("24"));
    let final_body = resized
        .edit(cut(2, "1.5", true), &mut Decider::default())
        .unwrap();
    assert_eq!(volume(final_body.body().definition()), q("18"));
    assert_eq!(final_body.body().counts(), original.body().counts());
    assert_eq!(final_body.body().evidence().tolerance_mm(), &q("0.002"));
    assert_eq!(
        final_body.body().evidence().repairs(),
        original.body().evidence().repairs()
    );
    let history = final_body.body().history();
    assert_eq!(history.len(), 2);
    assert_eq!(history[0].evidence, *original.body().evidence());
    assert_eq!(history[0].options.tolerance_mm, q("0.01"));
    assert!(history[0].evidence.checks().iter().any(
        |c| matches!(c,NumericalCheck::Passed{name,..} if name=="sewing.explicit-vertex-pair")
    ));
    for (key, origins) in original.face_origins() {
        assert!(origins.is_subset(&final_body.face_origins()[key]));
    }
    assert!(final_body.face_origins()["f0"].contains("edit:trim/cap"));
    assert_eq!(final_body.body().source_definition().vertices.len(), 17);
    let saved = final_body.recipe().to_json().unwrap();
    let parsed = NumericalRecipe::parse(&saved, &mut Decider::default()).unwrap();
    assert_eq!(parsed.to_json().unwrap(), saved);
    let replay = evaluate_numerical(&parsed, &mut Decider::default()).unwrap();
    assert_eq!(replay.face_origins(), final_body.face_origins());
    assert_eq!(replay.body().evidence(), final_body.body().evidence());
    for (a, b) in replay.body().history().iter().zip(history) {
        assert_eq!(a.evidence, b.evidence);
    }
    assert_eq!(step(replay.body()), step(final_body.body()));
    assert_eq!(volume(original.body().definition()), q("12"));
    if let Ok(dir) = std::env::var("CATALLAX_N04_CAPTURE") {
        std::fs::write(std::path::Path::new(&dir).join("recipe.json"), saved).unwrap();
        std::fs::write(
            std::path::Path::new(&dir).join("edited.step"),
            step(replay.body()),
        )
        .unwrap();
    }
    let mut expected_work = None;
    for _ in 0..3 {
        let (answer, work) = with_operation_budget(20_000_000, || {
            evaluate_numerical(&parsed, &mut Decider::default())
        });
        assert_eq!(answer.unwrap().body().evidence(), replay.body().evidence());
        assert!(!work.exhausted);
        if let Some(expected) = expected_work {
            assert_eq!(work.used, expected);
        } else {
            expected_work = Some(work.used);
        }
    }
    println!(
        "N04 volume=18 hole=retained origins=retained history=2 repairs=1 work={}",
        expected_work.unwrap()
    );
}

#[test]
fn refused_edits_and_replay_are_atomic_and_share_the_request_limit() {
    let original = evaluate_numerical(&recipe(), &mut Decider::default()).unwrap();
    let before = original.recipe().to_json().unwrap();
    let bytes = step(original.body());
    for operation in [
        resize(3),
        cut(2, "0", true),
        cut(2, "2", true),
        edit(
            "zero",
            NumericalOperation::Resize {
                axis: 2,
                length_mm: q("0"),
            },
        ),
    ] {
        assert!(original.edit(operation, &mut Decider::default()).is_err());
        assert_eq!(original.recipe().to_json().unwrap(), before);
        assert_eq!(step(original.body()), bytes);
    }
    let mut raw = recipe();
    let Surface::Plane(p) = &mut raw.source.faces[0].surface else {
        panic!("plane")
    };
    p.origin.z = p.origin.z.clone() + q("0.003").into();
    let noisy = evaluate_numerical(&raw, &mut Decider::default()).unwrap();
    let error = noisy.edit(resize(2), &mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "numbers.numerical-tolerance-exceeded");
    assert!(noisy.body().history().is_empty());
    let mut raw = recipe();
    raw.edits = vec![resize(2), cut(2, "1.5", true)];
    let (answer, work) =
        with_operation_budget(1, || evaluate_numerical(&raw, &mut Decider::default()));
    assert!(answer.is_err());
    assert!(work.exhausted);
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        ..BudgetPolicy::default()
    });
    assert!(evaluate_numerical(&raw, &mut dec).is_err());
    raw.edits.push(resize(2));
    assert_eq!(
        evaluate_numerical(&raw, &mut Decider::default())
            .unwrap_err()
            .code,
        "features.numerical-definition"
    );
    let mut saved: serde_json::Value = serde_json::from_str(&before).unwrap();
    saved["source"]["vertices"][16]["position"][0] = "0.02".into();
    let parsed = NumericalRecipe::parse(&saved.to_string(), &mut Decider::default()).unwrap();
    assert!(evaluate_numerical(&parsed, &mut Decider::default()).is_err());
    assert_eq!(original.recipe().to_json().unwrap(), before);
}

#[test]
fn lower_cut_and_cyclic_coordinate_frame_preserve_the_hole_and_material() {
    // Cycle x,y,z to y,z,x, a proper rotation, then translate all source carriers.
    let mut r = recipe();
    for v in &mut r.source.vertices {
        v.position.rotate_right(1);
        v.position[0] = &v.position[0] + q("10");
    }
    for e in &mut r.source.edges {
        let Curve::Line(l) = &mut e.curve else {
            panic!("line")
        };
        (l.origin.x, l.origin.y, l.origin.z) = (
            l.origin.z.clone() + q("10").into(),
            l.origin.x.clone(),
            l.origin.y.clone(),
        );
        (l.dir.x, l.dir.y, l.dir.z) = (l.dir.z.clone(), l.dir.x.clone(), l.dir.y.clone());
    }
    for f in &mut r.source.faces {
        let Surface::Plane(p) = &mut f.surface else {
            panic!("plane")
        };
        (p.origin.x, p.origin.y, p.origin.z) = (
            p.origin.z.clone() + q("10").into(),
            p.origin.x.clone(),
            p.origin.y.clone(),
        );
        for v in [&mut p.normal, &mut p.ref_dir] {
            (v.x, v.y, v.z) = (v.z.clone(), v.x.clone(), v.y.clone());
        }
    }
    r.edits = vec![resize(0), cut(0, "10.5", false)];
    let body = evaluate_numerical(&r, &mut Decider::default()).unwrap();
    assert_eq!(volume(body.body().definition()), q("18"));
    assert!(body.face_origins()["f1"].contains("edit:trim/cap"));
    assert_eq!(body.body().definition().faces[0].loops.len(), 2);
    assert_eq!(body.body().definition().faces[1].loops.len(), 2);
    let parsed = NumericalRecipe::parse(&r.to_json().unwrap(), &mut Decider::default()).unwrap();
    assert_eq!(
        step(
            evaluate_numerical(&parsed, &mut Decider::default())
                .unwrap()
                .body()
        ),
        step(body.body())
    );
}
