//! Preview delivery, revision identity and conserved deferred work.
use catallax_app::{
    Error,
    lifecycle::{Controller, ResolvedSelection, Status, WorkAllocation},
};
use catallax_mesh::{MeshOptions, PreviewOptions};
use catallax_numbers::{Decider, current_operation_work, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::json;

fn recipe(width: u32) -> Recipe {
    catallax_app::soliddraft::convert(&json!({
        "format":"soliddraft-recipe/1", "name":"edit", "steps":[
            {"id":"profile","op":"profile","plane":"XY","rectangle":{"center":[0,0],"width":width,"height":2}},
            {"id":"body","op":"extrude","profile":"profile","distance":3}
        ], "result":"body"
    }).to_string()).unwrap().recipe
}

fn start(
    controller: &Controller,
    width: u32,
    allocation: WorkAllocation,
    proof: &MeshOptions,
) -> catallax_app::lifecycle::ProofJob {
    controller
        .begin()
        .preview(
            &recipe(width),
            &PreviewOptions::default(),
            proof,
            allocation,
            &mut Decider::default(),
        )
        .unwrap()
        .unwrap()
}

fn allocation() -> WorkAllocation {
    WorkAllocation::application(1_000_000).unwrap()
}

#[test]
fn preview_then_confirm_resolves_names_only_on_the_proven_result() {
    let controller = Controller::default();
    let job = start(&controller, 4, allocation(), &MeshOptions::default());
    let preview = controller.latest().unwrap();
    assert!(matches!(preview.status(), Status::Unproved));
    assert!(preview.work().proof.is_none());
    assert_eq!(preview.preview().mesh.positions.len(), 8);
    assert_eq!(preview.preview().mesh.triangles.len(), 12);
    let face = preview.preview().mesh.faces[0].face;
    let edge = preview.preview().mesh.edges[0].edge;
    let vertex = *preview.preview().names.vertices.keys().next().unwrap();
    let selections = [
        preview.select_face(face).unwrap(),
        preview.select_edge(edge).unwrap(),
        preview.select_vertex(vertex).unwrap(),
    ];
    assert!(
        selections
            .iter()
            .all(|pick| preview.resolve(pick).is_none())
    );
    let completed = job.prove(&mut Decider::default());
    // Computing proof does not publish an unexamined or queued completion.
    assert!(matches!(
        controller.latest().unwrap().status(),
        Status::Unproved
    ));
    assert!(completed.publish());
    let confirmed = controller.latest().unwrap();
    assert_eq!(confirmed.revision(), preview.revision());
    let Status::Confirmed(body) = confirmed.status() else {
        panic!("expected confirmed result")
    };
    assert_eq!(body.built.recipe(), &recipe(4));
    assert_eq!(body.names, preview.preview().names);
    assert_eq!(
        confirmed.resolve(&selections[0]),
        Some(ResolvedSelection::Faces(vec![face]))
    );
    assert_eq!(
        confirmed.resolve(&selections[1]),
        Some(ResolvedSelection::Edges(vec![edge]))
    );
    assert_eq!(
        confirmed.resolve(&selections[2]),
        Some(ResolvedSelection::Vertices(vec![vertex]))
    );
    assert!(!confirmed.work().proof.unwrap().exhausted);
    assert!(matches!(preview.status(), Status::Unproved));
}

#[test]
fn typed_proof_failure_marks_the_visible_edit_failed() {
    let controller = Controller::default();
    let job = start(
        &controller,
        4,
        WorkAllocation {
            preview: 1_000_000,
            proof: 0,
        },
        &MeshOptions::default(),
    );
    let before = controller.latest().unwrap();
    assert!(matches!(before.status(), Status::Unproved));
    assert!(job.prove(&mut Decider::default()).publish());
    let failed = controller.latest().unwrap();
    let Status::Failed(Error::Display(refusal)) = failed.status() else {
        panic!("typed numeric refusal required")
    };
    assert_eq!(refusal.code, "numbers.operation-budget");
    assert_eq!(failed.revision(), before.revision());
    assert_eq!(failed.preview().mesh, before.preview().mesh);
    assert!(failed.work().proof.unwrap().exhausted);
    assert_eq!(failed.work().proof.unwrap().used, 0);
}

#[test]
fn late_success_and_failure_cannot_overwrite_a_newer_edit() {
    for fail in [false, true] {
        let controller = Controller::default();
        let mut caps = allocation();
        if fail {
            caps.proof = 0;
        }
        let old_job = start(&controller, 4, caps, &MeshOptions::default());
        let old = controller.latest().unwrap();
        let pick = old.select_face(old.preview().mesh.faces[0].face).unwrap();
        let completion = old_job.prove(&mut Decider::default());
        let new_job = start(&controller, 8, allocation(), &MeshOptions::default());
        let new = controller.latest().unwrap();
        assert!(new.revision() > old.revision());
        assert!(!completion.publish());
        assert_eq!(controller.latest().unwrap().revision(), new.revision());
        assert!(matches!(
            controller.latest().unwrap().status(),
            Status::Unproved
        ));
        assert!(new_job.prove(&mut Decider::default()).publish());
        assert!(controller.latest().unwrap().resolve(&pick).is_none());
    }
}

#[test]
fn invalidation_precedes_construction_and_cancels_unstarted_proof() {
    let controller = Controller::default();
    let old_job = start(&controller, 4, allocation(), &MeshOptions::default());
    let next = controller.begin();
    assert!(controller.latest().is_none());
    assert!(!old_job.prove(&mut Decider::default()).publish());
    let invalid = PreviewOptions {
        max_vertices: 0,
        ..PreviewOptions::default()
    };
    assert!(
        next.preview(
            &recipe(8),
            &invalid,
            &MeshOptions::default(),
            allocation(),
            &mut Decider::default()
        )
        .is_err()
    );
    assert!(controller.latest().is_none());
}

#[test]
fn equal_revision_numbers_from_different_controllers_do_not_share_picks() {
    let first = Controller::default();
    let second = Controller::default();
    let first_job = start(&first, 4, allocation(), &MeshOptions::default());
    let first_preview = first.latest().unwrap();
    let pick = first_preview
        .select_face(first_preview.preview().mesh.faces[0].face)
        .unwrap();
    let second_job = start(&second, 4, allocation(), &MeshOptions::default());
    assert!(first_job.prove(&mut Decider::default()).publish());
    assert!(second_job.prove(&mut Decider::default()).publish());
    assert_eq!(
        first.latest().unwrap().revision(),
        second.latest().unwrap().revision()
    );
    assert!(first.latest().unwrap().resolve(&pick).is_some());
    assert!(second.latest().unwrap().resolve(&pick).is_none());
}

#[test]
fn caller_reserves_once_and_deferred_proof_cannot_take_a_fresh_allowance() {
    let controller = Controller::default();
    let caps = WorkAllocation {
        preview: 1_000_000,
        proof: 1_000_000,
    };
    let (job, work) = with_operation_budget(2_000_000, || {
        let job = start(&controller, 4, caps, &MeshOptions::default());
        assert_eq!(current_operation_work().used, 2_000_000);
        job
    });
    assert_eq!(work.used, 2_000_000);
    assert!(!work.exhausted);
    // This job already owns its reservation; the originating scope can end.
    assert!(job.prove(&mut Decider::default()).publish());
    let report = controller.latest().unwrap().work();
    assert!(report.preview.used <= caps.preview);
    assert!(report.proof.unwrap().used <= caps.proof);
    assert!(!report.proof.unwrap().exhausted);
    let (_, duplicate) = with_operation_budget(2_000_000, || {
        let controller = Controller::default();
        let _job = start(&controller, 4, caps, &MeshOptions::default());
        let second = Controller::default().begin().preview(
            &recipe(4),
            &PreviewOptions::default(),
            &MeshOptions::default(),
            caps,
            &mut Decider::default(),
        );
        assert!(
            matches!(second, Err(Error::Allocation(refusal)) if refusal.code == "numbers.operation-budget")
        );
    });
    assert_eq!(duplicate.used, 2_000_000);
    assert!(duplicate.exhausted);
}

#[test]
fn exact_phase_caps_pass_and_one_unit_short_fails_without_borrowing() {
    let controller = Controller::default();
    let job = start(&controller, 4, allocation(), &MeshOptions::default());
    assert!(job.prove(&mut Decider::default()).publish());
    let report = controller.latest().unwrap().work();
    let caps = WorkAllocation {
        preview: report.preview.used,
        proof: report.proof.unwrap().used,
    };
    let exact = Controller::default();
    let job = start(&exact, 4, caps, &MeshOptions::default());
    assert!(job.prove(&mut Decider::default()).publish());
    assert!(matches!(
        exact.latest().unwrap().status(),
        Status::Confirmed(_)
    ));
    let short = Controller::default();
    let job = start(
        &short,
        4,
        WorkAllocation {
            proof: caps.proof - 1,
            ..allocation()
        },
        &MeshOptions::default(),
    );
    assert!(job.prove(&mut Decider::default()).publish());
    assert!(matches!(
        short.latest().unwrap().status(),
        Status::Failed(_)
    ));
    let preview_short = Controller::default();
    let answer = preview_short.begin().preview(
        &recipe(4),
        &PreviewOptions::default(),
        &MeshOptions::default(),
        WorkAllocation {
            preview: caps.preview - 1,
            proof: 19_000_000,
        },
        &mut Decider::default(),
    );
    assert!(answer.is_err());
    assert!(preview_short.latest().is_none());
}

#[test]
fn zero_caller_allowance_and_overflow_never_publish() {
    for caps in [
        allocation(),
        WorkAllocation {
            preview: u64::MAX,
            proof: 1,
        },
    ] {
        let controller = Controller::default();
        let (answer, work) = with_operation_budget(0, || {
            controller.begin().preview(
                &recipe(4),
                &PreviewOptions::default(),
                &MeshOptions::default(),
                caps,
                &mut Decider::default(),
            )
        });
        assert!(matches!(answer, Err(Error::Allocation(_))));
        assert_eq!(work.used, 0);
        assert!(controller.latest().is_none());
    }
}
