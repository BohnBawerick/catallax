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
    let display = preview.preview().brep().unwrap();
    let face = display.mesh.faces[0].face;
    let edge = display.mesh.edges[0].edge;
    let vertex = *display.names.vertices.keys().next().unwrap();
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
    assert_eq!(body.checkpoint.built().recipe(), &recipe(4));
    assert_eq!(body.names, display.names);
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
    assert_eq!(
        failed.preview().brep().unwrap().mesh,
        before.preview().brep().unwrap().mesh,
    );
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
        let pick = old
            .select_face(old.preview().brep().unwrap().mesh.faces[0].face)
            .unwrap();
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
        .select_face(first_preview.preview().brep().unwrap().mesh.faces[0].face)
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

#[test]
fn failed_edit_restores_last_proved_display_but_cannot_export_it_as_current() {
    let controller = Controller::default();
    let first = start(&controller, 4, allocation(), &MeshOptions::default());
    assert!(first.prove(&mut Decider::default()).publish());
    let proved = controller.displayed().unwrap();
    assert_eq!(
        controller.confirmed().unwrap().checkpoint.built().recipe(),
        &recipe(4)
    );

    let failing = start(
        &controller,
        8,
        WorkAllocation {
            preview: 1_000_000,
            proof: 0,
        },
        &MeshOptions::default(),
    );
    let requested = controller.latest().unwrap().revision();
    assert_eq!(controller.displayed().unwrap().revision(), requested);
    assert!(controller.confirmed().is_none());
    assert!(failing.prove(&mut Decider::default()).publish());
    assert!(matches!(
        controller.latest().unwrap().status(),
        Status::Failed(_)
    ));
    assert_eq!(controller.latest().unwrap().revision(), requested);
    assert_eq!(
        controller.displayed().unwrap().revision(),
        proved.revision()
    );
    assert!(controller.confirmed().is_none());

    let restored = controller.cancel().unwrap();
    assert_eq!(restored.revision(), proved.revision());
    assert_eq!(
        controller.confirmed().unwrap().checkpoint.built().recipe(),
        &recipe(4)
    );
}

#[test]
fn cancellation_rejects_queued_success_and_failure_without_losing_last_proof() {
    for fail in [false, true] {
        let controller = Controller::default();
        let first = start(&controller, 4, allocation(), &MeshOptions::default());
        assert!(first.prove(&mut Decider::default()).publish());
        let proved = controller.latest().unwrap();
        let mut caps = allocation();
        if fail {
            caps.proof = 0;
        }
        let pending = start(&controller, 8, caps, &MeshOptions::default());
        let completion = pending.prove(&mut Decider::default());
        assert_eq!(controller.cancel().unwrap().revision(), proved.revision());
        assert!(!completion.publish());
        assert_eq!(
            controller.confirmed().unwrap().checkpoint.built().recipe(),
            &recipe(4)
        );

        let next = controller.begin();
        assert!(controller.confirmed().is_none());
        assert_eq!(
            controller.displayed().unwrap().revision(),
            proved.revision()
        );
        drop(next);
        controller.cancel();
        let newest = start(&controller, 12, allocation(), &MeshOptions::default());
        assert!(newest.prove(&mut Decider::default()).publish());
        assert_eq!(
            controller.confirmed().unwrap().checkpoint.built().recipe(),
            &recipe(12)
        );
    }
    let empty = Controller::default();
    let unproved = start(&empty, 4, allocation(), &MeshOptions::default());
    assert!(empty.cancel().is_none());
    assert!(!unproved.prove(&mut Decider::default()).publish());
    assert!(empty.displayed().is_none());
    assert!(empty.confirmed().is_none());
}

#[test]
fn early_construction_delivers_edited_geometry_and_origins_before_exact_proof() {
    let controller = Controller::default();
    let mut input = recipe(4);
    let job = controller
        .begin()
        .construct(
            &input,
            &catallax_field::display::PreviewOptions::default(),
            &MeshOptions::default(),
            allocation(),
            &mut Decider::default(),
        )
        .unwrap()
        .unwrap();
    let pending = controller.latest().unwrap();
    let display = pending.preview().construction().unwrap();
    let bounds = display.mesh.positions.iter().fold(
        ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
        |(mut lo, mut hi), point| {
            for i in 0..3 {
                lo[i] = lo[i].min(point[i]);
                hi[i] = hi[i].max(point[i]);
            }
            (lo, hi)
        },
    );
    assert_eq!(bounds, ([-2.0, -1.0, 0.0], [2.0, 1.0, 3.0]));
    let origin = display.origins.values().next().unwrap();
    let selected = pending.select_origin(origin).unwrap();
    assert!(pending.resolve(&selected).is_none());
    assert!(controller.confirmed().is_none());
    input = recipe(12);
    assert!(job.prove(&mut Decider::default()).publish());
    let proved = controller.latest().unwrap();
    let Status::Confirmed(body) = proved.status() else {
        panic!("exact model must prove")
    };
    assert_eq!(body.checkpoint.built().recipe(), &recipe(4));
    assert_ne!(body.checkpoint.built().recipe(), &input);
    let expected = body
        .names
        .faces
        .iter()
        .filter_map(|(&face, origins)| origins.contains(origin).then_some(face))
        .collect::<Vec<_>>();
    assert_eq!(
        proved.resolve(&selected),
        Some(ResolvedSelection::Faces(expected))
    );
}

#[test]
fn incomplete_construction_preview_never_replaces_a_proved_shape() {
    let controller = Controller::default();
    let first = start(&controller, 4, allocation(), &MeshOptions::default());
    assert!(first.prove(&mut Decider::default()).publish());
    let proved = controller.latest().unwrap();
    let refused = controller.begin().construct(
        &recipe(8),
        &catallax_field::display::PreviewOptions {
            max_cells: 1,
            ..Default::default()
        },
        &MeshOptions::default(),
        allocation(),
        &mut Decider::default(),
    );
    assert!(matches!(
        refused,
        Err(Error::ConstructionPreview(refusal)) if refusal.code == "application.preview-incomplete"
    ));
    assert!(controller.latest().is_none());
    assert!(controller.confirmed().is_none());
    assert_eq!(
        controller.displayed().unwrap().revision(),
        proved.revision()
    );
}
