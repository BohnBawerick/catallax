//! Check once (roadmap item 3.0A step 5; Checkpoint 2 decisions, section 4.1
//! item 5): admission of a definition read back from a file that records its
//! tolerance. "Reopening trusts the tolerance recorded at export."
#[path = "common/analytic.rs"]
mod fixture;
use catallax_numbers::{Decider, NumericalCheck, with_operation_budget};
use catallax_shapes::Surface;
use catallax_topology::{
    BodyKind,
    numerical::{
        self as n, NumericalOptions,
        tolerance::{RECORDED_USE, RECORDED_VERTEX},
    },
};
use fixture::{q, quarter};

fn options(cap: &str) -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q(cap),
        vertex_tolerance_mm: None,
    }
}

/// The quarter cylinder with its wall this far off its edges.
fn gapped(radius: &str) -> n::NumericalDefinition {
    let mut definition = quarter();
    let Surface::Cylinder(wall) = &mut definition.faces[4].surface else {
        unreachable!()
    };
    wall.radius_sq = (&q(radius) * &q(radius)).into();
    definition
}

fn sample_names(body: &n::NumericalBody) -> Vec<&str> {
    body.evidence()
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples { name, .. } => Some(name.as_str()),
            NumericalCheck::Passed { .. } => None,
        })
        .collect()
}

#[test]
fn a_recorded_tolerance_stands_for_every_distance() {
    // The wall is 0.004 mm off its edges. Measured, its edges store 0.0044.
    let (measured, measuring) = with_operation_budget(20_000_000, || {
        n::build(gapped("1.004"), options("0.01"), &mut Decider::default())
    });
    let measured = measured.unwrap();
    assert!(measured.tolerances().largest_mm() > q("0.004"));

    // The file it came from records 0.005. Nothing is measured: a wall 0.004
    // off, or any other, is taken at the record.
    let recorded = q("0.005");
    let (body, trusting) = with_operation_budget(20_000_000, || {
        n::admit_recorded(
            BodyKind::Solid,
            gapped("1.004"),
            Vec::new(),
            options("0.01"),
            &recorded,
            &mut Decider::default(),
        )
    });
    let body = body.unwrap();
    let names = sample_names(&body);
    let uses: usize = body
        .definition()
        .faces
        .iter()
        .map(|face| face.loops.iter().map(Vec::len).sum::<usize>())
        .sum();
    assert_eq!(
        names.iter().filter(|name| **name == RECORDED_USE).count(),
        uses
    );
    assert_eq!(
        names
            .iter()
            .filter(|name| **name == RECORDED_VERTEX)
            .count(),
        body.definition().edges.len()
    );
    assert_eq!(names.len(), uses + body.definition().edges.len());
    // Rule T4's 1.1 belongs to a fresh measurement. The record has it
    // already: an edge holds 0.005 and so does a vertex, not 0.0055.
    for entity in &body.tolerances().entities {
        assert_eq!(entity.tolerance_mm, recorded, "{entity:?}");
    }
    assert!(body.tolerances().warnings.is_empty());
    assert_eq!(body.options().tolerance_mm, q("0.01"));
    assert!(
        trusting.used < measuring.used,
        "recorded {} and measured {}",
        trusting.used,
        measuring.used
    );
}

#[test]
fn a_record_over_the_callers_cap_raises_the_bodys_limit() {
    // Measured under a cap of 0.001 mm this wall refuses. Its file records
    // 0.005 mm, so the reopened body's limit is 0.005 mm: a later operation
    // measures again and is not refused for what the file already held.
    let refusal = n::build(gapped("1.004"), options("0.001"), &mut Decider::default()).unwrap_err();
    assert_eq!(refusal.code, "numbers.numerical-tolerance-exceeded");
    let body = n::admit_recorded(
        BodyKind::Solid,
        gapped("1.004"),
        Vec::new(),
        options("0.001"),
        &q("0.005"),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(body.options().tolerance_mm, q("0.005"));
    let again = body
        .rebuild(
            body.definition().clone(),
            body.options().clone(),
            &mut Decider::default(),
        )
        .unwrap();
    // The continuation measured: its records are distances again.
    assert!(sample_names(&again).contains(&"trim.edge-surface-distance"));
    assert!(!sample_names(&again).contains(&RECORDED_USE));
}

#[test]
fn a_record_under_rule_t1s_floor_leaves_the_floor() {
    let body = n::admit_recorded(
        BodyKind::Solid,
        quarter(),
        Vec::new(),
        options("0.01"),
        &q("0.00000001"),
        &mut Decider::default(),
    )
    .unwrap();
    for entity in &body.tolerances().entities {
        assert_eq!(entity.tolerance_mm, body.tolerances().floor_mm);
    }
}

#[test]
fn the_structure_is_still_checked_and_a_record_must_be_positive() {
    // A missing use leaves the shell open. No record makes that a body.
    let mut open = quarter();
    open.faces[0].loops[0].pop();
    let refusal = n::admit_recorded(
        BodyKind::Solid,
        open,
        Vec::new(),
        options("0.01"),
        &q("0.005"),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_ne!(refusal.code, "numbers.numerical-tolerance-exceeded");
    for record in ["0", "-0.001"] {
        let refusal = n::admit_recorded(
            BodyKind::Solid,
            quarter(),
            Vec::new(),
            options("0.01"),
            &q(record),
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(refusal.code, "numerical.invalid-definition");
    }
}
