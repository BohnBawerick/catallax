//! Gap measurement and the policy table, through `propose`.

mod common;

use catallax_import::{Disposition, GapKind, ImportOptions, Proposal, propose};
use common::*;

fn shifted(size: [&str; 3], face: usize, by: &str) -> String {
    let mut b = BoxFixture::new(size);
    b.shift = Some((face, [r(by), r("0"), r("0")]));
    box_step(b)
}

fn only(text: &str) -> Proposal {
    let set = propose(text.as_bytes(), &ImportOptions::default()).unwrap();
    assert_eq!(set.format, "catallax-import-proposals/1");
    assert_eq!(set.proposals.len(), 1, "{:?}", set.proposals);
    set.proposals.into_iter().next().unwrap()
}

#[test]
fn a_clean_file_has_no_proposals() {
    let set = propose(
        box_step(BoxFixture::new(["10", "10", "10"])).as_bytes(),
        &ImportOptions::default(),
    )
    .unwrap();
    assert!(set.proposals.is_empty());
    assert_eq!(set.source_sha256.len(), 64);
}

#[test]
fn a_small_seam_is_measured_and_automatic() {
    let p = only(&shifted(["10", "10", "10"], 1, "0.004"));
    assert_eq!(p.id, "g0");
    assert_eq!(p.kind, GapKind::SeamClose);
    assert_eq!(p.disposition, Disposition::Auto);
    assert_eq!(
        (p.distance_min_mm.as_str(), p.distance_max_mm.as_str()),
        ("0.004", "0.004")
    );
    assert_eq!(p.tolerance_mm, "0.01");
    assert!(p.question.starts_with("repair:") && p.question.ends_with(":g0"));
    // One gap for the whole displaced face, naming it and its neighbours' edges.
    assert!(p.entities.len() >= 5, "{:?}", p.entities);
}

#[test]
fn the_cap_is_inclusive_and_just_above_it_asks() {
    let at = only(&shifted(["10", "10", "10"], 1, "0.01"));
    assert_eq!(at.disposition, Disposition::Auto);
    let above = only(&shifted(["10", "10", "10"], 1, "0.010001"));
    assert_eq!(above.disposition, Disposition::Ask);
    assert!(above.reason.contains("cap"));
    assert_eq!(above.tolerance_mm, "0.010001");
}

#[test]
fn a_move_over_a_tenth_of_the_nearest_feature_asks() {
    // The box is 0.05 mm thick in x: 0.009 mm is within the cap but over a tenth of it.
    let p = only(&shifted(["0.05", "10", "10"], 1, "0.009"));
    assert_eq!(p.disposition, Disposition::Ask);
    assert!(p.reason.contains("tenth"), "{}", p.reason);
    let ok = only(&shifted(["0.5", "10", "10"], 1, "0.009"));
    assert_eq!(ok.disposition, Disposition::Auto);
}

#[test]
fn separate_bodies_a_few_microns_apart_have_no_proposal() {
    let mut b = BoxFixture::new(["1", "1", "1"]);
    b.origin = [r("1.005"), r("0"), r("0")];
    let text = boxes_step(&[BoxFixture::new(["1", "1", "1"]), b]);
    assert!(
        propose(text.as_bytes(), &ImportOptions::default())
            .unwrap()
            .proposals
            .is_empty()
    );
}

#[test]
fn a_thin_open_edge_is_an_opening_and_asks() {
    // A plate 0.005 mm thick with its +x edge face missing.
    let mut b = BoxFixture::new(["10", "10", "0.005"]);
    b.omit = Some(1);
    let p = only(&box_step(b));
    assert_eq!(p.kind, GapKind::Opening);
    assert_eq!(p.disposition, Disposition::Ask);
    assert_eq!(p.distance_max_mm, "0.005");
}

#[test]
fn a_far_open_boundary_is_one_opening_with_no_measured_distance() {
    let mut b = BoxFixture::new(["10", "10", "10"]);
    b.omit = Some(1);
    let p = only(&box_step(b));
    assert_eq!(p.kind, GapKind::Opening);
    assert_eq!(p.disposition, Disposition::Ask);
}

#[test]
fn a_gap_wider_than_the_search_radius_is_an_opening_not_a_seam() {
    let text = shifted(["10", "10", "10"], 1, "2");
    let set = propose(text.as_bytes(), &ImportOptions::default()).unwrap();
    assert!(set.proposals.iter().all(|p| p.kind == GapKind::Opening));
}
