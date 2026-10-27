//! Exact self-intersection checking of source and repaired solids.

mod common;

use catallax_import::{
    ImportOptions, ImportRefusal, apply, check_self_intersection, import_step_bodies, read_step,
    sha256_hex,
};
use catallax_recipe::StepId;
use catallax_topology::Model;
use common::*;

fn check(solids: &[Vec<Ring>]) -> Result<(), ImportRefusal> {
    check_self_intersection(&read_step(polyhedra_step(solids).as_bytes()).unwrap())
}

fn crossing(r: Result<(), ImportRefusal>) -> bool {
    matches!(r, Err(ImportRefusal::SelfIntersecting { .. }))
}

/// An L-shaped prism: a 2 x 2 square with the corner [1,2] x [1,2] removed, 1 thick.
fn l_prism() -> Vec<Ring> {
    let outline = [
        ("0", "0"),
        ("2", "0"),
        ("2", "1"),
        ("1", "1"),
        ("1", "2"),
        ("0", "2"),
    ];
    let bottom: Ring = outline.iter().rev().map(|(x, y)| [*x, *y, "0"]).collect();
    let top: Ring = outline.iter().map(|(x, y)| [*x, *y, "1"]).collect();
    let mut faces = vec![bottom, top];
    for k in 0..outline.len() {
        let ((x0, y0), (x1, y1)) = (outline[k], outline[(k + 1) % outline.len()]);
        faces.push(vec![
            [x0, y0, "0"],
            [x1, y1, "0"],
            [x1, y1, "1"],
            [x0, y0, "1"],
        ]);
    }
    faces
}

#[test]
fn clean_solids_do_not_false_positive() {
    assert!(check(&[box_rings(["0", "0", "0"], ["1", "2", "3"])]).is_ok());
    assert!(check(&[l_prism()]).is_ok());
}

#[test]
fn coplanar_faces_sharing_only_an_edge_are_not_an_intersection() {
    // A box whose top is split into two coplanar faces along x = 1.
    let mut faces = box_rings(["0", "0", "0"], ["2", "1", "1"]);
    faces.pop();
    faces.push(vec![
        ["0", "0", "1"],
        ["1", "0", "1"],
        ["1", "1", "1"],
        ["0", "1", "1"],
    ]);
    faces.push(vec![
        ["1", "0", "1"],
        ["2", "0", "1"],
        ["2", "1", "1"],
        ["1", "1", "1"],
    ]);
    assert!(check(&[faces]).is_ok());
}

#[test]
fn touching_at_one_vertex_is_not_an_intersection() {
    let mut faces = box_rings(["0", "0", "0"], ["1", "1", "1"]);
    faces.extend(box_rings(["1", "1", "1"], ["2", "2", "2"]));
    assert!(check(&[faces]).is_ok());
}

#[test]
fn touching_along_one_edge_is_not_an_intersection() {
    let mut faces = box_rings(["0", "0", "0"], ["1", "1", "1"]);
    faces.extend(box_rings(["1", "1", "0"], ["2", "2", "1"]));
    assert!(check(&[faces]).is_ok());
}

#[test]
fn overlapping_boxes_in_one_solid_intersect() {
    let mut faces = box_rings(["0", "0", "0"], ["2", "2", "2"]);
    faces.extend(box_rings(["1", "1", "1"], ["3", "3", "3"]));
    let r = check(&[faces]);
    assert!(crossing(r.clone()), "{r:?}");
}

#[test]
fn a_box_piercing_through_another_intersects() {
    // A thin post through the middle of a slab, passing all the way through.
    let mut faces = box_rings(["0", "0", "0"], ["4", "4", "1"]);
    faces.extend(box_rings(["1.5", "1.5", "-1"], ["2.5", "2.5", "2"]));
    assert!(crossing(check(&[faces])));
}

#[test]
fn an_edge_poking_into_a_face_interior_intersects() {
    // The second box's bottom face lies below the first's top face, side by side overlap in x only.
    let mut faces = box_rings(["0", "0", "0"], ["2", "2", "2"]);
    faces.extend(box_rings(["1", "0.5", "1.5"], ["3", "1.5", "3"]));
    assert!(crossing(check(&[faces])));
}

#[test]
fn coplanar_overlapping_faces_intersect_and_exact_ties_do_not() {
    // Two boxes standing on the same plane z = 0 with overlapping footprints.
    let mut faces = box_rings(["0", "0", "0"], ["2", "2", "1"]);
    faces.extend(box_rings(["1", "1", "0"], ["3", "3", "1"]));
    assert!(crossing(check(&[faces])));
    // Boxes sharing a whole face, back to back, touch and do not intersect.
    let mut faces = box_rings(["0", "0", "0"], ["2", "2", "1"]);
    faces.extend(box_rings(["2", "0", "0"], ["4", "2", "1"]));
    assert!(check(&[faces]).is_ok());
}

#[test]
fn separate_solids_are_never_compared() {
    let a = box_rings(["0", "0", "0"], ["2", "2", "2"]);
    let b = box_rings(["1", "1", "1"], ["3", "3", "3"]);
    assert!(check(&[a, b]).is_ok());
}

#[test]
fn the_import_entries_refuse_a_self_intersecting_source_with_the_faces_named() {
    let mut faces = box_rings(["0", "0", "0"], ["2", "2", "2"]);
    faces.extend(box_rings(["1", "1", "1"], ["3", "3", "3"]));
    let text = polyhedra_step(&[faces]);
    let step = StepId::new("imp").unwrap();
    let refusal = import_step_bodies(
        &mut Model::new(),
        text.as_bytes(),
        &sha256_hex(text.as_bytes()),
        &step,
    )
    .unwrap_err();
    let ImportRefusal::SelfIntersecting { solid, faces } = &refusal else {
        panic!("{refusal:?}");
    };
    assert!(solid.starts_with('#') && faces[0].starts_with('#') && faces[1].starts_with('#'));
    let applied = apply(
        &mut Model::new(),
        text.as_bytes(),
        "p.step",
        &ImportOptions::default(),
        &[],
        None,
    )
    .unwrap_err();
    assert_eq!(applied, refusal);
    let json = serde_json::to_string(&refusal).unwrap();
    assert!(json.starts_with(r#"{"code":"self_intersecting""#));
}

#[test]
fn repaired_and_clean_box_imports_still_pass() {
    let text = box_step(BoxFixture::new(["10", "10", "10"]));
    let step = StepId::new("imp").unwrap();
    assert!(
        import_step_bodies(
            &mut Model::new(),
            text.as_bytes(),
            &sha256_hex(text.as_bytes()),
            &step
        )
        .is_ok()
    );
}
