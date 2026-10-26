//! The real application entry consumes the exact open-cap shell and its edits.

use std::{collections::BTreeMap, time::Instant};

use catallax_features::{evaluate, evaluate_incremental};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Operation, Recipe, Role};

fn recipe(thickness: &str) -> Recipe {
    let mut recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/open-spherical-cap/base-open.json"
    ))
    .unwrap();
    let Operation::Shell { thickness: t, .. } = &mut recipe.steps[3].op else {
        unreachable!()
    };
    *t = Rational::parse(thickness).unwrap().into();
    recipe
}

fn closed_display(mesh: &Mesh) {
    assert!(!mesh.triangles.is_empty());
    assert!(mesh.unavailable_controls.is_empty());
    let mut uses = BTreeMap::<_, (usize, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let count = uses.entry((a.min(b), a.max(b))).or_default();
            count.0 += 1;
            count.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(uses.values().all(|&count| count == (2, 0)));
    for face in &mesh.faces {
        assert!(face.max_chord_error <= MeshOptions::default().chord_tolerance);
    }
}

#[test]
fn application_and_suffix_edit_publish_closed_named_display_with_default_limits() {
    let before_recipe = recipe("1/2");
    let start = Instant::now();
    let before =
        catallax_app::build(&before_recipe, &Default::default(), &mut Decider::default()).unwrap();
    eprintln!(
        "application build plus requested display: {:?}",
        start.elapsed()
    );
    assert!(before.display.marks.is_empty());
    closed_display(&before.display.value);
    assert_eq!(before.names.faces.len(), 5);
    assert!(before.names.edges.values().all(Option::is_some));
    assert_eq!(
        before
            .names
            .faces
            .values()
            .flatten()
            .filter(|o| matches!(o.role, Role::Inner { .. }))
            .count(),
        2
    );

    let edited_recipe = recipe("3/4");
    let mut prefix_recipe = before_recipe.clone();
    prefix_recipe.steps.truncate(3);
    prefix_recipe.result = prefix_recipe.steps[2].id.clone();
    let prefix = evaluate_incremental(&prefix_recipe, None, &mut Decider::default()).unwrap();
    let start = Instant::now();
    let edited = catallax_app::build_incremental(
        &edited_recipe,
        Some(prefix),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    eprintln!(
        "suffix shell edit plus requested display: {:?}",
        start.elapsed()
    );
    assert!(edited.display.marks.is_empty());
    closed_display(&edited.display.value);
    assert_eq!(
        before.names.faces.values().collect::<Vec<_>>(),
        edited.names.faces.values().collect::<Vec<_>>()
    );
    let fresh =
        catallax_app::build(&edited_recipe, &Default::default(), &mut Decider::default()).unwrap();
    assert_eq!(fresh.display.value, edited.display.value);
    assert_eq!(fresh.display.marks, edited.display.marks);
    assert_eq!(fresh.names, edited.names);

    let start = Instant::now();
    let mass = before
        .mass_properties(&Default::default(), &mut Decider::default())
        .unwrap();
    eprintln!("requested properties: {:?}", start.elapsed());
    assert!(mass.marks.is_empty());
    assert!(mass.value.volume.lo > Rational::from_i64(0));
}

#[test]
fn construction_display_and_output_are_accounted_at_the_application_allowance() {
    let (result, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        let start = Instant::now();
        let built = evaluate(&recipe("1/2"), &mut Decider::default()).unwrap();
        eprintln!("construction: {:?}", start.elapsed());
        let start = Instant::now();
        let display = mesh_body(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        eprintln!("requested B-rep display: {:?}", start.elapsed());
        closed_display(&display.value);
        let start = Instant::now();
        let mut bytes = Vec::new();
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        eprintln!("requested STEP output: {:?}", start.elapsed());
        assert!(bytes.starts_with(b"ISO-10303-21;"));
        bytes
    });
    assert!(!work.exhausted, "{work:?}");
    assert!(!result.is_empty());
    eprintln!(
        "construction/display/STEP work at {} allowance: {work:?}",
        catallax_app::DEFAULT_BUILD_WORK
    );
}
