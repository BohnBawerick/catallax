//! The application preserves exact signed offsets, selectable names and properties.
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{Operation, Recipe, Role};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&r(n), 192)
}
fn recipe(distance: i64) -> Recipe {
    let mut recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/algebraic-extrusions/rectangle.json"
    ))
    .unwrap();
    let Operation::Extrude { distance: d, .. } = &mut recipe.steps[1].op else {
        unreachable!()
    };
    *d = distance.into();
    recipe
}
fn contains(value: &catallax_measure::Enclosure, expected: Ball) {
    assert!(
        value.lo <= expected.lo() && value.hi >= expected.hi(),
        "{value:?} versus {expected:?}"
    );
    assert!(value.width() <= Rational::parse("1e-9").unwrap());
}
#[test]
fn signed_oblique_display_and_requested_mass_match_independent_affine_integrals() {
    for distance in [5, -5] {
        let app = catallax_app::build(
            &recipe(distance),
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(app.display.marks.is_empty());
        assert!(!app.display.value.triangles.is_empty());
        for role in [Role::Start, Role::End] {
            assert!(
                app.display
                    .value
                    .triangles
                    .iter()
                    .any(|t| app.names.faces[&t.face].iter().any(|o| o.role == role))
            );
        }
        assert!(app.names.edges.values().all(Option::is_some));
        assert!(app.names.vertices.values().all(Option::is_some));
        let mass = app
            .mass_properties(&Default::default(), &mut Decider::default())
            .unwrap();
        assert!(mass.marks.is_empty());
        let root = b(2).sqrt().unwrap();
        let a = b(distance) * root.recip().unwrap();
        let volume = b(40) * root.recip().unwrap();
        contains(&mass.value.volume, volume.clone());
        contains(&mass.value.area, b(36) + b(20) * root);
        let centroid = mass.value.centroid.unwrap();
        for (value, expected) in centroid.iter().zip([
            b(2) + &a * &b(2).recip().unwrap(),
            b(1),
            &a * &b(2).recip().unwrap(),
        ]) {
            contains(value, expected);
        }
        // Uniform rectangle plus independent uniform sweep parameter:
        // Cov = diag(16/12,4/12,0) + vv^T/12, I=V*(tr(Cov)Id-Cov).
        let inertia = app
            .inertia(&Default::default(), &mut Decider::default())
            .unwrap()
            .value
            .unwrap();
        for (i, row) in inertia.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let factor = match (i, j) {
                    (0, 0) => (11, 8),
                    (1, 1) => (41, 12),
                    (2, 2) => (65, 24),
                    (0, 2) | (2, 0) => (-25, 24),
                    _ => (0, 1),
                };
                contains(
                    value,
                    &volume
                        * &Ball::from_rational(
                            &Rational::from_ratio(factor.0, factor.1).unwrap(),
                            192,
                        ),
                );
            }
        }
    }
}
#[test]
fn distance_edit_retains_all_face_and_edge_origins() {
    let before = catallax_app::build_incremental(
        &recipe(5),
        None,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let face_origins: Vec<_> = before.names.faces.values().cloned().collect();
    let edge_names: Vec<_> = before.names.edges.values().cloned().collect();
    let after = catallax_app::build_incremental(
        &recipe(-7),
        Some(before.checkpoint),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        face_origins,
        after.names.faces.values().cloned().collect::<Vec<_>>()
    );
    assert_eq!(
        edge_names,
        after.names.edges.values().cloned().collect::<Vec<_>>()
    );
    assert!(!after.display.value.triangles.is_empty());
    assert!(after.display.marks.is_empty());
}
