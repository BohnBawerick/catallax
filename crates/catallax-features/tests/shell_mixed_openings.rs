//! A tray with circular and polygonal through openings keeps rounded inner walls.
use catallax_features::evaluate;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Surface;

#[test]
fn top_open_shell_offsets_both_kinds_of_opening() {
    for source in [
        include_str!("fixtures/shell-mixed-openings.json"),
        include_str!("fixtures/shell-regular-openings.json"),
    ] {
        let recipe = Recipe::parse(source).unwrap();
        let (built, work) = with_operation_budget(20_000_000, || {
            evaluate(&recipe, &mut Decider::default()).unwrap()
        });
        eprintln!("{work:?}");
        assert!(!work.exhausted);
        assert!(built.model.check(built.result).is_empty());
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 33);
        assert_eq!(
            faces
                .iter()
                .filter(|face| matches!(
                    built.model.geom().surface(built.model.face(**face).surface),
                    Surface::Cylinder(_)
                ))
                .count(),
            8
        );
    }
}

#[test]
fn regular_opening_faces_and_mass_are_available_to_the_editor() {
    let recipe = Recipe::parse(include_str!("fixtures/shell-regular-openings.json")).unwrap();
    let (built, work) = with_operation_budget(20_000_000, || {
        evaluate(&recipe, &mut Decider::default()).unwrap()
    });
    assert!(!work.exhausted);
    let mut area = 0.0;
    let mut max_work = 0;
    let tolerance = Rational::parse("1e-6").unwrap();
    for face in built.model.faces_of_body(built.result) {
        let (answer, work) = with_operation_budget(100_000_000, || {
            catallax_measure::face_area(
                &built.model,
                face,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
        });
        assert!(!work.exhausted);
        assert!(answer.value.lo > Rational::zero());
        assert!(answer.value.width() < Rational::parse("1e-6").unwrap());
        area += ((answer.value.lo + answer.value.hi) / Rational::from_i64(2)).to_f64_nearest();
        max_work = max_work.max(work.used);
        let (center, work) = with_operation_budget(100_000_000, || {
            catallax_measure::face_centroid(
                &built.model,
                face,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
        });
        assert!(!work.exhausted);
        max_work = max_work.max(work.used);
        let (bounds, work) = with_operation_budget(100_000_000, || {
            catallax_measure::face_bounding_box(
                &built.model,
                face,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
        });
        assert!(!work.exhausted);
        max_work = max_work.max(work.used);
        for (i, coordinate) in center.value.iter().enumerate() {
            assert!(coordinate.lo >= &bounds.value.min[i] - &tolerance);
            assert!(coordinate.hi <= &bounds.value.max[i] + &tolerance);
        }
    }
    for edge in built.model.edges_of_body(built.result) {
        let (length, work) = with_operation_budget(100_000_000, || {
            catallax_measure::edge_length(
                &built.model,
                edge,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
        });
        assert!(!work.exhausted);
        assert!(length.value.lo > Rational::zero());
        max_work = max_work.max(work.used);
        let (middle, work) = with_operation_budget(100_000_000, || {
            catallax_measure::edge_midpoint(
                &built.model,
                edge,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
        });
        assert!(!work.exhausted);
        max_work = max_work.max(work.used);
        for (coordinate, bound) in middle.value.iter().zip([40, 30, 12]) {
            assert!(coordinate.lo >= -tolerance.clone());
            assert!(coordinate.hi <= Rational::from_i64(bound) + tolerance.clone());
        }
    }
    let pi = std::f64::consts::PI;
    let hex_area = 54.0 * 3.0_f64.sqrt();
    let expected_area = 2.0 * (1200.0 - hex_area - 9.0 * pi)
        + 1680.0
        + 432.0
        + 72.0 * pi
        + 2.0 * 11.0 * 66.0
        + 11.0 * (36.0 + 10.0 * pi);
    assert!(
        (area - expected_area).abs() < 1e-6,
        "{area} vs {expected_area}"
    );
    let (answer, work) = with_operation_budget(100_000_000, || {
        catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
    });
    assert!(!work.exhausted);
    let volume = answer.value.volume;
    let measured = ((volume.lo + volume.hi) / Rational::from_i64(2)).to_f64_nearest();
    let expected = 12.0 * (1200.0 - hex_area - 9.0 * pi)
        - 11.0 * (38.0 * 28.0 - (hex_area + 36.0 + pi) - 16.0 * pi);
    assert!(
        (measured - expected).abs() < 1e-6,
        "{measured} vs {expected}"
    );
    eprintln!("regular topology query max work={max_work}, mass work={work:?}");
}
