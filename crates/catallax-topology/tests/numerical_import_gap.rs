//! The gap rule for imports (decision Q9, and section 4.2 step 2 of the
//! Checkpoint 2 decisions): own builds store the same-parameter gap; an
//! import whose copies run at different speeds is reparametrised, and its
//! tolerance is not raised.
//!
//! The quarter cylinder is the fixture of the earlier parts. One circle use
//! of its wall gets a u/v copy that lies on the right parameter line and
//! runs at another speed (the stretch of known-answer edge 8).
#[path = "common/analytic.rs"]
mod fixture;
use catallax_numbers::{Decider, NumericalCheck, Rational, with_operation_budget};
use catallax_shapes::Surface;
use catallax_topology::numerical::{
    self as n, ImportGapOutcome, NumericalOptions,
    tolerance::{self, EntityKind},
};
use fixture::{bowed_copy, q, quarter, stretched_copy};

/// The quarter cylinder's curved face.
const WALL: usize = 4;

fn options(cap: &str) -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q(cap),
        vertex_tolerance_mm: None,
    }
}

fn with_radius(value: &str) -> n::NumericalDefinition {
    let mut source = quarter();
    let Surface::Cylinder(cylinder) = &mut source.faces[WALL].surface else {
        unreachable!()
    };
    cylinder.radius_sq = (&q(value) * &q(value)).into();
    source
}

fn edge_mm(body: &n::NumericalBody, edge: &str) -> Rational {
    body.tolerances()
        .of(EntityKind::Edge, edge)
        .expect("edge")
        .clone()
}

/// The stored record of the use of `edge` on the wall.
fn wall_record(body: &n::NumericalBody, edge: &str) -> Rational {
    let scope = format!("face \"f4\", edge {edge:?}");
    body.evidence()
        .checks()
        .iter()
        .find_map(|check| match check {
            NumericalCheck::Samples {
                name,
                scope: at,
                samples,
            } if name == "trim.edge-surface-distance" && at.starts_with(&scope) => {
                Some(samples[0].distance_upper_mm.clone())
            }
            _ => None,
        })
        .expect("a measured record")
}

/// "Own builds: same-parameter gap." The stretched copy runs up to 1e-3 of
/// the quarter arc ahead of the circle, 1.57e-3 mm at radius 1. An own build
/// stores rule T4's tolerance of that same-parameter gap and reads nothing
/// geometrically.
///
/// "Imports: geometric gap first; if the parametric gap is much larger,
/// reparametrise instead of raising the tolerance." The same definition as
/// an import: the geometric gap is zero, the copy is given the circle's
/// parameter, and the edge stores rule T1's floor.
#[test]
fn an_own_build_stores_the_same_parameter_gap_and_an_import_is_reparametrised() {
    let mut dec = Decider::default();
    let source = quarter();
    let (stretched, edge) = stretched_copy(&source, "0.001");
    let chord = 2.0 * (std::f64::consts::FRAC_PI_2 * 0.001 / 2.0).sin();

    let own = n::build_with_native_trims(
        source.clone(),
        vec![stretched.clone()],
        options("0.01"),
        &mut dec,
    )
    .unwrap();
    let own_mm = wall_record(&own, &edge).to_f64_nearest();
    assert!(
        own_mm >= 1.1 * chord && own_mm <= 1.1 * 1.1 * chord,
        "{own_mm}"
    );
    assert!(own.import_gaps().is_empty());

    let import =
        n::admit_import(source, vec![stretched.clone()], options("0.01"), &mut dec).unwrap();
    let [row] = import.import_gaps() else {
        panic!("{:?}", import.import_gaps())
    };
    println!("{row:?}");
    assert_eq!(row.edge, edge);
    assert_eq!(row.face, "f4");
    assert_eq!(row.outcome, ImportGapOutcome::Reparametrised);
    assert_eq!(row.parametric_mm, wall_record(&own, &edge));
    assert!(row.geometric_mm < q("1e-12"), "{}", row.geometric_mm);
    assert!(row.reparametrised_mm.as_ref().unwrap() <= &tolerance::floor_minimum_mm());
    // Admission continued with the same-parameter gap of the new copy.
    let stored = wall_record(&import, &edge);
    assert!(stored < q("1e-9"), "{stored}");
    assert_eq!(edge_mm(&import, &edge), tolerance::floor_minimum_mm());
    // The new copy is the retained u/v curve of the use, in place of the old.
    let [kept] = import.native_trims() else {
        panic!("one retained curve")
    };
    assert_eq!(
        (kept.edge.as_str(), kept.use_index),
        (edge.as_str(), stretched.use_index)
    );
    assert!(!kept.curve.same_definition(&stretched.curve));
}

/// "Reparametrise instead of raising the tolerance." A copy that runs up to
/// 2% of the arc ahead stores 0.0345 mm at the same parameter, over the
/// automatic 0.01 mm cap, and an own build refuses it. As an import it is
/// admitted at the automatic cap with nothing raised: no opt-in, and no
/// tolerance above the floor.
#[test]
fn a_difference_of_speed_over_the_cap_is_not_a_reason_to_raise_it() {
    let mut dec = Decider::default();
    let source = quarter();
    let (stretched, edge) = stretched_copy(&source, "0.02");
    assert_eq!(
        n::build_with_native_trims(
            source.clone(),
            vec![stretched.clone()],
            options("0.01"),
            &mut dec
        )
        .unwrap_err()
        .code,
        "numbers.numerical-tolerance-exceeded"
    );
    let import = n::admit_import(source, vec![stretched], options("0.01"), &mut dec).unwrap();
    assert_eq!(import.options().tolerance_mm, q("0.01"));
    let [row] = import.import_gaps() else {
        panic!("{:?}", import.import_gaps())
    };
    println!("{row:?}");
    assert_eq!(row.outcome, ImportGapOutcome::Reparametrised);
    assert!(row.parametric_mm > q("0.03"), "{}", row.parametric_mm);
    assert_eq!(
        import.tolerances().largest_mm(),
        tolerance::floor_minimum_mm()
    );
    assert_eq!(edge_mm(&import, &edge), tolerance::floor_minimum_mm());
}

/// A difference of speed on top of a gap in the geometry. The wall is 1e-4
/// mm wider than its circles, and one copy also runs up to 1e-3 of the arc
/// ahead. The geometric gap is the 1e-4 mm. The reparametrised copy stores
/// rule T4's tolerance of that, 1.1e-4 mm, where the copy as imported would
/// store 1.7e-3 mm.
#[test]
fn a_reparametrised_copy_keeps_the_gap_that_is_in_the_geometry() {
    let mut dec = Decider::default();
    let source = with_radius("1.0001");
    let (stretched, edge) = stretched_copy(&source, "0.001");
    let import = n::admit_import(source, vec![stretched], options("0.01"), &mut dec).unwrap();
    let rows = import.import_gaps();
    for row in rows {
        println!("{row:?}");
    }
    // Every use of the wall's two circles is off its wall and is read. Only
    // the stretched one is at another speed.
    let moved: Vec<_> = rows
        .iter()
        .filter(|row| row.outcome == ImportGapOutcome::Reparametrised)
        .collect();
    let [row] = moved.as_slice() else {
        panic!("one reparametrised use")
    };
    assert_eq!(row.edge, edge);
    assert!(rows.len() > 1);
    for other in rows.iter().filter(|other| other != row) {
        assert_eq!(other.outcome, ImportGapOutcome::Geometric);
        assert!(other.reparametrised_mm.is_none());
    }
    let gap = q("0.0001");
    assert!(row.geometric_mm >= gap && row.geometric_mm <= &gap * &q("1.1"));
    assert!(row.parametric_mm > q("0.0017"), "{}", row.parametric_mm);
    let stored = wall_record(&import, &edge);
    assert!(stored >= &gap * &q("1.1"), "{stored}");
    assert!(stored <= &gap * &q("1.21"), "{stored}");
}

/// A gap in the geometry is not a difference of speed. The wall 0.004 mm
/// wider than its circles is read geometrically, the two gaps agree, no copy
/// changes, and the body is the one admission built before the rule looked.
#[test]
fn a_gap_in_the_geometry_is_read_and_kept() {
    let mut dec = Decider::default();
    let built = n::build(with_radius("1.004"), options("0.01"), &mut dec).unwrap();
    let import =
        n::admit_import(with_radius("1.004"), Vec::new(), options("0.01"), &mut dec).unwrap();
    assert!(!import.import_gaps().is_empty());
    for row in import.import_gaps() {
        assert_eq!(row.outcome, ImportGapOutcome::Geometric, "{row:?}");
        assert!(row.geometric_mm >= q("0.004") && row.geometric_mm <= q("0.0044"));
        assert!(row.parametric_mm >= q("0.0044") && row.parametric_mm <= q("0.00484"));
    }
    assert!(import.native_trims().is_empty());
    assert_eq!(
        import.tolerances().largest_mm(),
        built.tolerances().largest_mm()
    );
}

/// A face copy bowed off its 3D curve is a gap in the geometry of the copy:
/// the rule reads it and leaves it, and rule T3's refit, not a change of
/// parameter, is what lowers it under the cap.
#[test]
fn a_bowed_copy_is_left_to_the_refit_of_rule_t3() {
    let mut dec = Decider::default();
    let source = quarter();
    let (bowed, edge) = bowed_copy(&source, "0.03");
    let import = n::admit_import(source, vec![bowed], options("0.01"), &mut dec).unwrap();
    let [row] = import.import_gaps() else {
        panic!("{:?}", import.import_gaps())
    };
    println!("{row:?}");
    assert_eq!(row.edge, edge);
    assert_eq!(row.outcome, ImportGapOutcome::Geometric);
    assert!(row.parametric_mm > q("0.01"));
    assert!(import.tolerances().largest_mm() <= q("0.01"));
}

/// An import whose copies agree is not read geometrically: no use is over
/// rule T1's floor. Its admission charges the work of the one build it
/// always made and, per edge use, one unit and one deadline poll for looking
/// at the use's stored record. The same request gives the same rows and the
/// same work on a second run.
#[test]
fn a_clean_import_is_not_read_geometrically_and_costs_two_units_per_use() {
    let admit = || {
        with_operation_budget(n::DEFAULT_WORK, || {
            n::admit_import(
                quarter(),
                Vec::new(),
                options("0.01"),
                &mut Decider::default(),
            )
        })
    };
    let (built, plain) = with_operation_budget(n::DEFAULT_WORK, || {
        n::build_with_native_trims(
            quarter(),
            Vec::new(),
            tolerance::import_options(q("0.01")),
            &mut Decider::default(),
        )
    });
    built.unwrap();
    let (import, work) = admit();
    let import = import.unwrap();
    assert!(import.import_gaps().is_empty());
    let uses: u64 = quarter()
        .faces
        .iter()
        .flat_map(|face| &face.loops)
        .map(|row| row.len() as u64)
        .sum();
    assert_eq!(work.used, plain.used + 2 * uses);

    let stretched = || {
        with_operation_budget(n::DEFAULT_WORK, || {
            let source = quarter();
            let (copy, _) = stretched_copy(&source, "0.001");
            n::admit_import(source, vec![copy], options("0.01"), &mut Decider::default())
        })
    };
    let (first, first_work) = stretched();
    let (second, second_work) = stretched();
    assert_eq!(first.unwrap().import_gaps(), second.unwrap().import_gaps());
    assert_eq!(first_work.used, second_work.used);
    println!(
        "clean import {} work, plain build {} work, stretched import {} work",
        work.used, plain.used, first_work.used
    );
}
