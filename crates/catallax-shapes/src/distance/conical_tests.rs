use super::*;
use crate::Circle;

fn cone() -> Cone {
    Cone {
        apex: p(0, 0, 0),
        axis: v(0, 0, 3),
        slope: Real::one(),
        ref_dir: v(2, 0, 0),
    }
}
fn rim(height: i64) -> Circle {
    Circle {
        center: p(0, 0, height),
        normal: v(0, 0, 2),
        radius_sq: Real::from_i64(height * height),
        ref_dir: v(3, 0, 0),
    }
}
#[test]
fn conical_meridian_branches_exclude_support_continuations_and_opposite_nappe() {
    let mut dec = Decider::default();
    let (lo, hi) = (rim(2), rim(6));
    for tip in [false, true] {
        let apex = p(0, 0, 0);
        let rings = if tip {
            vec![(&hi, false)]
        } else {
            vec![(&hi, false), (&lo, true)]
        };
        let finite = Finite::recognize_conical(&cone(), &rings, tip.then_some(&apex), &mut dec)
            .unwrap()
            .unwrap();
        for (at, squared, branch) in [
            (p(4, 0, 2), 2, Branch::Interior),
            (p(2, 0, 2), 0, Branch::Interior),
            (p(6, 0, 6), 0, Branch::Interior),
            (p(8, 0, 8), 8, Branch::Rim),
            (p(0, 0, 4), 8, Branch::Axis),
            (p(0, 0, -3), if tip { 9 } else { 29 }, Branch::Axis),
            (
                p(-2, 0, -2),
                if tip { 8 } else { 16 },
                if tip { Branch::Interior } else { Branch::Rim },
            ),
        ] {
            let (bounds, found) = finite
                .point_distance(&at, &options(1), &mut dec)
                .unwrap()
                .unwrap();
            contains(&bounds, &Real::from_i64(squared).sqrt(), &mut dec);
            assert_eq!(found, branch);
            assert!(bounds.converged && bounds.cells == 0);
            let [_, w] = bounds.witnesses.unwrap();
            assert_eq!(
                math::sign(&(&w.x.square() + &w.y.square() - w.z.square()), &mut dec).unwrap(),
                Sign::Zero
            );
            for value in [
                &w.z - &Real::from_i64(if tip { 0 } else { 2 }),
                &Real::from_i64(6) - &w.z,
            ] {
                assert_ne!(math::sign(&value, &mut dec).unwrap(), Sign::Negative);
            }
        }
    }
}
#[test]
fn conical_recognition_requires_exact_complete_positive_nappe_and_retained_sides() {
    let mut dec = Decider::default();
    let (lo, hi) = (rim(2), rim(6));
    let negative = rim(-2);
    let wrong_radius = Circle {
        radius_sq: Real::from_i64(3),
        ..lo.clone()
    };
    let off_axis = Circle {
        center: p(1, 0, 2),
        ..lo.clone()
    };
    let oblique = Circle {
        normal: v(0, 1, 1),
        ..lo.clone()
    };
    let bad_reference = Circle {
        ref_dir: v(0, 0, 1),
        ..lo.clone()
    };
    for rings in [
        vec![],
        vec![(&hi, false)],
        vec![(&lo, true), (&hi, true)],
        vec![(&lo, false), (&hi, true)],
        vec![(&lo, true), (&lo, false)],
        vec![(&negative, true), (&hi, false)],
        vec![(&wrong_radius, true), (&hi, false)],
        vec![(&off_axis, true), (&hi, false)],
        vec![(&oblique, true), (&hi, false)],
        vec![(&bad_reference, true), (&hi, false)],
        vec![(&lo, true), (&hi, false), (&hi, false)],
    ] {
        assert!(
            Finite::recognize_conical(&cone(), &rings, None, &mut dec)
                .unwrap()
                .is_none()
        );
    }
    assert!(
        Finite::recognize(&Surface::Cone(cone()), &[&lo, &hi], &mut dec)
            .unwrap()
            .is_none()
    );
    for (rings, apex) in [
        (vec![(&hi, true)], p(0, 0, 0)),
        (vec![(&hi, false)], p(0, 0, 1)),
        (vec![(&lo, true), (&hi, false)], p(0, 0, 0)),
    ] {
        assert!(
            Finite::recognize_conical(&cone(), &rings, Some(&apex), &mut dec)
                .unwrap()
                .is_none()
        );
    }
    let flipped = Circle {
        normal: v(0, 0, -5),
        ..lo
    };
    assert!(
        Finite::recognize_conical(&cone(), &[(&flipped, false), (&hi, false)], None, &mut dec)
            .unwrap()
            .is_some()
    );
}
#[test]
fn conical_direct_controls_remain_width_cell_and_work_limited() {
    use catallax_numbers::BudgetPolicy;
    let mut dec = Decider::default();
    let (lo, hi) = (rim(2), rim(6));
    let finite = Finite::recognize_conical(&cone(), &[(&lo, true), (&hi, false)], None, &mut dec)
        .unwrap()
        .unwrap();
    for bad in [
        Options {
            width: q(0),
            max_cells: 1,
        },
        Options {
            width: q(1),
            max_cells: 0,
        },
    ] {
        assert_eq!(
            finite
                .point_distance(&p(4, 0, 2), &bad, &mut dec)
                .unwrap_err()
                .code,
            "shapes.distance-domain"
        );
    }
    let tiny = Options {
        width: q(1).mul_pow2(-512),
        max_cells: 1,
    };
    assert!(
        finite
            .point_distance(&p(4, 0, 2), &tiny, &mut dec)
            .unwrap()
            .is_none()
    );
    let mut exhausted = Decider::new(BudgetPolicy {
        max_work: 0,
        ..BudgetPolicy::default()
    });
    assert_eq!(
        finite
            .point_distance(&p(4, 0, 2), &options(1), &mut exhausted)
            .unwrap_err()
            .code,
        "numbers.budget-exhausted"
    );
}

#[test]
fn conical_algebraic_slope_retains_exact_incidence_and_axis_projection() {
    let mut dec = Decider::default();
    let cone = Cone {
        slope: Real::from_i64(2).sqrt(),
        ..cone()
    };
    let lo = Circle {
        radius_sq: Real::from_i64(8),
        ..rim(2)
    };
    let hi = Circle {
        radius_sq: Real::from_i64(72),
        ..rim(6)
    };
    let finite = Finite::recognize_conical(&cone, &[(&lo, true), (&hi, false)], None, &mut dec)
        .unwrap()
        .unwrap();
    let (answer, branch) = finite
        .point_distance(&p(0, 0, 9), &options(1), &mut dec)
        .unwrap()
        .unwrap();
    contains(&answer, &Real::from_i64(54).sqrt(), &mut dec);
    assert_eq!(branch, Branch::Axis);
    let [_, witness] = answer.witnesses.unwrap();
    assert_eq!(
        math::sign(&(&witness.z - &Real::from_i64(3)), &mut dec).unwrap(),
        Sign::Zero
    );
    assert_eq!(
        math::sign(
            &(witness.x.square() + witness.y.square() - Real::from_i64(18)),
            &mut dec
        )
        .unwrap(),
        Sign::Zero
    );
}
