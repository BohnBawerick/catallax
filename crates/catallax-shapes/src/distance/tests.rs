use super::*;
use crate::{Cone, Cylinder, EllipticCylinder, Plane, RPoint3, RVector3, Sphere, Torus};
fn p(x: i64, y: i64, z: i64) -> Point3 {
    RPoint3::new(q(x), q(y), q(z)).into()
}
fn v(x: i64, y: i64, z: i64) -> Vector3 {
    RVector3::new(q(x), q(y), q(z)).into()
}
fn options(cells: usize) -> Options {
    Options {
        width: Rational::from_ratio(1, 1000).unwrap(),
        max_cells: cells,
    }
}
fn domain(surface: Surface, dec: &mut Decider) -> Domain {
    Domain::new(&surface, std::array::from_fn(|_| range(&q(-8), &q(8))), dec).unwrap()
}
fn all(_: &Point3, _: &[Ball; 3], _: &mut Decider) -> Result<Membership, Refusal> {
    Ok(Membership {
        witness: true,
        excluded: false,
    })
}
fn contains(bounds: &Bounds, value: &Real, dec: &mut Decider) {
    assert_ne!(
        math::sign(&(value - &Real::from(bounds.lower.clone())), dec).unwrap(),
        Sign::Negative,
        "{bounds:?}"
    );
    assert_ne!(
        math::sign(
            &(&Real::from(bounds.upper.clone().expect("finite upper")) - value),
            dec
        )
        .unwrap(),
        Sign::Negative,
        "{bounds:?}"
    );
}
#[test]
fn support_distances_cover_all_six_carriers_and_nonunit_frames() {
    let mut dec = Decider::default();
    let cases = [
        (
            Surface::Plane(Plane::new(p(0, 0, 0), v(0, 0, 2), v(3, 0, 0), &mut dec).unwrap()),
            p(1, 1, 1),
            Real::one(),
        ),
        (
            Surface::Cylinder(
                Cylinder::new(
                    p(0, 0, 0),
                    v(0, 0, 2),
                    Real::from_i64(2),
                    v(3, 0, 0),
                    &mut dec,
                )
                .unwrap(),
            ),
            p(3, 0, 1),
            Real::one(),
        ),
        (
            Surface::Cone(
                Cone::new(p(0, 0, 0), v(0, 0, 2), Real::one(), v(3, 0, 0), &mut dec).unwrap(),
            ),
            p(2, 0, 0),
            Real::from_i64(2).sqrt(),
        ),
        (
            Surface::Sphere(
                Sphere::new(
                    p(0, 0, 0),
                    Real::from_i64(2),
                    v(0, 0, 2),
                    v(3, 0, 0),
                    &mut dec,
                )
                .unwrap(),
            ),
            p(3, 0, 0),
            Real::one(),
        ),
        (
            Surface::Torus(
                Torus::new(
                    p(0, 0, 0),
                    v(0, 0, 2),
                    Real::from_i64(3),
                    Real::one(),
                    v(3, 0, 0),
                    &mut dec,
                )
                .unwrap(),
            ),
            p(5, 0, 0),
            Real::one(),
        ),
        (
            Surface::EllipticCylinder(
                EllipticCylinder::new(
                    p(0, 0, 0),
                    v(0, 0, 2),
                    v(3, 0, 0),
                    Real::from_i64(4),
                    Real::one(),
                    &mut dec,
                )
                .unwrap(),
            ),
            p(0, 2, 1),
            Real::one(),
        ),
    ];
    for (surface, point, want) in cases {
        let domain = domain(surface, &mut dec);
        let answer = point_distance(&domain, &point, &options(4096), &mut all, &mut dec).unwrap();
        contains(&answer, &want, &mut dec);
        assert!(answer.converged, "{answer:?}");
    }
}
#[test]
fn a_support_hit_inside_a_hole_is_not_a_zero_distance_witness() {
    let mut dec = Decider::default();
    let surface = Surface::Plane(Plane::new(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0), &mut dec).unwrap());
    let domain = Domain::new(
        &surface,
        [range(&q(-2), &q(2)), range(&q(-2), &q(2)), ball(&q(0))],
        &mut dec,
    )
    .unwrap();
    let mut trim = |p: &Point3, b: &[Ball; 3], d: &mut Decider| {
        let radius = p.x.square() + p.y.square();
        let witness = math::sign(&(&radius - &Real::one()), d)? != Sign::Negative
            && math::sign(&(radius - Real::from_i64(4)), d)? != Sign::Positive;
        let radius = b[0].square() + b[1].square();
        Ok(Membership {
            witness,
            excluded: radius.hi() < q(1) || radius.lo() > q(4),
        })
    };
    let coarse = point_distance(&domain, &p(0, 0, 0), &options(32), &mut trim, &mut dec).unwrap();
    let fine = point_distance(&domain, &p(0, 0, 0), &options(4096), &mut trim, &mut dec).unwrap();
    contains(&fine, &Real::one(), &mut dec);
    assert!(fine.lower > q(0), "{fine:?}");
    assert!(fine.lower >= coarse.lower);
    if let Some(upper) = coarse.upper {
        assert!(fine.upper.as_ref().unwrap() <= &upper);
    }
    assert!(fine.cells <= 4096);
    let repeat = point_distance(&domain, &p(0, 0, 0), &options(4096), &mut trim, &mut dec).unwrap();
    assert_eq!(
        (fine.lower, fine.upper, fine.cells),
        (repeat.lower, repeat.upper, repeat.cells)
    );
}
#[test]
fn face_pair_is_a_minimum_distance_and_empty_trims_donate_no_witness() {
    let mut dec = Decider::default();
    let a = Surface::Plane(Plane::new(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0), &mut dec).unwrap());
    let b = Surface::Plane(Plane::new(p(0, 0, 3), v(0, 0, 1), v(1, 0, 0), &mut dec).unwrap());
    let a = Domain::new(
        &a,
        [range(&q(-1), &q(1)), range(&q(-1), &q(1)), ball(&q(0))],
        &mut dec,
    )
    .unwrap();
    let b = Domain::new(
        &b,
        [range(&q(-2), &q(2)), range(&q(-2), &q(2)), ball(&q(3))],
        &mut dec,
    )
    .unwrap();
    let answer = face_distance(&a, &b, &options(64), &mut all, &mut all, &mut dec).unwrap();
    contains(&answer, &Real::from_i64(3), &mut dec);
    assert!(answer.converged);
    let mut absent = |_: &Point3, _: &[Ball; 3], _: &mut Decider| {
        Ok(Membership {
            witness: false,
            excluded: true,
        })
    };
    let answer = point_distance(&a, &p(0, 0, 1), &options(1), &mut absent, &mut dec).unwrap();
    assert!(answer.upper.is_none());
    assert!(!answer.converged);
    assert!(answer.cells <= 1);
}
#[test]
fn a_cone_query_behind_the_nappe_uses_the_apex_not_the_other_nappe() {
    let mut dec = Decider::default();
    let cone = domain(
        Surface::Cone(
            Cone::new(p(0, 0, 0), v(0, 0, 1), Real::one(), v(1, 0, 0), &mut dec).unwrap(),
        ),
        &mut dec,
    );
    let answer = point_distance(&cone, &p(0, 0, -2), &options(64), &mut all, &mut dec).unwrap();
    contains(&answer, &Real::from_i64(2), &mut dec);
    assert!(answer.converged);
}

#[test]
fn algebraic_frames_and_ellipse_axes_keep_the_euclidean_metric() {
    let mut dec = Decider::default();
    let cylinder = domain(
        Surface::Cylinder(
            Cylinder::new(
                p(0, 0, 0),
                v(0, 0, 2),
                Real::from_i64(2).sqrt(),
                v(1, 1, 0),
                &mut dec,
            )
            .unwrap(),
        ),
        &mut dec,
    );
    let answer =
        point_distance(&cylinder, &p(2, 2, 1), &options(4096), &mut all, &mut dec).unwrap();
    contains(&answer, &Real::from_i64(2).sqrt(), &mut dec);
    assert!(answer.converged);
    let ellipse = domain(
        Surface::EllipticCylinder(
            EllipticCylinder::new(
                p(0, 0, 0),
                v(0, 0, 2),
                v(1, 1, 0),
                Real::from_i64(3),
                Real::from_i64(2),
                &mut dec,
            )
            .unwrap(),
        ),
        &mut dec,
    );
    let answer = point_distance(&ellipse, &p(0, 0, 0), &options(4096), &mut all, &mut dec).unwrap();
    contains(&answer, &Real::from_i64(2).sqrt(), &mut dec);
    assert!(answer.converged);
}

#[test]
fn trim_classification_receives_the_full_cell_not_only_its_bound_intersection() {
    let mut dec = Decider::default();
    let sphere = Surface::Sphere(
        Sphere::new(
            p(0, 0, 0),
            Real::from_i64(2),
            v(0, 0, 1),
            v(1, 0, 0),
            &mut dec,
        )
        .unwrap(),
    );
    let cap = Domain::new(
        &sphere,
        [
            range(&q(-2), &q(2)),
            range(&q(-2), &q(2)),
            range(&q(1), &q(2)),
        ],
        &mut dec,
    )
    .unwrap();
    let mut trim = |point: &Point3, cover: &[Ball; 3], dec: &mut Decider| {
        if let Some(point) = point.as_rational() {
            for (coordinate, bound) in [&point.x, &point.y, &point.z].into_iter().zip(cover) {
                assert!(
                    bound.contains(coordinate),
                    "{bound:?} excludes its exact chart centre {coordinate}"
                );
            }
        }
        Ok(Membership {
            witness: math::sign(&(&point.z - &Real::one()), dec)? != Sign::Negative,
            excluded: cover[2].hi() < q(1),
        })
    };
    let answer = point_distance(&cap, &p(0, 0, -3), &options(64), &mut trim, &mut dec).unwrap();
    contains(&answer, &Real::from_i64(19).sqrt(), &mut dec);
}
mod finite {
    use super::*;
    use crate::Circle;
    fn rim(z: i64, radius_sq: i64, normal: Vector3) -> Circle {
        Circle {
            center: p(0, 0, z),
            normal,
            radius_sq: Real::from_i64(radius_sq),
            ref_dir: v(1, 0, 0),
        }
    }
    fn band() -> Surface {
        Surface::Cylinder(Cylinder {
            origin: p(0, 0, 0),
            axis: v(0, 0, 2),
            radius_sq: Real::from_i64(4),
            ref_dir: v(1, 0, 0),
        })
    }
    fn query(finite: &Finite, at: Point3, dec: &mut Decider) -> (Bounds, Branch) {
        finite
            .point_distance(&at, &options(1), dec)
            .unwrap()
            .unwrap()
    }
    #[test]
    fn band_branches_and_boundary_ties() {
        let mut dec = Decider::default();
        let (a, b) = (rim(0, 4, v(0, 0, 1)), rim(3, 4, v(0, 0, -1)));
        let finite = Finite::recognize(&band(), &[&b, &a], &mut dec)
            .unwrap()
            .unwrap();
        for (at, want, branch) in [
            (p(5, 0, 1), 3, Branch::Interior),
            (p(2, 0, 3), 0, Branch::Interior),
            (p(2, 0, 7), 4, Branch::Rim),
            (p(5, 0, 7), 5, Branch::Rim),
            (p(0, 0, 1), 2, Branch::Axis),
            (p(0, 0, 3), 2, Branch::Axis),
        ] {
            let (bounds, got) = query(&finite, at, &mut dec);
            assert_eq!(got, branch);
            assert!(bounds.lower <= q(want) && bounds.upper.unwrap() >= q(want));
            assert!(bounds.converged && bounds.cells == 0);
        }
    }
    #[test]
    fn cap_branches_and_rim_tie() {
        let mut dec = Decider::default();
        let plane = Surface::Plane(Plane {
            origin: p(0, 0, 3),
            normal: v(0, 0, 1),
            ref_dir: v(1, 0, 0),
        });
        let r = rim(3, 4, v(0, 0, 1));
        let finite = Finite::recognize(&plane, &[&r], &mut dec).unwrap().unwrap();
        for (at, want, branch) in [
            (p(1, 0, 7), 4, Branch::Interior),
            (p(2, 0, 4), 1, Branch::Interior),
            (p(5, 0, 7), 5, Branch::Rim),
        ] {
            let (bounds, got) = query(&finite, at, &mut dec);
            assert_eq!(got, branch);
            assert!(bounds.lower <= q(want) && bounds.upper.unwrap() >= q(want));
        }
    }
    #[test]
    fn other_trims_are_not_recognized() {
        let mut dec = Decider::default();
        let good = rim(0, 4, v(0, 0, 1));
        let other = rim(3, 4, v(0, 0, 1));
        let narrow = rim(3, 1, v(0, 0, 1));
        let tilted = rim(3, 4, v(0, 1, 1));
        let same = rim(0, 4, v(0, 0, 1));
        let off_axis = Circle {
            center: p(1, 0, 3),
            ..rim(3, 4, v(0, 0, 1))
        };
        for rims in [
            vec![&good],
            vec![&good, &narrow],
            vec![&good, &tilted],
            vec![&good, &same],
            vec![&good, &off_axis],
            vec![&good, &other, &other],
        ] {
            assert!(
                Finite::recognize(&band(), &rims, &mut dec)
                    .unwrap()
                    .is_none()
            );
        }
        let plane = Surface::Plane(Plane {
            origin: p(0, 0, 1),
            normal: v(0, 0, 1),
            ref_dir: v(1, 0, 0),
        });
        for rims in [vec![&good], vec![&tilted], vec![&good, &other]] {
            assert!(
                Finite::recognize(&plane, &rims, &mut dec)
                    .unwrap()
                    .is_none()
            );
        }
    }
}

mod spherical_finite {
    use super::*;
    use crate::Circle;
    fn sphere() -> Surface {
        Surface::Sphere(Sphere {
            center: p(0, 0, 0),
            radius: Real::from_i64(5),
            axis: v(2, 0, 0),
            ref_dir: v(0, 3, 0),
        })
    }
    fn rim(z: i64) -> Circle {
        Circle {
            center: p(0, 0, z),
            normal: v(0, 0, 2),
            radius_sq: Real::from_i64(25 - z * z),
            ref_dir: v(3, 0, 0),
        }
    }
    fn check(finite: &Finite, at: Point3, squared: i64, branch: Branch, dec: &mut Decider) {
        let (bounds, found) = finite
            .point_distance(&at, &options(1), dec)
            .unwrap()
            .unwrap();
        contains(&bounds, &Real::from_i64(squared).sqrt(), dec);
        assert_eq!(found, branch);
        assert!(bounds.cells == 0 && bounds.converged);
        let [_, w] = bounds.witnesses.unwrap();
        // Independent sphere equation and retained-domain checks, not a chart.
        assert_eq!(
            math::sign(
                &(w.minus(&p(0, 0, 0)).norm_squared() - Real::from_i64(25)),
                dec
            )
            .unwrap(),
            Sign::Zero
        );
    }
    #[test]
    fn branches_include_poles_centre_both_rims_and_radial_boundary_ties() {
        let mut dec = Decider::default();
        let (lo, hi) = (rim(-3), rim(3));
        let cap = Finite::recognize_oriented(&sphere(), &[(&hi, true)], &mut dec)
            .unwrap()
            .unwrap();
        for (at, squared, branch) in [
            (p(6, 0, 8), 25, Branch::Interior),
            (p(3, 0, 4), 0, Branch::Interior),
            (p(4, 0, 3), 0, Branch::Interior),
            (p(0, 0, 7), 4, Branch::Axis),
            (p(0, 0, -5), 80, Branch::Axis),
            (p(0, 0, 0), 25, Branch::Axis),
            (p(7, 0, -8), 130, Branch::Rim),
        ] {
            check(&cap, at, squared, branch, &mut dec);
        }
        let band = Finite::recognize_oriented(&sphere(), &[(&hi, false), (&lo, true)], &mut dec)
            .unwrap()
            .unwrap();
        for (at, squared, branch) in [
            (p(6, 0, 8), 29, Branch::Rim),
            (p(7, 0, -8), 34, Branch::Rim),
            (p(10, 0, 0), 25, Branch::Interior),
            (p(0, 0, 7), 32, Branch::Axis),
            (p(0, 0, -7), 32, Branch::Axis),
            (p(0, 0, 0), 25, Branch::Axis),
        ] {
            check(&band, at, squared, branch, &mut dec);
        }
    }
    #[test]
    fn hemisphere_major_cap_and_reversed_circle_normal_have_exact_retained_sides() {
        let mut dec = Decider::default();
        for (ring, forward, south_squared, north_squared) in [
            (rim(0), true, 50, 4),
            (rim(-3), true, 20, 4),
            (rim(3), false, 0, 32),
            (
                Circle {
                    normal: v(0, 0, -3),
                    ..rim(3)
                },
                false,
                80,
                4,
            ),
        ] {
            let finite = Finite::recognize_oriented(&sphere(), &[(&ring, forward)], &mut dec)
                .unwrap()
                .unwrap();
            check(&finite, p(0, 0, -5), south_squared, Branch::Axis, &mut dec);
            check(&finite, p(0, 0, 7), north_squared, Branch::Axis, &mut dec);
        }
    }
    #[test]
    fn irrational_radial_and_rim_witnesses_are_exact_and_admitted() {
        let mut dec = Decider::default();
        let ring = rim(3);
        let cap = Finite::recognize_oriented(&sphere(), &[(&ring, true)], &mut dec)
            .unwrap()
            .unwrap();
        for (at, expected, branch) in [
            (
                p(1, 1, 4),
                Real::from_i64(5) - Real::from_i64(18).sqrt(),
                Branch::Interior,
            ),
            (
                p(1, 1, -1),
                (Real::from_i64(34) - Real::from_i64(8) * Real::from_i64(2).sqrt()).sqrt(),
                Branch::Rim,
            ),
        ] {
            let (bounds, found) = cap
                .point_distance(&at, &options(1), &mut dec)
                .unwrap()
                .unwrap();
            contains(&bounds, &expected, &mut dec);
            assert_eq!(found, branch);
            let [_, w] = bounds.witnesses.unwrap();
            assert_ne!(
                math::sign(&(&w.z - &Real::from_i64(3)), &mut dec).unwrap(),
                Sign::Negative
            );
            assert_eq!(
                math::sign(
                    &(w.minus(&p(0, 0, 0)).norm_squared() - Real::from_i64(25)),
                    &mut dec
                )
                .unwrap(),
                Sign::Zero
            );
        }
    }
    #[test]
    fn missing_orientation_incompatible_bands_and_nonincident_or_polar_rings_decline() {
        let mut dec = Decider::default();
        let (a, b) = (rim(-3), rim(3));
        assert!(
            Finite::recognize(&sphere(), &[&b], &mut dec)
                .unwrap()
                .is_none()
        );
        let displaced = Circle {
            center: p(1, 0, 3),
            ..rim(3)
        };
        let wrong_radius = Circle {
            radius_sq: Real::one(),
            ..rim(3)
        };
        let skew = Circle {
            normal: v(0, 1, 1),
            ..rim(3)
        };
        let bad_reference = Circle {
            ref_dir: v(0, 0, 1),
            ..rim(3)
        };
        let polar = rim(5);
        for rings in [
            vec![],
            vec![(&a, true), (&b, true)],
            vec![(&a, false), (&b, true)],
            vec![(&a, true), (&a, false)],
            vec![(&a, true), (&b, false), (&b, false)],
            vec![(&displaced, true)],
            vec![(&wrong_radius, true)],
            vec![(&skew, true)],
            vec![(&bad_reference, true)],
            vec![(&polar, true)],
        ] {
            assert!(
                Finite::recognize_oriented(&sphere(), &rings, &mut dec)
                    .unwrap()
                    .is_none()
            );
        }
    }
    #[test]
    fn spherical_direct_route_preserves_validation_width_and_typed_budget_refusals() {
        use catallax_numbers::BudgetPolicy;
        let mut dec = Decider::default();
        let ring = rim(3);
        let cap = Finite::recognize_oriented(&sphere(), &[(&ring, true)], &mut dec)
            .unwrap()
            .unwrap();
        let at = p(0, 0, -5);
        assert_eq!(
            cap.point_distance(&at, &options(0), &mut dec)
                .unwrap_err()
                .code,
            "shapes.distance-domain"
        );
        let mut invalid = options(1);
        invalid.width = q(0);
        assert_eq!(
            cap.point_distance(&at, &invalid, &mut dec)
                .unwrap_err()
                .code,
            "shapes.distance-domain"
        );
        let mut tight = options(1);
        tight.width = q(1).mul_pow2(-512);
        assert!(cap.point_distance(&at, &tight, &mut dec).unwrap().is_none());
        let mut exhausted = Decider::new(BudgetPolicy {
            max_work: 0,
            ..BudgetPolicy::default()
        });
        assert_eq!(
            cap.point_distance(&at, &options(1), &mut exhausted)
                .unwrap_err()
                .code,
            "numbers.budget-exhausted"
        );
    }
}
