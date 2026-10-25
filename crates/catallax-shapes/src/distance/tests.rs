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
