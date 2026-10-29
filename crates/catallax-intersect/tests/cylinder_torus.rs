//! Offset parallel axes: topology, branch identity, singular limits and finite traces.
mod common;
use common::*;

fn pair(
    major: &str,
    minor: &str,
    offset: &str,
    radius: &str,
    reversed: bool,
) -> (GeomStore, [SurfaceId; 2], SurfacePair) {
    let cylinder = cyl(p([offset, "0", "7"]), vi([0, 0, -2]), r(radius));
    let torus = Surface::Torus(
        Torus::new(
            pi([0, 0, 0]).into(),
            vi([0, 0, 3]).into(),
            r(major).into(),
            r(minor).into(),
            vi([1, 0, 0]).into(),
            &mut Decider::default(),
        )
        .unwrap(),
    );
    let mut geom = GeomStore::new();
    let mut supports = [geom.add_surface(cylinder), geom.add_surface(torus)];
    if reversed {
        supports.reverse();
    }
    let answer = exact(intersect_surfaces(
        &mut geom,
        supports[0],
        supports[1],
        &mut Decider::default(),
    ));
    (geom, supports, answer)
}

fn overlap(a: &[Ball; 3], b: &[Ball; 3]) {
    for (a, b) in a.iter().zip(b) {
        assert!(
            a.lo() <= b.hi() && b.lo() <= a.hi(),
            "disconnected trace: {a:?}, {b:?}"
        );
    }
}

fn check_trace(
    geom: &GeomStore,
    curve: CurveId,
    major: &str,
    minor: &str,
    offset: &str,
    reversed: bool,
) {
    let traced = exact(trace(geom, curve, None, None, &mut Decider::default()));
    let ball = |q: Rational| Ball::from_rational(&q, 128);
    let equations = geom.curve(curve).equations(geom).unwrap();
    let mut previous = None;
    let mut first = None;
    for piece in traced.pieces() {
        let (lo, hi) = piece.range();
        let start = piece
            .evaluate(&ball(lo.clone()), &Decider::default())
            .unwrap();
        let end = piece
            .evaluate(&ball(hi.clone()), &Decider::default())
            .unwrap();
        if let Some(previous) = &previous {
            overlap(previous, &start);
        }
        if first.is_none() {
            first = Some(start);
        }
        previous = Some(end);
        for parameter in [lo.clone(), (&lo + &hi) / n(2), hi.clone()] {
            let parameter = ball(parameter);
            let position = piece.evaluate(&parameter, &Decider::default()).unwrap();
            let derivative = piece.derivative(&parameter, &Decider::default()).unwrap();
            assert!(
                derivative
                    .iter()
                    .any(|coordinate| !coordinate.contains_zero()),
                "endpoint regularization lost the tangent"
            );
            for equation in &equations {
                assert!(
                    equation
                        .eval_ball(&position, &Decider::default())
                        .unwrap()
                        .contains_zero()
                );
            }
        }
        let parameter = ball((&lo + &hi) / n(2));
        let [x, y, z] = piece.evaluate(&parameter, &Decider::default()).unwrap();
        let [dx, dy, dz] = piece.derivative(&parameter, &Decider::default()).unwrap();
        // Independent implicit gradients, not the library's tangent helper.
        let major2 = ball(r(major) * r(major));
        let q = &((&x * &x) + (&y * &y) + (&z * &z)) + &major2 - ball(r(minor) * r(minor));
        let radial = &q - &(&ball(n(2)) * &major2);
        let d = ball(r(offset));
        let cross_x = &y * &(&z * &q);
        let cross_y = -(&(&x - &d) * &(&z * &q));
        let cross_z = -(&(&d * &y) * &radial);
        let forward = (dx * cross_x) + (dy * cross_y) + (dz * cross_z);
        if reversed {
            assert!(forward.hi() < n(0));
        } else {
            assert!(forward.lo() > n(0));
        }
    }
    overlap(previous.as_ref().unwrap(), first.as_ref().unwrap());
}

#[test]
fn radial_configurations_keep_components_nodes_and_orientation() {
    for (offset, radius, curves, nodes, points) in [
        ("6", "1", 0, 0, 0),
        ("1/2", "1/2", 0, 0, 0),
        ("5", "1", 0, 0, 1),
        ("1/2", "3/2", 0, 0, 1),
        ("5/2", "1", 1, 0, 0),
        ("4", "1", 1, 0, 0),
        ("3", "1/2", 2, 0, 0),
        ("3", "2", 2, 0, 0),
        ("7/2", "1/2", 2, 1, 0),
        ("5/2", "1/2", 2, 1, 0),
        ("3/2", "5/2", 2, 1, 0),
        ("9/2", "5/2", 2, 1, 0),
        ("3", "1", 2, 2, 0),
        ("1/2", "3", 2, 0, 0),
    ] {
        for reversed in [false, true] {
            let (mut geom, supports, answer) = pair("3", "1", offset, radius, reversed);
            if curves == 0 && points == 0 {
                assert_eq!(answer, SurfacePair::Disjoint);
                continue;
            }
            let SurfacePair::Meets(parts) = answer else {
                panic!("{offset}/{radius}: {answer:?}")
            };
            let mut found_curves = Vec::new();
            let mut found_points = 0;
            for part in parts {
                match part {
                    Component::Point { point } => {
                        found_points += 1;
                        for surface in supports {
                            coincident(geom.side(surface, point, &mut Decider::default()));
                        }
                    }
                    Component::Curve {
                        curve,
                        contact,
                        singular,
                    } => {
                        assert_eq!(contact, Contact::Crossing);
                        assert_eq!(singular.len(), nodes, "{offset}/{radius}");
                        assert!(curve_is_closed(&geom, curve, &mut Decider::default()).unwrap());
                        check_trace(&geom, curve, "3", "1", offset, reversed);
                        for node in singular {
                            coincident(on_curve(&geom, curve, node, &mut Decider::default()));
                            let [arriving, leaving] =
                                tangent_at(&geom, curve, node, &mut Decider::default()).unwrap();
                            for tangent in [arriving, leaving] {
                                assert!(matches!(
                                    Decider::default().sign(&tangent.norm_squared()),
                                    Outcome::Proven(StrictSign::Positive)
                                ));
                            }
                        }
                        found_curves.push(curve);
                    }
                }
            }
            assert_eq!(found_points, points, "{offset}/{radius}");
            assert_eq!(found_curves.len(), curves, "{offset}/{radius}");
            for &curve in &found_curves {
                for &other in &found_curves {
                    if curve == other {
                        continue;
                    }
                    let witness = geom.curve(other).parameter_origin(&geom);
                    let point = closed(&mut geom, witness, &supports);
                    assert_eq!(
                        on_curve(&geom, curve, point, &mut Decider::default()),
                        Outcome::Proven(Apart)
                    );
                    assert_eq!(
                        same_curve(&geom, curve, other, &mut Decider::default()),
                        Outcome::Proven(Apart)
                    );
                }
            }
        }
    }
}

#[test]
fn p5_drill_support_trims_keep_the_selected_half_loop() {
    let (mut geom, supports, answer) = pair("14", "2", "17", "2", false);
    let SurfacePair::Meets(parts) = answer else {
        panic!("P5 support must meet")
    };
    let [
        Component::Curve {
            curve, singular, ..
        },
    ] = parts.as_slice()
    else {
        panic!("one closed branch")
    };
    let curve = *curve;
    assert!(singular.is_empty());
    check_trace(&geom, curve, "14", "2", "17", false);
    let height = root(&n(3));
    let a = closed(
        &mut geom,
        Point3::new(rr(15), rr(0), height.clone()),
        &supports,
    );
    let b = closed(&mut geom, Point3::new(rr(15), rr(0), -height), &supports);
    for (from, to, sign) in [(a, b, 1), (b, a, -1)] {
        let traced = exact(trace(
            &geom,
            curve,
            Some(from),
            Some(to),
            &mut Decider::default(),
        ));
        let first = traced.pieces().first().unwrap();
        let last = traced.pieces().last().unwrap();
        for (piece, parameter, point) in
            [(first, first.range().0, from), (last, last.range().1, to)]
        {
            let actual = piece
                .evaluate(&Ball::from_rational(&parameter, 128), &Decider::default())
                .unwrap();
            let expected = xyz(&geom.position(point))
                .map(|x| Decider::default().enclose(&x, &r("1e-30")).unwrap());
            overlap(&actual, &expected);
        }
        for piece in traced.pieces() {
            let (lo, hi) = piece.range();
            let middle = Ball::from_rational(&((&lo + &hi) / n(2)), 128);
            let position = piece.evaluate(&middle, &Decider::default()).unwrap();
            if sign == 1 {
                assert!(position[1].lo() > n(0));
            } else {
                assert!(position[1].hi() < n(0));
            }
        }
    }
    let plane = geom.add_surface(axis_plane(2, n(0)));
    let CurveSurface::Meets(hits) = exact(intersect_curve_surface(
        &mut geom,
        curve,
        plane,
        &mut Decider::default(),
    )) else {
        panic!("equatorial plane must cross both halves");
    };
    assert_eq!(hits.len(), 2);
    let mut signs = Vec::new();
    for hit in hits {
        let position = geom.position(hit.point);
        real_eq(&position.x, &Real::from(r("541/34")));
        real_eq(&position.y.square(), &Real::from(r("3255/1156")));
        real_eq(&position.z, &rr(0));
        signs.push(match Decider::default().sign(&position.y) {
            Outcome::Proven(sign) => sign,
            other => panic!("expected a strict branch side: {other:?}"),
        });
    }
    assert!(signs.contains(&StrictSign::Positive) && signs.contains(&StrictSign::Negative));
}

#[test]
fn plane_exclusion_retains_nodes_and_tangent_contacts() {
    for (radius, plane_axis, expected_hits) in [("1/2", 2, 0), ("2", 1, 0), ("1", 2, 2)] {
        let (mut geom, _, answer) = pair("3", "1", "3", radius, false);
        let SurfacePair::Meets(parts) = answer else {
            panic!("two closed branches");
        };
        let plane = geom.add_surface(axis_plane(plane_axis, n(0)));
        for part in parts {
            let Component::Curve { curve, .. } = part else {
                panic!("expected a branch");
            };
            let answer = exact(intersect_curve_surface(
                &mut geom,
                curve,
                plane,
                &mut Decider::default(),
            ));
            if expected_hits == 0 {
                assert_eq!(answer, CurveSurface::Disjoint);
            } else {
                let CurveSurface::Meets(hits) = answer else {
                    panic!("the equatorial plane must retain both shared nodes");
                };
                assert_eq!(hits.len(), expected_hits);
                for hit in hits {
                    let point = geom.position(hit.point);
                    real_eq(&point.y, &rr(0));
                    real_eq(&point.z, &rr(0));
                    real_eq(&(&point.x - &rr(3)).square(), &rr(1));
                }
            }
        }
    }
    let (mut geom, _, answer) = pair("3", "1", "3", "1/2", false);
    let SurfacePair::Meets(parts) = answer else {
        panic!("two axial branches");
    };
    for part in parts {
        let Component::Curve { curve, .. } = part else {
            panic!("expected an axial branch");
        };
        let witness = geom.curve(curve).parameter_origin(&geom);
        let sign = match Decider::default().sign(&witness.z) {
            Outcome::Proven(StrictSign::Positive) => 1,
            Outcome::Proven(StrictSign::Negative) => -1,
            other => panic!("strict axial branch: {other:?}"),
        };
        for height in [n(sign), r("11/10") * n(sign)] {
            let plane = geom.add_surface(axis_plane(2, height.clone()));
            let answer = exact(intersect_curve_surface(
                &mut geom,
                curve,
                plane,
                &mut Decider::default(),
            ));
            if height.abs() > n(1) {
                assert_eq!(answer, CurveSurface::Disjoint);
                continue;
            }
            let CurveSurface::Meets(hits) = answer else {
                panic!("the top and bottom planes each touch twice");
            };
            assert_eq!(hits.len(), 2);
            for hit in hits {
                assert_eq!(hit.contact, Contact::Tangent);
                let point = geom.position(hit.point);
                real_eq(&point.x, &Real::from(r("71/24")));
                real_eq(&point.y.square(), &Real::from(r("143/576")));
                real_eq(&point.z, &Real::from(height.clone()));
            }
        }
    }
}

#[test]
fn algebraic_radii_and_rigid_placements_keep_certified_branch_bounds() {
    let major = root(&n(10));
    let cylinder = Surface::Cylinder(
        Cylinder::new(
            Point3::new(major.clone(), rr(0), rr(7)),
            vi([0, 0, -2]).into(),
            r("1/2").into(),
            vi([1, 0, 0]).into(),
            &mut Decider::default(),
        )
        .unwrap(),
    );
    let torus = Surface::Torus(
        Torus::new(
            pi([0, 0, 0]).into(),
            vi([0, 0, 3]).into(),
            major,
            rr(1),
            vi([1, 0, 0]).into(),
            &mut Decider::default(),
        )
        .unwrap(),
    );
    let pose = Pose::new(
        Rotation::from_quaternion(n(1), n(2), n(3), n(4)).unwrap(),
        vi([7, 11, -5]).into(),
        &mut Decider::default(),
    )
    .unwrap();
    let mut geom = GeomStore::new();
    let supports = [cylinder, torus].map(|surface| geom.add_surface(surface.transformed(&pose)));
    let SurfacePair::Meets(parts) = exact(intersect_surfaces(
        &mut geom,
        supports[0],
        supports[1],
        &mut Decider::default(),
    )) else {
        panic!("two displaced algebraic branches");
    };
    assert_eq!(parts.len(), 2);
    let middle_plane = geom.add_surface(axis_plane(2, n(0)).transformed(&pose));
    for part in parts {
        let Component::Curve {
            curve, singular, ..
        } = part
        else {
            panic!("regular closed branch");
        };
        assert!(singular.is_empty());
        let traced = exact(trace(&geom, curve, None, None, &mut Decider::default()));
        for piece in traced.pieces() {
            let exact_point = piece
                .evaluate_exact(&r("1/2"), &Decider::default())
                .unwrap();
            for surface in supports {
                coincident(geom.side_of(surface, &exact_point, &mut Decider::default()));
            }
            let domain = Ball::from_rational(&n(0), 128).hull(&Ball::from_rational(&n(1), 128));
            let bounds = piece.taylor(&domain, 4, &Decider::default()).unwrap();
            let at_middle = piece
                .taylor(&Ball::from_rational(&r("1/2"), 128), 4, &Decider::default())
                .unwrap();
            for (whole, middle) in bounds.iter().zip(&at_middle) {
                for (whole, middle) in whole.iter().zip(middle) {
                    assert!(whole.lo() <= middle.lo() && middle.hi() <= whole.hi());
                }
            }
        }
        assert_eq!(
            exact(intersect_curve_surface(
                &mut geom,
                curve,
                middle_plane,
                &mut Decider::default(),
            )),
            CurveSurface::Disjoint
        );
    }
}

#[test]
fn skew_axes_are_not_misclassified_as_parallel() {
    let mut geom = GeomStore::new();
    let cylinder = geom.add_surface(cyl(pi([3, 0, 0]), vi([1, 0, 1]), n(1)));
    let torus = geom.add_surface(Surface::Torus(
        Torus::new(
            pi([0, 0, 0]).into(),
            vi([0, 0, 1]).into(),
            rr(3),
            rr(1),
            vi([1, 0, 0]).into(),
            &mut Decider::default(),
        )
        .unwrap(),
    ));
    assert_eq!(
        intersect_surfaces(&mut geom, cylinder, torus, &mut Decider::default())
            .unwrap_err()
            .code,
        "intersect.traced-pair-required"
    );
}
