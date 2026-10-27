//! Minor-angle interval sets, retained nappes, finite trims and singular limits.
mod common;
use common::*;

fn supports(major: Real, minor: Real, offset: Real, base: Real, slope: Real) -> [Surface; 2] {
    [
        Surface::Cone(
            Cone::new(
                Point3::new(offset, rr(0), -base / slope.clone()),
                vi([0, 0, 1]).into(),
                slope,
                vi([1, 0, 0]).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ),
        Surface::Torus(
            Torus::new(
                pi([0, 0, 0]).into(),
                vi([0, 0, 1]).into(),
                major,
                minor,
                vi([1, 0, 0]).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ),
    ]
}
fn construct(
    surfaces: [Surface; 2],
    reversed: bool,
) -> (GeomStore, [SurfaceId; 2], Vec<Component>) {
    let mut geom = GeomStore::new();
    let mut ids = surfaces.map(|surface| geom.add_surface(surface));
    if reversed {
        ids.reverse();
    }
    let SurfacePair::Meets(parts) = exact(intersect_surfaces(
        &mut geom,
        ids[0],
        ids[1],
        &mut Decider::default(),
    )) else {
        panic!("expected retained components");
    };
    (geom, ids, parts)
}
fn overlap(a: &[Ball; 3], b: &[Ball; 3]) {
    for (a, b) in a.iter().zip(b) {
        assert!(a.lo() <= b.hi() && b.lo() <= a.hi(), "disconnected trace");
    }
}
fn check_trace(geom: &GeomStore, curve: CurveId) {
    assert!(curve_is_closed(geom, curve, &mut Decider::default()).unwrap());
    let traced = exact(trace(geom, curve, None, None, &mut Decider::default()));
    let domain = Ball::from_rational(&n(0), 128).hull(&Ball::from_rational(&n(1), 128));
    let mut first = None;
    let mut previous = None;
    let equations = geom.curve(curve).equations(geom).unwrap();
    for piece in traced.pieces() {
        let bounds = piece.taylor(&domain, 4, &Decider::default()).unwrap();
        let (lo, hi) = piece.range();
        let start = piece
            .evaluate(&Ball::from_rational(&lo, 128), &Decider::default())
            .unwrap();
        let end = piece
            .evaluate(&Ball::from_rational(&hi, 128), &Decider::default())
            .unwrap();
        if let Some(previous) = &previous {
            overlap(previous, &start);
        }
        if first.is_none() {
            first = Some(start);
        }
        previous = Some(end);
        for parameter in [lo.clone(), (&lo + &hi) / n(2), hi.clone()] {
            let position = piece
                .evaluate(&Ball::from_rational(&parameter, 128), &Decider::default())
                .unwrap();
            let exact_position = piece
                .evaluate_exact(&parameter, &Decider::default())
                .unwrap();
            for (coordinate, bound) in xyz(&exact_position).iter().zip(&bounds[0]) {
                for separation in [
                    coordinate - &Real::from(bound.lo()),
                    &Real::from(bound.hi()) - coordinate,
                ] {
                    assert!(
                        matches!(
                            Decider::default().sign(&separation),
                            Outcome::Proven(StrictSign::Positive) | Outcome::ProvenCoincident(_)
                        ),
                        "exact point outside whole-piece bounds at {parameter}"
                    );
                }
            }
            for equation in &equations {
                assert!(
                    equation
                        .eval_ball(&position, &Decider::default())
                        .unwrap()
                        .contains_zero()
                );
            }
        }
        let middle = (&lo + &hi) / n(2);
        let point = piece.evaluate_exact(&middle, &Decider::default()).unwrap();
        let Curve::Intersection(definition) = geom.curve(curve) else {
            panic!("branch")
        };
        let expected = geom
            .surface(definition.surfaces[0])
            .normal_at(&point, &mut Decider::default())
            .unwrap()
            .cross(
                &geom
                    .surface(definition.surfaces[1])
                    .normal_at(&point, &mut Decider::default())
                    .unwrap(),
            );
        let tangent = piece
            .derivative(&Ball::from_rational(&middle, 128), &Decider::default())
            .unwrap();
        let expected = [&expected.x, &expected.y, &expected.z]
            .map(|value| Decider::default().enclose_at(value, 128).unwrap());
        let dot = &(&expected[0] * &tangent[0])
            + &(&(&expected[1] * &tangent[1]) + &(&expected[2] * &tangent[2]));
        assert_eq!(
            dot.sign(),
            Some(StrictSign::Positive),
            "natural branch orientation"
        );
    }
    overlap(previous.as_ref().unwrap(), first.as_ref().unwrap());
}

#[test]
fn three_minor_angle_intervals_do_not_alias_components() {
    for reversed in [false, true] {
        let (mut geom, ids, parts) = construct(
            supports(
                r("101/100").into(),
                rr(1),
                rr(1),
                r("101/100").into(),
                rr(1),
            ),
            reversed,
        );
        assert_eq!(parts.len(), 3);
        let curves: Vec<_> = parts
            .into_iter()
            .map(|part| {
                let Component::Curve {
                    curve, singular, ..
                } = part
                else {
                    panic!("three regular loops")
                };
                assert!(singular.is_empty());
                check_trace(&geom, curve);
                curve
            })
            .collect();
        let mut selected_components = Vec::new();
        for (cosine, sine) in [
            ("3/5", "4/5"),
            ("-39999/40001", "-400/40001"),
            ("-400/40001", "-39999/40001"),
        ] {
            let rho = Real::from(r("101/100") + r(cosine));
            let radius = Real::from(r("101/100") + r(sine));
            let x = (rho.square() + rr(1) - radius.square()) / rr(2);
            let y = (rho.square() - x.square()).sqrt();
            for side in [1, -1] {
                let position = Point3::new(x.clone(), &y * &rr(side), r(sine).into());
                for surface in ids {
                    coincident(geom.side_of(surface, &position, &mut Decider::default()));
                }
                let point = geom.add_point(PointDef::Explicit(position));
                let mut members = 0;
                let mut selected = None;
                for &curve in &curves {
                    match on_curve(&geom, curve, point, &mut Decider::default()) {
                        Outcome::ProvenCoincident(_) => {
                            members += 1;
                            selected = Some(curve);
                        }
                        Outcome::Proven(Apart) => {}
                        other => panic!("component membership: {other:?}"),
                    }
                }
                assert_eq!(members, 1);
                if side > 0 {
                    selected_components.push(selected.unwrap());
                } else {
                    assert_eq!(selected, selected_components.last().copied());
                }
            }
        }
        let plane = geom.add_surface(axis_plane(2, n(0)));
        for (curve, expected_x) in
            selected_components
                .into_iter()
                .zip([Some("201/100"), Some("-1/100"), None])
        {
            let answer = exact(intersect_curve_surface(
                &mut geom,
                curve,
                plane,
                &mut Decider::default(),
            ));
            if let Some(expected_x) = expected_x {
                let CurveSurface::Meets(hits) = answer else {
                    panic!("one equatorial contact on this interval")
                };
                assert_eq!(hits.len(), 1);
                assert_eq!(hits[0].contact, Contact::Tangent);
                let point = geom.position(hits[0].point);
                real_eq(&point.x, &Real::from(r(expected_x)));
                real_eq(&point.y, &rr(0));
                real_eq(&point.z, &rr(0));
            } else {
                assert_eq!(answer, CurveSurface::Disjoint);
            }
        }
        for (i, &curve) in curves.iter().enumerate() {
            for &other in &curves[i + 1..] {
                assert_eq!(
                    same_curve(&geom, curve, other, &mut Decider::default()),
                    Outcome::Proven(Apart)
                );
            }
        }
    }
}

#[test]
fn housing_cone_trims_retain_both_contact_rails_and_plane_hits() {
    let pose = Pose::new(
        Rotation::from_quaternion(n(1), n(-1), n(0), n(0)).unwrap(),
        vi([0, 7, 0]).into(),
        &mut Decider::default(),
    )
    .unwrap();
    let (mut geom, ids, parts) = construct(
        supports(rr(14), rr(2), rr(17), rr(6), rr(1)).map(|surface| surface.transformed(&pose)),
        false,
    );
    assert_eq!(parts.len(), 2);
    let plane = geom.add_surface(axis_plane(1, n(5)));
    for part in parts {
        let Component::Curve {
            curve, singular, ..
        } = part
        else {
            panic!("two closed support loops")
        };
        assert!(singular.is_empty());
        check_trace(&geom, curve);
        let witness = geom.curve(curve).parameter_origin(&geom);
        let sign = match Decider::default().sign(&witness.z) {
            Outcome::Proven(StrictSign::Positive) => 1,
            Outcome::Proven(StrictSign::Negative) => -1,
            other => panic!("regular transverse side: {other:?}"),
        };
        let mouth = Point3::new(
            r("469/34").into(),
            rr(5),
            rr(sign) * root(&n(15)) * Real::from(r("21/34")),
        );
        let wall = Point3::new(
            r("397/34").into(),
            rr(7),
            rr(sign) * root(&n(8855)) / rr(34),
        );
        let a = closed(&mut geom, mouth.clone(), &ids);
        let b = closed(&mut geom, wall, &ids);
        let (from, to) = if sign < 0 { (a, b) } else { (b, a) };
        let selected = exact(trace(
            &geom,
            curve,
            Some(from),
            Some(to),
            &mut Decider::default(),
        ));
        for piece in selected.pieces() {
            let (lo, hi) = piece.range();
            let point = piece
                .evaluate(
                    &Ball::from_rational(&((&lo + &hi) / n(2)), 128),
                    &Decider::default(),
                )
                .unwrap();
            assert!(point[1].lo() > n(5) && point[1].hi() < n(7));
            let radial = &point[0].square() + &point[2].square();
            assert!(radial.lo() > n(144) && radial.hi() < n(196));
        }
        let CurveSurface::Meets(hits) = exact(intersect_curve_surface(
            &mut geom,
            curve,
            plane,
            &mut Decider::default(),
        )) else {
            panic!("mouth contact")
        };
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].contact, Contact::Tangent);
        for (actual, expected) in xyz(&geom.position(hits[0].point)).iter().zip(xyz(&mouth)) {
            real_eq(actual, &expected);
        }
    }
}

#[test]
fn cone_apex_isolated_contacts_nodes_and_cusp_limits_are_distinct() {
    let half_root_two = root(&r("1/2"));
    for (offset, base, slope, expected_points, expected_nodes) in [
        (rr(3), rr(1), rr(1), 1, 0),
        (rr(4), rr(0), rr(2), 0, 1),
        (rr(3) - half_root_two.clone(), half_root_two, rr(1), 0, 1),
    ] {
        let (geom, _, parts) = construct(supports(rr(3), rr(1), offset, base, slope), false);
        let mut points = 0;
        let mut curves = 0;
        for part in parts {
            match part {
                Component::Point { .. } => points += 1,
                Component::Curve {
                    curve, singular, ..
                } => {
                    curves += 1;
                    assert_eq!(singular.len(), expected_nodes);
                    check_trace(&geom, curve);
                    for node in singular {
                        let tangents =
                            tangent_at(&geom, curve, node, &mut Decider::default()).unwrap();
                        for tangent in tangents {
                            assert!(matches!(
                                Decider::default().sign(&tangent.norm_squared()),
                                Outcome::Proven(StrictSign::Positive)
                            ));
                        }
                    }
                }
            }
        }
        assert_eq!(points, expected_points);
        assert_eq!(curves, 1);
    }
}

#[test]
fn full_circle_sheets_and_shared_nodes_keep_their_domains() {
    let root_two = root(&n(2));
    for (offset, base, nodes) in [
        (rr(2), rr(3), 0),
        (rr(1) + root_two.clone(), rr(4), 1),
        (root_two, rr(3), 2),
    ] {
        let (geom, _, parts) = construct(supports(rr(3), rr(1), offset, base, rr(1)), false);
        assert_eq!(parts.len(), 2);
        let mut curves = Vec::new();
        for part in parts {
            let Component::Curve {
                curve, singular, ..
            } = part
            else {
                panic!("two closed loops")
            };
            assert_eq!(singular.len(), nodes);
            check_trace(&geom, curve);
            for node in singular {
                for tangent in tangent_at(&geom, curve, node, &mut Decider::default()).unwrap() {
                    assert!(matches!(
                        Decider::default().sign(&tangent.norm_squared()),
                        Outcome::Proven(StrictSign::Positive)
                    ));
                }
            }
            curves.push(curve);
        }
        assert_eq!(
            same_curve(&geom, curves[0], curves[1], &mut Decider::default()),
            Outcome::Proven(Apart)
        );
    }
}

#[test]
fn empty_nappe_and_isolated_external_tangency_are_not_closed_curves() {
    for surfaces in [
        supports(rr(3), rr(1), rr(8), rr(3), rr(1)),
        supports(rr(3), rr(1), rr(2), rr(-2), rr(1)),
    ] {
        let mut geom = GeomStore::new();
        let [a, b] = surfaces.map(|surface| geom.add_surface(surface));
        assert_eq!(
            exact(intersect_surfaces(&mut geom, a, b, &mut Decider::default())),
            SurfacePair::Disjoint
        );
    }
    let (geom, ids, parts) = construct(
        supports(rr(3), rr(1), rr(6) + root(&n(2)), rr(3), rr(1)),
        false,
    );
    let [Component::Point { point }] = parts.as_slice() else {
        panic!("one isolated tangency")
    };
    let position = geom.position(*point);
    for surface in ids {
        coincident(geom.side_of(surface, &position, &mut Decider::default()));
    }
    real_eq(&position.x, &(rr(3) + root(&r("1/2"))));
    real_eq(&position.y, &rr(0));
    real_eq(&position.z, &root(&r("1/2")));
}

#[test]
fn algebraic_placement_and_opposed_axes_keep_certified_consumers() {
    let major = root(&n(10));
    let pose = Pose::new(
        Rotation::from_quaternion(n(1), n(2), n(3), n(4)).unwrap(),
        vi([7, 11, -5]).into(),
        &mut Decider::default(),
    )
    .unwrap();
    let surfaces = [
        Surface::Cone(
            Cone::new(
                Point3::new(rr(2), rr(0), major.clone()),
                vi([0, 0, -3]).into(),
                rr(1),
                vi([1, 0, 0]).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ),
        Surface::Torus(
            Torus::new(
                pi([0, 0, 0]).into(),
                vi([0, 0, 5]).into(),
                major,
                rr(1),
                vi([1, 0, 0]).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ),
    ]
    .map(|surface| surface.transformed(&pose));
    let (mut geom, _, parts) = construct(surfaces, false);
    assert_eq!(parts.len(), 2);
    let plane = geom.add_surface(axis_plane(2, n(2)).transformed(&pose));
    for part in parts {
        let Component::Curve {
            curve, singular, ..
        } = part
        else {
            panic!("two regular algebraic loops")
        };
        assert!(singular.is_empty());
        check_trace(&geom, curve);
        assert_eq!(
            exact(intersect_curve_surface(
                &mut geom,
                curve,
                plane,
                &mut Decider::default()
            )),
            CurveSurface::Disjoint
        );
    }
}
