//! Conic sections, nappe/component selection and finite open-curve traces.
mod common;
use common::*;

fn cone(apex: RPoint3, axis: RVector3, slope: Rational) -> Surface {
    let reference = perp(&axis);
    Surface::Cone(
        Cone::new(
            apex.into(),
            axis.into(),
            slope.into(),
            reference.into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
fn torus(center: RPoint3, major: Rational, minor: Rational) -> Surface {
    Surface::Torus(
        Torus::new(
            center.into(),
            vi([0, 0, 1]).into(),
            major.into(),
            minor.into(),
            vi([1, 0, 0]).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
fn section(a: Surface, b: Surface) -> (GeomStore, [SurfaceId; 2], SurfacePair) {
    let mut geom = GeomStore::new();
    let supports = [geom.add_surface(a), geom.add_surface(b)];
    let answer = exact(intersect_surfaces(
        &mut geom,
        supports[0],
        supports[1],
        &mut Decider::default(),
    ));
    (geom, supports, answer)
}
fn circles(a: Surface, b: Surface, expected: &[(Point3, Real, Contact)], points: &[Point3]) {
    let (mut geom, supports, answer) = section(a, b);
    let SurfacePair::Meets(parts) = answer else {
        panic!("expected circles: {answer:?}")
    };
    let mut curves = vec![];
    let mut isolated = vec![];
    for component in parts {
        match component {
            Component::Curve {
                curve,
                contact,
                singular,
            } => {
                assert!(singular.is_empty());
                curves.push((curve, contact));
            }
            Component::Point { point } => isolated.push(point),
        }
    }
    assert_eq!(isolated.len(), points.len());
    for (&point, expected) in isolated.iter().zip(points) {
        point_eq(&geom.position(point), expected);
    }
    assert_eq!(curves.len(), expected.len());
    for (&(curve, contact), (center, radius_sq, expected_contact)) in curves.iter().zip(expected) {
        assert_eq!(contact, *expected_contact);
        let form = curve_form(&geom, curve, &mut Decider::default()).unwrap();
        let CurveForm::Circle {
            center: actual,
            radius_sq: radius,
            ..
        } = &form
        else {
            panic!("expected circle, got {form:?}")
        };
        point_eq(actual, center);
        real_eq(radius, radius_sq);
        assert!(curve_is_closed(&geom, curve, &mut Decider::default()).unwrap());
        let point = exact(arc_point(&geom, curve, None, None, &mut Decider::default()));
        for support in supports {
            coincident(geom.side_of(support, &point, &mut Decider::default()));
        }
        let point = closed(&mut geom, point, &supports);
        coincident(on_curve(&geom, curve, point, &mut Decider::default()));
        if contact == Contact::Crossing {
            let cross = geom
                .normal(supports[0], point)
                .cross(&geom.normal(supports[1], point));
            positive_multiple(
                &tangent_at(&geom, curve, point, &mut Decider::default()).unwrap()[1],
                &cross,
            );
        }
        let trace = exact(trace(&geom, curve, None, None, &mut Decider::default()));
        let equations = supports.map(|surface| geom.implicit(surface));
        for piece in trace.pieces() {
            for t in [r("0"), r("1/3"), r("1")] {
                let position = piece
                    .evaluate(&Ball::from_rational(&t, 128), &Decider::default())
                    .unwrap();
                for equation in &equations {
                    assert!(
                        equation
                            .eval_ball(&position, &catallax_numbers::Decider::default())
                            .unwrap()
                            .contains_zero()
                    );
                }
            }
        }
    }
    for (i, &(curve, _)) in curves.iter().enumerate() {
        for (j, &(other, _)) in curves.iter().enumerate() {
            if i == j {
                continue;
            }
            let point = geom.curve(other).parameter_origin(&geom);
            let point = closed(&mut geom, point, &supports);
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

#[test]
fn plane_cone_closed_sections_exclude_the_other_nappe() {
    let c = cone(pi([0, 0, 0]), vi([0, 0, 1]), n(1));
    circles(
        axis_plane(2, n(2)),
        c.clone(),
        &[(pi([0, 0, 2]).to_real(), rr(4), Contact::Crossing)],
        &[],
    );
    let (_, _, result) = section(axis_plane(2, n(-2)), c.clone());
    assert_eq!(result, SurfacePair::Disjoint);
    let (g, _, result) = section(axis_plane(2, n(0)), c.clone());
    let SurfacePair::Meets(parts) = result else {
        panic!("expected apex")
    };
    let [Component::Point { point }] = parts.as_slice() else {
        panic!("expected isolated apex")
    };
    point_eq(&g.position(*point), &pi([0, 0, 0]).to_real());
    let (g, _, result) = section(plane(pi([0, 0, 3]), vi([1, 0, 2])), c);
    let SurfacePair::Meets(parts) = result else {
        panic!("expected ellipse")
    };
    let [Component::Curve { curve, .. }] = parts.as_slice() else {
        panic!("expected one ellipse")
    };
    let CurveForm::Ellipse {
        center,
        major_sq,
        minor_sq,
        ..
    } = curve_form(&g, *curve, &mut Decider::default()).unwrap()
    else {
        panic!("expected ellipse")
    };
    point_eq(&center, &pi([-2, 0, 4]).to_real());
    real_eq(&major_sq, &rr(20));
    real_eq(&minor_sq, &rr(12));
}

#[test]
fn open_conics_have_ordered_finite_traces_and_correct_contacts() {
    for (plane, endpoints, vertex, cap, second) in [
        (
            plane(pi([-1, 0, 0]), vi([1, 0, -1])),
            [p(["0", "-1", "1"]), p(["0", "1", "1"])],
            p(["-1/2", "0", "1/2"]),
            r("1"),
            r("2"),
        ),
        (
            axis_plane(0, n(1)),
            [p(["1", "-3/4", "5/4"]), p(["1", "3/4", "5/4"])],
            pi([1, 0, 1]),
            r("5/4"),
            r("9/8"),
        ),
    ] {
        let (mut g, supports, result) = section(plane, cone(pi([0, 0, 0]), vi([0, 0, 1]), n(1)));
        let SurfacePair::Meets(parts) = result else {
            panic!("expected open conic")
        };
        let [Component::Curve { curve, .. }] = parts.as_slice() else {
            panic!("expected one component")
        };
        let curve = *curve;
        assert!(!curve_is_closed(&g, curve, &mut Decider::default()).unwrap());
        let [a, b] = endpoints.map(|p| explicit(&mut g, p));
        coincident(on_curve(&g, curve, a, &mut Decider::default()));
        coincident(on_curve(&g, curve, b, &mut Decider::default()));
        assert_eq!(
            order_on_curve(&g, curve, a, b, &mut Decider::default()),
            Outcome::Proven(StrictOrder::Less)
        );
        point_eq(
            &exact(arc_point(
                &g,
                curve,
                Some(a),
                Some(b),
                &mut Decider::default(),
            )),
            &vertex.to_real(),
        );
        let trace = exact(trace(&g, curve, Some(a), Some(b), &mut Decider::default()));
        assert_eq!(trace.pieces().len(), 1);
        let jets = trace.pieces()[0]
            .taylor(&Ball::from_rational(&r("1/2"), 128), 4, &Decider::default())
            .unwrap();
        assert!(jets[0][2].contains(&vertex.z));
        assert!(jets[2][2].contains(&second));
        let equations = supports.map(|s| g.implicit(s));
        for t in [r("0"), r("1/4"), r("1/2"), r("3/4"), r("1")] {
            let position = trace.pieces()[0]
                .evaluate(&Ball::from_rational(&t, 128), &Decider::default())
                .unwrap();
            for equation in &equations {
                assert!(
                    equation
                        .eval_ball(&position, &catallax_numbers::Decider::default())
                        .unwrap()
                        .contains_zero()
                );
            }
        }
        let cutting = g.add_surface(axis_plane(2, cap));
        let CurveSurface::Meets(hits) = exact(intersect_curve_surface(
            &mut g,
            curve,
            cutting,
            &mut Decider::default(),
        )) else {
            panic!("expected two cap hits")
        };
        assert_eq!(hits.len(), 2);
        for (hit, expected) in hits.iter().zip([a, b]) {
            point_eq(&g.position(hit.point), &g.position(expected));
            assert_eq!(hit.contact, Contact::Crossing);
        }
        let touching = g.add_surface(axis_plane(2, vertex.z.clone()));
        let CurveSurface::Meets(hits) = exact(intersect_curve_surface(
            &mut g,
            curve,
            touching,
            &mut Decider::default(),
        )) else {
            panic!("expected vertex tangency")
        };
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].contact, Contact::Tangent);
        point_eq(&g.position(hits[0].point), &vertex.to_real());
        if matches!(
            curve_form(&g, curve, &mut Decider::default()).unwrap(),
            CurveForm::Hyperbola { .. }
        ) {
            let other_nappe = explicit(&mut g, pi([1, 0, -1]));
            assert_eq!(
                on_curve(&g, curve, other_nappe, &mut Decider::default()),
                Outcome::Proven(Apart)
            );
        }
    }
}

#[test]
fn generator_components_are_rays_not_whole_lines() {
    let (mut g, supports, result) = section(
        axis_plane(0, n(0)),
        cone(pi([0, 0, 0]), vi([0, 0, 1]), n(1)),
    );
    let SurfacePair::Meets(parts) = result else {
        panic!("expected generators")
    };
    assert_eq!(parts.len(), 2);
    for part in parts {
        let Component::Curve {
            curve, singular, ..
        } = part
        else {
            panic!("expected ray")
        };
        assert_eq!(singular.len(), 1);
        let CurveForm::Line { dir, .. } = curve_form(&g, curve, &mut Decider::default()).unwrap()
        else {
            panic!("expected line carrier")
        };
        let whole = g.add_curve(Curve::Line(
            Line::new(
                pi([0, 0, 0]).into(),
                dir.clone(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ));
        assert_eq!(
            same_curve(&g, curve, whole, &mut Decider::default()),
            Outcome::Proven(Apart)
        );
        let position = g.curve(curve).parameter_origin(&g);
        let behind = position.plus(&position.minus(&pi([0, 0, 0]).to_real()).scaled(&rr(-2)));
        let point = closed(&mut g, behind, &supports);
        assert_eq!(
            on_curve(&g, curve, point, &mut Decider::default()),
            Outcome::Proven(Apart)
        );
        coincident(on_curve(&g, curve, singular[0], &mut Decider::default()));
    }
    let (g, _, result) = section(
        plane(pi([0, 0, 0]), vi([1, 0, -1])),
        cone(pi([0, 0, 0]), vi([0, 0, 1]), n(1)),
    );
    let SurfacePair::Meets(parts) = result else {
        panic!("expected tangent generator")
    };
    let [
        Component::Curve {
            curve,
            contact: Contact::Tangent,
            ..
        },
    ] = parts.as_slice()
    else {
        panic!("expected one tangent ray")
    };
    let curve = *curve;
    let Curve::Intersection(definition) = g.curve(curve) else {
        panic!("a ray retains its nappe definition")
    };
    let tangent = tangent_at(&g, curve, definition.witness, &mut Decider::default()).unwrap();
    positive_multiple(&tangent[1], &vi([1, 0, 1]).to_real());
}

#[test]
fn torus_sections_keep_concentric_and_axial_components_distinct() {
    let t = torus(pi([0, 0, 0]), n(3), n(1));
    circles(
        axis_plane(2, r("1/2")),
        t.clone(),
        &[
            (
                p(["0", "0", "1/2"]).to_real(),
                Real::from(r("39/4")) + rr(3) * root(&n(3)),
                Contact::Crossing,
            ),
            (
                p(["0", "0", "1/2"]).to_real(),
                Real::from(r("39/4")) - rr(3) * root(&n(3)),
                Contact::Crossing,
            ),
        ],
        &[],
    );
    circles(
        axis_plane(1, n(0)),
        t.clone(),
        &[
            (pi([3, 0, 0]).to_real(), rr(1), Contact::Crossing),
            (pi([-3, 0, 0]).to_real(), rr(1), Contact::Crossing),
        ],
        &[],
    );
    circles(
        axis_plane(2, n(1)),
        t.clone(),
        &[(pi([0, 0, 1]).to_real(), rr(9), Contact::Tangent)],
        &[],
    );
    let (_, _, result) = section(axis_plane(2, n(2)), t.clone());
    assert_eq!(result, SurfacePair::Disjoint);
    let (g, _, result) = section(axis_plane(0, n(4)), t);
    let SurfacePair::Meets(parts) = result else {
        panic!("expected support tangency")
    };
    let [Component::Point { point }] = parts.as_slice() else {
        panic!("expected isolated point")
    };
    point_eq(&g.position(*point), &pi([4, 0, 0]).to_real());
}

#[test]
fn coaxial_cone_and_torus_pairs_reduce_to_exact_circles() {
    let c = cone(pi([0, 0, 0]), vi([0, 0, 1]), n(1));
    let t = torus(pi([0, 0, 0]), n(3), n(1));
    circles(
        c.clone(),
        zc(2),
        &[(pi([0, 0, 2]).to_real(), rr(4), Contact::Crossing)],
        &[],
    );
    circles(
        c.clone(),
        sphere(pi([0, 0, 3]), n(3)),
        &[(pi([0, 0, 3]).to_real(), rr(9), Contact::Crossing)],
        &[pi([0, 0, 0]).to_real()],
    );
    circles(
        c.clone(),
        cone(pi([0, 0, 4]), vi([0, 0, -1]), n(1)),
        &[(pi([0, 0, 2]).to_real(), rr(4), Contact::Crossing)],
        &[],
    );
    circles(
        c,
        torus(pi([0, 0, 3]), n(3), n(1)),
        &[
            (
                Point3::new(rr(0), rr(0), rr(3) + root(&n(2)) / rr(2)),
                Real::from(r("19/2")) + rr(3) * root(&n(2)),
                Contact::Crossing,
            ),
            (
                Point3::new(rr(0), rr(0), rr(3) - root(&n(2)) / rr(2)),
                Real::from(r("19/2")) - rr(3) * root(&n(2)),
                Contact::Crossing,
            ),
        ],
        &[],
    );
    circles(
        cyl(pi([0, 0, 0]), vi([0, 0, 1]), r("7/2")),
        t.clone(),
        &[
            (
                Point3::new(rr(0), rr(0), root(&n(3)) / rr(2)),
                Real::from(r("49/4")),
                Contact::Crossing,
            ),
            (
                Point3::new(rr(0), rr(0), -root(&n(3)) / rr(2)),
                Real::from(r("49/4")),
                Contact::Crossing,
            ),
        ],
        &[],
    );
    circles(
        sphere(pi([0, 0, 0]), n(3)),
        t.clone(),
        &[
            (
                Point3::new(rr(0), rr(0), root(&n(35)) / rr(6)),
                Real::from(r("289/36")),
                Contact::Crossing,
            ),
            (
                Point3::new(rr(0), rr(0), -root(&n(35)) / rr(6)),
                Real::from(r("289/36")),
                Contact::Crossing,
            ),
        ],
        &[],
    );
    circles(
        t.clone(),
        torus(pi([0, 0, 1]), n(3), n(1)),
        &[
            (
                p(["0", "0", "1/2"]).to_real(),
                Real::from(r("39/4")) + rr(3) * root(&n(3)),
                Contact::Crossing,
            ),
            (
                p(["0", "0", "1/2"]).to_real(),
                Real::from(r("39/4")) - rr(3) * root(&n(3)),
                Contact::Crossing,
            ),
        ],
        &[],
    );
    circles(
        t,
        torus(pi([0, 0, 0]), n(4), n(1)),
        &[
            (
                Point3::new(rr(0), rr(0), root(&n(3)) / rr(2)),
                Real::from(r("49/4")),
                Contact::Crossing,
            ),
            (
                Point3::new(rr(0), rr(0), -root(&n(3)) / rr(2)),
                Real::from(r("49/4")),
                Contact::Crossing,
            ),
        ],
        &[],
    );
}

#[test]
fn sphere_inside_tube_and_meridian_tangency_are_certified() {
    let t = torus(pi([0, 0, 0]), n(3), n(1));
    for center in [pi([3, 0, 0]), pi([0, 0, 5]), pi([9, 0, 0])] {
        let (_, _, result) = section(sphere(center, r("1/2")), t.clone());
        assert_eq!(result, SurfacePair::Disjoint);
    }
    circles(
        sphere(pi([3, 0, 0]), n(1)),
        t,
        &[(pi([3, 0, 0]).to_real(), rr(1), Contact::Tangent)],
        &[],
    );
}

#[test]
fn torus_line_contact_uses_the_first_nonzero_derivative() {
    let mut g = GeomStore::new();
    let origin = p(["12/5", "0", "4/5"]);
    let direction = vi([-4, 10, -3]);
    let curve = g.add_curve(Curve::Line(
        Line::new(
            origin.clone().into(),
            direction.clone().into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let surface = g.add_surface(torus(pi([0, 0, 0]), n(3), n(1)));
    // F(origin+t*direction) = t^3 (15625 t - 6000): both roots cross.
    let CurveSurface::Meets(hits) = exact(intersect_curve_surface(
        &mut g,
        curve,
        surface,
        &mut Decider::default(),
    )) else {
        panic!("expected two hits")
    };
    assert_eq!(hits.len(), 2);
    for (hit, t) in hits.iter().zip([r("0"), r("48/125")]) {
        assert_eq!(hit.contact, Contact::Crossing);
        point_eq(
            &g.position(hit.point),
            &origin.plus(&direction.scaled(&t)).to_real(),
        );
    }
}

#[test]
fn circle_results_are_covariant_under_generated_rational_poses() {
    let mut rng = Rng(0x434f4e4943);
    for _ in 0..16 {
        let pose = rng.pose();
        let a = cone(pi([0, 0, 0]), vi([0, 0, 1]), n(1)).transformed(&pose);
        let b = cyl(pi([0, 0, 0]), vi([0, 0, -1]), n(2)).transformed(&pose);
        let expected = pose.apply_point(&pi([0, 0, 2]).to_real());
        circles(
            a.clone(),
            b.clone(),
            &[(expected.clone(), rr(4), Contact::Crossing)],
            &[],
        );
        circles(b, a, &[(expected, rr(4), Contact::Crossing)], &[]);
    }
}
