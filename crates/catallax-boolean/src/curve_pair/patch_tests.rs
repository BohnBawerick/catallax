use super::*;
use catallax_numbers::{Rational, ZeroProof, with_operation_budget};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Cylinder, Line, Plane, RationalPatch, SectionSurface, Surface,
    Vector3,
};

fn q(n: i64, d: i64) -> Real {
    Rational::from_ratio(n, d).unwrap().into()
}

#[test]
fn finite_line_covers_keep_native_domains_and_do_not_mark_whole_curves_coincident() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        for reverse in [false, true] {
            for direction in [-2, 2] {
                let mut model = Model::new();
                let mut controls = vec![
                    Point3::new(q(1, 1), q(0, 1), q(0, 1)),
                    Point3::new(q(3, 1), q(0, 1), q(0, 1)),
                ];
                if reverse {
                    controls.reverse();
                }
                let spline = BSplineCurve::new(
                    1,
                    [2, 2, 5, 5].map(Rational::from_i64).to_vec(),
                    controls.clone(),
                    vec![Real::one(); 2],
                    &mut dec,
                )
                .unwrap();
                let finite = model.geom_mut().add_curve(Curve::BSpline(spline));
                let line = model.geom_mut().add_curve(Curve::Line(Line {
                    origin: Point3::new(q(0, 1), q(0, 1), q(0, 1)),
                    dir: Vector3::new(q(direction, 1), q(0, 1), q(0, 1)),
                }));
                for (a, b) in [(finite, line), (line, finite)] {
                    let cover = finite_line_cover(&mut model, a, b, &mut dec)
                        .unwrap()
                        .unwrap();
                    assert_eq!((cover.finite, cover.carrier), (finite, line));
                    for (point, expected) in cover.ends.into_iter().zip(&controls) {
                        exact_zero(
                            &model.geom().position(point).minus(expected).norm_squared(),
                            &mut dec,
                        );
                    }
                    let Curve::BSpline(source) = model.geom().curve(finite) else {
                        unreachable!()
                    };
                    assert_eq!(
                        source.domain(),
                        [&Rational::from_i64(2), &Rational::from_i64(5)]
                    );
                    assert!(matches!(
                        same_curve(model.geom(), finite, line, &mut dec),
                        Outcome::Proven(false)
                    ));
                    assert!(!model.geom().facts().holds(
                        FactKind::Coincident,
                        GeomKey::Curve(finite),
                        GeomKey::Curve(line)
                    ));
                }
                let off = model.geom_mut().add_curve(Curve::Line(Line {
                    origin: Point3::new(q(0, 1), q(1, 1), q(0, 1)),
                    dir: Vector3::new(q(1, 1), q(0, 1), q(0, 1)),
                }));
                assert!(
                    finite_line_cover(&mut model, finite, off, &mut dec)
                        .unwrap()
                        .is_none()
                );
            }
        }
    });
    assert!(!work.exhausted, "{work:?}");
}

#[test]
fn native_straight_spline_meetings_keep_finite_ends_and_operand_senses() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        for reverse in [false, true] {
            for (z, expected) in [(q(1, 2), 1), (q(2, 1), 0)] {
                let mut model = Model::new();
                let mut ends = vec![
                    Point3::new(q(0, 1), q(0, 1), q(0, 1)),
                    Point3::new(q(0, 1), q(0, 1), q(1, 1)),
                ];
                if reverse {
                    ends.reverse();
                }
                let line = BSplineCurve::new(
                    1,
                    [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
                    ends,
                    vec![Real::one(); 2],
                    &mut dec,
                )
                .unwrap();
                let curved = BSplineCurve::new(
                    2,
                    [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
                    vec![
                        Point3::new(q(-1, 1), q(1, 1), z.clone()),
                        Point3::new(q(0, 1), q(-1, 1), z.clone()),
                        Point3::new(q(1, 1), q(1, 1), z.clone()),
                    ],
                    vec![Real::one(); 3],
                    &mut dec,
                )
                .unwrap();
                let a = model.geom_mut().add_curve(Curve::BSpline(line));
                let b = model.geom_mut().add_curve(Curve::BSpline(curved));
                for (a, b) in [(a, b), (b, a)] {
                    let points = positions(&mut model, a, b, &mut dec);
                    assert_eq!(points.len(), expected);
                    for point in points {
                        exact_zero(&point.x, &mut dec);
                        exact_zero(&point.y, &mut dec);
                        exact_zero(&(&point.z - &z), &mut dec);
                    }
                }
            }
        }
    });
    assert!(!work.exhausted, "{} units", work.used);
}

fn plane() -> SectionSurface {
    SectionSurface::Plane(Plane {
        origin: Point3::new(Real::zero(), Real::zero(), q(5, 2)),
        normal: Vector3::new(Real::one(), Real::one(), Real::one()),
        ref_dir: Vector3::new(Real::one(), Real::from_i64(-1), Real::zero()),
    })
}

fn cylinder() -> SectionSurface {
    SectionSurface::Cylinder(Cylinder {
        origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
        axis: Vector3::new(Real::zero(), Real::one(), Real::zero()),
        radius_sq: q(25, 4),
        ref_dir: Vector3::new(Real::zero(), Real::zero(), Real::one()),
    })
}

// The positive-weight loft-like patch of the section oracle tests.
fn patch(shift: i64, dec: &mut Decider) -> RationalPatch {
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    let bulge = [q(1, 1), q(2, 1), q(1, 2), q(3, 2)];
    for i in 0..3 {
        for (j, height) in bulge.iter().enumerate() {
            let y = if i == 1 { height.clone() } else { Real::zero() };
            controls.push(Point3::new(
                Real::from_i64(i),
                y + Real::from_i64(shift),
                Real::from_i64(j as i64),
            ));
            weights.push(if i == 1 { q(3, 4) } else { Real::one() });
        }
    }
    let knots = |values: &[i64]| values.iter().map(|n| Rational::from_i64(*n)).collect();
    let carrier = BSplineSurface::new(
        [2, 3],
        [knots(&[0, 0, 0, 1, 1, 1]), knots(&[0, 0, 0, 0, 1, 1, 1, 1])],
        [3, 4],
        controls,
        weights,
        dec,
    )
    .unwrap();
    RationalPatch::new(carrier, dec).unwrap()
}

fn exact_zero(value: &Real, dec: &mut Decider) {
    let Outcome::ProvenCoincident(proof) = dec.sign(value) else {
        panic!("expected an exact zero")
    };
    assert_ne!(proof.proof, ZeroProof::SeparationBound);
}

fn positions(model: &mut Model, a: CurveId, b: CurveId, dec: &mut Decider) -> Vec<Point3> {
    let Meetings::Points(points) = meetings(model, a, b, dec).unwrap() else {
        panic!("these sections have only isolated meetings")
    };
    points
        .into_iter()
        .map(|p| model.geom().position(p))
        .collect()
}

fn check_ends(points: &[Point3], dec: &mut Decider) {
    assert_eq!(points.len(), 2);
    for point in points {
        exact_zero(&point.y, dec);
        exact_zero(&(&point.x * &(&point.x - &Real::from_i64(2))), dec);
        exact_zero(&(&(&point.x + &point.z) - &q(5, 2)), dec);
    }
    assert!(matches!(
        dec.sign(&points[0].minus(&points[1]).norm_squared()),
        Outcome::Proven(_)
    ));
}

#[test]
fn patch_sections_meet_on_the_selected_components_in_both_orders() {
    for separate in [false, true] {
        for reverse in [false, true] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let mut model = Model::new();
                let source = patch(0, &mut dec);
                let section = source.sections(&plane(), &mut dec).unwrap();
                assert_eq!(section.branches().len(), 1);
                let a = model
                    .geom_mut()
                    .add_curve(Curve::PatchSection(section.branches()[0].clone()));
                let other = if separate {
                    patch(10, &mut dec)
                } else {
                    source
                };
                let section = other.sections(&cylinder(), &mut dec).unwrap();
                assert_eq!(section.branches().len(), 1);
                let b = model
                    .geom_mut()
                    .add_curve(Curve::PatchSection(section.branches()[0].clone()));
                let (a, b) = if reverse { (b, a) } else { (a, b) };
                let points = positions(&mut model, a, b, &mut dec);
                if separate {
                    assert!(points.is_empty());
                } else {
                    assert_eq!(points.len(), 1);
                    exact_zero(
                        &points[0]
                            .minus(&Point3::new(Real::zero(), Real::zero(), q(5, 2)))
                            .norm_squared(),
                        &mut dec,
                    );
                }
            });
            assert!(!work.exhausted, "{} units", work.used);
        }
    }
}

#[test]
fn patch_section_meets_a_line_contained_in_its_cutting_plane() {
    for reverse in [false, true] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let mut model = Model::new();
            let section = patch(0, &mut dec).sections(&plane(), &mut dec).unwrap();
            assert_eq!(section.branches().len(), 1);
            let a = model
                .geom_mut()
                .add_curve(Curve::PatchSection(section.branches()[0].clone()));
            let line = Line::new(
                Point3::new(Real::zero(), Real::zero(), q(5, 2)),
                Vector3::new(Real::one(), Real::zero(), Real::from_i64(-1)),
                &mut dec,
            )
            .unwrap();
            let b = model.geom_mut().add_curve(Curve::Line(line));
            let (a, b) = if reverse { (b, a) } else { (a, b) };
            check_ends(&positions(&mut model, a, b, &mut dec), &mut dec);
        });
        assert!(!work.exhausted, "{} units", work.used);
    }
}

#[test]
fn a_reversed_section_keeps_its_component_identity() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut model = Model::new();
        let section = patch(0, &mut dec).sections(&plane(), &mut dec).unwrap();
        let branch = section.branches()[0].clone();
        let a = model
            .geom_mut()
            .add_curve(Curve::PatchSection(branch.clone()));
        let b = model
            .geom_mut()
            .add_curve(Curve::PatchSection(branch.reversed()));
        assert!(matches!(
            meetings(&mut model, a, b, &mut dec),
            Ok(Meetings::Coincident)
        ));
    });
    assert!(!work.exhausted, "{} units", work.used);
}

#[test]
fn a_missing_common_surface_solver_is_never_an_empty_meeting() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut model = Model::new();
        let section = patch(0, &mut dec).sections(&plane(), &mut dec).unwrap();
        let a = model
            .geom_mut()
            .add_curve(Curve::PatchSection(section.branches()[0].clone()));
        let curve = BSplineCurve::new(
            1,
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
            vec![
                Point3::new(Real::zero(), Real::zero(), q(5, 2)),
                Point3::new(Real::from_i64(2), Real::zero(), q(1, 2)),
            ],
            vec![Real::one(); 2],
            &mut dec,
        )
        .unwrap();
        let b = model.geom_mut().add_curve(Curve::BSpline(curve));
        let Err(error) = meetings(&mut model, a, b, &mut dec) else {
            panic!("the unsupported finite curve query must refuse")
        };
        assert_eq!(error.code, "boolean.patch-section-meetings");
        assert!(error.detail.contains("BSpline"));
        assert!(error.detail.contains("solver"));
    });
    assert!(!work.exhausted, "{} units", work.used);
}

#[test]
fn spline_span_meetings_keep_a_shared_knot_once_and_refuse_partial_containment() {
    for partial in [false, true] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let mut model = Model::new();
            let section = patch(0, &mut dec).sections(&plane(), &mut dec).unwrap();
            let a = model
                .geom_mut()
                .add_curve(Curve::PatchSection(section.branches()[0].clone()));
            let controls = if partial {
                vec![
                    Point3::new(Real::zero(), Real::zero(), q(5, 2)),
                    Point3::new(Real::one(), Real::zero(), q(3, 2)),
                    Point3::new(Real::from_i64(2), Real::zero(), q(5, 2)),
                ]
            } else {
                [-1, 0, 1]
                    .map(|x| Point3::new(Real::from_i64(x), Real::zero(), q(5, 2)))
                    .to_vec()
            };
            let curve = BSplineCurve::new(
                1,
                [0, 0, 1, 2, 2].map(Rational::from_i64).to_vec(),
                controls,
                vec![Real::one(); 3],
                &mut dec,
            )
            .unwrap();
            let b = model.geom_mut().add_curve(Curve::BSpline(curve));
            if partial {
                let Err(error) = meetings(&mut model, a, b, &mut dec) else {
                    panic!("partial containment must refuse")
                };
                assert_eq!(error.code, "boolean.patch-section-meetings");
                assert!(error.detail.contains("some spans"));
            } else {
                let points = positions(&mut model, a, b, &mut dec);
                assert_eq!(points.len(), 1);
                exact_zero(
                    &points[0]
                        .minus(&Point3::new(Real::zero(), Real::zero(), q(5, 2)))
                        .norm_squared(),
                    &mut dec,
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }
}

#[test]
fn finite_atlas_incidence_uses_parameters_without_a_global_side() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut model = Model::new();
        let source = patch(0, &mut dec);
        let on = source.point_at(&[q(1, 2), q(1, 3)], &mut dec).unwrap();
        let mut atlases = catallax_shapes::RationalAtlas::partition(
            source.carrier().clone(),
            [false, false],
            &mut dec,
        )
        .unwrap();
        assert_eq!(atlases.len(), 1);
        let support = atlases.pop().unwrap();
        let off = Point3::new(Real::one(), q(1, 10), Real::one());
        assert_eq!(
            support.side_at(&off, &mut dec).unwrap_err().code,
            "shapes.parametric-query"
        );
        let surface = model
            .geom_mut()
            .add_surface(Surface::BSpline(support.into()));
        let face = model.add_face(catallax_topology::Face {
            surface,
            sense: catallax_topology::Sense::Same,
            loops: Vec::new(),
        });
        for (position, expected) in [(on, true), (off, false)] {
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit(position.clone()));
            assert_eq!(
                crate::geometry::on_surface(&model, surface, point, &mut dec, "test.incidence")
                    .unwrap(),
                expected
            );
            if !expected {
                assert!(!crate::face_contains(&mut model, face, &position, &mut dec).unwrap());
            }
        }
    });
    assert!(!work.exhausted, "{} units", work.used);
}
