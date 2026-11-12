use super::*;
use catallax_numbers::{Real, with_operation_budget};
use catallax_shapes::{
    BSplineSurface, Curve, Plane, PointDef, RationalPatch, SectionSurface, Vector3,
};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn q(n: i64, d: i64) -> Rational {
    Rational::from_ratio(n, d).unwrap()
}
fn xyz(p: &Point3) -> [&Real; 3] {
    [&p.x, &p.y, &p.z]
}
fn zero(value: &Real, dec: &mut Decider) {
    assert!(dec.sign(value).is_coincident(), "{value:?}");
}

fn fixture(dec: &mut Decider) -> (GeomStore, CurveId) {
    let mut controls = Vec::new();
    for (x, zx) in [(-1, 1), (0, -1), (1, 1)] {
        for (y, zy) in [(-1, 1), (0, -1), (1, 1)] {
            controls.push(Point3::new(x.into(), y.into(), (zx + zy).into()));
        }
    }
    let knots = vec![r(-1), r(-1), r(-1), r(1), r(1), r(1)];
    let patch = RationalPatch::new(
        BSplineSurface::new(
            [2, 2],
            [knots.clone(), knots],
            [3, 3],
            controls,
            vec![Real::one(); 9],
            dec,
        )
        .unwrap(),
        dec,
    )
    .unwrap();
    // 3/8 has no rational two-square representation, so cell corners are off the cut.
    let cut = SectionSurface::Plane(Plane {
        origin: Point3::new(Real::zero(), Real::zero(), q(3, 8).into()),
        normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
    });
    let section = patch.sections(&cut, dec).unwrap();
    let branch = section
        .branches()
        .first()
        .expect("a closed regular section");
    assert!(branch.is_closed());
    let mut geom = GeomStore::new();
    let curve = geom.add_curve(Curve::PatchSection(branch.clone()));
    (geom, curve)
}

#[test]
fn retained_closed_section_traces_enclose_exact_points_and_implicit_jets() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let (geom, curve) = fixture(&mut dec);
        let traced = trace(&geom, curve, None, None, &mut dec).unwrap();
        assert!(traced.marks.is_empty());
        let mut first = None;
        let mut previous: Option<Point3> = None;
        for piece in traced.value.pieces() {
            let entire = Ball::from_rational(&r(0), 128).hull(&Ball::from_rational(&r(1), 128));
            let enclosure = piece.enclose(&entire, &dec).unwrap();
            let start = piece.evaluate_exact(&r(0), &mut dec).unwrap();
            if let Some(previous) = &previous {
                zero(&start.minus(previous).norm_squared(), &mut dec);
            }
            if first.is_none() {
                first = Some(start);
            }
            for t in [r(0), q(1, 2), r(1)] {
                let exact = piece.evaluate_exact(&t, &mut dec).unwrap();
                zero(&(&exact.z - &Real::from(q(3, 8))), &mut dec);
                zero(
                    &(&(&exact.x.square() + &exact.y.square()) - &Real::from(q(3, 8))),
                    &mut dec,
                );
                for (coordinate, bound) in xyz(&exact).into_iter().zip(&enclosure) {
                    for margin in [
                        coordinate - &Real::from(bound.lo()),
                        &Real::from(bound.hi()) - coordinate,
                    ] {
                        assert!(matches!(
                            dec.sign(&margin),
                            catallax_numbers::Outcome::Proven(
                                catallax_numbers::StrictSign::Positive
                            ) | catallax_numbers::Outcome::ProvenCoincident(_)
                        ));
                    }
                }
            }
            let jets = piece
                .taylor(&Ball::from_rational(&q(1, 2), 128), 3, &dec)
                .unwrap();
            assert!(jets[1][2].contains_zero());
            let tangent_contact = &(&jets[0][0] * &jets[1][0]) + &(&jets[0][1] * &jets[1][1]);
            assert!(tangent_contact.contains_zero());
            let speed = &jets[1][0].square() + &jets[1][1].square();
            assert!(speed.lo() > r(0));
            let filtered = piece
                .taylor_interval(catallax_numbers::Interval::point(0.5).unwrap(), 3, &dec)
                .unwrap()
                .unwrap();
            assert!(filtered[1][2].lo() <= 0.0 && filtered[1][2].hi() >= 0.0);
            previous = Some(piece.evaluate_exact(&r(1), &mut dec).unwrap());
        }
        zero(
            &first.unwrap().minus(&previous.unwrap()).norm_squared(),
            &mut dec,
        );
    });
    assert!(!work.exhausted);
}

#[test]
fn retained_section_trace_clips_interior_endpoints_instead_of_replaying_the_component() {
    let mut dec = Decider::default();
    let (mut geom, curve) = fixture(&mut dec);
    let whole = trace(&geom, curve, None, None, &mut dec).unwrap();
    let source = &whole.value.pieces()[0];
    let a = source.evaluate_exact(&q(1, 4), &mut dec).unwrap();
    let b = source.evaluate_exact(&q(3, 4), &mut dec).unwrap();
    let expected_middle = source.evaluate_exact(&q(1, 2), &mut dec).unwrap();
    let from = geom.add_point(PointDef::Explicit(a.clone()));
    let to = geom.add_point(PointDef::Explicit(b.clone()));
    let clipped = trace(&geom, curve, Some(from), Some(to), &mut dec).unwrap();
    let first = clipped.value.pieces().first().unwrap();
    let last = clipped.value.pieces().last().unwrap();
    zero(
        &first
            .evaluate_exact(&r(0), &mut dec)
            .unwrap()
            .minus(&a)
            .norm_squared(),
        &mut dec,
    );
    zero(
        &last
            .evaluate_exact(&r(1), &mut dec)
            .unwrap()
            .minus(&b)
            .norm_squared(),
        &mut dec,
    );
    zero(
        &first
            .evaluate_exact(&q(1, 2), &mut dec)
            .unwrap()
            .minus(&expected_middle)
            .norm_squared(),
        &mut dec,
    );
    assert_eq!(dec.stats().coincident_separation, 0);
}

#[test]
fn retained_section_meetings_keep_algebraic_cut_coefficients_exact() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let (geom, curve) = fixture(&mut dec);
        let Curve::PatchSection(branch) = geom.curve(curve) else {
            unreachable!()
        };
        let offset = Real::from(q(3, 8)).sqrt() * Real::from(q(1, 2));
        let cutting = SectionSurface::Plane(Plane {
            origin: Point3::new(offset.clone(), Real::zero(), Real::zero()),
            normal: Vector3::new(Real::one(), Real::zero(), Real::zero()),
            ref_dir: Vector3::new(Real::zero(), Real::one(), Real::zero()),
        });
        let catallax_shapes::SectionHits::Points(points) =
            branch.surface_hits(&cutting, &mut dec).unwrap()
        else {
            panic!("a transverse plane has finite circle meetings")
        };
        assert_eq!(points.len(), 2);
        let mut signs = Vec::new();
        for point in points {
            zero(&(&point.x - &offset), &mut dec);
            zero(&(&point.z - &Real::from(q(3, 8))), &mut dec);
            zero(
                &(&(&point.x.square() + &point.y.square()) - &Real::from(q(3, 8))),
                &mut dec,
            );
            signs.push(match dec.sign(&point.y) {
                catallax_numbers::Outcome::Proven(sign) => sign,
                other => panic!("a transverse meeting needs a strict side: {other:?}"),
            });
        }
        assert!(signs.contains(&catallax_numbers::StrictSign::Positive));
        assert!(signs.contains(&catallax_numbers::StrictSign::Negative));
        assert_eq!(dec.stats().coincident_separation, 0);
    });
    assert!(!work.exhausted);
}

#[test]
fn a_section_and_its_endpoint_chord_are_not_coincident() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut controls = Vec::new();
        for (x, y) in [(r(0), r(0)), (q(1, 2), q(1, 2)), (r(1), r(0))] {
            for z in [r(0), r(1)] {
                controls.push(Point3::new(x.clone().into(), y.clone().into(), z.into()));
            }
        }
        let carrier = BSplineSurface::new(
            [2, 1],
            [
                vec![r(0), r(0), r(0), r(1), r(1), r(1)],
                vec![r(0), r(0), r(1), r(1)],
            ],
            [3, 2],
            controls,
            vec![Real::one(); 6],
            &mut dec,
        )
        .unwrap();
        let patch = RationalPatch::new(carrier, &mut dec).unwrap();
        let section = patch
            .sections(
                &SectionSurface::Plane(Plane {
                    origin: Point3::new(Real::zero(), Real::zero(), q(1, 2).into()),
                    normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
                    ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
                }),
                &mut dec,
            )
            .unwrap();
        assert_eq!(section.branches().len(), 1);
        let branch = section.branches()[0].clone();
        for endpoint in branch.endpoints() {
            zero(&endpoint.y, &mut dec);
            zero(&(&endpoint.z - &Real::from(q(1, 2))), &mut dec);
        }
        let mut geom = GeomStore::new();
        let curve = geom.add_curve(Curve::PatchSection(branch));
        let chord = geom.add_curve(Curve::Line(catallax_shapes::Line {
            origin: Point3::new(Real::zero(), Real::zero(), q(1, 2).into()),
            dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        }));
        for (a, b) in [(curve, chord), (chord, curve)] {
            assert!(matches!(
                same_curve(&geom, a, b, &mut dec),
                catallax_numbers::Outcome::Proven(Apart)
            ));
        }
        assert_eq!(dec.stats().coincident_separation, 0);
    });
    assert!(!work.exhausted);
}
