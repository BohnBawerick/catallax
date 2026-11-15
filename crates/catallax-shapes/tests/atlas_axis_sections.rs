//! Complete chart rows retain physical ends and shared seam ownership.

use catallax_numbers::{Decider, Outcome, Rational, Real, ZeroProof, with_operation_budget};
use catallax_shapes::{BSplineSurface, Plane, Point3, RationalAtlas, SectionSurface, Vector3};

fn q(value: &str) -> Real {
    Rational::parse(value).unwrap().into()
}

fn point(u: i64, v: i64, scale: i64, transpose: bool) -> Point3 {
    let x = Real::from_i64(3 + u * scale);
    let y = Real::from_i64(-5 + v * scale);
    if transpose {
        Point3::new(y, x, Real::from_i64(7))
    } else {
        Point3::new(x, y, Real::from_i64(7))
    }
}

fn sheet(scale: i64, transpose: bool, dec: &mut Decider) -> RationalAtlas {
    let controls = (0..3)
        .flat_map(|u| (0..2).map(move |v| point(u, v, scale, transpose)))
        .collect();
    let source = BSplineSurface::new(
        [1, 1],
        [
            [0, 0, 1, 2, 2].map(Rational::from_i64).to_vec(),
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        ],
        [3, 2],
        controls,
        vec![Real::one(); 6],
        dec,
    )
    .unwrap();
    let mut atlases = RationalAtlas::partition(source, [false, false], dec).unwrap();
    assert_eq!(atlases.len(), 1);
    let atlas = atlases.pop().unwrap();
    assert_eq!(atlas.charts().len(), 2);
    atlas
}

fn plane(axis: usize, at: Real, reversed: bool) -> SectionSurface {
    let mut origin = [Real::zero(), Real::zero(), Real::zero()];
    origin[axis] = at;
    let mut normal = [Real::zero(), Real::zero(), Real::zero()];
    normal[axis] = Real::from_i64(if reversed { -1 } else { 1 });
    let mut reference = [Real::zero(), Real::zero(), Real::zero()];
    reference[(axis + 1) % 3] = Real::one();
    let [x, y, z] = origin;
    let [nx, ny, nz] = normal;
    let [rx, ry, rz] = reference;
    SectionSurface::Plane(Plane {
        origin: Point3::new(x, y, z),
        normal: Vector3::new(nx, ny, nz),
        ref_dir: Vector3::new(rx, ry, rz),
    })
}

fn equal(a: &Point3, b: &Point3, dec: &mut Decider) -> bool {
    matches!(
        dec.sign(&a.minus(b).norm_squared()),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound
    )
}

#[test]
fn physical_rows_join_and_interior_seam_rows_have_one_owner() {
    for scale in [1, 7] {
        for transpose in [false, true] {
            for reversed in [false, true] {
                let ((), work) = with_operation_budget(20_000_000, || {
                    let mut dec = Decider::default();
                    let atlas = sheet(scale, transpose, &mut dec);
                    for seam in [false, true] {
                        let axis = usize::from(if seam { transpose } else { !transpose });
                        let at = Real::from_i64(if seam { 3 + scale } else { -5 });
                        let section = atlas
                            .sections(&plane(axis, at, reversed), &mut dec)
                            .unwrap();
                        assert!(section.points().is_empty());
                        assert!(section.branches().is_empty());
                        let curves: Vec<_> = section.curves().collect();
                        assert_eq!(curves.len(), 1);
                        let (curve, crossing) = curves[0];
                        assert!(crossing);
                        let ends = curve
                            .domain()
                            .map(|t| curve.point_at(&t.clone().into(), &mut dec).unwrap());
                        let expected = if seam {
                            [point(1, 0, scale, transpose), point(1, 1, scale, transpose)]
                        } else {
                            [point(0, 0, scale, transpose), point(2, 0, scale, transpose)]
                        };
                        assert!(
                            expected
                                .iter()
                                .all(|p| ends.iter().any(|e| equal(p, e, &mut dec)))
                        );
                        assert_eq!(curve.spans().count(), if seam { 1 } else { 2 });
                    }
                    for (offset, count) in [(q("-1/1000000"), 0), (q("1/1000000"), 1)] {
                        let section = atlas
                            .sections(
                                &plane(
                                    usize::from(!transpose),
                                    Real::from_i64(-5) + offset,
                                    reversed,
                                ),
                                &mut dec,
                            )
                            .unwrap();
                        assert_eq!(section.curves().count() + section.branches().len(), count);
                    }
                });
                assert!(!work.exhausted, "{} work units", work.used);
            }
        }
    }
}

fn ruled_row(extra_crossing: bool, dec: &mut Decider) -> RationalAtlas {
    let mut controls = Vec::new();
    for u in 0..2 {
        for v in 0..3 {
            let height = if extra_crossing {
                Real::from_i64(v) - q("1/2")
            } else {
                Real::from_i64(2 + v)
            };
            controls.push(Point3::new(u.into(), Real::from_i64(u) * height, v.into()));
        }
    }
    let source = BSplineSurface::new(
        [1, 1],
        [
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
            [0, 0, 1, 2, 2].map(Rational::from_i64).to_vec(),
        ],
        [2, 3],
        controls,
        vec![Real::one(); 6],
        dec,
    )
    .unwrap();
    let mut atlases = RationalAtlas::partition(source, [false, false], dec).unwrap();
    assert_eq!(atlases.len(), 1);
    atlases.pop().unwrap()
}

#[test]
fn row_factors_need_a_nonzero_remainder_over_the_whole_chart() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let atlas = ruled_row(false, &mut dec);
        let cut = plane(1, Real::zero(), false);
        let section = atlas.sections(&cut, &mut dec).unwrap();
        let curves: Vec<_> = section.curves().collect();
        assert_eq!(curves.len(), 1);
        let (curve, crossing) = curves[0];
        assert!(crossing);
        assert_eq!(curve.spans().count(), 2);
        let ends = curve
            .domain()
            .map(|t| curve.point_at(&t.clone().into(), &mut dec).unwrap());
        assert!(equal(
            &ends[0],
            &Point3::new(0.into(), 0.into(), 0.into()),
            &mut dec
        ));
        assert!(equal(
            &ends[1],
            &Point3::new(0.into(), 0.into(), 2.into()),
            &mut dec
        ));
        let crossing = ruled_row(true, &mut dec);
        let error = crossing.sections(&cut, &mut dec).unwrap_err();
        assert_eq!(error.code, "intersect.atlas-section-row");
    });
    assert!(!work.exhausted, "{} work units", work.used);
}

#[test]
fn a_cubic_loft_row_keeps_its_finite_span_after_cancellation() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut controls = Vec::new();
        for (u, row) in [(q("0"), q("0")), (q("1/8"), q("1/8")), (q("1/4"), q("1/8"))] {
            for (height, v) in [(24, q("1")), (24, q("4/3")), (20, q("5/3")), (16, q("2"))] {
                controls.push(Point3::new(u.clone(), &row * &Real::from_i64(height), v));
            }
        }
        let source = BSplineSurface::new(
            [2, 3],
            [
                ["0", "0", "0", "1/4", "1/4", "1/4"]
                    .map(|t| Rational::parse(t).unwrap())
                    .to_vec(),
                [1, 1, 1, 1, 2, 2, 2, 2].map(Rational::from_i64).to_vec(),
            ],
            [3, 4],
            controls,
            vec![Real::one(); 12],
            &mut dec,
        )
        .unwrap();
        let mut atlases = RationalAtlas::partition(source, [false, false], &mut dec).unwrap();
        assert_eq!(atlases.len(), 1);
        let section = atlases
            .pop()
            .unwrap()
            .sections(&plane(1, Real::zero(), false), &mut dec)
            .unwrap();
        let curves: Vec<_> = section.curves().collect();
        assert_eq!(curves.len(), 1);
        let (curve, crossing) = curves[0];
        assert!(crossing);
        let ends = curve
            .domain()
            .map(|t| curve.point_at(&t.clone().into(), &mut dec).unwrap());
        for z in [1, 2] {
            assert!(ends.iter().any(|p| equal(
                p,
                &Point3::new(0.into(), 0.into(), z.into()),
                &mut dec
            )));
        }
    });
    assert!(!work.exhausted, "{} work units", work.used);
}
