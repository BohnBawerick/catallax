//! Independent polar-law and cap-intersection evidence, without publishing a body.
use catallax_features::draft::{
    DraftDefinition,
    annular::{AnnularSeparation, AnnularWall, ObliqueAnnularBound},
};
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real};
use catallax_shapes::{Cylinder, Plane, Point3, Vector3};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn r(s: &str) -> Real {
    q(s).into()
}
fn p(x: &str, y: &str, z: &str) -> Point3 {
    Point3::new(r(x), r(y), r(z))
}
fn v(x: &str, y: &str, z: &str) -> Vector3 {
    Vector3::new(r(x), r(y), r(z))
}
fn fixture(rate: &str, reversed: bool) -> ObliqueAnnularBound {
    let mut walls = [
        AnnularWall {
            cylinder: Cylinder {
                origin: p("0", "0", "0"),
                axis: v("0", "0", "2"),
                ref_dir: v("3", "0", "0"),
                radius_sq: r("64"),
            },
            origins: vec![],
            same_sense: true,
            selected: true,
        },
        AnnularWall {
            cylinder: Cylinder {
                // Axial origin, axis sense and azimuth do not change the cylinder.
                origin: p("0", "0", "7"),
                axis: v("0", "0", "-1"),
                ref_dir: v("0", "2", "0"),
                radius_sq: r("16"),
            },
            origins: vec![],
            same_sense: false,
            selected: true,
        },
    ];
    let mut caps = [
        (
            Plane {
                origin: p("0", "0", "0"),
                normal: v("0", "0", "-2"),
                ref_dir: v("1", "0", "0"),
            },
            vec![],
        ),
        (
            Plane {
                origin: p("0", "0", "2"),
                normal: v("0", "0", "3"),
                ref_dir: v("1", "0", "0"),
            },
            vec![],
        ),
    ];
    if reversed {
        walls.swap(0, 1);
        caps.swap(0, 1);
    }
    let neutral = Plane {
        origin: p("0", "0", "1"),
        normal: v("1", "0", rate),
        ref_dir: v("0", "1", "0"),
    };
    let definition =
        DraftDefinition::new(v("0", "0", "1"), neutral, &q("45"), &mut Decider::default()).unwrap();
    ObliqueAnnularBound::new(walls, caps, &definition, &mut Decider::default()).unwrap()
}

#[test]
fn complete_finite_gap_has_independently_attained_extrema_and_exact_cap_curves() {
    for reversed in [false, true] {
        let bound = fixture("8", reversed);
        // Independently solve x+8(z-1)=0 on each source cylinder and
        // advance the authored generators. Outer rho=9-z-cos(u),
        // inner rho=3+z+cos(u)/2; gap=6-2z-3cos(u)/2.
        assert_eq!(bound.separation(), AnnularSeparation::Separated);
        assert_eq!(bound.radial_gap(), &[q("1/2"), q("15/2")]);
        for z in ["0", "1/3", "1", "5/3", "2"].map(q) {
            for [c, s] in [["1", "0"], ["-1", "0"], ["3/5", "4/5"]] {
                let (c, s) = (q(c), q(s));
                let radii = [&(&q("9") - &z) - &c, &(&q("3") + &z) + &(&c / &q("2"))];
                let gap = &radii[0] - &radii[1];
                assert!(bound.radial_gap()[0] <= gap && gap <= bound.radial_gap()[1]);
                for j in 0..2 {
                    let point = p_from(&radii[j] * &c, &radii[j] * &s, z.clone());
                    assert_eq!(
                        bound.walls()[j]
                            .support
                            .side_at(&point, &mut Decider::default())
                            .unwrap(),
                        catallax_numbers::Sign::Zero
                    );
                }
            }
        }
        for (i, (plane, rings)) in bound.cap_intersections().iter().enumerate() {
            let z = Rational::from_i64(2 * i as i64);
            for quarter in 0..4 {
                for tau in ["0", "1/3", "1"] {
                    let tau = q(tau);
                    let den = Rational::one() + &tau * &tau;
                    let c = (Rational::one() - &tau * &tau) / den.clone();
                    let s = (&Rational::from_i64(2) * &tau) / den;
                    let [c, s] = match quarter {
                        0 => [c, s],
                        1 => [-s, c],
                        2 => [-c, -s],
                        _ => [s, -c],
                    };
                    let radii = [&(&q("9") - &z) - &c, &(&q("3") + &z) + &(&c / &q("2"))];
                    let gap = &radii[0] - &radii[1];
                    assert!(bound.radial_gap()[0] <= gap && gap <= bound.radial_gap()[1]);
                    let u = Rational::from_i64(quarter) + tau;
                    for j in 0..2 {
                        let want = [&radii[j] * &c, &radii[j] * &s, z.clone()];
                        let got = rings[j]
                            .evaluate(&Ball::from_rational(&u, 192), &Decider::default())
                            .unwrap();
                        for (got, want) in got.iter().zip(want) {
                            assert!(got.lo() <= want && want <= got.hi());
                            assert!(got.hi() - got.lo() < q("1e-40"));
                        }
                        assert_eq!(rings[j].is_closed(), Some(true));
                        let point = p_from(&radii[j] * &c, &radii[j] * &s, z.clone());
                        assert!(matches!(
                            Decider::default().sign(&point.minus(&plane.origin).dot(&plane.normal)),
                            Outcome::ProvenCoincident(_)
                        ));
                        assert_eq!(
                            bound.walls()[j]
                                .support
                                .side_at(&point, &mut Decider::default())
                                .unwrap(),
                            catallax_numbers::Sign::Zero
                        );
                    }
                }
            }
        }
        // Both exact extrema occur at cap intersections, including the
        // periodic seam (the minimum at z=2, cos=1).
        for (z, c, gap) in [("2", "1", "1/2"), ("0", "-1", "15/2")] {
            assert_eq!(q("6") - q("2") * q(z) - q("3/2") * q(c), q(gap));
        }
    }
}
fn p_from(x: Rational, y: Rational, z: Rational) -> Point3 {
    Point3::new(x.into(), y.into(), z.into())
}

#[test]
fn touching_and_crossing_are_exact_relations_not_numeric_coverage_failures() {
    let touching = fixture("6", false);
    assert_eq!(touching.separation(), AnnularSeparation::Touching);
    assert_eq!(touching.radial_gap()[0], q("0"));
    for wall in touching.walls() {
        assert_eq!(
            wall.support
                .side_at(&p("17/3", "0", "2"), &mut Decider::default())
                .unwrap(),
            catallax_numbers::Sign::Zero
        );
    }
    let crossing = fixture("4", false);
    assert_eq!(crossing.separation(), AnnularSeparation::Crossing);
    assert_eq!(crossing.radial_gap(), &[q("-1"), q("9")]);
    // At the upper cap, cos(u)=2/3 and sin(u)=sqrt(5)/3 give
    // equal radii 17/3. This is an actual cap-curve crossing.
    let point = Point3::new(r("34/9"), &r("17/9") * &r("5").sqrt(), r("2"));
    for wall in crossing.walls() {
        assert_eq!(
            wall.support
                .side_at(&point, &mut Decider::default())
                .unwrap(),
            catallax_numbers::Sign::Zero
        );
    }
}
