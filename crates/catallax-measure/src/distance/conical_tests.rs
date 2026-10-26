use super::*;
use catallax_numbers::{Real, Sign};
use catallax_shapes::{Circle, Cone, PointDef, RVector3, Rotation};
use catallax_topology::{Coedge, Edge, Face, Vertex};

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    use catallax_numbers::{Outcome, StrictSign};
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        outcome => panic!("unproved witness control: {outcome:?}"),
    }
}

// A public synthetic full-band fixture, independent of the stock solid builder.
// Reversal reflects the axis; the three senses select the same finite region
// when effective travel is lower=increasing / upper=decreasing.
fn fixture(
    pose: &Pose,
    scale: &Rational,
    reverse: bool,
    senses: [Sense; 3],
    tip: bool,
) -> (Model, FaceId) {
    let mut model = Model::new();
    let sign = if reverse { -1 } else { 1 };
    let pt = |x, y, z| {
        pose.apply_point(
            &RPoint3::new(&q(x) * scale, &q(y) * scale, &q(z * sign) * scale).to_real(),
        )
    };
    let axis = pose.apply_vector(&RVector3::new(q(0), q(0), q(3 * sign)).to_real());
    let reference = pose.apply_vector(&RVector3::new(q(2), q(0), q(0)).to_real());
    let surface = model.geom_mut().add_surface(Surface::Cone(Cone {
        apex: pt(0, 0, 0),
        axis: axis.clone(),
        slope: Real::one(),
        ref_dir: reference.clone(),
    }));
    let mut loops = Vec::new();
    for height in if tip { vec![6] } else { vec![2, 6] } {
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: pt(0, 0, height),
            normal: axis.clone(),
            radius_sq: (&q(height * height) * &(scale * scale)).into(),
            ref_dir: reference.clone(),
        }));
        let edge = model.add_edge(Edge {
            curve,
            sense: senses[1],
            start: None,
            end: None,
        });
        let sense = if height == 2 {
            senses[2]
        } else {
            if senses[2] == Sense::Same {
                Sense::Reversed
            } else {
                Sense::Same
            }
        };
        let coedge = model.add_coedge(Coedge { edge, sense });
        loops.push(model.add_loop(Loop::Coedges(vec![coedge])));
    }
    if tip {
        let point = model.geom_mut().add_point(PointDef::Explicit(pt(0, 0, 0)));
        let vertex = model.add_vertex(Vertex { point });
        loops.push(model.add_loop(Loop::Vertex(vertex)));
    }
    let face = model.add_face(Face {
        surface,
        sense: senses[0],
        loops,
    });
    (model, face)
}
fn check(model: &mut Model, face: FaceId, at: Point3, expected: [Rational; 2]) -> Bounds {
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(model, face, &mut dec).unwrap();
    assert!(prepared.finite.is_some());
    let answer = prepared
        .point_to_face(model, &at, &controls(), &mut dec)
        .unwrap();
    assert!(answer.converged && answer.cells == 0, "{answer:?}");
    assert!(
        answer.lower <= expected[0] && answer.upper.as_ref().unwrap() >= &expected[1],
        "{answer:?} oracle={expected:?}"
    );
    assert!(answer.upper.as_ref().unwrap() - &answer.lower <= controls().width);
    let [_, witness] = answer.witnesses.as_ref().unwrap();
    assert!(
        catallax_boolean::face_contains(model, face, witness, &mut dec).unwrap(),
        "{witness:?}"
    );
    // Check the carrier equation and positive nappe independently of face_contains.
    let Surface::Cone(c) = model.geom().surface(model.face(face).surface) else {
        panic!()
    };
    let delta = witness.minus(&c.apex);
    let axial = delta.dot(&c.axis);
    let axial_sq = axial.square() / c.axis.norm_squared();
    assert_eq!(
        sign(
            &(delta.norm_squared() - axial_sq.clone() - c.slope.square() * axial_sq),
            &mut dec
        )
        .unwrap(),
        Sign::Zero
    );
    assert_ne!(sign(&axial, &mut dec).unwrap(), Sign::Negative);
    answer
}
#[test]
fn conical_exact_and_mpfr_enclosures_survive_pose_scale_axis_and_all_effective_senses() {
    let moved = Pose::from_rational(
        Rotation::from_quaternion(q(1), q(2), q(3), q(4)).unwrap(),
        RVector3::new(q(7), q(-11), q(13)),
    );
    for pose in [Pose::identity(), moved] {
        for scale in [Rational::from_ratio(1, 1000).unwrap(), q(1), q(1000)] {
            for reverse in [false, true] {
                for tip in [false, true] {
                    // Four compatible combinations and four complementary ones.
                    for fs in [Sense::Same, Sense::Reversed] {
                        for es in [Sense::Same, Sense::Reversed] {
                            for cs in [Sense::Same, Sense::Reversed] {
                                let (mut model, face) =
                                    fixture(&pose, &scale, reverse, [fs, es, cs], tip);
                                let retained = (es == cs) == (fs == Sense::Same);
                                if !retained {
                                    assert!(
                                        PreparedFace::recognize(
                                            &model,
                                            face,
                                            &mut Decider::default()
                                        )
                                        .unwrap()
                                        .is_none()
                                    );
                                    continue;
                                }
                                for (xyz, squared) in [
                                    ([4, 0, 2], 2),
                                    ([1, 0, 1], if tip { 0 } else { 2 }),
                                    ([8, 0, 8], 8),
                                    ([0, 0, 4], 8),
                                    ([0, 0, -3], if tip { 9 } else { 29 }),
                                    ([-2, 0, -2], if tip { 8 } else { 16 }),
                                    ([2, 0, 2], 0),
                                    ([6, 0, 6], 0),
                                    ([0, 0, 0], if tip { 0 } else { 8 }),
                                ] {
                                    let at = pose.apply_point(
                                        &RPoint3::new(
                                            &q(xyz[0]) * &scale,
                                            &q(xyz[1]) * &scale,
                                            &q(xyz[2] * if reverse { -1 } else { 1 }) * &scale,
                                        )
                                        .to_real(),
                                    );
                                    check(
                                        &mut model,
                                        face,
                                        at,
                                        independent_sqrt(squared).map(|x| &x * &scale),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn conical_irrational_azimuths_match_independent_directed_mpfr_enclosures() {
    let (mut model, face) = fixture(&Pose::identity(), &q(1), false, [Sense::Same; 3], false);
    for (at, oracle) in [
        (p(1, 1, 4), independent_radical(9, 4, 2, true)),
        (p(1, 1, 0), independent_radical(10, 4, 2, true)),
    ] {
        check(&mut model, face, at, oracle);
    }
}
#[test]
fn stock_tips_growing_and_shrinking_frusta_have_direct_admitted_distances() {
    for (base, top, height, point, squared) in [
        (6, 0, 6, p(0, 0, 8), 4),
        (6, 2, 4, p(8, 0, -2), 8),
        (2, 6, 4, p(8, 0, 6), 8),
    ] {
        let mut model = Model::new();
        let body =
            catallax_topology::cone(&mut model, q(base), q(top), q(height), &Pose::identity())
                .unwrap();
        let face = model
            .faces_of_body(body)
            .into_iter()
            .find(|f| {
                matches!(
                    model.geom().surface(model.face(*f).surface),
                    Surface::Cone(_)
                )
            })
            .unwrap();
        check(&mut model, face, point, independent_sqrt(squared));
    }
}
#[test]
fn conical_narrow_width_uses_the_same_complete_fallback_and_cell_cap() {
    let (mut model, face) = fixture(&Pose::identity(), &q(1), false, [Sense::Same; 3], false);
    let mut prepared = PreparedFace::new(&model, face, &mut Decider::default()).unwrap();
    let opts = Options {
        width: q(1).mul_pow2(-512),
        max_cells: 1,
    };
    let at = p(4, 0, 2);
    assert!(
        prepared
            .finite
            .as_ref()
            .unwrap()
            .point_distance(&at, &opts, &mut Decider::default())
            .unwrap()
            .is_none()
    );
    let a = prepared
        .point_to_face(&mut model, &at, &opts, &mut Decider::default())
        .unwrap();
    prepared.finite = None;
    let b = prepared
        .point_to_face(&mut model, &at, &opts, &mut Decider::default())
        .unwrap();
    assert_eq!(
        (a.lower, a.upper, a.cells, a.converged),
        (b.lower, b.upper, b.cells, b.converged)
    );
    // One cell cannot initialize the cone's complete multi-chart cover.
    assert_eq!(b.cells, 0);
    assert!(!b.converged);
}

#[test]
fn conical_collapsed_boundary_inventory_must_be_the_single_exact_apex() {
    let (mut model, face) = fixture(&Pose::identity(), &q(1), false, [Sense::Same; 3], true);
    let original = model.face(face).clone();
    let point = model.geom_mut().add_point(PointDef::Explicit(p(0, 0, 1)));
    let vertex = model.add_vertex(Vertex { point });
    let wrong = model.add_loop(Loop::Vertex(vertex));
    let mut bad = original.clone();
    bad.loops[1] = wrong;
    let face = model.add_face(bad);
    assert!(
        PreparedFace::recognize(&model, face, &mut Decider::default())
            .unwrap()
            .is_none()
    );
    let mut duplicate = original.clone();
    duplicate.loops.push(original.loops[1]);
    let face = model.add_face(duplicate);
    assert!(
        PreparedFace::recognize(&model, face, &mut Decider::default())
            .unwrap()
            .is_none()
    );
    let mut missing = original;
    missing.loops.pop();
    let face = model.add_face(missing);
    assert!(
        PreparedFace::recognize(&model, face, &mut Decider::default())
            .unwrap()
            .is_none()
    );
}

fn half_cone() -> (Model, FaceId) {
    use catallax_shapes::Line;
    let mut model = Model::new();
    let vector = |x, y, z| RVector3::new(q(x), q(y), q(z)).to_real();
    let surface = model.geom_mut().add_surface(Surface::Cone(Cone {
        apex: p(0, 0, 0),
        axis: vector(0, 0, 1),
        slope: Real::one(),
        ref_dir: vector(1, 0, 0),
    }));
    let vertices: Vec<_> = [p(0, 0, 0), p(-2, 0, 2), p(2, 0, 2)]
        .into_iter()
        .map(|at| {
            let point = model.geom_mut().add_point(PointDef::Explicit(at));
            model.add_vertex(Vertex { point })
        })
        .collect();
    let curves = [
        Curve::Line(Line {
            origin: p(0, 0, 0),
            dir: vector(-1, 0, 1),
        }),
        Curve::Circle(Circle {
            center: p(0, 0, 2),
            normal: vector(0, 0, 1),
            radius_sq: q(4).into(),
            ref_dir: vector(1, 0, 0),
        }),
        Curve::Line(Line {
            origin: p(2, 0, 2),
            dir: vector(-1, 0, -1),
        }),
    ];
    let coedges = curves
        .into_iter()
        .enumerate()
        .map(|(i, curve)| {
            let curve = model.geom_mut().add_curve(curve);
            let edge = model.add_edge(Edge {
                curve,
                sense: if i == 1 { Sense::Reversed } else { Sense::Same },
                start: Some(vertices[i]),
                end: Some(vertices[(i + 1) % 3]),
            });
            model.add_coedge(Coedge {
                edge,
                sense: Sense::Same,
            })
        })
        .collect();
    let lp = model.add_loop(Loop::Coedges(coedges));
    let face = model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: vec![lp],
    });
    (model, face)
}
#[test]
fn conical_excluded_angular_sector_cannot_supply_a_witness_and_keeps_fallback() {
    let (mut model, face) = half_cone();
    let mut dec = Decider::default();
    let mut prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    assert!(prepared.finite.is_none());
    assert!(!catallax_boolean::face_contains(&mut model, face, &p(0, -2, 2), &mut dec).unwrap());
    for (at, squared) in [(p(2, 0, 2), 0), (p(0, -2, 2), 6)] {
        let options = Options {
            width: controls().width,
            max_cells: 128,
        };
        let answer = prepared
            .point_to_face(&mut model, &at, &options, &mut dec)
            .unwrap();
        let oracle = independent_sqrt(squared);
        assert!(answer.lower <= oracle[0], "{answer:?}");
        assert!(answer.upper.as_ref().unwrap() >= &oracle[1], "{answer:?}");
        if squared > 0 {
            assert!(answer.lower > q(0), "{answer:?}");
        }
        let [_, witness] = answer.witnesses.as_ref().unwrap();
        assert!(catallax_boolean::face_contains(&mut model, face, witness, &mut dec).unwrap());
        assert_ne!(sign(&witness.y, &mut dec).unwrap(), Sign::Negative);
        prepared.finite = None;
        let fallback = distance::point_distance(
            &prepared.domain,
            &at,
            &options,
            &mut |p: &Point3, b: &[Ball; 3], d: &mut Decider| {
                prepared.classify(&mut model, p, b, d)
            },
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(
            (answer.lower, answer.upper, answer.cells, answer.converged),
            (
                fallback.lower,
                fallback.upper,
                fallback.cells,
                fallback.converged
            )
        );
    }
}

#[test]
fn conical_nonunit_algebraic_frame_and_slope_keep_the_euclidean_metric() {
    // Axis (1,2,2) has length 3. Reference (2,-1,0) has irrational length sqrt5.
    // This is independent of the cone's azimuthal parametrisation.
    let mut model = Model::new();
    let axis = RVector3::new(q(1), q(2), q(2)).to_real();
    let reference = RVector3::new(q(2), q(-1), q(0)).to_real();
    let surface = model.geom_mut().add_surface(Surface::Cone(Cone {
        apex: p(0, 0, 0),
        axis: axis.clone(),
        slope: Rational::from_ratio(3, 4).unwrap().into(),
        ref_dir: reference.clone(),
    }));
    let mut loops = Vec::new();
    for (height, sense) in [(4, Sense::Same), (8, Sense::Reversed)] {
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: p(height / 4, height / 2, height / 2), // physical height is 3*height/4
            normal: axis.clone(),
            radius_sq: Rational::from_ratio(81 * height * height, 256)
                .unwrap()
                .into(),
            ref_dir: reference.clone(),
        }));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        let coedge = model.add_coedge(Coedge { edge, sense });
        loops.push(model.add_loop(Loop::Coedges(vec![coedge])));
    }
    let face = model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops,
    });
    // Query on axis at height 6: min distance 6*k/sqrt(1+k²) = 18/5,
    // attained at height 96/25, strictly between heights 3 and 6.
    let expected = Rational::from_ratio(18, 5).unwrap();
    check(&mut model, face, p(2, 4, 4), [expected.clone(), expected]);
}
