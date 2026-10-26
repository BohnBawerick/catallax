use super::*;
use catallax_shapes::{Pose, RPoint3, Surface};
use catallax_topology::FaceId;
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(x: i64, y: i64, z: i64) -> Point3 {
    RPoint3::new(q(x), q(y), q(z)).into()
}
fn controls() -> Options {
    Options {
        width: Rational::from_ratio(1, 1000).unwrap(),
        max_cells: 4096,
    }
}
#[test]
fn cylinder_arc_trims_exclude_the_nearest_untrimmed_support_point() {
    let (mut model, body) = crate::fixtures::half_cylinder(&Pose::identity());
    let wall = model
        .faces_of_body(body)
        .into_iter()
        .find(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Cylinder(_)
            )
        })
        .unwrap();
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, wall, &mut dec).unwrap();
    let answer = prepared
        .point_to_face(&mut model, &p(0, -2, 1), &controls(), &mut dec)
        .unwrap();
    let want = Ball::from_rational(&q(8), 256).sqrt().unwrap();
    assert!(answer.lower <= want.lo(), "{answer:?}");
    assert!(answer.upper.as_ref().unwrap() >= &want.hi(), "{answer:?}");
    assert!(answer.lower > q(2), "{answer:?}");
    assert!(answer.upper.as_ref().unwrap() < &q(3), "{answer:?}");
    let [_, witness] = answer.witnesses.unwrap();
    assert!(catallax_boolean::face_contains(&mut model, wall, &witness, &mut dec).unwrap());
}
#[test]
fn cap_pair_minimum_distance_and_cylinder_height_boundary_are_certified() {
    let mut model = Model::new();
    let body = catallax_topology::cylinder(&mut model, q(1), q(2), &Pose::identity()).unwrap();
    let faces = model.faces_of_body(body);
    let wall = faces
        .iter()
        .copied()
        .find(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Cylinder(_)
            )
        })
        .unwrap();
    let caps: Vec<_> = faces
        .into_iter()
        .filter(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Plane(_)
            )
        })
        .collect();
    let mut dec = Decider::default();
    let wall = PreparedFace::new(&model, wall, &mut dec).unwrap();
    let answer = wall
        .point_to_face(&mut model, &p(1, 0, 3), &controls(), &mut dec)
        .unwrap();
    assert!(
        answer.lower <= q(1) && answer.lower >= Rational::parse("0.999").unwrap(),
        "{answer:?}"
    );
    assert!(
        answer.upper.as_ref().unwrap() >= &q(1)
            && answer.upper.as_ref().unwrap() <= &Rational::parse("1.001").unwrap(),
        "{answer:?}"
    );
    let a = PreparedFace::new(&model, caps[0], &mut dec).unwrap();
    let b = PreparedFace::new(&model, caps[1], &mut dec).unwrap();
    let answer = a
        .face_to_face(&b, &mut model, &controls(), &mut dec)
        .unwrap();
    assert_eq!(answer.lower, q(2));
    assert_eq!(answer.upper, Some(q(2)));
}
#[test]
fn spherical_cap_does_not_use_the_removed_south_pole() {
    let (mut model, body) = crate::fixtures::cap(5, 3);
    let face = model
        .faces_of_body(body)
        .into_iter()
        .find(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Sphere(_)
            )
        })
        .unwrap();
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    let answer = prepared
        .point_to_face(&mut model, &p(0, 0, -5), &controls(), &mut dec)
        .unwrap();
    assert_eq!(
        answer.cells, 0,
        "complete spherical caps should avoid subdivision: {answer:?}"
    );
    // The boundary circle has radius four and z=3, so the minimum is sqrt(80).
    let want = Ball::from_rational(&q(80), 256).sqrt().unwrap();
    assert!(answer.lower <= want.lo(), "{answer:?}");
    assert!(answer.upper.as_ref().unwrap() >= &want.hi(), "{answer:?}");
    assert!(answer.lower > q(8), "{answer:?}");
}

#[test]
fn elliptic_cylinder_preserves_euclidean_distance_and_oblique_circle_trims() {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Circle, Curve, EllipticCylinder, GeomKey, RVector3};
    use catallax_topology::{Coedge, Edge, Face, Loop, Sense};

    let mut model = Model::new();
    let surface = model
        .geom_mut()
        .add_surface(Surface::EllipticCylinder(EllipticCylinder {
            origin: p(0, 0, 0),
            axis: RVector3::new(q(3), q(0), q(4)).into(),
            major_dir: RVector3::new(q(0), q(1), q(0)).into(),
            major_sq: q(1).into(),
            minor_sq: Rational::from_ratio(16, 25).unwrap().into(),
        }));
    // The unit circles at (0,0,0) and (3,0,4), both parallel to xy,
    // lie on this carrier. They bound an oblique circular sweep.
    let mut loops = Vec::new();
    for (center, sense) in [(p(0, 0, 0), Sense::Same), (p(3, 0, 4), Sense::Reversed)] {
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center,
            normal: RVector3::new(q(0), q(0), q(1)).into(),
            radius_sq: q(1).into(),
            ref_dir: RVector3::new(q(1), q(0), q(0)).into(),
        }));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Curve(curve),
            b: GeomKey::Surface(surface),
            at: None,
            origin: FactOrigin::Construction("measure.distance-elliptic-cylinder"),
        });
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
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    for (xyz, expected) in [
        (["29/10", "0", "17/10"], "1/2"),
        (["23/5", "0", "24/5"], "1"),
    ] {
        let point = RPoint3::new(
            Rational::parse(xyz[0]).unwrap(),
            Rational::parse(xyz[1]).unwrap(),
            Rational::parse(xyz[2]).unwrap(),
        )
        .to_real();
        let answer = prepared
            .point_to_face(&mut model, &point, &controls(), &mut dec)
            .unwrap();
        let expected = Rational::parse(expected).unwrap();
        let upper = answer.upper.as_ref().unwrap();
        assert!(answer.lower <= expected && upper >= &expected, "{answer:?}");
        assert!(
            answer.converged && upper - &answer.lower <= controls().width,
            "{answer:?}"
        );
        let [_, witness] = answer.witnesses.unwrap();
        assert!(catallax_boolean::face_contains(&mut model, face, &witness, &mut dec).unwrap());
    }
}

fn placed() -> Pose {
    use catallax_shapes::{RVector3, Rotation};
    Pose::from_rational(
        Rotation::from_quaternion(q(1), q(1), q(1), q(1)).unwrap(),
        RVector3::new(q(1), q(-2), q(3)),
    )
}
fn world(pose: &Pose, local: [i64; 3]) -> Point3 {
    pose.apply_point(&p(local[0], local[1], local[2]))
}
fn solid_faces(pose: &Pose) -> (Model, FaceId, [FaceId; 2]) {
    let mut model = Model::new();
    let body = catallax_topology::cylinder(&mut model, q(2), q(3), pose).unwrap();
    let faces = model.faces_of_body(body);
    let is = |model: &Model, f: FaceId, wall: bool| {
        matches!(
            model.geom().surface(model.face(f).surface),
            Surface::Cylinder(_)
        ) == wall
    };
    let wall = faces
        .iter()
        .copied()
        .find(|f| is(&model, *f, true))
        .unwrap();
    let caps: Vec<_> = faces
        .into_iter()
        .filter(|f| is(&model, *f, false))
        .collect();
    (model, wall, [caps[0], caps[1]])
}
fn exact(
    model: &mut Model,
    face: FaceId,
    point: &Point3,
    want: &Rational,
    dec: &mut Decider,
) -> Bounds {
    let prepared = PreparedFace::new(model, face, dec).unwrap();
    let answer = prepared
        .point_to_face(model, point, &controls(), dec)
        .unwrap();
    let upper = answer.upper.clone().unwrap();
    assert!(answer.lower <= *want && upper >= *want, "{answer:?}");
    assert!(answer.converged && answer.cells == 0, "{answer:?}");
    assert!(
        &upper - &answer.lower <= Rational::from_ratio(1, 1_000_000_000).unwrap(),
        "{answer:?}"
    );
    let [origin, witness] = answer.witnesses.clone().unwrap();
    assert_eq!(origin.as_rational(), point.as_rational());
    assert!(catallax_boolean::face_contains(model, face, &witness, dec).unwrap());
    answer
}
#[test]
fn placed_cylinder_band_matches_independent_distances_for_every_nearest_case() {
    for pose in [Pose::identity(), placed()] {
        let (mut model, wall, _) = solid_faces(&pose);
        let mut dec = Decider::default();
        // Local points and exact distances to the radius 2, height 3 wall.
        for (local, want) in [
            ([5, 0, 1], q(3)),
            ([0, 0, 1], q(2)),
            ([1, 0, 1], q(1)),
            ([5, 0, 7], q(5)),
            ([5, 0, -4], q(5)),
            ([2, 0, 7], q(4)),
            ([0, 0, 3], q(2)),
            ([2, 0, 3], q(0)),
            ([3, 4, 0], q(3)),
        ] {
            exact(&mut model, wall, &world(&pose, local), &want, &mut dec);
        }
        // On the axis beyond the end the nearest set is a rim circle at sqrt(2^2+4^2).
        let axis = world(&pose, [0, 0, 7]);
        let answer = PreparedFace::new(&model, wall, &mut dec)
            .unwrap()
            .point_to_face(&mut model, &axis, &controls(), &mut dec)
            .unwrap();
        let want = Ball::from_rational(&q(20), 256).sqrt().unwrap();
        assert!(answer.lower <= want.lo() && answer.upper.unwrap() >= want.hi());
    }
}
#[test]
fn finite_band_is_not_its_infinite_support() {
    let pose = placed();
    let (mut model, wall, _) = solid_faces(&pose);
    let mut dec = Decider::default();
    // The point lies on the cylinder's infinite support but 4 above the band.
    let answer = exact(&mut model, wall, &world(&pose, [2, 0, 7]), &q(4), &mut dec);
    assert!(answer.lower > q(3));
}
#[test]
fn placed_caps_match_independent_distances_including_the_rim_tie() {
    for pose in [Pose::identity(), placed()] {
        let (mut model, _, caps) = solid_faces(&pose);
        let mut dec = Decider::default();
        // Faces are built bottom (z = 0) then top (z = 3); the zero-distance
        // rows below fail if that order ever changes.
        for (cap, base) in caps.into_iter().zip([0, 3]) {
            for (local, want) in [
                ([1, 1, base + 4], q(4)),
                ([0, 0, base - 2], q(2)),
                ([5, 0, base + 4], q(5)),
                ([2, 0, base + 1], q(1)),
                ([2, 0, base], q(0)),
                ([3, 4, base], q(3)),
                ([1, 0, base], q(0)),
            ] {
                exact(&mut model, cap, &world(&pose, local), &want, &mut dec);
            }
        }
    }
}

fn sqrt2() -> Ball {
    Ball::from_rational(&q(2), 256).sqrt().unwrap()
}
/// Checks one irrational nearest-point case on the radius 2, height 3 cylinder.
fn irrational_case(wall: bool, at: Point3, expected: &Ball, witness: [Ball; 3]) {
    let (mut model, side, caps) = solid_faces(&Pose::identity());
    let face = if wall { side } else { caps[1] };
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    let answer = prepared
        .point_to_face(&mut model, &at, &controls(), &mut dec)
        .unwrap();
    let upper = answer.upper.clone().unwrap();
    let tight = Rational::from_ratio(1, 1_000_000_000).unwrap();
    assert!(answer.converged && answer.cells == 0, "{answer:?}");
    assert!(
        answer.lower <= expected.hi() && upper >= expected.lo(),
        "{answer:?}"
    );
    assert!(&upper - &answer.lower <= tight, "{answer:?}");
    assert!(&expected.hi() - &answer.lower <= tight, "{answer:?}");
    assert!(&upper - &expected.lo() <= tight, "{answer:?}");
    let [_, found] = answer.witnesses.unwrap();
    for (coordinate, want) in [&found.x, &found.y, &found.z].into_iter().zip(&witness) {
        let got = dec.enclose_budgeted_at(coordinate, 128).unwrap();
        assert!(got.lo() <= want.hi() && want.lo() <= got.hi(), "{found:?}");
        assert!(&got.hi() - &got.lo() <= tight, "{found:?}");
    }
    assert!(catallax_boolean::face_contains(&mut model, face, &found, &mut dec).unwrap());
}
#[test]
fn wall_point_with_an_irrational_witness_is_enclosed_without_subdivision() {
    let r = sqrt2();
    let one = Ball::from_rational(&q(1), 256);
    // The point (1, 1, 1) is inside the wall, so the nearest point is
    // (sqrt2, sqrt2, 1) at distance 2 - sqrt2.
    let expected = Ball::from_rational(&q(2), 256) - r.clone();
    irrational_case(true, p(1, 1, 1), &expected, [r.clone(), r, one]);
}
#[test]
fn cap_point_with_an_irrational_rim_witness_is_enclosed_without_subdivision() {
    let r = sqrt2();
    // The point (3, 3, 7) projects outside the top cap, so the nearest point is
    // the rim point (sqrt2, sqrt2, 3) at distance sqrt(38 - 12 sqrt2).
    let squared = Ball::from_rational(&q(38), 256) - Ball::from_rational(&q(12), 256) * r.clone();
    let expected = squared.sqrt().unwrap();
    let three = Ball::from_rational(&q(3), 256);
    irrational_case(false, p(3, 3, 7), &expected, [r.clone(), r, three]);
}

// Public synthetic finite spherical faces. Trim axes deliberately differ from
// the carrier frame, so recognition cannot infer a cap from its parametrisation.
fn spherical_face(
    pose: &Pose,
    scale: &Rational,
    face_sense: Sense,
    edge_sense: Sense,
    rims: &[(i64, Sense)],
) -> (Model, FaceId) {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Circle, GeomKey, RVector3, Sphere};
    use catallax_topology::{Coedge, Edge, Face};
    let mut model = Model::new();
    let mut dec = Decider::default();
    let vector = |x, y, z| pose.apply_vector(&RVector3::new(q(x), q(y), q(z)).to_real());
    let surface = model.geom_mut().add_surface(Surface::Sphere(
        Sphere::new(
            pose.apply_point(&p(0, 0, 0)),
            (&q(5) * scale).into(),
            vector(2, 0, 0),
            vector(0, 3, 0),
            &mut dec,
        )
        .unwrap(),
    ));
    let loops = rims
        .iter()
        .map(|(z, sense)| {
            let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                center: pose.apply_point(&RPoint3::new(q(0), q(0), &q(*z) * scale).to_real()),
                normal: vector(0, 0, 2),
                radius_sq: (&q(25 - z * z) * &(scale * scale)).into(),
                ref_dir: vector(3, 0, 0),
            }));
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(curve),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("public-spherical-distance-tests"),
            });
            let edge = model.add_edge(Edge {
                curve,
                sense: edge_sense,
                start: None,
                end: None,
            });
            let coedge = model.add_coedge(Coedge {
                edge,
                sense: *sense,
            });
            model.add_loop(Loop::Coedges(vec![coedge]))
        })
        .collect();
    let face = model.add_face(Face {
        surface,
        sense: face_sense,
        loops,
    });
    (model, face)
}

// MPFR directed rounding is independent of kernel Real/Ball and sphere charts.
fn independent_sqrt(n: i64) -> [Rational; 2] {
    use rug::{Float, float::Round};
    [Round::Down, Round::Up].map(|round| {
        let mut f = Float::with_val(256, n);
        f.sqrt_round(round);
        Rational::parse(&f.to_rational().unwrap().to_string()).unwrap()
    })
}
fn check_spherical(
    model: &mut Model,
    face: FaceId,
    at: &Point3,
    expected: [Rational; 2],
) -> Bounds {
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(model, face, &mut dec).unwrap();
    assert!(prepared.finite.is_some());
    let answer = prepared
        .point_to_face(model, at, &controls(), &mut dec)
        .unwrap();
    let upper = answer.upper.as_ref().unwrap();
    assert!(answer.converged && answer.cells == 0, "{answer:?}");
    assert!(
        answer.lower <= expected[0] && *upper >= expected[1],
        "{answer:?}, {expected:?}"
    );
    assert!(upper - &answer.lower <= controls().width, "{answer:?}");
    let [_, witness] = answer.witnesses.as_ref().unwrap();
    assert!(
        catallax_boolean::face_contains(model, face, witness, &mut dec).unwrap(),
        "{answer:?}"
    );
    answer
}
#[test]
fn spherical_caps_and_bands_match_independent_oracles_under_pose_and_scale() {
    for pose in [Pose::identity(), placed()] {
        for scale in [Rational::from_ratio(1, 1000).unwrap(), q(1), q(1000)] {
            for (rims, cases) in [
                (
                    vec![(3, Sense::Same)],
                    vec![
                        ([0, 0, -5], 80),
                        ([0, 0, 7], 4),
                        ([0, 0, 0], 25),
                        ([3, 0, 4], 0),
                        ([6, 0, 8], 25),
                        ([4, 0, 3], 0),
                    ],
                ),
                (
                    vec![(-3, Sense::Same), (3, Sense::Reversed)],
                    vec![
                        ([6, 0, 8], 29),
                        ([7, 0, -8], 34),
                        ([0, 0, 0], 25),
                        ([0, 0, 7], 32),
                        ([0, 0, -7], 32),
                        ([10, 0, 0], 25),
                        ([4, 0, 3], 0),
                    ],
                ),
                (
                    vec![(0, Sense::Same), (3, Sense::Reversed)],
                    vec![([7, 0, -8], 68), ([6, 0, 8], 29), ([5, 0, 0], 0)],
                ),
                (
                    vec![(3, Sense::Same), (4, Sense::Reversed)],
                    vec![([0, 0, 7], 18), ([0, 0, -5], 80), ([3, 0, 4], 0)],
                ),
            ] {
                let (mut model, face) =
                    spherical_face(&pose, &scale, Sense::Same, Sense::Same, &rims);
                for (xyz, squared) in cases {
                    let at = pose.apply_point(
                        &RPoint3::new(
                            &q(xyz[0]) * &scale,
                            &q(xyz[1]) * &scale,
                            &q(xyz[2]) * &scale,
                        )
                        .to_real(),
                    );
                    check_spherical(
                        &mut model,
                        face,
                        &at,
                        independent_sqrt(squared).map(|x| &x * &scale),
                    );
                }
            }
        }
    }
}
#[test]
fn complementary_caps_and_effective_senses_keep_opposite_regions() {
    // Inverting any one sense selects the complement. Inverting two preserves it.
    for (face_sense, edge_sense, coedge_sense, north) in [
        (Sense::Same, Sense::Same, Sense::Same, true),
        (Sense::Same, Sense::Same, Sense::Reversed, false),
        (Sense::Same, Sense::Reversed, Sense::Same, false),
        (Sense::Reversed, Sense::Same, Sense::Same, false),
        (Sense::Reversed, Sense::Reversed, Sense::Same, true),
        (Sense::Reversed, Sense::Same, Sense::Reversed, true),
        (Sense::Same, Sense::Reversed, Sense::Reversed, true),
        (Sense::Reversed, Sense::Reversed, Sense::Reversed, false),
    ] {
        let (mut model, face) = spherical_face(
            &Pose::identity(),
            &q(1),
            face_sense,
            edge_sense,
            &[(3, coedge_sense)],
        );
        check_spherical(
            &mut model,
            face,
            &p(0, 0, -5),
            independent_sqrt(if north { 80 } else { 0 }),
        );
        check_spherical(
            &mut model,
            face,
            &p(0, 0, 7),
            independent_sqrt(if north { 4 } else { 32 }),
        );
    }
}
#[test]
fn spherical_arc_trims_keep_the_complete_search_and_admitted_witnesses() {
    let (mut model, body) = crate::fixtures::octant(&Pose::identity());
    let face = model
        .faces_of_body(body)
        .into_iter()
        .find(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Sphere(_)
            )
        })
        .unwrap();
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    assert!(prepared.finite.is_none());
    let at = p(5, 0, 0);
    let answer = prepared
        .point_to_face(&mut model, &at, &controls(), &mut dec)
        .unwrap();
    let fallback = distance::point_distance(
        &prepared.domain,
        &at,
        &controls(),
        &mut |p: &Point3, b: &[Ball; 3], d: &mut Decider| prepared.classify(&mut model, p, b, d),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(answer.lower, fallback.lower);
    assert_eq!(answer.upper, fallback.upper);
    assert_eq!(answer.cells, fallback.cells);
    assert!(answer.converged);
    // The radius-3 octant admits its positive-x vertex, at distance 2.
    assert!(answer.lower <= q(2) && answer.upper.as_ref().unwrap() >= &q(2));
    let [_, witness] = answer.witnesses.unwrap();
    assert!(catallax_boolean::face_contains(&mut model, face, &witness, &mut dec).unwrap());
    assert!(!catallax_boolean::face_contains(&mut model, face, &p(-3, 0, 0), &mut dec).unwrap());
}
#[test]
fn insufficient_direct_width_retains_the_original_fallback_controls() {
    let (mut model, face) = spherical_face(
        &Pose::identity(),
        &q(1),
        Sense::Same,
        Sense::Same,
        &[(3, Sense::Same)],
    );
    let mut dec = Decider::default();
    let mut prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    let options = Options {
        width: Rational::from_ratio(1, 10).unwrap().mul_pow2(-512),
        max_cells: 1,
    };
    assert!(
        prepared
            .finite
            .as_ref()
            .unwrap()
            .point_distance(&p(0, 0, -5), &options, &mut dec)
            .unwrap()
            .is_none()
    );
    let direct_fallback = prepared
        .point_to_face(&mut model, &p(0, 0, -5), &options, &mut dec)
        .unwrap();
    prepared.finite = None;
    let fallback = prepared
        .point_to_face(&mut model, &p(0, 0, -5), &options, &mut Decider::default())
        .unwrap();
    assert_eq!(direct_fallback.lower, fallback.lower);
    assert_eq!(direct_fallback.upper, fallback.upper);
    assert_eq!(direct_fallback.cells, fallback.cells);
    assert_eq!(direct_fallback.converged, fallback.converged);
    assert!(!fallback.converged);
}

fn independent_radical(a: i32, b: i32, c: i32, outer_root: bool) -> [Rational; 2] {
    use rug::{Float, float::Round};
    [(Round::Down, Round::Up), (Round::Up, Round::Down)].map(|(outward, inward)| {
        let mut root = Float::with_val(256, c);
        root.sqrt_round(inward);
        let product = Float::with_val_round(256, &root * b, inward).0;
        let mut value = Float::with_val_round(256, a - &product, outward).0;
        if outer_root {
            value.sqrt_round(outward);
        }
        Rational::parse(&value.to_rational().unwrap().to_string()).unwrap()
    })
}
#[test]
fn irrational_spherical_witnesses_match_independent_mpfr_oracles() {
    use catallax_shapes::{RVector3, Rotation};
    let pose = Pose::from_rational(
        Rotation::from_quaternion(q(1), q(2), q(3), q(4)).unwrap(),
        RVector3::new(q(7), q(-11), q(13)),
    );
    for pose in [Pose::identity(), pose] {
        for (rims, xyz, oracle) in [
            (
                vec![(3, Sense::Same)],
                [1, 1, 4],
                independent_radical(5, 1, 18, false),
            ),
            (
                vec![(3, Sense::Same)],
                [1, 1, -1],
                independent_radical(34, 8, 2, true),
            ),
            (
                vec![(-3, Sense::Same), (3, Sense::Reversed)],
                [1, 1, 0],
                independent_radical(5, 1, 2, false),
            ),
        ] {
            let (mut model, face) = spherical_face(&pose, &q(1), Sense::Same, Sense::Same, &rims);
            check_spherical(&mut model, face, &world(&pose, xyz), oracle);
        }
    }
}

#[path = "conical_tests.rs"]
mod conical;
