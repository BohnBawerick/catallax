//! Complete conic boolean paths: closed shells, exact-reference properties,
//! rational-pose covariance, analytic STEP and boolean volume conservation.
use catallax_boolean::{BoolOp, boolean};
use catallax_measure::{Accuracy, Enclosure, MassProperties, mass_properties};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation, Surface};
use catallax_step::{StepOptions, write_step};
use catallax_topology::{self as topology, BodyId, Loop, Model, Sense};
use std::collections::{BTreeMap, BTreeSet};

const BITS: u32 = 256;
fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&r(n), BITS)
}
fn fraction(n: i64, d: i64) -> Ball {
    Ball::from_rational(&Rational::from_ratio(n, d).unwrap(), BITS)
}
fn divide(a: &Ball, b: &Ball) -> Ball {
    a.checked_div(b).unwrap()
}
fn shifted(pose: &Pose, xyz: [i64; 3]) -> Pose {
    Pose::new(
        Rotation::identity(),
        RVector3::new(r(xyz[0]), r(xyz[1]), r(xyz[2])),
    )
    .then(pose)
}
fn poses() -> [Pose; 3] {
    fn next(state: &mut u64) -> i64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 32) % 11) as i64 - 5
    }
    let mut state = 0x434f4e45544f5255;
    std::array::from_fn(|i| {
        if i == 0 {
            return Pose::identity();
        }
        let rotation = Rotation::from_quaternion(
            r(next(&mut state).abs() + 1),
            r(next(&mut state)),
            r(next(&mut state)),
            r(next(&mut state)),
        )
        .unwrap();
        Pose::new(
            rotation,
            RVector3::new(
                r(next(&mut state)),
                r(next(&mut state)),
                r(next(&mut state)),
            ),
        )
    })
}
fn world(pose: &Pose, point: &[Ball; 3]) -> [Ball; 3] {
    let origin = pose.apply_point(&RPoint3::origin());
    let axes = [
        RVector3::new(r(1), r(0), r(0)),
        RVector3::new(r(0), r(1), r(0)),
        RVector3::new(r(0), r(0), r(1)),
    ]
    .map(|v| pose.apply_vector(&v));
    let translations = [&origin.x, &origin.y, &origin.z];
    std::array::from_fn(|i| {
        let components = axes.each_ref().map(|v| [&v.x, &v.y, &v.z][i]);
        point
            .iter()
            .zip(components)
            .fold(Ball::from_rational(translations[i], BITS), |sum, (p, c)| {
                &sum + &(p * &Ball::from_rational(c, BITS))
            })
    })
}
fn agrees(actual: &Enclosure, expected: &Ball) {
    // Reference arithmetic is enclosed too; an exact answer can be narrower.
    assert!(
        actual.lo <= expected.hi() && actual.hi >= expected.lo(),
        "{actual:?} is disjoint from reference [{}, {}]",
        expected.lo(),
        expected.hi()
    );
}
fn check_shells(model: &Model, body: BodyId, expected: &[i64]) {
    let faults = model.check(body);
    assert!(faults.is_empty(), "{faults:?}");
    let mut euler = Vec::new();
    for &lump in &model.body(body).lumps {
        for &shell in &model.lump(lump).shells {
            let mut edges = BTreeSet::new();
            let mut vertices = BTreeSet::new();
            let mut uses: BTreeMap<_, (usize, i64)> = BTreeMap::new();
            let mut faces = 0i64;
            for &face in &model.shell(shell).faces {
                let face = model.face(face);
                faces += if face.loops.is_empty()
                    && matches!(model.geom().surface(face.surface), Surface::Torus(_))
                {
                    0
                } else {
                    2 - face.loops.len() as i64
                };
                for &boundary in &face.loops {
                    match model.loop_(boundary) {
                        Loop::Vertex(vertex) => {
                            vertices.insert(*vertex);
                        }
                        Loop::Coedges(coedges) => {
                            for &coedge in coedges {
                                let coedge = model.coedge(coedge);
                                edges.insert(coedge.edge);
                                let entry = uses.entry(coedge.edge).or_default();
                                entry.0 += 1;
                                entry.1 += if coedge.sense == Sense::Same { 1 } else { -1 };
                                let edge = model.edge(coedge.edge);
                                vertices.extend([edge.start, edge.end].into_iter().flatten());
                            }
                        }
                    }
                }
            }
            for (edge, use_count) in uses {
                assert_eq!(use_count, (2, 0), "open or misoriented edge {edge:?}");
            }
            let intervals = edges
                .iter()
                .filter(|edge| model.edge(**edge).start.is_some())
                .count();
            euler.push(vertices.len() as i64 - intervals as i64 + faces);
        }
    }
    euler.sort_unstable();
    assert_eq!(euler, expected);
}
fn measured(model: &Model, body: BodyId, chi: &[i64], dec: &mut Decider) -> MassProperties {
    check_shells(model, body, chi);
    let accuracy = Accuracy::default();
    let properties = mass_properties(model, body, &accuracy, dec).unwrap();
    assert!(properties.marks.is_empty(), "{:?}", properties.marks);
    for enclosure in [&properties.value.volume, &properties.value.area]
        .into_iter()
        .chain(properties.value.centroid.iter().flatten())
    {
        let magnitude = enclosure.lo.abs().min(enclosure.hi.abs());
        let goal = (&accuracy.relative * &magnitude).max(accuracy.absolute.clone());
        assert!(&enclosure.hi - &enclosure.lo <= goal, "{enclosure:?}");
    }
    let mut bytes = Vec::new();
    let report = write_step(model, &[body], &StepOptions::default(), dec, &mut bytes).unwrap();
    assert!(report.marks.is_empty(), "{:?}", report.marks);
    assert_eq!(report.value.approximated_curves, 0);
    assert!(report.value.max_deviation < q("1e-8"));
    properties.value
}
struct Expected {
    volume: Ball,
    area: Ball,
    centroid: [Ball; 3],
    chi: &'static [i64],
}
fn properties_match(actual: &MassProperties, expected: &Expected, pose: &Pose) {
    agrees(&actual.volume, &expected.volume);
    agrees(&actual.area, &expected.area);
    for (actual, expected) in actual
        .centroid
        .as_ref()
        .unwrap()
        .iter()
        .zip(world(pose, &expected.centroid))
    {
        agrees(actual, &expected);
    }
}

#[derive(Clone, Copy, Debug)]
enum ClosedCase {
    ConeCylinder,
    ConeHalf,
    TorusUpper,
    TorusVertical,
    SphereTube,
}
fn closed_operands(case: ClosedCase, model: &mut Model, pose: &Pose) -> (BodyId, BodyId) {
    match case {
        ClosedCase::ConeCylinder => (
            topology::cone(model, r(2), r(0), r(2), pose).unwrap(),
            topology::cylinder(model, r(1), r(4), &shifted(pose, [0, 0, -1])).unwrap(),
        ),
        ClosedCase::ConeHalf => (
            topology::cone(model, r(2), r(0), r(2), pose).unwrap(),
            topology::block(model, [r(4), r(6), r(4)], &shifted(pose, [0, -3, -1])).unwrap(),
        ),
        ClosedCase::TorusUpper => (
            topology::torus(model, r(3), r(1), pose).unwrap(),
            topology::block(model, [r(10), r(10), r(4)], &shifted(pose, [-5, -5, 0])).unwrap(),
        ),
        ClosedCase::TorusVertical => (
            topology::torus(model, r(3), r(1), pose).unwrap(),
            topology::block(model, [r(5), r(10), r(4)], &shifted(pose, [0, -5, -2])).unwrap(),
        ),
        ClosedCase::SphereTube => (
            topology::torus(model, r(3), r(1), pose).unwrap(),
            topology::sphere(model, q("1/2"), &shifted(pose, [3, 0, 0])).unwrap(),
        ),
    }
}
fn closed_expected(case: ClosedCase, op: BoolOp) -> Expected {
    let pi = Ball::pi(BITS);
    let pi2 = pi.square();
    let sqrt2 = b(2).sqrt().unwrap();
    let zero = || [b(0), b(0), b(0)];
    match case {
        ClosedCase::ConeCylinder => {
            let (volume, factor, z, chi) = match op {
                BoolOp::Unite => (
                    fraction(16, 3),
                    &b(11) + &(&b(3) * &sqrt2),
                    fraction(53, 64),
                    &[2][..],
                ),
                BoolOp::Subtract => (
                    fraction(4, 3),
                    &b(5) + &(&b(3) * &sqrt2),
                    fraction(5, 16),
                    &[0][..],
                ),
                BoolOp::Intersect => (fraction(4, 3), &b(3) + &sqrt2, fraction(11, 16), &[2][..]),
                _ => unreachable!(),
            };
            Expected {
                volume: &volume * &pi,
                area: &factor * &pi,
                centroid: [b(0), b(0), z],
                chi,
            }
        }
        ClosedCase::ConeHalf => {
            let half_volume = &fraction(4, 3) * &pi;
            let curved = &(&b(2) * &sqrt2) * &pi;
            let area = &(&b(4) + &(&b(2) * &pi)) + &curved;
            match op {
                BoolOp::Unite => {
                    let volume = &b(96) + &half_volume;
                    Expected {
                        area: &area + &b(120),
                        centroid: [
                            divide(&fraction(568, 3), &volume),
                            b(0),
                            divide(&(&b(96) + &(&fraction(2, 3) * &pi)), &volume),
                        ],
                        volume,
                        chi: &[2],
                    }
                }
                BoolOp::Subtract | BoolOp::Intersect => Expected {
                    volume: half_volume,
                    area,
                    centroid: [
                        if op == BoolOp::Subtract {
                            divide(&b(-2), &pi)
                        } else {
                            divide(&b(2), &pi)
                        },
                        b(0),
                        fraction(1, 2),
                    ],
                    chi: &[2],
                },
                _ => unreachable!(),
            }
        }
        ClosedCase::TorusUpper | ClosedCase::TorusVertical => {
            let vertical = matches!(case, ClosedCase::TorusVertical);
            let half_volume = &b(3) * &pi2;
            let curved = &b(6) * &pi2;
            let caps = &b(if vertical { 2 } else { 12 }) * &pi;
            match op {
                BoolOp::Unite => {
                    let volume = &b(if vertical { 200 } else { 400 }) + &half_volume;
                    let mut centroid = zero();
                    centroid[if vertical { 0 } else { 2 }] = if vertical {
                        divide(&(&b(500) - &(&fraction(37, 2) * &pi)), &volume)
                    } else {
                        divide(&(&b(800) - &(&b(4) * &pi)), &volume)
                    };
                    Expected {
                        volume,
                        area: &(&b(if vertical { 220 } else { 360 }) + &curved) - &caps,
                        centroid,
                        chi: if vertical { &[0] } else { &[2] },
                    }
                }
                BoolOp::Subtract | BoolOp::Intersect => {
                    let mut centroid = zero();
                    let value = if vertical {
                        divide(&fraction(37, 6), &pi)
                    } else {
                        divide(&fraction(4, 3), &pi)
                    };
                    centroid[if vertical { 0 } else { 2 }] = if op == BoolOp::Subtract {
                        -value
                    } else {
                        value
                    };
                    Expected {
                        volume: half_volume,
                        area: &curved + &caps,
                        centroid,
                        chi: if vertical { &[2] } else { &[0] },
                    }
                }
                _ => unreachable!(),
            }
        }
        ClosedCase::SphereTube => match op {
            BoolOp::Unite => Expected {
                volume: &b(6) * &pi2,
                area: &b(12) * &pi2,
                centroid: zero(),
                chi: &[0],
            },
            BoolOp::Subtract => {
                let volume = &(&b(6) * &pi2) - &divide(&pi, &b(6));
                Expected {
                    centroid: [divide(&(-&pi), &(&b(2) * &volume)), b(0), b(0)],
                    volume,
                    area: &(&b(12) * &pi2) + &pi,
                    chi: &[0, 2],
                }
            }
            BoolOp::Intersect => Expected {
                volume: divide(&pi, &b(6)),
                area: pi,
                centroid: [b(3), b(0), b(0)],
                chi: &[2],
            },
            _ => unreachable!(),
        },
    }
}

#[test]
fn closed_form_conic_booleans_preserve_shells_and_pose_covariance() {
    for pose in poses() {
        for case in [
            ClosedCase::ConeCylinder,
            ClosedCase::ConeHalf,
            ClosedCase::TorusUpper,
            ClosedCase::TorusVertical,
            ClosedCase::SphereTube,
        ] {
            let mut model = Model::new();
            let (a, b) = closed_operands(case, &mut model, &pose);
            for op in [BoolOp::Unite, BoolOp::Subtract, BoolOp::Intersect] {
                let mut dec = Decider::default();
                let answer = boolean(&mut model, a, b, op, &mut dec)
                    .unwrap_or_else(|e| panic!("{case:?} {op:?}: {e:?}"));
                assert!(answer.marks.is_empty());
                let expected = closed_expected(case, op);
                let actual = measured(&model, answer.body, expected.chi, &mut dec);
                properties_match(&actual, &expected, &pose);
            }
        }
    }
}

#[test]
fn countersink_with_cylindrical_bore_has_exact_mass_and_one_closed_handle() {
    for pose in poses() {
        let mut model = Model::new();
        let block = topology::block(&mut model, [r(10), r(10), r(10)], &pose).unwrap();
        let bore =
            topology::cylinder(&mut model, r(1), r(12), &shifted(&pose, [5, 5, -1])).unwrap();
        let sink =
            topology::cone(&mut model, r(1), r(3), r(2), &shifted(&pose, [5, 5, 8])).unwrap();
        let mut dec = Decider::default();
        let drilled = boolean(&mut model, block, bore, BoolOp::Subtract, &mut dec)
            .unwrap()
            .body;
        let answer = boolean(&mut model, drilled, sink, BoolOp::Subtract, &mut dec).unwrap();
        assert!(answer.marks.is_empty());
        let pi = Ball::pi(BITS);
        let volume = &b(1000) - &(&fraction(50, 3) * &pi);
        let root = b(2).sqrt().unwrap();
        let area = &(&b(600) + &(&b(6) * &pi)) + &(&(&b(8) * &root) * &pi);
        let z = divide(&(&b(5000) - &(&fraction(338, 3) * &pi)), &volume);
        let expected = Expected {
            volume,
            area,
            centroid: [b(5), b(5), z],
            chi: &[0],
        };
        properties_match(
            &measured(&model, answer.body, expected.chi, &mut dec),
            &expected,
            &pose,
        );
    }
}

#[derive(Clone, Copy, Debug)]
enum Family {
    ConeSphere,
    ConeCone,
    ConeTorus,
    CylinderTorus,
    SphereTorus,
    TorusTorus,
}
fn family_operands(family: Family, model: &mut Model, pose: &Pose) -> (BodyId, BodyId) {
    match family {
        Family::ConeSphere => (
            topology::cone(model, r(2), r(0), r(2), pose).unwrap(),
            topology::sphere(model, q("3/4"), &shifted(pose, [0, 0, 1])).unwrap(),
        ),
        Family::ConeCone => (
            topology::cone(model, r(3), r(1), r(4), pose).unwrap(),
            topology::cone(model, r(1), r(3), r(4), pose).unwrap(),
        ),
        Family::ConeTorus => (
            topology::cone(model, r(1), r(5), r(4), &shifted(pose, [0, 0, -2])).unwrap(),
            topology::torus(model, r(3), r(1), pose).unwrap(),
        ),
        Family::CylinderTorus => (
            topology::cylinder(model, r(3), r(4), &shifted(pose, [0, 0, -2])).unwrap(),
            topology::torus(model, r(3), r(1), pose).unwrap(),
        ),
        Family::SphereTorus => (
            topology::sphere(model, r(3), pose).unwrap(),
            topology::torus(model, r(3), r(1), pose).unwrap(),
        ),
        Family::TorusTorus => (
            topology::torus(model, r(3), r(1), pose).unwrap(),
            topology::torus(model, r(3), r(1), &shifted(pose, [0, 0, 1])).unwrap(),
        ),
    }
}
fn family_volumes(family: Family) -> [Ball; 2] {
    let pi = Ball::pi(BITS);
    let t = &b(6) * &pi.square();
    match family {
        Family::ConeSphere => [&fraction(8, 3) * &pi, &fraction(9, 16) * &pi],
        Family::ConeCone => [&fraction(52, 3) * &pi, &fraction(52, 3) * &pi],
        Family::ConeTorus => [&fraction(124, 3) * &pi, t],
        Family::CylinderTorus | Family::SphereTorus => [&b(36) * &pi, t],
        Family::TorusTorus => [t.clone(), t],
    }
}
fn family_chi(family: Family, op: BoolOp) -> &'static [i64] {
    match (family, op) {
        (Family::ConeSphere, BoolOp::Subtract) => &[2, 2],
        (Family::ConeSphere, _) => &[2],
        (Family::ConeCone, BoolOp::Subtract) => &[0],
        (Family::ConeCone, _) => &[2],
        (Family::TorusTorus, _) | (_, BoolOp::Intersect) => &[0],
        _ => &[2],
    }
}
#[test]
fn coaxial_pair_families_conserve_volume_in_generated_poses() {
    for pose in poses().into_iter().take(2) {
        for family in [
            Family::ConeSphere,
            Family::ConeCone,
            Family::ConeTorus,
            Family::CylinderTorus,
            Family::SphereTorus,
            Family::TorusTorus,
        ] {
            let mut model = Model::new();
            let (a, b) = family_operands(family, &mut model, &pose);
            let properties = [BoolOp::Unite, BoolOp::Subtract, BoolOp::Intersect].map(|op| {
                let mut dec = Decider::default();
                let answer = boolean(&mut model, a, b, op, &mut dec)
                    .unwrap_or_else(|e| panic!("{family:?} {op:?}: {e:?}"));
                assert!(answer.marks.is_empty());
                measured(&model, answer.body, family_chi(family, op), &mut dec)
            });
            let [va, vb] = family_volumes(family);
            let sum = Enclosure {
                lo: &properties[0].volume.lo + &properties[2].volume.lo,
                hi: &properties[0].volume.hi + &properties[2].volume.hi,
            };
            agrees(&sum, &(&va + &vb));
            let difference = Enclosure {
                lo: &properties[1].volume.lo + &properties[2].volume.lo,
                hi: &properties[1].volume.hi + &properties[2].volume.hi,
            };
            agrees(&difference, &va);
            assert!(properties[2].volume.lo > r(0));
            assert!(properties[2].volume.hi < va.lo().min(vb.lo()));
        }
    }
}
