//! An authored quarter cylinder, with retained circles and shared 3D vertices.
use catallax_numbers::{Ball, Decider, NumericalCheck, Rational};
use catallax_shapes::{Circle, Curve, Cylinder, Line, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{
    Sense,
    numerical::{
        self as n, NumericalDefinition, NumericalEdge, NumericalFace, NumericalTrim,
        NumericalVertex,
    },
};
use std::collections::BTreeMap;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn point(p: [i64; 3]) -> catallax_shapes::Point3 {
    RPoint3::new(p[0].into(), p[1].into(), p[2].into()).to_real()
}
fn vector(p: [i64; 3]) -> catallax_shapes::Vector3 {
    RVector3::new(p[0].into(), p[1].into(), p[2].into()).to_real()
}
fn quarter() -> NumericalDefinition {
    let positions = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
    ];
    let mut result = NumericalDefinition {
        name: "quarter-cylinder".into(),
        vertices: positions
            .iter()
            .enumerate()
            .map(|(i, p)| NumericalVertex {
                key: format!("v{i}"),
                position: p.map(Rational::from_i64),
            })
            .collect(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let mut indices = BTreeMap::new();
    let half_pi = Ball::pi(128).mid() / q("2");
    let polygons = [
        vec![0, 2, 1],
        vec![3, 4, 5],
        vec![0, 1, 4, 3],
        vec![0, 3, 5, 2],
        vec![1, 2, 5, 4],
    ];
    for (f, polygon) in polygons.into_iter().enumerate() {
        let uv = |i: usize| {
            let p = positions[i];
            match f {
                0 | 1 => [p[0].into(), p[1].into()],
                2 => [p[0].into(), p[2].into()],
                3 => [p[2].into(), p[1].into()],
                _ => [
                    if p[0] == 1 { q("0") } else { half_pi.clone() },
                    p[2].into(),
                ],
            }
        };
        let surface = match f {
            0 | 1 => Surface::Plane(Plane {
                origin: point([0, 0, f as i64]),
                normal: vector([0, 0, 1]),
                ref_dir: vector([1, 0, 0]),
            }),
            2 => Surface::Plane(Plane {
                origin: point([0, 0, 0]),
                normal: vector([0, -1, 0]),
                ref_dir: vector([1, 0, 0]),
            }),
            3 => Surface::Plane(Plane {
                origin: point([0, 0, 0]),
                normal: vector([-1, 0, 0]),
                ref_dir: vector([0, 0, 1]),
            }),
            _ => Surface::Cylinder(Cylinder {
                origin: point([0, 0, 0]),
                axis: vector([0, 0, 1]),
                radius_sq: 1.into(),
                ref_dir: vector([1, 0, 0]),
            }),
        };
        let mut trims = Vec::new();
        for i in 0..polygon.len() {
            let ends = [polygon[i], polygon[(i + 1) % polygon.len()]];
            let mut pair = ends;
            pair.sort();
            let index = *indices.entry(pair).or_insert_with(|| {
                let index = result.edges.len();
                let a = positions[ends[0]];
                let b = positions[ends[1]];
                let circular = pair == [1, 2] || pair == [4, 5];
                let (curve, parameters) = if circular {
                    (
                        Curve::Circle(Circle {
                            center: point([0, 0, a[2]]),
                            normal: vector([0, 0, 1]),
                            radius_sq: 1.into(),
                            ref_dir: vector([1, 0, 0]),
                        }),
                        ends.map(|i| {
                            if positions[i][0] == 1 {
                                q("0")
                            } else {
                                half_pi.clone()
                            }
                        }),
                    )
                } else {
                    (
                        Curve::Line(Line {
                            origin: point(a),
                            dir: vector(std::array::from_fn(|k| b[k] - a[k])),
                        }),
                        [q("0"), q("1")],
                    )
                };
                result.edges.push(NumericalEdge {
                    key: format!("e{index}"),
                    curve,
                    parameters,
                    vertices: ends,
                });
                index
            });
            trims.push(NumericalTrim {
                edge: index,
                sense: if result.edges[index].vertices == ends {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                uv: ends.map(uv),
            });
        }
        result.faces.push(NumericalFace {
            key: format!("f{f}"),
            surface,
            sense: if f == 0 { Sense::Reversed } else { Sense::Same },
            loops: vec![trims],
        });
    }
    result
}

#[test]
fn finite_circular_caps_and_cylinder_keep_numerical_checks() {
    let body = n::build(
        quarter(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let count = body.counts();
    assert_eq!((count.vertices, count.edges, count.faces), (6, 9, 5));
    assert!(body.evidence().checks().iter().any(|check| matches!(check,NumericalCheck::Passed { name,.. } if name=="trim.sampled-simple-oriented-parameter-loop")));
    assert!(
        body.evidence()
            .checks()
            .iter()
            .filter_map(|check| match check {
                NumericalCheck::Samples { samples, .. } => Some(samples),
                _ => None,
            })
            .flatten()
            .all(|sample| sample.distance_upper_mm <= q("0.01"))
    );
    assert!(!body.exact_topology_faults().is_empty());
}

#[test]
fn a_cylinder_residual_over_the_source_tolerance_refuses() {
    for (radius, pass) in [("1.0099", true), ("1.0101", false)] {
        let mut source = quarter();
        let Surface::Cylinder(cylinder) = &mut source.faces[4].surface else {
            unreachable!()
        };
        cylinder.radius_sq = (&q(radius) * &q(radius)).into();
        let result = n::build(
            source,
            n::NumericalOptions::default(),
            &mut Decider::default(),
        );
        assert_eq!(result.is_ok(), pass);
        if !pass {
            assert_eq!(
                result.unwrap_err().code,
                "numbers.numerical-tolerance-exceeded"
            );
        }
    }
}

#[test]
fn reversed_cylinder_winding_and_unresolved_periodic_cover_refuse() {
    let mut source = quarter();
    source.faces[4].sense = Sense::Reversed;
    assert!(
        n::build(
            source,
            n::NumericalOptions::default(),
            &mut Decider::default()
        )
        .is_err()
    );
    let mut source = quarter();
    source.faces[4].loops[0][0].uv[1][0] = Ball::pi(128).hi() * q("2");
    assert!(
        n::build(
            source,
            n::NumericalOptions::default(),
            &mut Decider::default()
        )
        .is_err()
    );
}
