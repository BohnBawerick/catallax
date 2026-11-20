//! An authored quarter cylinder, with retained circles and shared 3D vertices.
use catallax_numbers::{Ball, Rational};
use catallax_shapes::{Circle, Curve, Cylinder, Line, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{
    Sense,
    numerical::{
        NumericalDefinition, NumericalEdge, NumericalFace, NumericalNativeTrim, NumericalTrim,
        NumericalVertex,
    },
};
use std::collections::BTreeMap;

pub fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn point(p: [i64; 3]) -> catallax_shapes::Point3 {
    RPoint3::new(p[0].into(), p[1].into(), p[2].into()).to_real()
}
fn vector(p: [i64; 3]) -> catallax_shapes::Vector3 {
    RVector3::new(p[0].into(), p[1].into(), p[2].into()).to_real()
}
pub fn quarter() -> NumericalDefinition {
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
                tag: catallax_topology::numerical::CopyTag::Independent,
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

/// A native u/v copy of one circle use on the wall (face 4): the right u at the
/// right parameter, but v bowed by three quarters of `bow` in the middle, so
/// its ends still close the loop. Returns the binding and the edge's key.
#[allow(dead_code)]
pub fn bowed_copy(source: &NumericalDefinition, bow: &str) -> (NumericalNativeTrim, String) {
    use catallax_shapes::{BSplineCurve, Point3};
    let (use_index, trim) = source.faces[4].loops[0]
        .iter()
        .enumerate()
        .find(|(_, trim)| matches!(source.edges[trim.edge].curve, Curve::Circle(_)))
        .unwrap();
    let edge = &source.edges[trim.edge];
    let (lo, hi) = (
        edge.parameters[0].clone().min(edge.parameters[1].clone()),
        edge.parameters[0].clone().max(edge.parameters[1].clone()),
    );
    // The wall's u is the circle's angle; its v is the circle's height.
    let height = trim.uv[0][1].clone();
    let third = (&hi - &lo) / q("3");
    let us = [lo.clone(), &lo + &third, &hi - &third, hi.clone()];
    let vs = [
        height.clone(),
        &height + &q(bow),
        &height + &q(bow),
        height.clone(),
    ];
    let controls = us
        .iter()
        .zip(&vs)
        .map(|(u, v)| Point3::new(u.clone().into(), v.clone().into(), 0.into()))
        .collect();
    let knots = [
        lo.clone(),
        lo.clone(),
        lo.clone(),
        lo,
        hi.clone(),
        hi.clone(),
        hi.clone(),
        hi,
    ]
    .to_vec();
    let curve = BSplineCurve::new(
        3,
        knots,
        controls,
        vec![catallax_numbers::Real::one(); 4],
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    (
        NumericalNativeTrim {
            face: source.faces[4].key.clone(),
            boundary: 0,
            use_index,
            edge: edge.key.clone(),
            curve,
            tag: catallax_topology::numerical::CopyTag::Independent,
        },
        edge.key.clone(),
    )
}

/// The quarter cylinder scaled by `k` about the origin.
#[allow(dead_code)]
pub fn scaled(mut d: NumericalDefinition, k: &Rational) -> NumericalDefinition {
    use catallax_shapes::Point3;
    let real = |r: &Rational| -> catallax_numbers::Real { r.clone().into() };
    let point = |p: &Point3| {
        Point3::new(
            p.x.clone() * real(k),
            p.y.clone() * real(k),
            p.z.clone() * real(k),
        )
    };
    for vertex in &mut d.vertices {
        for c in &mut vertex.position {
            *c = &*c * k;
        }
    }
    for edge in &mut d.edges {
        match &mut edge.curve {
            Curve::Line(line) => {
                line.origin = point(&line.origin);
                line.dir = line.dir.scaled(&real(k));
            }
            Curve::Circle(circle) => {
                circle.center = point(&circle.center);
                circle.radius_sq = circle.radius_sq.clone() * real(&(k * k));
            }
            _ => unreachable!(),
        }
    }
    for face in &mut d.faces {
        let planar = matches!(face.surface, Surface::Plane(_));
        match &mut face.surface {
            Surface::Plane(plane) => plane.origin = point(&plane.origin),
            Surface::Cylinder(cylinder) => {
                cylinder.origin = point(&cylinder.origin);
                cylinder.radius_sq = cylinder.radius_sq.clone() * real(&(k * k));
            }
            _ => unreachable!(),
        }
        for trim in face.loops.iter_mut().flatten() {
            for uv in &mut trim.uv {
                if planar {
                    uv[0] = &uv[0] * k;
                }
                uv[1] = &uv[1] * k;
            }
        }
    }
    d
}
