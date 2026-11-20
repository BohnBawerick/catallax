//! A tent: a 100 mm square base at z = 0 and an apex above its centre, with
//! the apex vertex lifted off the ends of its four edges. The apex is near
//! the base's interior while every edge and vertex of the base is tens of
//! millimetres away.
use catallax_numbers::Rational;
use catallax_shapes::{Curve, Line, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{
    Sense,
    numerical::{
        CopyTag, NumericalDefinition, NumericalEdge, NumericalFace, NumericalTrim, NumericalVertex,
    },
};
use std::collections::BTreeMap;

/// The base face's key.
#[allow(dead_code)]
pub const BASE: &str = "base";
/// The apex vertex's key.
#[allow(dead_code)]
pub const APEX: &str = "apex";

fn rational(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn sub(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|k| &a[k] - &b[k])
}

fn cross(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

fn dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, k| sum + &a[k] * &b[k])
}

fn point(p: &[Rational; 3]) -> catallax_shapes::Point3 {
    RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()).to_real()
}

fn vector(p: &[Rational; 3]) -> catallax_shapes::Vector3 {
    RVector3::new(p[0].clone(), p[1].clone(), p[2].clone()).to_real()
}

/// The tent whose apex lies at `height` above the base's centre, its vertex
/// a further `lift` above that, so the vertex is `height + lift` above the
/// base and `lift` off its edges' ends.
#[allow(dead_code)]
pub fn tent(height: &str, lift: &str) -> NumericalDefinition {
    let h = rational(height);
    let corner = |x: i64, y: i64| {
        [
            Rational::from_i64(x),
            Rational::from_i64(y),
            Rational::zero(),
        ]
    };
    let corners = [
        corner(0, 0),
        corner(100, 0),
        corner(100, 100),
        corner(0, 100),
    ];
    let apex = [Rational::from_i64(50), Rational::from_i64(50), h.clone()];
    let mut lifted = apex.clone();
    lifted[2] = &h + &rational(lift);
    let mut points: Vec<[Rational; 3]> = corners.to_vec();
    points.push(apex.clone());
    let mut d = NumericalDefinition {
        name: "tent".into(),
        vertices: corners
            .iter()
            .enumerate()
            .map(|(i, p)| NumericalVertex {
                key: format!("b{i}"),
                position: p.clone(),
            })
            .chain(std::iter::once(NumericalVertex {
                key: APEX.into(),
                position: lifted,
            }))
            .collect(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    // Each face, counter-clockwise around its outward normal.
    let faces: Vec<(&str, Vec<usize>)> = vec![
        (BASE, vec![0, 3, 2, 1]),
        ("front", vec![0, 1, 4]),
        ("right", vec![1, 2, 4]),
        ("back", vec![2, 3, 4]),
        ("left", vec![3, 0, 4]),
    ];
    let mut edges = BTreeMap::new();
    for (key, polygon) in faces {
        let (p0, p1, p2) = (
            &points[polygon[0]],
            &points[polygon[1]],
            &points[polygon[2]],
        );
        let reference = sub(p1, p0);
        let normal = cross(&reference, &sub(p2, p0));
        let in_plane = cross(&normal, &reference);
        let scale = |v: &[Rational; 3]| dot(v, v).to_f64_nearest().sqrt();
        let (r_len, n_len) = (scale(&reference), scale(&in_plane));
        let uv = |p: &[Rational; 3]| -> [Rational; 2] {
            let offset = sub(p, p0);
            [
                Rational::from_f64_exact(dot(&offset, &reference).to_f64_nearest() / r_len)
                    .unwrap(),
                Rational::from_f64_exact(dot(&offset, &in_plane).to_f64_nearest() / n_len).unwrap(),
            ]
        };
        let mut trims = Vec::new();
        for k in 0..polygon.len() {
            let ends = [polygon[k], polygon[(k + 1) % polygon.len()]];
            let mut pair = ends;
            pair.sort();
            let index = *edges.entry(pair).or_insert_with(|| {
                let (a, b) = (&points[pair[0]], &points[pair[1]]);
                d.edges.push(NumericalEdge {
                    key: format!("e{}{}", pair[0], pair[1]),
                    curve: Curve::Line(Line {
                        origin: point(a),
                        dir: vector(&sub(b, a)),
                    }),
                    vertices: pair,
                    parameters: [Rational::zero(), Rational::one()],
                });
                d.edges.len() - 1
            });
            trims.push(NumericalTrim {
                edge: index,
                sense: if d.edges[index].vertices == ends {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                uv: ends.map(|i| uv(&points[i])),
                tag: CopyTag::Independent,
            });
        }
        d.faces.push(NumericalFace {
            key: key.into(),
            surface: Surface::Plane(Plane {
                origin: point(p0),
                normal: vector(&normal),
                ref_dir: vector(&reference),
            }),
            sense: Sense::Same,
            loops: vec![trims],
        });
    }
    d
}
