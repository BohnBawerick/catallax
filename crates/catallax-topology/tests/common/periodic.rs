//! Authored square stock with a through bore and one explicit periodic seam.
use catallax_numbers::{Ball, Rational};
use catallax_shapes::{Circle, Curve, Cylinder, Line, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{
    Sense,
    numerical::{
        NumericalDefinition, NumericalEdge, NumericalFace, NumericalTrim, NumericalVertex,
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
pub fn bore() -> NumericalDefinition {
    let positions = [
        [-2, -2, 0],
        [2, -2, 0],
        [2, 2, 0],
        [-2, 2, 0],
        [-2, -2, 1],
        [2, -2, 1],
        [2, 2, 1],
        [-2, 2, 1],
        [1, 0, 0],
        [1, 0, 1],
    ];
    let mut d = NumericalDefinition {
        name: "authored-periodic-bore".into(),
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
    let rings = [
        vec![0, 3, 2, 1],
        vec![4, 5, 6, 7],
        vec![0, 1, 5, 4],
        vec![1, 2, 6, 5],
        vec![2, 3, 7, 6],
        vec![3, 0, 4, 7],
    ];
    let origins = [
        [0, 0, 0],
        [0, 0, 1],
        [0, -2, 0],
        [2, 0, 0],
        [0, 2, 0],
        [-2, 0, 0],
    ];
    let normals = [
        [0, 0, 1],
        [0, 0, 1],
        [0, -1, 0],
        [1, 0, 0],
        [0, 1, 0],
        [-1, 0, 0],
    ];
    let references = [
        [1, 0, 0],
        [1, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [-1, 0, 0],
        [0, -1, 0],
    ];
    for (f, ring) in rings.into_iter().enumerate() {
        let uv = |i: usize| {
            let p = positions[i];
            match f {
                0 | 1 => [p[0].into(), p[1].into()],
                2 => [p[0].into(), p[2].into()],
                3 => [p[1].into(), p[2].into()],
                4 => [(-p[0]).into(), p[2].into()],
                _ => [(-p[1]).into(), p[2].into()],
            }
        };
        let mut trims = Vec::new();
        for i in 0..ring.len() {
            let ends = [ring[i], ring[(i + 1) % ring.len()]];
            let mut pair = ends;
            pair.sort();
            let edge = *indices.entry(pair).or_insert_with(|| {
                let edge = d.edges.len();
                let a = positions[ends[0]];
                let b = positions[ends[1]];
                d.edges.push(NumericalEdge {
                    key: format!("e{edge}"),
                    curve: Curve::Line(Line {
                        origin: point(a),
                        dir: vector(std::array::from_fn(|k| b[k] - a[k])),
                    }),
                    parameters: [q("0"), q("1")],
                    vertices: ends,
                });
                edge
            });
            trims.push(NumericalTrim {
                edge,
                sense: if d.edges[edge].vertices == ends {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                uv: ends.map(uv),
            });
        }
        d.faces.push(NumericalFace {
            key: format!("f{f}"),
            surface: Surface::Plane(Plane {
                origin: point(origins[f]),
                normal: vector(normals[f]),
                ref_dir: vector(references[f]),
            }),
            sense: if f == 0 { Sense::Reversed } else { Sense::Same },
            loops: vec![trims],
        });
    }
    let period = Ball::pi(128).mid() * q("2");
    let circles = [d.edges.len(), d.edges.len() + 1];
    for (height, vertex) in [(0, 8), (1, 9)] {
        d.edges.push(NumericalEdge {
            key: format!("rim-{height}"),
            curve: Curve::Circle(Circle {
                center: point([0, 0, height]),
                normal: vector([0, 0, 1]),
                radius_sq: 1.into(),
                ref_dir: vector([1, 0, 0]),
            }),
            vertices: [vertex, vertex],
            parameters: [q("0"), period.clone()],
        });
        d.faces[height as usize].loops.push(vec![NumericalTrim {
            edge: circles[height as usize],
            sense: if height == 0 {
                Sense::Same
            } else {
                Sense::Reversed
            },
            uv: [[q("1"), q("0")], [q("1"), q("0")]],
        }]);
    }
    let seam = d.edges.len();
    d.edges.push(NumericalEdge {
        key: "cylinder-u-seam".into(),
        curve: Curve::Line(Line {
            origin: point([1, 0, 0]),
            dir: vector([0, 0, 1]),
        }),
        vertices: [8, 9],
        parameters: [q("0"), q("1")],
    });
    d.faces.push(NumericalFace {
        key: "bore-wall".into(),
        surface: Surface::Cylinder(Cylinder {
            origin: point([0, 0, 0]),
            axis: vector([0, 0, 1]),
            radius_sq: 1.into(),
            ref_dir: vector([1, 0, 0]),
        }),
        sense: Sense::Reversed,
        loops: vec![vec![
            NumericalTrim {
                edge: circles[0],
                sense: Sense::Reversed,
                uv: [[period.clone(), q("0")], [q("0"), q("0")]],
            },
            NumericalTrim {
                edge: seam,
                sense: Sense::Same,
                uv: [[q("0"), q("0")], [q("0"), q("1")]],
            },
            NumericalTrim {
                edge: circles[1],
                sense: Sense::Same,
                uv: [[q("0"), q("1")], [period.clone(), q("1")]],
            },
            NumericalTrim {
                edge: seam,
                sense: Sense::Reversed,
                uv: [[period.clone(), q("1")], [period, q("0")]],
            },
        ]],
    });
    d
}
