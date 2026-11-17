//! Circular iso-curves on a ring torus, including a wrapped-contact counterexample.
use catallax_numbers::{Ball, Rational};
use catallax_shapes::{Circle, Curve, RPoint3, RVector3, Surface, Torus};
use catallax_topology::{Sense, numerical as n};
use std::collections::BTreeMap;

fn q(value: f64) -> Rational {
    Rational::from_f64_exact(value).unwrap()
}
fn point(value: [f64; 3]) -> catallax_shapes::Point3 {
    RPoint3::new(q(value[0]), q(value[1]), q(value[2])).to_real()
}
fn vector(value: [f64; 3]) -> catallax_shapes::Vector3 {
    RVector3::new(q(value[0]), q(value[1]), q(value[2])).to_real()
}
pub fn period() -> Rational {
    Ball::pi(128).mid() * Rational::from_i64(2)
}

fn sin_cos(angle: &Rational) -> (f64, f64) {
    // Cartesian quadrant points are exact; their stored angles remain numerical.
    let quarter = period() / Rational::from_i64(4);
    for (index, pair) in [
        (0, (0.0, 1.0)),
        (1, (1.0, 0.0)),
        (2, (0.0, -1.0)),
        (3, (-1.0, 0.0)),
        (4, (0.0, 1.0)),
    ] {
        if angle == &(&quarter * &Rational::from_i64(index)) {
            return pair;
        }
    }
    angle.to_f64_nearest().sin_cos()
}

fn patch(corners: Vec<[Rational; 2]>, seam: Option<[usize; 2]>) -> n::NumericalDefinition {
    let mut result = n::NumericalDefinition {
        name: "torus-sheet".into(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let count = corners.len();
    let mut canonical: Vec<_> = (0..count).collect();
    if let Some([a, b]) = seam {
        canonical[b] = (a + 1) % count;
        canonical[(b + 1) % count] = a;
    }
    let mut indices = BTreeMap::new();
    let vertices: Vec<_> = canonical
        .iter()
        .map(|&corner| {
            *indices.entry(corner).or_insert_with(|| {
                let (s, c) = sin_cos(&corners[corner][0]);
                let (sv, cv) = sin_cos(&corners[corner][1]);
                let index = result.vertices.len();
                result.vertices.push(n::NumericalVertex {
                    key: format!("v{index}"),
                    position: [(3.0 + cv) * c, (3.0 + cv) * s, sv].map(q),
                });
                index
            })
        })
        .collect();
    let mut trims: Vec<n::NumericalTrim> = Vec::new();
    for i in 0..count {
        let uv = [corners[i].clone(), corners[(i + 1) % count].clone()];
        if let Some([a, b]) = seam
            && i == b
        {
            trims.push(n::NumericalTrim {
                edge: trims[a].edge,
                sense: if trims[a].sense == Sense::Same {
                    Sense::Reversed
                } else {
                    Sense::Same
                },
                uv,
            });
            continue;
        }
        let axis = usize::from(uv[0][0] == uv[1][0]);
        let fixed = &uv[0][1 - axis];
        let (s, c) = sin_cos(fixed);
        let (center, normal, reference, radius) = if axis == 0 {
            ([0.0, 0.0, s], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 3.0 + c)
        } else {
            ([3.0 * c, 3.0 * s, 0.0], [s, -c, 0.0], [c, s, 0.0], 1.0)
        };
        let forward = uv[0][axis] < uv[1][axis];
        let ends = if forward {
            [i, (i + 1) % count]
        } else {
            [(i + 1) % count, i]
        };
        let edge = result.edges.len();
        let radius = q(radius);
        result.edges.push(n::NumericalEdge {
            key: format!("e{edge}"),
            curve: Curve::Circle(Circle {
                center: point(center),
                normal: vector(normal),
                ref_dir: vector(reference),
                radius_sq: (&radius * &radius).into(),
            }),
            parameters: ends.map(|end| corners[end][axis].clone()),
            vertices: ends.map(|end| vertices[end]),
        });
        trims.push(n::NumericalTrim {
            edge,
            sense: if forward {
                Sense::Same
            } else {
                Sense::Reversed
            },
            uv,
        });
    }
    result.faces.push(n::NumericalFace {
        key: "torus".into(),
        surface: Surface::Torus(Torus {
            center: point([0.0; 3]),
            axis: vector([0.0, 0.0, 1.0]),
            major: 3.into(),
            minor: 1.into(),
            ref_dir: vector([1.0, 0.0, 0.0]),
        }),
        sense: Sense::Same,
        loops: vec![trims],
    });
    result
}

pub fn quarter() -> n::NumericalDefinition {
    let h = period() / Rational::from_i64(4);
    patch(
        vec![
            [q(0.0), q(0.0)],
            [h.clone(), q(0.0)],
            [h.clone(), h.clone()],
            [q(0.0), h],
        ],
        None,
    )
}
pub fn band() -> n::NumericalDefinition {
    let t = period();
    let h = &t / &Rational::from_i64(4);
    patch(
        vec![
            [q(0.0), q(0.0)],
            [t.clone(), q(0.0)],
            [t, h.clone()],
            [q(0.0), h],
        ],
        Some([1, 3]),
    )
}
pub fn unpaired_wrapped_contact() -> n::NumericalDefinition {
    let t = period();
    let inset = &t - &q(1.0);
    let half = &t / &Rational::from_i64(2);
    patch(
        vec![
            [q(0.0), q(0.0)],
            [t.clone(), q(0.0)],
            [t.clone(), q(1.0)],
            [inset.clone(), q(1.0)],
            [inset, q(2.0)],
            [t.clone(), q(2.0)],
            [t, q(3.0)],
            [half, q(3.0)],
            [q(0.0), q(3.0)],
            [q(0.0), q(2.0)],
            [q(1.0), q(2.0)],
            [q(1.0), q(1.0)],
            [q(0.0), q(1.0)],
        ],
        Some([1, 12]),
    )
}
