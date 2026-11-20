//! N02 square tube construction copied without geometry changes from
//! catallax-topology/tests/numerical_trim.rs at N02 merge 3e2a929e.
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Curve, Line, Plane, Point3, Surface, Vector3};
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
fn point(p: [i64; 3]) -> Point3 {
    Point3::new(p[0].into(), p[1].into(), p[2].into())
}
fn vector(p: [i64; 3]) -> Vector3 {
    Vector3::new(p[0].into(), p[1].into(), p[2].into())
}
fn sub(a: [i64; 3], b: [i64; 3]) -> [i64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [i64; 3], b: [i64; 3]) -> [i64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: [i64; 3]) -> [i64; 3] {
    a.map(|x| x.signum())
}
fn dot(a: [i64; 3], b: [i64; 3]) -> i64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

pub fn tube() -> NumericalDefinition {
    let outer = [[0, 0], [4, 0], [4, 4], [0, 4]];
    let hole = [[1, 1], [1, 3], [3, 3], [3, 1]];
    let at = |p: [i64; 2], z| [p[0], p[1], z];
    let mut polygons = vec![
        vec![
            outer.map(|p| at(p, 1)).to_vec(),
            hole.map(|p| at(p, 1)).to_vec(),
        ],
        vec![
            outer.into_iter().rev().map(|p| at(p, 0)).collect(),
            hole.into_iter().rev().map(|p| at(p, 0)).collect(),
        ],
    ];
    for ring in [outer, hole] {
        for i in 0..4 {
            let (a, b) = (ring[i], ring[(i + 1) % 4]);
            polygons.push(vec![vec![at(a, 0), at(b, 0), at(b, 1), at(a, 1)]]);
        }
    }
    let mut result = NumericalDefinition {
        name: "square-through-hole".into(),
        vertices: vec![],
        edges: vec![],
        faces: vec![],
    };
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    let mut dec = Decider::default();
    for (f, loops) in polygons.into_iter().enumerate() {
        let p = &loops[0];
        let origin = p[0];
        let x = unit(sub(p[1], origin));
        let normal = unit(cross(x, sub(p[2], p[1])));
        let y = cross(normal, x);
        let surface =
            Surface::Plane(Plane::new(point(origin), vector(normal), vector(x), &mut dec).unwrap());
        let mut trims = vec![];
        for polygon in loops {
            let mut boundary = vec![];
            for i in 0..polygon.len() {
                let ends = [polygon[i], polygon[(i + 1) % polygon.len()]];
                let ids = ends.map(|p| {
                    *vertices.entry(p).or_insert_with(|| {
                        let index = result.vertices.len();
                        result.vertices.push(NumericalVertex {
                            key: format!("v{index}"),
                            position: p.map(Rational::from),
                        });
                        index
                    })
                });
                let mut key = ids;
                key.sort();
                let edge = *edges.entry(key).or_insert_with(|| {
                    let index = result.edges.len();
                    result.edges.push(NumericalEdge {
                        key: format!("e{index}"),
                        curve: Curve::Line(
                            Line::new(point(ends[0]), vector(sub(ends[1], ends[0])), &mut dec)
                                .unwrap(),
                        ),
                        vertices: ids,
                        parameters: [q("0"), q("1")],
                    });
                    index
                });
                boundary.push(NumericalTrim {
                    tag: catallax_topology::numerical::CopyTag::Independent,
                    edge,
                    sense: if result.edges[edge].vertices == ids {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                    uv: ends.map(|p| {
                        [
                            Rational::from(dot(sub(p, origin), x)),
                            Rational::from(dot(sub(p, origin), y)),
                        ]
                    }),
                });
            }
            trims.push(boundary);
        }
        result.faces.push(NumericalFace {
            key: format!("f{f}"),
            surface,
            sense: Sense::Same,
            loops: trims,
        });
    }
    result
}
