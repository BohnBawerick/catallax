//! Shared real B-rep fixtures for numerical admission and delivery tests.
#![allow(dead_code)]

use catallax_app::{
    NumericalDefinition, NumericalEdge, NumericalFace, NumericalTrim, NumericalVertex,
};
use catallax_numbers::{Decider, Rational, Real};
use catallax_shapes::{BSplineSurface, Curve, Point3, Pose, RationalPatch, Surface, Vector3};
use catallax_topology::{Loop, Model, Sense};
use std::collections::BTreeMap;

pub fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn axis_length(v: &Vector3) -> Rational {
    let coordinates = [&v.x, &v.y, &v.z].map(|x| x.as_rational().unwrap());
    assert_eq!(coordinates.iter().filter(|x| !x.is_zero()).count(), 1);
    coordinates.into_iter().fold(q("0"), |sum, x| sum + x.abs())
}

pub fn roof(delta: &str) -> Surface {
    let mut controls = Vec::new();
    for (i, x) in ["-1", "0.5", "2"].iter().enumerate() {
        for (j, y) in ["-1", "0.5", "2"].iter().enumerate() {
            let z = q("1") + if i == 1 && j == 1 { q(delta) } else { q("0") };
            controls.push(Point3::new(q(x).into(), q(y).into(), z.into()));
        }
    }
    let knots = ["0", "0", "0", "1", "1", "1"].map(q).to_vec();
    let mut dec = Decider::default();
    Surface::BSpline(
        RationalPatch::new(
            BSplineSurface::new(
                [2, 2],
                [knots.clone(), knots],
                [3, 3],
                controls,
                vec![Real::one(); 9],
                &mut dec,
            )
            .unwrap(),
            &mut dec,
        )
        .unwrap()
        .into(),
    )
}

pub fn fixture(delta: &str) -> NumericalDefinition {
    let mut model = Model::new();
    let body =
        catallax_topology::block(&mut model, [q("1"), q("1"), q("1")], &Pose::identity()).unwrap();
    from_model(&model, body, delta)
}

pub fn from_model(
    model: &Model,
    body: catallax_topology::BodyId,
    delta: &str,
) -> NumericalDefinition {
    let vertex_ids = model.vertices_of_body(body);
    let vertex_map: BTreeMap<_, _> = vertex_ids
        .iter()
        .enumerate()
        .map(|(i, id)| (*id, i))
        .collect();
    let vertices: Vec<_> = vertex_ids
        .iter()
        .enumerate()
        .map(|(i, id)| {
            let point = model.geom().position(model.vertex(*id).point);
            NumericalVertex {
                key: format!("v{i}"),
                position: [&point.x, &point.y, &point.z].map(|v| v.as_rational().unwrap()),
            }
        })
        .collect();
    let edge_ids = model.edges_of_body(body);
    let edge_map: BTreeMap<_, _> = edge_ids
        .iter()
        .enumerate()
        .map(|(i, id)| (*id, i))
        .collect();
    let edges = edge_ids
        .iter()
        .enumerate()
        .map(|(i, id)| {
            let edge = model.edge(*id);
            let curve = model.geom().curve(edge.curve).clone();
            let Curve::Line(line) = &curve else {
                panic!("line")
            };
            let ends = [edge.start.unwrap(), edge.end.unwrap()];
            let parameters = ends.map(|v| {
                model
                    .geom()
                    .position(model.vertex(v).point)
                    .minus(&line.origin)
                    .dot(&line.dir)
                    .as_rational()
                    .unwrap()
                    / line.dir.norm_squared().as_rational().unwrap()
            });
            NumericalEdge {
                key: format!("e{i}"),
                curve,
                vertices: ends.map(|v| vertex_map[&v]),
                parameters,
            }
        })
        .collect();
    let faces = model
        .faces_of_body(body)
        .into_iter()
        .enumerate()
        .map(|(i, id)| {
            let face = model.face(id);
            let Surface::Plane(plane) = model.geom().surface(face.surface) else {
                panic!("plane")
            };
            let top = plane.normal.x.as_rational() == Some(q("0"))
                && plane.normal.y.as_rational() == Some(q("0"))
                && plane.normal.z.as_rational().unwrap() > q("0");
            let Loop::Coedges(uses) = model.loop_(face.loops[0]) else {
                panic!("coedges")
            };
            let boundary = uses
                .iter()
                .map(|id| {
                    let c = model.coedge(*id);
                    let edge = model.edge(c.edge);
                    let mut ends = [edge.start.unwrap(), edge.end.unwrap()];
                    if c.sense == Sense::Reversed {
                        ends.reverse();
                    }
                    let uv = ends.map(|v| {
                        let p = model.geom().position(model.vertex(v).point);
                        if top {
                            [
                                (p.x.as_rational().unwrap() + q("1")) / q("3"),
                                (p.y.as_rational().unwrap() + q("1")) / q("3"),
                            ]
                        } else {
                            let d = p.minus(&plane.origin);
                            [
                                d.dot(&plane.ref_dir).as_rational().unwrap()
                                    / axis_length(&plane.ref_dir),
                                d.dot(&plane.normal.cross(&plane.ref_dir))
                                    .as_rational()
                                    .unwrap()
                                    / axis_length(&plane.normal.cross(&plane.ref_dir)),
                            ]
                        }
                    });
                    NumericalTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        edge: edge_map[&c.edge],
                        sense: c.sense,
                        uv,
                    }
                })
                .collect();
            NumericalFace {
                key: if top {
                    "roof".into()
                } else {
                    format!("side{i}")
                },
                surface: if top {
                    roof(delta)
                } else {
                    Surface::Plane(plane.clone())
                },
                sense: face.sense,
                loops: vec![boundary],
            }
        })
        .collect();
    NumericalDefinition {
        name: "general-trim-roof".into(),
        vertices,
        edges,
        faces,
    }
}
