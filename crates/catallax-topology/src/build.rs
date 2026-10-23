//! Primitive solid builders.

use std::collections::BTreeMap;

use catallax_numbers::{Fact, FactKind, FactOrigin, MarkLog, Rational, Refusal, Sign};
use catallax_shapes::{
    Circle, Cone, Curve, Cylinder, GeomKey, Line, Plane, Point3, PointDef, Pose, RPoint3, RVector3,
    Sphere, Surface, SurfaceId, Torus, Vector3,
};

use crate::{
    Body, BodyId, BodyKind, Coedge, Edge, EdgeId, Face, FaceId, Loop, LoopId, Lump, Model, Sense,
    Shell, Vertex,
};

/// The box `[0, sx] x [0, sy] x [0, sz]`, moved by `pose`.
pub fn block(model: &mut Model, size: [Rational; 3], pose: &Pose) -> Result<BodyId, Refusal> {
    if size.iter().any(|s| s.sign() != Sign::Positive) {
        return Err(non_positive());
    }
    let origin = "topology.block";
    let corners: [Point3; 8] = std::array::from_fn(|i| {
        let c: [Rational; 3] = std::array::from_fn(|axis| {
            if i & (4 >> axis) == 0 {
                Rational::zero()
            } else {
                size[axis].clone()
            }
        });
        let [x, y, z] = c;
        RPoint3::new(x, y, z).to_real()
    });
    let vertices: [_; 8] = std::array::from_fn(|i| {
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(pose.apply_point(&corners[i])));
        model.add_vertex(Vertex { point })
    });
    let mut edges = BTreeMap::new();
    for (axis, mask) in [4, 2, 1].into_iter().enumerate() {
        for i in 0..8 {
            if i & mask == 0 {
                let curve = Curve::Line(Line {
                    origin: corners[i].clone(),
                    dir: unit(axis, 1),
                });
                let curve = model
                    .geom_mut()
                    .add_curve(curve.transformed(pose).expect("analytic curve"));
                let edge = model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: Some(vertices[i]),
                    end: Some(vertices[i | mask]),
                });
                edges.insert((i, i | mask), edge);
            }
        }
    }
    let cycles = [
        [0, 1, 3, 2],
        [4, 6, 7, 5],
        [0, 4, 5, 1],
        [2, 3, 7, 6],
        [0, 2, 6, 4],
        [1, 5, 7, 3],
    ];
    let mut faces = Vec::with_capacity(6);
    let mut surfaces = Vec::with_capacity(6);
    for (i, cycle) in cycles.iter().enumerate() {
        let axis = i / 2;
        let plane = Surface::Plane(Plane {
            origin: corners[if i % 2 == 0 { 0 } else { 4 >> axis }].clone(),
            normal: unit(axis, if i % 2 == 0 { -1 } else { 1 }),
            ref_dir: unit((axis + 1) % 3, 1),
        });
        let surface = model.geom_mut().add_surface(plane.transformed(pose));
        surfaces.push(surface);
        let mut coedges = Vec::with_capacity(4);
        for j in 0..4 {
            let (a, b) = (cycle[j], cycle[(j + 1) % 4]);
            let edge = edges[&(a.min(b), a.max(b))];
            coedges.push(model.add_coedge(Coedge {
                edge,
                sense: if a < b { Sense::Same } else { Sense::Reversed },
            }));
        }
        let loop_ = model.add_loop(Loop::Coedges(coedges));
        faces.push(face(model, surface, vec![loop_], origin));
    }
    for pair in surfaces.as_chunks::<2>().0 {
        record(
            model,
            FactKind::Parallel,
            GeomKey::Surface(pair[0]),
            GeomKey::Surface(pair[1]),
            origin,
        );
    }
    Ok(solid(model, faces))
}

/// A cylinder on the local `z` axis from 0 to `height`, moved by `pose`.
pub fn cylinder(
    model: &mut Model,
    radius: Rational,
    height: Rational,
    pose: &Pose,
) -> Result<BodyId, Refusal> {
    if radius.sign() != Sign::Positive || height.sign() != Sign::Positive {
        return Err(non_positive());
    }
    let side = Surface::Cylinder(Cylinder {
        origin: Point3::origin(),
        axis: unit(2, 1),
        radius_sq: catallax_numbers::Real::from(radius.clone()).square(),
        ref_dir: unit(0, 1),
    });
    Ok(round_solid(
        model,
        [radius.clone(), radius],
        height,
        side,
        pose,
        "topology.cylinder",
    ))
}

/// A sphere centred at the local origin, moved by `pose`.
pub fn sphere(model: &mut Model, radius: Rational, pose: &Pose) -> Result<BodyId, Refusal> {
    if radius.sign() != Sign::Positive {
        return Err(non_positive());
    }
    let surface = Surface::Sphere(Sphere {
        center: Point3::origin(),
        radius: radius.into(),
        axis: unit(2, 1),
        ref_dir: unit(0, 1),
    });
    Ok(closed_solid(model, surface, pose))
}

/// A cone frustum on the local `z` axis, `base_radius` at 0 and `top_radius`
/// at `height` (0 for an apex), moved by `pose`.
pub fn cone(
    model: &mut Model,
    base_radius: Rational,
    top_radius: Rational,
    height: Rational,
    pose: &Pose,
) -> Result<BodyId, Refusal> {
    if base_radius.sign() != Sign::Positive
        || height.sign() != Sign::Positive
        || top_radius.sign() == Sign::Negative
    {
        return Err(non_positive());
    }
    if top_radius == base_radius {
        return Err(Refusal::new(
            "topology.cone-is-cylinder",
            "equal cone radii define a cylinder",
        ));
    }
    let difference = &base_radius - &top_radius;
    let apex_z = (&height * &base_radius)
        .checked_div(&difference)
        .expect("unequal radii");
    let slope = difference
        .abs()
        .checked_div(&height)
        .expect("positive height");
    let side = Surface::Cone(Cone {
        apex: at_z(apex_z),
        axis: unit(
            2,
            if difference.sign() == Sign::Positive {
                -1
            } else {
                1
            },
        ),
        slope: slope.into(),
        ref_dir: unit(0, 1),
    });
    Ok(round_solid(
        model,
        [base_radius, top_radius],
        height,
        side,
        pose,
        "topology.cone",
    ))
}

/// A torus centred at the local origin about the local `z` axis, moved by `pose`.
pub fn torus(
    model: &mut Model,
    major: Rational,
    minor: Rational,
    pose: &Pose,
) -> Result<BodyId, Refusal> {
    if major.sign() != Sign::Positive || minor.sign() != Sign::Positive {
        return Err(non_positive());
    }
    if minor >= major {
        return Err(Refusal::new(
            "topology.torus-self-intersecting",
            "minor radius must be less than major radius",
        ));
    }
    let surface = Surface::Torus(Torus {
        center: Point3::origin(),
        axis: unit(2, 1),
        major: major.into(),
        minor: minor.into(),
        ref_dir: unit(0, 1),
    });
    Ok(closed_solid(model, surface, pose))
}

/// A new empty body.
pub fn empty(model: &mut Model) -> BodyId {
    model.add_body(Body {
        kind: BodyKind::Empty,
        lumps: Vec::new(),
        marks: MarkLog::new(),
    })
}

fn non_positive() -> Refusal {
    Refusal::new(
        "topology.non-positive-size",
        "sizes must be positive, except that a cone's top radius may be zero",
    )
}

fn unit(axis: usize, sign: i64) -> Vector3 {
    let [x, y, z] = std::array::from_fn(|i| Rational::from_i64(if i == axis { sign } else { 0 }));
    RVector3::new(x, y, z).to_real()
}

fn at_z(z: Rational) -> Point3 {
    RPoint3::new(Rational::zero(), Rational::zero(), z).to_real()
}

fn solid(model: &mut Model, faces: Vec<FaceId>) -> BodyId {
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    })
}

fn closed_solid(model: &mut Model, surface: Surface, pose: &Pose) -> BodyId {
    let surface = model.geom_mut().add_surface(surface.transformed(pose));
    let face = model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: Vec::new(),
    });
    solid(model, vec![face])
}

fn record(model: &mut Model, kind: FactKind, a: GeomKey, b: GeomKey, origin: &'static str) {
    model.geom_mut().record_fact(Fact {
        kind,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction(origin),
    });
}

fn face(model: &mut Model, surface: SurfaceId, loops: Vec<LoopId>, origin: &'static str) -> FaceId {
    for l in &loops {
        match model.loop_(*l) {
            Loop::Vertex(v) => {
                let point = model.vertex(*v).point;
                record(
                    model,
                    FactKind::Incident,
                    GeomKey::Point(point),
                    GeomKey::Surface(surface),
                    origin,
                );
            }
            Loop::Coedges(cs) => {
                let count = cs.len();
                for i in 0..count {
                    let Loop::Coedges(cs) = model.loop_(*l) else {
                        unreachable!()
                    };
                    let edge = model.edge(model.coedge(cs[i]).edge);
                    let (curve, start, end) = (edge.curve, edge.start, edge.end);
                    record(
                        model,
                        FactKind::Incident,
                        GeomKey::Curve(curve),
                        GeomKey::Surface(surface),
                        origin,
                    );
                    for v in [start, end].into_iter().flatten() {
                        record(
                            model,
                            FactKind::Incident,
                            GeomKey::Point(model.vertex(v).point),
                            GeomKey::Curve(curve),
                            origin,
                        );
                    }
                }
            }
        }
    }
    model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops,
    })
}

fn ring(model: &mut Model, radius: &Rational, z: Rational, pose: &Pose) -> EdgeId {
    let curve = Curve::Circle(Circle {
        center: at_z(z),
        normal: unit(2, 1),
        radius_sq: (radius * radius).into(),
        ref_dir: unit(0, 1),
    });
    let curve = model
        .geom_mut()
        .add_curve(curve.transformed(pose).expect("analytic curve"));
    model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    })
}

fn ring_loop(model: &mut Model, edge: EdgeId, sense: Sense) -> LoopId {
    let coedge = model.add_coedge(Coedge { edge, sense });
    model.add_loop(Loop::Coedges(vec![coedge]))
}

fn cap(model: &mut Model, z: Rational, sign: i64, pose: &Pose) -> SurfaceId {
    let plane = Surface::Plane(Plane {
        origin: at_z(z),
        normal: unit(2, sign),
        ref_dir: unit(0, 1),
    });
    model.geom_mut().add_surface(plane.transformed(pose))
}

fn round_solid(
    model: &mut Model,
    radii: [Rational; 2],
    height: Rational,
    side: Surface,
    pose: &Pose,
    origin: &'static str,
) -> BodyId {
    let bottom_surface = cap(model, Rational::zero(), -1, pose);
    let top_surface = if radii[1].is_zero() {
        None
    } else {
        Some(cap(model, height.clone(), 1, pose))
    };
    let side_surface = model.geom_mut().add_surface(side.transformed(pose));
    let bottom = ring(model, &radii[0], Rational::zero(), pose);
    let bottom_loop = ring_loop(model, bottom, Sense::Reversed);
    let mut faces = vec![face(model, bottom_surface, vec![bottom_loop], origin)];
    let top = top_surface.map(|surface| {
        let edge = ring(model, &radii[1], height.clone(), pose);
        let loop_ = ring_loop(model, edge, Sense::Same);
        faces.push(face(model, surface, vec![loop_], origin));
        record(
            model,
            FactKind::Parallel,
            GeomKey::Surface(bottom_surface),
            GeomKey::Surface(surface),
            origin,
        );
        edge
    });
    let mut side_loops = vec![ring_loop(model, bottom, Sense::Same)];
    side_loops.push(if let Some(top) = top {
        ring_loop(model, top, Sense::Reversed)
    } else {
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(pose.apply_point(&at_z(height))));
        let vertex = model.add_vertex(Vertex { point });
        model.add_loop(Loop::Vertex(vertex))
    });
    faces.push(face(model, side_surface, side_loops, origin));
    solid(model, faces)
}
