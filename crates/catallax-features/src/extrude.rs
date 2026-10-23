//! Closed profile prisms with shared boundary edges and sweep origins.

use catallax_numbers::{Fact, FactKind, FactOrigin, Integer, MarkLog, Rational, Refusal, Sign};
use catallax_recipe::{Origin, Role, StepId, Turn};
use catallax_shapes::{
    Circle, Curve, Cylinder, EllipticCylinder, GeomKey, Line, Plane, PointDef, RPoint3, RVector3,
    Surface,
};
use catallax_topology::{
    Body, BodyId, BodyKind, Coedge, Edge, EdgeId, Face, FaceId, Loop, Lump, Model, Sense, Shell,
    Vertex, VertexId,
};

use crate::Profile;

pub(crate) fn sqrt(n: &Rational) -> Option<Rational> {
    fn root(n: &Integer) -> Option<Integer> {
        if n.sign() == Sign::Negative {
            return None;
        }
        if n.is_zero() {
            return Some(Integer::zero());
        }
        let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        loop {
            let q = n.div_rem_floor(&x)?.0;
            let y = (&x + &q).div_rem_floor(&Integer::from_i64(2))?.0;
            if y >= x {
                return (&x * &x == *n).then_some(x);
            }
            x = y;
        }
    }
    Rational::new(root(n.numerator())?, root(n.denominator())?)
}

pub(crate) fn reverse(sense: Sense) -> Sense {
    match sense {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    }
}

pub(crate) fn face(
    model: &mut Model,
    surface: Surface,
    sense: Sense,
    paths: Vec<Vec<(EdgeId, Sense)>>,
    apex: Option<VertexId>,
    origins: Vec<Origin>,
) -> FaceId {
    let surface = model.geom_mut().add_surface(surface);
    let mut loops = Vec::with_capacity(paths.len() + usize::from(apex.is_some()));
    for path in paths {
        let mut coedges = Vec::with_capacity(path.len());
        for (edge, direction) in path {
            let e = model.edge(edge).clone();
            incident(model, GeomKey::Curve(e.curve), GeomKey::Surface(surface));
            for vertex in [e.start, e.end].into_iter().flatten() {
                incident(
                    model,
                    GeomKey::Point(model.vertex(vertex).point),
                    GeomKey::Curve(e.curve),
                );
            }
            coedges.push(model.add_coedge(Coedge {
                edge,
                sense: direction,
            }));
        }
        loops.push(model.add_loop(Loop::Coedges(coedges)));
    }
    if let Some(vertex) = apex {
        incident(
            model,
            GeomKey::Point(model.vertex(vertex).point),
            GeomKey::Surface(surface),
        );
        loops.push(model.add_loop(Loop::Vertex(vertex)));
    }
    model.add_named_face(
        Face {
            surface,
            sense,
            loops,
        },
        origins,
    )
}

fn incident(model: &mut Model, a: GeomKey, b: GeomKey) {
    model.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("features.sweep"),
    });
}

/// Sweeps a profile by a signed distance. Directions use model-space coordinates.
/// The profile's normal is the default; changing the sign preserves its origins.
pub fn extrude(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    distance: &Rational,
    direction: Option<&RVector3>,
) -> Result<BodyId, Refusal> {
    if distance.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-distance",
            "extrusion distance must not be zero",
        ));
    }
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    let direction = direction.unwrap_or(&normal);
    if direction.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-direction",
            "extrusion direction must not be zero",
        ));
    }
    let height = direction.dot(&normal);
    if height.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-thickness",
            "extrusion direction lies in the sketch plane",
        ));
    }
    if profile.is_empty() {
        return Ok(catallax_topology::empty(model));
    }
    let length = sqrt(&direction.norm_squared()).ok_or_else(|| {
        Refusal::new(
            "features.unsupported-algebraic-placement",
            "the analytic geometry store needs an exact algebraic placement for this direction",
        )
    })?;
    let vector = direction.scaled(&(distance / &length));
    extrude_vector(model, profile, step, &vector)
}

pub(crate) fn extrude_vector(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    vector: &RVector3,
) -> Result<BodyId, Refusal> {
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    if vector.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-direction",
            "extrusion vector is zero",
        ));
    }
    let height = vector.dot(&normal);
    if height.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-thickness",
            "extrusion direction lies in the sketch plane",
        ));
    }
    let positive = height.sign() == Sign::Positive;
    let forward = if positive {
        Sense::Same
    } else {
        Sense::Reversed
    };
    let zero = RVector3::zero();
    let mut lumps = Vec::new();
    for (outer, boundary) in profile
        .loops
        .iter()
        .enumerate()
        .filter(|(_, b)| b.depth % 2 == 0)
    {
        let indices: Vec<_> = std::iter::once(outer)
            .chain(profile.loops.iter().enumerate().filter_map(|(i, b)| {
                (b.parent == Some(outer) && b.depth == boundary.depth + 1).then_some(i)
            }))
            .collect();
        let mut bottom = Vec::new();
        let mut top = Vec::new();
        let mut sides = Vec::new();
        for index in indices {
            let boundary = &profile.loops[index];
            let segments = &boundary.segments;
            let vertices: Vec<_> = if segments[0].full {
                Vec::new()
            } else {
                segments
                    .iter()
                    .map(|s| {
                        [zero.clone(), vector.clone()].map(|offset| {
                            let point = model.geom_mut().add_point(PointDef::Explicit(
                                profile.world(&s.start, &offset).into(),
                            ));
                            model.add_vertex(Vertex { point })
                        })
                    })
                    .collect()
            };
            let vertical: Vec<_> = vertices
                .iter()
                .enumerate()
                .map(|(i, pair)| {
                    let curve = model.geom_mut().add_curve(Curve::Line(Line {
                        origin: profile.world(&segments[i].start, &zero).into(),
                        dir: vector.clone().into(),
                    }));
                    model.add_edge(Edge {
                        curve,
                        sense: Sense::Same,
                        start: Some(pair[0]),
                        end: Some(pair[1]),
                    })
                })
                .collect();
            let mut start_path = Vec::new();
            let mut end_path = Vec::new();
            for (i, s) in segments.iter().enumerate() {
                let mut edges = Vec::with_capacity(2);
                for (layer, offset) in [&zero, vector].into_iter().enumerate() {
                    let curve = if let Some((c, radius_sq, _)) = &s.circle {
                        Curve::Circle(Circle {
                            center: profile.world(c, offset).into(),
                            normal: normal.clone().into(),
                            radius_sq: radius_sq.clone().into(),
                            ref_dir: profile.pose.apply_vector(
                                &RVector3::new(
                                    &s.start[0] - &c[0],
                                    &s.start[1] - &c[1],
                                    Rational::zero(),
                                )
                                .to_real(),
                            ),
                        })
                    } else {
                        Curve::Line(Line {
                            origin: profile.world(&s.start, offset).into(),
                            dir: profile
                                .world(&s.end, offset)
                                .minus(&profile.world(&s.start, offset))
                                .into(),
                        })
                    };
                    let curve = model.geom_mut().add_curve(curve);
                    let sense = if matches!(&s.circle, Some((_, _, Turn::Cw))) {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    };
                    let (start, end) = if s.full {
                        (None, None)
                    } else {
                        (
                            Some(vertices[i][layer]),
                            Some(vertices[(i + 1) % segments.len()][layer]),
                        )
                    };
                    edges.push(model.add_edge(Edge {
                        curve,
                        sense,
                        start,
                        end,
                    }));
                }
                start_path.push((edges[0], reverse(forward)));
                end_path.push((edges[1], forward));
                let mut path = vec![(edges[0], Sense::Same)];
                let paths = if s.full {
                    vec![path, vec![(edges[1], Sense::Reversed)]]
                } else {
                    path.extend([
                        (vertical[(i + 1) % segments.len()], Sense::Same),
                        (edges[1], Sense::Reversed),
                        (vertical[i], Sense::Reversed),
                    ]);
                    vec![path]
                };
                let paths = if positive {
                    paths
                } else {
                    paths
                        .into_iter()
                        .map(|p| p.into_iter().rev().map(|(e, s)| (e, reverse(s))).collect())
                        .collect()
                };
                let (surface, sense) = if let Some((c, radius_sq, turn)) = &s.circle {
                    let major_dir = normal.cross(vector);
                    let surface = if major_dir.is_zero() {
                        Surface::Cylinder(Cylinder {
                            origin: profile.world(c, &zero).into(),
                            axis: normal.clone().into(),
                            radius_sq: radius_sq.clone().into(),
                            ref_dir: profile.pose.apply_vector(
                                &RVector3::new(
                                    &s.start[0] - &c[0],
                                    &s.start[1] - &c[1],
                                    Rational::zero(),
                                )
                                .to_real(),
                            ),
                        })
                    } else {
                        // The profile normal is a rigid image of the unit z
                        // axis. Nonzero height and major_dir prove 0 < cos² < 1.
                        // Keep the original circular end curves and oblique caps.
                        Surface::EllipticCylinder(EllipticCylinder {
                            origin: profile.world(c, &zero).into(),
                            axis: vector.clone().into(),
                            major_dir: major_dir.into(),
                            major_sq: radius_sq.clone().into(),
                            minor_sq: (radius_sq * &(&height * &height) / vector.norm_squared())
                                .into(),
                        })
                    };
                    (
                        surface,
                        if *turn == Turn::Ccw {
                            Sense::Same
                        } else {
                            Sense::Reversed
                        },
                    )
                } else {
                    let tangent = profile
                        .world(&s.end, &zero)
                        .minus(&profile.world(&s.start, &zero));
                    let outward = tangent
                        .cross(vector)
                        .scaled(&Rational::from_i64(if positive { 1 } else { -1 }));
                    (
                        Surface::Plane(Plane {
                            origin: profile.world(&s.start, &zero).into(),
                            normal: outward.into(),
                            ref_dir: vector.clone().into(),
                        }),
                        Sense::Same,
                    )
                };
                let origins = s
                    .ids
                    .iter()
                    .map(|segment| Origin {
                        step: step.clone(),
                        role: Role::Side {
                            segment: segment.clone(),
                            path: None,
                        },
                    })
                    .collect();
                sides.push(face(model, surface, sense, paths, None, origins));
            }
            if positive {
                start_path.reverse();
            } else {
                end_path.reverse();
            }
            bottom.push(start_path);
            top.push(end_path);
        }
        let reference = profile.pose.apply_vector(
            &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
        );
        let origin = profile.pose.apply_point(&RPoint3::origin().to_real());
        let start = face(
            model,
            Surface::Plane(Plane {
                origin: origin.clone(),
                normal: normal.clone().into(),
                ref_dir: reference.clone(),
            }),
            reverse(forward),
            bottom,
            None,
            vec![Origin {
                step: step.clone(),
                role: Role::Start,
            }],
        );
        let end = face(
            model,
            Surface::Plane(Plane {
                origin: origin.plus(&vector.to_real()),
                normal: normal.clone().into(),
                ref_dir: reference,
            }),
            forward,
            top,
            None,
            vec![Origin {
                step: step.clone(),
                role: Role::End,
            }],
        );
        let mut faces = vec![start, end];
        faces.extend(sides);
        let shell = model.add_shell(Shell { faces });
        lumps.push(model.add_lump(Lump {
            shells: vec![shell],
        }));
    }
    Ok(model.add_body(Body {
        kind: if lumps.is_empty() {
            BodyKind::Empty
        } else {
            BodyKind::Solid
        },
        lumps,
        marks: MarkLog::new(),
    }))
}
