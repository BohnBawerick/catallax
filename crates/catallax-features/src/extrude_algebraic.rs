//! Finite signed sweeps with exact-real offsets and rational analytic carriers.
//!
//! Keep the rational construction path unchanged. Carrier orientation and the
//! elliptic reduction use the input direction, so no radical is introduced into
//! a side surface merely to describe the finite distance between its caps.

use super::{face, reverse};
use crate::Profile;
use catallax_numbers::{MarkLog, Rational, Real, Refusal, Sign};
use catallax_recipe::{Origin, Role, StepId, Turn};
use catallax_shapes::{
    Circle, Curve, Cylinder, EllipticCylinder, Line, Plane, Point3, PointDef, RPoint3, RVector3,
    Surface, Vector3,
};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Lump, Model, Sense, Shell, Vertex};

fn world(profile: &Profile, p: &[Rational; 2], offset: &Vector3) -> Point3 {
    profile.world(p, &RVector3::zero()).to_real().plus(offset)
}

pub(super) fn extrude(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    distance: &Rational,
    direction: &RVector3,
) -> Result<BodyId, Refusal> {
    // The caller has proved nonzero distance, direction and normal height.
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    let height = direction.dot(&normal);
    let positive = (&height * distance).sign() == Sign::Positive;
    let forward = if positive {
        Sense::Same
    } else {
        Sense::Reversed
    };
    let scale = &Real::from(distance.clone()) / &Real::from(direction.norm_squared()).sqrt();
    let vector = direction.to_real().scaled(&scale);
    let zero = Vector3::zero();
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
                            let point = model
                                .geom_mut()
                                .add_point(PointDef::Explicit(world(profile, &s.start, &offset)));
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
                        origin: world(profile, &segments[i].start, &zero),
                        dir: vector.clone(),
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
                for (layer, offset) in [&zero, &vector].into_iter().enumerate() {
                    let curve = if let Some((c, radius_sq, _)) = &s.circle {
                        Curve::Circle(Circle {
                            center: world(profile, c, offset),
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
                            origin: world(profile, &s.start, offset),
                            dir: world(profile, &s.end, offset)
                                .minus(&world(profile, &s.start, offset)),
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
                    let major_dir = normal.cross(direction);
                    let surface = if major_dir.is_zero() {
                        Surface::Cylinder(Cylinder {
                            origin: world(profile, c, &zero),
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
                            origin: world(profile, c, &zero),
                            axis: direction.clone().into(),
                            major_dir: major_dir.into(),
                            major_sq: radius_sq.clone().into(),
                            minor_sq: (radius_sq * &(&height * &height) / direction.norm_squared())
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
                        .world(&s.end, &RVector3::zero())
                        .minus(&profile.world(&s.start, &RVector3::zero()));
                    let outward = tangent.cross(direction).scaled(&Rational::from_i64(
                        if height.sign() == Sign::Positive {
                            1
                        } else {
                            -1
                        },
                    ));
                    (
                        Surface::Plane(Plane {
                            origin: world(profile, &s.start, &zero),
                            normal: outward.into(),
                            ref_dir: direction.clone().into(),
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
            [Role::Start]
                .iter()
                .map(|role| Origin {
                    step: step.clone(),
                    role: role.clone(),
                })
                .collect(),
        );
        let end = face(
            model,
            Surface::Plane(Plane {
                origin: origin.plus(&vector),
                normal: normal.clone().into(),
                ref_dir: reference,
            }),
            forward,
            top,
            None,
            [Role::End]
                .iter()
                .map(|role| Origin {
                    step: step.clone(),
                    role: role.clone(),
                })
                .collect(),
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
