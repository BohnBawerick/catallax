//! Revolutions of planar line and circle boundaries about in-plane axes.

use catallax_numbers::{Decider, MarkLog, Rational, Real, Refusal, Sign};
use catallax_recipe::{Axis, Origin, Role, StepId, Turn};
use catallax_shapes::{
    Circle, Cone, Curve, Cylinder, GeomKey, Line, Plane, Point3, PointDef, RPoint3, RVector3,
    Sphere, Surface, Torus,
};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Loop, Lump, Model, Sense, Shell, Vertex};

use crate::{
    Profile,
    extrude::{face, reverse, sqrt},
    profile::{Arc, sign},
};

fn bad(code: &'static str, detail: &str) -> Refusal {
    Refusal::new(code, detail)
}
fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    RVector3::new(
        v.0[0].value().clone(),
        v.0[1].value().clone(),
        v.0[2].value().clone(),
    )
}
fn qpoint(p: &[Rational; 2]) -> RPoint3 {
    RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero())
}
fn unit_z() -> RVector3 {
    RVector3::new(Rational::zero(), Rational::zero(), Rational::one())
}
fn rational_radius(radius_sq: &Rational) -> Result<Rational, Refusal> {
    sqrt(radius_sq).ok_or_else(|| {
        bad(
            "features.algebraic-radius",
            "the analytic surface store needs an exact squared-radius surface definition",
        )
    })
}

struct Frame {
    origin: RPoint3,
    axis: RVector3,
    norm: Rational,
}
impl Frame {
    fn project(&self, p: &RPoint3) -> RPoint3 {
        self.origin.plus(
            &self
                .axis
                .scaled(&(p.minus(&self.origin).dot(&self.axis) / self.norm.clone())),
        )
    }
    fn signed(&self, p: &RPoint3) -> Rational {
        self.axis.cross(&p.minus(&self.origin)).z
    }
    fn rotate(&self, p: &RPoint3, cosine: &Rational, sine: &Rational) -> Result<RPoint3, Refusal> {
        let projection = self.project(p);
        let radial = p.minus(&projection);
        let transverse = if sine.is_zero() {
            RVector3::zero()
        } else {
            let length = sqrt(&self.norm).ok_or_else(|| {
                bad(
                    "features.algebraic-frame",
                    "the analytic geometry store needs an algebraic rotation frame",
                )
            })?;
            self.axis.cross(&radial).scaled(&(sine / &length))
        };
        Ok(projection.plus(&radial.scaled(cosine)).plus(&transverse))
    }
}

fn carrier(segment: &Arc, frame: &Frame) -> Result<Option<Surface>, Refusal> {
    let start = qpoint(&segment.start);
    let end = qpoint(&segment.end);
    let radial = start.minus(&frame.project(&start));
    if let Some((centre, radius_sq, _)) = &segment.circle {
        let centre = qpoint(centre);
        let projection = frame.project(&centre);
        let offset = centre.minus(&projection);
        let radius = rational_radius(radius_sq)?;
        return Ok(Some(if offset.is_zero() {
            Surface::Sphere(Sphere {
                center: centre,
                radius,
                axis: frame.axis.clone(),
                ref_dir: frame.axis.cross(&unit_z()),
            })
        } else {
            let major = rational_radius(&offset.norm_squared())?;
            Surface::Torus(Torus {
                center: projection,
                axis: frame.axis.clone(),
                major,
                minor: radius,
                ref_dir: offset,
            })
        }));
    }
    let tangent = end.minus(&start);
    if radial.is_zero() && end.minus(&frame.project(&end)).is_zero() {
        return Ok(None);
    }
    let axial = tangent.dot(&frame.axis);
    if axial.is_zero() {
        let reference = if radial.is_zero() {
            end.minus(&frame.project(&end))
        } else {
            radial
        };
        return Ok(Some(Surface::Plane(Plane {
            origin: frame.project(&start),
            normal: frame.axis.clone(),
            ref_dir: reference,
        })));
    }
    if tangent.cross(&frame.axis).is_zero() {
        return Ok(Some(Surface::Cylinder(Cylinder {
            origin: frame.project(&start),
            axis: frame.axis.clone(),
            radius: rational_radius(&radial.norm_squared())?,
            ref_dir: radial,
        })));
    }
    let initial = frame.signed(&start);
    let change = frame.axis.cross(&tangent).z;
    let apex = start.plus(&tangent.scaled(&(-initial / change.clone())));
    let direction = if start.minus(&apex).dot(&frame.axis).sign() == Sign::Positive
        || (start == apex && end.minus(&apex).dot(&frame.axis).sign() == Sign::Positive)
    {
        frame.axis.clone()
    } else {
        frame.axis.scaled(&Rational::from_i64(-1))
    };
    Ok(Some(Surface::Cone(Cone {
        apex,
        axis: direction,
        slope: (change / axial).abs(),
        ref_dir: if radial.is_zero() {
            end.minus(&frame.project(&end))
        } else {
            radial
        },
    })))
}

fn material_side(boundary: &[Arc], frame: &Frame, dec: &mut Decider) -> Result<Sign, Refusal> {
    let mut positive = false;
    let mut negative = false;
    for segment in boundary {
        let mut values = vec![
            Real::from(frame.signed(&qpoint(&segment.start))),
            Real::from(frame.signed(&qpoint(&segment.end))),
        ];
        if let Some((centre, radius_sq, _)) = &segment.circle {
            if segment.full {
                let distance = frame.signed(&qpoint(centre));
                if &distance * &distance == radius_sq * &frame.norm {
                    return Err(bad(
                        "boolean.self-touching-shell",
                        "a circle tangent to the revolution axis makes a pinched shell",
                    ));
                }
            }
            let radius = (Real::from(radius_sq.clone()) / Real::from(frame.norm.clone())).sqrt();
            for flip in [-1, 1] {
                let scale = &Real::from_i64(flip) * &radius;
                let point = [
                    Real::from(centre[0].clone()) - &Real::from(frame.axis.y.clone()) * &scale,
                    Real::from(centre[1].clone()) + &Real::from(frame.axis.x.clone()) * &scale,
                ];
                if segment.contains(&point, dec)? {
                    let x = &point[0] - &Real::from(frame.origin.x.clone());
                    let y = &point[1] - &Real::from(frame.origin.y.clone());
                    values.push(
                        &Real::from(frame.axis.x.clone()) * &y
                            - &Real::from(frame.axis.y.clone()) * &x,
                    );
                }
            }
        }
        for value in values {
            match sign(&value, dec)? {
                Sign::Positive => positive = true,
                Sign::Negative => negative = true,
                Sign::Zero => {}
            }
        }
    }
    match (positive, negative) {
        (true, true) => Err(bad(
            "revolve.axis-interior",
            "revolution axis meets the interior of the profile",
        )),
        (true, false) => Ok(Sign::Positive),
        (false, true) => Ok(Sign::Negative),
        _ => Err(bad(
            "revolve.axis-interior",
            "profile has no area away from the revolution axis",
        )),
    }
}

/// Turns a profile about an in-plane model-space axis, preserving segment origins.
pub fn revolve(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    axis: &Axis,
    angle: &Rational,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    if angle.sign() != Sign::Positive || *angle > Rational::from_i64(360) {
        return Err(bad(
            "revolve.angle",
            "revolution angle must be in (0, 360] degrees",
        ));
    }
    let world_axis = vector(&axis.direction);
    if world_axis.is_zero() {
        return Err(bad(
            "revolve.zero-axis",
            "revolution axis direction is zero",
        ));
    }
    let inverse = profile.pose.inverse();
    let origin = inverse.apply_point(&RPoint3::new(
        axis.point.0[0].value().clone(),
        axis.point.0[1].value().clone(),
        axis.point.0[2].value().clone(),
    ));
    let direction = inverse.apply_vector(&world_axis);
    if !origin.z.is_zero() || !direction.z.is_zero() {
        return Err(bad(
            "revolve.axis-plane",
            "revolution axis must lie in the sketch plane",
        ));
    }
    let frame = Frame {
        origin,
        norm: direction.norm_squared(),
        axis: direction,
    };
    if profile.is_empty() {
        return Ok(catallax_topology::empty(model));
    }
    let full = *angle == Rational::from_i64(360);
    let (cosine, sine) = match angle.to_string().as_str() {
        "90" => (Rational::zero(), Rational::one()),
        "180" => (-Rational::one(), Rational::zero()),
        "270" => (Rational::zero(), -Rational::one()),
        "360" => (Rational::one(), Rational::zero()),
        _ => {
            return Err(bad(
                "features.algebraic-angle",
                "the analytic geometry store needs exact algebraic cap planes for this rational-degree angle",
            ));
        }
    };
    let normal = profile.pose.apply_vector(&unit_z());
    let mut bodies = Vec::new();
    for (outer, boundary) in profile
        .loops
        .iter()
        .enumerate()
        .filter(|(_, b)| b.depth % 2 == 0)
    {
        let side = material_side(&boundary.segments, &frame, dec)?;
        let forward = if side == Sign::Positive {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let indices: Vec<_> = std::iter::once(outer)
            .chain(profile.loops.iter().enumerate().filter_map(|(i, b)| {
                (b.parent == Some(outer) && b.depth == boundary.depth + 1).then_some(i)
            }))
            .collect();
        let mut faces = Vec::new();
        let mut start_paths = Vec::new();
        let mut end_paths = Vec::new();
        for index in indices {
            let segments = &profile.loops[index].segments;
            let mut vertices = Vec::new();
            let mut orbit = Vec::new();
            for s in segments {
                let p = qpoint(&s.start);
                let on_axis = frame.signed(&p).is_zero();
                let first = if s.full || (full && !on_axis) {
                    None
                } else {
                    let point = model
                        .geom_mut()
                        .add_point(PointDef::Explicit(profile.pose.apply_point(&p)));
                    Some(model.add_vertex(Vertex { point }))
                };
                let second = if s.full || on_axis || full {
                    first
                } else {
                    let point = model.geom_mut().add_point(PointDef::Explicit(
                        profile.pose.apply_point(&frame.rotate(&p, &cosine, &sine)?),
                    ));
                    Some(model.add_vertex(Vertex { point }))
                };
                vertices.push([first, second]);
                orbit.push(if on_axis || s.full {
                    None
                } else {
                    let centre = frame.project(&p);
                    let radial = p.minus(&centre);
                    let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                        center: profile.pose.apply_point(&centre),
                        normal: world_axis.clone(),
                        radius_sq: radial.norm_squared(),
                        ref_dir: profile.pose.apply_vector(&radial),
                    }));
                    Some(model.add_edge(Edge {
                        curve,
                        sense: Sense::Same,
                        start: first,
                        end: second,
                    }))
                });
            }
            let mut start_path = Vec::new();
            let mut end_path = Vec::new();
            for (i, s) in segments.iter().enumerate() {
                let next = (i + 1) % segments.len();
                let surface = carrier(s, &frame)?;
                let mut profile_edges = Vec::new();
                if !full {
                    for layer in 0..2 {
                        if layer == 1 && surface.is_none() {
                            profile_edges.push(profile_edges[0]);
                            continue;
                        }
                        let curve = if let Some((c, radius_sq, _)) = &s.circle {
                            let centre = qpoint(c);
                            let centre = if layer == 0 {
                                centre
                            } else {
                                frame.rotate(&centre, &cosine, &sine)?
                            };
                            let start = if layer == 0 {
                                qpoint(&s.start)
                            } else {
                                frame.rotate(&qpoint(&s.start), &cosine, &sine)?
                            };
                            let rotated_normal = if layer == 0 {
                                unit_z()
                            } else {
                                frame
                                    .rotate(&frame.origin.plus(&unit_z()), &cosine, &sine)?
                                    .minus(&frame.origin)
                            };
                            Curve::Circle(Circle {
                                center: profile.pose.apply_point(&centre),
                                normal: profile.pose.apply_vector(&rotated_normal),
                                radius_sq: radius_sq.clone(),
                                ref_dir: profile.pose.apply_vector(&start.minus(&centre)),
                            })
                        } else {
                            let start = if layer == 0 {
                                qpoint(&s.start)
                            } else {
                                frame.rotate(&qpoint(&s.start), &cosine, &sine)?
                            };
                            let end = if layer == 0 {
                                qpoint(&s.end)
                            } else {
                                frame.rotate(&qpoint(&s.end), &cosine, &sine)?
                            };
                            Curve::Line(Line {
                                origin: profile.pose.apply_point(&start),
                                dir: profile.pose.apply_vector(&end.minus(&start)),
                            })
                        };
                        let curve = model.geom_mut().add_curve(curve);
                        let edge_sense = if matches!(s.circle, Some((_, _, Turn::Cw))) {
                            Sense::Reversed
                        } else {
                            Sense::Same
                        };
                        profile_edges.push(model.add_edge(Edge {
                            curve,
                            sense: edge_sense,
                            start: vertices[i][layer],
                            end: vertices[next][layer],
                        }));
                    }
                    start_path.push((profile_edges[0], reverse(forward)));
                    end_path.push((profile_edges[1], forward));
                }
                let Some(local_surface) = surface else {
                    continue;
                };
                let surface = local_surface.transformed(&profile.pose);
                let mid = s.midpoint();
                let point = profile.pose.apply_real_point(&Point3::new(
                    mid[0].clone(),
                    mid[1].clone(),
                    Real::zero(),
                ));
                let outward = if let Some((c, _, turn)) = &s.circle {
                    let radial = Point3::new(mid[0].clone(), mid[1].clone(), Real::zero())
                        .minus(&qpoint(c).to_real());
                    profile.pose.rotation.apply_real(
                        &radial.scaled(&Real::from_i64(if *turn == Turn::Ccw { 1 } else { -1 })),
                    )
                } else {
                    profile
                        .pose
                        .apply_vector(&qpoint(&s.end).minus(&qpoint(&s.start)).cross(&unit_z()))
                        .to_real()
                };
                let sense = match sign(&surface.gradient(&point).dot(&outward), dec)? {
                    Sign::Positive => Sense::Same,
                    Sign::Negative => Sense::Reversed,
                    Sign::Zero => {
                        return Err(bad(
                            "features.singular-face",
                            "revolved boundary has a singular normal",
                        ));
                    }
                };
                let mut paths = Vec::new();
                if full {
                    if let Some(edge) = orbit[next] {
                        paths.push(vec![(edge, forward)]);
                    }
                    if let Some(edge) = orbit[i] {
                        paths.push(vec![(edge, reverse(forward))]);
                    }
                } else if s.full {
                    paths.push(vec![(profile_edges[0], forward)]);
                    paths.push(vec![(profile_edges[1], reverse(forward))]);
                } else {
                    let mut path = vec![(profile_edges[0], Sense::Same)];
                    if let Some(edge) = orbit[next] {
                        path.push((edge, Sense::Same));
                    }
                    path.push((profile_edges[1], Sense::Reversed));
                    if let Some(edge) = orbit[i] {
                        path.push((edge, Sense::Reversed));
                    }
                    if forward == Sense::Reversed {
                        path.reverse();
                        for (_, s) in &mut path {
                            *s = reverse(*s);
                        }
                    }
                    paths.push(path);
                }
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
                let stored = face(model, surface, sense, paths, origins);
                // A full cone retains its singular apex as a point-only loop.
                let stored = if full && matches!(local_surface, Surface::Cone(_)) {
                    let apex = vertices[i][0].or(vertices[next][0]);
                    if let Some(apex) = apex {
                        let point_loop = model.add_loop(Loop::Vertex(apex));
                        let mut f = model.face(stored).clone();
                        f.loops.push(point_loop);
                        let origins = model.face_origins(stored).to_vec();
                        let surface = f.surface;
                        let point = model.vertex(apex).point;
                        model.geom_mut().record_fact(catallax_numbers::Fact {
                            kind: catallax_numbers::FactKind::Incident,
                            a: GeomKey::Point(point),
                            b: GeomKey::Surface(surface),
                            at: None,
                            origin: catallax_numbers::FactOrigin::Construction("features.revolve"),
                        });
                        model.add_named_face(f, origins)
                    } else {
                        stored
                    }
                } else {
                    stored
                };
                faces.push(stored);
            }
            if forward == Sense::Same {
                start_path.reverse();
            } else {
                end_path.reverse();
            }
            start_paths.push(start_path);
            end_paths.push(end_path);
        }
        if !full {
            let origin = profile.pose.apply_point(&RPoint3::origin());
            let reference = profile.pose.apply_vector(&frame.axis);
            faces.push(face(
                model,
                Surface::Plane(Plane {
                    origin,
                    normal: normal.clone(),
                    ref_dir: reference.clone(),
                }),
                reverse(forward),
                start_paths,
                vec![Origin {
                    step: step.clone(),
                    role: Role::Start,
                }],
            ));
            let origin =
                profile
                    .pose
                    .apply_point(&frame.rotate(&RPoint3::origin(), &cosine, &sine)?);
            let end_normal = frame
                .rotate(&frame.origin.plus(&unit_z()), &cosine, &sine)?
                .minus(&frame.origin);
            faces.push(face(
                model,
                Surface::Plane(Plane {
                    origin,
                    normal: profile.pose.apply_vector(&end_normal),
                    ref_dir: reference,
                }),
                forward,
                end_paths,
                vec![Origin {
                    step: step.clone(),
                    role: Role::End,
                }],
            ));
        }
        let shell = model.add_shell(Shell { faces });
        let lump = model.add_lump(Lump {
            shells: vec![shell],
        });
        bodies.push(model.add_body(Body {
            kind: BodyKind::Solid,
            lumps: vec![lump],
            marks: MarkLog::new(),
        }));
    }
    let mut bodies = bodies.into_iter();
    let Some(mut result) = bodies.next() else {
        return Ok(catallax_topology::empty(model));
    };
    for body in bodies {
        result =
            catallax_boolean::boolean(model, result, body, catallax_boolean::BoolOp::Unite, dec)
                .map_err(|e| e.refusal)?
                .body;
    }
    Ok(result)
}
