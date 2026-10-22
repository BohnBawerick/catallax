//! Revolutions of planar line and circle boundaries about in-plane axes.

use catallax_numbers::{Decider, MarkLog, Rational, Real, Refusal, Sign};
use catallax_recipe::{Axis, Origin, Role, StepId, Turn};
use catallax_shapes::{
    Circle, Cone, Curve, Cylinder, Line, Plane, Point3, PointDef, RPoint3, RVector3, Sphere,
    Surface, Torus,
};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Lump, Model, Sense, Shell, Vertex};

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
fn rational_radius(radius_sq: &Rational, code: &'static str) -> Result<Rational, Refusal> {
    sqrt(radius_sq).ok_or_else(|| {
        bad(
            code,
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
                    "features.unsupported-algebraic-placement",
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
        let radius = rational_radius(
            radius_sq,
            if offset.is_zero() {
                "features.unsupported-algebraic-sphere-radius"
            } else {
                "features.unsupported-algebraic-torus-radius"
            },
        )?;
        return Ok(Some(if offset.is_zero() {
            Surface::Sphere(Sphere {
                center: centre.into(),
                radius: radius.into(),
                axis: frame.axis.clone().into(),
                ref_dir: frame.axis.cross(&unit_z()).into(),
            })
        } else {
            let major = rational_radius(
                &offset.norm_squared(),
                "features.unsupported-algebraic-torus-radius",
            )?;
            if major <= radius {
                return Err(bad(
                    "features.unsupported-spindle-torus",
                    "revolving this arc needs a horn or spindle torus sheet, not the ring-torus representation",
                ));
            }
            Surface::Torus(Torus {
                center: projection.into(),
                axis: frame.axis.clone().into(),
                major: major.into(),
                minor: radius.into(),
                ref_dir: offset.into(),
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
            origin: frame.project(&start).into(),
            normal: frame.axis.clone().into(),
            ref_dir: reference.into(),
        })));
    }
    if tangent.cross(&frame.axis).is_zero() {
        return Ok(Some(Surface::Cylinder(Cylinder {
            origin: frame.project(&start).into(),
            axis: frame.axis.clone().into(),
            radius: rational_radius(
                &radial.norm_squared(),
                "features.unsupported-algebraic-cylinder-radius",
            )?
            .into(),
            ref_dir: radial.into(),
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
        apex: apex.into(),
        axis: direction.into(),
        slope: (change / axial).abs().into(),
        ref_dir: if radial.is_zero() {
            end.minus(&frame.project(&end)).into()
        } else {
            radial.into()
        },
    })))
}

fn full_turn_axis_contacts(
    boundary: &[Arc],
    frame: &Frame,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let on_axis = |s: &Arc| {
        s.circle.is_none()
            && frame.signed(&qpoint(&s.start)).is_zero()
            && frame.signed(&qpoint(&s.end)).is_zero()
    };
    for (i, segment) in boundary.iter().enumerate() {
        let previous = &boundary[(i + boundary.len() - 1) % boundary.len()];
        if !segment.full
            && frame.signed(&qpoint(&segment.start)).is_zero()
            && !on_axis(segment)
            && !on_axis(previous)
        {
            return Err(bad(
                "boolean.self-touching-shell",
                "an isolated profile vertex on the axis pinches the full-turn shell",
            ));
        }
        if let Some((centre, radius_sq, _)) = &segment.circle {
            let distance = frame.signed(&qpoint(centre));
            if &distance * &distance == radius_sq * &frame.norm {
                let touch = frame.project(&qpoint(centre));
                let point = [touch.x, touch.y];
                if (segment.full || (point != segment.start && point != segment.end))
                    && segment.contains(&point.each_ref().map(|r| Real::from(r.clone())), dec)?
                {
                    return Err(bad(
                        "boolean.self-touching-shell",
                        "an interior arc tangency to the axis pinches the full-turn shell",
                    ));
                }
            }
        }
    }
    Ok(())
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
    let origin = crate::profile::rational_point(
        &inverse,
        &RPoint3::new(
            axis.point.0[0].value().clone(),
            axis.point.0[1].value().clone(),
            axis.point.0[2].value().clone(),
        ),
    );
    let direction = crate::profile::rational_vector(&inverse, &world_axis);
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
    let outer_sides: Vec<_> = profile
        .loops
        .iter()
        .filter(|b| b.depth % 2 == 0)
        .map(|b| {
            let side = material_side(&b.segments, &frame, dec)?;
            if full {
                full_turn_axis_contacts(&b.segments, &frame, dec)?;
            }
            Ok(side)
        })
        .collect::<Result<_, Refusal>>()?;
    let (cosine, sine) = if *angle == Rational::from_i64(90) {
        (Rational::zero(), Rational::one())
    } else if *angle == Rational::from_i64(180) {
        (-Rational::one(), Rational::zero())
    } else if *angle == Rational::from_i64(270) {
        (Rational::zero(), -Rational::one())
    } else if full {
        (Rational::one(), Rational::zero())
    } else {
        return Err(bad(
            "features.unsupported-algebraic-revolve-frame",
            "this rational-degree angle needs algebraic rotated cap and trimming-curve definitions",
        ));
    };
    let half_turn = !full && sine.is_zero();
    let normal = profile.pose.apply_vector(&unit_z().to_real());
    let mut bodies = Vec::new();
    for ((outer, boundary), side) in profile
        .loops
        .iter()
        .enumerate()
        .filter(|(_, b)| b.depth % 2 == 0)
        .zip(outer_sides)
    {
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
        let mut shells = Vec::new();
        let mut start_paths = Vec::new();
        let mut end_paths = Vec::new();
        let mut merge_caps = false;
        for index in indices {
            let segments = &profile.loops[index].segments;
            let mut vertices = Vec::new();
            let mut orbit = Vec::new();
            for s in segments {
                let p = qpoint(&s.start);
                let on_axis = frame.signed(&p).is_zero();
                let first = if s.full || full {
                    None
                } else {
                    let point = model
                        .geom_mut()
                        .add_point(PointDef::Explicit(profile.pose.apply_point(&p.to_real())));
                    Some(model.add_vertex(Vertex { point }))
                };
                let second = if s.full || on_axis || full {
                    first
                } else {
                    let point = model.geom_mut().add_point(PointDef::Explicit(
                        profile
                            .pose
                            .apply_point(&frame.rotate(&p, &cosine, &sine)?.into()),
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
                        center: profile.pose.apply_point(&centre.to_real()),
                        normal: world_axis.clone().into(),
                        radius_sq: radial.norm_squared().into(),
                        ref_dir: profile.pose.apply_vector(&radial.to_real()),
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
                merge_caps |= half_turn && surface.is_none();
                let mut profile_edges = Vec::new();
                if !full {
                    for (layer, start_vertex) in vertices[i].iter().copied().enumerate() {
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
                                center: profile.pose.apply_point(&centre.to_real()),
                                normal: profile.pose.apply_vector(&rotated_normal.to_real()),
                                radius_sq: radius_sq.clone().into(),
                                ref_dir: profile.pose.apply_vector(&start.minus(&centre).to_real()),
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
                                origin: profile.pose.apply_point(&start.to_real()),
                                dir: profile.pose.apply_vector(&end.minus(&start).to_real()),
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
                            start: start_vertex,
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
                    profile.pose.apply_vector(
                        &qpoint(&s.end)
                            .minus(&qpoint(&s.start))
                            .cross(&unit_z())
                            .to_real(),
                    )
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
                let apex = if full {
                    if let Surface::Cone(cone) = &local_surface
                        && {
                            let [x, y, z] = cone.apex.rational_components()
                                .expect("the local meridian is constructed from rational profile coordinates");
                            z.is_zero()
                                && [&s.start, &s.end]
                                    .into_iter()
                                    .any(|p| x == &p[0] && y == &p[1])
                        }
                    {
                        let point = model
                            .geom_mut()
                            .add_point(PointDef::Explicit(profile.pose.apply_point(&cone.apex)));
                        Some(model.add_vertex(Vertex { point }))
                    } else {
                        None
                    }
                } else {
                    None
                };
                let stored = face(model, surface, sense, paths, apex, origins);
                faces.push(stored);
            }
            if forward == Sense::Same {
                start_path.reverse();
            } else {
                end_path.reverse();
            }
            start_paths.push(start_path);
            end_paths.push(end_path);
            if full {
                shells.push(model.add_shell(Shell {
                    faces: std::mem::take(&mut faces),
                }));
            }
        }
        if !full {
            let origin = profile.pose.apply_point(&RPoint3::origin().to_real());
            let reference = profile.pose.apply_vector(&frame.axis.to_real());
            faces.push(face(
                model,
                Surface::Plane(Plane {
                    origin,
                    normal: normal.clone(),
                    ref_dir: reference.clone(),
                }),
                reverse(forward),
                start_paths,
                None,
                vec![Origin {
                    step: step.clone(),
                    role: Role::Start,
                }],
            ));
            let origin = profile
                .pose
                .apply_point(&frame.rotate(&RPoint3::origin(), &cosine, &sine)?.into());
            let end_normal = frame
                .rotate(&frame.origin.plus(&unit_z()), &cosine, &sine)?
                .minus(&frame.origin);
            faces.push(face(
                model,
                Surface::Plane(Plane {
                    origin,
                    normal: profile.pose.apply_vector(&end_normal.to_real()),
                    ref_dir: reference,
                }),
                forward,
                end_paths,
                None,
                vec![Origin {
                    step: step.clone(),
                    role: Role::End,
                }],
            ));
        }
        if !full {
            shells.push(model.add_shell(Shell { faces }));
        }
        let lump = model.add_lump(Lump { shells });
        let body = model.add_body(Body {
            kind: BodyKind::Solid,
            lumps: vec![lump],
            marks: MarkLog::new(),
        });
        let body = if merge_caps {
            catallax_boolean::make_maximal(model, body, dec)?.value
        } else {
            body
        };
        bodies.push(body);
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
