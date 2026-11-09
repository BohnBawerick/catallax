//! Exact line/cardinal-arc motions, stationary turns and regularised analytic unions.

use std::sync::Arc;

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{
    Axis, ElementId, Origin, PathSegment, Point3, Role, StepId, SweepOrientation,
};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};
use catallax_topology::{BodyId, BodyKind, Model};

use crate::Profile;

#[derive(Clone, Debug)]
pub(crate) struct Path(Arc<Definition>);

#[derive(Debug)]
struct Definition {
    start: RPoint3,
    segments: Vec<Segment>,
    extent: Rational,
}

#[derive(Debug)]
struct Segment {
    id: ElementId,
    end: RPoint3,
    incoming: RVector3,
    outgoing: RVector3,
    motion: Motion,
    corner_axis: Option<RVector3>,
}

#[derive(Debug)]
enum Motion {
    Line,
    Arc(Box<Turn>),
}

#[derive(Debug)]
struct Turn {
    centre: RPoint3,
    axis: RVector3,
    angle: i64,
}

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("features.unsupported-path", detail)
}
fn point(p: &Point3) -> RPoint3 {
    RPoint3::new(
        p.0[0].value().clone(),
        p.0[1].value().clone(),
        p.0[2].value().clone(),
    )
}
fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    RVector3::new(
        v.0[0].value().clone(),
        v.0[1].value().clone(),
        v.0[2].value().clone(),
    )
}
fn length_bound(v: &RVector3) -> Rational {
    v.x.abs() + v.y.abs() + v.z.abs()
}

impl Turn {
    fn rotation(&self) -> Result<Rotation, Refusal> {
        let axis = if self.angle == 180 {
            self.axis.clone()
        } else {
            let length = crate::extrude::sqrt(&self.axis.norm_squared()).ok_or_else(|| {
                unsupported("quarter-turn transport needs a rational unit rotation axis")
            })?;
            self.axis.scaled(&(Rational::one() / length))
        };
        let w = match self.angle {
            90 => 1,
            180 => 0,
            270 => -1,
            _ => unreachable!("checked cardinal turn"),
        };
        Rotation::from_quaternion(Rational::from_i64(w), axis.x, axis.y, axis.z)
    }

    fn pose(&self) -> Result<Pose, Refusal> {
        let rotation = self.rotation()?;
        let centre = self.centre.minus(&RPoint3::origin());
        let rotated = rotation
            .apply(&centre.to_real())
            .as_rational()
            .ok_or_else(|| {
                unsupported("analytic sweep continuation needs a rational transported frame")
            })?;
        Ok(Pose::from_rational(
            rotation,
            RVector3::new(
                centre.x - rotated.x,
                centre.y - rotated.y,
                centre.z - rotated.z,
            ),
        ))
    }

    fn build(
        &self,
        model: &mut Model,
        profile: &Profile,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        let normal = profile.pose.apply_vector(
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
        );
        if crate::profile::sign(&normal.cross(&self.axis.to_real()).norm_squared(), dec)?
            == Sign::Zero
        {
            // Rotation about the profile normal stays in one plane.
            return Ok(catallax_topology::empty(model));
        }
        let axis = Axis {
            point: Point3([
                self.centre.x.clone().into(),
                self.centre.y.clone().into(),
                self.centre.z.clone().into(),
            ]),
            direction: catallax_recipe::Vector3([
                self.axis.x.clone().into(),
                self.axis.y.clone().into(),
                self.axis.z.clone().into(),
            ]),
        };
        crate::revolve::revolve(
            model,
            profile,
            step,
            &axis,
            &Rational::from_i64(self.angle),
            dec,
        )
        .map_err(|e| match e.code {
            "revolve.axis-plane" | "revolve.axis-interior" => {
                unsupported("this rotating profile needs folded or non-meridional sweep carriers")
            }
            _ => e,
        })
    }
}

impl Path {
    pub(crate) fn new(start: &Point3, segments: &[PathSegment]) -> Result<Self, Refusal> {
        let start = point(start);
        let mut previous = start.clone();
        let mut checked: Vec<Segment> = Vec::with_capacity(segments.len());
        let mut extent = Rational::zero();
        let mut rotates = false;
        for segment in segments {
            let (end, incoming, outgoing, motion) = match segment {
                PathSegment::Line { to, .. } => {
                    let end = point(to);
                    let delta = end.minus(&previous);
                    if delta.is_zero() {
                        return Err(Refusal::new(
                            "path.zero-length",
                            "a line path segment must have distinct endpoints",
                        ));
                    }
                    extent = extent + length_bound(&delta);
                    (end, delta.clone(), delta, Motion::Line)
                }
                PathSegment::Arc {
                    centre, to, axis, ..
                } => {
                    let centre = point(centre);
                    let end = point(to);
                    let axis = vector(axis);
                    let a = previous.minus(&centre);
                    let b = end.minus(&centre);
                    if axis.is_zero()
                        || a.is_zero()
                        || a == b
                        || a.norm_squared() != b.norm_squared()
                        || !a.dot(&axis).is_zero()
                        || !b.dot(&axis).is_zero()
                    {
                        return Err(Refusal::new(
                            "path.arc",
                            "an arc needs distinct endpoints, a nonzero axis and equal nonzero perpendicular radii",
                        ));
                    }
                    let dot = a.dot(&b);
                    let angle = if a.cross(&b).is_zero() {
                        180
                    } else if dot.is_zero() {
                        if axis.dot(&a.cross(&b)).sign() == Sign::Positive {
                            90
                        } else {
                            270
                        }
                    } else {
                        return Err(unsupported(
                            "analytic arc sweeps currently need quarter, half or three-quarter circular turns",
                        ));
                    };
                    extent = extent + Rational::from_i64(2) * length_bound(&a);
                    rotates = true;
                    let incoming = axis.cross(&a);
                    let outgoing = axis.cross(&b);
                    (
                        end,
                        incoming,
                        outgoing,
                        Motion::Arc(Box::new(Turn {
                            centre,
                            axis,
                            angle,
                        })),
                    )
                }
                PathSegment::Spline { .. } => {
                    return Err(unsupported(
                        "spline path definitions are retained, but checked sweep carriers are not implemented",
                    ));
                }
            };
            let corner_axis = segment.corner_axis().map(vector);
            if let Some(axis) = &corner_axis {
                let reversing = checked.last().is_some_and(|s| {
                    s.outgoing.cross(&incoming).is_zero()
                        && s.outgoing.dot(&incoming).sign() == Sign::Negative
                });
                if !reversing || axis.is_zero() || !axis.dot(&incoming).is_zero() {
                    return Err(Refusal::new(
                        "path.corner-axis",
                        "corner_axis requires a reversing join and a nonzero perpendicular axis",
                    ));
                }
            }
            checked.push(Segment {
                id: segment.id().clone(),
                end: end.clone(),
                incoming,
                outgoing,
                motion,
                corner_axis,
            });
            previous = end;
        }
        // |gamma(s)-Q gamma(0)| <= path displacement + 2 |gamma(0)|.
        // Include this bound for turns at line joins as well as circular spans.
        if rotates || checked.len() > 1 {
            extent =
                extent + Rational::from_i64(2) * length_bound(&start.minus(&RPoint3::origin()));
        }
        Ok(Self(Arc::new(Definition {
            start,
            segments: checked,
            extent,
        })))
    }

    pub(crate) fn extent(&self) -> Rational {
        self.0.extent.clone()
    }

    pub(crate) fn build(
        &self,
        model: &mut Model,
        profile: &Profile,
        step: &StepId,
        orientation: SweepOrientation,
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        profile.require_circular("sweep.ellipse-profile")?;
        if profile.is_empty() || self.0.segments.is_empty() {
            return Ok(catallax_topology::empty(model));
        }
        if let [
            Segment {
                id,
                incoming: delta,
                motion: Motion::Line,
                ..
            },
        ] = self.0.segments.as_slice()
        {
            return line(model, profile, step, id, delta, dec);
        }
        let mut carried = profile.clone();
        let mut previous = self.0.start.clone();
        let mut tangent: Option<&RVector3> = None;
        let mut result = catallax_topology::empty(model);
        for (i, segment) in self.0.segments.iter().enumerate() {
            if orientation == SweepOrientation::RotationMinimising
                && let Some(incoming) = tangent
            {
                let cross = incoming.cross(&segment.incoming);
                let dot = incoming.dot(&segment.incoming);
                let turn = if cross.is_zero() {
                    if dot.sign() == Sign::Negative {
                        let axis = segment.corner_axis.as_ref().ok_or_else(|| Refusal::new(
                                "sweep.missing-corner-axis", "a reversing rotation-minimising path needs the outgoing corner_axis",
                            ))?;
                        Some(Turn {
                            centre: previous.clone(),
                            axis: axis.clone(),
                            angle: 180,
                        })
                    } else {
                        None
                    }
                } else if dot.is_zero() {
                    Some(Turn {
                        centre: previous.clone(),
                        axis: cross,
                        angle: 90,
                    })
                } else {
                    return Err(unsupported(
                        "stationary turns currently need agreeing, perpendicular or reversing tangents",
                    ));
                };
                if let Some(turn) = turn {
                    let body = turn.build(model, &carried, step, dec)?;
                    remap(model, body, step, &segment.id, false, false, &segment.id);
                    result = unite(model, result, body, dec)?;
                    carried.pose = carried.pose.then(&turn.pose()?);
                }
            }
            let body = match &segment.motion {
                Motion::Line => {
                    let delta = &segment.incoming;
                    let body = line(model, &carried, step, &segment.id, delta, dec)?;
                    carried.pose = carried
                        .pose
                        .then(&Pose::from_rational(Rotation::identity(), delta.clone()));
                    body
                }
                Motion::Arc(turn) => {
                    if orientation == SweepOrientation::Fixed {
                        return Err(unsupported(
                            "fixed-orientation circular motion needs its translated-profile envelope carriers",
                        ));
                    }
                    let body = turn.build(model, &carried, step, dec)?;
                    carried.pose = carried.pose.then(&turn.pose()?);
                    body
                }
            };
            let next = self.0.segments.get(i + 1).map_or(&segment.id, |s| &s.id);
            remap(
                model,
                body,
                step,
                &segment.id,
                i == 0,
                i + 1 == self.0.segments.len(),
                next,
            );
            result = unite(model, result, body, dec)?;
            tangent = Some(&segment.outgoing);
            previous = segment.end.clone();
        }
        Ok(result)
    }
}

fn remap(
    model: &mut Model,
    body: BodyId,
    step: &StepId,
    path: &ElementId,
    start: bool,
    end: bool,
    next: &ElementId,
) {
    model.map_face_origins(body, |origin| Origin {
        step: step.clone(),
        role: match &origin.role {
            Role::Side { segment, .. } => Role::Side {
                segment: segment.clone(),
                path: Some(path.clone()),
            },
            Role::Start if !start => Role::Envelope { path: path.clone() },
            Role::End if !end => Role::Envelope { path: next.clone() },
            role => role.clone(),
        },
    });
}

fn unite(model: &mut Model, a: BodyId, b: BodyId, dec: &mut Decider) -> Result<BodyId, Refusal> {
    if model.body(a).kind == BodyKind::Empty {
        return Ok(b);
    }
    if model.body(b).kind == BodyKind::Empty {
        return Ok(a);
    }
    Ok(
        catallax_boolean::boolean(model, a, b, catallax_boolean::BoolOp::Unite, dec)
            .map_err(|e| e.refusal)?
            .body,
    )
}

fn line(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    id: &ElementId,
    delta: &RVector3,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    if delta.dot(&normal).is_zero() {
        return Ok(catallax_topology::empty(model));
    }
    crate::extrude::extrude_vector(
        model,
        profile,
        step,
        delta,
        [Some(Role::Start), Some(Role::End)],
        Some(id),
        dec,
    )
}
