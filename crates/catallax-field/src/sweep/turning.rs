//! Field transport of a profile by meridional sectors and straight translations.

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{Axis, BooleanKind, PathSegment, Role, StepId};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};

use crate::{Plan, features::Feature, profile::Profile};

use super::{circular_sectors, point};

pub(super) fn needs_turn(origin: &RPoint3, segments: &[PathSegment]) -> bool {
    let mut previous = origin.clone();
    let mut direction: Option<RVector3> = None;
    for segment in segments {
        let PathSegment::Line { to, .. } = segment else {
            return matches!(segment, PathSegment::Arc { .. });
        };
        let end = point(to);
        let next = end.minus(&previous);
        if direction
            .as_ref()
            .is_some_and(|d| !d.cross(&next).is_zero() || d.dot(&next).sign() == Sign::Negative)
        {
            return true;
        }
        direction = Some(next);
        previous = end;
    }
    false
}

fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    let [x, y, z] = v.0.each_ref().map(|n| n.value().clone());
    RVector3::new(x, y, z)
}

struct Turn {
    centre: RPoint3,
    axis: RVector3,
    angle: i64,
}

impl Turn {
    fn pose(&self) -> Result<Pose, Refusal> {
        let w = if self.angle == 180 {
            Rational::zero()
        } else {
            let length = crate::profile::sqrt(&self.axis.norm_squared()).ok_or_else(|| {
                Refusal::new(
                    "field.sweep-motion-unsupported",
                    "quarter turn requires a rational unit axis",
                )
            })?;
            if self.angle == 270 { -length } else { length }
        };
        let RVector3 { x, y, z } = self.axis.clone();
        let rotation = Rotation::from_quaternion(w, x, y, z)?;
        let centre = self.centre.minus(&RPoint3::origin());
        let rotated = rotation
            .apply(&centre.clone().into())
            .as_rational()
            .expect("rational turn");
        Ok(Pose::from_rational(rotation, centre.minus(&rotated)))
    }

    fn sweep(
        &self,
        plan: &mut Plan,
        profile: &Profile,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let normal = crate::profile::rational_vector(
            &profile.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
        // Rotation about the sketch normal stays in a single plane.
        if self.axis.cross(&normal).is_zero() {
            return Ok(0);
        }
        let axis = Axis {
            point: catallax_recipe::Point3(
                [
                    self.centre.x.clone(),
                    self.centre.y.clone(),
                    self.centre.z.clone(),
                ]
                .map(Into::into),
            ),
            direction: catallax_recipe::Vector3(
                [
                    self.axis.x.clone(),
                    self.axis.y.clone(),
                    self.axis.z.clone(),
                ]
                .map(Into::into),
            ),
        };
        let built = Feature::revolve(profile, step, &axis, &Rational::from_i64(self.angle), dec)?;
        if let Some(sectors) = built
            .as_ref()
            .and_then(|(_, faces)| circular_sectors(profile, faces, self.angle))
        {
            let mut result = 0;
            for (cell, faces) in sectors {
                let node = plan.leaf(crate::geometry::Leaf::Cell(cell), faces, dec)?;
                result = plan.boolean(BooleanKind::Union, result, node, dec)?;
            }
            Ok(result)
        } else {
            plan.motion_feature(built, dec)
        }
    }
}

impl Plan {
    pub(super) fn turning_sweep(
        &mut self,
        profile: &Profile,
        mut previous: RPoint3,
        segments: &[PathSegment],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut carried = profile.clone();
        let mut tangent: Option<RVector3> = None;
        let mut result = 0;
        let normal = RVector3::new(Rational::zero(), Rational::zero(), Rational::one());
        for (index, segment) in segments.iter().enumerate() {
            let (end, incoming, outgoing, arc) = match segment {
                PathSegment::Line { to, .. } => {
                    let end = point(to);
                    let delta = end.minus(&previous);
                    if delta.is_zero() {
                        return Err(Refusal::new(
                            "features.zero-path-segment",
                            format!("segment {index}: zero displacement"),
                        ));
                    }
                    (end, delta.clone(), delta, None)
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
                        || !axis.dot(&a).is_zero()
                        || !axis.dot(&b).is_zero()
                    {
                        return Err(Refusal::new(
                            "path.arc",
                            format!("invalid circular path segment {}", segment.id()),
                        ));
                    }
                    let angle = match (a.dot(&b).sign(), a.cross(&b).dot(&axis).sign()) {
                        (Sign::Zero, Sign::Positive) => 90,
                        (Sign::Negative, Sign::Zero) => 180,
                        (Sign::Zero, Sign::Negative) => 270,
                        _ => {
                            return Err(Refusal::new(
                                "field.sweep-motion-unsupported",
                                format!("arc {} is outside the cardinal angles", segment.id()),
                            ));
                        }
                    };
                    let incoming = axis.cross(&a);
                    let outgoing = axis.cross(&b);
                    (
                        end,
                        incoming,
                        outgoing,
                        Some(Turn {
                            centre,
                            axis,
                            angle,
                        }),
                    )
                }
                _ => {
                    return Err(Refusal::new(
                        "field.sweep-motion-unsupported",
                        format!(
                            "segment {} is outside line and circular transport",
                            segment.id()
                        ),
                    ));
                }
            };
            let saved = segment.corner_axis().map(vector);
            let reversing = tangent.as_ref().is_some_and(|t| {
                t.cross(&incoming).is_zero() && t.dot(&incoming).sign() == Sign::Negative
            });
            if saved
                .as_ref()
                .is_some_and(|axis| !reversing || axis.is_zero() || !axis.dot(&incoming).is_zero())
            {
                return Err(Refusal::new(
                    "path.corner-axis",
                    "corner_axis requires a reversing join and a nonzero perpendicular axis",
                ));
            }
            if let Some(tangent) = tangent {
                let axis = tangent.cross(&incoming);
                let turn = if reversing {
                    Some(Turn {
                        centre: previous.clone(),
                        axis: saved.ok_or_else(|| {
                            Refusal::new(
                                "sweep.missing-corner-axis",
                                "a reversing path needs corner_axis",
                            )
                        })?,
                        angle: 180,
                    })
                } else if axis.is_zero() {
                    None
                } else if tangent.dot(&incoming).is_zero() {
                    Some(Turn {
                        centre: previous.clone(),
                        axis,
                        angle: 90,
                    })
                } else {
                    return Err(Refusal::new(
                        "field.sweep-corner-unsupported",
                        format!(
                            "join {index} of {} is outside the cardinal turns",
                            segments.len()
                        ),
                    ));
                };
                if let Some(turn) = turn {
                    let pose = turn.pose()?;
                    let first = self.faces.len();
                    let node = turn.sweep(self, &carried, step, dec)?;
                    self.name_motion_faces(first, segment, Some(segment), false);
                    result = self.boolean(BooleanKind::Union, result, node, dec)?;
                    carried.pose = carried.pose.then(&pose);
                }
            }

            let first = self.faces.len();
            let node = if let Some(arc) = arc {
                let pose = arc.pose()?;
                let node = arc.sweep(self, &carried, step, dec)?;
                carried.pose = carried.pose.then(&pose);
                node
            } else {
                let carried_normal = crate::profile::rational_vector(&carried.pose, &normal);
                let node = if incoming.dot(&carried_normal).is_zero() {
                    0
                } else {
                    let built =
                        Feature::extrude_vector(&carried, step, incoming.clone(), false, dec)?;
                    self.motion_feature(built, dec)?
                };
                carried.pose = carried
                    .pose
                    .then(&Pose::from_rational(Rotation::identity(), incoming));
                node
            };
            self.name_motion_faces(first, segment, segments.get(index + 1), index == 0);
            result = self.boolean(BooleanKind::Union, result, node, dec)?;
            tangent = Some(outgoing);
            previous = end;
        }
        Ok(result)
    }

    fn name_motion_faces(
        &mut self,
        first: usize,
        segment: &PathSegment,
        next: Option<&PathSegment>,
        start: bool,
    ) {
        for i in first..self.faces.len() {
            let original = self.faces[i].origin.role.clone();
            let mapped = match &original {
                Role::Side {
                    segment: profile, ..
                } => Role::Side {
                    segment: profile.clone(),
                    path: Some(segment.id().clone()),
                },
                Role::Start if !start => Role::Envelope {
                    path: segment.id().clone(),
                },
                Role::End if next.is_some() => Role::Envelope {
                    path: next.expect("next segment").id().clone(),
                },
                _ => continue,
            };
            self.definition_roles.insert(i, original);
            self.faces[i].origin.role = mapped;
        }
    }
}
