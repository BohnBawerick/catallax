//! Exact translational sweep and loft reductions over field profile definitions.

use crate::{Plan, features::Feature, profile::Profile};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{
    BooleanKind, LoftSection, LoftShape, Operation, Recipe, Role, StepId, SweepOrientation,
};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};
use std::collections::BTreeMap;

fn point(p: &catallax_recipe::Point3) -> RPoint3 {
    let [x, y, z] = p.0.each_ref().map(|n| n.value().clone());
    RPoint3::new(x, y, z)
}

fn authored_loops<'a>(recipe: &'a Recipe, id: &StepId) -> &'a [catallax_recipe::Loop] {
    let Operation::Sketch { loops, .. } = &recipe.step(id).expect("checked section").op else {
        unreachable!("checked sketch reference")
    };
    loops
}

fn same_segments(a: &catallax_recipe::Segment, b: &catallax_recipe::Segment) -> bool {
    use catallax_recipe::Segment;
    let mut b = b.clone();
    let id = match &mut b {
        Segment::Line { id, .. }
        | Segment::Arc { id, .. }
        | Segment::Circle { id, .. }
        | Segment::Bezier { id, .. }
        | Segment::Ellipse { id, .. }
        | Segment::EllipticArc { id, .. } => id,
    };
    *id = a.id().clone();
    a == &b
}

// A full circular profile rotated through a cardinal angle is a ring torus
// clipped by exact meridional halfspaces. For 180 degrees, split at the
// perpendicular meridian so Start and End keep their separate cap domains.
fn circular_sectors(
    profile: &Profile,
    faces: &[crate::Face],
    angle: i64,
) -> Option<Vec<(crate::geometry::Cell, Vec<crate::Face>)>> {
    use catallax_shapes::{Plane, Surface};
    let [edges] = profile.loops.as_slice() else {
        return None;
    };
    let [edge] = edges.as_slice() else {
        return None;
    };
    let (centre, radius_sq, _) = edge.circle.as_ref()?;
    if !edge.full || faces.len() != 3 {
        return None;
    }
    let side = faces
        .iter()
        .find(|f| matches!(f.origin.role, Role::Side { .. }))?;
    let start = faces.iter().find(|f| f.origin.role == Role::Start)?;
    let end = faces.iter().find(|f| f.origin.role == Role::End)?;
    let Surface::Torus(torus) = &side.surface else {
        return None;
    };
    let center = torus.center.as_rational()?;
    let axis = torus.axis.as_rational()?;
    let major = torus.major.as_rational()?;
    let minor = torus.minor.as_rational()?;
    let radial = profile.world(centre).minus(&center);
    if minor.sign() != Sign::Positive
        || major <= minor
        || &minor * &minor != *radius_sq
        || radial.norm_squared() != &major * &major
        || !axis.dot(&radial).is_zero()
        || axis.is_zero()
    {
        return None;
    }
    let tangent = axis.cross(&radial);
    let witness = match angle {
        90 => radial.plus(&tangent),
        180 => tangent.clone(),
        270 => tangent.minus(&radial),
        _ => return None,
    };
    let orient = |face: &crate::Face| -> Option<Plane> {
        let Surface::Plane(plane) = &face.surface else {
            return None;
        };
        let normal = plane.normal.as_rational()?;
        let origin = plane.origin.as_rational()?;
        if !center.minus(&origin).dot(&normal).is_zero() || !axis.dot(&normal).is_zero() {
            return None;
        }
        let sign = normal.dot(&witness).sign();
        if sign == Sign::Zero {
            return None;
        }
        let mut plane = plane.clone();
        if sign == Sign::Positive {
            plane.normal = normal.scaled(&Rational::from_i64(-1)).into();
        }
        Some(plane)
    };
    let start_plane = orient(start)?;
    let end_plane = orient(end)?;
    if !start_plane.normal.as_rational()?.dot(&radial).is_zero()
        || !end_plane
            .normal
            .as_rational()?
            .dot(if angle == 180 { &radial } else { &tangent })
            .is_zero()
    {
        return None;
    }
    let radius = center.x.abs() + center.y.abs() + center.z.abs() + major + minor;
    let make = |planes: Vec<Plane>, faces: Vec<crate::Face>| {
        let mut surfaces = vec![side.surface.clone()];
        surfaces.extend(planes.into_iter().map(Surface::Plane));
        crate::geometry::Cell::new(surfaces, Some(radius.clone()))
            .ok()
            .map(|cell| (cell, faces))
    };
    match angle {
        90 => Some(vec![make(vec![start_plane, end_plane], faces.to_vec())?]),
        180 => {
            let split = |sign| Plane {
                origin: center.clone().into(),
                normal: radial.scaled(&Rational::from_i64(sign)).into(),
                ref_dir: axis.clone().into(),
            };
            Some(vec![
                make(
                    vec![start_plane, split(-1)],
                    vec![side.clone(), start.clone()],
                )?,
                make(vec![end_plane, split(1)], vec![side.clone(), end.clone()])?,
            ])
        }
        270 => Some(vec![
            make(vec![start_plane], vec![side.clone(), start.clone()])?,
            make(vec![end_plane], vec![side.clone(), end.clone()])?,
        ]),
        _ => None,
    }
}

impl Plan {
    fn motion_feature(
        &mut self,
        built: Option<(Feature, Vec<crate::Face>)>,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let Some((feature, faces)) = built else {
            return Ok(0);
        };
        if let Some(planes) = feature.convex_planes(dec)? {
            let radius = feature.display_radius(dec)?;
            let cell = crate::geometry::Cell::new(
                planes
                    .into_iter()
                    .map(catallax_shapes::Surface::Plane)
                    .collect(),
                radius,
            )?;
            return self.leaf(crate::geometry::Leaf::Cell(cell), faces, dec);
        }
        self.feature(Some((feature, faces)), dec)
    }

    pub(crate) fn sweep(
        &mut self,
        profile: &Profile,
        path: &Operation,
        orientation: SweepOrientation,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let Operation::Path { start, segments } = path else {
            unreachable!("checked path reference")
        };
        let origin = point(start);
        if let [
            catallax_recipe::PathSegment::Arc {
                id,
                centre,
                to,
                axis,
                ..
            },
        ] = segments.as_slice()
            && orientation == SweepOrientation::RotationMinimising
        {
            let centre_point = point(centre);
            let a = origin.minus(&centre_point);
            let b = point(to).minus(&centre_point);
            let [x, y, z] = axis.0.each_ref().map(|n| n.value().clone());
            let axis_vector = RVector3::new(x, y, z);
            if axis_vector.is_zero()
                || !a.dot(&axis_vector).is_zero()
                || !b.dot(&axis_vector).is_zero()
                || a.is_zero()
                || a.norm_squared() != b.norm_squared()
            {
                return Err(Refusal::new(
                    "features.invalid-path-arc",
                    "arc radii and axis do not define a circular span",
                ));
            }
            let dot = a.dot(&b);
            let cross = a.cross(&b).dot(&axis_vector);
            let angle = if dot.is_zero() && cross.sign() == Sign::Positive {
                90
            } else if dot.is_zero() && cross.sign() == Sign::Negative {
                270
            } else if cross.is_zero() && dot.sign() == Sign::Negative {
                180
            } else {
                return Err(Refusal::new(
                    "field.sweep-motion-unsupported",
                    "1 arc span: angle is outside 90, 180 and 270 degrees",
                ));
            };
            let first = self.faces.len();
            let built = Feature::revolve(
                profile,
                step,
                &catallax_recipe::Axis {
                    point: centre.clone(),
                    direction: axis.clone(),
                },
                &Rational::from_i64(angle),
                dec,
            )?;
            let sectors = built
                .as_ref()
                .and_then(|(_, faces)| circular_sectors(profile, faces, angle));
            let node = if let Some(sectors) = sectors {
                let mut result = 0;
                for (cell, faces) in sectors {
                    let node = self.leaf(crate::geometry::Leaf::Cell(cell), faces, dec)?;
                    result = self.boolean(BooleanKind::Union, result, node, dec)?;
                }
                result
            } else {
                self.motion_feature(built, dec)?
            };
            for i in first..self.faces.len() {
                if let Role::Side { segment, .. } = &self.faces[i].origin.role {
                    self.faces[i].origin.role = Role::Side {
                        segment: segment.clone(),
                        path: Some(id.clone()),
                    };
                }
            }
            return Ok(node);
        }
        let mut previous = origin.clone();
        let mut incoming: Option<RVector3> = None;
        let mut result = 0;
        for (index, segment) in segments.iter().enumerate() {
            let catallax_recipe::PathSegment::Line { id, to, .. } = segment else {
                return Err(Refusal::new(
                    "field.sweep-motion-unsupported",
                    format!(
                        "segment {index} of {}: {} needs a certified non-translational field; 0 such spans admitted",
                        segments.len(),
                        segment.id()
                    ),
                ));
            };
            let end = point(to);
            let vector = end.minus(&previous);
            if vector.is_zero() {
                return Err(Refusal::new(
                    "features.zero-path-segment",
                    format!("segment {index}: zero displacement"),
                ));
            }
            if orientation == SweepOrientation::RotationMinimising
                && incoming.as_ref().is_some_and(|v| {
                    !v.cross(&vector).is_zero() || v.dot(&vector).sign() != Sign::Positive
                })
            {
                return Err(Refusal::new(
                    "field.sweep-corner-unsupported",
                    format!(
                        "join {index} of {} needs a stationary rotation; 0 corner rotations admitted",
                        segments.len()
                    ),
                ));
            }
            let mut span = profile.clone();
            // M(0) is identity: the path start never relocates the sketch.
            span.pose = profile.pose.then(&Pose::from_rational(
                Rotation::identity(),
                previous.minus(&origin),
            ));
            let first = self.faces.len();
            let node = self.motion_feature(
                Feature::extrude_vector(&span, step, vector.clone(), false, dec)?,
                dec,
            )?;
            for i in first..self.faces.len() {
                let original = self.faces[i].origin.role.clone();
                self.faces[i].origin.role = match &original {
                    Role::Side { segment, .. } => Role::Side {
                        segment: segment.clone(),
                        path: Some(id.clone()),
                    },
                    Role::Start if index != 0 => Role::Envelope { path: id.clone() },
                    Role::End if index + 1 != segments.len() => Role::Envelope { path: id.clone() },
                    _ => original.clone(),
                };
                self.definition_roles.insert(i, original);
            }
            result = self.boolean(BooleanKind::Union, result, node, dec)?;
            previous = end;
            incoming = Some(vector);
        }
        Ok(result)
    }

    pub(crate) fn loft(
        &mut self,
        recipe: &Recipe,
        sections: &[LoftSection],
        shape: LoftShape,
        profiles: &BTreeMap<StepId, Profile>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if recipe.loft_parameterization().ok()
            != Some(catallax_recipe::LoftParameterization::RationalCharts)
        {
            return Err(Refusal::new(
                "loft.legacy-parameterization",
                "legacy arc-length lofts are not rational-chart lofts",
            ));
        }
        let first = &profiles[&sections[0].profile];
        let last = &profiles[&sections.last().expect("checked sections").profile];
        if let Some(node) = self.circular_loft(recipe, sections, profiles, step, dec)? {
            return Ok(node);
        }
        let vector = last
            .pose
            .translation
            .minus(&first.pose.translation)
            .as_rational()
            .ok_or_else(|| {
                Refusal::new(
                    "field.loft-translation-unsupported",
                    "end displacement is not rational",
                )
            })?;
        let fail = |index: usize| {
            Refusal::new(
                "field.loft-correspondence-unsupported",
                format!(
                    "section {index} of {} is not a certified straight translated correspondence; 0 general loft patches admitted",
                    sections.len()
                ),
            )
        };
        let normal = crate::profile::rational_vector(
            &first.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
        let height = vector.dot(&normal);
        if height.is_zero() {
            return Err(fail(sections.len() - 1));
        }
        let mut previous = Rational::zero();
        for (index, section) in sections.iter().enumerate() {
            let p = &profiles[&section.profile];
            let original = authored_loops(recipe, &sections[0].profile);
            let current = authored_loops(recipe, &section.profile);
            if original.len() != current.len()
                || original.iter().zip(current).any(|(a, b)| {
                    a.0.len() != b.0.len()
                        || a.0.iter().zip(&b.0).any(|(a, b)| !same_segments(a, b))
                })
            {
                return Err(fail(index));
            }
            if !p.pose.rotation.same_rotation(&first.pose.rotation, dec)?
                || p.loops.len() != first.loops.len()
            {
                return Err(fail(index));
            }
            for (a, b) in first.loops.iter().zip(&p.loops) {
                if a.len() != b.len()
                    || a.iter().zip(b).any(|(a, b)| {
                        a.start != b.start
                            || a.end != b.end
                            || a.circle != b.circle
                            || a.full != b.full
                    })
                {
                    return Err(fail(index));
                }
            }
            for (a, b) in sections[0].start.iter().zip(&section.start) {
                let locate = |loops: &[catallax_recipe::Loop], id: &catallax_recipe::ElementId| {
                    loops.iter().enumerate().find_map(|(l, edges)| {
                        edges.0.iter().position(|e| e.id() == id).map(|e| (l, e))
                    })
                };
                if locate(original, a) != locate(current, b) {
                    return Err(fail(index));
                }
            }
            let shift = p
                .pose
                .translation
                .minus(&first.pose.translation)
                .as_rational()
                .ok_or_else(|| fail(index))?;
            let t = shift.dot(&normal) / height.clone();
            if !shift.cross(&vector).is_zero() || t < previous || t > Rational::one() {
                return Err(fail(index));
            }
            if shape == LoftShape::Smooth
                && t != Rational::from_ratio(index as i64, (sections.len() - 1) as i64)
                    .expect("checked sections")
            {
                return Err(fail(index));
            }
            previous = t;
        }
        self.motion_feature(
            Feature::extrude_vector(first, step, vector, false, dec)?,
            dec,
        )
    }

    fn circular_loft(
        &mut self,
        recipe: &Recipe,
        sections: &[LoftSection],
        profiles: &BTreeMap<StepId, Profile>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        if sections.len() != 2 {
            return Ok(None);
        }
        let first = &profiles[&sections[0].profile];
        let last = &profiles[&sections[1].profile];
        let ([a], [b]) = (first.loops.as_slice(), last.loops.as_slice()) else {
            return Ok(None);
        };
        let ([a], [b]) = (a.as_slice(), b.as_slice()) else {
            return Ok(None);
        };
        let (Some((ca, ra, _)), Some((cb, rb, _))) = (&a.circle, &b.circle) else {
            return Ok(None);
        };
        if !a.full
            || !b.full
            || ca != cb
            || !first
                .pose
                .rotation
                .same_rotation(&last.pose.rotation, dec)?
        {
            return Ok(None);
        }
        // Circle charts start at the same local +x point in both sections.
        if sections[0].start[0] != a.ids[0] || sections[1].start[0] != b.ids[0] {
            return Ok(None);
        }
        let normal = crate::profile::rational_vector(
            &first.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
        let shift = last
            .pose
            .translation
            .minus(&first.pose.translation)
            .as_rational()
            .ok_or_else(|| {
                Refusal::new(
                    "field.loft-translation-unsupported",
                    "circular loft displacement is not rational",
                )
            })?;
        let height = shift.dot(&normal);
        if height.sign() != Sign::Positive || !shift.cross(&normal).is_zero() {
            return Ok(None);
        }
        let Some(ra) = crate::profile::sqrt(ra) else {
            return Ok(None);
        };
        let Some(rb) = crate::profile::sqrt(rb) else {
            return Ok(None);
        };
        let Operation::Sketch { plane, .. } = &recipe
            .step(&sections[0].profile)
            .expect("checked section")
            .op
        else {
            unreachable!()
        };
        let mut pose = plane.clone();
        let centre = crate::profile::rational_point(
            &first.pose,
            &RPoint3::new(ca[0].clone(), ca[1].clone(), Rational::zero()),
        );
        pose.translation = catallax_recipe::Vector3([centre.x, centre.y, centre.z].map(Into::into));
        let op = Operation::Cone {
            base_radius: ra.into(),
            top_radius: rb.into(),
            height: height.into(),
            pose,
        };
        let (leaf, mut faces) = crate::geometry::compile(&op, step, dec)?;
        let offset = self.faces.len();
        for (i, f) in faces.iter_mut().enumerate() {
            if matches!(f.origin.role, Role::Side { .. }) {
                self.definition_roles
                    .insert(offset + i, f.origin.role.clone());
                f.origin.role = Role::Side {
                    segment: a.ids[0].clone(),
                    path: None,
                };
            }
        }
        self.leaf(leaf, faces, dec).map(Some)
    }
}
