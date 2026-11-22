//! Exact sweep and loft constructions over field profile definitions.

mod affine_loft;
mod bezier_loft;
mod spline;
mod stations;
mod turning;

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
        if let [segment @ catallax_recipe::PathSegment::Spline { .. }] = segments.as_slice()
            && orientation == SweepOrientation::Fixed
        {
            return self.spline_sweep(profile, segment, dec);
        }
        if orientation == SweepOrientation::RotationMinimising
            && segments.len() > 1
            && turning::needs_turn(&origin, segments)
        {
            return self.turning_sweep(profile, origin, segments, step, dec);
        }
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
            let normal = crate::profile::rational_vector(
                &profile.pose,
                &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
            );
            if axis_vector.cross(&normal).is_zero() {
                return Ok(0);
            }
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
        let normal = crate::profile::rational_vector(
            &profile.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
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
            if vector.dot(&normal).is_zero() {
                previous = end;
                incoming = Some(vector);
                continue;
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
            if sections.iter().all(|section| {
                let pose = &profiles[&section.profile].pose;
                let other_normal = crate::profile::rational_vector(
                    pose,
                    &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
                );
                other_normal.cross(&normal).is_zero()
                    && pose
                        .translation
                        .minus(&first.pose.translation)
                        .as_rational()
                        .is_some_and(|shift| shift.dot(&normal).is_zero())
            }) {
                return Ok(0);
            }
            return self.patch_or(
                recipe,
                fail(sections.len() - 1),
                sections,
                shape,
                profiles,
                step,
                dec,
            );
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
                return self.patch_or(recipe, fail(index), sections, shape, profiles, step, dec);
            }
            if !p.pose.rotation.same_rotation(&first.pose.rotation, dec)?
                || p.loops.len() != first.loops.len()
            {
                return self.patch_or(recipe, fail(index), sections, shape, profiles, step, dec);
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
                    return self.patch_or(recipe, fail(index), sections, shape, profiles, step, dec);
                }
            }
            for (a, b) in sections[0].start.iter().zip(&section.start) {
                let locate = |loops: &[catallax_recipe::Loop], id: &catallax_recipe::ElementId| {
                    loops.iter().enumerate().find_map(|(l, edges)| {
                        edges.0.iter().position(|e| e.id() == id).map(|e| (l, e))
                    })
                };
                if locate(original, a) != locate(current, b) {
                    return self.patch_or(recipe, fail(index), sections, shape, profiles, step, dec);
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
                return self.patch_or(recipe, fail(index), sections, shape, profiles, step, dec);
            }
            if shape == LoftShape::Smooth
                && t != Rational::from_ratio(index as i64, (sections.len() - 1) as i64)
                    .expect("checked sections")
            {
                return self.patch_or(recipe, fail(index), sections, shape, profiles, step, dec);
            }
            previous = t;
        }
        self.motion_feature(
            Feature::extrude_vector(first, step, vector, false, dec)?,
            dec,
        )
    }

    /// Tries the patch-shell lofts before returning the correspondence refusal.
    fn patch_or(
        &mut self,
        recipe: &Recipe,
        refusal: Refusal,
        sections: &[LoftSection],
        shape: LoftShape,
        profiles: &BTreeMap<StepId, Profile>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if let Some(node) = self.circular_station_loft(recipe, sections, shape, profiles, dec)? {
            return Ok(node);
        }
        let node = self.patch_loft(sections, profiles, step, dec)?;
        if let Some(node) = node {
            return Ok(node);
        }
        if let Some(node) = self.affine_loft(sections, profiles, dec)? {
            return Ok(node);
        }
        self.polygon_station_loft(recipe, sections, shape, profiles, dec)?
            .ok_or(refusal)
    }

    /// A loft between two parallel polygons, every slice of which is a
    /// convex polygon, held as a patch shell (4.2a step 7.2): a bilinear side
    /// for each pair of corresponding segments, and the two sections as
    /// polygon caps. With two stations natural interpolation is linear, so
    /// ruled and smooth use the same sides. `None` for any other loft.
    ///
    /// The sections share a rotation, so the slice at `v` in `[0, 1]` of the
    /// height is the polygon whose vertices are `(1 - v) a_i + v b_i`, with
    /// `b_i` moved by the in-plane offset of the second section. The turn at
    /// a vertex, the cross product of its two edges, is then a quadratic in
    /// `v` with Bernstein coefficients `c0 = a' x a''`, `c1 = (a' x b'' + b' x
    /// a'') / 2` and `c2 = b' x b''` for the edges `'` before and `''` after
    /// it. When every coefficient of every vertex has one strict sign, every
    /// slice turns strictly one way at every vertex, so it is a strictly
    /// convex polygon (it starts as the simple first section and no edge
    /// vanishes on the way). Then the sides meet only along the lines through
    /// corresponding vertices and the caps only at the sections' edges: the
    /// faces form a closed shell with no crossings, which is what the leaf
    /// needs. The sides are exactly the B-rep's ruled sides (the M11 loft's
    /// bilinear patches): straight lines between points at equal parameters
    /// of corresponding line segments.
    fn patch_loft(
        &mut self,
        sections: &[LoftSection],
        profiles: &BTreeMap<StepId, Profile>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        use crate::patch_shell::{Face, PatchFace, Polygon, Shell};
        use catallax_shapes::{BSplineSurface, Point3, RationalPatch};
        if sections.len() != 2 {
            return Ok(None);
        }
        let first = &profiles[&sections[0].profile];
        let last = &profiles[&sections[1].profile];
        let ([a], [b]) = (first.loops.as_slice(), last.loops.as_slice()) else {
            return Ok(None);
        };
        let n = a.len();
        if n < 3
            || b.len() != n
            || a.iter()
                .chain(b)
                .any(|e| e.circle.is_some() || e.full || e.ids.len() != 1)
            || !last
                .pose
                .rotation
                .same_rotation(&first.pose.rotation, dec)?
        {
            return Ok(None);
        }
        let start = |edges: &[crate::profile::Edge], section: &LoftSection| {
            let id = section.start.first()?;
            edges.iter().position(|e| &e.ids[0] == id)
        };
        let (Some(sa), Some(sb)) = (start(a, &sections[0]), start(b, &sections[1])) else {
            return Ok(None);
        };
        let Some(shift) = last
            .pose
            .translation
            .minus(&first.pose.translation)
            .as_rational()
        else {
            return Ok(None);
        };
        let axis = |k: usize| {
            let mut v = [Rational::zero(), Rational::zero(), Rational::zero()];
            v[k] = Rational::one();
            crate::profile::rational_vector(
                &first.pose,
                &RVector3::new(v[0].clone(), v[1].clone(), v[2].clone()),
            )
        };
        let offset = [shift.dot(&axis(0)), shift.dot(&axis(1))];
        let height = shift.dot(&axis(2));
        if height.is_zero() {
            return Ok(None);
        }
        let lower: Vec<[Rational; 2]> = (0..n).map(|k| a[(sa + k) % n].start.clone()).collect();
        let upper: Vec<[Rational; 2]> = (0..n)
            .map(|k| {
                let p = &b[(sb + k) % n].start;
                [&p[0] + &offset[0], &p[1] + &offset[1]]
            })
            .collect();
        let edge = |p: &[[Rational; 2]], i: usize| {
            let (s, t) = (&p[i % n], &p[(i + 1) % n]);
            [&t[0] - &s[0], &t[1] - &s[1]]
        };
        let cross = |u: &[Rational; 2], v: &[Rational; 2]| &(&u[0] * &v[1]) - &(&u[1] * &v[0]);
        let mut sense = None;
        for i in 0..n {
            let (before_a, after_a) = (edge(&lower, i + n - 1), edge(&lower, i));
            let (before_b, after_b) = (edge(&upper, i + n - 1), edge(&upper, i));
            let middle = (&cross(&before_a, &after_b) + &cross(&before_b, &after_a)).mul_pow2(-1);
            for c in [
                cross(&before_a, &after_a),
                middle,
                cross(&before_b, &after_b),
            ] {
                let sign = c.sign();
                if sign == Sign::Zero || *sense.get_or_insert(sign) != sign {
                    return Ok(None);
                }
            }
        }
        let counterclockwise = sense == Some(Sign::Positive);
        let up = height > Rational::zero();
        let world = |pose: &Pose, p: &[Rational; 2]| {
            crate::profile::rational_point(
                pose,
                &RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero()),
            )
        };
        let bottom: Vec<RPoint3> = (0..n)
            .map(|k| world(&first.pose, &a[(sa + k) % n].start))
            .collect();
        let top: Vec<RPoint3> = (0..n)
            .map(|k| world(&last.pose, &b[(sb + k) % n].start))
            .collect();
        let control =
            |p: &RPoint3| Point3::new(p.x.clone().into(), p.y.clone().into(), p.z.clone().into());
        let knots = vec![
            Rational::zero(),
            Rational::zero(),
            Rational::one(),
            Rational::one(),
        ];
        let mut faces = Vec::with_capacity(n + 2);
        for k in 0..n {
            let next = (k + 1) % n;
            let carrier = BSplineSurface::new(
                [1, 1],
                [knots.clone(), knots.clone()],
                [2, 2],
                vec![
                    control(&bottom[k]),
                    control(&top[k]),
                    control(&bottom[next]),
                    control(&top[next]),
                ],
                vec![catallax_numbers::Real::one(); 4],
                dec,
            )?;
            let patch = match RationalPatch::new(carrier, dec) {
                Ok(patch) => patch,
                Err(refusal) if refusal.code == "shapes.rational-patch-unproven" => {
                    return Ok(None);
                }
                Err(refusal) => return Err(refusal),
            };
            // `u` runs along the segments in the loop's sense and `v` toward
            // the second section, so the natural normal points out of a
            // counterclockwise loft that rises along the frame's normal.
            faces.push(Face::Patch(Box::new(PatchFace::new(
                &patch,
                counterclockwise == up,
                dec,
            )?)));
        }
        // Each cap's vertices run counterclockwise seen from outside: the
        // first section is seen from below the loft, the second from above.
        let mut bottom = bottom;
        let mut top = top;
        if counterclockwise == up {
            bottom.reverse();
        } else {
            top.reverse();
        }
        faces.push(Face::Polygon(Box::new(Polygon::new(bottom)?)));
        faces.push(Face::Polygon(Box::new(Polygon::new(top)?)));
        let sides = (0..n)
            .map(|k| catallax_recipe::Origin {
                step: step.clone(),
                role: Role::Side {
                    segment: a[(sa + k) % n].ids[0].clone(),
                    path: None,
                },
            })
            .collect();
        let shell = Shell::new(faces)?.named_loft(sides);
        self.leaf(crate::geometry::Leaf::Shell(Box::new(shell)), vec![], dec)
            .map(Some)
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
        // A whole-circle segment starts at local +x. Merged arc loops can
        // carry a different start and need the station loft's phase check.
        if sections.iter().any(|section| {
            !matches!(authored_loops(recipe, &section.profile), [lp]
                if matches!(lp.0.as_slice(), [catallax_recipe::Segment::Circle { id, .. }]
                    if id == &section.start[0]))
        }) {
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
