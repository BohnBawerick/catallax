//! Definition-based profile sweeps, independent of B-rep construction.

#[path = "extrude_algebraic.rs"]
mod extrude_algebraic;

use crate::{
    Face, Location, geometry,
    metric::Scalar,
    profile::{self, Edge, P, Profile},
    support::Support,
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{Axis, Origin, Role, StepId};
use catallax_shapes::{
    Cone, Cylinder, EllipticCylinder, Isometry, Plane, Point3, RPoint3, RVector3, Sphere, Surface,
    Torus, Vector3,
};

mod revolution;
use revolution::{Angle, radius};

#[derive(Clone, Debug)]
enum Sweep {
    Extrude {
        vector: RVector3,
        through: bool,
    },
    Revolve {
        origin: RPoint3,
        axis: RVector3,
        angle: Box<Angle>,
    },
}
#[derive(Clone, Debug)]
pub(crate) struct Feature {
    loops: Vec<Vec<Edge>>,
    convex_revolve: Option<(Sign, i64)>,
    inverse: Isometry,
    sweep: Sweep,
    // Irrational finite end parameter over a rational extrusion carrier.
    algebraic_end: Option<Real>,
    axis_norm: Rational,
    merge_caps: Vec<usize>,
    supports: Vec<Support>,
    pub(crate) constraints: Vec<geometry::Constraint>,
    fast: FastProfile,
    pub(crate) fillet: Option<Box<catallax_features::fillet::FilletDefinition>>,
    pub(crate) chamfer: Option<Box<catallax_features::chamfer::ChamferDefinition>>,
}
fn bad(code: &'static str) -> Refusal {
    Refusal::new(
        code,
        "feature definition requires an unsupported exact surface or frame",
    )
}
fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    let [x, y, z] = v.0.each_ref().map(|n| n.value().clone());
    RVector3::new(x, y, z)
}
fn point(p: &P) -> RPoint3 {
    RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero())
}
fn z() -> RVector3 {
    RVector3::new(Rational::zero(), Rational::zero(), Rational::one())
}
fn face(surface: Surface, step: &StepId, role: Role, dec: &mut Decider) -> Result<Face, Refusal> {
    geometry::face(
        surface,
        Origin {
            step: step.clone(),
            role,
        },
        dec,
    )
}
fn side(segment: &catallax_recipe::ElementId) -> Role {
    Role::Side {
        segment: segment.clone(),
        path: None,
    }
}

impl Feature {
    /// Collects the finite planar sweep skin from field definitions, without
    /// constructing or consulting B-rep topology.
    pub(crate) fn planar_skin(
        &self,
        origins: &[&Origin],
        through_bound: Option<&Rational>,
        dec: &mut Decider,
    ) -> Result<catallax_features::draft::arrangement::Skin, Refusal> {
        use catallax_features::draft::arrangement::{Patch, Skin};
        if self.algebraic_end.is_some() {
            return Err(bad("draft.unsupported-arrangement"));
        }
        let Sweep::Extrude { vector, through } = &self.sweep else {
            return Err(bad("draft.unsupported-arrangement"));
        };
        let (offset, extent) = if *through {
            let bound = through_bound.ok_or_else(|| bad("draft.unsupported-arrangement"))?;
            let mut height = Rational::zero();
            for corner in 0..8 {
                let [x, y, z] = [0, 1, 2].map(|axis| {
                    if corner & (1 << axis) == 0 {
                        -bound
                    } else {
                        bound.clone()
                    }
                });
                let local = self.inverse.apply_point(&RPoint3::new(x, y, z).to_real());
                let ball = geometry::enclose(&local.z, 128, dec)?;
                height = height.max(ball.lo().abs().max(ball.hi().abs()));
            }
            let k = (height + Rational::one()) / vector.z.abs();
            (vector.scaled(&(-&k)), vector.scaled(&(&k + &k)))
        } else {
            (RVector3::zero(), vector.clone())
        };
        let vector = &extent;
        if self.fillet.is_some() || self.loops.iter().flatten().any(|e| e.circle.is_some()) {
            return Err(bad("draft.unsupported-arrangement"));
        }
        fn source_role(origin: &Origin) -> &Role {
            let mut role = &origin.role;
            while let Role::Copy { of, .. } = role {
                role = &of.role;
            }
            role
        }
        let named = |role: &Role| -> Vec<Origin> {
            origins
                .iter()
                .filter(|o| source_role(o) == role)
                .map(|o| (**o).clone())
                .collect()
        };
        let placement = self.inverse.inverse();
        let sweep = placement.apply_vector(&vector.to_real());
        let mut vertices = Vec::new();
        let mut patches = Vec::new();
        let mut caps = Vec::new();
        for (index, edges) in self.loops.iter().enumerate() {
            let mut depth = 0;
            for (j, other) in self.loops.iter().enumerate() {
                if index != j
                    && profile::loop_location(other, &profile::real(&edges[0].start), dec)?
                        == Location::Inside
                {
                    depth += 1;
                }
            }
            let area = edges.iter().fold(Rational::zero(), |sum, e| {
                sum + profile::cross(&e.start, &e.end)
            });
            let desired = if depth % 2 == 0 {
                Sign::Positive
            } else {
                Sign::Negative
            };
            let follows = (area.sign() == desired) == (vector.z.sign() == Sign::Positive);
            let start = vertices.len();
            for edge in edges {
                vertices.push(placement.apply_point(&point(&edge.start).plus(&offset).to_real()));
            }
            for edge in edges {
                vertices.push(
                    placement.apply_point(&point(&edge.start).plus(&offset).plus(vector).to_real()),
                );
            }
            let mut cap: Vec<_> = (start..start + edges.len()).collect();
            if area.sign() != desired {
                cap.reverse();
            }
            caps.push((depth, index, cap));
            for (i, edge) in edges.iter().enumerate() {
                let j = (i + 1) % edges.len();
                let mut normal = point(&edge.end)
                    .minus(&point(&edge.start))
                    .cross(vector)
                    .to_real();
                if !follows {
                    normal = normal.scaled(&Real::from_i64(-1));
                }
                let names: Vec<Origin> = origins.iter().filter(|o|matches!(source_role(o),Role::Side{segment,..} if edge.ids.contains(segment))).map(|o|(**o).clone()).collect();
                if names.is_empty() {
                    return Err(bad("draft.missing-source-origin"));
                }
                patches.push(Patch {
                    plane: Plane {
                        origin: vertices[start + i].clone(),
                        normal: placement.apply_vector(&normal),
                        ref_dir: sweep.clone(),
                    },
                    loops: vec![vec![
                        start + i,
                        start + j,
                        start + edges.len() + j,
                        start + edges.len() + i,
                    ]],
                    origins: names,
                });
            }
        }
        caps.sort_by_key(|(depth, index, _)| (*depth, *index));
        let normal = placement.apply_vector(&Vector3::new(
            Real::zero(),
            Real::zero(),
            Real::from_i64(if vector.z.sign() == Sign::Positive {
                1
            } else {
                -1
            }),
        ));
        let reference =
            placement.apply_vector(&Vector3::new(Real::one(), Real::zero(), Real::zero()));
        let origin = placement.apply_point(&Point3::origin().plus(&offset.to_real()));
        for (depth, index, outer) in &caps {
            if depth % 2 != 0 {
                continue;
            }
            let mut start_loops = vec![outer.clone()];
            for (hole_depth, hole_index, hole) in &caps {
                if *hole_depth == depth + 1
                    && profile::loop_location(
                        &self.loops[*index],
                        &profile::real(&self.loops[*hole_index][0].start),
                        dec,
                    )? == Location::Inside
                {
                    start_loops.push(hole.clone());
                }
            }
            let end_loops = start_loops
                .iter()
                .map(|cycle| {
                    let count = cycle.len();
                    cycle.iter().map(|i| i + count).collect()
                })
                .collect();
            patches.push(Patch {
                plane: Plane {
                    origin: origin.clone(),
                    normal: normal.scaled(&Real::from_i64(-1)),
                    ref_dir: reference.clone(),
                },
                loops: start_loops,
                origins: named(&Role::Start),
            });
            patches.push(Patch {
                plane: Plane {
                    origin: origin.plus(&sweep),
                    normal: normal.clone(),
                    ref_dir: reference.clone(),
                },
                loops: end_loops,
                origins: named(&Role::End),
            });
        }
        if !through && patches.iter().any(|p| p.origins.is_empty()) {
            return Err(bad("draft.missing-source-origin"));
        }
        Ok(Skin { vertices, patches })
    }

    pub(crate) fn draft_witness(&self) -> Result<Point3, Refusal> {
        let Sweep::Extrude {
            vector,
            through: false,
        } = &self.sweep
        else {
            return Err(bad("draft.unsupported-arrangement"));
        };
        let [edges] = self.loops.as_slice() else {
            return Err(bad("draft.unsupported-arrangement"));
        };
        let centre = if let [edge] = edges.as_slice() {
            edge.circle.as_ref().map(|(centre, _, _)| centre.clone())
        } else {
            None
        };
        let centre = centre.unwrap_or_else(|| {
            let count = Rational::from_i64(edges.len() as i64);
            std::array::from_fn(|i| {
                edges
                    .iter()
                    .fold(Rational::zero(), |sum, edge| sum + edge.start[i].clone())
                    / count.clone()
            })
        });
        let half = Rational::from_ratio(1, 2).expect("nonzero denominator");
        let local = point(&centre)
            .to_real()
            .plus(&self.extrusion_vector(vector).scaled(&Real::from(half)));
        Ok(self.inverse.inverse().apply_point(&local))
    }

    pub(crate) fn certify_planar_convexity(&self, _dec: &mut Decider) -> Result<(), Refusal> {
        let Sweep::Extrude { through: false, .. } = self.sweep else {
            return Err(bad("draft.unsupported-arrangement"));
        };
        let [edges] = self.loops.as_slice() else {
            return Err(bad("draft.unsupported-arrangement"));
        };
        for edge in edges {
            if edge.circle.is_some() {
                return Err(bad("draft.unsupported-arrangement"));
            }
            let direction = [&edge.end[0] - &edge.start[0], &edge.end[1] - &edge.start[1]];
            let mut side = Sign::Zero;
            for other in edges {
                let cross = &direction[0] * &(&other.start[1] - &edge.start[1])
                    - &direction[1] * &(&other.start[0] - &edge.start[0]);
                let sign = cross.sign();
                if sign == Sign::Zero {
                    continue;
                }
                if side != Sign::Zero && side != sign {
                    return Err(bad("draft.unsupported-arrangement"));
                }
                side = sign;
            }
        }
        Ok(())
    }

    pub(crate) fn convex_planes(&self, dec: &mut Decider) -> Result<Option<Vec<Plane>>, Refusal> {
        let Sweep::Extrude {
            vector,
            through: false,
        } = &self.sweep
        else {
            return Ok(None);
        };
        if self.fillet.is_some()
            || self.loops.len() != 1
            || self.loops[0].iter().any(|e| e.circle.is_some())
        {
            return Ok(None);
        }
        let map = self.inverse.inverse();
        let mut vertices = Vec::with_capacity(self.loops[0].len() * 2);
        let mut centre = Vector3::zero();
        for edge in &self.loops[0] {
            let p = RPoint3::new(
                edge.start[0].clone(),
                edge.start[1].clone(),
                Rational::zero(),
            );
            let p = p.to_real();
            let end = p.plus(&self.extrusion_vector(vector));
            for p in [p, end] {
                let p = map.apply_point(&p);
                centre = centre.plus(&Vector3::new(p.x.clone(), p.y.clone(), p.z.clone()));
                vertices.push(p);
            }
        }
        let centre = centre.scaled(&(&Real::one() / &Real::from_i64(vertices.len() as i64)));
        let centre = Point3::new(centre.x, centre.y, centre.z);
        let mut planes = Vec::with_capacity(self.constraints.len());
        for constraint in &self.constraints {
            let Surface::Plane(plane) = constraint.surface() else {
                return Ok(None);
            };
            let mut plane = plane.clone();
            match geometry::sign(&centre.minus(&plane.origin).dot(&plane.normal), dec)? {
                Sign::Positive => plane.normal = plane.normal.scaled(&Real::from_i64(-1)),
                Sign::Zero => return Ok(None),
                Sign::Negative => (),
            }
            for p in &vertices {
                if geometry::sign(&p.minus(&plane.origin).dot(&plane.normal), dec)?
                    == Sign::Positive
                {
                    return Ok(None);
                }
            }
            planes.push(plane);
        }
        Ok(Some(planes))
    }
    pub(crate) fn display_box(&self, dec: &mut Decider) -> Result<Option<[Ball; 3]>, Refusal> {
        if self.algebraic_end.is_some() {
            return self.algebraic_display_box(dec);
        }
        let Sweep::Extrude { vector, through } = &self.sweep else {
            let radius = self.display_radius(dec)?.expect("finite revolve");
            return Ok(Some(std::array::from_fn(|_| {
                Ball::from_rational(&(-&radius), 128).hull(&Ball::from_rational(&radius, 128))
            })));
        };
        if *through {
            return Ok(None);
        }
        let mut bounds: Option<[[Rational; 2]; 3]> = None;
        for edge in self.loops.iter().flatten() {
            for point in [&edge.start, &edge.end] {
                let b = bounds.get_or_insert_with(|| {
                    [
                        [point[0].clone(), point[0].clone()],
                        [point[1].clone(), point[1].clone()],
                        [Rational::zero(), Rational::zero()],
                    ]
                });
                for i in 0..2 {
                    b[i][0] = b[i][0].clone().min(point[i].clone());
                    b[i][1] = b[i][1].clone().max(point[i].clone());
                }
            }
            if let Some((centre, radius, _)) = &edge.circle {
                let radius = Ball::from_rational(radius, 128)
                    .sqrt()
                    .expect("positive profile radius")
                    .hi();
                let b = bounds.as_mut().expect("edge endpoints");
                for i in 0..2 {
                    b[i][0] = b[i][0].clone().min(&centre[i] - &radius);
                    b[i][1] = b[i][1].clone().max(&centre[i] + &radius);
                }
            }
        }
        let Some(mut bounds) = bounds else {
            return Ok(None);
        };
        for (i, v) in [&vector.x, &vector.y, &vector.z].into_iter().enumerate() {
            bounds[i][0] = &bounds[i][0] + &v.clone().min(Rational::zero());
            bounds[i][1] = &bounds[i][1] + &v.clone().max(Rational::zero());
        }
        Ok(Some(geometry::display_box(
            &self.inverse.inverse(),
            &bounds,
            dec,
        )?))
    }
    pub(crate) fn display_radius(&self, dec: &mut Decider) -> Result<Option<Rational>, Refusal> {
        if self.algebraic_end.is_some() {
            return Ok(self.algebraic_display_box(dec)?.map(|bounds| {
                bounds
                    .iter()
                    .fold(Rational::zero(), |r, b| r + b.lo().abs().max(b.hi().abs()))
            }));
        }
        let mut radius = Rational::zero();
        for edge in self.loops.iter().flatten() {
            for p in [&edge.start, &edge.end] {
                radius = radius.max(p[0].abs() + p[1].abs());
            }
            if let Some((c, r, _)) = &edge.circle {
                let r = catallax_numbers::Ball::from_rational(r, 128)
                    .sqrt()
                    .expect("positive profile radius")
                    .hi();
                radius = radius.max(c[0].abs() + c[1].abs() + r);
            }
        }
        radius = match &self.sweep {
            Sweep::Extrude { through: true, .. } => return Ok(None),
            Sweep::Extrude { vector, .. } => {
                radius + vector.x.abs() + vector.y.abs() + vector.z.abs()
            }
            Sweep::Revolve { origin, .. } => {
                radius + Rational::from_i64(2) * (origin.x.abs() + origin.y.abs() + origin.z.abs())
            }
        };
        let translation = self.inverse.inverse().translation;
        for x in [&translation.x, &translation.y, &translation.z] {
            let b = geometry::enclose(x, 128, dec)?;
            radius = radius + b.lo().abs().max(b.hi().abs());
        }
        Ok(Some(radius))
    }
    pub fn extrude(
        profile: &Profile,
        step: &StepId,
        distance: Option<&Rational>,
        direction: Option<&catallax_recipe::Vector3>,
        empty_target: bool,
        dec: &mut Decider,
    ) -> Result<Option<(Self, Vec<Face>)>, Refusal> {
        if distance.is_some_and(Rational::is_zero) {
            return Err(bad("extrude.zero-distance"));
        }
        let normal = profile::rational_vector(&profile.pose, &z());
        let direction = direction.map(vector).unwrap_or(normal.clone());
        if direction.is_zero() {
            return Err(bad("extrude.zero-direction"));
        }
        if direction.dot(&normal).is_zero() {
            return Err(bad("extrude.zero-thickness"));
        }
        if profile.loops.is_empty() || empty_target {
            return Ok(None);
        }
        let through = distance.is_none();
        let v = if let Some(distance) = distance {
            let Some(length) = profile::sqrt(&direction.norm_squared()) else {
                return Self::algebraic_extrude(profile, step, distance, &direction, dec);
            };
            direction.scaled(&(distance / &length))
        } else {
            direction
        };
        let local = profile::rational_vector(&profile.pose.inverse(), &v);
        let mut faces = Vec::new();
        for edge in profile.loops.iter().flatten() {
            let surface = if let Some((c, r, _)) = &edge.circle {
                let major_dir = normal.cross(&v);
                if major_dir.is_zero() {
                    Surface::Cylinder(Cylinder {
                        origin: profile.world(c).into(),
                        axis: normal.clone().into(),
                        radius_sq: r.clone().into(),
                        ref_dir: profile.world(&edge.start).minus(&profile.world(c)).into(),
                    })
                } else {
                    let height = v.dot(&normal);
                    Surface::EllipticCylinder(EllipticCylinder {
                        origin: profile.world(c).into(),
                        axis: v.clone().into(),
                        major_dir: major_dir.into(),
                        major_sq: r.clone().into(),
                        minor_sq: (r * &(&height * &height) / v.norm_squared()).into(),
                    })
                }
            } else {
                Surface::Plane(Plane {
                    origin: profile.world(&edge.start).into(),
                    normal: profile
                        .world(&edge.end)
                        .minus(&profile.world(&edge.start))
                        .cross(&v)
                        .into(),
                    ref_dir: v.clone().into(),
                })
            };
            for id in &edge.ids {
                faces.push(face(surface.clone(), step, side(id), dec)?);
            }
        }
        if !through {
            let origin = profile::rational_point(&profile.pose, &RPoint3::origin());
            let reference = profile::rational_vector(
                &profile.pose,
                &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()),
            );
            for (origin, role) in [(origin.clone(), Role::Start), (origin.plus(&v), Role::End)] {
                faces.push(face(
                    Surface::Plane(Plane {
                        origin: origin.into(),
                        normal: normal.clone().into(),
                        ref_dir: reference.clone().into(),
                    }),
                    step,
                    role,
                    dec,
                )?);
            }
        }
        let sweep = Sweep::Extrude {
            vector: local,
            through,
        };
        Ok(Some((
            Self::new(profile, sweep, &faces, Rational::zero(), Vec::new())?,
            faces,
        )))
    }
    pub fn revolve(
        profile: &Profile,
        step: &StepId,
        axis: &Axis,
        angle: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<(Self, Vec<Face>)>, Refusal> {
        if angle.sign() != Sign::Positive || *angle > Rational::from_i64(360) {
            return Err(bad("revolve.angle"));
        }
        let world_axis = vector(&axis.direction);
        if world_axis.is_zero() {
            return Err(bad("revolve.zero-axis"));
        }
        let inverse = profile.pose.inverse();
        let origin = profile::rational_point(
            &inverse,
            &RPoint3::new(
                axis.point.0[0].value().clone(),
                axis.point.0[1].value().clone(),
                axis.point.0[2].value().clone(),
            ),
        );
        let axis = profile::rational_vector(&inverse, &world_axis);
        if !origin.z.is_zero() || !axis.z.is_zero() {
            return Err(bad("revolve.axis-plane"));
        }
        if profile.loops.is_empty() {
            return Ok(None);
        }
        let norm = axis.norm_squared();
        let project =
            |p: &RPoint3| origin.plus(&axis.scaled(&(p.minus(&origin).dot(&axis) / norm.clone())));
        let signed = |p: &RPoint3| axis.cross(&p.minus(&origin)).z;
        let mut depths = Vec::new();
        for (i, edges) in profile.loops.iter().enumerate() {
            let mut depth = 0;
            for (j, other) in profile.loops.iter().enumerate() {
                if i != j
                    && profile::loop_location(other, &profile::real(&edges[0].start), dec)?
                        == Location::Inside
                {
                    depth += 1;
                }
            }
            depths.push(depth);
        }
        let full = *angle == Rational::from_i64(360);
        let on_axis = |e: &Edge| {
            e.circle.is_none()
                && signed(&point(&e.start)).is_zero()
                && signed(&point(&e.end)).is_zero()
        };
        for (_, edges) in profile
            .loops
            .iter()
            .enumerate()
            .filter(|(i, _)| depths[*i] % 2 == 0)
        {
            let (mut positive, mut negative) = (false, false);
            let mut pinches_axis = false;
            for (j, e) in edges.iter().enumerate() {
                let mut values = vec![
                    Real::from(signed(&point(&e.start))),
                    Real::from(signed(&point(&e.end))),
                ];
                if let Some((c, r, _)) = &e.circle {
                    let scale = profile::positive_sqrt(&(r / &norm), dec)?;
                    for flip in [-1, 1] {
                        let scale = &Real::from_i64(flip) * &scale;
                        let p = [
                            Real::from(c[0].clone()) - &Real::from(axis.y.clone()) * &scale,
                            Real::from(c[1].clone()) + &Real::from(axis.x.clone()) * &scale,
                        ];
                        if e.contains(&p, dec)? {
                            values.push(
                                &Real::from(axis.x.clone())
                                    * &(&p[1] - &Real::from(origin.y.clone()))
                                    - &Real::from(axis.y.clone())
                                        * &(&p[0] - &Real::from(origin.x.clone())),
                            );
                        }
                    }
                    let distance = signed(&point(c));
                    if full && &distance * &distance == r * &norm {
                        let touch = project(&point(c));
                        let p = [touch.x, touch.y];
                        pinches_axis |= (e.full || (p != e.start && p != e.end))
                            && e.contains(&profile::real(&p), dec)?;
                    }
                }
                for value in values {
                    match profile::sign(&value, dec)? {
                        Sign::Positive => positive = true,
                        Sign::Negative => negative = true,
                        _ => (),
                    }
                }
                pinches_axis |= full
                    && !e.full
                    && signed(&point(&e.start)).is_zero()
                    && !on_axis(e)
                    && !on_axis(&edges[(j + edges.len() - 1) % edges.len()]);
            }
            if positive == negative {
                return Err(bad("revolve.axis-interior"));
            }
            if pinches_axis {
                return Err(bad("boolean.self-touching-shell"));
            }
        }
        let angle = Angle::new(angle, dec)?;
        let mut faces = Vec::new();
        let mut merge_caps = Vec::new();
        for (index, edges) in profile.loops.iter().enumerate() {
            for edge in edges {
                let start = point(&edge.start);
                let end = point(&edge.end);
                let radial = start.minus(&project(&start));
                let surface = if let Some((c, r, _)) = &edge.circle {
                    let centre = point(c);
                    let projection = project(&centre);
                    let offset = centre.minus(&projection);
                    let radius = radius(r);
                    if offset.is_zero() {
                        Surface::Sphere(Sphere {
                            center: centre.into(),
                            radius,
                            axis: axis.clone().into(),
                            ref_dir: axis.cross(&z()).into(),
                        })
                    } else {
                        let major_sq = offset.norm_squared();
                        if major_sq <= *r {
                            return Err(bad("features.unsupported-spindle-torus"));
                        }
                        Surface::Torus(Torus {
                            center: projection.into(),
                            axis: axis.clone().into(),
                            major: revolution::radius(&major_sq),
                            minor: radius,
                            ref_dir: offset.into(),
                        })
                    }
                } else {
                    let tangent = end.minus(&start);
                    if on_axis(edge) {
                        if angle.is(180) && depths[index] % 2 == 0 {
                            merge_caps.push(index);
                        }
                        continue;
                    }
                    let axial = tangent.dot(&axis);
                    if axial.is_zero() {
                        Surface::Plane(Plane {
                            origin: project(&start).into(),
                            normal: axis.clone().into(),
                            ref_dir: if radial.is_zero() {
                                end.minus(&project(&end)).into()
                            } else {
                                radial.into()
                            },
                        })
                    } else if tangent.cross(&axis).is_zero() {
                        Surface::Cylinder(Cylinder {
                            origin: project(&start).into(),
                            axis: axis.clone().into(),
                            radius_sq: radial.norm_squared().into(),
                            ref_dir: radial.into(),
                        })
                    } else {
                        let change = axis.cross(&tangent).z;
                        let apex = start.plus(&tangent.scaled(&(-signed(&start) / change.clone())));
                        let direction = if start.minus(&apex).dot(&axis).sign() == Sign::Positive
                            || (start == apex
                                && end.minus(&apex).dot(&axis).sign() == Sign::Positive)
                        {
                            axis.clone()
                        } else {
                            axis.scaled(&Rational::from_i64(-1))
                        };
                        Surface::Cone(Cone {
                            apex: apex.into(),
                            axis: direction.into(),
                            slope: (change / axial).abs().into(),
                            ref_dir: if radial.is_zero() {
                                end.minus(&project(&end)).into()
                            } else {
                                radial.into()
                            },
                        })
                    }
                };
                let surface = surface.transformed(&profile.pose);
                for id in &edge.ids {
                    faces.push(face(surface.clone(), step, side(id), dec)?);
                }
            }
        }
        if !full {
            let normal = profile::rational_vector(&profile.pose, &z());
            let end_normal = profile.pose.apply_vector(&angle.end_normal(&axis));
            let reference = world_axis.clone();
            let world_origin = profile::rational_point(&profile.pose, &origin);
            faces.push(face(
                Surface::Plane(Plane {
                    origin: world_origin.clone().into(),
                    normal: normal.into(),
                    ref_dir: reference.clone().into(),
                }),
                step,
                Role::Start,
                dec,
            )?);
            faces.push(face(
                Surface::Plane(Plane {
                    origin: world_origin.into(),
                    normal: end_normal,
                    ref_dir: reference.into(),
                }),
                step,
                Role::End,
                dec,
            )?);
        }
        let sweep = Sweep::Revolve {
            origin,
            axis,
            angle: Box::new(angle),
        };
        Ok(Some((
            Self::new(profile, sweep, &faces, norm, merge_caps)?,
            faces,
        )))
    }
    fn new(
        profile: &Profile,
        sweep: Sweep,
        faces: &[Face],
        axis_norm: Rational,
        merge_caps: Vec<usize>,
    ) -> Result<Self, Refusal> {
        let mut constraints = Vec::new();
        let mut supports = Vec::new();
        for face in faces {
            if !constraints
                .iter()
                .any(|c: &geometry::Constraint| c.polynomial.same_definition(&face.polynomial))
            {
                constraints.push(geometry::constraint(face.surface.clone())?);
                supports.push(Support::new(&face.surface));
            }
        }
        Ok(Self {
            loops: profile.loops.clone(),
            convex_revolve: convex_revolve(profile, &sweep),
            inverse: Isometry::from(&profile.pose.inverse()),
            fast: FastProfile::new(profile, &sweep, &axis_norm),
            axis_norm,
            merge_caps,
            sweep,
            algebraic_end: None,
            constraints,
            supports,
            fillet: None,
            chamfer: None,
        })
    }
    pub fn transform(&mut self, pose: &Isometry) -> Result<(), Refusal> {
        if let Some(definition) = &mut self.fillet {
            definition.transform(pose)?;
        }
        if let Some(definition) = &mut self.chamfer {
            definition.transform(pose)?;
        }
        self.inverse = pose.inverse().then(&self.inverse);
        self.fast.matrix = geometry::fast_matrix(&self.inverse);
        self.fast.translation = geometry::fast_translation(&self.inverse);
        self.fast.valid = self.fast.finite();
        for support in &mut self.supports {
            support.transform(pose);
        }
        for constraint in &mut self.constraints {
            constraint.transform(pose)?;
        }
        Ok(())
    }
    pub(crate) fn display_classifier(&self) -> Option<crate::display::membership::Classifier<'_>> {
        if self.algebraic_end.is_some() {
            return None;
        }
        let Sweep::Extrude { vector, through } = &self.sweep else {
            return None;
        };
        crate::display::membership::Classifier::new(&self.inverse, vector, *through, &self.loops)
    }
    // Exact stadium section of a plain extrusion, for the preview's obstruction
    // proof. Blended, chamfered and other swept definitions have none.
    pub(crate) fn display_capsule(&self) -> Option<crate::display::obstruction::Capsule<'_>> {
        if self.algebraic_end.is_some() || self.fillet.is_some() || self.chamfer.is_some() {
            return None;
        }
        let Sweep::Extrude { vector, through } = &self.sweep else {
            return None;
        };
        crate::display::obstruction::Capsule::new(&self.inverse, vector, *through, &self.loops)
    }
    // Query-only reuse of certified elementary fillet sectors. Unsupported
    // sectors and non-rational ray directions retain the existing evaluator.
    fn fillet_query_classifier(&self) -> Option<crate::display::membership::Classifier<'_>> {
        self.fillet.as_deref().and_then(|definition| {
            crate::display::membership::Classifier::fillet(definition)
                .or_else(|| self.display_classifier())
        })
    }
    pub fn location(&self, p: &[Rational; 3], dec: &mut Decider) -> Result<Location, Refusal> {
        if let Some(classifier) = self.fillet_query_classifier() {
            return Ok(classifier.location(p));
        }
        self.classify(
            &Point3::new(
                p[0].clone().into(),
                p[1].clone().into(),
                p[2].clone().into(),
            ),
            dec,
        )
    }
    pub(crate) fn classify(&self, p: &Point3, dec: &mut Decider) -> Result<Location, Refusal> {
        if let Some(p) = p.as_rational()
            && let Some(classifier) = self.fillet_query_classifier()
        {
            return Ok(classifier.location(&[p.x, p.y, p.z]));
        }
        let local = self.inverse.apply_point(p);
        if let Some((location, _)) = self.convex_revolve_query(&local, dec)? {
            return Ok(location);
        }
        match &self.sweep {
            Sweep::Extrude { vector, through } => {
                let t = &local.z / &Real::from(vector.z.clone());
                let cap_boundary = if *through {
                    false
                } else {
                    let start = geometry::sign(&t, dec)?;
                    if start == Sign::Negative {
                        return Ok(Location::Outside);
                    }
                    let end = geometry::sign(&(&t - &self.extrusion_end()), dec)?;
                    if end == Sign::Positive {
                        return Ok(Location::Outside);
                    }
                    start == Sign::Zero || end == Sign::Zero
                };
                let q = [
                    &local.x - &(&t * &Real::from(vector.x.clone())),
                    &local.y - &(&t * &Real::from(vector.y.clone())),
                ];
                let loc = profile::location(&self.loops, &q, dec)?;
                Ok(if loc == Location::Inside && cap_boundary {
                    Location::Boundary
                } else {
                    loc
                })
            }
            Sweep::Revolve {
                origin,
                axis,
                angle,
            } => {
                let axis_real = axis.to_real();
                let origin_real = origin.to_real();
                let norm = Real::from(self.axis_norm.clone());
                let v = local.minus(&origin_real);
                let t = &v.dot(&axis_real) / &norm;
                let b = &axis_real.cross(&v).z / &norm;
                let k = sqrt_real(&(&b * &b + &(&v.z * &v.z) / &norm), dec)?;
                let b_sign = geometry::sign(&b, dec)?;
                let z_sign = geometry::sign(&v.z, dec)?;
                let mut boundary = false;
                for flip in [-1, 1] {
                    let angular =
                        angle.location(flip, [&b, &v.z], &self.axis_norm, [b_sign, z_sign], dec)?;
                    if angular == Location::Outside {
                        continue;
                    }
                    let k = &Real::from_i64(flip) * &k;
                    let q = [
                        &origin_real.x + &(&axis_real.x * &t) - &axis_real.y * &k,
                        &origin_real.y + &(&axis_real.y * &t) + &axis_real.x * &k,
                    ];
                    let loc = profile::location(&self.loops, &q, dec)?;
                    if loc == Location::Inside && angular == Location::Inside {
                        return Ok(Location::Inside);
                    }
                    boundary |= loc != Location::Outside;
                }
                if boundary && angle.is(360) && b_sign == Sign::Zero && z_sign == Sign::Zero {
                    let q = [local.x, local.y];
                    for edge in self.loops.iter().flatten().filter(|e| e.circle.is_none()) {
                        if axis.cross(&point(&edge.start).minus(origin)).z.is_zero()
                            && axis.cross(&point(&edge.end).minus(origin)).z.is_zero()
                            && !same_profile_point(&q, &edge.start, dec)?
                            && !same_profile_point(&q, &edge.end, dec)?
                            && edge.contains(&q, dec)?
                        {
                            return Ok(Location::Inside);
                        }
                    }
                }
                Ok(if boundary {
                    Location::Boundary
                } else {
                    Location::Outside
                })
            }
        }
    }
    // A full revolution of a convex line polygon is an intersection of
    // oriented section half-planes evaluated at its unique radial preimage.
    // Compare A*sqrt(D)+B by sign and squaring, avoiding root isolation at
    // every rounded display vertex. The convexity and one-sided-axis proofs
    // are cached from the immutable profile, and other profiles use the
    // general winding evaluator.
    fn convex_revolve_query(
        &self,
        local: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<(Location, Vec<usize>)>, Refusal> {
        let Some((winding, radial_sign)) = self.convex_revolve else {
            return Ok(None);
        };
        dec.check_deadline()?;
        let Sweep::Revolve { origin, axis, .. } = &self.sweep else {
            return Ok(None);
        };
        let origin_real = origin.to_real();
        let axis_real = axis.to_real();
        let norm = Real::from(self.axis_norm.clone());
        let v = local.minus(&origin_real);
        let t = v.dot(&axis_real) / norm.clone();
        let base = origin_real.plus(&axis_real.scaled(&t));
        let cross = axis_real.cross(&v).z;
        let radicand = cross.square() / norm.square() + local.z.square() / norm;
        let mut boundary = Vec::new();
        for (i, edge) in self.loops[0].iter().enumerate() {
            let d = profile::sub(&edge.end, &edge.start);
            let a =
                Real::from((&d[0] * &axis.x + &d[1] * &axis.y) * Rational::from_i64(radial_sign));
            let b = &Real::from(d[0].clone()) * &(&base.y - &Real::from(edge.start[1].clone()))
                - &Real::from(d[1].clone()) * &(&base.x - &Real::from(edge.start[0].clone()));
            let side = radical_sign(&a, &b, &radicand, dec)?;
            if side != winding && side != Sign::Zero {
                return Ok(Some((Location::Outside, Vec::new())));
            }
            if side == Sign::Zero {
                let axis_edge = [&edge.start, &edge.end]
                    .iter()
                    .all(|p| axis.cross(&point(p).minus(origin)).z.is_zero());
                if !axis_edge {
                    boundary.push(i);
                }
            }
        }
        Ok(Some((
            if boundary.is_empty() {
                Location::Inside
            } else {
                Location::Boundary
            },
            boundary,
        )))
    }

    pub fn on_face(
        &self,
        role: &Role,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut base_role = role;
        while let Role::Copy { of, .. } = base_role {
            base_role = &of.role;
        }
        if self
            .fillet
            .as_ref()
            .is_some_and(|definition| &definition.origin.role == base_role)
            && let Some(classifier) = self.fillet_query_classifier()
        {
            return Ok(classifier.on_face(base_role, p));
        }
        self.on_face_point(
            role,
            &Point3::new(
                p[0].clone().into(),
                p[1].clone().into(),
                p[2].clone().into(),
            ),
            dec,
        )
    }
    pub(crate) fn on_face_point(
        &self,
        role: &Role,
        p: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut role = role;
        while let Role::Copy { of, .. } = role {
            role = &of.role;
        }
        if self.fillet.as_ref().is_some_and(|f| &f.origin.role == role)
            && let Some(p) = p.as_rational()
            && let Some(classifier) = self.fillet_query_classifier()
        {
            return Ok(classifier.on_face(role, &[p.x, p.y, p.z]));
        }
        let alias = self
            .fillet
            .as_ref()
            .filter(|f| &f.origin.role == role)
            .map(|_| Role::Side {
                segment: catallax_recipe::ElementId::new("blend").expect("fixed section id"),
                path: None,
            });
        let alias = alias.or_else(|| {
            self.chamfer
                .as_ref()
                .filter(|f| &f.origin.role == role)
                .map(|_| Role::Side {
                    segment: catallax_recipe::ElementId::new("bevel").expect("fixed section id"),
                    path: None,
                })
        });
        let role = alias.as_ref().unwrap_or(role);
        let local = self.inverse.apply_point(p);
        if let Role::Side { segment, .. } = role
            && let Some((_, boundary)) = self.convex_revolve_query(&local, dec)?
        {
            return Ok(boundary
                .into_iter()
                .any(|i| self.loops[0][i].ids.contains(segment)));
        }
        let norm = &self.axis_norm;
        let merge_caps = &self.merge_caps;
        match &self.sweep {
            Sweep::Extrude { vector, through } => {
                let t = &local.z / &Real::from(vector.z.clone());
                if !through
                    && (geometry::sign(&t, dec)? == Sign::Negative
                        || geometry::sign(&(&t - &self.extrusion_end()), dec)? == Sign::Positive)
                {
                    return Ok(false);
                }
                let q = [
                    &local.x - &(&t * &Real::from(vector.x.clone())),
                    &local.y - &(&t * &Real::from(vector.y.clone())),
                ];
                match role {
                    Role::Side { segment, .. } => {
                        for e in self
                            .loops
                            .iter()
                            .flatten()
                            .filter(|e| e.ids.contains(segment))
                        {
                            if e.contains(&q, dec)? {
                                return Ok(true);
                            }
                        }
                        Ok(false)
                    }
                    Role::Start | Role::End => {
                        Ok(profile::location(&self.loops, &q, dec)? != Location::Outside)
                    }
                    _ => Ok(false),
                }
            }
            Sweep::Revolve {
                origin,
                axis,
                angle,
            } => {
                let origin = origin.to_real();
                let axis = axis.to_real();
                let norm_real = Real::from(norm.clone());
                let v = local.minus(&origin);
                let t = &v.dot(&axis) / &norm_real;
                let b = &axis.cross(&v).z / &norm_real;
                if let Role::Side { segment, .. } = role {
                    let k = sqrt_real(&(&b * &b + &(&v.z * &v.z) / &norm_real), dec)?;
                    let b_sign = geometry::sign(&b, dec)?;
                    let z_sign = geometry::sign(&v.z, dec)?;
                    for flip in [-1, 1] {
                        if angle.location(flip, [&b, &v.z], norm, [b_sign, z_sign], dec)?
                            == Location::Outside
                        {
                            continue;
                        }
                        let k = &Real::from_i64(flip) * &k;
                        let q = [
                            &origin.x + &(&axis.x * &t) - &axis.y * &k,
                            &origin.y + &(&axis.y * &t) + &axis.x * &k,
                        ];
                        for e in self
                            .loops
                            .iter()
                            .flatten()
                            .filter(|e| e.ids.contains(segment))
                        {
                            if e.contains(&q, dec)? {
                                return Ok(true);
                            }
                        }
                    }
                    return Ok(false);
                }
                let radial = v.minus(&axis.scaled(&t));
                let rotated = if matches!(role, Role::Start) {
                    radial
                } else {
                    angle.unrotate_radial(&radial, &axis, norm)
                };
                if geometry::sign(&rotated.z, dec)? != Sign::Zero {
                    return Ok(false);
                }
                let q = origin.plus(&axis.scaled(&t)).plus(&rotated);
                let q = [q.x, q.y];
                if profile::location(&self.loops, &q, dec)? != Location::Outside {
                    return Ok(true);
                }
                if angle.is(180) {
                    let reflected = origin
                        .plus(&axis.scaled(&t))
                        .plus(&rotated.scaled(&Real::from_i64(-1)));
                    let reflected = [reflected.x, reflected.y];
                    if profile::location(&self.loops, &reflected, dec)? != Location::Outside {
                        for &i in merge_caps {
                            if profile::loop_location(&self.loops[i], &reflected, dec)?
                                != Location::Outside
                            {
                                return Ok(true);
                            }
                        }
                    }
                }
                Ok(false)
            }
        }
    }
    pub fn proven_distance(
        &self,
        p: &RPoint3,
        location: Location,
        dec: &mut Decider,
    ) -> Result<Ball, Refusal> {
        if location == Location::Boundary {
            return Ok(Ball::integer(0));
        }
        let p = [&p.x, &p.y, &p.z].map(|q| Ball::from_rational(q, crate::metric::PRECISION));
        let mut supports = self.supports.iter();
        let mut residual = supports
            .next()
            .expect("nonempty sweep supports")
            .proven(&p, dec)?;
        for support in supports {
            residual = residual.min(&support.proven(&p, dec)?);
        }
        Ok(match location {
            Location::Inside => -residual,
            Location::Indeterminate => residual.hull(&(-&residual)),
            _ => residual,
        })
    }
    pub fn ray(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        if let Some(classifier) = self.fillet_query_classifier()
            && let Some(location) = classifier.ray(&self.constraints, p, d)
        {
            return Ok(location);
        }
        // Stay within every support's initial sign cell; this is not a tolerance.
        let mut epsilon = Rational::one();
        for constraint in &self.constraints {
            if let Some(coefficients) = geometry::native_ray_coefficients(constraint, p, d) {
                if let Some((i, c)) = coefficients.iter().enumerate().find(|(_, c)| !c.is_zero()) {
                    let tail = coefficients[i + 1..]
                        .iter()
                        .fold(Rational::zero(), |a, c| a + c.abs());
                    let bound = c.abs() / (c.abs() + tail + Rational::one());
                    epsilon = epsilon.min(bound);
                }
            } else {
                let coefficients = geometry::algebraic_ray_coefficients(constraint, p, d);
                for (i, c) in coefficients.iter().enumerate() {
                    let sign = geometry::sign(c, dec)?;
                    if sign == Sign::Zero {
                        continue;
                    }
                    let magnitude = if sign == Sign::Negative {
                        -c
                    } else {
                        c.clone()
                    };
                    let mut denominator = &magnitude + &Real::one();
                    for tail in &coefficients[i + 1..] {
                        // |x| <= x^2 + 1, without a second sign decision.
                        denominator = denominator
                            + if let Some(q) = tail.rational_ref() {
                                Real::from(q.abs())
                            } else {
                                tail * tail + Real::one()
                            };
                    }
                    let bound = magnitude / denominator;
                    let mut precision = crate::metric::PRECISION;
                    let lower = loop {
                        let lower = geometry::enclose(&bound, precision, dec)?.lo();
                        if lower.sign() == Sign::Positive {
                            break lower;
                        }
                        precision = precision.checked_mul(2).ok_or_else(|| {
                            Refusal::new(
                                "numbers.refinement-budget",
                                "ray-germ separation exceeds the precision index",
                            )
                        })?;
                    };
                    epsilon = epsilon.min(lower);
                    break;
                }
            }
        }
        let p = Point3::new(
            p[0].clone().into(),
            p[1].clone().into(),
            p[2].clone().into(),
        );
        self.classify(&p.plus(&d.scaled(&Real::from(epsilon))), dec)
    }
    pub fn fast_distance(&self, p: [f64; 3]) -> f64 {
        if !self.fast.valid {
            return f64::NAN;
        }
        let inside = self.fast.inside(p);
        let mut residual = f64::INFINITY;
        for s in &self.supports {
            let distance = s.fast(&p);
            if !distance.is_finite() {
                return f64::NAN;
            }
            residual = residual.min(distance);
        }
        if inside { -residual } else { residual }
    }
}

// The constructor has already validated the loops. A single straight loop
// is convex exactly when every vertex lies in every oriented edge half-plane.
fn convex_revolve(profile: &Profile, sweep: &Sweep) -> Option<(Sign, i64)> {
    let Sweep::Revolve {
        origin,
        axis,
        angle,
    } = sweep
    else {
        return None;
    };
    if !angle.is(360) {
        return None;
    }
    let [edges] = profile.loops.as_slice() else {
        return None;
    };
    if edges.len() < 3 || edges.iter().any(|edge| edge.circle.is_some()) {
        return None;
    }
    let area = edges.iter().fold(Rational::zero(), |area, e| {
        area + profile::cross(&e.start, &e.end)
    });
    let winding = area.sign();
    if winding == Sign::Zero {
        return None;
    }
    let mut radial = None;
    for edge in edges {
        let d = profile::sub(&edge.end, &edge.start);
        for other in edges {
            let s = profile::cross(&d, &profile::sub(&other.start, &edge.start)).sign();
            if s != Sign::Zero && s != winding {
                return None;
            }
        }
        let s = axis.cross(&point(&edge.start).minus(origin)).z.sign();
        if s != Sign::Zero {
            if radial.is_some_and(|r| r != s) {
                return None;
            }
            radial = Some(s);
        }
    }
    Some((winding, if radial? == Sign::Positive { 1 } else { -1 }))
}

fn radical_sign(a: &Real, b: &Real, d: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    let sa = geometry::sign(a, dec)?;
    let sb = geometry::sign(b, dec)?;
    if sa == Sign::Zero {
        return Ok(sb);
    }
    if sb == Sign::Zero {
        return Ok(if geometry::sign(d, dec)? == Sign::Zero {
            Sign::Zero
        } else {
            sa
        });
    }
    if sa == sb {
        return Ok(sa);
    }
    Ok(
        match geometry::sign(&(&a.square() * d - b.square()), dec)? {
            Sign::Positive => sa,
            Sign::Negative => sb,
            Sign::Zero => Sign::Zero,
        },
    )
}

fn sqrt_real(n: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    if let Some(q) = n.rational_ref() {
        profile::positive_sqrt(q, dec)
    } else {
        Ok(n.sqrt())
    }
}

fn same_profile_point(a: &[Real; 2], b: &P, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        geometry::sign(&(&a[0] - &Real::from(b[0].clone())), dec)? == Sign::Zero
            && geometry::sign(&(&a[1] - &Real::from(b[1].clone())), dec)? == Sign::Zero,
    )
}

fn angular(angle: u16, flip: i64, b: Sign, z: Sign) -> Location {
    if angle == 360 {
        return Location::Inside;
    }
    let oriented = |s| {
        if flip == 1 {
            s
        } else {
            match s {
                Sign::Positive => Sign::Negative,
                Sign::Negative => Sign::Positive,
                Sign::Zero => Sign::Zero,
            }
        }
    };
    let (b, z) = (oriented(b), oriented(z));
    let outside = match angle {
        90 => b == Sign::Negative || z == Sign::Negative,
        180 => z == Sign::Negative,
        270 => b == Sign::Positive && z == Sign::Negative,
        _ => unreachable!(),
    };
    if outside {
        Location::Outside
    } else if z == Sign::Zero && (angle == 180 || b != Sign::Negative)
        || b == Sign::Zero
            && (angle != 180
                && z == if angle == 90 {
                    Sign::Positive
                } else {
                    Sign::Negative
                })
    {
        Location::Boundary
    } else {
        Location::Inside
    }
}

#[derive(Clone, Debug)]
struct FastEdge {
    start: [f64; 2],
    end: [f64; 2],
    circle: Option<([f64; 2], f64, bool)>,
    full: bool,
}
#[derive(Clone, Debug)]
struct FastProfile {
    matrix: [[f64; 3]; 3],
    translation: [f64; 3],
    edges: Vec<FastEdge>,
    sweep: [f64; 5],
    angle: u16,
    rotation: [f64; 2],
    through: bool,
    valid: bool,
}
impl FastProfile {
    fn new(profile: &Profile, sweep: &Sweep, norm: &Rational) -> Self {
        let inverse = profile.pose.inverse();
        let matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| inverse.rotation.matrix()[i][j].to_f64())
        });
        let translation = [
            &inverse.translation.x,
            &inverse.translation.y,
            &inverse.translation.z,
        ]
        .map(Real::to_f64);
        let edges = profile
            .loops
            .iter()
            .flatten()
            .map(|e| FastEdge {
                start: e.start.each_ref().map(Rational::to_f64_nearest),
                end: e.end.each_ref().map(Rational::to_f64_nearest),
                circle: e.circle.as_ref().map(|(c, r, t)| {
                    (
                        c.each_ref().map(Rational::to_f64_nearest),
                        r.to_f64_nearest(),
                        *t == catallax_recipe::Turn::Ccw,
                    )
                }),
                full: e.full,
            })
            .collect::<Vec<_>>();
        let (parameters, angle, rotation) = match sweep {
            Sweep::Extrude { vector, .. } => (
                [
                    vector.x.to_f64_nearest(),
                    vector.y.to_f64_nearest(),
                    vector.z.to_f64_nearest(),
                    0.0,
                    0.0,
                ],
                0,
                [1.0, 0.0],
            ),
            Sweep::Revolve {
                origin,
                axis,
                angle,
                ..
            } => (
                [
                    origin.x.to_f64_nearest(),
                    origin.y.to_f64_nearest(),
                    axis.x.to_f64_nearest(),
                    axis.y.to_f64_nearest(),
                    norm.to_f64_nearest(),
                ],
                angle.fast().0,
                angle.fast().1,
            ),
        };
        let mut result = Self {
            matrix,
            translation,
            edges,
            sweep: parameters,
            angle,
            rotation,
            through: matches!(sweep, Sweep::Extrude { through: true, .. }),
            valid: false,
        };
        result.valid = result.finite();
        result
    }
    fn finite(&self) -> bool {
        self.matrix
            .iter()
            .flatten()
            .chain(self.translation.iter())
            .chain(self.sweep.iter())
            .chain(self.rotation.iter())
            .all(|v| v.is_finite())
            && self.edges.iter().all(|e| {
                e.start.iter().chain(e.end.iter()).all(|v| v.is_finite())
                    && e.circle
                        .as_ref()
                        .is_none_or(|(c, r, _)| c.iter().all(|v| v.is_finite()) && r.is_finite())
            })
    }
    fn inside(&self, p: [f64; 3]) -> bool {
        if !self.valid {
            return false;
        }
        let p: [f64; 3] = std::array::from_fn(|i| {
            ((self.matrix[i][0] * p[0] + self.matrix[i][1] * p[1]) + self.matrix[i][2] * p[2])
                + self.translation[i]
        });
        let a = self.sweep;
        if self.angle == 0 {
            let t = p[2] / a[2];
            return (self.through || (0.0..=1.0).contains(&t))
                && self.region([p[0] - t * a[0], p[1] - t * a[1]]);
        }
        let v = [p[0] - a[0], p[1] - a[1]];
        let t = (v[0] * a[2] + v[1] * a[3]) / a[4];
        let b = (a[2] * v[1] - a[3] * v[0]) / a[4];
        let k = Scalar::sqrt(&(b * b + p[2] * p[2] / a[4]));
        for flip in [-1.0, 1.0] {
            let angular = match self.angle {
                90 => b * flip >= 0.0 && p[2] * flip >= 0.0,
                180 => p[2] * flip >= 0.0,
                270 => b * flip <= 0.0 || p[2] * flip >= 0.0,
                360 => true,
                sector => {
                    let start = p[2] * flip;
                    let end = (b * Scalar::sqrt(&a[4]) * self.rotation[1]
                        - p[2] * self.rotation[0])
                        * flip;
                    if sector == 1 {
                        start >= 0.0 && end >= 0.0
                    } else {
                        start >= 0.0 || end >= 0.0
                    }
                }
            };
            if angular
                && self.region([
                    a[0] + a[2] * t - a[3] * k * flip,
                    a[1] + a[3] * t + a[2] * k * flip,
                ])
            {
                return true;
            }
        }
        false
    }
    fn region(&self, p: [f64; 2]) -> bool {
        let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
        let sub = |a: [f64; 2], b: [f64; 2]| [a[0] - b[0], a[1] - b[1]];
        let mut crossings = 0;
        for e in &self.edges {
            if let Some((c, r, ccw)) = e.circle {
                let dy = p[1] - c[1];
                let discriminant = r - dy * dy;
                if discriminant <= 0.0 {
                    continue;
                }
                let root = Scalar::sqrt(&discriminant);
                for x in [c[0] - root, c[0] + root] {
                    if x <= p[0] {
                        continue;
                    }
                    let v = [x - c[0], dy];
                    let (start, end) = if ccw {
                        (e.start, e.end)
                    } else {
                        (e.end, e.start)
                    };
                    let a = sub(start, c);
                    let b = sub(end, c);
                    let av = cross(a, v);
                    let vb = cross(v, b);
                    let contains = e.full
                        || if cross(a, b) > 0.0 {
                            av >= 0.0 && vb >= 0.0
                        } else if cross(a, b) < 0.0 {
                            av >= 0.0 || vb >= 0.0
                        } else {
                            av >= 0.0
                        };
                    if contains {
                        // Half-open endpoints count a crossing once across a join.
                        let at_start = av == 0.0 && a[0] * v[0] + a[1] * v[1] > 0.0;
                        let at_end = vb == 0.0 && b[0] * v[0] + b[1] * v[1] > 0.0;
                        if e.full || (!at_start || a[0] > 0.0) && (!at_end || b[0] < 0.0) {
                            crossings += 1;
                        }
                    }
                }
            } else if (e.start[1] > p[1]) != (e.end[1] > p[1]) {
                let x = e.start[0]
                    + (p[1] - e.start[1]) * (e.end[0] - e.start[0]) / (e.end[1] - e.start[1]);
                if x > p[0] {
                    crossings += 1;
                }
            }
        }
        crossings % 2 == 1
    }
}

impl Feature {
    fn sample_profile_bounds(&self) -> [[Rational; 2]; 2] {
        let first = self
            .loops
            .iter()
            .flatten()
            .next()
            .expect("nonempty compiled feature");
        let mut bounds = first.start.each_ref().map(|x| [x.clone(), x.clone()]);
        for edge in self.loops.iter().flatten() {
            let limits = if let Some((center, radius, _)) = &edge.circle {
                let radius = Ball::from_rational(radius, 192).sqrt().unwrap().hi();
                center.each_ref().map(|x| [x - &radius, x + &radius])
            } else {
                std::array::from_fn(|i| {
                    [
                        edge.start[i].clone().min(edge.end[i].clone()),
                        edge.start[i].clone().max(edge.end[i].clone()),
                    ]
                })
            };
            for i in 0..2 {
                bounds[i][0] = bounds[i][0].clone().min(limits[i][0].clone());
                bounds[i][1] = bounds[i][1].clone().max(limits[i][1].clone());
            }
        }
        bounds
    }

    pub(crate) fn sample_bounds(&self) -> Result<Option<crate::samples::Bounds>, Refusal> {
        if self.algebraic_end.is_some() {
            return Ok(self
                .algebraic_display_box(&mut Decider::default())?
                .map(|b| b.map(|v| [v.lo(), v.hi()])));
        }
        let profile = self.sample_profile_bounds();
        let zero = Rational::zero();
        let local = match &self.sweep {
            Sweep::Extrude { vector, through } => {
                if *through {
                    return Ok(None);
                }
                let [x, y] = std::array::from_fn(|i| {
                    let v = [&vector.x, &vector.y][i];
                    [
                        &profile[i][0] + &zero.clone().min(v.clone()),
                        &profile[i][1] + &zero.clone().max(v.clone()),
                    ]
                });
                [
                    x,
                    y,
                    [
                        zero.clone().min(vector.z.clone()),
                        zero.max(vector.z.clone()),
                    ],
                ]
            }
            Sweep::Revolve { origin, .. } => {
                let radius = profile
                    .iter()
                    .map(|b| b[0].abs().max(b[1].abs()))
                    .fold(zero, |a, b| a + b)
                    + origin.x.abs()
                    + origin.y.abs()
                    + origin.z.abs();
                [&origin.x, &origin.y, &origin.z].map(|x| [x - &radius, x + &radius])
            }
        };
        Ok(Some(crate::samples::transform(
            &local,
            &self.inverse.inverse(),
        )?))
    }

    pub(crate) fn sample_candidates(
        &self,
        role: &Role,
        n: u32,
        result: &crate::samples::Bounds,
        visit: &mut impl FnMut(RPoint3) -> bool,
    ) -> Result<bool, Refusal> {
        let mut role = role;
        while let Role::Copy { of, .. } = role {
            role = &of.role;
        }
        let alias = self
            .fillet
            .as_ref()
            .filter(|f| &f.origin.role == role)
            .map(|_| Role::Side {
                segment: catallax_recipe::ElementId::new("blend").expect("fixed section id"),
                path: None,
            });
        let alias = alias.or_else(|| {
            self.chamfer
                .as_ref()
                .filter(|f| &f.origin.role == role)
                .map(|_| Role::Side {
                    segment: catallax_recipe::ElementId::new("bevel").expect("fixed section id"),
                    path: None,
                })
        });
        let role = alias.as_ref().unwrap_or(role);
        use crate::samples::{circle, fraction};
        if self.algebraic_end.is_some() {
            return Err(crate::samples::placement_error());
        }
        let pose = self.inverse.inverse();
        if pose
            .matrix()
            .iter()
            .flatten()
            .any(|v| v.rational_ref().is_none())
            || pose.translation.rational_components().is_none()
        {
            return Err(crate::samples::placement_error());
        }
        let mut emit = |p: RPoint3| {
            visit(
                pose.apply_rational_point(&p)
                    .expect("the native rational sampling frame was checked"),
            )
        };
        if matches!(role, Role::Start | Role::End) {
            let bounds = self.sample_profile_bounds();
            for i in 0..n {
                for j in 0..n {
                    let [x, y] = std::array::from_fn(|k| {
                        let f = fraction([i, j][k], n);
                        &bounds[k][0] + &((&bounds[k][1] - &bounds[k][0]) * f)
                    });
                    let p = RPoint3::new(x, y, Rational::zero());
                    let p = match &self.sweep {
                        Sweep::Extrude { vector, through } => {
                            if *through {
                                return Err(Refusal::new(
                                    "field.sample-role",
                                    "through-all tools have no cap faces",
                                ));
                            }
                            if matches!(role, Role::End) {
                                p.plus(vector)
                            } else {
                                p
                            }
                        }
                        Sweep::Revolve {
                            origin,
                            axis,
                            angle,
                        } => {
                            if matches!(role, Role::Start) {
                                p
                            } else {
                                if angle.is(360) {
                                    return Err(Refusal::new(
                                        "field.sample-role",
                                        "whole revolutions have no end cap",
                                    ));
                                }
                                let (cosine, sine) = angle
                                    .rational_rotation(&axis.norm_squared())
                                    .ok_or_else(crate::samples::placement_error)?;
                                sample_rotate(&p, origin, axis, &cosine, &sine)
                            }
                        }
                    };
                    if !emit(p) {
                        return Ok(false);
                    }
                }
            }
            return Ok(true);
        }
        let Role::Side { segment, .. } = role else {
            return Err(Refusal::new(
                "field.sample-role",
                "feature face role has no sampling definition",
            ));
        };
        let mut dec = Decider::default();
        for edge in self
            .loops
            .iter()
            .flatten()
            .filter(|e| e.ids.contains(segment))
        {
            let mut points = Vec::new();
            if let Some((center, radius, _)) = &edge.circle {
                let radial = profile::sqrt(radius)
                    .map(|r| [r, Rational::zero()])
                    .unwrap_or_else(|| profile::sub(&edge.start, center));
                for [x, y] in circle(n) {
                    let q = [
                        &center[0] + &(&(&radial[0] * &x) - &(&radial[1] * &y)),
                        &center[1] + &(&(&radial[0] * &y) + &(&radial[1] * &x)),
                    ];
                    if edge.contains(&profile::real(&q), &mut dec)? {
                        points.push(q);
                    }
                }
            } else {
                for i in 0..n {
                    let f = fraction(i, n);
                    points.push(std::array::from_fn(|k| {
                        &edge.start[k] + &((&edge.end[k] - &edge.start[k]) * f.clone())
                    }));
                }
            }
            match &self.sweep {
                Sweep::Extrude { vector, through } => {
                    let parameters = if *through {
                        let local = crate::samples::transform(result, &self.inverse)?;
                        let mut range = local[2].each_ref().map(|x| x / &vector.z);
                        range.sort();
                        range
                    } else {
                        [Rational::zero(), Rational::one()]
                    };
                    for q in points {
                        let p = point(&q);
                        for i in 0..n {
                            let t = &parameters[0]
                                + &((&parameters[1] - &parameters[0]) * fraction(i, n));
                            if !emit(p.plus(&vector.scaled(&t))) {
                                return Ok(false);
                            }
                        }
                    }
                }
                Sweep::Revolve { origin, axis, .. } => {
                    // Rational circle parametrisation works even when the axis
                    // length is irrational. Face queries apply the angular trim.
                    let norm = axis.norm_squared();
                    let scale = Ball::from_rational(&norm, 192)
                        .sqrt()
                        .unwrap()
                        .recip()
                        .unwrap()
                        .hi();
                    let mut rotations = Vec::with_capacity((4 * n) as usize);
                    for sign in [1, -1] {
                        for i in 0..2 * n {
                            let t = &(Rational::from_i64(2) * fraction(i, 2 * n) - Rational::one())
                                * &scale;
                            let squared = &norm * &(&t * &t);
                            let denominator = &Rational::one() + &squared;
                            let cosine = &(Rational::one() - squared) / &denominator
                                * Rational::from_i64(sign);
                            let sine =
                                Rational::from_i64(2) * t / denominator * Rational::from_i64(sign);
                            rotations.push((cosine, sine));
                        }
                    }
                    for q in points {
                        let p = point(&q);
                        let height = &p.minus(origin).dot(axis) / &norm;
                        let projection = origin.plus(&axis.scaled(&height));
                        let radial = p.minus(&projection);
                        let tangent = axis.cross(&radial);
                        for (cosine, sine) in &rotations {
                            if !emit(
                                projection
                                    .plus(&radial.scaled(cosine))
                                    .plus(&tangent.scaled(sine)),
                            ) {
                                return Ok(false);
                            }
                        }
                    }
                }
            }
        }
        Ok(true)
    }
}

fn sample_rotate(
    p: &RPoint3,
    origin: &RPoint3,
    axis: &RVector3,
    cosine: &Rational,
    sine_over_axis_length: &Rational,
) -> RPoint3 {
    let height = p.minus(origin).dot(axis) / axis.norm_squared();
    let projection = origin.plus(&axis.scaled(&height));
    let radial = p.minus(&projection);
    projection
        .plus(&radial.scaled(cosine))
        .plus(&axis.cross(&radial).scaled(sine_over_axis_length))
}

#[cfg(test)]
#[path = "features/distance_tests.rs"]
mod distance_tests;
mod planar_thicken;
