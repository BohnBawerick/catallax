//! Convex coaxial plane, cylinder and sphere cells with complete circular trims.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, FaceRef, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Circle, Curve, Point3, Surface, Torus, Vector3};

use crate::{
    cone_skin::{ConeBoundary, ConeFace, reduced},
    profile::sign,
};

/// One finite meridian interval of a physical curved face.
#[derive(Clone, Debug)]
pub struct CircularBand {
    /// Index in the cell's physical faces.
    pub face: usize,
    /// Lower and upper axial coordinates relative to the cell origin.
    pub heights: [Real; 2],
}

/// Intersection of convex rotational supports, with its complete active boundary.
#[derive(Clone, Debug)]
pub struct CircularSkin {
    /// Physical faces, including disconnected bands of the same support.
    pub faces: Vec<ConeFace>,
    /// Complete circular sections, ordered from lower to upper face.
    pub boundaries: Vec<ConeBoundary>,
    /// Physical curved bands. Poles do not create artificial boundary circles.
    pub bands: Vec<CircularBand>,
    /// Axial coordinate origin.
    pub origin: Point3,
    /// Unit axis shared by every support.
    pub axis: Vector3,
    /// Perpendicular reference direction.
    pub reference: Vector3,
}

pub(crate) fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.cylinder-domain", detail)
}
fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Zero)
}
fn insert(values: &mut Vec<Real>, value: Real, dec: &mut Decider) -> Result<(), Refusal> {
    for i in 0..values.len() {
        match sign(&(&value - &values[i]), dec)? {
            Sign::Zero => return Ok(()),
            Sign::Negative => {
                values.insert(i, value);
                return Ok(());
            }
            Sign::Positive => (),
        }
    }
    values.push(value);
    Ok(())
}
fn radial_sq(surface: &Surface, height: &Real, origin: &Point3, axis: &Vector3) -> Real {
    match surface {
        Surface::Cylinder(c) => c.radius_sq.clone(),
        Surface::Sphere(s) => {
            &s.radius.square() - &(height - &s.center.minus(origin).dot(axis)).square()
        }
        _ => unreachable!("a meridian band is curved"),
    }
}
fn equivalent(a: &Surface, b: &Surface, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(match (a, b) {
        (Surface::Cylinder(a), Surface::Cylinder(b)) => zero(&(&a.radius_sq - &b.radius_sq), dec)?,
        (Surface::Sphere(a), Surface::Sphere(b)) => {
            zero(&a.center.minus(&b.center).norm_squared(), dec)?
                && zero(&(&a.radius - &b.radius), dec)?
        }
        (Surface::Plane(a), Surface::Plane(b)) => {
            sign(&a.normal.dot(&b.normal), dec)? == Sign::Positive
                && zero(&a.origin.minus(&b.origin).dot(&a.normal), dec)?
        }
        _ => false,
    })
}

impl CircularSkin {
    /// Recognises a bounded convex intersection. Source trim completeness must
    /// be checked separately by a B-rep caller, not inferred from these supports.
    pub fn new(
        surfaces: &[(Surface, Vec<Origin>)],
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let Some(cylinder) = surfaces.iter().find_map(|(s, _)| {
            if let Surface::Cylinder(c) = s {
                Some(c)
            } else {
                None
            }
        }) else {
            return Ok(None);
        };
        if surfaces.iter().any(|(s, _)| {
            !matches!(
                s,
                Surface::Cylinder(_) | Surface::Plane(_) | Surface::Sphere(_)
            )
        }) {
            return Ok(None);
        }
        let length = reduced(cylinder.axis.norm_squared().sqrt(), dec)?;
        let axis = cylinder.axis.scaled(&(&Real::one() / &length));
        let origin = cylinder.origin.clone();
        let reference = cylinder.ref_dir.clone();
        let mut supports: Vec<(Surface, Vec<Origin>)> = Vec::new();
        for (surface, names) in surfaces {
            let coaxial = match surface {
                Surface::Cylinder(c) => {
                    zero(&c.axis.cross(&axis).norm_squared(), dec)?
                        && zero(&c.origin.minus(&origin).cross(&axis).norm_squared(), dec)?
                }
                Surface::Sphere(s) => {
                    zero(&s.center.minus(&origin).cross(&axis).norm_squared(), dec)?
                }
                Surface::Plane(p) => zero(&p.normal.cross(&axis).norm_squared(), dec)?,
                _ => unreachable!(),
            };
            if !coaxial {
                return Err(bad(
                    "noncoaxial or oblique cylindrical joins need a finite edge-pipe trim arrangement",
                ));
            }
            let mut duplicate = false;
            for (previous, origins) in &mut supports {
                if equivalent(previous, surface, dec)? {
                    origins.extend(names.iter().cloned());
                    origins.sort();
                    origins.dedup();
                    duplicate = true;
                    break;
                }
            }
            if !duplicate {
                supports.push((surface.clone(), names.clone()));
            }
        }
        let mut lower: Option<Real> = None;
        let mut upper: Option<Real> = None;
        let mut curved = Vec::new();
        for (i, (surface, _)) in supports.iter().enumerate() {
            let limits = match surface {
                Surface::Sphere(s) => {
                    curved.push(i);
                    let center = s.center.minus(&origin).dot(&axis);
                    [Some(&center - &s.radius), Some(&center + &s.radius)]
                }
                Surface::Cylinder(_) => {
                    curved.push(i);
                    [None, None]
                }
                Surface::Plane(p) => {
                    let h = p.origin.minus(&origin).dot(&axis);
                    if sign(&p.normal.dot(&axis), dec)? == Sign::Negative {
                        [Some(h), None]
                    } else {
                        [None, Some(h)]
                    }
                }
                _ => unreachable!(),
            };
            if let Some(h) = limits[0].clone()
                && (lower.is_none()
                    || sign(&(&h - lower.as_ref().unwrap()), dec)? == Sign::Positive)
            {
                lower = Some(h);
            }
            if let Some(h) = limits[1].clone()
                && (upper.is_none()
                    || sign(&(&h - upper.as_ref().unwrap()), dec)? == Sign::Negative)
            {
                upper = Some(h);
            }
        }
        let (Some(lower), Some(upper)) = (lower, upper) else {
            return Err(bad("a cylindrical cell needs finite axial limits"));
        };
        if sign(&(&upper - &lower), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "distance.empty-cylinder",
                "the convex cylindrical supports consume the cell",
            ));
        }
        let mut heights = vec![lower.clone(), upper.clone()];
        for (k, i) in curved.iter().enumerate() {
            for j in &curved[k + 1..] {
                let a = &supports[*i].0;
                let b = &supports[*j].0;
                let roots = match (a, b) {
                    (Surface::Cylinder(c), Surface::Sphere(s))
                    | (Surface::Sphere(s), Surface::Cylinder(c)) => {
                        let q = &s.radius.square() - &c.radius_sq;
                        if sign(&q, dec)? == Sign::Negative {
                            Vec::new()
                        } else {
                            let dz = reduced(q.sqrt(), dec)?;
                            let z = s.center.minus(&origin).dot(&axis);
                            vec![&z - &dz, &z + &dz]
                        }
                    }
                    (Surface::Sphere(a), Surface::Sphere(b)) => {
                        let za = a.center.minus(&origin).dot(&axis);
                        let zb = b.center.minus(&origin).dot(&axis);
                        if zero(&(&za - &zb), dec)? {
                            Vec::new()
                        } else {
                            vec![reduced(
                                (&b.radius.square() - &a.radius.square() + za.square()
                                    - zb.square())
                                    / ((&za - &zb) * Real::from_i64(2)),
                                dec,
                            )?]
                        }
                    }
                    _ => Vec::new(),
                };
                for h in roots {
                    if sign(&(&h - &lower), dec)? == Sign::Positive
                        && sign(&(&upper - &h), dec)? == Sign::Positive
                    {
                        insert(&mut heights, h, dec)?;
                    }
                }
            }
        }
        let mut intervals: Vec<(usize, [Real; 2])> = Vec::new();
        for pair in heights.windows(2) {
            let middle = (&pair[0] + &pair[1]) / Real::from_i64(2);
            let mut active = curved[0];
            let mut radius = radial_sq(&supports[active].0, &middle, &origin, &axis);
            for i in &curved[1..] {
                let candidate = radial_sq(&supports[*i].0, &middle, &origin, &axis);
                if sign(&(&candidate - &radius), dec)? == Sign::Negative {
                    active = *i;
                    radius = candidate;
                }
            }
            if sign(&radius, dec)? != Sign::Positive {
                return Err(bad(
                    "a cylindrical meridian interval has no positive interior",
                ));
            }
            if let Some((previous, range)) = intervals.last_mut()
                && *previous == active
            {
                range[1] = pair[1].clone();
            } else {
                intervals.push((active, [pair[0].clone(), pair[1].clone()]));
            }
        }
        let mut result = Self {
            faces: Vec::new(),
            boundaries: Vec::new(),
            bands: Vec::new(),
            origin,
            axis,
            reference,
        };
        for (source, range) in intervals {
            let face = result.faces.len();
            result.faces.push(ConeFace {
                surface: supports[source].0.clone(),
                origins: supports[source].1.clone(),
            });
            result.bands.push(CircularBand {
                face,
                heights: range,
            });
        }
        for i in 0..result.bands.len() - 1 {
            let height = result.bands[i].heights[1].clone();
            let circle = result.circle(result.bands[i].face, &height, dec)?;
            result.boundaries.push(ConeBoundary {
                curve: Curve::Circle(circle),
                faces: [result.bands[i].face, result.bands[i + 1].face],
            });
        }
        for (end, height) in [(0, lower), (1, upper)] {
            let band = if end == 0 {
                &result.bands[0]
            } else {
                result.bands.last().unwrap()
            };
            let curved_face = band.face;
            let circle = result.circle(curved_face, &height, dec)?;
            if zero(&circle.radius_sq, dec)? {
                continue;
            }
            let mut source = None;
            for (surface, origins) in &supports {
                if let Surface::Plane(p) = surface
                    && zero(
                        &(&p.origin.minus(&result.origin).dot(&result.axis) - &height),
                        dec,
                    )?
                    && (sign(&p.normal.dot(&result.axis), dec)? == Sign::Negative) == (end == 0)
                {
                    source = Some((surface.clone(), origins.clone()));
                    break;
                }
            }
            let Some((surface, origins)) = source else {
                return Err(bad("a finite positive-radius end has no active plane"));
            };
            let face = result.faces.len();
            result.faces.push(ConeFace { surface, origins });
            result.boundaries.push(ConeBoundary {
                curve: Curve::Circle(circle),
                faces: if end == 0 {
                    [face, curved_face]
                } else {
                    [curved_face, face]
                },
            });
        }
        Ok(Some(result))
    }

    /// Exact complete section at an axial coordinate of a curved face.
    pub fn circle(&self, face: usize, height: &Real, dec: &mut Decider) -> Result<Circle, Refusal> {
        Ok(Circle {
            center: self.origin.plus(&self.axis.scaled(height)),
            normal: self.axis.clone(),
            radius_sq: reduced(
                radial_sq(&self.faces[face].surface, height, &self.origin, &self.axis),
                dec,
            )?,
            ref_dir: self.reference.clone(),
        })
    }

    /// Radial and axial components of an incident outward unit normal.
    pub fn normal(
        &self,
        face: usize,
        circle: &Circle,
        dec: &mut Decider,
    ) -> Result<[Real; 2], Refusal> {
        Ok(match &self.faces[face].surface {
            Surface::Plane(p) => [
                Real::zero(),
                Real::from_i64(if sign(&p.normal.dot(&self.axis), dec)? == Sign::Negative {
                    -1
                } else {
                    1
                }),
            ],
            Surface::Cylinder(_) => [Real::one(), Real::zero()],
            Surface::Sphere(s) => [
                reduced(&circle.radius_sq.sqrt() / &s.radius, dec)?,
                reduced(
                    &circle.center.minus(&s.center).dot(&self.axis) / &s.radius,
                    dec,
                )?,
            ],
            _ => unreachable!(),
        })
    }

    /// Complete edge-distance tube. Its exposed normal sector is selected by
    /// the offset boundary in B-rep and by regularised union in the field.
    pub fn edge_round(
        &self,
        boundary: &ConeBoundary,
        distance: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<ConeFace>, Refusal> {
        let Curve::Circle(circle) = &boundary.curve else {
            unreachable!()
        };
        let [a, b] = boundary.faces;
        let na = self.normal(a, circle, dec)?;
        let nb = self.normal(b, circle, dec)?;
        if zero(&(&na[0] - &nb[0]), dec)? && zero(&(&na[1] - &nb[1]), dec)? {
            return Ok(None);
        }
        if sign(&(&na[0] * &nb[1] - &na[1] * &nb[0]), dec)? != Sign::Positive {
            return Err(bad(
                "a nonconvex circular join needs its exposed envelope arrangement",
            ));
        }
        let radius = reduced(circle.radius_sq.sqrt(), dec)?;
        if sign(&(&radius - distance), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "offset.spindle-carrier",
                "the circular-edge distance tube reaches a horn or spindle torus",
            ));
        }
        let left = self.faces[a]
            .origins
            .first()
            .ok_or_else(|| bad("a circular edge needs named incident faces"))?;
        let right = self.faces[b]
            .origins
            .first()
            .ok_or_else(|| bad("a circular edge needs named incident faces"))?;
        let edge = EdgeName::new(left.clone(), right.clone()).map_err(|e| bad(&e.to_string()))?;
        Ok(Some(ConeFace {
            surface: Surface::Torus(Torus {
                center: circle.center.clone(),
                axis: self.axis.clone(),
                major: radius,
                minor: distance.clone(),
                ref_dir: circle.ref_dir.clone(),
            }),
            origins: vec![Origin {
                step: step.clone(),
                role: Role::Round { edge },
            }],
        }))
    }

    /// Parallel faces and convex circular normal sectors. A complete edge tube
    /// is admitted only when its existing torus carrier remains a ring torus.
    pub fn outward(
        &self,
        distance: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<(Vec<ConeFace>, Vec<ConeBoundary>), Refusal> {
        let mut faces =
            self.faces
                .iter()
                .map(|f| {
                    Ok(ConeFace {
                        surface: f.surface.parallel_support(distance, dec)?.ok_or_else(|| {
                            bad("a positive offset unexpectedly collapsed a support")
                        })?,
                        origins: f.origins.clone(),
                    })
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
        let mut boundaries = Vec::new();
        for boundary in &self.boundaries {
            let Curve::Circle(circle) = &boundary.curve else {
                unreachable!()
            };
            let radius = reduced(circle.radius_sq.sqrt(), dec)?;
            let [a, b] = boundary.faces;
            let na = self.normal(a, circle, dec)?;
            let nb = self.normal(b, circle, dec)?;
            let shifted = |normal: &[Real; 2]| Circle {
                center: circle
                    .center
                    .plus(&self.axis.scaled(&(distance * &normal[1]))),
                normal: self.axis.clone(),
                radius_sq: (&radius + &(distance * &normal[0])).square(),
                ref_dir: circle.ref_dir.clone(),
            };
            let Some(tube) = self.edge_round(boundary, distance, step, dec)? else {
                boundaries.push(ConeBoundary {
                    curve: Curve::Circle(shifted(&na)),
                    faces: [a, b],
                });
                continue;
            };
            let round = faces.len();
            faces.push(tube);
            boundaries.push(ConeBoundary {
                curve: Curve::Circle(shifted(&na)),
                faces: [a, round],
            });
            boundaries.push(ConeBoundary {
                curve: Curve::Circle(shifted(&nb)),
                faces: [round, b],
            });
        }
        Ok((faces, boundaries))
    }

    /// Closed-shell cavity from the intersection of inward convex parallels.
    /// Open spherical trims need a separate kept-face edge-distance proof.
    pub fn cavity(
        &self,
        open: &[bool],
        thickness: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if open.iter().any(|v| *v) {
            return Err(bad(
                "open curved circular cells need their kept-face edge-distance envelope",
            ));
        }
        let mut surfaces = Vec::new();
        for face in &self.faces {
            let Some(surface) = face.surface.parallel_support(&(-thickness), dec)? else {
                return Ok(None);
            };
            let names = face
                .origins
                .iter()
                .map(|of| Origin {
                    step: step.clone(),
                    role: Role::Inner {
                        of: Box::new(of.clone()),
                    },
                })
                .collect();
            surfaces.push((surface, names));
        }
        match Self::new(&surfaces, dec) {
            Err(e) if e.code == "distance.empty-cylinder" => Ok(None),
            other => other,
        }
    }

    /// Resolves complete physical bands and their circular neighbours.
    pub fn select(&self, references: &[FaceRef]) -> Result<Vec<bool>, Refusal> {
        let mut selected = vec![false; self.faces.len()];
        for reference in references {
            let mut candidates = Vec::new();
            for (i, face) in self.faces.iter().enumerate() {
                if face.origins.iter().any(|o| reference.face.matches(o))
                    && reference.touching.iter().all(|pattern| {
                        self.boundaries
                            .iter()
                            .filter(|b| b.faces.contains(&i))
                            .any(|b| {
                                b.faces.iter().any(|j| {
                                    *j != i
                                        && self.faces[*j].origins.iter().any(|o| pattern.matches(o))
                                })
                            })
                    })
                {
                    candidates.push(i);
                }
            }
            if candidates.is_empty() {
                return Err(Refusal::new(
                    "reference.lost-face",
                    format!("no circular face resolves {reference:?}"),
                ));
            }
            if reference.near.is_some() && candidates.len() > 1 {
                return Err(bad(
                    "nearest circular-band selection needs a certified trimmed-face distance comparison",
                ));
            }
            for i in candidates {
                selected[i] = true;
            }
        }
        Ok(selected)
    }

    /// Normal image of a single smooth cylindrical or spherical band. The
    /// spherical rails follow constant polar angle, not constant axial height.
    pub fn face_sweep(
        &self,
        face: usize,
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<crate::cone_sweep::ConeSweep, Refusal> {
        let band = self.bands.iter().find(|b| b.face == face).ok_or_else(|| {
            bad("a circular planar disk sweep needs its separate cap construction")
        })?;
        let source = &self.faces[face];
        let front = matches!(side, ThickenSide::Front);
        let delta = if front { thickness.clone() } else { -thickness };
        let far = source
            .surface
            .parallel_support(&delta, dec)?
            .ok_or_else(|| {
                bad("a folded finite radial band needs its exposed normal-image arrangement")
            })?;
        let mut original = Vec::new();
        let mut shifted = Vec::new();
        let mut rims = Vec::new();
        for height in &band.heights {
            let circle = self.circle(face, height, dec)?;
            let radius = reduced(circle.radius_sq.sqrt(), dec)?;
            let normal = self.normal(face, &circle, dec)?;
            original.push([radius.clone(), height.clone()]);
            shifted.push([
                &radius + &(&delta * &normal[0]),
                height + &(&delta * &normal[1]),
            ]);
            if zero(&radius, dec)? {
                rims.push(Vec::new());
                continue;
            }
            let mut boundary = None;
            for b in &self.boundaries {
                if b.faces.contains(&face)
                    && let Curve::Circle(c) = &b.curve
                    && zero(&c.center.minus(&circle.center).norm_squared(), dec)?
                {
                    boundary = Some(b);
                    break;
                }
            }
            let boundary = boundary.ok_or_else(|| bad("a band end has no physical neighbour"))?;
            let other = boundary.faces.iter().copied().find(|i| *i != face).unwrap();
            let a = source
                .origins
                .first()
                .ok_or_else(|| bad("a radial band needs named source faces"))?;
            let b = self.faces[other]
                .origins
                .first()
                .ok_or_else(|| bad("a radial band needs named rim neighbours"))?;
            let edge = EdgeName::new(a.clone(), b.clone()).map_err(|e| bad(&e.to_string()))?;
            rims.push(vec![Origin {
                step: step.clone(),
                role: Role::Rim { edge },
            }]);
        }
        let far_names = source
            .origins
            .iter()
            .map(|of| Origin {
                step: step.clone(),
                role: Role::Offset {
                    of: Box::new(of.clone()),
                },
            })
            .collect();
        let (vertices, names, source_ends, far_ends) = if front {
            (
                vec![
                    original[0].clone(),
                    shifted[0].clone(),
                    shifted[1].clone(),
                    original[1].clone(),
                ],
                vec![
                    rims[0].clone(),
                    far_names,
                    rims[1].clone(),
                    source.origins.clone(),
                ],
                [3, 0],
                [1, 2],
            )
        } else {
            (
                vec![
                    shifted[0].clone(),
                    original[0].clone(),
                    original[1].clone(),
                    shifted[1].clone(),
                ],
                vec![
                    rims[0].clone(),
                    source.origins.clone(),
                    rims[1].clone(),
                    far_names,
                ],
                [1, 2],
                [3, 0],
            )
        };
        let mut sweep = crate::cone_sweep::ConeSweep::new(
            self.origin.clone(),
            self.axis.clone(),
            self.reference.clone(),
            vertices,
            names,
            dec,
        )?;
        for patch in &mut sweep.faces {
            if patch.ends == source_ends {
                patch.surface = source.surface.clone();
            }
            if patch.ends == far_ends {
                patch.surface = far.clone();
            }
        }
        Ok(sweep)
    }
}
