//! Exact rotational normal sweeps with their finite meridian domains.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId, ThickenSide, VertexName};
use catallax_shapes::{Circle, Cone, Curve, Cylinder, Plane, Point3, Sphere, Surface, Vector3};

use crate::{
    cone_skin::{ConeSkin, reduced},
    profile::sign,
};

/// An oriented analytic patch of a closed rotational sweep.
#[derive(Clone, Debug)]
pub struct SweepFace {
    /// Native rotational support.
    pub surface: Surface,
    /// True for an inner radial support.
    pub reversed: bool,
    /// Exact meridian endpoint indices.
    pub ends: [usize; 2],
    /// Physical source, offset or rim origins.
    pub origins: Vec<Origin>,
}

/// A closed rotational sweep with exact meridian endpoints and oriented patches.
/// The polygon constructor uses straight meridians; circular cells can retain
/// spherical patches between the same exact endpoints.
#[derive(Clone, Debug)]
pub struct ConeSweep {
    /// Axis reference point.
    pub origin: Point3,
    /// Unit axial direction.
    pub axis: Vector3,
    /// Perpendicular reference direction.
    pub reference: Vector3,
    /// Counter-clockwise radius/height endpoints of the physical patches.
    pub vertices: Vec<[Real; 2]>,
    /// Physical patches. Segments on the rotation axis are not faces.
    pub faces: Vec<SweepFace>,
}
fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.cone-domain", detail)
}
fn named(step: &StepId, source: &[Origin]) -> Vec<Origin> {
    source
        .iter()
        .map(|of| Origin {
            step: step.clone(),
            role: Role::Offset {
                of: Box::new(of.clone()),
            },
        })
        .collect()
}
impl ConeSweep {
    /// Builds a proven convex, positive-radius meridian sweep. Axis segments are
    /// permitted for planar disk sweeps; radial folds are not silently clipped.
    pub fn new(
        origin: Point3,
        axis: Vector3,
        reference: Vector3,
        vertices: Vec<[Real; 2]>,
        names: Vec<Vec<Origin>>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if vertices.len() < 3 || vertices.len() != names.len() {
            return Err(bad(
                "a rotational sweep needs a closed named meridian polygon",
            ));
        }
        let mut faces = Vec::new();
        for (i, p) in vertices.iter().enumerate() {
            let j = (i + 1) % vertices.len();
            let q = &vertices[j];
            if sign(&p[0], dec)? == Sign::Negative {
                return Err(bad(
                    "a folded conical sweep needs its exposed normal-image arrangement",
                ));
            }
            let dr = reduced(&q[0] - &p[0], dec)?;
            let dh = reduced(&q[1] - &p[1], dec)?;
            for v in &vertices {
                let cross = &(&dr * &(&v[1] - &p[1])) - &(&dh * &(&v[0] - &p[0]));
                if sign(&cross, dec)? == Sign::Negative {
                    return Err(bad(
                        "the rotational sweep meridian must be convex and counter-clockwise",
                    ));
                }
            }
            if sign(&p[0], dec)? == Sign::Zero && sign(&q[0], dec)? == Sign::Zero {
                continue;
            }
            let height_sign = sign(&dh, dec)?;
            let (surface, reversed) = if height_sign == Sign::Zero {
                let orientation = match sign(&dr, dec)? {
                    Sign::Positive => -1,
                    Sign::Negative => 1,
                    Sign::Zero => return Err(bad("a conical sweep has a collapsed meridian edge")),
                };
                (
                    Surface::Plane(Plane {
                        origin: origin.plus(&axis.scaled(&p[1])),
                        normal: axis.scaled(&Real::from_i64(orientation)),
                        ref_dir: reference.clone(),
                    }),
                    false,
                )
            } else if sign(&dr, dec)? == Sign::Zero {
                (
                    Surface::Cylinder(Cylinder {
                        origin: origin.clone(),
                        axis: axis.clone(),
                        radius_sq: p[0].square(),
                        ref_dir: reference.clone(),
                    }),
                    height_sign == Sign::Negative,
                )
            } else {
                let slope = reduced(&dr / &dh, dec)?;
                let increasing = sign(&slope, dec)? == Sign::Positive;
                let shift = reduced(&p[1] - &(&p[0] / &slope), dec)?;
                let apex = origin.plus(&axis.scaled(&shift));
                (
                    Surface::Cone(Cone {
                        apex,
                        axis: axis.scaled(&Real::from_i64(if increasing { 1 } else { -1 })),
                        slope: if increasing { slope } else { -slope },
                        ref_dir: reference.clone(),
                    }),
                    height_sign == Sign::Negative,
                )
            };
            faces.push(SweepFace {
                surface,
                reversed,
                ends: [i, j],
                origins: names[i].clone(),
            });
        }
        Ok(Self {
            origin,
            axis,
            reference,
            vertices,
            faces,
        })
    }

    /// Exact circle at one meridian vertex. A zero radius is not a ring edge.
    pub fn circle(&self, i: usize) -> Circle {
        Circle {
            center: self.origin.plus(&self.axis.scaled(&self.vertices[i][1])),
            normal: self.axis.clone(),
            radius_sq: self.vertices[i][0].square(),
            ref_dir: self.reference.clone(),
        }
    }
}
impl ConeSkin {
    /// Sweeps one complete conical band or perpendicular circular cap along its
    /// own normal. Sharp multi-face joins and cone-tip normal sectors refuse.
    pub fn face_sweep(
        &self,
        selected: usize,
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<ConeSweep, Refusal> {
        let f = &self.faces[selected];
        let front = matches!(side, ThickenSide::Front);
        let delta = if front { thickness.clone() } else { -thickness };
        let boundaries: Vec<_> = self
            .boundaries
            .iter()
            .filter(|b| b.faces.contains(&selected))
            .collect();
        let rim = |i: usize| -> Result<Vec<Origin>, Refusal> {
            let adjacent = &boundaries[i].faces;
            let a = self.faces[adjacent[0]]
                .origins
                .first()
                .ok_or_else(|| bad("a conical rim needs named incident faces"))?;
            let b = self.faces[adjacent[1]]
                .origins
                .first()
                .ok_or_else(|| bad("a conical rim needs named incident faces"))?;
            Ok(vec![Origin {
                step: step.clone(),
                role: Role::Rim {
                    edge: EdgeName::new(a.clone(), b.clone()).map_err(|e| bad(&e.to_string()))?,
                },
            }])
        };
        match &f.surface {
            Surface::Cone(c) => {
                if boundaries.len() != 2 || self.apex.as_ref().is_some_and(|a| a.1 == selected) {
                    return Err(bad("a cone-tip sweep needs its spherical normal sector"));
                }
                let mut points = Vec::new();
                for b in &boundaries {
                    let Curve::Circle(circle) = &b.curve else {
                        return Err(bad(
                            "oblique conical normal-sweep rims need their analytic envelope trims",
                        ));
                    };
                    if sign(&circle.normal.cross(&self.axis).norm_squared(), dec)? != Sign::Zero {
                        return Err(bad("a conical sweep requires perpendicular circular trims"));
                    }
                    points.push([
                        reduced(circle.radius_sq.sqrt(), dec)?,
                        reduced(circle.center.minus(&c.apex).dot(&self.axis), dec)?,
                    ]);
                }
                if sign(&(&points[1][1] - &points[0][1]), dec)? != Sign::Positive {
                    return Err(bad("conical sweep sections are not ordered"));
                }
                let factor = reduced((&Real::one() + &c.slope.square()).sqrt(), dec)?;
                let shift = [&delta / &factor, -(&(&delta * &c.slope) / &factor)];
                let far = points
                    .iter()
                    .map(|p| [&p[0] + &shift[0], &p[1] + &shift[1]])
                    .collect::<Vec<_>>();
                for p in &far {
                    if sign(&p[0], dec)? != Sign::Positive {
                        return Err(bad("a conical normal sweep reaches its radial fold"));
                    }
                }
                let lower_rim = rim(0)?;
                let upper_rim = rim(1)?;
                let (vertices, names) = if front {
                    (
                        vec![
                            points[0].clone(),
                            far[0].clone(),
                            far[1].clone(),
                            points[1].clone(),
                        ],
                        vec![
                            lower_rim,
                            named(step, &f.origins),
                            upper_rim,
                            f.origins.clone(),
                        ],
                    )
                } else {
                    (
                        vec![
                            far[0].clone(),
                            points[0].clone(),
                            points[1].clone(),
                            far[1].clone(),
                        ],
                        vec![
                            lower_rim,
                            f.origins.clone(),
                            upper_rim,
                            named(step, &f.origins),
                        ],
                    )
                };
                ConeSweep::new(
                    c.apex.clone(),
                    self.axis.clone(),
                    c.ref_dir.clone(),
                    vertices,
                    names,
                    dec,
                )
            }
            Surface::Plane(p) => {
                let [boundary] = boundaries.as_slice() else {
                    return Err(bad("a conical cap sweep needs one complete disk"));
                };
                let Curve::Circle(circle) = &boundary.curve else {
                    return Err(bad(
                        "an elliptic cap sweep needs a finite elliptic-cylinder definition",
                    ));
                };
                if sign(&p.normal.cross(&self.axis).norm_squared(), dec)? != Sign::Zero {
                    return Err(bad("a conical cap sweep needs a perpendicular disk"));
                }
                let orientation = if sign(&p.normal.dot(&self.axis), dec)? == Sign::Positive {
                    1
                } else {
                    -1
                };
                let h = &delta * &Real::from_i64(orientation);
                let radius = reduced(circle.radius_sq.sqrt(), dec)?;
                let (low, high, low_names, high_names) = if sign(&h, dec)? == Sign::Positive {
                    (Real::zero(), h, f.origins.clone(), named(step, &f.origins))
                } else {
                    (h, Real::zero(), named(step, &f.origins), f.origins.clone())
                };
                ConeSweep::new(
                    circle.center.clone(),
                    self.axis.clone(),
                    circle.ref_dir.clone(),
                    vec![
                        [Real::zero(), low.clone()],
                        [radius.clone(), low],
                        [radius, high.clone()],
                        [Real::zero(), high],
                    ],
                    vec![low_names, rim(0)?, high_names, Vec::new()],
                    dec,
                )
            }
            _ => Err(bad(
                "a selected cylindrical neighbour uses the complete cylinder sweep definition",
            )),
        }
    }
}

/// One exact nearest-point stratum of an outward conical distance body.
#[derive(Clone, Debug)]
pub enum ConeDistancePiece {
    /// Finite positive normal image of an input face.
    Sweep(ConeSweep),
    /// Complete circular edge-distance tube, with exposed sectors selected by union.
    Round(Surface, Vec<Origin>),
}

impl ConeSkin {
    /// Decomposes outward conical distance by nearest face, circular edge and tip.
    /// The field unions these finite normal images independently of the B-rep's
    /// directly assembled boundary sectors.
    pub fn outward_pieces(
        &self,
        distance: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<ConeDistancePiece>, Refusal> {
        use catallax_shapes::Torus;
        if self
            .faces
            .iter()
            .any(|f| matches!(f.surface, Surface::Cylinder(_)))
        {
            return Err(bad(
                "outward cylindrical neighbours need their circular-edge join arrangement",
            ));
        }
        let mut pieces = Vec::new();
        for (i, f) in self.faces.iter().enumerate() {
            let mut sweep = if let Some((apex, incident)) = &self.apex
                && *incident == i
            {
                let Surface::Cone(cone) = &f.surface else {
                    unreachable!("a conical tip is incident to its conical wall")
                };
                let [boundary] = self.boundaries.as_slice() else {
                    return Err(bad("a pointed conical wall needs one complete cap"));
                };
                let Curve::Circle(circle) = &boundary.curve else {
                    return Err(bad("an oblique tip normal image needs its elliptic rim"));
                };
                if sign(&circle.normal.cross(&self.axis).norm_squared(), dec)? != Sign::Zero {
                    return Err(bad("a tip normal image needs a perpendicular circular cap"));
                }
                let radius = reduced(circle.radius_sq.sqrt(), dec)?;
                let height = reduced(circle.center.minus(apex).dot(&self.axis), dec)?;
                let factor_sq = &Real::one() + &cone.slope.square();
                let factor = reduced(factor_sq.sqrt(), dec)?;
                let radial = &(distance / &factor_sq) * &factor;
                let axial = -(&radial * &cone.slope);
                ConeSweep::new(
                    apex.clone(),
                    self.axis.clone(),
                    cone.ref_dir.clone(),
                    vec![
                        [Real::zero(), Real::zero()],
                        [radial.clone(), axial.clone()],
                        [&radius + &radial, &height + &axial],
                        [radius, height],
                    ],
                    vec![Vec::new(), f.origins.clone(), Vec::new(), Vec::new()],
                    dec,
                )?
            } else {
                self.face_sweep(i, distance, ThickenSide::Front, step, dec)?
            };
            for patch in &mut sweep.faces {
                if patch
                    .origins
                    .iter()
                    .any(|o| matches!(o.role, Role::Offset { .. }))
                {
                    patch.origins = f.origins.clone();
                } else if patch
                    .origins
                    .iter()
                    .any(|o| matches!(o.role, Role::Rim { .. }))
                {
                    patch.origins.clear();
                }
            }
            pieces.push(ConeDistancePiece::Sweep(sweep));
        }
        if let Some((apex, incident)) = &self.apex {
            let Surface::Cone(cone) = &self.faces[*incident].surface else {
                unreachable!("a conical tip is incident to its conical wall")
            };
            let origin = self.faces[*incident]
                .origins
                .first()
                .ok_or_else(|| bad("a conical tip needs a named incident face"))?;
            let vertex = VertexName::new(vec![origin.clone()]).map_err(|e| bad(&e.to_string()))?;
            pieces.push(ConeDistancePiece::Round(
                Surface::Sphere(Sphere {
                    center: apex.clone(),
                    radius: distance.clone(),
                    axis: self.axis.clone(),
                    ref_dir: cone.ref_dir.clone(),
                }),
                vec![Origin {
                    step: step.clone(),
                    role: Role::RoundCorner { vertex },
                }],
            ));
        }
        for boundary in &self.boundaries {
            let Curve::Circle(circle) = &boundary.curve else {
                return Err(bad(
                    "elliptic conical edges need a certified non-toroidal distance carrier",
                ));
            };
            if sign(&circle.normal.cross(&self.axis).norm_squared(), dec)? != Sign::Zero {
                return Err(bad(
                    "an outward conical edge needs a perpendicular circular section",
                ));
            }
            let radius = reduced(circle.radius_sq.sqrt(), dec)?;
            if sign(&(&radius - distance), dec)? != Sign::Positive {
                return Err(Refusal::new(
                    "offset.spindle-carrier",
                    "the circular-edge distance tube reaches a horn or spindle torus",
                ));
            }
            let a = self.faces[boundary.faces[0]]
                .origins
                .first()
                .ok_or_else(|| bad("a conical edge needs named incident faces"))?;
            let b = self.faces[boundary.faces[1]]
                .origins
                .first()
                .ok_or_else(|| bad("a conical edge needs named incident faces"))?;
            let edge = EdgeName::new(a.clone(), b.clone()).map_err(|e| bad(&e.to_string()))?;
            pieces.push(ConeDistancePiece::Round(
                Surface::Torus(Torus {
                    center: circle.center.clone(),
                    axis: self.axis.clone(),
                    major: radius,
                    minor: distance.clone(),
                    ref_dir: circle.ref_dir.clone(),
                }),
                vec![Origin {
                    step: step.clone(),
                    role: Role::Round { edge },
                }],
            ));
        }
        Ok(pieces)
    }
}
