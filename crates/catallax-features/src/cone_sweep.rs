//! Exact rotational normal sweeps with their finite meridian domains.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Circle, Cone, Curve, Cylinder, Plane, Point3, Surface, Vector3};

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

/// A convex meridian polygon revolved about an exact axis.
/// Convexity is in radius and axial height, not in three-dimensional space.
#[derive(Clone, Debug)]
pub struct ConeSweep {
    /// Axis reference point.
    pub origin: Point3,
    /// Unit axial direction.
    pub axis: Vector3,
    /// Perpendicular reference direction.
    pub reference: Vector3,
    /// Counter-clockwise radius/height vertices.
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
            let dr = &q[0] - &p[0];
            let dh = &q[1] - &p[1];
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
                let slope = &dr / &dh;
                let increasing = sign(&slope, dec)? == Sign::Positive;
                let apex = origin.plus(&axis.scaled(&(&p[1] - &(&p[0] / &slope))));
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
        self.sweep_face(selected, thickness, side, step, false, dec)
    }

    fn sweep_face(
        &self,
        selected: usize,
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        tip_sector: bool,
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
                let has_tip = self.apex.as_ref().is_some_and(|a| a.1 == selected);
                if boundaries.len() + usize::from(has_tip) != 2 || (has_tip && !tip_sector) {
                    return Err(bad("a cone-tip sweep needs its spherical normal sector"));
                }
                let mut points = if has_tip {
                    vec![[Real::zero(), Real::zero()]]
                } else {
                    Vec::new()
                };
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
                let lower_rim = if has_tip { Vec::new() } else { rim(0)? };
                let upper_rim = rim(usize::from(!has_tip))?;
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
