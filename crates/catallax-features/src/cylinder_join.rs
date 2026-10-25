//! Bounded complete circular joins, with explicit finite source patches.
//!
//! Orthogonal meridians have exact parallel planes/cylinders. A convex turn
//! sweeps a ring-torus sector; a concave turn intersects the two parallels.
//! Clearances below keep these local strata disjoint. No fitted carrier is used.

use crate::{
    cone_skin::{ConeBoundary, ConeFace, reduced},
    cylinder_skin::CylinderSkin,
    polyhedron::ConvexCell,
    profile::sign,
};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, FaceRef, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Circle, Curve, Surface, Torus};

/// Complete source bands of either a two-radius shaft or a contained plate boss.
#[derive(Clone, Debug)]
pub struct CylinderJoin {
    /// Physical patches, in meridian order (plate top, wall, cap for a boss).
    pub faces: Vec<ConeFace>,
    /// Complete source rings, ordered along the material meridian.
    pub boundaries: Vec<ConeBoundary>,
    /// Complete cylindrical source bands and their wall patch indices.
    pub bands: Vec<(usize, CylinderSkin)>,
    /// A rectangular base and the index of its attachment plane.
    pub plate: Option<(ConvexCell, usize)>,
}
/// Stable refusal for domains outside the bounded analytic join proof.
pub fn unsupported(detail: &str) -> Refusal {
    Refusal::new("distance.unsupported-envelope", detail)
}
fn zero(x: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(x, dec)? == Sign::Zero)
}
fn same_axis(a: &CylinderSkin, b: &CylinderSkin, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        zero(&a.cylinder.axis.cross(&b.cylinder.axis).norm_squared(), dec)?
            && sign(&a.cylinder.axis.dot(&b.cylinder.axis), dec)? == Sign::Positive
            && zero(
                &a.ends[0]
                    .minus(&b.ends[0])
                    .cross(&a.cylinder.axis)
                    .norm_squared(),
                dec,
            )?,
    )
}
impl CylinderJoin {
    /// Checks two axially consecutive complete cylinders. Their shared cap is
    /// removed; the exposed larger cap is a finite annulus, never a full disk.
    pub fn stepped(
        mut a: CylinderSkin,
        mut b: CylinderSkin,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if !same_axis(&a, &b, dec)? {
            return Err(unsupported("cylinder joins require proven coaxial bands"));
        }
        if sign(&b.ends[0].minus(&a.ends[0]).dot(&a.cylinder.axis), dec)? == Sign::Negative {
            std::mem::swap(&mut a, &mut b);
        }
        if !zero(&a.ends[1].minus(&b.ends[0]).norm_squared(), dec)? {
            return Err(unsupported(
                "shaft bands must meet at one complete perpendicular section",
            ));
        }
        let order = sign(&(&a.cylinder.radius_sq - &b.cylinder.radius_sq), dec)?;
        if order == Sign::Zero {
            return Err(unsupported(
                "equal-radius bands use the existing complete-cylinder path",
            ));
        }
        let shoulder = if order == Sign::Positive {
            a.surfaces()[1].clone()
        } else {
            b.surfaces()[0].clone()
        };
        let shoulder_names = if order == Sign::Positive {
            a.origins[1].clone()
        } else {
            b.origins[0].clone()
        };
        a.origins[1] = shoulder_names.clone();
        b.origins[0] = shoulder_names.clone();
        let faces = vec![
            ConeFace {
                surface: a.surfaces()[0].clone(),
                origins: a.origins[0].clone(),
            },
            ConeFace {
                surface: Surface::Cylinder(a.cylinder.clone()),
                origins: a.origins[2].clone(),
            },
            ConeFace {
                surface: shoulder,
                origins: shoulder_names,
            },
            ConeFace {
                surface: Surface::Cylinder(b.cylinder.clone()),
                origins: b.origins[2].clone(),
            },
            ConeFace {
                surface: b.surfaces()[1].clone(),
                origins: b.origins[1].clone(),
            },
        ];
        let boundaries = vec![
            ring(&a, 0, [0, 1]),
            ring(&a, 1, [1, 2]),
            ring(&b, 0, [2, 3]),
            ring(&b, 1, [3, 4]),
        ];
        Ok(Self {
            faces,
            boundaries,
            bands: vec![(1, a), (3, b)],
            plate: None,
        })
    }
    /// Proves a perpendicular boss attached wholly inside one rectangular face.
    pub fn boss(
        cell: ConvexCell,
        mut cylinder: CylinderSkin,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if cell.faces.len() != 6
            || cell.vertices().len() != 8
            || cell.polygons.iter().any(|p| p.len() != 4)
        {
            return Err(unsupported(
                "a plate boss requires one complete rectangular base",
            ));
        }
        for (i, a) in cell.faces.iter().enumerate() {
            for b in &cell.faces[i + 1..] {
                if !zero(&a.plane.normal.cross(&b.plane.normal).norm_squared(), dec)?
                    && !zero(&a.plane.normal.dot(&b.plane.normal), dec)?
                {
                    return Err(unsupported(
                        "plate normals must be exactly orthogonal or parallel",
                    ));
                }
            }
        }
        let mut attached = None;
        for (i, f) in cell.faces.iter().enumerate() {
            if zero(
                &f.plane.normal.cross(&cylinder.cylinder.axis).norm_squared(),
                dec,
            )? && sign(&f.plane.normal.dot(&cylinder.cylinder.axis), dec)? == Sign::Positive
                && zero(
                    &cylinder.ends[0].minus(&f.plane.origin).dot(&f.plane.normal),
                    dec,
                )?
            {
                attached = Some(i);
            }
        }
        let top = attached
            .ok_or_else(|| unsupported("boss must start on an outward perpendicular plate face"))?;
        cylinder.origins[0] = cell.faces[top].origins.clone();
        let result = Self {
            faces: vec![
                ConeFace {
                    surface: Surface::Plane(cell.faces[top].plane.clone()),
                    origins: cell.faces[top].origins.clone(),
                },
                ConeFace {
                    surface: Surface::Cylinder(cylinder.cylinder.clone()),
                    origins: cylinder.origins[2].clone(),
                },
                ConeFace {
                    surface: cylinder.surfaces()[1].clone(),
                    origins: cylinder.origins[1].clone(),
                },
            ],
            boundaries: vec![ring(&cylinder, 0, [0, 1]), ring(&cylinder, 1, [1, 2])],
            bands: vec![(1, cylinder)],
            plate: Some((cell, top)),
        };
        result.clearance(&Real::zero(), dec)?;
        Ok(result)
    }
    /// Exact bounds preventing a new contact between disjoint source strata.
    pub fn clearance(&self, d: &Real, dec: &mut Decider) -> Result<(), Refusal> {
        for (_, band) in &self.bands {
            let radius = reduced(band.cylinder.radius(), dec)?;
            let height = band.ends[1].minus(&band.ends[0]).norm_squared();
            if sign(&(&radius - d), dec)? != Sign::Positive
                || sign(&(height - (d * &Real::from_i64(2)).square()), dec)? != Sign::Positive
            {
                return Err(unsupported(
                    "offset reaches a cylinder radius or a finite band clearance",
                ));
            }
        }
        if let Some((cell, _)) = &self.plate {
            let band = &self.bands[0].1;
            let radius = reduced(band.cylinder.radius(), dec)? + d.clone();
            for (i, a) in cell.faces.iter().enumerate() {
                for b in &cell.faces[i + 1..] {
                    if zero(&a.plane.normal.cross(&b.plane.normal).norm_squared(), dec)?
                        && sign(&a.plane.normal.dot(&b.plane.normal), dec)? == Sign::Negative
                    {
                        let gap = a.plane.origin.minus(&b.plane.origin).dot(&a.plane.normal);
                        if sign(
                            &(gap.square()
                                - (d * &Real::from_i64(2)).square()
                                    * a.plane.normal.norm_squared()),
                            dec,
                        )? != Sign::Positive
                        {
                            return Err(unsupported(
                                "offset reaches a finite rectangular plate clearance",
                            ));
                        }
                    }
                }
            }
            for f in &cell.faces {
                if !zero(&f.plane.normal.dot(&band.cylinder.axis), dec)? {
                    continue;
                }
                let gap = -band.ends[0].minus(&f.plane.origin).dot(&f.plane.normal);
                if sign(&gap, dec)? != Sign::Positive
                    || sign(
                        &(gap.square() - radius.square() * f.plane.normal.norm_squared()),
                        dec,
                    )? != Sign::Positive
                {
                    return Err(unsupported(
                        "the expanded boss ring must remain strictly inside the rectangular plate",
                    ));
                }
            }
        } else {
            let a = reduced(self.bands[0].1.cylinder.radius(), dec)?;
            let b = reduced(self.bands[1].1.cylinder.radius(), dec)?;
            if sign(&((a - b).square() - (d * &Real::from_i64(2)).square()), dec)? != Sign::Positive
            {
                return Err(unsupported(
                    "offset reaches the finite annular shoulder clearance",
                ));
            }
        }
        Ok(())
    }
    /// Exact displaced faces and convex torus sectors, with shared contact rings.
    pub fn outward(
        &self,
        d: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<(Vec<ConeFace>, Vec<ConeBoundary>), Refusal> {
        self.clearance(d, dec)?;
        let band = &self.bands[0].1;
        let axis = band
            .cylinder
            .axis
            .scaled(&(&Real::one() / &reduced(band.cylinder.axis.norm_squared().sqrt(), dec)?));
        let mut faces = self
            .faces
            .iter()
            .map(|f| {
                Ok(ConeFace {
                    surface: f
                        .surface
                        .parallel_support(d, dec)?
                        .ok_or_else(|| unsupported("positive join parallel vanished"))?,
                    origins: f.origins.clone(),
                })
            })
            .collect::<Result<Vec<_>, Refusal>>()?;
        let mut boundaries = Vec::new();
        for boundary in &self.boundaries {
            let Curve::Circle(circle) = &boundary.curve else {
                unreachable!()
            };
            let normal = |f: usize, dec: &mut Decider| -> Result<[Real; 2], Refusal> {
                Ok(match &self.faces[f].surface {
                    Surface::Cylinder(_) => [Real::one(), Real::zero()],
                    Surface::Plane(p) => [
                        Real::zero(),
                        Real::from_i64(if sign(&p.normal.dot(&axis), dec)? == Sign::Positive {
                            1
                        } else {
                            -1
                        }),
                    ],
                    _ => unreachable!(),
                })
            };
            let [a, b] = boundary.faces;
            let na = normal(a, dec)?;
            let nb = normal(b, dec)?;
            let turn = sign(&(&na[0] * &nb[1] - &na[1] * &nb[0]), dec)?;
            let radius = reduced(circle.radius_sq.sqrt(), dec)?;
            let shifted = |n: &[Real; 2]| {
                Curve::Circle(Circle {
                    center: circle.center.plus(&axis.scaled(&(d * &n[1]))),
                    normal: axis.clone(),
                    radius_sq: (&radius + &(d * &n[0])).square(),
                    ref_dir: circle.ref_dir.clone(),
                })
            };
            if turn == Sign::Positive {
                let edge = EdgeName::new(
                    self.faces[a]
                        .origins
                        .first()
                        .cloned()
                        .ok_or_else(|| unsupported("a join needs named source faces"))?,
                    self.faces[b]
                        .origins
                        .first()
                        .cloned()
                        .ok_or_else(|| unsupported("a join needs named source faces"))?,
                )
                .map_err(|e| unsupported(&e.to_string()))?;
                let index = faces.len();
                faces.push(ConeFace {
                    surface: Surface::Torus(Torus::new(
                        circle.center.clone(),
                        axis.clone(),
                        radius.clone(),
                        d.clone(),
                        circle.ref_dir.clone(),
                        dec,
                    )?),
                    origins: vec![Origin {
                        step: step.clone(),
                        role: Role::Round { edge },
                    }],
                });
                boundaries.push(ConeBoundary {
                    curve: shifted(&na),
                    faces: [a, index],
                });
                boundaries.push(ConeBoundary {
                    curve: shifted(&nb),
                    faces: [index, b],
                });
            } else if turn == Sign::Negative {
                boundaries.push(ConeBoundary {
                    curve: shifted(&[&na[0] + &nb[0], &na[1] + &nb[1]]),
                    faces: [a, b],
                });
            } else {
                return Err(unsupported(
                    "a cylinder join needs a regular orthogonal normal sector",
                ));
            }
        }
        Ok((faces, boundaries))
    }
}
fn ring(skin: &CylinderSkin, end: usize, faces: [usize; 2]) -> ConeBoundary {
    ConeBoundary {
        curve: Curve::Circle(Circle {
            center: skin.ends[end].clone(),
            normal: skin.cylinder.axis.clone(),
            radius_sq: skin.cylinder.radius_sq.clone(),
            ref_dir: skin.cylinder.ref_dir.clone(),
        }),
        faces,
    }
}

/// Exact equality of the oriented plane/cylinder laws used by this lane.
pub fn same_support(a: &Surface, b: &Surface, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(match (a, b) {
        (Surface::Plane(a), Surface::Plane(b)) => {
            zero(&a.normal.cross(&b.normal).norm_squared(), dec)?
                && sign(&a.normal.dot(&b.normal), dec)? == Sign::Positive
                && zero(&a.origin.minus(&b.origin).dot(&a.normal), dec)?
        }
        (Surface::Cylinder(a), Surface::Cylinder(b)) => {
            zero(&a.axis.cross(&b.axis).norm_squared(), dec)?
                && zero(
                    &a.origin.minus(&b.origin).cross(&a.axis).norm_squared(),
                    dec,
                )?
                && zero(&(&a.radius_sq - &b.radius_sq), dec)?
        }
        _ => false,
    })
}
impl CylinderJoin {
    /// Single selected smooth circular patch. Sharp multi-face selections are
    /// deliberately outside this lane. Circular rims sweep finite cylinders.
    pub fn face_band(
        &self,
        index: usize,
        t: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<crate::cylinder_skin::CylinderBand, Refusal> {
        if let Some((_, band)) = self.bands.iter().find(|(i, _)| *i == index) {
            return band.thicken(2, t, side, step, dec);
        }
        if self.plate.is_some() && index == 0 {
            return Err(unsupported(
                "polygonal plate holes require their selected planar sweep arrangement",
            ));
        }
        let Surface::Plane(plane) = &self.faces[index].surface else {
            return Err(unsupported(
                "selected join patch is not a complete circular band",
            ));
        };
        let mut rings: Vec<_> = self
            .boundaries
            .iter()
            .filter(|b| b.faces.contains(&index))
            .collect();
        if rings.len() == 2 {
            let radius = |b: &ConeBoundary| {
                let Curve::Circle(c) = &b.curve else {
                    unreachable!()
                };
                c.radius_sq.clone()
            };
            if sign(&(radius(rings[0]) - radius(rings[1])), dec)? == Sign::Negative {
                rings.swap(0, 1);
            }
        }
        let source = &self.bands[0].1;
        let axis = &source.cylinder.axis;
        let delta = plane.normal.scaled(
            &(t / &reduced(plane.normal.norm_squared().sqrt(), dec)?
                * Real::from_i64(if matches!(side, ThickenSide::Front) {
                    1
                } else {
                    -1
                })),
        );
        let make = |b: &ConeBoundary, dec: &mut Decider| -> Result<CylinderSkin, Refusal> {
            let Curve::Circle(c) = &b.curve else {
                unreachable!()
            };
            let far = c.center.plus(&delta);
            let forward = sign(&delta.dot(axis), dec)? == Sign::Positive;
            let old = self.faces[index].origins.clone();
            let shifted = old
                .iter()
                .map(|of| Origin {
                    step: step.clone(),
                    role: Role::Offset {
                        of: Box::new(of.clone()),
                    },
                })
                .collect::<Vec<_>>();
            let other = *b.faces.iter().find(|f| **f != index).unwrap();
            let edge = EdgeName::new(
                old.first()
                    .cloned()
                    .ok_or_else(|| unsupported("a rim needs a named patch"))?,
                self.faces[other]
                    .origins
                    .first()
                    .cloned()
                    .ok_or_else(|| unsupported("a rim needs a named neighbour"))?,
            )
            .map_err(|e| unsupported(&e.to_string()))?;
            Ok(CylinderSkin {
                cylinder: catallax_shapes::Cylinder {
                    origin: c.center.clone(),
                    axis: axis.clone(),
                    radius_sq: c.radius_sq.clone(),
                    ref_dir: c.ref_dir.clone(),
                },
                ends: if forward {
                    [c.center.clone(), far]
                } else {
                    [far, c.center.clone()]
                },
                origins: if forward {
                    [
                        old,
                        shifted,
                        vec![Origin {
                            step: step.clone(),
                            role: Role::Rim { edge },
                        }],
                    ]
                } else {
                    [
                        shifted,
                        old,
                        vec![Origin {
                            step: step.clone(),
                            role: Role::Rim { edge },
                        }],
                    ]
                },
            })
        };
        Ok(crate::cylinder_skin::CylinderBand {
            outer: make(rings[0], dec)?,
            inner: if rings.len() == 2 {
                Some(make(rings[1], dec)?)
            } else {
                None
            },
        })
    }
}

impl CylinderJoin {
    /// Exact adjacency selection is shared; nearest disambiguation is not inferred.
    pub fn select(&self, refs: &[FaceRef]) -> Result<usize, Refusal> {
        let mut chosen = None;
        for reference in refs {
            let candidates: Vec<_> = self
                .faces
                .iter()
                .enumerate()
                .filter(|(i, f)| {
                    f.origins.iter().any(|o| reference.face.matches(o))
                        && reference.touching.iter().all(|p| {
                            self.boundaries.iter().any(|b| {
                                b.faces.contains(i)
                                    && b.faces.iter().any(|j| {
                                        *j != *i
                                            && self.faces[*j].origins.iter().any(|o| p.matches(o))
                                    })
                            })
                        })
                })
                .map(|(i, _)| i)
                .collect();
            match candidates.as_slice() {
                [] => {
                    return Err(Refusal::new(
                        "reference.lost-face",
                        format!("no bounded join patch resolves {reference:?}"),
                    ));
                }
                [one] if chosen.is_none_or(|c| c == *one) => chosen = Some(*one),
                _ => {
                    return Err(unsupported(
                        "sharp multi-face join thickening is outside the bounded lane",
                    ));
                }
            }
            if reference.near.is_some() && candidates.len() != 1 {
                return Err(unsupported(
                    "nearest finite join disambiguation is unsupported",
                ));
            }
        }
        chosen.ok_or_else(|| {
            Refusal::new(
                "reference.lost-face",
                "thicken needs a surviving join patch",
            )
        })
    }
}

impl CylinderJoin {
    /// Finite normal sweep of the rectangular attachment face with its circular
    /// hole. The hole cylinder spans exactly the same two cap planes.
    pub fn plate_sweep(
        &self,
        index: usize,
        t: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<(ConvexCell, CylinderSkin)>, Refusal> {
        let Some((cell, top)) = &self.plate else {
            return Ok(None);
        };
        if index != 0 {
            return Ok(None);
        }
        let prism = cell.face_sweep(*top, t, side, step, dec)?;
        let source = &self.bands[0].1;
        let direction = source.cylinder.axis.scaled(
            &(t / &reduced(source.cylinder.axis.norm_squared().sqrt(), dec)?
                * Real::from_i64(if matches!(side, ThickenSide::Front) {
                    1
                } else {
                    -1
                })),
        );
        let at = source.ends[0].clone();
        let far = at.plus(&direction);
        let edge = EdgeName::new(
            self.faces[0]
                .origins
                .first()
                .cloned()
                .ok_or_else(|| unsupported("a plate rim requires a named face"))?,
            self.faces[1]
                .origins
                .first()
                .cloned()
                .ok_or_else(|| unsupported("a plate rim requires a named neighbour"))?,
        )
        .map_err(|e| unsupported(&e.to_string()))?;
        let hole = CylinderSkin {
            cylinder: source.cylinder.clone(),
            ends: if matches!(side, ThickenSide::Front) {
                [at, far]
            } else {
                [far, at]
            },
            origins: [
                Vec::new(),
                Vec::new(),
                vec![Origin {
                    step: step.clone(),
                    role: Role::Rim { edge },
                }],
            ],
        };
        Ok(Some((prism, hole)))
    }
}
