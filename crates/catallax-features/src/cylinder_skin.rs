//! Exact complete capped-cylinder definitions for distance features.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, FaceRef, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Cylinder, Plane, Point3, Surface};

use crate::profile::sign;

/// A complete cylindrical skin, ordered lower cap, upper cap, wall.
#[derive(Clone, Debug)]
pub struct CylinderSkin {
    /// Exact cylindrical support.
    pub cylinder: Cylinder,
    /// Axis points on the lower and upper caps.
    pub ends: [Point3; 2],
    /// Original cap and wall generators.
    pub origins: [Vec<Origin>; 3],
}

/// An exact radial cylinder band with planar normal-sweep rims.
#[derive(Clone, Debug)]
pub struct CylinderBand {
    /// Complete outer cylinder.
    pub outer: CylinderSkin,
    /// Inner cylinder, sharing the cap coordinates when the ends are rims.
    pub inner: Option<CylinderSkin>,
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.cylinder-domain", detail)
}
fn names(step: &StepId, originals: &[Origin], inner: bool) -> Vec<Origin> {
    originals
        .iter()
        .map(|of| Origin {
            step: step.clone(),
            role: if inner {
                Role::Inner {
                    of: Box::new(of.clone()),
                }
            } else {
                Role::Offset {
                    of: Box::new(of.clone()),
                }
            },
        })
        .collect()
}
impl CylinderSkin {
    /// Recognises two complete perpendicular caps and a complete cylindrical wall.
    /// Completeness of the trim inventory is a caller obligation.
    pub fn new(
        surfaces: &[(Surface, Vec<Origin>)],
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if surfaces.len() != 3 {
            return Ok(None);
        }
        let mut wall = None;
        let mut caps = [None, None];
        for (i, (surface, _)) in surfaces.iter().enumerate() {
            match surface {
                Surface::Cylinder(c) if wall.is_none() => wall = Some((i, c)),
                Surface::Plane(_) => (),
                _ => return Ok(None),
            }
        }
        let Some((wall_id, cylinder)) = wall else {
            return Ok(None);
        };
        for (i, (surface, _)) in surfaces.iter().enumerate() {
            let Surface::Plane(plane) = surface else {
                continue;
            };
            if sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero {
                return Ok(None);
            }
            let orientation = sign(&plane.normal.dot(&cylinder.axis), dec)?;
            let slot = match orientation {
                Sign::Negative => 0,
                Sign::Positive => 1,
                Sign::Zero => return Ok(None),
            };
            if caps[slot].is_some() {
                return Ok(None);
            }
            let h = plane.origin.minus(&cylinder.origin).dot(&cylinder.axis)
                / cylinder.axis.norm_squared();
            caps[slot] = Some((i, cylinder.origin.plus(&cylinder.axis.scaled(&h))));
        }
        let [Some((lower, a)), Some((upper, b))] = caps else {
            return Ok(None);
        };
        if sign(&b.minus(&a).dot(&cylinder.axis), dec)? != Sign::Positive {
            return Ok(None);
        }
        Ok(Some(Self {
            cylinder: cylinder.clone(),
            ends: [a, b],
            origins: [
                surfaces[lower].1.clone(),
                surfaces[upper].1.clone(),
                surfaces[wall_id].1.clone(),
            ],
        }))
    }
    /// Complete outward-oriented cap and wall equations.
    pub fn surfaces(&self) -> [Surface; 3] {
        [
            Surface::Plane(Plane {
                origin: self.ends[0].clone(),
                normal: self.cylinder.axis.scaled(&Real::from_i64(-1)),
                ref_dir: self.cylinder.ref_dir.clone(),
            }),
            Surface::Plane(Plane {
                origin: self.ends[1].clone(),
                normal: self.cylinder.axis.clone(),
                ref_dir: self.cylinder.ref_dir.clone(),
            }),
            Surface::Cylinder(self.cylinder.clone()),
        ]
    }
    /// Resolves the complete caps and wall, including their circular adjacency.
    pub fn select(&self, refs: &[FaceRef]) -> Result<[bool; 3], Refusal> {
        let mut selected = [false; 3];
        for reference in refs {
            let candidates: Vec<_> = (0..3)
                .filter(|i| {
                    self.origins[*i].iter().any(|o| reference.face.matches(o))
                        && reference.touching.iter().all(|p| {
                            (0..3).any(|j| {
                                j != *i
                                    && (*i == 2 || j == 2)
                                    && self.origins[j].iter().any(|o| p.matches(o))
                            })
                        })
                })
                .collect();
            if candidates.is_empty() {
                return Err(Refusal::new(
                    "reference.lost-face",
                    format!("no cylinder face resolves {reference:?}"),
                ));
            }
            if reference.near.is_some() && candidates.len() > 1 {
                return Err(bad(
                    "nearest disambiguation among cylinder trims needs exact disk and band distance",
                ));
            }
            for i in candidates {
                selected[i] = true;
            }
        }
        Ok(selected)
    }
    /// Erodes exactly by radius and cap clearances. An overlarge erosion is empty.
    pub fn eroded(&self, distance: &Real, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        self.cavity(&[false; 3], distance, None, dec)
    }
    /// Distance from kept complete caps and wall, with open faces left unshifted.
    pub fn cavity(
        &self,
        open: &[bool; 3],
        t: &Real,
        step: Option<&StepId>,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let mut result = self.clone();
        if !open[2] {
            let Some(Surface::Cylinder(cylinder)) =
                Surface::Cylinder(self.cylinder.clone()).parallel_support(&(-t), dec)?
            else {
                return Ok(None);
            };
            result.cylinder = cylinder;
        }
        let direction = self
            .cylinder
            .axis
            .scaled(&(&Real::one() / &self.cylinder.axis.norm_squared().sqrt()));
        for (i, opened) in open.iter().take(2).enumerate() {
            if !opened {
                result.ends[i] = result.ends[i]
                    .plus(&direction.scaled(&(t * &Real::from_i64(if i == 0 { 1 } else { -1 }))));
            }
        }
        if sign(
            &result.ends[1]
                .minus(&result.ends[0])
                .dot(&self.cylinder.axis),
            dec,
        )? != Sign::Positive
        {
            return Ok(None);
        }
        if let Some(step) = step {
            for (i, opened) in open.iter().enumerate() {
                if !opened {
                    result.origins[i] = names(step, &self.origins[i], true);
                }
            }
        }
        Ok(Some(result))
    }
    /// Normal sweep of one complete smooth face, retaining its exact circular rims.
    pub fn thicken(
        &self,
        selected: usize,
        t: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<CylinderBand, Refusal> {
        let front = matches!(side, ThickenSide::Front);
        let mut outer = self.clone();
        let mut inner = None;
        if selected < 2 {
            let index = selected;
            let signum = if (index == 1) == front { 1 } else { -1 };
            let direction = self.cylinder.axis.scaled(&(t * &Real::from_i64(signum)));
            let direction =
                direction.scaled(&(&Real::one() / &self.cylinder.axis.norm_squared().sqrt()));
            let far = self.ends[index].plus(&direction);
            let far_index = if signum > 0 { 1 } else { 0 };
            outer.ends[far_index] = far;
            outer.ends[1 - far_index] = self.ends[index].clone();
            outer.origins[far_index] = names(step, &self.origins[index], false);
            outer.origins[1 - far_index] = self.origins[index].clone();
            outer.origins[2] = self.rim(index, step)?;
        } else {
            let r = self.cylinder.radius();
            let delta = if front { t.clone() } else { -t };
            let far_origins = names(step, &self.origins[2], false);
            if front {
                let Some(Surface::Cylinder(cylinder)) =
                    Surface::Cylinder(self.cylinder.clone()).parallel_support(&delta, dec)?
                else {
                    unreachable!("positive radial growth")
                };
                outer.cylinder = cylinder;
                outer.origins[2] = far_origins;
                inner = Some(self.clone());
            } else {
                let fold = t - &(&r * &Real::from_i64(2));
                match sign(&fold, dec)? {
                    Sign::Positive => {
                        let Some(Surface::Cylinder(cylinder)) =
                            Surface::Cylinder(self.cylinder.clone())
                                .parallel_support(&fold, dec)?
                        else {
                            unreachable!("positive folded growth")
                        };
                        outer.cylinder = cylinder;
                        outer.origins[2] = far_origins;
                    }
                    Sign::Zero => {
                        outer.origins[2].extend(far_origins);
                        outer.origins[2].sort();
                        outer.origins[2].dedup();
                    }
                    Sign::Negative => {
                        if let Some(Surface::Cylinder(cylinder)) =
                            Surface::Cylinder(self.cylinder.clone())
                                .parallel_support(&delta, dec)?
                        {
                            let mut hole = self.clone();
                            hole.cylinder = cylinder;
                            hole.origins[2] = far_origins;
                            inner = Some(hole);
                        }
                    }
                }
            }
            for i in 0..2 {
                outer.origins[i] = self.rim(i, step)?;
                if let Some(hole) = &mut inner {
                    hole.origins[i] = outer.origins[i].clone();
                }
            }
        }
        Ok(CylinderBand { outer, inner })
    }
    fn rim(&self, cap: usize, step: &StepId) -> Result<Vec<Origin>, Refusal> {
        let a = self.origins[cap]
            .first()
            .ok_or_else(|| bad("a cylinder rim requires named caps"))?;
        let b = self.origins[2]
            .first()
            .ok_or_else(|| bad("a cylinder rim requires a named wall"))?;
        Ok(vec![Origin {
            step: step.clone(),
            role: Role::Rim {
                edge: EdgeName::new(a.clone(), b.clone()).map_err(|e| bad(&e.to_string()))?,
            },
        }])
    }
}
