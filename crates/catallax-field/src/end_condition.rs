//! Independent sign queries for a certified single-branch curved extrusion.
//!
//! Queries read exact support equations and the retained ray chart. They never
//! build a B-rep, call B-rep containment, or treat a polynomial as true distance.

use catallax_features::end_condition::curved::{AuthoredCap, Branch, Location};
use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign, StrictSign};
use catallax_recipe::{Origin, Role, StepId};
use catallax_shapes::{Point3, RPoint3};

/// A field definition retaining the certified profile and analytic cap.
#[derive(Clone, Debug)]
pub struct CurvedField {
    cap: AuthoredCap,
    step: StepId,
}

fn sign(x: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    dec.check_deadline()?;
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::Refused(r) => Err(r),
        Outcome::NotExact(_) => Err(Refusal::new(
            "field.curved-limit-unresolved",
            "curved cap sign is not certified",
        )),
    }
}

impl CurvedField {
    /// Retains an independently compiled field definition. No topology is read.
    pub fn new(cap: AuthoredCap, step: StepId) -> Self {
        Self { cap, step }
    }

    fn coordinates(&self, point: &Point3, dec: &mut Decider) -> Result<[Real; 3], Refusal> {
        let f = self.cap.definition().frame();
        let d = point.minus(&f.origin);
        let determinant = f.u.cross(&f.v).dot(&f.direction);
        // The definition constructor proves this determinant nonzero.
        dec.check_deadline()?;
        Ok([
            d.dot(&f.v.cross(&f.direction)) / determinant.clone(),
            d.dot(&f.direction.cross(&f.u)) / determinant.clone(),
            d.dot(&f.u.cross(&f.v)) / determinant,
        ])
    }

    /// Exact set classification using the polynomial and its ray derivative,
    /// independently of the B-rep builder and the definition's height evaluator.
    pub fn location(&self, point: &RPoint3, dec: &mut Decider) -> Result<Location, Refusal> {
        let point = point.to_real();
        let [u, v, t] = self.coordinates(&point, dec)?;
        let mut boundary = false;
        if let Some(radius) = self.cap.definition().circular_radius() {
            let radial = sign(&(&u * &u + &v * &v - Real::from(radius * radius)), dec)?;
            if radial == Sign::Positive {
                return Ok(Location::Outside);
            }
            boundary = radial == Sign::Zero;
        }
        for [a, b, c] in self.cap.definition().side_constraints() {
            let side = sign(
                &(&Real::from(a.clone()) * &u
                    + &Real::from(b.clone()) * &v
                    + Real::from(c.clone())),
                dec,
            )?;
            if side == Sign::Positive {
                return Ok(Location::Outside);
            }
            boundary |= side == Sign::Zero;
        }
        let start = sign(&t, dec)?;
        if start == Sign::Negative {
            return Ok(Location::Outside);
        }
        let surface = self.cap.definition().scene().surface();
        let value = surface
            .implicit()
            .ok_or_else(|| {
                Refusal::new(
                    "field.curved-limit-unresolved",
                    "limiting carrier has no implicit polynomial",
                )
            })?
            .eval_real(&[point.x.clone(), point.y.clone(), point.z.clone()]);
        let value = sign(&value, dec)?;
        match self.cap.definition().branch() {
            Branch::Exit if value == Sign::Positive => return Ok(Location::Outside),
            Branch::Entry => {
                if value == Sign::Negative {
                    return Ok(Location::Outside);
                }
                let derivative = surface
                    .gradient(&point)
                    .ok_or_else(|| {
                        Refusal::new(
                            "field.curved-limit-unresolved",
                            "limiting carrier has no analytic gradient",
                        )
                    })?
                    .dot(&self.cap.definition().frame().direction);
                if sign(&derivative, dec)? != Sign::Negative {
                    return Ok(Location::Outside);
                }
            }
            _ => (),
        }
        Ok(if boundary || start == Sign::Zero || value == Sign::Zero {
            Location::Boundary
        } else {
            Location::Inside
        })
    }

    /// Origins of the actual boundary contributions at a certified boundary
    /// point. Edge queries return the union of their incident contributions.
    pub fn origins(&self, point: &RPoint3, dec: &mut Decider) -> Result<Vec<Origin>, Refusal> {
        if self.location(point, dec)? != Location::Boundary {
            return Ok(vec![]);
        }
        let p = point.to_real();
        let [u, v, t] = self.coordinates(&p, dec)?;
        let mut roles = Vec::new();
        if let Some(radius) = self.cap.definition().circular_radius()
            && sign(&(&u * &u + &v * &v - Real::from(radius * radius)), dec)? == Sign::Zero
        {
            roles.push(Role::Side {
                segment: self.cap.side_ids()[0].clone(),
                path: None,
            });
        }
        if sign(&t, dec)? == Sign::Zero {
            roles.push(Role::Start);
        }
        let surface = self.cap.definition().scene().surface();
        let polynomial = surface.implicit().ok_or_else(|| {
            Refusal::new(
                "field.curved-limit-unresolved",
                "limiting carrier has no implicit polynomial",
            )
        })?;
        if sign(&polynomial.eval_real(&[p.x, p.y, p.z]), dec)? == Sign::Zero {
            roles.push(Role::End);
        }
        for ([a, b, c], id) in self
            .cap
            .definition()
            .side_constraints()
            .iter()
            .zip(self.cap.side_ids())
        {
            if sign(
                &(&Real::from(a.clone()) * &u
                    + &Real::from(b.clone()) * &v
                    + Real::from(c.clone())),
                dec,
            )? == Sign::Zero
            {
                roles.push(Role::Side {
                    segment: id.clone(),
                    path: None,
                });
            }
        }
        let mut origins: Vec<_> = roles
            .into_iter()
            .map(|role| Origin {
                step: self.step.clone(),
                role,
            })
            .collect();
        origins.sort();
        origins.dedup();
        Ok(origins)
    }
}
