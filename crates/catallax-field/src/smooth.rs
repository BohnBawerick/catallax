//! Circular distance blends. Operand residuals never substitute for signed distance.

use crate::{
    Face, Instruction, Location, Plan,
    geometry::{self, Leaf, Primitive},
    metric::{Arithmetic, Kind, PRECISION},
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, Origin, Role, StepId};
use catallax_shapes::{Isometry, RPoint3};

#[derive(Clone, Debug)]
pub(crate) struct Smooth {
    operands: [Primitive; 2],
    faces: [Vec<Face>; 2],
    kind: BooleanKind,
    radius: Rational,
    origin: Origin,
}

fn min(a: Real, b: Real, dec: &mut Decider) -> Result<Real, Refusal> {
    Ok(if geometry::sign(&(&a - &b), dec)? == Sign::Positive {
        b
    } else {
        a
    })
}
fn max(a: Real, b: Real, dec: &mut Decider) -> Result<Real, Refusal> {
    min(-a, -b, dec).map(|x| -x)
}

fn exact_distance(a: &Primitive, p: &RPoint3, dec: &mut Decider) -> Result<Real, Refusal> {
    let q = a.inverse_placement().apply_point(&p.to_real());
    let [x, y, z] = [q.x, q.y, q.z];
    let n = a.parameters.each_ref().map(|q| Real::from(q.clone()));
    let zero = Real::zero();
    Ok(match a.kind {
        Kind::Sphere => (x.square() + y.square() + z.square()).sqrt() - n[0].clone(),
        Kind::Torus => {
            (((x.square() + y.square()).sqrt() - n[0].clone()).square() + z.square()).sqrt()
                - n[1].clone()
        }
        Kind::Block => {
            let q = [
                max(-&x, &x - &n[0], dec)?,
                max(-&y, &y - &n[1], dec)?,
                max(-&z, &z - &n[2], dec)?,
            ];
            let outside = q
                .iter()
                .try_fold(Real::zero(), |sum, v| {
                    Ok::<_, Refusal>(sum + max(v.clone(), zero.clone(), dec)?.square())
                })?
                .sqrt();
            outside
                + min(
                    max(max(q[0].clone(), q[1].clone(), dec)?, q[2].clone(), dec)?,
                    zero,
                    dec,
                )?
        }
        Kind::Round => {
            let rho = (x.square() + y.square()).sqrt();
            let dr = &n[1] - &n[0];
            let t = min(
                max(
                    ((&rho - &n[0]) * dr.clone() + &z * &n[2]) / (dr.square() + n[2].square()),
                    zero.clone(),
                    dec,
                )?,
                Real::one(),
                dec,
            )?;
            let side = ((rho.clone() - (n[0].clone() + &t * &dr)).square()
                + (z.clone() - &t * &n[2]).square())
            .sqrt();
            let bottom = (max(&rho - &n[0], zero.clone(), dec)?.square() + z.square()).sqrt();
            let top = (max(&rho - &n[1], zero, dec)?.square() + (&z - &n[2]).square()).sqrt();
            let d = min(min(side, bottom, dec)?, top, dec)?;
            if a.location(&[p.x.clone(), p.y.clone(), p.z.clone()], dec)? == Location::Inside {
                -d
            } else {
                d
            }
        }
    })
}

impl Plan {
    pub(crate) fn require_analytic_boundary(&self) -> Result<(), Refusal> {
        if self.leaves.iter().enumerate().any(|(i, l)| {
            self.leaf_registers.get(i).is_some_and(Option::is_some) && matches!(l, Leaf::Smooth(_))
        }) {
            return Err(Refusal::new(
                "smooth.boundary-unsupported",
                "exact membership and terminal boundary origins are available; 0 certified smooth display/print patches admitted",
            ));
        }
        Ok(())
    }

    pub(crate) fn smooth(
        &mut self,
        kind: BooleanKind,
        nodes: [usize; 2],
        faces: [&Vec<usize>; 2],
        radius: &Rational,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if radius.sign() != Sign::Positive {
            return Err(Refusal::new(
                "smooth.non-positive-radius",
                format!("radius {radius} mm must exceed 0"),
            ));
        }
        let operand = |index: usize| -> Result<Primitive, Refusal> {
            if let Instruction::Leaf(leaf) = self.instructions[nodes[index]]
                && let Leaf::Primitive(p) = &self.leaves[leaf]
            {
                return Ok(p.clone());
            }
            Err(Refusal::new(
                "smooth.distance-unsupported",
                format!(
                    "operand {index}: {} tape instructions; exact stock distance required, residual substitution forbidden",
                    self.instructions.len()
                ),
            ))
        };
        let smooth = Smooth {
            operands: [operand(0)?, operand(1)?],
            faces: faces.map(|ids| ids.iter().map(|&i| self.faces[i].clone()).collect()),
            kind,
            radius: radius.clone(),
            origin: Origin {
                step: step.clone(),
                role: Role::Smooth,
            },
        };
        self.leaf(Leaf::Smooth(Box::new(smooth)), vec![], dec)
    }

    /// Exact boundary names for a terminal circular smooth blend. The analytic
    /// face registry has no carrier for its new implicit patch. Other roots refuse.
    pub fn smooth_origins_checked(
        &self,
        point: &RPoint3,
        dec: &mut Decider,
    ) -> Result<Vec<Origin>, Refusal> {
        if let Instruction::Leaf(i) = self.instructions[self.root]
            && let Leaf::Smooth(s) = &self.leaves[i]
        {
            return s.origins(point, dec);
        }
        Err(Refusal::new(
            "smooth.not-terminal",
            "boundary-name query requires one terminal smooth blend",
        ))
    }
}

impl Smooth {
    fn values(&self, p: &RPoint3, dec: &mut Decider) -> Result<(Real, Real, Real), Refusal> {
        let mut a = exact_distance(&self.operands[0], p, dec)?;
        let mut b = exact_distance(&self.operands[1], p, dec)?;
        if self.kind == BooleanKind::Subtract {
            b = -b;
        }
        let union = self.kind == BooleanKind::Union;
        if !union {
            a = -a;
            b = -b;
        }
        let k = Real::from(self.radius.clone());
        let x = max(&k - &a, Real::zero(), dec)?;
        let y = max(&k - &b, Real::zero(), dec)?;
        let f = max(k, min(a.clone(), b.clone(), dec)?, dec)? - (x.square() + y.square()).sqrt();
        Ok((if union { f } else { -f }, a, b))
    }
    pub(crate) fn location(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let (f, _, _) =
            self.values(&RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()), dec)?;
        Ok(match geometry::sign(&f, dec)? {
            Sign::Negative => Location::Inside,
            Sign::Positive => Location::Outside,
            Sign::Zero => Location::Boundary,
        })
    }
    pub(crate) fn residual(&self, p: &RPoint3, dec: &mut Decider) -> Result<Ball, Refusal> {
        // The circular formula is sqrt(2)-Lipschitz in model space. Dividing
        // by 2 preserves the zero set and certifies the distance lower bound.
        let (f, _, _) = self.values(p, dec)?;
        geometry::enclose(&(f / Real::from_i64(2)), PRECISION, dec)
    }
    pub(crate) fn origins(&self, p: &RPoint3, dec: &mut Decider) -> Result<Vec<Origin>, Refusal> {
        let (f, a, b) = self.values(p, dec)?;
        if geometry::sign(&f, dec)? != Sign::Zero {
            return Ok(vec![]);
        }
        let k = Real::from(self.radius.clone());
        if geometry::sign(&(&a - &k), dec)? == Sign::Negative
            && geometry::sign(&(&b - &k), dec)? == Sign::Negative
        {
            return Ok(vec![self.origin.clone()]);
        }
        let xyz = [p.x.clone(), p.y.clone(), p.z.clone()];
        let mut origins = Vec::new();
        for (i, faces) in self.faces.iter().enumerate() {
            if geometry::sign(if i == 0 { &a } else { &b }, dec)? != Sign::Zero {
                continue;
            }
            for face in faces {
                if crate::polynomial_sign(&face.polynomial, &xyz, dec)? == Sign::Zero {
                    origins.push(face.origin.clone());
                }
            }
        }
        origins.sort();
        origins.dedup();
        Ok(origins)
    }
    pub(crate) fn transform(&mut self, _pose: &Isometry) -> Result<(), Refusal> {
        Err(Refusal::new(
            "smooth.copy-names-unsupported",
            "smooth copies need transformed implicit boundary names; 0 copied smooth patches admitted",
        ))
    }
    pub(crate) fn radius(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        let mut radius = Rational::zero();
        for p in &self.operands {
            radius = radius.max(
                Leaf::Primitive(p.clone())
                    .display_radius(dec)?
                    .expect("finite stock"),
            );
        }
        Ok(radius + &self.radius * &Rational::from_i64(2))
    }
    pub(crate) fn fast<A: Arithmetic>(&self, p: [f64; 3]) -> f64 {
        let mut a = self.operands[0].fast_distance::<A>(p);
        let mut b = self.operands[1].fast_distance::<A>(p);
        if self.kind == BooleanKind::Subtract {
            b = -b;
        }
        let union = self.kind == BooleanKind::Union;
        if !union {
            a = -a;
            b = -b;
        }
        let k = self.radius.to_f64_nearest();
        let x = (k - a).max(0.0);
        let y = (k - b).max(0.0);
        let f = k.max(a.min(b)) - A::sqrt(x * x + y * y);
        if union { f / 2.0 } else { -f / 2.0 }
    }
}
