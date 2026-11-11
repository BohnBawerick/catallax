//! Analytic cutting surfaces and their exact homogeneous patch substitutions.

use super::*;
use crate::{Cone, Cylinder, Sphere, Surface, Torus};

mod torus;
use torus::torus_net;

/// An analytic surface whose section of a rational patch has a complete
/// parameter-cell search. Its implicit equation is substituted exactly into
/// the patch's homogeneous coordinates.
#[derive(Clone, Debug)]
pub enum SectionSurface {
    /// Plane, with the cut equation linear in homogeneous coordinates.
    Plane(Plane),
    /// Circular cylinder, with the cut equation quadratic.
    Cylinder(Cylinder),
    /// Sphere, with the cut equation quadratic.
    Sphere(Sphere),
    /// The retained nappe of a circular cone, with a quadratic cut equation.
    Cone(Cone),
    /// Ring torus, with a quartic cut equation.
    Torus(Torus),
}

/// Stable refusal codes for one cutting-surface family.
pub(super) struct Codes {
    pub(super) tangent: &'static str,
    pub(super) candidate_domain: &'static str,
    pub(super) candidate_endpoint: &'static str,
    pub(super) candidate_cap: &'static str,
}

impl SectionSurface {
    /// The admitted cutting surface for an analytic carrier, if any.
    pub fn of(surface: &Surface) -> Option<Self> {
        Some(match surface {
            Surface::Plane(s) => Self::Plane(s.clone()),
            Surface::Cylinder(s) => Self::Cylinder(s.clone()),
            Surface::Sphere(s) => Self::Sphere(s.clone()),
            Surface::Cone(s) => Self::Cone(s.clone()),
            Surface::Torus(s) => Self::Torus(s.clone()),
            _ => return None,
        })
    }
    /// The same carrier as a general surface.
    pub fn surface(&self) -> Surface {
        match self {
            Self::Plane(s) => Surface::Plane(s.clone()),
            Self::Cylinder(s) => Surface::Cylinder(s.clone()),
            Self::Sphere(s) => Surface::Sphere(s.clone()),
            Self::Cone(s) => Surface::Cone(s.clone()),
            Self::Torus(s) => Surface::Torus(s.clone()),
        }
    }
    pub(super) fn codes(&self) -> Codes {
        match self {
            Self::Plane(_) => Codes {
                tangent: "intersect.patch-plane-tangent",
                candidate_domain: "intersect.patch-plane-candidate-domain",
                candidate_endpoint: "intersect.patch-plane-candidate-endpoint",
                candidate_cap: "intersect.patch-plane-candidate-cap",
            },
            Self::Cylinder(_) => Codes {
                tangent: "intersect.patch-cylinder-tangent",
                candidate_domain: "intersect.patch-cylinder-candidate-domain",
                candidate_endpoint: "intersect.patch-cylinder-candidate-endpoint",
                candidate_cap: "intersect.patch-cylinder-candidate-cap",
            },
            Self::Sphere(_) => Codes {
                tangent: "intersect.patch-sphere-tangent",
                candidate_domain: "intersect.patch-sphere-candidate-domain",
                candidate_endpoint: "intersect.patch-sphere-candidate-endpoint",
                candidate_cap: "intersect.patch-sphere-candidate-cap",
            },
            Self::Cone(_) => Codes {
                tangent: "intersect.patch-cone-tangent",
                candidate_domain: "intersect.patch-cone-candidate-domain",
                candidate_endpoint: "intersect.patch-cone-candidate-endpoint",
                candidate_cap: "intersect.patch-cone-candidate-cap",
            },
            Self::Torus(_) => Codes {
                tangent: "intersect.patch-torus-tangent",
                candidate_domain: "intersect.patch-torus-candidate-domain",
                candidate_endpoint: "intersect.patch-torus-candidate-endpoint",
                candidate_cap: "intersect.patch-torus-candidate-cap",
            },
        }
    }
    pub(super) fn regular_contact_code(&self) -> Option<&'static str> {
        match self {
            Self::Cone(_) | Self::Torus(_) => Some(self.codes().tangent),
            _ => None,
        }
    }
    pub(super) fn unresolved(&self, detail: &str) -> Refusal {
        Refusal::new(
            match self {
                Self::Cone(_) => "intersect.patch-cone-nonregular",
                Self::Torus(_) => "intersect.patch-torus-nonregular",
                _ => "intersect.rational-patch-section",
            },
            detail,
        )
    }
    pub(crate) fn validate(&self, dec: &mut Decider) -> Result<(), Refusal> {
        let positive = |value: &Real, detail: &str, dec: &mut Decider| {
            if exact_sign(value, dec)? == Sign::Positive {
                Ok(())
            } else {
                Err(failure(detail))
            }
        };
        match self {
            Self::Plane(plane) => {
                math::validate_point(&plane.origin, dec)?;
                math::validate_vector(&plane.normal, dec)?;
                positive(
                    &plane.normal.norm_squared(),
                    "a cutting plane needs a nonzero normal",
                    dec,
                )
            }
            Self::Cylinder(cylinder) => {
                math::validate_point(&cylinder.origin, dec)?;
                math::validate_vector(&cylinder.axis, dec)?;
                positive(
                    &cylinder.axis.norm_squared(),
                    "a cutting cylinder needs a nonzero axis",
                    dec,
                )?;
                positive(
                    &cylinder.radius_sq,
                    "a cutting cylinder needs a positive radius",
                    dec,
                )
            }
            Self::Sphere(sphere) => {
                math::validate_point(&sphere.center, dec)?;
                positive(
                    &sphere.radius,
                    "a cutting sphere needs a positive radius",
                    dec,
                )
            }
            Self::Cone(cone) => {
                math::validate_point(&cone.apex, dec)?;
                math::validate_vector(&cone.axis, dec)?;
                positive(
                    &cone.axis.norm_squared(),
                    "a cutting cone needs a nonzero axis",
                    dec,
                )?;
                positive(&cone.slope, "a cutting cone needs a positive slope", dec)
            }
            Self::Torus(torus) => {
                math::validate_point(&torus.center, dec)?;
                math::validate_vector(&torus.axis, dec)?;
                positive(
                    &torus.axis.norm_squared(),
                    "a cutting torus needs a nonzero axis",
                    dec,
                )?;
                positive(
                    &torus.minor,
                    "a cutting torus needs a positive minor radius",
                    dec,
                )?;
                positive(
                    &(&torus.major - &torus.minor),
                    "a cutting torus must be a ring torus",
                    dec,
                )
            }
        }
    }
    /// A cone's polynomial contains both nappes. The homogeneous height net
    /// proves this whole span retained or excluded before any reduction.
    pub(super) fn retains_span(
        &self,
        controls: &[[Real; 4]],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let Self::Cone(cone) = self else {
            return Ok(true);
        };
        let origin = [&cone.apex.x, &cone.apex.y, &cone.apex.z];
        let axis = [&cone.axis.x, &cone.axis.y, &cone.axis.z];
        let mut signs = Vec::with_capacity(controls.len());
        for h in controls {
            let mut height = Real::zero();
            for i in 0..3 {
                let offset = exact(&h[i], &(-exact(origin[i], &h[3], true, dec)?), false, dec)?;
                let term = exact(&offset, axis[i], true, dec)?;
                height = exact(&height, &term, false, dec)?;
            }
            signs.push(exact_sign(&height, dec)?);
        }
        if signs.iter().all(|s| *s == Sign::Positive) {
            return Ok(true);
        }
        if signs.iter().all(|s| *s == Sign::Negative) {
            return Ok(false);
        }
        Err(Refusal::new(
            "intersect.patch-cone-nappe",
            "the span height net does not prove one strict nappe; apex and nappe-boundary contacts need a separate certificate",
        ))
    }
    /// The implicit equation at a point, zero exactly on the surface.
    pub(super) fn value_at(&self, point: &Point3) -> Real {
        self.surface()
            .value(point)
            .expect("admitted cutting surfaces have equations")
    }
    /// Homogeneous quadratic form `sum c X_i X_j` of a quadric, indexing
    /// `X, Y, Z, W`. Each monomial of the affine equation keeps its exact
    /// coefficient; missing powers of `W` complete it to degree two.
    fn quadratic_form(&self) -> Vec<([usize; 2], Real)> {
        let implicit = self
            .surface()
            .implicit()
            .expect("admitted cutting surfaces have equations");
        let mut form = Vec::new();
        for (monomial, coefficient) in implicit.terms() {
            if coefficient.rational_ref().is_some_and(Rational::is_zero) {
                continue;
            }
            let mut indices = [3, 3];
            let mut slot = 0;
            for (axis, &power) in monomial.exponents().iter().enumerate() {
                for _ in 0..power {
                    indices[slot] = axis;
                    slot += 1;
                }
            }
            form.push((indices, coefficient));
        }
        form
    }
    /// The cut equation in a span's absolute knot parameters, homogeneous of
    /// this surface's degree in the patch coordinates. Its zero set in the
    /// span is the section, and its sign is the side on positive weights.
    pub(super) fn equation(&self, span: &Span) -> RealPoly {
        match self {
            Self::Plane(plane) => {
                let point = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
                let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
                (0..3).fold(RealPoly::zero(2), |sum, i| {
                    sum + (&span.coordinates[i] - &span.coordinates[3].scale(point[i]))
                        .scale(normal[i])
                })
            }
            Self::Cylinder(_) | Self::Sphere(_) | Self::Cone(_) => self
                .quadratic_form()
                .into_iter()
                .fold(RealPoly::zero(2), |sum, ([i, j], coefficient)| {
                    sum + (&span.coordinates[i] * &span.coordinates[j]).scale(&coefficient)
                }),
            Self::Torus(torus) => {
                let d: [RealPoly; 3] = std::array::from_fn(|i| {
                    &span.coordinates[i]
                        - &span.coordinates[3]
                            .scale([&torus.center.x, &torus.center.y, &torus.center.z][i])
                });
                let q = d.iter().fold(RealPoly::zero(2), |sum, p| sum + p * p);
                let h = d
                    .iter()
                    .zip([&torus.axis.x, &torus.axis.y, &torus.axis.z])
                    .fold(RealPoly::zero(2), |sum, (p, a)| sum + p.scale(a));
                let w2 = &span.coordinates[3] * &span.coordinates[3];
                let a = &q + &w2.scale(&(&torus.major.square() - &torus.minor.square()));
                let radial = q - (&h * &h).scale(&(&Real::one() / &torus.axis.norm_squared()));
                &a * &a - (&radial * &w2).scale(&(&Real::from_i64(4) * &torus.major.square()))
            }
        }
    }
    /// Tensor Bernstein coefficients of the cut equation over the rectangle
    /// whose homogeneous Bezier controls are `controls`, of bidegree
    /// `degrees`. Returns the coefficients and their count in each direction.
    pub(super) fn net(
        &self,
        controls: &[[Real; 4]],
        degrees: [usize; 2],
        dec: &mut Decider,
    ) -> Result<(Vec<Real>, [usize; 2]), Refusal> {
        match self {
            Self::Plane(plane) => {
                let point = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
                let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
                let values = controls
                    .iter()
                    .map(|h| {
                        (0..3).fold(Real::zero(), |s, i| {
                            s + &(&h[i] - &(point[i] * &h[3])) * normal[i]
                        })
                    })
                    .collect();
                Ok((values, degrees.map(|d| d + 1)))
            }
            Self::Cylinder(_) | Self::Sphere(_) | Self::Cone(_) => {
                quadric_net(&self.quadratic_form(), controls, degrees, dec)
            }
            Self::Torus(torus) => torus_net(torus, controls, degrees, dec),
        }
    }
}

/// Bernstein product of a quadratic form of four Bernstein nets. For nets of
/// bidegree `(m, n)` the product has bidegree `(2m, 2n)` with coefficient
/// `sum C(m,i1) C(m,i2) C(n,j1) C(n,j2) B(h1, h2) / (C(2m,i) C(2n,j))`
/// over `i1 + i2 = i` and `j1 + j2 = j`, where `B` is the symmetric form.
fn quadric_net(
    form: &[([usize; 2], Real)],
    controls: &[[Real; 4]],
    degrees: [usize; 2],
    dec: &mut Decider,
) -> Result<(Vec<Real>, [usize; 2]), Refusal> {
    let [m, n] = degrees;
    let count = [2 * m + 1, 2 * n + 1];
    let half = Real::from_rational(Rational::from_ratio(1, 2).expect("nonzero"));
    // Symmetric bilinear coefficients: off-diagonal monomials split in half.
    let mut bilinear: Vec<([usize; 2], Real)> = Vec::new();
    for ([i, j], c) in form {
        if i == j {
            bilinear.push(([*i, *j], c.clone()));
        } else {
            let c = exact(&half, c, true, dec)?;
            bilinear.push(([*i, *j], c.clone()));
            bilinear.push(([*j, *i], c));
        }
    }
    let binomial = |n: usize, k: usize| {
        let mut value = Rational::one();
        for i in 0..k {
            value = value * Rational::from_i64((n - i) as i64) / Rational::from_i64((i + 1) as i64);
        }
        value
    };
    let mut values = vec![Real::zero(); count[0] * count[1]];
    for i1 in 0..=m {
        for j1 in 0..=n {
            let a = &controls[i1 * (n + 1) + j1];
            for i2 in 0..=m {
                for j2 in 0..=n {
                    let b = &controls[i2 * (n + 1) + j2];
                    let mut pair = Real::zero();
                    for ([i, j], c) in &bilinear {
                        let term = exact(&a[*i], &b[*j], true, dec)?;
                        let term = exact(c, &term, true, dec)?;
                        pair = exact(&pair, &term, false, dec)?;
                    }
                    let weight =
                        binomial(m, i1) * binomial(m, i2) * binomial(n, j1) * binomial(n, j2)
                            / (binomial(2 * m, i1 + i2) * binomial(2 * n, j1 + j2));
                    let term = exact(&Real::from_rational(weight), &pair, true, dec)?;
                    let k = (i1 + i2) * count[1] + j1 + j2;
                    values[k] = exact(&values[k], &term, false, dec)?;
                }
            }
        }
    }
    Ok((values, count))
}

/// One reserved exact operation. Rational results become fresh leaves so the
/// net does not retain a deep expression of rational arithmetic.
fn exact(a: &Real, b: &Real, multiply: bool, dec: &mut Decider) -> Result<Real, Refusal> {
    if let (Some(x), Some(y)) = (a.rational_ref(), b.rational_ref()) {
        dec.reserve_exact_arithmetic(x, y, multiply)?;
        return Ok(Real::from_rational(if multiply { x * y } else { x + y }));
    }
    Ok(if multiply { a * b } else { a + b })
}
