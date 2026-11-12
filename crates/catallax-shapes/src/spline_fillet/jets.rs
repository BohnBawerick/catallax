//! The exact rolling-ball definition on interval Taylor jets in `tau`.
//!
//! A jet's coefficients are Taylor derivatives divided by factorial over its
//! whole base ball, so one evaluation over a span of `tau` encloses every
//! derivative there up to fourth order.

use super::*;

/// Working precision of every ball.
pub(super) const BITS: u32 = 64;
const LIMBS: u64 = (BITS as u64).div_ceil(64);
/// Charges of one ball product and one ball sum, as the decider prices them.
const BALL_MUL: u64 = LIMBS * LIMBS * 8;
const BALL_ADD: u64 = LIMBS * 4;

pub(super) fn ball(value: &Rational) -> Ball {
    Ball::from_rational(value, BITS)
}

pub(super) fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

fn zero() -> Ball {
    ball(&Rational::zero())
}

/// Taylor coefficients in `tau` to some order, each a ball over the whole
/// base ball: coefficient `k` encloses the `k`-th derivative over `k!`.
#[derive(Clone, Debug)]
pub(super) struct Jet(pub(super) Vec<Ball>);

impl Jet {
    pub(super) fn constant(value: Ball, order: usize) -> Self {
        let mut c = vec![zero(); order + 1];
        c[0] = value;
        Self(c)
    }
    pub(super) fn linear(value: Ball, slope: Ball, order: usize) -> Self {
        let mut jet = Self::constant(value, order);
        if order > 0 {
            jet.0[1] = slope;
        }
        jet
    }
    pub(super) fn order(&self) -> usize {
        self.0.len() - 1
    }
    /// The first `order + 1` coefficients.
    pub(super) fn truncated(&self, order: usize) -> Self {
        Self(self.0[..=order.min(self.order())].to_vec())
    }
    /// Value and first derivative.
    pub(super) fn data(&self) -> [Ball; 2] {
        [
            self.0[0].clone(),
            self.0.get(1).cloned().unwrap_or_else(zero),
        ]
    }
}

/// Jet arithmetic, each ball operation counted; the caller charges the count.
#[derive(Default)]
pub(super) struct Meter(u64);

fn unresolved(detail: &str) -> Refusal {
    Refusal::new(SPLINE_FILLET_SPINE, detail.to_string())
}

impl Meter {
    pub(super) fn charge(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        let units = std::mem::take(&mut self.0);
        dec.charge_work(units)
    }
    pub(super) fn balls(&mut self, count: u64, mul: bool) {
        self.0 += count * if mul { BALL_MUL } else { BALL_ADD };
    }
    pub(super) fn mul(&mut self, a: &Jet, c: &Jet) -> Jet {
        let order = a.order().min(c.order());
        let pairs = ((order + 1) * (order + 2) / 2) as u64;
        self.balls(pairs, true);
        self.balls(pairs, false);
        Jet((0..=order)
            .map(|n| (0..=n).fold(zero(), |sum, i| sum + &a.0[i] * &c.0[n - i]))
            .collect())
    }
    pub(super) fn add(&mut self, a: &Jet, c: &Jet) -> Jet {
        let order = a.order().min(c.order());
        self.balls(order as u64 + 1, false);
        Jet((0..=order).map(|n| &a.0[n] + &c.0[n]).collect())
    }
    pub(super) fn sub(&mut self, a: &Jet, c: &Jet) -> Jet {
        let order = a.order().min(c.order());
        self.balls(order as u64 + 1, false);
        Jet((0..=order).map(|n| &a.0[n] - &c.0[n]).collect())
    }
    pub(super) fn scale(&mut self, a: &Jet, k: &Ball) -> Jet {
        self.balls(a.order() as u64 + 1, true);
        Jet(a.0.iter().map(|x| x * k).collect())
    }
    fn reciprocal(&mut self, c: &Jet) -> Result<Jet, Refusal> {
        let order = c.order();
        let mut r = vec![zero(); order + 1];
        r[0] = ball(&Rational::one())
            .checked_div(&c.0[0])
            .ok_or_else(|| unresolved("a divisor is not separated from zero"))?;
        for n in 1..=order {
            let sum = (1..=n).fold(zero(), |s, i| s + &c.0[i] * &r[n - i]);
            r[n] = (-sum)
                .checked_div(&c.0[0])
                .ok_or_else(|| unresolved("a divisor is not separated from zero"))?;
        }
        let pairs = ((order + 1) * (order + 2) / 2) as u64;
        self.balls(pairs + order as u64 + 1, true);
        Ok(Jet(r))
    }
    pub(super) fn div(&mut self, a: &Jet, c: &Jet) -> Result<Jet, Refusal> {
        let r = self.reciprocal(c)?;
        Ok(self.mul(a, &r))
    }
    pub(super) fn sqrt(&mut self, a: &Jet) -> Result<Jet, Refusal> {
        let order = a.order();
        if a.0[0].lo().sign() != Sign::Positive {
            return Err(unresolved(
                "a square root's argument is not separated from zero",
            ));
        }
        let mut r = vec![zero(); order + 1];
        r[0] = a.0[0]
            .sqrt()
            .ok_or_else(|| unresolved("an undefined square root"))?;
        let twice = &ball(&Rational::from_i64(2)) * &r[0];
        for n in 1..=order {
            let sum = (1..n).fold(zero(), |s, i| s + &r[i] * &r[n - i]);
            r[n] = (&a.0[n] - &sum)
                .checked_div(&twice)
                .ok_or_else(|| unresolved("a square root's derivative is unresolved"))?;
        }
        let pairs = ((order + 1) * (order + 2) / 2) as u64;
        self.balls(pairs + order as u64 + 2, true);
        Ok(Jet(r))
    }
    pub(super) fn dot(&mut self, a: &[Jet; 3], c: &[Jet; 3]) -> Jet {
        let mut total = self.mul(&a[0], &c[0]);
        for i in 1..3 {
            let term = self.mul(&a[i], &c[i]);
            total = self.add(&total, &term);
        }
        total
    }
}

pub(super) fn constant(value: &Rational, order: usize) -> Jet {
    Jet::constant(ball(value), order)
}

/// A polynomial in two parameters with ball coefficients.
pub(super) struct Poly {
    terms: Vec<(usize, usize, Ball, f64)>,
    pub(super) degree: [usize; 2],
}

impl Poly {
    pub(super) fn new(poly: &MPoly) -> Self {
        let mut degree = [0, 0];
        let terms = poly
            .terms()
            .into_iter()
            .map(|(monomial, coefficient)| {
                let [i, j] = [
                    monomial.exponents()[0] as usize,
                    monomial.exponents()[1] as usize,
                ];
                degree = [degree[0].max(i), degree[1].max(j)];
                (i, j, ball(&coefficient), coefficient.to_f64_nearest())
            })
            .collect();
        Self { terms, degree }
    }

    pub(super) fn at(&self, powers: &[Vec<Jet>; 2], meter: &mut Meter) -> Jet {
        let order = powers[0][0].order();
        let mut total = Jet::constant(zero(), order);
        for (i, j, coefficient, _) in &self.terms {
            let monomial = meter.mul(&powers[0][*i], &powers[1][*j]);
            let term = meter.scale(&monomial, coefficient);
            total = meter.add(&total, &term);
        }
        total
    }

    pub(super) fn f64(&self, at: [f64; 2]) -> f64 {
        self.terms
            .iter()
            .map(|(i, j, _, c)| c * at[0].powi(*i as i32) * at[1].powi(*j as i32))
            .sum()
    }
}

/// Powers of two jets up to the given degrees.
pub(super) fn powers(at: &[Jet; 2], degree: [usize; 2], meter: &mut Meter) -> [Vec<Jet>; 2] {
    let order = at[0].order().min(at[1].order());
    std::array::from_fn(|axis| {
        let mut out = vec![Jet::constant(ball(&Rational::one()), order)];
        for k in 0..degree[axis] {
            let next = meter.mul(&out[k], &at[axis]);
            out.push(next);
        }
        out
    })
}

/// One support: its homogeneous coordinates and normal `N`, a positive
/// multiple of `P_u x P_v`, with the side of the ball's centre.
pub(super) struct Support {
    pub(super) exact: [MPoly; 4],
    h: [Poly; 4],
    normal: [Poly; 3],
    degree: [usize; 2],
    /// `+1` when the centre is `P + r N/|N|`, `-1` for `P - r N/|N|`.
    pub(super) side: Rational,
    pub(super) domain: [[Rational; 2]; 2],
}

fn size(poly: &MPoly) -> u64 {
    poly.terms().iter().map(|(_, c)| limbs(c)).sum()
}

fn limbs(value: &Rational) -> u64 {
    value.numerator().bit_length().div_ceil(64).max(1)
        + value.denominator().bit_length().div_ceil(64).max(1)
}

fn product(a: &MPoly, c: &MPoly, dec: &mut Decider) -> Result<MPoly, Refusal> {
    dec.charge_work(size(a).saturating_mul(size(c)))?;
    Ok(a * c)
}

fn difference(a: &MPoly, c: &MPoly, dec: &mut Decider) -> Result<MPoly, Refusal> {
    dec.charge_work(size(a).saturating_add(size(c)))?;
    Ok(a - c)
}

impl Support {
    pub(super) fn new(
        patch: &RationalPatch,
        side: Sign,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let spans = patch.rational_span_polynomials().ok_or_else(|| {
            definition("a spline fillet needs supports with rational coefficients")
        })?;
        let [(domain, h)] = <[_; 1]>::try_from(spans).map_err(|_| {
            Refusal::new(
                SUPPORT_SPANS,
                "a spline fillet on a support of more than one knot span is not covered yet",
            )
        })?;
        let mut tangents: Vec<[MPoly; 3]> = Vec::with_capacity(2);
        for axis in 0..2 {
            let dw = h[3].partial(axis);
            let mut parts = Vec::with_capacity(3);
            for coordinate in &h[..3] {
                let left = product(&coordinate.partial(axis), &h[3], dec)?;
                let right = product(coordinate, &dw, dec)?;
                parts.push(difference(&left, &right, dec)?);
            }
            tangents.push(parts.try_into().expect("three components"));
        }
        let mut normal = Vec::with_capacity(3);
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let left = product(&tangents[0][j], &tangents[1][k], dec)?;
            let right = product(&tangents[0][k], &tangents[1][j], dec)?;
            normal.push(difference(&left, &right, dec)?);
        }
        let normal: [MPoly; 3] = normal.try_into().expect("three components");
        let h_polys = h.each_ref().map(Poly::new);
        let normal_polys = normal.each_ref().map(Poly::new);
        let degree = h_polys.iter().chain(&normal_polys).fold([0, 0], |d, p| {
            [d[0].max(p.degree[0]), d[1].max(p.degree[1])]
        });
        Ok(Self {
            exact: h,
            h: h_polys,
            normal: normal_polys,
            degree,
            side: Rational::from_i64(match side {
                Sign::Positive => 1,
                Sign::Negative => -1,
                Sign::Zero => return Err(definition("a ball side must be positive or negative")),
            }),
            domain,
        })
    }

    /// The point, unit normal and homogeneous weight at parameter jets.
    pub(super) fn contact(&self, at: &[Jet; 2], meter: &mut Meter) -> Result<Contact, Refusal> {
        let powers = powers(at, self.degree, meter);
        let h: [Jet; 4] = std::array::from_fn(|i| self.h[i].at(&powers, meter));
        let normal: [Jet; 3] = std::array::from_fn(|i| self.normal[i].at(&powers, meter));
        let mut point = Vec::with_capacity(3);
        for coordinate in &h[..3] {
            point.push(meter.div(coordinate, &h[3])?);
        }
        let square = meter.dot(&normal, &normal);
        let length = meter.sqrt(&square)?;
        let mut unit = Vec::with_capacity(3);
        for component in &normal {
            unit.push(meter.div(component, &length)?);
        }
        Ok(Contact {
            point: point.try_into().expect("three"),
            unit: unit.try_into().expect("three"),
            weight: h[3].clone(),
        })
    }

    /// The point alone, for exclusion: fewer operations.
    pub(super) fn point(&self, at: &[Jet; 2], meter: &mut Meter) -> Result<[Jet; 3], Refusal> {
        let powers = powers(at, self.degree, meter);
        let h: [Jet; 4] = std::array::from_fn(|i| self.h[i].at(&powers, meter));
        let mut point = Vec::with_capacity(3);
        for coordinate in &h[..3] {
            point.push(meter.div(coordinate, &h[3])?);
        }
        Ok(point.try_into().expect("three"))
    }

    /// The centre `P + side r N/|N|` in binary64, for proposals only.
    fn centre_f64(&self, at: [f64; 2], radius: f64) -> [f64; 3] {
        let h: [f64; 4] = std::array::from_fn(|i| self.h[i].f64(at));
        let n: [f64; 3] = std::array::from_fn(|i| self.normal[i].f64(at));
        let length = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let side = self.side.to_f64_nearest();
        std::array::from_fn(|i| h[i] / h[3] + side * radius * n[i] / length)
    }
}

/// A support's point, unit normal and weight at some parameters.
#[derive(Clone)]
pub(super) struct Contact {
    pub(super) point: [Jet; 3],
    pub(super) unit: [Jet; 3],
    pub(super) weight: Jet,
}

/// The exact definition: both supports, the radius and the parameter `tau`.
pub(super) struct Exact {
    pub(super) supports: [Support; 2],
    pub(super) radius: Rational,
    pub(super) axis: usize,
}

/// The spine over a ball of `tau`: jets of the first support's other
/// parameter and both of the second's, with the proved box round each.
#[derive(Clone)]
pub(super) struct Solution {
    pub(super) x: [Jet; 3],
    pub(super) tau: Jet,
}

/// Everything the face reads at a spine point: both contacts, the centre,
/// the section's middle point and weights, and `cos(phi/2)`.
pub(super) struct Rows {
    pub(super) contacts: [Contact; 2],
    pub(super) centre: [Jet; 3],
    pub(super) middle: [Jet; 3],
    /// The middle weight `cos(phi/2) sqrt(W_A W_B)`.
    pub(super) weight: Jet,
    /// `cos(phi/2)`, from the contact directions.
    pub(super) half: Jet,
}

fn matrix_f64(a: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det == 0.0 || !det.is_finite() {
        return None;
    }
    Some(std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let (r0, r1) = ((j + 1) % 3, (j + 2) % 3);
            let (c0, c1) = ((i + 1) % 3, (i + 2) % 3);
            (a[r0][c0] * a[r1][c1] - a[r0][c1] * a[r1][c0]) / det
        })
    }))
}

impl Exact {
    /// The exact definition's supports, radius and spine parameter.
    pub(super) fn new(
        definition: &SplineFilletDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Ok(Self {
            supports: [
                Support::new(&definition.supports[0], definition.sides[0], dec)?,
                Support::new(&definition.supports[1], definition.sides[1], dec)?,
            ],
            radius: definition.radius.clone(),
            axis: definition.axis,
        })
    }

    /// The two supports' parameters from the unknowns and `tau`.
    pub(super) fn params<T: Clone>(&self, x: &[T; 3], tau: &T) -> [[T; 2]; 2] {
        let first = if self.axis == 0 {
            [tau.clone(), x[0].clone()]
        } else {
            [x[0].clone(), tau.clone()]
        };
        [first, [x[1].clone(), x[2].clone()]]
    }

    fn centre(&self, contact: &Contact, side: &Rational, meter: &mut Meter) -> [Jet; 3] {
        let reach = ball(&(&self.radius * side));
        std::array::from_fn(|i| {
            let along = meter.scale(&contact.unit[i], &reach);
            meter.add(&contact.point[i], &along)
        })
    }

    /// The centre difference `C_A - C_B`, zero on the spine.
    pub(super) fn residual(
        &self,
        x: &[Jet; 3],
        tau: &Jet,
        meter: &mut Meter,
    ) -> Result<[Jet; 3], Refusal> {
        let [pa, pb] = self.params(x, tau);
        let [a, b] = &self.supports;
        let ca = self.centre(&a.contact(&pa, meter)?, &a.side, meter);
        let cb = self.centre(&b.contact(&pb, meter)?, &b.side, meter);
        Ok(std::array::from_fn(|i| meter.sub(&ca[i], &cb[i])))
    }

    fn residual_f64(&self, x: [f64; 3], tau: f64) -> [f64; 3] {
        let [pa, pb] = self.params(&x, &tau);
        let r = self.radius.to_f64_nearest();
        let ca = self.supports[0].centre_f64(pa, r);
        let cb = self.supports[1].centre_f64(pb, r);
        std::array::from_fn(|i| ca[i] - cb[i])
    }

    /// A binary64 Newton estimate of the spine point at `tau`.
    pub(super) fn newton(&self, tau: f64, mut x: [f64; 3]) -> Option<[f64; 3]> {
        for _ in 0..40 {
            let f = self.residual_f64(x, tau);
            let mut jacobian = [[0.0; 3]; 3];
            for j in 0..3 {
                let h = 1e-7 * x[j].abs().max(1.0);
                let mut up = x;
                let mut down = x;
                up[j] += h;
                down[j] -= h;
                let (fu, fd) = (self.residual_f64(up, tau), self.residual_f64(down, tau));
                for i in 0..3 {
                    jacobian[i][j] = (fu[i] - fd[i]) / (2.0 * h);
                }
            }
            let inverse = matrix_f64(jacobian)?;
            let mut size = 0.0f64;
            for i in 0..3 {
                let step: f64 = (0..3).map(|j| inverse[i][j] * f[j]).sum();
                x[i] -= step;
                size = size.max(step.abs());
            }
            if !x.iter().all(|v| v.is_finite()) {
                return None;
            }
            if size < 1e-15 {
                break;
            }
        }
        Some(x)
    }

    /// The residual and its Jacobian in the unknowns over a box and a ball
    /// of `tau`: the value at the box's centre and the partials over the box.
    fn linear(
        &self,
        centre: &[Ball; 3],
        boxes: &[Ball; 3],
        tau: &Ball,
        meter: &mut Meter,
    ) -> Result<([Ball; 3], [[Ball; 3]; 3]), Refusal> {
        let at: [Jet; 3] = std::array::from_fn(|i| Jet::constant(centre[i].clone(), 0));
        let value = self.residual(&at, &Jet::constant(tau.clone(), 0), meter)?;
        let t = Jet::constant(tau.clone(), 1);
        let mut columns = Vec::with_capacity(3);
        for j in 0..3 {
            let x: [Jet; 3] = std::array::from_fn(|i| {
                if i == j {
                    Jet::linear(boxes[i].clone(), ball(&Rational::one()), 1)
                } else {
                    Jet::constant(boxes[i].clone(), 1)
                }
            });
            let g = self.residual(&x, &t, meter)?;
            columns.push(g.map(|jet| jet.0[1].clone()));
        }
        let jacobian = std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i].clone()));
        Ok((value.map(|jet| jet.0[0].clone()), jacobian))
    }

    /// Proves a unique spine point in a box for every `tau` in `span`, then
    /// encloses its Taylor coefficients to `order` over the span. The box
    /// comes from binary64 estimates at both ends.
    pub(super) fn solve(
        &self,
        span: [&Rational; 2],
        estimates: [[f64; 3]; 2],
        order: usize,
        dec: &mut Decider,
    ) -> Result<Solution, Refusal> {
        let mut meter = Meter::default();
        let tau = ball(span[0]).hull(&ball(span[1]));
        // Proposed box: the hull of both estimates, widened by a quarter of
        // its width and a little more.
        let mut boxes: [Ball; 3] = std::array::from_fn(|i| {
            let (a, c) = (estimates[0][i], estimates[1][i]);
            let width = (a - c).abs();
            let pad = width / 4.0 + 1e-9 * a.abs().max(1.0);
            let lo = Rational::from_f64_exact(a.min(c) - pad).unwrap_or_else(Rational::zero);
            let hi = Rational::from_f64_exact(a.max(c) + pad).unwrap_or_else(Rational::zero);
            ball(&lo).hull(&ball(&hi))
        });
        let mut proved = None;
        for _ in 0..8 {
            let centre: [Ball; 3] = std::array::from_fn(|i| ball(&boxes[i].mid()));
            let (value, jacobian) = self.linear(&centre, &boxes, &tau, &mut meter)?;
            meter.charge(dec)?;
            let middle: [[f64; 3]; 3] = std::array::from_fn(|i| {
                std::array::from_fn(|j| jacobian[i][j].mid().to_f64_nearest())
            });
            let inverse =
                matrix_f64(middle).ok_or_else(|| spine("the spine Jacobian is singular"))?;
            let y: [[Ball; 3]; 3] = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    ball(&Rational::from_f64_exact(inverse[i][j]).unwrap_or_else(Rational::zero))
                })
            });
            // I - Y J, and the Krawczyk image m - Y g(m) + (I - Y J)(X - m).
            let mut residue: [[Ball; 3]; 3] =
                std::array::from_fn(|_| std::array::from_fn(|_| ball(&Rational::zero())));
            for i in 0..3 {
                for j in 0..3 {
                    let mut total = ball(&Rational::from_i64(i64::from(i == j)));
                    for k in 0..3 {
                        total = &total - &(&y[i][k] * &jacobian[k][j]);
                    }
                    residue[i][j] = total;
                }
            }
            meter.balls(81, true);
            let mut image: Vec<Ball> = Vec::with_capacity(3);
            for i in 0..3 {
                let mut total = centre[i].clone();
                for k in 0..3 {
                    total = &total - &(&y[i][k] * &value[k]);
                    total = &total + &(&residue[i][k] * &(&boxes[k] - &centre[k]));
                }
                image.push(total);
            }
            meter.balls(36, true);
            meter.charge(dec)?;
            let inside =
                (0..3).all(|i| image[i].lo() > boxes[i].lo() && image[i].hi() < boxes[i].hi());
            let rho = (0..3)
                .map(|i| (0..3).fold(Rational::zero(), |s, j| s + magnitude(&residue[i][j])))
                .max()
                .expect("three rows");
            if inside && rho < Rational::one() {
                proved = Some((boxes.clone(), y, residue, rho));
                break;
            }
            // Epsilon inflation round the image.
            boxes = std::array::from_fn(|i| {
                let hull = image[i].hull(&boxes[i]);
                let width = hull.hi() - hull.lo();
                let pad = &width / &Rational::from_i64(8) + Rational::one().mul_pow2(-60);
                ball(&(hull.lo() - pad.clone())).hull(&ball(&(hull.hi() + pad)))
            });
        }
        let Some((boxes, y, residue, rho)) = proved else {
            return Err(spine(format!(
                "no unique spine point was proved for tau in [{}, {}]",
                span[0], span[1]
            )));
        };
        // Higher Taylor coefficients: J x_k = -R_k, with R_k the k-th
        // coefficient of the residual while x_k is zero. With the
        // preconditioner Y and |I - Y J| <= rho < 1 over the box,
        // x_k lies in -Y R_k + (I - Y J) [-beta, beta], beta = |Y R_k| / (1 - rho).
        let mut x: [Jet; 3] = std::array::from_fn(|i| Jet::constant(boxes[i].clone(), order));
        let one_minus = Rational::one() - rho;
        for k in 1..=order {
            let t = Jet::linear(tau.clone(), ball(&Rational::one()), k);
            let at: [Jet; 3] = std::array::from_fn(|i| x[i].truncated(k));
            let g = self.residual(&at, &t, &mut meter)?;
            meter.charge(dec)?;
            let r: [Ball; 3] = std::array::from_fn(|i| g[i].0[k].clone());
            let yr: [Ball; 3] = std::array::from_fn(|i| {
                (0..3).fold(ball(&Rational::zero()), |s, j| &s + &(&y[i][j] * &r[j]))
            });
            let beta = yr.iter().map(magnitude).max().expect("three") / one_minus.clone();
            let spread = ball(&-beta.clone()).hull(&ball(&beta));
            for (i, xi) in x.iter_mut().enumerate() {
                xi.0[k] = residue[i]
                    .iter()
                    .fold(-&yr[i], |value, entry| &value + &(entry * &spread));
            }
            meter.balls(30, true);
        }
        meter.charge(dec)?;
        let tau = Jet::linear(tau, ball(&Rational::one()), order);
        Ok(Solution { x, tau })
    }

    /// The face's rows at a solved spine point.
    pub(super) fn rows(&self, solution: &Solution, meter: &mut Meter) -> Result<Rows, Refusal> {
        let [pa, pb] = self.params(&solution.x, &solution.tau);
        let [a, b] = &self.supports;
        let ca = a.contact(&pa, meter)?;
        let cb = b.contact(&pb, meter)?;
        let centre = self.centre(&ca, &a.side, meter);
        // Unit directions from the centre to each contact: -side N/|N|.
        let toward = |contact: &Contact, side: &Rational, meter: &mut Meter| -> [Jet; 3] {
            let s = ball(&-side.clone());
            std::array::from_fn(|i| meter.scale(&contact.unit[i], &s))
        };
        let da = toward(&ca, &a.side, meter);
        let db = toward(&cb, &b.side, meter);
        let cosine = meter.dot(&da, &db);
        let one = constant(&Rational::one(), cosine.order());
        let plus = meter.add(&one, &cosine);
        let halved = meter.scale(&plus, &ball(&Rational::from_ratio(1, 2).expect("nonzero")));
        let half = meter.sqrt(&halved)?;
        let r = ball(&self.radius);
        let mut middle = Vec::with_capacity(3);
        for i in 0..3 {
            let sum = meter.add(&da[i], &db[i]);
            let scaled = meter.scale(&sum, &r);
            let along = meter.div(&scaled, &plus)?;
            middle.push(meter.add(&centre[i], &along));
        }
        let weights = meter.mul(&ca.weight, &cb.weight);
        let root = meter.sqrt(&weights)?;
        let weight = meter.mul(&half, &root);
        Ok(Rows {
            contacts: [ca, cb],
            centre,
            middle: middle.try_into().expect("three"),
            weight,
            half,
        })
    }
}
