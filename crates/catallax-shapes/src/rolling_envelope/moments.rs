//! Material added by the junction blend, on the exact definition.
//!
//! With the ball centre C(s), unit contact directions a (to the cylinder) and c
//! (to the mouth edge), and b the unit vector of the contact plane orthogonal to
//! a, the envelope is S = C + r (cos p a + sin p b) for p in [0, alpha]. The
//! contact plane is orthogonal to C', so with T = C'/|C'| the area element is
//! A r dp ds for A = |C'| + r (cos p a'.T + sin p b'.T). Above the mouth plane
//! the sector's meridian walls and the cylinder have no axial normal, so the
//! divergence theorem turns each moment into one integral over S. Each integral
//! over p is a closed form, leaving one over s, enclosed panel by panel. About a
//! panel's midpoint with half-width w, the Taylor terms give 2 w c0 + 2/3 w^3 c2.
//! The mean of c3 cancels by symmetry, so the third-order tail adds only the
//! width of c3 over the panel times w^4/2.

use std::{cell::RefCell, collections::BTreeMap};

use super::*;

/// What one mouth sector adds relative to the unfilleted body, enclosed on the
/// exact definition, independent of the fitted carrier.
#[derive(Clone, Debug)]
pub struct ContactMoments {
    /// Half the sector's angle at the cylinder axis.
    pub half_angle: Ball,
    /// Added material, in mm^3.
    pub volume: Ball,
    /// First moment of that material along the axis, from the mouth plane.
    pub axial_moment: Ball,
    /// First moment along `toward`, from the cylinder axis. Its moment along
    /// `across` is zero by the definition's mirror symmetry.
    pub radial_moment: Ball,
    /// Area of the exact envelope.
    pub envelope_area: Ball,
    /// Plane area the envelope covers.
    pub plane_area: Ball,
    /// Cylinder area the envelope covers.
    pub cylinder_area: Ball,
    /// Upper bound on the integral of the envelope's non-tube principal
    /// curvature magnitude over its area.
    pub turning: Rational,
}

/// Per-radian constants of the untouched circular section between the plane,
/// the cylinder and the ball.
#[derive(Clone, Debug)]
pub struct JunctionSection {
    /// Added material.
    pub volume: Ball,
    /// First moment along the axis, from the plane.
    pub axial_moment: Ball,
    /// First moment along the radial direction, from the axis.
    pub radial_moment: Ball,
    /// Added torus area minus the plane and cylinder area it covers.
    pub area_change: Ball,
}

// A trigonometric polynomial in p with jet coefficients: cos^i sin^j terms.
struct Trig(BTreeMap<(usize, usize), Jet>);

impl Trig {
    fn term(i: usize, j: usize, value: Jet) -> Self {
        Self(BTreeMap::from([((i, j), value)]))
    }
    fn add(&mut self, key: (usize, usize), value: Jet) {
        match self.0.get_mut(&key) {
            Some(sum) => *sum = sum.add(&value),
            None => {
                self.0.insert(key, value);
            }
        }
    }
    fn mul(&self, other: &Self) -> Self {
        let mut out = Self(BTreeMap::new());
        for ((i, j), a) in &self.0 {
            for ((k, l), c) in &other.0 {
                out.add((i + k, j + l), a.mul(c));
            }
        }
        out
    }
    fn plus(mut self, other: Self) -> Self {
        for (key, value) in other.0 {
            self.add(key, value);
        }
        self
    }
}

// Integrals of cos^i sin^j over [0, alpha], by the reduction formulas.
struct Powers<'a> {
    alpha: &'a Jet,
    sin: &'a Jet,
    cos: &'a Jet,
    one: Jet,
    known: RefCell<BTreeMap<(usize, usize), Jet>>,
}

impl Powers<'_> {
    fn power(base: &Jet, exponent: usize, one: &Jet) -> Jet {
        (0..exponent).fold(one.clone(), |product, _| product.mul(base))
    }
    fn integral(&self, i: usize, j: usize, k: &Constants) -> Jet {
        if let Some(known) = self.known.borrow().get(&(i, j)) {
            return known.clone();
        }
        let value = self.reduced(i, j, k);
        self.known.borrow_mut().insert((i, j), value.clone());
        value
    }
    fn reduced(&self, i: usize, j: usize, k: &Constants) -> Jet {
        let bits = k.one.precision();
        let ratio = |a: usize, c: usize| {
            b(
                &Rational::from_ratio(a as i64, c as i64).expect("nonzero"),
                bits,
            )
        };
        match (i, j) {
            (0, 0) => self.alpha.clone(),
            (1, 0) => self.sin.clone(),
            (0, 1) => self.one.sub(self.cos),
            (1, 1) => self.sin.square().scale(&ratio(1, 2)),
            (_, j) if j >= 2 => {
                // int c^i s^j = (-s^(j-1) c^(i+1) + (j-1) int c^i s^(j-2)) / (i+j)
                let boundary = Self::power(self.sin, j - 1, &self.one)
                    .mul(&Self::power(self.cos, i + 1, &self.one))
                    .neg();
                boundary
                    .add(&self.integral(i, j - 2, k).scale(&ratio(j - 1, 1)))
                    .scale(&ratio(1, i + j))
            }
            (i, j) => {
                // int c^i s^j = (c^(i-1) s^(j+1) + (i-1) int c^(i-2) s^j) / (i+j)
                let boundary = Self::power(self.cos, i - 1, &self.one).mul(&Self::power(
                    self.sin,
                    j + 1,
                    &self.one,
                ));
                boundary
                    .add(&self.integral(i - 2, j, k).scale(&ratio(i - 1, 1)))
                    .scale(&ratio(1, i + j))
            }
        }
    }
    fn total(&self, trig: &Trig, k: &Constants) -> Jet {
        trig.0.iter().fold(
            Jet::constant(b(&Rational::zero(), k.one.precision())),
            |sum, ((i, j), c)| sum.add(&c.mul(&self.integral(*i, *j, k))),
        )
    }
}

fn shifted(f: &Jet) -> Jet {
    let bits = f.0[0].precision();
    Jet(std::array::from_fn(|i| {
        if i == 4 {
            b(&Rational::zero(), bits)
        } else {
            &f.0[i + 1] * &b(&q(i as i64 + 1), bits)
        }
    }))
}
fn dot3(a: &[Jet; 3], c: &[Jet; 3]) -> Jet {
    a[0].mul(&c[0]).add(&a[1].mul(&c[1])).add(&a[2].mul(&c[2]))
}

// Six densities over the chart parameter, valid through their third Taylor
// coefficient: volume, axial and radial moment, envelope, plane, cylinder area.
// With them, a bound on the density of the curvature integral: the curvature is
// (cos p a'.T + sin p b'.T) / A against the area element A r dp ds.
fn densities(
    k: &Constants,
    body: &Ball,
    s: &Jet,
    dec: &Decider,
) -> Result<([Jet; 6], Rational), Refusal> {
    let station = station_jets(k, s)?;
    let zero = Jet::constant(b(&Rational::zero(), k.one.precision()));
    let r = &k.radius;
    let centre = [
        station.cx.clone(),
        station.height.clone(),
        station.cz.clone(),
    ];
    let a = [station.cos.neg(), zero.clone(), station.sin.neg()];
    let r_inv = k
        .one
        .checked_div(r)
        .ok_or_else(|| unresolved("ball radius is not separated from zero"))?;
    let c = [
        station.rx.sub(&centre[0]).scale(&r_inv),
        centre[1].neg().scale(&r_inv),
        station.rz.sub(&centre[2]).scale(&r_inv),
    ];
    let cosine = dot3(&a, &c);
    let sine = Jet::constant(k.one.clone()).sub(&cosine.square()).sqrt()?;
    let inverse = sine.reciprocal()?;
    let frame: [Jet; 3] = std::array::from_fn(|i| c[i].sub(&cosine.mul(&a[i])).mul(&inverse));
    let velocity: [Jet; 3] = std::array::from_fn(|i| shifted(&centre[i]));
    let speed = dot3(&velocity, &velocity).sqrt()?;
    let unit_inverse = speed.reciprocal()?;
    let tangent: [Jet; 3] = std::array::from_fn(|i| velocity[i].mul(&unit_inverse));
    let turn_a = dot3(&a.clone().map(|x| shifted(&x)), &tangent);
    let turn_b = dot3(&frame.clone().map(|x| shifted(&x)), &tangent);
    // cos(alpha) = cosine gives alpha' = -cosine' / sin(alpha); its value needs
    // the arctangent, and higher coefficients integrate that derivative.
    let rate = shifted(&cosine).neg().mul(&inverse);
    let value = Ball::atan2(&sine.0[0], &cosine.0[0])
        .ok_or_else(|| unresolved("contact angle is not separated from zero"))?;
    let alpha = Jet(std::array::from_fn(|i| {
        if i == 0 {
            value.clone()
        } else {
            &rate.0[i - 1]
                * &b(
                    &Rational::from_ratio(1, i as i64).expect("nonzero"),
                    k.one.precision(),
                )
        }
    }));
    let powers = Powers {
        alpha: &alpha,
        sin: &sine,
        cos: &cosine,
        one: Jet::constant(k.one.clone()),
        known: RefCell::new(BTreeMap::new()),
    };
    let drop = frame[1].neg().scale(r);
    let term = Trig::term;
    let area = term(0, 0, speed.clone())
        .plus(term(1, 0, turn_a.scale(r)))
        .plus(term(0, 1, turn_b.scale(r)));
    let height = term(0, 0, station.height.clone()).plus(term(0, 1, drop.neg()));
    let radial = term(0, 0, station.cx.clone())
        .plus(term(1, 0, a[0].scale(r)))
        .plus(term(0, 1, frame[0].scale(r)));
    let lifted = term(0, 1, drop.clone()).mul(&area);
    let half = b(&Rational::from_ratio(1, 2).unwrap(), k.one.precision());
    let volume = powers.total(&lifted.mul(&height), k);
    let axial = powers
        .total(&lifted.mul(&height).mul(&height), k)
        .scale(&half);
    let moment = powers.total(&lifted.mul(&height).mul(&radial), k);
    let envelope = powers.total(&area, k).scale(r);
    let plane = powers.total(&lifted, k);
    // Cylinder strip: body * height * d(theta)/ds, theta = 2 atan(t s).
    let u = s.scale(&k.t);
    let dtheta =
        Jet::constant(&k.two * &k.t).div(&Jet::constant(k.one.clone()).add(&u.square()))?;
    let cylinder = station.height.mul(&dtheta).scale(body);
    dec.check_deadline()?;
    let turning =
        (magnitude(&turn_a.0[0]) + magnitude(&turn_b.0[0])) * magnitude(&alpha.0[0]) * magnitude(r);
    Ok(([volume, axial, moment, envelope, plane, cylinder], turning))
}

impl RollingEnvelope {
    /// The definition's sector moments, shared by every rigid copy of this fit.
    pub fn contact_moments(
        &self,
        width: &Rational,
        dec: &Decider,
    ) -> Result<Arc<ContactMoments>, Refusal> {
        let mut known = self
            .data
            .net
            .moments
            .lock()
            .expect("no panic holds this lock");
        if let Some((key, moments)) = known.as_ref()
            && key == width
        {
            return Ok(moments.clone());
        }
        let moments = Arc::new(self.data.definition.contact_moments(width, dec)?);
        *known = Some((width.clone(), moments.clone()));
        Ok(moments)
    }
}

impl RollingEnvelopeDefinition {
    /// Encloses one sector's moments with volume and area widths at most
    /// `width` and first-moment widths at most eight times that. The integrals
    /// are even in the chart parameter, so only its upper half is summed.
    pub fn contact_moments(
        &self,
        width: &Rational,
        dec: &Decider,
    ) -> Result<ContactMoments, Refusal> {
        let bits = 128.min(dec.policy().max_precision_bits);
        let k = Constants::new(self, bits, dec)?;
        let body = dec.enclose_at(&self.body, bits)?;
        let twelfth = b(&Rational::from_ratio(1, 12).unwrap(), bits);
        let limits: [Rational; 6] = std::array::from_fn(|i| {
            if (1..=2).contains(&i) {
                width.mul_pow2(3)
            } else {
                width.clone()
            }
        });
        let mut panels = 64i64;
        loop {
            dec.check_deadline()?;
            let h = q(1) / q(panels);
            let hb = b(&h, bits);
            let curvature = &(&hb * &(&hb * &hb)) * &twelfth;
            let half = &h * &Rational::from_ratio(1, 2).unwrap();
            let tail =
                &(&(&half * &half) * &(&half * &half)) * &Rational::from_ratio(1, 2).unwrap();
            let mut sums: [Ball; 6] = std::array::from_fn(|_| b(&Rational::zero(), bits));
            let mut turning = Rational::zero();
            for panel in 0..panels {
                dec.check_deadline()?;
                let lo = &h * &q(panel);
                let hi = &lo + &h;
                let middle = (&lo + &hi) * Rational::from_ratio(1, 2).unwrap();
                let point = Jet::linear(b(&middle, bits), k.one.clone());
                let whole = Jet::linear(b(&lo, bits).hull(&b(&hi, bits)), k.one.clone());
                let (at, _) = densities(&k, &body, &point, dec)?;
                let (over, density) = densities(&k, &body, &whole, dec)?;
                turning = turning + (&h * &density).mul_pow2(1);
                for (sum, (mid, full)) in sums.iter_mut().zip(at.iter().zip(&over)) {
                    let extent = &full.0[3].rad() * &tail;
                    let rest = b(&(-&extent), bits).hull(&b(&extent, bits));
                    let panel = &(&(&hb * &mid.0[0]) + &(&curvature * &mid.0[2])) + &rest;
                    *sum = &*sum + &(&k.two * &panel);
                }
            }
            // Each halving of a panel gains a factor of sixteen.
            let mut ratio = sums
                .iter()
                .zip(&limits)
                .map(|(sum, limit)| sum.rad().mul_pow2(1) / limit.clone())
                .max()
                .expect("six moments");
            if ratio <= Rational::one() {
                let t = dec.enclose_at(&self.transition().t, bits)?;
                let half_angle = &k.two * &t.atan_checked(dec)?;
                let [
                    volume,
                    axial_moment,
                    radial_moment,
                    envelope_area,
                    plane_area,
                    cylinder_area,
                ] = sums;
                return Ok(ContactMoments {
                    half_angle,
                    volume,
                    axial_moment,
                    radial_moment,
                    envelope_area,
                    plane_area,
                    cylinder_area,
                    turning,
                });
            }
            while ratio > Rational::one() {
                ratio = ratio.mul_pow2(-4);
                panels *= 2;
            }
            if panels > 1 << 16 {
                return Err(unresolved(
                    "envelope moments do not reach the requested width",
                ));
            }
        }
    }

    /// Exact per-radian section constants, for the ring's untouched arcs.
    pub fn section(&self, dec: &Decider) -> Result<JunctionSection, Refusal> {
        let bits = 128.min(dec.policy().max_precision_bits);
        let r = dec.enclose_at(&self.radius, bits)?;
        let body = dec.enclose_at(&self.body, bits)?;
        let big = &r + &body;
        let pi = Ball::pi_checked(bits, dec)?;
        let ratio = |a: i64, c: i64| b(&Rational::from_ratio(a, c).expect("nonzero"), bits);
        let r2 = r.square();
        let r3 = &r2 * &r;
        let r4 = r2.square();
        let quarter = &pi * &ratio(1, 4);
        // With u = R - rho and w = sqrt(r^2 - u^2), the section's height is
        // r - w; int w = pi r^2/4, int u w = r^3/3, int u^2 w = pi r^4/16.
        let volume = &(&big * &(&r2 * &(&ratio(1, 1) - &quarter))) - &(&r3 * &ratio(1, 6));
        let radial_moment = &(&r
            * &(&(&(&big.square() * &r) - &(&big * &r2)) + &(&r3 * &ratio(1, 3))))
            - &(&(&(&big.square() * &(&quarter * &r2)) - &(&(&big * &r3) * &ratio(2, 3)))
                + &(&(&pi * &r4) * &ratio(1, 16)));
        let axial_moment = &(&(&(&(&r2 * &ratio(2, 1)) * &(&(&big * &r) - &(&r2 * &ratio(1, 2))))
            - &(&(&big * &r3) * &ratio(1, 3)))
            + &(&(&r4 * &ratio(1, 4))
                - &(&(&r * &ratio(2, 1))
                    * &(&(&big * &(&quarter * &r2)) - &(&r3 * &ratio(1, 3))))))
            * &ratio(1, 2);
        let torus = &r * &(&(&(&big * &pi) * &ratio(1, 2)) - &r);
        let plane = &(&big.square() - &body.square()) * &ratio(1, 2);
        let cylinder = &body * &r;
        let area_change = &(&torus - &plane) - &cylinder;
        Ok(JunctionSection {
            volume,
            axial_moment,
            radial_moment,
            area_change,
        })
    }
}
