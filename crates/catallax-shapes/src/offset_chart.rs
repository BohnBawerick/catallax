//! Algebraic quarter charts of an elliptic-cylinder normal offset.

use catallax_numbers::{Ball, Decider, EvidenceBound, EvidenceUnit, Rational, Real, Refusal};

use crate::{EllipticCylinder, EllipticOffsetDefinition, OffsetSide, Point3, Surface, Vector3};

fn unresolved(detail: impl Into<String>) -> Refusal {
    Refusal::new("numbers.geometric-gap-unproven", detail)
}
fn b(value: &Rational, bits: u32) -> Ball {
    Ball::from_rational(value, bits)
}
fn zero(bits: u32) -> Ball {
    b(&Rational::zero(), bits)
}

// Coefficients are Taylor derivatives divided by factorial, on the entire
// supplied parameter interval. Fixed order avoids allocating expression jets.
#[derive(Clone)]
struct Jet([Ball; 5]);
impl Jet {
    fn constant(value: Ball) -> Self {
        let bits = value.precision();
        Self(std::array::from_fn(|i| {
            if i == 0 { value.clone() } else { zero(bits) }
        }))
    }
    fn variable(value: Ball) -> Self {
        let bits = value.precision();
        let mut result = Self::constant(value);
        result.0[1] = b(&Rational::one(), bits);
        result
    }
    fn add(&self, other: &Self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] + &other.0[i]))
    }
    fn neg(&self) -> Self {
        Self(std::array::from_fn(|i| -&self.0[i]))
    }
    fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }
    fn mul(&self, other: &Self) -> Self {
        Self(std::array::from_fn(|n| {
            (0..=n).fold(zero(self.0[0].precision()), |sum, i| {
                sum + &self.0[i] * &other.0[n - i]
            })
        }))
    }
    fn square(&self) -> Self {
        let mut result = self.mul(self);
        result.0[0] = self.0[0].square();
        result
    }
    fn reciprocal(&self) -> Result<Self, Refusal> {
        let bits = self.0[0].precision();
        let mut result: [Ball; 5] = std::array::from_fn(|_| zero(bits));
        result[0] = b(&Rational::one(), bits)
            .checked_div(&self.0[0])
            .ok_or_else(|| unresolved("offset chart divisor is not separated from zero"))?;
        for n in 1..5 {
            let sum = (1..=n).fold(zero(bits), |sum, i| sum + &self.0[i] * &result[n - i]);
            result[n] = (-sum)
                .checked_div(&self.0[0])
                .ok_or_else(|| unresolved("offset chart derivative divisor is unresolved"))?;
        }
        Ok(Self(result))
    }
    fn div(&self, other: &Self) -> Result<Self, Refusal> {
        Ok(self.mul(&other.reciprocal()?))
    }
    fn sqrt(&self) -> Result<Self, Refusal> {
        let bits = self.0[0].precision();
        if self.0[0].lo() <= Rational::zero() {
            return Err(unresolved(
                "offset chart requires a positive square-root lower bound",
            ));
        }
        let mut result: [Ball; 5] = std::array::from_fn(|_| zero(bits));
        result[0] = self.0[0]
            .sqrt()
            .ok_or_else(|| unresolved("undefined offset chart square root"))?;
        let denominator = &b(&Rational::from_i64(2), bits) * &result[0];
        for n in 1..5 {
            let sum = (1..n).fold(zero(bits), |sum, i| sum + &result[i] * &result[n - i]);
            result[n] = (&self.0[n] - &sum)
                .checked_div(&denominator)
                .ok_or_else(|| {
                    unresolved("offset square-root derivative lower bound unresolved")
                })?;
        }
        Ok(Self(result))
    }
}

fn original(definition: &EllipticOffsetDefinition) -> &EllipticCylinder {
    match definition.original() {
        Surface::EllipticCylinder(original) => original,
        _ => unreachable!("checked elliptic-offset definition invariant"),
    }
}
fn frame(original: &EllipticCylinder) -> [Vector3; 3] {
    let major = original
        .major_dir
        .scaled(&(&Real::one() / &original.major_dir.norm_squared().sqrt()));
    let axis = original
        .axis
        .scaled(&(&Real::one() / &original.axis.norm_squared().sqrt()));
    let minor = axis.cross(&major);
    [major, minor, axis]
}
fn signed(definition: &EllipticOffsetDefinition) -> Real {
    match definition.side() {
        OffsetSide::Outward => definition.distance().clone(),
        OffsetSide::Inward => -definition.distance(),
    }
}
fn rotate<T>(quarter: usize, x: T, y: T, negate: impl Fn(&T) -> T) -> (T, T) {
    match quarter {
        0 => (x, y),
        1 => (negate(&y), x),
        2 => (negate(&x), negate(&y)),
        3 => (y, negate(&x)),
        _ => unreachable!("quarter range checked by caller"),
    }
}

/// Compiled exact quarter charts shared by one construction or verification pass.
/// Frames, semi-axes and the algebraic chart scale are built only once per pass.
pub struct EllipticOffsetChart<'a> {
    definition: &'a EllipticOffsetDefinition,
    axes: [Vector3; 3],
    semiaxes: [Real; 2],
    scale: Real,
    distance: Real,
}

impl EllipticOffsetDefinition {
    /// Compiles the exact algebraic chart constants without fitting geometry.
    pub fn chart(&self) -> EllipticOffsetChart<'_> {
        let original = original(self);
        EllipticOffsetChart {
            definition: self,
            axes: frame(original),
            semiaxes: [original.major_sq.sqrt(), original.minor_sq.sqrt()],
            scale: &Real::from_i64(2).sqrt() - &Real::one(),
            distance: signed(self),
        }
    }
}

impl EllipticOffsetChart<'_> {
    /// Exact unit direction of the cylinder's linear sweep.
    pub fn axial_direction(&self) -> &Vector3 {
        &self.axes[2]
    }

    /// Exact point and derivative in a bounded algebraic quarter chart.
    /// Each quarter uses tan(theta/2)=(sqrt(2)-1)*(2t-1), 0 <= t <= 1.
    /// Adjacent quarters have equal endpoint parameter speed, including closure.
    /// An inward chart is not itself proof of a valid erosion branch.
    pub fn point_and_tangent(
        &self,
        quarter: usize,
        parameter: &Rational,
        axial: &Real,
        dec: &mut Decider,
    ) -> Result<(Point3, Vector3), Refusal> {
        dec.check_deadline()?;
        crate::math::validate(axial, dec)?;
        if quarter >= 4 || parameter < &Rational::zero() || parameter > &Rational::one() {
            return Err(unresolved(
                "offset quarter chart parameter is outside its finite domain",
            ));
        }
        let original = original(self.definition);
        let r = &self.scale;
        let t = r * &Real::from_rational(parameter * &Rational::from_i64(2) - Rational::one());
        let denominator = &Real::one() + &t.square();
        let c = &(&Real::one() - &t.square()) / &denominator;
        let s = &(&Real::from_i64(2) * &t) / &denominator;
        let dc = &(&(&Real::from_i64(-8) * r) * &t) / &denominator.square();
        let ds =
            &(&(&Real::from_i64(4) * r) * &(&Real::one() - &t.square())) / &denominator.square();
        let (c, s) = rotate(quarter, c, s, |x| -x);
        let (dc, ds) = rotate(quarter, dc, ds, |x| -x);
        let [a, b] = &self.semiaxes;
        let rho_sq = &(&original.minor_sq * &c.square()) + &(&original.major_sq * &s.square());
        let rho = rho_sq.sqrt();
        let product = &(&(&original.minor_sq * &c) * &dc) + &(&(&original.major_sq * &s) * &ds);
        let distance = &self.distance;
        let x = &(a * &c) + &(&(distance * &(b * &c)) / &rho);
        let y = &(b * &s) + &(&(distance * &(a * &s)) / &rho);
        let nx = &(b * &dc) / &rho - &(&(b * &c) * &product) / &(&rho_sq * &rho);
        let ny = &(a * &ds) / &rho - &(&(a * &s) * &product) / &(&rho_sq * &rho);
        let dx = &(a * &dc) + &(distance * &nx);
        let dy = &(b * &ds) + &(distance * &ny);
        let axes = &self.axes;
        let point = original.origin.plus(
            &axes[0]
                .scaled(&x)
                .plus(&axes[1].scaled(&y))
                .plus(&axes[2].scaled(axial)),
        );
        let tangent = axes[0].scaled(&dx).plus(&axes[1].scaled(&dy));
        Ok((point, tangent))
    }

    /// Cubic Hermite interpolation's uniform position-error upper bound over
    /// one quarter-chart interval, in millimetres. It uses interval fourth
    /// derivatives over the full interval, not point samples or fit residuals.
    /// A consumer must separately verify endpoint interpolation and full coverage.
    pub fn hermite_error_bound(
        &self,
        quarter: usize,
        interval: [&Rational; 2],
        dec: &Decider,
    ) -> Result<EvidenceBound, Refusal> {
        dec.check_deadline()?;
        if quarter >= 4
            || interval[0] < &Rational::zero()
            || interval[1] > &Rational::one()
            || interval[0] >= interval[1]
        {
            return Err(unresolved("invalid offset Hermite chart interval"));
        }
        let bits = 128.min(dec.policy().max_precision_bits);
        if bits == 0 {
            return Err(unresolved("no offset chart proof precision available"));
        }
        let original = original(self.definition);
        let constant = |value: &Real| -> Result<Jet, Refusal> {
            Ok(Jet::constant(dec.enclose_at(value, bits)?))
        };
        let one = constant(&Real::one())?;
        let two = constant(&Real::from_i64(2))?;
        let radius = constant(&self.scale)?;
        let parameter = Jet::variable(b(interval[0], bits).hull(&b(interval[1], bits)));
        let t = radius.mul(&two.mul(&parameter).sub(&one));
        let denominator = one.add(&t.square());
        let c = one.sub(&t.square()).div(&denominator)?;
        let s = two.mul(&t).div(&denominator)?;
        let (c, s) = rotate(quarter, c, s, Jet::neg);
        let a = constant(&self.semiaxes[0])?;
        let b = constant(&self.semiaxes[1])?;
        let rho = constant(&original.minor_sq)?
            .mul(&c.square())
            .add(&constant(&original.major_sq)?.mul(&s.square()))
            .sqrt()?;
        let distance = constant(&self.distance)?;
        let x = a.mul(&c).add(&distance.mul(&b.mul(&c)).div(&rho)?);
        let y = b.mul(&s).add(&distance.mul(&a.mul(&s)).div(&rho)?);
        let axes = &self.axes;
        let mut coefficient = Rational::zero();
        for (u, v) in [
            (&axes[0].x, &axes[1].x),
            (&axes[0].y, &axes[1].y),
            (&axes[0].z, &axes[1].z),
        ] {
            let world = &dec.enclose_at(u, bits)? * &x.0[4] + &dec.enclose_at(v, bits)? * &y.0[4];
            coefficient = coefficient + world.lo().abs().max(world.hi().abs());
        }
        let width = interval[1] - interval[0];
        // Taylor coefficient four is f''''/24; Hermite's remainder constant is 1/384.
        let upper = coefficient * width.pow(4).unwrap() / Rational::from_i64(16);
        Ok(EvidenceBound {
            upper,
            unit: EvidenceUnit::Millimetres,
        })
    }
}
