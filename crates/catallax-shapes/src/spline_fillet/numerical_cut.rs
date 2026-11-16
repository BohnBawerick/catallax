//! Sampled end-plane cuts of numerical fillet carriers.

use super::*;
use super::super::carrier::cubic;

/// Numerical world and native-UV curves for one end-plane cut.
#[derive(Clone, Debug)]
pub struct NumericalFilletCut {
    curve: BSplineCurve,
    parameters: BSplineCurve,
    tau_ends: [Rational; 2],
    sampled_error_mm: Rational,
}
impl NumericalFilletCut {
    /// Three-dimensional curve, parameterized from the first rail to the second.
    pub fn curve(&self) -> &BSplineCurve { &self.curve }
    /// Native fillet parameters at the same curve parameter.
    pub fn parameters(&self) -> &BSplineCurve { &self.parameters }
    /// First-support spine parameters where the cut meets both rails.
    pub fn tau_ends(&self) -> &[Rational; 2] { &self.tau_ends }
    /// Largest checked plane or world/UV residual, in millimetres.
    pub fn sampled_error_mm(&self) -> &Rational { &self.sampled_error_mm }
}
fn plane_value(point: &[Ball; 3], normal: &[Rational; 3], offset: &Rational) -> Ball {
    (0..3).fold(-&Ball::from_rational(offset, 128), |sum, i| {
        &sum + &(&point[i] * &Ball::from_rational(&normal[i], 128))
    })
}
fn magnitude(value: &Ball) -> Rational { value.lo().abs().max(value.hi().abs()) }

impl NumericalSplineFillet {
    /// Intersects an end plane numerically as one graph `tau(s)` across the
    /// carrier. Bisection and finite residual checks have no exact authority.
    /// A non-bracketed or unresolved cut returns its cause rather than a curve.
    pub fn plane_cut(
        &self,
        normal: [Rational; 3],
        offset: Rational,
        dec: &mut Decider,
    ) -> Result<NumericalFilletCut, Refusal> {
        let norm_square = normal.iter().fold(Rational::zero(), |sum, n| sum + n * n);
        if norm_square.sign() != Sign::Positive { return Err(bad("fillet end plane has zero normal")); }
        let length = Ball::from_rational(&norm_square, 128).sqrt()
            .ok_or_else(|| bad("fillet end plane normal length is unresolved"))?;
        let lower = length.lo();
        if lower.sign() != Sign::Positive { return Err(bad("fillet end plane normal is not separated from zero")); }
        let tolerance = self.spine().tolerance_mm().mul_pow2(-2);
        let root_tolerance = &tolerance * &lower;
        let source = self.surface().source();
        let domain = source.domain()[0];
        let root = |s: &Rational, dec: &mut Decider| -> Result<(Rational, [Rational; 3]), Refusal> {
            let mut lo = domain[0].clone();
            let mut hi = domain[1].clone();
            let at = |tau: &Rational, dec: &Decider| source.evaluate_box([tau, tau], [s, s], 128, dec);
            let a = plane_value(&at(&lo, dec)?, &normal, &offset).mid();
            let b = plane_value(&at(&hi, dec)?, &normal, &offset).mid();
            if a.sign() == b.sign() && a.sign() != Sign::Zero {
                return Err(bad("fillet end plane does not bracket a transverse cut across the carrier"));
            }
            for _ in 0..80 {
                dec.check_deadline()?;
                let tau = (&lo + &hi).mul_pow2(-1);
                let point = at(&tau, dec)?;
                let value = plane_value(&point, &normal, &offset);
                dec.charge_work(32)?;
                if magnitude(&value) <= root_tolerance.mul_pow2(-4) {
                    return Ok((tau, point.each_ref().map(Ball::mid)));
                }
                if value.mid().sign() == a.sign() { lo = tau; } else { hi = tau; }
            }
            Err(bad("fillet end-plane bisection cannot meet the stated tolerance"))
        };
        for density in [8usize, 16, 32, 64, 128, 256] {
            dec.charge_work(density as u64 * 128)?;
            let mut samples = Vec::with_capacity(density + 1);
            for i in 0..=density {
                let s = Rational::from_ratio(i as i64, density as i64).unwrap();
                let (tau, point) = root(&s, dec)?;
                samples.push([point[0].clone(), point[1].clone(), point[2].clone(), tau, s]);
            }
            let mut slopes: Vec<[Rational; 5]> = Vec::with_capacity(samples.len());
            for i in 0..=density {
                let a = i.saturating_sub(1);
                let b = (i + 1).min(density);
                let width = &samples[b][4] - &samples[a][4];
                slopes.push(std::array::from_fn(|j| (&samples[b][j] - &samples[a][j]) / width.clone()));
            }
            let mut controls = vec![samples[0].clone()];
            let mut knots = vec![Rational::zero(); 4];
            for i in 0..density {
                let span = [&samples[i][4], &samples[i + 1][4]];
                let rows: [Vec<Rational>; 5] = std::array::from_fn(|j| cubic(span, [&samples[i][j], &samples[i + 1][j]], [&slopes[i][j], &slopes[i + 1][j]]));
                for k in 1..4 { controls.push(std::array::from_fn(|j| rows[j][k].clone())); }
                if i + 1 < density { knots.extend(std::iter::repeat_n(samples[i + 1][4].clone(), 3)); }
            }
            knots.extend(std::iter::repeat_n(Rational::one(), 4));
            let world = controls.iter().map(|p| Point3::new(p[0].clone().into(), p[1].clone().into(), p[2].clone().into())).collect();
            let uv = controls.iter().map(|p| Point3::new(p[3].clone().into(), p[4].clone().into(), Real::zero())).collect();
            let curve = BSplineCurve::new(3, knots.clone(), world, vec![Real::one(); controls.len()], dec)?;
            let parameters = BSplineCurve::new(3, knots, uv, vec![Real::one(); controls.len()], dec)?;
            let mut upper = Rational::zero();
            for i in 0..density {
                for j in 0..=8 {
                    let s = Rational::from_ratio((8 * i + j) as i64, (8 * density) as i64).unwrap();
                    let point = curve.evaluate_box([&s, &s], 128, dec)?;
                    let uv = parameters.evaluate_box([&s, &s], 128, dec)?;
                    let tau = uv[0].mid();
                    let across = uv[1].mid();
                    let composed = source.evaluate_box([&tau, &tau], [&across, &across], 128, dec)?;
                    let gap = (0..3).fold(Rational::zero(), |sum, k| {
                        sum + (&point[k].lo() - &composed[k].hi()).abs().max((&point[k].hi() - &composed[k].lo()).abs())
                    });
                    let plane_error = magnitude(&plane_value(&point, &normal, &offset)) / lower.clone();
                    upper = upper.max(gap).max(plane_error);
                }
            }
            if upper <= tolerance {
                return Ok(NumericalFilletCut {
                    curve, parameters,
                    tau_ends: [samples[0][3].clone(), samples[density][3].clone()],
                    sampled_error_mm: upper,
                });
            }
        }
        Err(bad("fillet end-plane curves cannot meet the stated sampled tolerance"))
    }
}
