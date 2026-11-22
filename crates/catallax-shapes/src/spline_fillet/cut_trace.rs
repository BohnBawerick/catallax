//! Enclosures of a proved cap cut, without replacing it by its STEP fit.

use super::*;
use catallax_numbers::Interval;

/// A proved plane cut, parametrised by the carrier's second parameter.
/// Its private coefficients retain the monotonicity proof of `plane_cut`.
#[derive(Clone, Debug)]
pub struct FilletCutTrace {
    fillet: SplineFillet,
    cut: FilletCut,
    coefficients: Vec<[Vec<Rational>; 3]>,
}

fn evaluate(mut row: Vec<Ball>, t: &Ball) -> Ball {
    for level in 1..row.len() {
        for i in 0..row.len() - level {
            row[i] = &row[i] + &(t * &(&row[i + 1] - &row[i]));
        }
    }
    row[0].clone()
}

fn difference(row: &[Interval]) -> Vec<Interval> {
    let degree = Interval::point(row.len().saturating_sub(1) as f64).unwrap();
    row.windows(2)
        .map(|pair| degree * (pair[1] - pair[0]))
        .collect()
}

fn hull(values: impl IntoIterator<Item = Interval>) -> Interval {
    values
        .into_iter()
        .reduce(|a, b| a.hull(&b))
        .unwrap_or_else(|| Interval::point(0.0).unwrap())
}

impl SplineFillet {
    /// Proves a plane cut and prepares its exact equation for point queries.
    pub fn cut_trace(
        &self,
        normal: [Rational; 3],
        offset: Rational,
        dec: &mut Decider,
    ) -> Result<Option<FilletCutTrace>, Refusal> {
        let Some(cut) = self.plane_cut(normal, offset, dec)? else {
            return Ok(None);
        };
        let (normal, offset) = surface::local_plane(&self.data.frame, &cut.normal, &cut.offset)
            .expect("plane_cut checked the rational frame");
        let coefficients = (cut.spans[0]..=cut.spans[1])
            .map(|span| cut::coefficients(&self.data.fit.rows[span], &normal, &offset, dec))
            .collect::<Result<_, _>>()?;
        Ok(Some(FilletCutTrace {
            fillet: self.clone(),
            cut,
            coefficients,
        }))
    }
}

impl FilletCutTrace {
    /// The checked cut, including its exact rail endpoints.
    pub fn cut(&self) -> &FilletCut {
        &self.cut
    }

    /// The exact endpoint on rail `side`.
    pub fn endpoint(&self, side: usize) -> Point3 {
        self.fillet
            .rail_point(side, &self.cut.ends[side])
            .expect("a cut ends at rail nodes")
    }

    /// Upper bound on the spatial acceleration with respect to theta.
    /// The display sampler uses this to plan chords of the proved plane cut.
    pub fn display_acceleration(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        let carrier = self.fillet.carrier();
        let [degree, across] = carrier.degrees();
        let fit = &self.fillet.data.fit;
        let bits = dec.policy().max_precision_bits.min(128);
        let zero = Interval::point(0.0).unwrap();
        let twice = Interval::point(2.0).unwrap();
        let mut maximum = 0.0_f64;
        let mut workspace = [Vec::new(), Vec::new()];
        for (span, rows) in (self.cut.spans[0]..).zip(&self.coefficients) {
            let g = rows
                .each_ref()
                .map(|row| row.iter().map(Interval::from_rational).collect::<Vec<_>>());
            let gs = g.each_ref().map(|row| difference(row));
            let gt = std::array::from_fn::<_, 2, _>(|i| {
                g[i + 1]
                    .iter()
                    .zip(&g[i])
                    .map(|(&b, &a)| twice * (b - a))
                    .collect::<Vec<_>>()
            });
            dec.charge_work(
                (3 * degree + 2 * (degree + 1) + 3 * (degree - 1) + 2 * degree) as u64,
            )?;
            let gss = hull(gs.iter().flat_map(|row| difference(row)));
            let gst = hull(gt.iter().flat_map(|row| difference(row)));
            let gtt = hull(gt[1].iter().zip(&gt[0]).map(|(&b, &a)| b - a));
            let gs = hull(gs.into_iter().flatten());
            if gs.lo() <= 0.0 && gs.hi() >= 0.0 {
                return Err(Refusal::new(
                    SPLINE_FILLET_CUT,
                    "binary64 could not separate the plane cut's monotone derivative",
                ));
            }
            // Implicit differentiation of g(s(theta), theta) = 0.
            // The carrier composes Hermite spans with shared first jets, so
            // crossing a knot does not introduce a jump in the curve tangent.
            let rate = -hull(gt.into_iter().flatten()) / gs;
            let rate_sq = rate.square();
            let accel = -(gtt + twice * gst * rate + gss * rate_sq) / gs;
            let domain = [
                [fit.nodes[span].tau.clone(), fit.nodes[span + 1].tau.clone()],
                [Rational::zero(), Rational::one()],
            ];
            let width = Interval::from_rational(&(&domain[0][1] - &domain[0][0]));
            let jet = carrier.jet_on_span(
                [(span + 1) * degree, across],
                &domain,
                bits,
                dec,
                &mut workspace,
            )?;
            let tau_rate = rate * width;
            let tau_accel = accel * width;
            let tau_rate_sq = tau_rate.square();
            dec.charge_work(2 + 4 + 1 + 3 * 4)?;
            let bound = (0..3).fold(zero, |sum, i| {
                let acceleration = jet.second[2][i].to_interval()
                    + twice * jet.second[1][i].to_interval() * tau_rate
                    + jet.second[0][i].to_interval() * tau_rate_sq
                    + jet.first[0][i].to_interval() * tau_accel;
                sum + Interval::point(acceleration.lo().abs().max(acceleration.hi().abs())).unwrap()
            });
            maximum = maximum.max(bound.hi());
        }
        Rational::from_f64_exact(maximum).ok_or_else(|| {
            Refusal::new(SPLINE_FILLET_CUT, "plane cut acceleration exceeds binary64")
        })
    }

    /// Encloses every cut point for `theta` in `[0,1]`. Bisection only
    /// tightens a box; an undecided sign keeps it, and never discards a root.
    pub fn enclose(&self, theta: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        let bits = theta.precision().min(dec.policy().max_precision_bits);
        let zero = Rational::zero();
        let one = Rational::one();
        if bits == 0 || theta.lo() < zero || theta.hi() > one {
            return Err(Refusal::new(
                SPLINE_FILLET_CUT,
                "cut parameter is outside [0,1] or has no precision",
            ));
        }
        if theta.lo() == theta.hi() && (theta.lo() == zero || theta.lo() == one) {
            let endpoint = self.endpoint(usize::from(theta.lo() == one));
            let [x, y, z] =
                [&endpoint.x, &endpoint.y, &endpoint.z].map(|v| dec.enclose_at(v, bits));
            return Ok([x?, y?, z?]);
        }
        let fit = &self.fillet.data.fit;
        let [first, last] = self.cut.spans;
        let ball = |q: &Rational| Ball::from_rational(q, bits);
        let theta = theta.with_precision(bits);
        let side = |tau: &Rational| -> Result<Option<Sign>, Refusal> {
            dec.check_deadline()?;
            let span = (first..=last)
                .find(|&i| tau <= &fit.nodes[i + 1].tau)
                .unwrap_or(last);
            let start = &fit.nodes[span].tau;
            let width = &fit.nodes[span + 1].tau - start;
            let s = ball(&((tau - start) / width));
            let rows = self.coefficients[span - first]
                .each_ref()
                .map(|row| evaluate(row.iter().map(&ball).collect(), &s));
            let value = evaluate(rows.to_vec(), &theta);
            Ok(if value.lo() > zero {
                Some(Sign::Positive)
            } else if value.hi() < zero {
                Some(Sign::Negative)
            } else {
                None
            })
        };
        let range = [
            fit.nodes[first].tau.clone(),
            fit.nodes[last + 1].tau.clone(),
        ];
        // The lower bracket stays before every root; the upper stays after
        // every root, even when theta is an interval with many root values.
        let mut lower = range.clone();
        let mut upper = range;
        for _ in 0..24 {
            let middle = (&lower[0] + &lower[1]).mul_pow2(-1);
            if side(&middle)? == Some(self.cut.before) {
                lower[0] = middle;
            } else {
                lower[1] = middle;
            }
            let middle = (&upper[0] + &upper[1]).mul_pow2(-1);
            if side(&middle)?.is_some_and(|sign| sign != self.cut.before) {
                upper[1] = middle;
            } else {
                upper[0] = middle;
            }
        }
        let tau = ball(&lower[0]).hull(&ball(&upper[1]));
        self.fillet.carrier().evaluate(&tau, &theta, dec)
    }
}
