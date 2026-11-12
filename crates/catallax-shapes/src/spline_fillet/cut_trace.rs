//! Enclosures of a proved cap cut, without replacing it by its STEP fit.

use super::*;

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
