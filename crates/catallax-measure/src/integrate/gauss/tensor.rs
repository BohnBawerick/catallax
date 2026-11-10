//! Product Gauss integration with complete pure-partial remainder bounds.

use super::*;

impl Program {
    pub(crate) fn parameter_tensor_rectangle(
        &self,
        u: [&Rational; 2],
        v: [&Rational; 2],
        work: &mut Work,
    ) -> Checked<Vec<Ball>> {
        let ub = work.interval(u[0], u[1])?;
        let vb = work.interval(v[0], v[1])?;
        let area = work.rational(&((u[1] - u[0]) * (v[1] - v[0])))?;
        self.evaluate(&[[ub], [vb]], Coefficients::All, work)?
            .into_iter()
            .map(|[value]| work.mul(&value, &area))
            .collect()
    }

    /// Iu Iv - Qu Qv = Iv (Iu - Qu) + Qu (Iv - Qv). Positive Gauss
    /// weights sum to each interval's width, so complete-box pure-partial
    /// remainders bound both terms. Samples supply only Qu Qv, never the proof.
    pub(crate) fn parameter_tensor_panel(
        &self,
        u: [&Rational; 2],
        v: [&Rational; 2],
        rule: &Rule,
        work: &mut Work,
    ) -> Checked<Vec<Ball>> {
        if u[0] >= u[1] || v[0] >= v[1] {
            return Err(certificate());
        }
        let widths = [u[1] - u[0], v[1] - v[0]];
        let centers = [
            work.rational(&(u[0] + u[1]).mul_pow2(-1))?,
            work.rational(&(v[0] + v[1]).mul_pow2(-1))?,
        ];
        let halves = [
            work.rational(&widths[0].mul_pow2(-1))?,
            work.rational(&widths[1].mul_pow2(-1))?,
        ];
        let ranges = [work.interval(u[0], u[1])?, work.interval(v[0], v[1])?];
        let mut nodes = [Vec::with_capacity(8), Vec::with_capacity(8)];
        for axis in 0..2 {
            for (point, weight) in &rule.points {
                let offset = work.mul(&halves[axis], point)?;
                nodes[axis].push((work.add(&centers[axis], &offset)?, weight));
            }
        }
        let zero = work.rational(&Rational::zero())?;
        let mut sums = vec![zero; self.outputs.len()];
        for (u, wu) in &nodes[0] {
            for (v, wv) in &nodes[1] {
                let weight = work.mul(wu, wv)?;
                let values = self.evaluate(&[[u.clone()], [v.clone()]], Coefficients::All, work)?;
                for (sum, [value]) in sums.iter_mut().zip(values) {
                    let value = work.mul(&value, &weight)?;
                    *sum = work.add(sum, &value)?;
                }
            }
        }
        let scale = work.mul(&halves[0], &halves[1])?;
        for sum in &mut sums {
            *sum = work.mul(sum, &scale)?;
        }
        let one = work.rational(&Rational::one())?;
        for axis in 0..2 {
            let mut input: [[Ball; SIZE]; 2] = std::array::from_fn(|_| work.zeros());
            input[0][0] = ranges[0].clone();
            input[1][0] = ranges[1].clone();
            input[axis][1] = one.clone();
            let remainder = self.evaluate(&input, Coefficients::GaussRemainder, work)?;
            let factor = widths[axis].pow(17).expect("positive power") * widths[1 - axis].clone()
                / Rational::from_i64(17 * 12870 * 12870);
            let factor = work.rational(&factor)?;
            for (sum, remainder) in sums.iter_mut().zip(remainder) {
                let error = work.mul(&remainder[ORDER], &factor)?;
                *sum = work.add(sum, &error)?;
            }
        }
        Ok(sums)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrate::Builder;
    use catallax_numbers::{Decider, with_operation_budget};

    #[test]
    fn product_rule_bounds_a_nonzero_remainder_in_both_directions() {
        let ((), work_used) = with_operation_budget(10_000_000, || {
            let dec = Decider::default();
            let policy = *dec.policy();
            let mut work = Work::new(&policy, 128, policy.max_work, &dec).unwrap();
            let rule = work.gauss_rule().unwrap();
            let k = Builder::default();
            fn power<'a>(mut value: crate::integrate::Expr<'a>) -> crate::integrate::Expr<'a> {
                for _ in 0..4 {
                    value = value.square();
                }
                value
            }
            let program = k.finish(&[power(k.input(0)) * power(k.input(1))]);
            let zero = Rational::zero();
            let one = Rational::one();
            let value = program
                .parameter_tensor_panel([&zero, &one], [&zero, &one], &rule, &mut work)
                .unwrap();
            let exact = Rational::from_ratio(1, 17 * 17).unwrap();
            assert!(value[0].lo() <= exact && value[0].hi() >= exact);
            assert!(value[0].hi() - value[0].lo() < Rational::parse("1e-6").unwrap());
        });
        assert!(!work_used.exhausted);
    }
}
