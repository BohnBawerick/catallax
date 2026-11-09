use super::*;
use catallax_numbers::Ball;

fn interval(lo: &Rational, hi: &Rational, bits: u32) -> Ball {
    Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits))
}

impl BezierExtrusion {
    /// Encloses every preimage of a spatial box on the finite patch.
    /// The cover proves no incidence; monotone projection selects one interval.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        dec.check_deadline()?;
        let precision = point
            .iter()
            .map(Ball::precision)
            .max()
            .expect("three coordinates");
        if precision == 0 || precision > dec.policy().max_precision_bits {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "Bezier inverse precision exceeds the query policy",
            ));
        }
        let bits = precision
            .saturating_add(32)
            .min(dec.policy().max_precision_bits);
        let origin = math::point(&self.controls()[0], bits, dec)?;
        let offset: [Ball; 3] = std::array::from_fn(|i| &point[i] - &origin[i]);
        let project = |axis: &Vector3, norm: &Real| -> Result<Ball, Refusal> {
            let axis = [
                dec.enclose_at(&axis.x, bits)?,
                dec.enclose_at(&axis.y, bits)?,
                dec.enclose_at(&axis.z, bits)?,
            ];
            let dot = (0..3).fold(Ball::from_rational(&Rational::zero(), bits), |sum, i| {
                &sum + &(&offset[i] * &axis[i])
            });
            dot.checked_div(&dec.enclose_at(norm, bits)?)
                .ok_or_else(|| {
                    query_domain("the projection denominator is not separated from zero")
                })
        };
        let along = project(&self.data.chord, &self.data.chord_sq)?;
        let height = project(self.direction(), &self.data.axis_sq)?;
        let lo = along.lo().max(Rational::zero());
        let hi = along.hi().min(Rational::one());
        let v_lo = height.lo().max(Rational::zero());
        let v_hi = height.hi().min(Rational::one());
        if lo > hi || v_lo > v_hi {
            return Ok(None);
        }
        let u_lo = self.projection_inverse(&lo, bits, dec)?.0;
        let u_hi = self.projection_inverse(&hi, bits, dec)?.1;
        Ok(Some([
            interval(&u_lo, &u_hi, precision),
            interval(&v_lo, &v_hi, precision),
        ]))
    }

    fn projection_inverse(
        &self,
        target: &Rational,
        bits: u32,
        dec: &Decider,
    ) -> Result<(Rational, Rational), Refusal> {
        if target.is_zero() || *target == Rational::one() {
            return Ok((target.clone(), target.clone()));
        }
        let mut lo = Rational::zero();
        let mut hi = Rational::one();
        // An unfinished bounded search still returns a complete bracket.
        for _ in 0..bits.min(256) {
            dec.check_deadline()?;
            let mid = (&lo + &hi) / Rational::from_i64(2);
            let complement = &Rational::one() - &mid;
            let mut row = self.data.projection.clone();
            for count in (1..row.len()).rev() {
                for i in 0..count {
                    row[i] = &(&row[i] * &complement) + &(&row[i + 1] * &mid);
                }
            }
            match row[0].cmp(target) {
                std::cmp::Ordering::Less => lo = mid,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => return Ok((mid.clone(), mid)),
            }
        }
        Ok((lo, hi))
    }
}
