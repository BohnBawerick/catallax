//! Exact clipped graph pieces of a retained regular section.

use super::*;
use crate::rational_patch::trace_series::series;
use catallax_numbers::Ball;

/// One ordered subarc in a certified section cell. Its free chart coordinate
/// varies affinely on [0, 1]; its other coordinate is the cell's unique root.
#[derive(Clone, Debug)]
pub struct SectionTracePiece {
    branch: PatchSectionBranch,
    cell: usize,
    ends: [[Real; 2]; 2],
    equation: RealPoly,
    transverse: RealPoly,
    slope: Sign,
}

impl PatchSectionBranch {
    /// Complete ordered cell pieces of a finite natural arc. Equal endpoints
    /// select a whole closed component. An open component needs both endpoints.
    pub fn trace_arc(
        &self,
        from: Option<&Point3>,
        to: Option<&Point3>,
        dec: &mut Decider,
    ) -> Result<Vec<SectionTracePiece>, Refusal> {
        let before = dec.stats().coincident_separation;
        let endpoint = |point: &Point3, dec: &mut Decider| -> Result<(Real, [Real; 2]), Refusal> {
            let (cell, fraction, at) = self.locate(point, dec)?.ok_or_else(|| {
                failure("trace endpoint is outside its retained section component")
            })?;
            let at = match at {
                Some(at) => at,
                None => self.cells[cell].ends
                    [usize::from(exact_sign(&fraction, dec)? != Sign::Zero)]
                .clone(),
            };
            Ok((Real::from_i64(cell as i64) + fraction, at))
        };
        let count = Real::from_i64(self.cells.len() as i64);
        let (start, mut end, start_at, end_at) = match (from, to) {
            (None, None) if self.closed => (
                Real::zero(),
                count.clone(),
                self.cells[0].ends[0].clone(),
                self.cells.last().unwrap().ends[1].clone(),
            ),
            (Some(from), Some(to)) => {
                let (start, start_at) = endpoint(from, dec)?;
                let (end, end_at) = endpoint(to, dec)?;
                (start, end, start_at, end_at)
            }
            _ => return Err(failure("an open section trace needs both finite endpoints")),
        };
        let order = exact_sign(&(&end - &start), dec)?;
        if order != Sign::Positive {
            if !self.closed {
                return Err(failure(
                    "open section trace endpoints oppose its natural direction",
                ));
            }
            end = &end + &count;
        }
        let mut pieces = Vec::new();
        for turn in 0..=usize::from(self.closed) {
            for (index, cell) in self.cells.iter().enumerate() {
                dec.check_deadline()?;
                let low = Real::from_i64((turn * self.cells.len() + index) as i64);
                let high = &low + &Real::one();
                let a = if exact_sign(&(&start - &low), dec)? == Sign::Positive {
                    start.clone()
                } else {
                    low.clone()
                };
                let b = if exact_sign(&(&end - &high), dec)? == Sign::Negative {
                    end.clone()
                } else {
                    high.clone()
                };
                if exact_sign(&(&b - &a), dec)? != Sign::Positive {
                    continue;
                }
                let first = if exact_sign(&(&a - &low), dec)? == Sign::Zero {
                    cell.ends[0].clone()
                } else {
                    start_at.clone()
                };
                let last = if exact_sign(&(&b - &high), dec)? == Sign::Zero {
                    cell.ends[1].clone()
                } else {
                    end_at.clone()
                };
                let equation = self.surface.equation(self.span(cell));
                let transverse = equation.partial(cell.solved);
                let middle = midpoint(&cell.domain).map(Real::from_rational);
                let slope = exact_sign(&transverse.eval_real(&middle), dec)?;
                if slope == Sign::Zero {
                    return Err(failure(
                        "a retained graph cell lost its strict transverse derivative",
                    ));
                }
                pieces.push(SectionTracePiece {
                    branch: self.clone(),
                    cell: index,
                    ends: [first, last],
                    equation,
                    transverse,
                    slope,
                });
            }
        }
        if pieces.is_empty() || dec.stats().coincident_separation != before {
            return Err(failure(
                "section trace clipping needs exact endpoint ownership",
            ));
        }
        Ok(pieces)
    }
}

impl SectionTracePiece {
    fn cell(&self) -> &ArcCell {
        &self.branch.cells[self.cell]
    }

    fn parameters(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 2], Refusal> {
        let bits = t.precision();
        if bits > dec.policy().max_precision_bits
            || t.lo() < Rational::zero()
            || t.hi() > Rational::one()
        {
            return Err(failure(
                "section trace enclosure is outside its finite parameter domain or precision policy",
            ));
        }
        let cell = self.cell();
        if t.lo() == t.hi() {
            for end in 0..2 {
                if t.lo() == Rational::from_i64(end as i64) {
                    return Ok([
                        dec.enclose_at(&self.ends[end][0], bits)?,
                        dec.enclose_at(&self.ends[end][1], bits)?,
                    ]);
                }
            }
        }
        let free = 1 - cell.solved;
        let start = dec.enclose_at(&self.ends[0][free], bits)?;
        let speed = dec.enclose_at(&(&self.ends[1][free] - &self.ends[0][free]), bits)?;
        let coordinate = &start + &(&speed * t);
        let mut parameters = box_balls(&cell.domain, bits);
        let free_low = coordinate.lo().max(cell.domain[free][0].clone());
        let free_high = coordinate.hi().min(cell.domain[free][1].clone());
        if free_low > free_high {
            return Err(failure(
                "section trace free-coordinate bounds lost their cell",
            ));
        }
        parameters[free] =
            Ball::from_rational(&free_low, bits).hull(&Ball::from_rational(&free_high, bits));
        let mut lo = cell.domain[cell.solved][0].clone();
        let mut hi = cell.domain[cell.solved][1].clone();
        let minimum = Rational::one().mul_pow2(-i64::from(bits));
        let target = (coordinate.hi() - coordinate.lo()).max(minimum);
        for _ in 0..bits {
            dec.check_deadline()?;
            if &hi - &lo <= target {
                break;
            }
            let old = [lo.clone(), hi.clone()];
            parameters[cell.solved] =
                Ball::from_rational(&lo, bits).hull(&Ball::from_rational(&hi, bits));
            let middle = (&lo + &hi).mul_pow2(-1);
            let derivative = self.transverse.eval_ball(&parameters, dec)?;
            let mut at_middle = parameters.clone();
            at_middle[cell.solved] = Ball::from_rational(&middle, bits);
            let residual = self.equation.eval_ball(&at_middle, dec)?;
            if let Some(delta) = residual.checked_div(&derivative) {
                let image = &Ball::from_rational(&middle, bits) - &delta;
                lo = lo.max(image.lo());
                hi = hi.min(image.hi());
                if lo > hi {
                    return Err(failure(
                        "section trace root contraction lost its certified graph",
                    ));
                }
            }
            if old == [lo.clone(), hi.clone()] {
                let width = &hi - &lo;
                for cut in [&lo + &width.mul_pow2(-2), &hi - &width.mul_pow2(-2)] {
                    let mut test = parameters.clone();
                    test[cell.solved] = Ball::from_rational(&cut, bits);
                    let value = self.equation.eval_ball(&test, dec)?;
                    if let Some(sign) = value.sign() {
                        let positive = sign == catallax_numbers::StrictSign::Positive;
                        if positive == (self.slope == Sign::Positive) {
                            hi = hi.min(cut);
                        } else {
                            lo = lo.max(cut);
                        }
                    }
                }
                if old == [lo.clone(), hi.clone()] {
                    break;
                }
            }
        }
        parameters[cell.solved] =
            Ball::from_rational(&lo, bits).hull(&Ball::from_rational(&hi, bits));
        Ok(parameters)
    }

    /// Exact point on the selected finite implicit graph, including clipped ends.
    pub fn point_at(&self, t: &Rational, dec: &mut Decider) -> Result<Point3, Refusal> {
        if t < &Rational::zero() || t > &Rational::one() {
            return Err(failure("section trace point is outside [0, 1]"));
        }
        let cell = self.cell();
        let span = self.branch.span(cell);
        for end in 0..2 {
            if t == &Rational::from_i64(end as i64) {
                return Ok(span.point(&self.ends[end]));
            }
        }
        let before_separation = dec.stats().coincident_separation;
        let free = 1 - cell.solved;
        let value = &self.ends[0][free]
            + &(&(&self.ends[1][free] - &self.ends[0][free]) * &Real::from(t.clone()));
        let restricted = self
            .equation
            .substitute(free, &RealPoly::constant(2, value.clone()));
        let polynomial = axis_polynomial(&restricted, cell.solved, dec)?;
        let mut found = None;
        for root in polynomial.real_roots(dec)? {
            if in_range(&root, &cell.domain[cell.solved], dec)? {
                if found.is_some() {
                    return Err(failure(
                        "section trace selected more than one transverse root",
                    ));
                }
                let mut at = [value.clone(), value.clone()];
                at[cell.solved] = root;
                found = Some(span.point(&at));
            }
        }
        if dec.stats().coincident_separation != before_separation {
            return Err(failure(
                "section trace root selection needs exact cell membership",
            ));
        }
        found.ok_or_else(|| failure("section trace selected no transverse root"))
    }

    /// Taylor coefficient enclosures valid at every parameter in `t`. The
    /// implicit-function recurrence divides only by a strict transverse bound.
    pub fn taylor(&self, t: &Ball, order: usize, dec: &Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        let at = self.parameters(t, dec)?;
        let bits = t.precision();
        let zero = || Ball::from_rational(&Rational::zero(), bits);
        let cell = self.cell();
        let free = 1 - cell.solved;
        let mut parameters: [Vec<Ball>; 2] = std::array::from_fn(|axis| {
            let mut v = vec![zero(); order + 1];
            v[0] = at[axis].clone();
            v
        });
        if order > 0 {
            parameters[free][1] =
                dec.enclose_at(&(&self.ends[1][free] - &self.ends[0][free]), bits)?;
            let derivative = self.transverse.eval_ball(&at, dec)?;
            for n in 1..=order {
                let residual = series(&self.equation, &parameters, n, bits, dec)?
                    .pop()
                    .unwrap();
                parameters[cell.solved][n] =
                    (-&residual).checked_div(&derivative).ok_or_else(|| {
                        failure("section trace jets need a strict transverse derivative enclosure")
                    })?;
            }
        }
        let span = self.branch.span(cell);
        let mut homogeneous = Vec::with_capacity(4);
        for coordinate in &span.coordinates {
            homogeneous.push(series(coordinate, &parameters, order, bits, dec)?);
        }
        let mut result: Vec<[Ball; 3]> = Vec::with_capacity(order + 1);
        for n in 0..=order {
            let mut xyz: [Ball; 3] = std::array::from_fn(|_| zero());
            for axis in 0..3 {
                let mut numerator = homogeneous[axis][n].clone();
                for k in 1..=n {
                    numerator = &numerator - &(&homogeneous[3][k] * &result[n - k][axis]);
                }
                xyz[axis] = numerator.checked_div(&homogeneous[3][0]).ok_or_else(|| {
                    failure("section trace homogeneous denominator needs a positive enclosure")
                })?;
            }
            result.push(xyz);
        }
        Ok(result)
    }
}
