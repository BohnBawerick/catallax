use super::*;
use std::{borrow::Cow, cell::RefCell};

mod cache;
mod unit;

#[cfg(test)]
mod unit_tests;

/// Derivative numerators in the original span's unit parameters. The atlas
/// cover uses these only for exact Bezier restrictions of that same span.
pub(crate) struct PolynomialBounds {
    domain: Domain,
    widths: [Rational; 2],
    weight: Grid,
    tangents: [[Grid; 3]; 2],
    cache: RefCell<cache::RestrictionCache>,
}

pub(super) fn quotient(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a / b)
}

pub(super) fn bernstein(poly: &Dense, dec: &mut Decider) -> Result<Grid, Refusal> {
    let unsupported = || {
        Refusal::new(
            "shapes.rational-atlas-polynomial-degree",
            "the exact Bernstein binomial recurrence exceeds its integer representation",
        )
    };
    for degree in poly.deg {
        let n = i64::try_from(degree).map_err(|_| unsupported())?;
        let mut coefficient = 1i64;
        for i in 0..n {
            coefficient = coefficient
                .checked_mul(n - i)
                .map(|product| product / (i + 1))
                .ok_or_else(unsupported)?;
        }
    }
    poly.bernstein(dec)
}

pub(super) fn restriction(
    domain: &Domain,
    widths: &[Rational; 2],
    cell: &Domain,
    dec: &mut Decider,
) -> Result<[[Rational; 2]; 2], Refusal> {
    let mut cuts = std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
    for (axis, bounds) in cell.iter().enumerate() {
        if less(&bounds[0], &domain[axis][0], dec)?
            || less(&domain[axis][1], &bounds[1], dec)?
            || !less(&bounds[0], &bounds[1], dec)?
        {
            return Err(Refusal::new(
                "shapes.rational-atlas-polynomial-domain",
                "a derivative cell must be a nonempty restriction of its source span",
            ));
        }
        let low = sub(&bounds[0], &domain[axis][0], dec)?;
        let high = sub(&bounds[1], &domain[axis][0], dec)?;
        cuts[axis] = [
            quotient(&low, &high, dec)?,
            quotient(&high, &widths[axis], dec)?,
        ];
    }
    Ok(cuts)
}

pub(super) fn restricted_net<'a>(
    grid: &'a Grid,
    cuts: &[[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<Cow<'a, Grid>, Refusal> {
    let mut grid = Cow::Borrowed(grid);
    for (axis, [low, high]) in cuts.iter().enumerate() {
        if *high != Rational::one() {
            let [next, _] = grid.split(axis, high, dec)?;
            grid = Cow::Owned(next);
        }
        if !low.is_zero() {
            let [_, next] = grid.split(axis, low, dec)?;
            grid = Cow::Owned(next);
        }
    }
    Ok(grid)
}

pub(super) fn restricted_range(
    grid: &Grid,
    cuts: &[[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    restricted_net(grid, cuts, dec)?.range(dec)
}

impl PolynomialBounds {
    /// Cancels derivative numerators from the defining span's exact unit net.
    pub(crate) fn from_span(
        source: &crate::BSplineSurface,
        span: [usize; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let local = unit::polynomials(source, span, dec)?;
        let knots = source.knots();
        let domain = std::array::from_fn(|axis| {
            [
                knots[axis][span[axis]].clone(),
                knots[axis][span[axis] + 1].clone(),
            ]
        });
        Self::from_local(local, domain, dec)
    }

    #[cfg(test)]
    pub(crate) fn new(
        homogeneous: &[MPoly; 4],
        domain: Domain,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let local = [
            Dense::local(&homogeneous[0], &domain, dec)?,
            Dense::local(&homogeneous[1], &domain, dec)?,
            Dense::local(&homogeneous[2], &domain, dec)?,
            Dense::local(&homogeneous[3], &domain, dec)?,
        ];
        Self::from_local(local, domain, dec)
    }

    fn from_local(local: [Dense; 4], domain: Domain, dec: &mut Decider) -> Result<Self, Refusal> {
        let weight = bernstein(&local[3], dec)?;
        if !weight.all_positive() {
            return Err(Refusal::new(
                "shapes.rational-atlas-polynomial-weights",
                "the exact span weight has no strictly positive Bernstein hull",
            ));
        }
        let mut tangent = |axis| -> Result<[Grid; 3], Refusal> {
            let dw = local[3].derivative(axis, dec)?;
            let mut component = |coordinate: &Dense| -> Result<Grid, Refusal> {
                let dn = coordinate.derivative(axis, dec)?;
                let left = local[3].mul(&dn, dec)?;
                let right = coordinate.mul(&dw, dec)?;
                bernstein(&left.sub(&right, dec)?, dec)
            };
            Ok([
                component(&local[0])?,
                component(&local[1])?,
                component(&local[2])?,
            ])
        };
        let tangents = [tangent(0)?, tangent(1)?];
        let widths = [
            sub(&domain[0][1], &domain[0][0], dec)?,
            sub(&domain[1][1], &domain[1][0], dec)?,
        ];
        Ok(Self {
            domain,
            widths,
            weight,
            tangents,
            cache: RefCell::default(),
        })
    }

    pub(crate) fn derivatives(
        &self,
        cell: &Domain,
        bits: u32,
        dec: &mut Decider,
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        let cuts = restriction(&self.domain, &self.widths, cell, dec)?;
        let mut cache = self.cache.borrow_mut();
        let prepared = cache.prepare(&self.domain, cell, cuts, dec)?;
        let mut range =
            |grid: &Grid, component: usize, dec: &mut Decider| -> Result<Ball, Refusal> {
                let [lo, hi] = cache.range(grid, component, &prepared, dec)?;
                Ok(Ball::from_rational(&lo, bits).hull(&Ball::from_rational(&hi, bits)))
            };
        let weight = range(&self.weight, 0, dec)?;
        if weight.lo() <= Rational::zero() {
            return Err(Refusal::new(
                "shapes.rational-atlas-polynomial-weights",
                "the arithmetic precision does not separate the positive span weight from zero",
            ));
        }
        let squared = weight.square();
        let mut derivative = |axis: usize| -> Result<[Ball; 3], Refusal> {
            let denominator = &squared * &Ball::from_rational(&self.widths[axis], bits);
            let mut component = |coordinate: usize| {
                range(
                    &self.tangents[axis][coordinate],
                    1 + axis * 3 + coordinate,
                    dec,
                )?
                .checked_div(&denominator)
                .ok_or_else(|| {
                    Refusal::new(
                        "shapes.rational-atlas-polynomial-weights",
                        "the derivative denominator has no strictly positive enclosure",
                    )
                })
            };
            Ok([component(0)?, component(1)?, component(2)?])
        };
        Ok([derivative(0)?, derivative(1)?])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BSplineSurface, Point3, RationalPatch};
    use catallax_numbers::with_operation_budget;

    fn q(n: i64, d: i64) -> Rational {
        Rational::from_ratio(n, d).unwrap()
    }

    #[test]
    fn derivative_hulls_cover_translated_nonunit_rectangles_after_cancellation() {
        let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
            let mut dec = Decider::default();
            let u = MPoly::var(2, 0);
            let v = MPoly::var(2, 1);
            let w = &MPoly::constant(2, Rational::one()) + &u;
            let homogeneous = [u, &w * &v, MPoly::zero(2), w];
            let bounds = PolynomialBounds::new(
                &homogeneous,
                [[q(2, 1), q(4, 1)], [q(-3, 1), q(1, 1)]],
                &mut dec,
            )?;
            for (cell, lo, hi) in [
                ([[q(2, 1), q(4, 1)], [q(-3, 1), q(1, 1)]], q(1, 25), q(1, 9)),
                (
                    [[q(5, 2), q(3, 1)], [q(-2, 1), q(0, 1)]],
                    q(1, 16),
                    q(4, 49),
                ),
            ] {
                let [du, dv] = bounds.derivatives(&cell, 128, &mut dec)?;
                // X=u/(1+u), Y=v, Z=0 throughout each rectangle.
                assert!(du[0].contains(&lo) && du[0].contains(&hi));
                assert!(du[0].lo() > Rational::zero() && du[0].hi() < q(1, 8));
                assert!(du[1].is_exact_zero() && du[2].is_exact_zero());
                assert!(dv[0].is_exact_zero() && dv[2].is_exact_zero());
                assert!(dv[1].contains(&Rational::one()));
            }
            Ok(())
        });
        assert!(!work.exhausted);
        result.unwrap();
    }

    #[test]
    fn polynomial_bounds_reject_an_interior_fold_away_from_the_frame_point() {
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let carrier = BSplineSurface::new(
                [2, 1],
                [
                    [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
                    [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
                ],
                [3, 2],
                // X=(u-1/4)^2, Y=v: regular at the centre, singular at u=1/4.
                [1, -3, 9]
                    .into_iter()
                    .flat_map(|x| {
                        [0, 1]
                            .map(move |y| Point3::new(Real::from(q(x, 16)), y.into(), Real::zero()))
                    })
                    .collect(),
                vec![Real::one(); 6],
                &mut dec,
            )?;
            let homogeneous = carrier
                .span_polynomials([2, 1], &dec)?
                .each_ref()
                .map(|p| p.as_rational().unwrap());
            let bounds = PolynomialBounds::new(
                &homogeneous,
                std::array::from_fn(|_| [Rational::zero(), Rational::one()]),
                &mut dec,
            )?;
            RationalPatch::new_cell_with_polynomial_bounds(carrier, &bounds, &mut dec)
        });
        assert!(!work.exhausted);
        assert_eq!(result.unwrap_err().code, "shapes.rational-patch-unproven");
    }

    #[test]
    fn large_binomial_recurrences_refuse_before_integer_overflow() {
        let mut polynomial = Dense::zero([67, 0]);
        *polynomial.at_mut(67, 0) = Rational::one();
        let error = bernstein(&polynomial, &mut Decider::default()).unwrap_err();
        assert_eq!(error.code, "shapes.rational-atlas-polynomial-degree");
    }
}
