//! Implicit graph pieces of the existing complete patch-pair cover.

use super::*;
use crate::rational_patch::trace_series::series;
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
struct Definition {
    first: RetainedSource,
    second: RetainedSource,
    branch: PatchPairBranch,
}

#[derive(Debug)]
struct Equations {
    polynomials: [RealPoly; 3],
    jacobian: [[RealPoly; 4]; 3],
}

/// One producer-owned regular graph in four patch parameters. Its fifth
/// polynomial variable is the unit trace parameter, not another unknown.
/// The selected branch and its complete cell-pair certificate remain retained.
#[derive(Clone, Debug)]
pub struct PairTracePiece {
    definition: Arc<Definition>,
    arc: usize,
    system: Arc<Equations>,
    parameter: RealPoly,
    parameter_jacobian: [RealPoly; 4],
}

fn failure(detail: &str) -> Refusal {
    Refusal::new("intersect.patch-pair-trace", detail)
}

fn lift(poly: &RealPoly, targets: &[usize], variables: usize) -> RealPoly {
    RealPoly::from_terms(
        variables,
        poly.terms()
            .into_iter()
            .map(|(monomial, coefficient)| {
                let mut powers = vec![0; variables];
                for (axis, target) in targets.iter().enumerate() {
                    powers[*target] = monomial.exponents()[axis];
                }
                (powers, coefficient)
            })
            .collect(),
    )
}

impl PatchPairIntersection {
    /// Ordered graph pieces of one complete rational-patch branch. No fitted
    /// curve is admitted. Exact-offset inputs require their offset-graph proof.
    pub fn trace_branch(
        &self,
        branch: usize,
        dec: &mut Decider,
    ) -> Result<Vec<PairTracePiece>, Refusal> {
        if self.distance.is_some() || self.second_distance.is_some() {
            return Err(Refusal::new(
                "intersect.patch-pair-offset-trace",
                "offset branches need a retained offset graph certificate",
            ));
        }
        let branch = self.branches.get(branch).ok_or_else(|| {
            Refusal::new(
                "intersect.patch-pair-trace-domain",
                "branch index is outside the complete patch-pair meeting",
            )
        })?;
        let definition = Arc::new(Definition {
            first: self.first.clone(),
            second: self.second.clone(),
            branch: branch.clone(),
        });
        let mut cache = BTreeMap::new();
        let mut pieces = Vec::with_capacity(branch.arcs.len());
        for (index, arc) in branch.arcs.iter().enumerate() {
            let first = definition.first.span(arc.first_chart, arc.first_span);
            let second = definition.second.span(arc.second_chart, arc.second_span);
            let key = (
                arc.first_chart,
                arc.first_span,
                arc.second_chart,
                arc.second_span,
            );
            let equalities: &Arc<Equations> = match cache.entry(key) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    let a = first.coordinates.each_ref().map(|p| {
                        embed(&p.as_rational().expect("pair input is rational"), [0, 1], 5)
                    });
                    let b = second.coordinates.each_ref().map(|p| {
                        embed(&p.as_rational().expect("pair input is rational"), [2, 3], 5)
                    });
                    let mut equations = Vec::with_capacity(3);
                    for k in 0..3 {
                        equations.push(RealPoly::from(difference(
                            &product(&a[k], &b[3], dec)?,
                            &product(&b[k], &a[3], dec)?,
                            dec,
                        )?));
                    }
                    let polynomials: [RealPoly; 3] =
                        equations.try_into().expect("three patch equations");
                    let jacobian = std::array::from_fn(|row| {
                        std::array::from_fn(|column| polynomials[row].partial(column))
                    });
                    entry.insert(Arc::new(Equations {
                        polynomials,
                        jacobian,
                    }))
                }
                std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
            };
            let (coordinate, weight) = match arc.certificate {
                PairCertificate::FieldFirst([dx, dy]) => (
                    &RealPoly::var(5, 1).scale(&Real::from_i64(i64::from(dx)))
                        - &RealPoly::var(5, 0).scale(&Real::from_i64(i64::from(dy))),
                    RealPoly::constant(5, Real::one()),
                ),
                PairCertificate::FieldSecond([dx, dy]) => (
                    &RealPoly::var(5, 2).scale(&Real::from_i64(i64::from(dy)))
                        - &RealPoly::var(5, 3).scale(&Real::from_i64(i64::from(dx))),
                    RealPoly::constant(5, Real::one()),
                ),
                PairCertificate::NormalCone(axis) => {
                    let mut coordinate = RealPoly::zero(5);
                    for k in 0..3 {
                        coordinate = &coordinate
                            + &lift(&first.coordinates[k], &[0, 1], 5)
                                .scale(&Real::from(self.rotation[usize::from(axis)][k].clone()));
                    }
                    (coordinate, lift(&first.coordinates[3], &[0, 1], 5))
                }
            };
            let at = |end: usize| {
                let mut p = arc.ends[end].parameters().to_vec();
                p.push(Real::zero());
                coordinate.eval_real(&p) / weight.eval_real(&p)
            };
            let start = at(0);
            let end = at(1);
            if exact_sign(&(&end - &start), dec)? == Sign::Zero {
                return Err(failure(
                    "the cell certificate lost its strict graph coordinate",
                ));
            }
            let parameter = &RealPoly::constant(5, start.clone())
                + &RealPoly::var(5, 4).scale(&(&end - &start));
            let parameter = &coordinate - &(&parameter * &weight);
            let parameter_jacobian = std::array::from_fn(|column| parameter.partial(column));
            pieces.push(PairTracePiece {
                definition: definition.clone(),
                arc: index,
                system: equalities.clone(),
                parameter,
                parameter_jacobian,
            });
        }
        Ok(pieces)
    }
}

impl PairTracePiece {
    fn equation(&self, row: usize) -> &RealPoly {
        if row < 3 {
            &self.system.polynomials[row]
        } else {
            &self.parameter
        }
    }
    fn derivative(&self, row: usize, column: usize) -> &RealPoly {
        if row < 3 {
            &self.system.jacobian[row][column]
        } else {
            &self.parameter_jacobian[column]
        }
    }
    fn jacobian_at(&self, at: &[Ball; 5], dec: &Decider) -> Result<[[Ball; 4]; 4], Refusal> {
        let row = |r| -> Result<[Ball; 4], Refusal> {
            Ok([
                self.derivative(r, 0).eval_ball(at, dec)?,
                self.derivative(r, 1).eval_ball(at, dec)?,
                self.derivative(r, 2).eval_ball(at, dec)?,
                self.derivative(r, 3).eval_ball(at, dec)?,
            ])
        };
        Ok([row(0)?, row(1)?, row(2)?, row(3)?])
    }
    fn arc(&self) -> &PairArc {
        &self.definition.branch.arcs[self.arc]
    }
    fn first(&self) -> &Span {
        let arc = self.arc();
        self.definition.first.span(arc.first_chart, arc.first_span)
    }
    fn parameters(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 4], Refusal> {
        let bits = t.precision();
        if bits > dec.policy().max_precision_bits
            || t.lo() < Rational::zero()
            || t.hi() > Rational::one()
        {
            return Err(failure(
                "the graph enclosure is outside its finite domain or precision policy",
            ));
        }
        let arc = self.arc();
        if t.lo() == t.hi() {
            for end in 0..2 {
                if t.lo() == Rational::from_i64(end as i64) {
                    let p = arc.ends[end].parameters();
                    return Ok([
                        dec.enclose_at(&p[0], bits)?,
                        dec.enclose_at(&p[1], bits)?,
                        dec.enclose_at(&p[2], bits)?,
                        dec.enclose_at(&p[3], bits)?,
                    ]);
                }
            }
        }
        let domain = [
            arc.first_cell[0].clone(),
            arc.first_cell[1].clone(),
            arc.second_cell[0].clone(),
            arc.second_cell[1].clone(),
        ];
        let mut at: [Ball; 4] = domain
            .each_ref()
            .map(|[lo, hi]| Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits)));
        for _ in 0..bits {
            dec.check_deadline()?;
            let middle: [Ball; 5] = std::array::from_fn(|axis| {
                if axis < 4 {
                    Ball::from_rational(&at[axis].mid(), bits)
                } else {
                    t.clone()
                }
            });
            let all: [Ball; 5] = std::array::from_fn(|axis| {
                if axis < 4 {
                    at[axis].clone()
                } else {
                    t.clone()
                }
            });
            let residual = [
                self.equation(0).eval_ball(&middle, dec)?,
                self.equation(1).eval_ball(&middle, dec)?,
                self.equation(2).eval_ball(&middle, dec)?,
                self.equation(3).eval_ball(&middle, dec)?,
            ];
            let Some(factor) = Factor::new(self.jacobian_at(&all, dec)?, dec)? else {
                break;
            };
            let Some(delta) = factor.solve(residual, dec)? else {
                break;
            };
            let mut next = at.clone();
            for axis in 0..4 {
                let image = &middle[axis] - &delta[axis];
                let lo = at[axis].lo().max(image.lo());
                let hi = at[axis].hi().min(image.hi());
                if lo > hi {
                    return Err(failure(
                        "interval Newton contraction lost its certified graph",
                    ));
                }
                next[axis] = Ball::from_rational(&lo, bits).hull(&Ball::from_rational(&hi, bits));
            }
            if next
                .iter()
                .zip(&at)
                .all(|(a, b)| a.lo() == b.lo() && a.hi() == b.hi())
            {
                break;
            }
            at = next;
        }
        Ok(at)
    }

    /// Exact point on the selected four-parameter graph, with shared system roots.
    pub fn point_at(&self, t: &Rational, dec: &mut Decider) -> Result<Point3, Refusal> {
        if t < &Rational::zero() || t > &Rational::one() {
            return Err(failure("graph point is outside [0, 1]"));
        }
        let arc = self.arc();
        for end in 0..2 {
            if t == &Rational::from_i64(end as i64) {
                let p = arc.ends[end].parameters();
                return Ok(self.first().point(&[p[0].clone(), p[1].clone()]));
            }
        }
        let bounds = [
            arc.first_cell[0].clone(),
            arc.first_cell[1].clone(),
            arc.second_cell[0].clone(),
            arc.second_cell[1].clone(),
        ];
        let equations: Vec<_> = (0..4)
            .map(|row| {
                let p = self.equation(row);
                let fixed = p.substitute(4, &RealPoly::constant(5, Real::from(t.clone())));
                RealPoly::from_terms(
                    4,
                    fixed
                        .terms()
                        .into_iter()
                        .map(|(m, c)| (m.exponents()[..4].to_vec(), c))
                        .collect(),
                )
            })
            .collect();
        let at = RealPoly::simple_algebraic_root_in(
            &equations,
            &bounds,
            128.min(dec.policy().max_precision_bits),
            dec,
        )?;
        Ok(self.first().point(&[at[0].clone(), at[1].clone()]))
    }

    /// Whole-panel Taylor enclosures of the graph, not derivatives of a fit.
    pub fn taylor(&self, t: &Ball, order: usize, dec: &Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        let at = self.parameters(t, dec)?;
        let bits = t.precision();
        let zero = || Ball::from_rational(&Rational::zero(), bits);
        let mut parameters: [Vec<Ball>; 5] = std::array::from_fn(|axis| {
            let mut row = vec![zero(); order + 1];
            row[0] = if axis < 4 {
                at[axis].clone()
            } else {
                t.clone()
            };
            row
        });
        if order > 0 {
            parameters[4][1] = Ball::from_rational(&Rational::one(), bits);
            let all: [Ball; 5] = std::array::from_fn(|axis| {
                if axis < 4 {
                    at[axis].clone()
                } else {
                    t.clone()
                }
            });
            let factor = Factor::new(self.jacobian_at(&all, dec)?, dec)?.ok_or_else(|| {
                Refusal::new(
                    "intersect.patch-pair-trace-jacobian",
                    "the graph panel needs a strict interval Jacobian enclosure",
                )
            })?;
            for n in 1..=order {
                let coefficient = |row| -> Result<Ball, Refusal> {
                    Ok(-series(self.equation(row), &parameters, n, bits, dec)?
                        .pop()
                        .expect("nonempty series"))
                };
                let residual = [
                    coefficient(0)?,
                    coefficient(1)?,
                    coefficient(2)?,
                    coefficient(3)?,
                ];
                let values = factor.solve(residual, dec)?.ok_or_else(|| {
                    Refusal::new(
                        "intersect.patch-pair-trace-jacobian",
                        "the graph panel needs a strict interval Jacobian enclosure",
                    )
                })?;
                for axis in 0..4 {
                    parameters[axis][n] = values[axis].clone();
                }
            }
        }
        let homogeneous: Vec<_> = self
            .first()
            .coordinates
            .iter()
            .map(|p| series(p, &parameters[..2], order, bits, dec))
            .collect::<Result<_, _>>()?;
        let mut result: Vec<[Ball; 3]> = Vec::with_capacity(order + 1);
        for n in 0..=order {
            let mut xyz = std::array::from_fn(|_| zero());
            for axis in 0..3 {
                let mut numerator = homogeneous[axis][n].clone();
                for k in 1..=n {
                    numerator = &numerator - &(&homogeneous[3][k] * &result[n - k][axis]);
                }
                xyz[axis] = numerator.checked_div(&homogeneous[3][0]).ok_or_else(|| {
                    failure("graph homogeneous denominator needs a positive enclosure")
                })?;
            }
            result.push(xyz);
        }
        Ok(result)
    }
}

/// One interval factorization for all Taylor right-hand sides on a panel.
/// Strict pivots cover every matrix in the original interval box.
struct Factor {
    upper: [[Ball; 4]; 4],
    multipliers: [[Option<Ball>; 4]; 4],
    pivots: [usize; 4],
}

impl Factor {
    // Prefer the narrowest relative interval, not the first nonzero entry.
    // In particular, keep a pinned graph-coordinate row ahead of a wide
    // nonlinear surface row. Cross products compare ratios exactly.
    fn pivot(a: &[[Ball; 4]; 4], column: usize) -> Option<usize> {
        let mut best: Option<(usize, Rational, Rational)> = None;
        for (row, values) in a.iter().enumerate().skip(column) {
            let value = &values[column];
            if value.contains_zero() {
                continue;
            }
            let lo = value.lo();
            let hi = value.hi();
            let width = &hi - &lo;
            let margin = if lo.sign() == Sign::Positive {
                lo
            } else {
                -&hi
            };
            if best
                .as_ref()
                .is_none_or(|(_, old_width, old_margin)| &width * old_margin < old_width * &margin)
            {
                best = Some((row, width, margin));
            }
        }
        best.map(|(row, _, _)| row)
    }

    fn new(mut a: [[Ball; 4]; 4], dec: &Decider) -> Result<Option<Self>, Refusal> {
        let mut multipliers: [[Option<Ball>; 4]; 4] =
            std::array::from_fn(|_| std::array::from_fn(|_| None));
        let mut pivots = [0; 4];
        for column in 0..4 {
            dec.check_deadline()?;
            let Some(pivot) = Self::pivot(&a, column) else {
                return Ok(None);
            };
            pivots[column] = pivot;
            a.swap(column, pivot);
            let (leading, trailing) = a.split_at_mut(column + 1);
            let pivot_row = &leading[column];
            for (offset, values) in trailing.iter_mut().enumerate() {
                let Some(scale) = values[column].checked_div(&pivot_row[column]) else {
                    return Ok(None);
                };
                for (target, source) in values.iter_mut().zip(pivot_row).skip(column + 1) {
                    *target = &*target - &(&scale * source);
                }
                multipliers[column + 1 + offset][column] = Some(scale);
            }
        }
        Ok(Some(Self {
            upper: a,
            multipliers,
            pivots,
        }))
    }

    fn solve(&self, mut b: [Ball; 4], dec: &Decider) -> Result<Option<[Ball; 4]>, Refusal> {
        for column in 0..4 {
            dec.check_deadline()?;
            b.swap(column, self.pivots[column]);
            for row in column + 1..4 {
                let scale = self.multipliers[row][column]
                    .as_ref()
                    .expect("retained elimination step");
                b[row] = &b[row] - &(scale * &b[column]);
            }
        }
        for row in (0..4).rev() {
            let mut value = b[row].clone();
            for (coefficient, solved) in self.upper[row].iter().zip(&b).skip(row + 1) {
                value = &value - &(coefficient * solved);
            }
            let Some(value) = value.checked_div(&self.upper[row][row]) else {
                return Ok(None);
            };
            b[row] = value;
        }
        Ok(Some(b))
    }
}
