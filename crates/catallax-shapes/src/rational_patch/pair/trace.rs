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

/// One producer-owned regular graph in four patch parameters. Its fifth
/// polynomial variable is the unit trace parameter, not another unknown.
/// The selected branch and its complete cell-pair certificate remain retained.
#[derive(Clone, Debug)]
pub struct PairTracePiece {
    definition: Arc<Definition>,
    arc: usize,
    equations: [RealPoly; 4],
    jacobian: [[RealPoly; 4]; 4],
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
            let equalities: &[RealPoly; 3] = match cache.entry(key) {
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
                    entry.insert(equations.try_into().expect("three patch equations"))
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
            let equations = [
                equalities[0].clone(),
                equalities[1].clone(),
                equalities[2].clone(),
                &coordinate - &(&parameter * &weight),
            ];
            let jacobian = std::array::from_fn(|row| {
                std::array::from_fn(|column| equations[row].partial(column))
            });
            pieces.push(PairTracePiece {
                definition: definition.clone(),
                arc: index,
                equations,
                jacobian,
            });
        }
        Ok(pieces)
    }
}

impl PairTracePiece {
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
            let mut middle: Vec<Ball> = at
                .iter()
                .map(|x| Ball::from_rational(&x.mid(), bits))
                .collect();
            middle.push(t.clone());
            let mut all = at.to_vec();
            all.push(t.clone());
            let mut matrix: [[Ball; 4]; 4] = std::array::from_fn(|_| {
                std::array::from_fn(|_| Ball::from_rational(&Rational::zero(), bits))
            });
            let mut residual: [Ball; 4] =
                std::array::from_fn(|_| Ball::from_rational(&Rational::zero(), bits));
            for row in 0..4 {
                residual[row] = self.equations[row].eval_ball(&middle, dec)?;
                for column in 0..4 {
                    matrix[row][column] = self.jacobian[row][column].eval_ball(&all, dec)?;
                }
            }
            let Some(delta) = solve(matrix, residual, dec)? else {
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
        let equations: Vec<_> = self
            .equations
            .iter()
            .map(|p| {
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
            let mut all = at.to_vec();
            all.push(t.clone());
            let mut matrix: [[Ball; 4]; 4] =
                std::array::from_fn(|_| std::array::from_fn(|_| zero()));
            for (row, line) in matrix.iter_mut().enumerate() {
                for (column, value) in line.iter_mut().enumerate() {
                    *value = self.jacobian[row][column].eval_ball(&all, dec)?;
                }
            }
            for n in 1..=order {
                let mut residual: [Ball; 4] = std::array::from_fn(|_| zero());
                for (row, value) in residual.iter_mut().enumerate() {
                    *value = -series(&self.equations[row], &parameters, n, bits, dec)?
                        .pop()
                        .unwrap();
                }
                let values = solve(matrix.clone(), residual, dec)?.ok_or_else(|| {
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

/// Interval Gaussian elimination encloses every system in the matrix box.
/// A pivot is divided only after its whole interval excludes zero.
fn solve(
    mut a: [[Ball; 4]; 4],
    mut b: [Ball; 4],
    dec: &Decider,
) -> Result<Option<[Ball; 4]>, Refusal> {
    for column in 0..4 {
        dec.check_deadline()?;
        let Some(pivot) = (column..4).find(|row| !a[*row][column].contains_zero()) else {
            return Ok(None);
        };
        a.swap(column, pivot);
        b.swap(column, pivot);
        for row in column + 1..4 {
            let Some(scale) = a[row][column].checked_div(&a[column][column]) else {
                return Ok(None);
            };
            for k in column + 1..4 {
                a[row][k] = &a[row][k] - &(&scale * &a[column][k]);
            }
            b[row] = &b[row] - &(&scale * &b[column]);
        }
    }
    let mut result = b.clone();
    for row in (0..4).rev() {
        let mut value = b[row].clone();
        for column in row + 1..4 {
            value = &value - &(&a[row][column] * &result[column]);
        }
        let Some(value) = value.checked_div(&a[row][row]) else {
            return Ok(None);
        };
        result[row] = value;
    }
    Ok(Some(result))
}
