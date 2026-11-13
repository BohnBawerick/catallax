//! A tensor grid between two certified cut graphs.

use super::super::rational_patch::Grid;
use super::*;

fn ball(q: &Rational) -> Ball {
    Ball::from_rational(q, BITS)
}
fn interval(a: &Rational, b: &Rational) -> Ball {
    ball(a).hull(&ball(b))
}
fn magnitude(b: &Ball) -> Rational {
    b.lo().abs().max(b.hi().abs())
}
fn l1(v: &[Ball; 3]) -> Rational {
    v.iter().fold(Rational::zero(), |sum, b| sum + magnitude(b))
}
#[derive(Clone)]
pub(in super::super) struct Mapping {
    position: Box<super::super::position::PositionSurface>,
    cuts: [CutCurve; 2],
    planes: [([Rational; 3], Rational); 2],
    boundaries: [BSplineCurve; 4],
    deviation: Rational,
}

impl Mapping {
    pub(super) fn new(
        fillet: SplineFillet,
        rails: [BSplineCurve; 2],
        cuts: [CutCurve; 2],
        planes: [([Rational; 3], Rational); 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let low = cuts[0].tau.iter().max().expect("cut nodes");
        let high = cuts[1].tau.iter().min().expect("cut nodes");
        let shift = cuts
            .iter()
            .map(|cut| cut.parameter_upper.clone())
            .max()
            .unwrap();
        if high - low <= shift.mul_pow2(1) {
            return Err(bad(format!(
                "cut graph separation {} does not exceed parameter uncertainty {}",
                high - low,
                shift.mul_pow2(1)
            )));
        }
        let normalized = |rail: &BSplineCurve, dec: &mut Decider| {
            let [lo, hi] = rail.domain();
            let width = hi - lo;
            BSplineCurve::new(
                rail.degree(),
                rail.knots().iter().map(|t| &(t - lo) / &width).collect(),
                rail.controls().cloned().collect(),
                rail.weights().cloned().collect(),
                dec,
            )
        };
        let cut_curve = |cut: &CutCurve, dec: &mut Decider| {
            let knots = std::iter::repeat_n(cut.theta[0].clone(), cut.degree + 1)
                .chain(
                    cut.theta[1..cut.theta.len() - 1]
                        .iter()
                        .flat_map(|t| std::iter::repeat_n(t.clone(), cut.degree)),
                )
                .chain(std::iter::repeat_n(
                    cut.theta.last().unwrap().clone(),
                    cut.degree + 1,
                ))
                .collect();
            BSplineCurve::new(
                cut.degree,
                knots,
                cut.controls.iter().map(|(p, _)| p.clone()).collect(),
                cut.controls
                    .iter()
                    .map(|(_, w)| Real::from(w.clone()))
                    .collect(),
                dec,
            )
        };
        let boundaries = [
            normalized(&rails[0], dec)?,
            normalized(&rails[1], dec)?,
            cut_curve(&cuts[0], dec)?,
            cut_curve(&cuts[1], dec)?,
        ];
        let carrier = fillet.carrier();
        let position = Box::new(super::super::position::PositionSurface::new(carrier, dec)?);
        let domain = carrier.domain().map(|a| a.map(Clone::clone));
        let reference = std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
        let origin = std::array::from_fn(|_| Rational::zero());
        let first = position.derivative_residual(&domain, &origin, &reference, dec)?;
        let deviation = (&shift * &l1(&first[0])) / (&domain[0][1] - &domain[0][0]);
        Ok(Self {
            position,
            cuts,
            planes,
            boundaries,
            deviation,
        })
    }

    pub(in super::super) fn boundary(&self, row: BoundaryRow) -> BSplineCurve {
        self.boundaries[match row {
            BoundaryRow::Rail => 0,
            BoundaryRow::Mouth => 1,
            BoundaryRow::Negative => 2,
            BoundaryRow::Positive => 3,
        }]
        .clone()
    }

    pub(in super::super) fn exact_row(
        &self,
        row: BoundaryRow,
        parameters: &[Real],
        dec: &mut Decider,
    ) -> Result<Vec<[Real; 3]>, Refusal> {
        let curve = self.boundary(row);
        parameters
            .iter()
            .map(|t| {
                let p = curve.point_at(t, dec)?;
                let mut xyz = [p.x, p.y, p.z];
                if let Some(end) = match row {
                    BoundaryRow::Negative => Some(0),
                    BoundaryRow::Positive => Some(1),
                    _ => None,
                } {
                    let (n, c) = &self.planes[end];
                    let n2 = n.iter().fold(Rational::zero(), |sum, n| sum + n * n);
                    if n2 <= Rational::zero() {
                        return Err(bad("a cap projection has zero normal"));
                    }
                    let distance = (0..3).fold(-Real::from(c.clone()), |sum, i| {
                        &sum + &(&Real::from(n[i].clone()) * &xyz[i])
                    });
                    xyz =
                        std::array::from_fn(|i| &xyz[i] - &(&distance * &Real::from(&n[i] / &n2)));
                }
                Ok(xyz)
            })
            .collect()
    }

    pub(in super::super) fn position_grid(
        &self,
        tolerance: &Rational,
        max_vertices: usize,
        required: [&[Real]; 2],
        dec: &mut Decider,
    ) -> Result<Grid, Refusal> {
        let surface = &self.position;
        let mut axes = [
            vec![Real::from_i64(0), Real::from_i64(1)],
            vec![Real::from_i64(0), Real::from_i64(1)],
        ];
        for curve in &self.boundaries[..2] {
            for knot in curve.knots() {
                super::super::rational_patch::insert(&mut axes[0], Real::from(knot.clone()), dec)?;
            }
        }
        for cut in &self.cuts {
            for knot in &cut.theta {
                super::super::rational_patch::insert(&mut axes[1], Real::from(knot.clone()), dec)?;
            }
        }
        surface.grid(
            axes,
            required,
            super::super::position::Controls {
                tolerance,
                max_vertices,
                displacement: &self.deviation,
                cartesian: false,
            },
            |at| {
                let graph = |cut: usize| -> Result<Ball, Refusal> {
                    let curve = &self.cuts[cut];
                    let mut result: Option<Ball> = None;
                    for i in 0..curve.theta.len() - 1 {
                        let lo = at[1].lo().max(curve.theta[i].clone());
                        let hi = at[1].hi().min(curve.theta[i + 1].clone());
                        if lo > hi {
                            continue;
                        }
                        let slope = (&curve.tau[i + 1] - &curve.tau[i])
                            / (&curve.theta[i + 1] - &curve.theta[i]);
                        let value = &ball(&curve.tau[i])
                            + &(&ball(&slope) * &(&interval(&lo, &hi) - &ball(&curve.theta[i])));
                        result = Some(result.map_or(value.clone(), |old| old.hull(&value)));
                    }
                    result.ok_or_else(|| bad("a measured cut station leaves its complete graph"))
                };
                let low = graph(0)?;
                let high = graph(1)?;
                Ok([&low + &(&at[0] * &(&high - &low)), at[1].clone()])
            },
            dec,
        )
    }
}
