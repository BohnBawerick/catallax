//! Numerical offset-spine continuation across retained atlas charts.
//!
//! Newton solves and adaptive quarter-point checks are numerical. The retained
//! tolerance describes those checks, not a complete branch or Hausdorff proof.

use super::jets::{Jet, Meter, Support, ball};
use crate::{BSplineCurve, Point3, RationalAtlas};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use std::sync::Arc;

#[path = "numerical_analytic.rs"]
mod analytic;
#[path = "numerical_corner.rs"]
mod corner;
pub use corner::NumericalFilletCorner;
#[path = "numerical_carrier.rs"]
mod carrier;
pub use carrier::NumericalSplineFillet;
#[path = "numerical_cut.rs"]
mod cut;
#[path = "numerical_extend.rs"]
mod extend;
#[path = "numerical_support.rs"]
mod support;
pub use support::NumericalFilletSupport;
#[path = "numerical_walk.rs"]
mod walk;
pub use cut::NumericalFilletCut;

#[path = "numerical_surface_cut.rs"]
mod surface_cut;
pub use surface_cut::NumericalSurfaceCut;

const MAX_DEPTH: u32 = 24;
const MAX_NODES: usize = 16_385;

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("fillet.numerical-spine", detail)
}
fn rational(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| bad("non-finite spine coordinate"))
}
fn norm(value: [f64; 3]) -> f64 {
    value.iter().map(|v| v * v).sum::<f64>().sqrt()
}
fn inverse<const N: usize>(mut a: [[f64; N]; N], mut b: [f64; N]) -> Option<[f64; N]> {
    for k in 0..N {
        let pivot = (k..N).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))?;
        if !a[pivot][k].is_finite() || a[pivot][k].abs() < f64::MIN_POSITIVE {
            return None;
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        let divisor = a[k][k];
        for j in k..N {
            a[k][j] /= divisor;
        }
        b[k] /= divisor;
        for i in 0..N {
            if i == k {
                continue;
            }
            let scale = a[i][k];
            for j in k..N {
                a[i][j] -= scale * a[k][j];
            }
            b[i] -= scale * b[k];
        }
    }
    b.iter().all(|v| v.is_finite()).then_some(b)
}

/// One checked numerical station. Chart indices refer to the input atlases.
#[derive(Clone, Debug)]
pub struct NumericalSpineNode {
    /// Parameter of the first support used to march the spine.
    pub tau: Rational,
    /// Contact parameters on each original support.
    pub parameters: [[Rational; 2]; 2],
    /// Retained chart owning each contact.
    pub charts: [usize; 2],
    /// Centre checked against both signed offsets at this station.
    pub centre: [Rational; 3],
}

/// A C1 cubic centre curve with sampled numerical checks, never exact authority.
#[derive(Clone, Debug)]
pub struct NumericalFilletSpine {
    problem: Arc<Problem>,
    nodes: Vec<NumericalSpineNode>,
    centre_curve: BSplineCurve,
    tolerance_mm: Rational,
    sampled_error_mm: Rational,
    samples: usize,
}

#[derive(Debug)]
enum SpineGuide {
    Coordinate(usize),
    Curve(BSplineCurve),
}

#[derive(Debug)]
struct Problem {
    supports: [NumericalFilletSupport; 2],
    guide: SpineGuide,
    tolerance: f64,
}
impl Problem {
    fn params<T: Clone>(&self, tau: T, x: &[T; 3]) -> [[T; 2]; 2] {
        let SpineGuide::Coordinate(axis) = self.guide else {
            unreachable!("coordinate solver needs a coordinate guide")
        };
        let first = if axis == 0 {
            [tau, x[0].clone()]
        } else {
            [x[0].clone(), tau]
        };
        [first, [x[1].clone(), x[2].clone()]]
    }
    fn residual(&self, tau: f64, x: [f64; 3], dec: &mut Decider) -> Result<[f64; 3], Refusal> {
        let at = self.params(tau, &x);
        let a = self.supports[0].centre(at[0], dec)?;
        let b = self.supports[1].centre(at[1], dec)?;
        Ok(std::array::from_fn(|i| a[i] - b[i]))
    }
    fn solve(
        &self,
        tau: &Rational,
        seed: [f64; 4],
        dec: &mut Decider,
    ) -> Result<NumericalSpineNode, Refusal> {
        match &self.guide {
            SpineGuide::Coordinate(axis) => {
                self.solve_axis(tau, [seed[1 - axis], seed[2], seed[3]], dec)
            }
            SpineGuide::Curve(curve) => self.solve_guided(curve, tau, seed, dec),
        }
    }

    fn solve_axis(
        &self,
        tau: &Rational,
        mut x: [f64; 3],
        dec: &mut Decider,
    ) -> Result<NumericalSpineNode, Refusal> {
        let t = tau.to_f64_nearest();
        for _ in 0..40 {
            dec.check_deadline()?;
            let residual = self.residual(t, x, dec)?;
            let error = norm(residual);
            if error <= self.tolerance / 64.0 {
                let unknowns = [rational(x[0])?, rational(x[1])?, rational(x[2])?];
                let parameters = self.params(tau.clone(), &unknowns);
                let (ca, a) = self.supports[0].checked(&parameters[0], dec)?;
                let (cb, b) = self.supports[1].checked(&parameters[1], dec)?;
                let centre = std::array::from_fn(|i| (&a[i].mid() + &b[i].mid()).mul_pow2(-1));
                let upper = offset_error(&centre, [&a, &b]);
                if upper.to_f64_nearest() > self.tolerance / 8.0 {
                    return Err(bad(
                        "offset residual enclosure exceeds the stated numerical tolerance",
                    ));
                }
                return Ok(NumericalSpineNode {
                    tau: tau.clone(),
                    parameters,
                    charts: [ca, cb],
                    centre,
                });
            }
            let mut jacobian = [[0.0; 3]; 3];
            for j in 0..3 {
                let h = 1e-6 * x[j].abs().max(1.0);
                let mut up = x;
                let mut down = x;
                up[j] += h;
                down[j] -= h;
                let plus = self.residual(t, up, dec);
                let minus = self.residual(t, down, dec);
                let (difference, width): ([f64; 3], f64) = match (plus, minus) {
                    (Ok(a), Ok(b)) => (std::array::from_fn(|i| a[i] - b[i]), 2.0 * h),
                    (Ok(a), Err(e)) if e.code == "fillet.numerical-spine" => {
                        (std::array::from_fn(|i| a[i] - residual[i]), h)
                    }
                    (Err(e), Ok(b)) if e.code == "fillet.numerical-spine" => {
                        (std::array::from_fn(|i| residual[i] - b[i]), h)
                    }
                    (Err(e), _) | (_, Err(e)) => return Err(e),
                };
                for i in 0..3 {
                    jacobian[i][j] = difference[i] / width;
                }
            }
            dec.charge_work(64)?;
            let step = inverse(jacobian, residual)
                .ok_or_else(|| bad("singular numerical offset-spine Jacobian"))?;
            let mut accepted = None;
            for power in 0..16 {
                let share = 2.0f64.powi(-power);
                let next = std::array::from_fn(|i| x[i] - share * step[i]);
                match self.residual(t, next, dec) {
                    Ok(value) if norm(value) < error => {
                        accepted = Some(next);
                        break;
                    }
                    Ok(_) => {}
                    Err(e) if e.code == "fillet.numerical-spine" => {}
                    Err(e) => return Err(e),
                }
            }
            x = accepted.ok_or_else(|| {
                bad("numerical offset-spine Newton step did not reduce the residual")
            })?;
        }
        Err(bad("numerical offset-spine Newton iteration limit reached"))
    }
    fn guess(&self, a: &NumericalSpineNode, b: &NumericalSpineNode, share: f64) -> [f64; 4] {
        let extract = |n: &NumericalSpineNode| {
            std::array::from_fn::<_, 4, _>(|i| n.parameters[i / 2][i % 2].to_f64_nearest())
        };
        let (a, b) = (extract(a), extract(b));
        std::array::from_fn(|i| a[i] + share * (b[i] - a[i]))
    }
}
fn offset_error(centre: &[Rational; 3], offsets: [&[Ball; 3]; 2]) -> Rational {
    offsets
        .into_iter()
        .map(|offset| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&offset[i].lo() - &centre[i])
                    .abs()
                    .max((&offset[i].hi() - &centre[i]).abs())
            })
        })
        .max()
        .expect("two offsets")
}

fn centre_fit(nodes: &[NumericalSpineNode], dec: &mut Decider) -> Result<BSplineCurve, Refusal> {
    let (knots, controls) = carrier::hermite(nodes, |i| &nodes[i].centre);
    let weights = vec![Real::one(); controls.len()];
    let points = controls
        .into_iter()
        .map(|[x, y, z]| {
            Point3::new(
                Real::from_rational(x),
                Real::from_rational(y),
                Real::from_rational(z),
            )
        })
        .collect();
    BSplineCurve::new(3, knots, points, weights, dec)
}

impl NumericalFilletSpine {
    /// Marches one selected branch as a graph over `axis` of the first support.
    /// The same support adapters serve offset walks and zero-offset sharp guides.
    /// The seed supplies its other parameter and both parameters of the second.
    /// Quarter-point position checks adaptively split each span. They do not
    /// establish branch completeness, uniqueness, or deviation between samples.
    pub fn trace(
        supports: [NumericalFilletSupport; 2],
        axis: usize,
        range: [Rational; 2],
        seed: [Rational; 3],
        tolerance_mm: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if axis > 1 || range[0] >= range[1] || tolerance_mm.sign() != Sign::Positive {
            return Err(bad(
                "spine needs a valid axis, increasing range and positive tolerance",
            ));
        }
        let tolerance = tolerance_mm.to_f64_nearest();
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(bad("spine tolerance is outside binary64 range"));
        }
        let problem = Problem {
            supports,
            guide: SpineGuide::Coordinate(axis),
            tolerance,
        };
        let parameters = problem.params(range[0].clone(), &seed);
        Self::trace_problem(
            problem,
            range,
            std::array::from_fn(|i| parameters[i / 2][i % 2].to_f64_nearest()),
            tolerance_mm,
            dec,
        )
    }

    /// Walks the intersection of signed offsets in planes normal to `guide`.
    /// Contact, guide-plane and fitted-quarter checks are numerical. A failed
    /// walk returns its cause rather than replacing the selected branch.
    /// Outside the guide's finite interval its endpoint tangent lines define
    /// the sections. Those tails select the walk, not an extrapolated support.
    pub fn walk(
        supports: [NumericalFilletSupport; 2],
        guide: BSplineCurve,
        range: [Rational; 2],
        seed: [[Rational; 2]; 2],
        tolerance_mm: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let tolerance = tolerance_mm.to_f64_nearest();
        if range[0] >= range[1]
            || tolerance_mm.sign() != Sign::Positive
            || !tolerance.is_finite()
            || tolerance <= 0.0
        {
            return Err(bad(
                "walking needs an increasing range and finite positive tolerance",
            ));
        }
        Self::trace_problem(
            Problem {
                supports,
                guide: SpineGuide::Curve(guide),
                tolerance,
            },
            range,
            std::array::from_fn(|i| seed[i / 2][i % 2].to_f64_nearest()),
            tolerance_mm,
            dec,
        )
    }

    fn trace_problem(
        problem: Problem,
        range: [Rational; 2],
        seed: [f64; 4],
        tolerance_mm: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let first = problem.solve(&range[0], seed, dec)?;
        let mut nodes = vec![first];
        for i in 1..=8 {
            let tau = &range[0] + &((&range[1] - &range[0]) * Rational::from_ratio(i, 8).unwrap());
            let previous = nodes.last().expect("first station");
            let guess = problem.guess(previous, previous, 0.0);
            nodes.push(problem.solve(&tau, guess, dec)?);
        }
        let centre_curve = centre_fit(&nodes, dec)?;
        let mut spine = Self {
            problem: Arc::new(problem),
            samples: nodes.len(),
            nodes,
            centre_curve,
            tolerance_mm,
            sampled_error_mm: Rational::zero(),
        };
        for round in 0..=MAX_DEPTH {
            let mut inserted = Vec::new();
            let mut upper = Rational::zero();
            for (index, pair) in spine.nodes.windows(2).enumerate() {
                let mut worst = Rational::zero();
                let mut middle = None;
                for i in 1..4 {
                    let share = Rational::from_ratio(i, 4).unwrap();
                    let tau = &pair[0].tau + &(&(&pair[1].tau - &pair[0].tau) * &share);
                    let checked = spine.problem.solve(
                        &tau,
                        spine.problem.guess(&pair[0], &pair[1], i as f64 / 4.0),
                        dec,
                    )?;
                    worst = worst.max(spine.centre_error(&checked, dec)?);
                    spine.samples += 1;
                    if i == 2 {
                        middle = Some(checked);
                    }
                }
                if worst > spine.tolerance_mm.mul_pow2(-2) {
                    inserted.push((index + 1, middle.expect("checked midpoint")));
                }
                upper = upper.max(worst);
            }
            if inserted.is_empty() {
                spine.sampled_error_mm = upper;
                return Ok(spine);
            }
            if round == MAX_DEPTH || spine.nodes.len() + inserted.len() > MAX_NODES {
                return Err(bad(
                    "numerical cubic spine cannot meet its sampled tolerance within the subdivision limit",
                ));
            }
            for (index, node) in inserted.into_iter().rev() {
                spine.nodes.insert(index, node);
            }
            spine.centre_curve = centre_fit(&spine.nodes, dec)?;
        }
        unreachable!("trace returns at its subdivision limit")
    }
    /// Ordered stations, including both ends and every accepted span boundary.
    pub fn nodes(&self) -> &[NumericalSpineNode] {
        &self.nodes
    }
    /// Declared numerical tolerance in millimetres, not a global bound.
    pub fn tolerance_mm(&self) -> &Rational {
        &self.tolerance_mm
    }
    /// Largest checked quarter-point error from the cubic centre curve, in L1 mm.
    pub fn sampled_error_mm(&self) -> &Rational {
        &self.sampled_error_mm
    }
    /// Number of Newton stations checked, including rejected coarse spans.
    pub fn samples(&self) -> usize {
        self.samples
    }
    /// The C1 cubic reconstruction checked against the offset-centre stations.
    pub fn centre_curve(&self) -> &BSplineCurve {
        &self.centre_curve
    }

    fn centre_error(
        &self,
        node: &NumericalSpineNode,
        dec: &mut Decider,
    ) -> Result<Rational, Refusal> {
        let fitted = self
            .centre_curve
            .evaluate_box([&node.tau, &node.tau], 128, dec)?;
        Ok((0..3).fold(Rational::zero(), |sum, axis| {
            sum + (&fitted[axis].lo() - &node.centre[axis])
                .abs()
                .max((&fitted[axis].hi() - &node.centre[axis]).abs())
        }))
    }
}
