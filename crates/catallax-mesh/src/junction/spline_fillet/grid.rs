//! A tensor grid between two certified cut graphs.

use super::super::rational_patch::{Embedding as GridEmbedding, Grid};
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
fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    (0..3).fold(ball(&Rational::zero()), |sum, i| &sum + &(&a[i] * &b[i]))
}

/// A fixed projection makes the symmetric part of every cell's independent
/// secant-column family positive. This proves the complete PL map is embedded.
pub(in super::super) struct Embedding {
    lower: Rational,
    projection_norms: [Rational; 2],
    steps: [Rational; 2],
}

impl Embedding {
    pub(in super::super) fn accepts_rounding(&self, error: &Rational) -> bool {
        if error < &Rational::zero() {
            return false;
        }
        let twice = error.mul_pow2(1);
        let p = self.projection_norms.each_ref().map(|n| &twice * n);
        let diagonal = (&p[0] / &self.steps[0]).max(&p[1] / &self.steps[1]);
        let off = (&p[0] / &self.steps[1] + &p[1] / &self.steps[0]).mul_pow2(-1);
        self.lower > diagonal + off
    }
}

#[derive(Clone)]
pub(in super::super) struct Mapping {
    fillet: SplineFillet,
    prepared: Vec<BSplineSurface>,
    position: Option<Box<super::super::position::PositionSurface>>,
    cuts: [CutCurve; 2],
    planes: [([Rational; 3], Rational); 2],
    boundaries: [BSplineCurve; 4],
    deviation: Rational,
    projection: [[Ball; 3]; 2],
}

impl Mapping {
    pub(super) fn new(
        fillet: SplineFillet,
        rails: [BSplineCurve; 2],
        cuts: [CutCurve; 2],
        planes: [([Rational; 3], Rational); 2],
        measured: bool,
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
        let mut prepared = Vec::new();
        if !measured {
            for span in carrier.span_indices() {
                let knots = carrier.knots();
                let domain: [[Rational; 2]; 2] = std::array::from_fn(|axis| {
                    [
                        knots[axis][span[axis]].clone(),
                        knots[axis][span[axis] + 1].clone(),
                    ]
                });
                let patch = carrier.bezier_patch(
                    span,
                    domain.each_ref().map(|bounds| bounds.each_ref()),
                    dec,
                )?;
                let first = &patch.controls()[0];
                let last = patch.controls().last().expect("span controls");
                let two = Real::from_i64(2);
                let origin = Point3::new(
                    &(&first.x + &last.x) / &two,
                    &(&first.y + &last.y) / &two,
                    &(&first.z + &last.z) / &two,
                );
                // Derivatives are translation invariant. Keep exact, span-local
                // coefficients so interval quotients do not subtract world origins.
                let controls = patch
                    .controls()
                    .iter()
                    .map(|point| {
                        Point3::new(
                            &point.x - &origin.x,
                            &point.y - &origin.y,
                            &point.z - &origin.z,
                        )
                    })
                    .collect();
                prepared.push(BSplineSurface::new(
                    patch.degrees(),
                    patch.knots().map(|knots| knots.to_vec()),
                    patch.control_count(),
                    controls,
                    patch.weights().to_vec(),
                    dec,
                )?);
            }
        }
        let position = if measured {
            Some(Box::new(super::super::position::PositionSurface::new(
                carrier, dec,
            )?))
        } else {
            None
        };
        let mut mapping = Self {
            fillet,
            prepared,
            position,
            cuts,
            planes,
            boundaries,
            deviation: Rational::zero(),
            projection: std::array::from_fn(|_| std::array::from_fn(|_| ball(&Rational::zero()))),
        };
        let domain = mapping.fillet.carrier().domain();
        if let Some(position) = &mapping.position {
            let domain = domain.map(|a| a.map(Clone::clone));
            let reference = std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
            let origin = std::array::from_fn(|_| Rational::zero());
            let first = position.derivative_residual(&domain, &origin, &reference, dec)?;
            mapping.deviation = (&shift * &l1(&first[0])) / (&domain[0][1] - &domain[0][0]);
            return Ok(mapping);
        }
        let derivatives = mapping.first_bounds(domain[0], domain[1], dec)?;
        mapping.deviation = &shift * &l1(&derivatives[0]);
        let half = Rational::from_ratio(1, 2).unwrap();
        let first = mapping.derivatives([&half, &half], [&half, &half], dec)?;
        // A binary rational proposal only. Every published cell proves the
        // resulting projection inequality with outward interval arithmetic.
        let a = first[0].each_ref().map(|b| b.mid());
        let b = first[1].each_ref().map(|b| b.mid());
        let qdot = |a: &[Rational; 3], b: &[Rational; 3]| {
            (0..3).fold(Rational::zero(), |sum, i| sum + &a[i] * &b[i])
        };
        let aa = qdot(&a, &a);
        let bb = qdot(&b, &b);
        let ab = qdot(&a, &b);
        let determinant = &aa * &bb - &ab * &ab;
        if determinant <= Rational::zero() {
            return Err(bad(format!(
                "mapped chart centre determinant {determinant} is not positive"
            )));
        }
        mapping.projection = [
            std::array::from_fn(|i| ball(&(&(&a[i] * &bb - &b[i] * &ab) / &determinant))),
            std::array::from_fn(|i| ball(&(&(&b[i] * &aa - &a[i] * &ab) / &determinant))),
        ];
        Ok(mapping)
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

    fn graph(&self, cut: usize, theta: &Rational) -> (Rational, Rational) {
        let curve = &self.cuts[cut];
        let i = curve
            .theta
            .partition_point(|t| t <= theta)
            .saturating_sub(1)
            .min(curve.theta.len() - 2);
        let slope = (&curve.tau[i + 1] - &curve.tau[i]) / (&curve.theta[i + 1] - &curve.theta[i]);
        (
            &curve.tau[i] + &(&slope * &(theta - &curve.theta[i])),
            slope,
        )
    }

    fn tau(&self, s: &Rational, theta: &Rational) -> Rational {
        let low = self.graph(0, theta).0;
        let high = self.graph(1, theta).0;
        &low + &(s * &(&high - &low))
    }

    fn first_bounds(
        &self,
        u: [&Rational; 2],
        v: [&Rational; 2],
        dec: &Decider,
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        let mut result = None::<[[Ball; 3]; 2]>;
        for patch in &self.prepared {
            let whole = patch.domain();
            let bounds: [[Rational; 2]; 2] = std::array::from_fn(|axis| {
                [
                    [u, v][axis][0].clone().max(whole[axis][0].clone()),
                    [u, v][axis][1].clone().min(whole[axis][1].clone()),
                ]
            });
            if bounds.iter().any(|b| b[0] > b[1]) {
                continue;
            }
            let first =
                patch.derivatives_box(bounds[0].each_ref(), bounds[1].each_ref(), BITS, dec)?;
            match &mut result {
                None => result = Some(first),
                Some(previous) => {
                    for axis in 0..2 {
                        for coordinate in 0..3 {
                            previous[axis][coordinate] =
                                previous[axis][coordinate].hull(&first[axis][coordinate]);
                        }
                    }
                }
            }
        }
        result.ok_or_else(|| bad("mapped derivative box meets no prepared source span"))
    }

    fn derivatives(
        &self,
        s: [&Rational; 2],
        theta: [&Rational; 2],
        dec: &Decider,
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        let a = theta.map(|t| self.graph(0, t));
        let b = theta.map(|t| self.graph(1, t));
        let taus = [
            self.tau(s[0], theta[0]),
            self.tau(s[0], theta[1]),
            self.tau(s[1], theta[0]),
            self.tau(s[1], theta[1]),
        ];
        let tau = interval(taus.iter().min().unwrap(), taus.iter().max().unwrap());
        let ds = interval(&(&b[0].0 - &a[0].0), &(&b[1].0 - &a[1].0));
        // A theta cell contains no graph breakpoint in its interior.
        // At its right endpoint graph() selects the next piece; use its
        // secant instead, which is this cell's exact slope.
        let slope = if theta[0] == theta[1] {
            [a[0].1.clone(), b[0].1.clone()]
        } else {
            [
                (&a[1].0 - &a[0].0) / (theta[1] - theta[0]),
                (&b[1].0 - &b[0].0) / (theta[1] - theta[0]),
            ]
        };
        let dtheta = interval(
            &(&slope[0] + &(s[0] * &(&slope[1] - &slope[0]))),
            &(&slope[0] + &(s[1] * &(&slope[1] - &slope[0]))),
        );
        let tau_bounds = [tau.lo(), tau.hi()];
        let raw = self.first_bounds(tau_bounds.each_ref(), theta, dec)?;
        Ok([
            std::array::from_fn(|i| &raw[0][i] * &ds),
            std::array::from_fn(|i| &(&raw[0][i] * &dtheta) + &raw[1][i]),
        ])
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
        let surface = self
            .position
            .as_ref()
            .ok_or_else(|| bad("positional evaluator is missing"))?;
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

    pub(in super::super) fn grid(
        &self,
        tolerance: &Rational,
        max_vertices: usize,
        required: [&[Real]; 2],
        dec: &mut Decider,
    ) -> Result<Grid, Refusal> {
        if self.deviation >= tolerance.mul_pow2(-2) {
            return Err(bad(format!(
                "cut-domain displacement {} mm exceeds grid reserve {} mm",
                self.deviation,
                tolerance.mul_pow2(-2)
            )));
        }
        let mut axes = [
            vec![Rational::zero(), Rational::one()],
            vec![Rational::zero(), Rational::one()],
        ];
        for (axis, values) in required.into_iter().enumerate() {
            for value in values {
                let q = value
                    .as_rational()
                    .ok_or_else(|| bad("mapped fillet stations require rational parameters"))?;
                if q < Rational::zero() || q > Rational::one() {
                    return Err(bad(format!("mapped station {q} leaves [0,1]")));
                }
                axes[axis].push(q);
            }
        }
        for curve in &self.boundaries[..2] {
            axes[0].extend(curve.knots().iter().cloned());
        }
        for cut in &self.cuts {
            axes[1].extend(cut.theta.iter().cloned());
        }
        for axis in &mut axes {
            axis.sort();
            axis.dedup();
        }
        let target = (tolerance - &self.deviation).mul_pow2(-1);
        let mut bounds = std::collections::BTreeMap::new();
        let mut certified = std::collections::BTreeMap::<_, (Rational, Rational)>::new();
        loop {
            dec.check_deadline()?;
            let count = axes[0]
                .len()
                .checked_mul(axes[1].len())
                .ok_or_else(|| bad("mapped grid vertex count overflow"))?;
            if count > max_vertices {
                return Err(bad(format!(
                    "mapped grid needs {count} vertices, cap {max_vertices}"
                )));
            }
            let mut worst = Rational::zero();
            let mut lower = None::<Rational>;
            let mut refine = None;
            'cells: for (i, s) in axes[0].windows(2).enumerate() {
                for (j, theta) in axes[1].windows(2).enumerate() {
                    let key = [
                        [s[0].clone(), s[1].clone()],
                        [theta[0].clone(), theta[1].clone()],
                    ];
                    dec.charge_work(1)?;
                    if let Some((chord, margin)) = certified.get(&key) {
                        worst = worst.max(chord.clone());
                        lower = Some(lower.map_or(margin.clone(), |old| old.min(margin.clone())));
                        continue;
                    }
                    if let std::collections::btree_map::Entry::Vacant(entry) =
                        bounds.entry(key.clone())
                    {
                        entry.insert(self.derivatives(
                            [&s[0], &s[1]],
                            [&theta[0], &theta[1]],
                            dec,
                        )?);
                    }
                    let first = &bounds[&key];
                    let widths = [&s[1] - &s[0], &theta[1] - &theta[0]];
                    // Subtract any constant affine map. Its barycentric residual
                    // vanishes. Integrating the remaining derivative boxes bounds
                    // the error in both directions, even across carrier knots.
                    let contributions: [Rational; 2] = std::array::from_fn(|axis| {
                        &first[axis]
                            .iter()
                            .fold(Rational::zero(), |sum, b| sum + b.rad())
                            * &widths[axis]
                    });
                    let chord = &contributions[0] + &contributions[1];
                    let diagonal = [
                        dot(&self.projection[0], &first[0]),
                        dot(&self.projection[1], &first[1]),
                    ];
                    let off = (&dot(&self.projection[0], &first[1])
                        + &dot(&self.projection[1], &first[0]))
                        * ball(&Rational::from_ratio(1, 2).unwrap());
                    let margin = diagonal[0].lo().min(diagonal[1].lo()) - magnitude(&off);
                    if chord > target || margin <= Rational::zero() {
                        let axis = if margin <= Rational::zero() {
                            usize::from(widths[1] > widths[0])
                        } else {
                            usize::from(contributions[1] > contributions[0])
                        };

                        refine = Some((axis, [i, j][axis]));
                        break 'cells;
                    }
                    certified.insert(key, (chord.clone(), margin.clone()));
                    worst = worst.max(chord);
                    lower = Some(lower.map_or(margin.clone(), |old| old.min(margin)));
                }
            }
            if let Some((axis, index)) = refine {
                let next = (axes[axis].len() + 1).checked_mul(axes[1 - axis].len());
                if next.is_none_or(|count| count > max_vertices) {
                    return Err(bad(format!(
                        "mapped grid refinement exceeds vertex cap {max_vertices}"
                    )));
                }
                let midpoint = (&axes[axis][index] + &axes[axis][index + 1]).mul_pow2(-1);
                bounds
                    .retain(|domain, _| domain[axis][1] <= midpoint || domain[axis][0] >= midpoint);
                certified
                    .retain(|domain, _| domain[axis][1] <= midpoint || domain[axis][0] >= midpoint);
                axes[axis].insert(index + 1, midpoint);
                continue;
            }
            let mut points = Vec::with_capacity(count);
            let mut uncertainty = Rational::zero();
            for s in &axes[0] {
                for theta in &axes[1] {
                    let tau = self.tau(s, theta);
                    let p = self
                        .fillet
                        .carrier()
                        .evaluate(&ball(&tau), &ball(theta), dec)?;
                    uncertainty = uncertainty.max(
                        p.iter()
                            .fold(Rational::zero(), |sum, b| sum + b.rad().mul_pow2(1)),
                    );
                    points.push(p);
                }
            }
            let chord = (&worst + &self.deviation) + uncertainty;
            if chord > *tolerance {
                return Err(bad(format!(
                    "mapped grid error {chord} mm exceeds {tolerance} mm"
                )));
            }
            let embedding = Embedding {
                lower: lower.expect("grid cells"),
                projection_norms: self.projection.each_ref().map(l1),
                steps: axes
                    .each_ref()
                    .map(|axis| axis.windows(2).map(|p| &p[1] - &p[0]).min().unwrap()),
            };
            return Ok(Grid {
                grid: EnvelopeGrid {
                    columns: axes[0].len(),
                    rows: axes[1].len(),
                    points,
                    chord,
                },
                parameters: axes.map(|axis| axis.into_iter().map(Real::from).collect()),
                normal_boxes: None,
                embedding: GridEmbedding::Fillet(embedding),
            });
        }
    }
}
