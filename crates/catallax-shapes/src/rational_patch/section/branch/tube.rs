//! Two-sided graph tubes, including their exact end caps.

use super::*;
use catallax_numbers::{GeometricProof, LocalToleranceSettings};

/// Untrusted polynomial cubic Bezier controls in the patch's knot parameters.
/// Success of the complete search never depends on this list being nonempty.
#[derive(Clone, Debug)]
pub struct SectionCandidate {
    /// Each piece has four exact (u,v) controls. Traced binary64 values can be
    /// retained exactly with `Rational::from_f64_exact`.
    pub pieces: Vec<[[Rational; 2]; 4]>,
}

/// A checked representation with a symmetric whole-component distance bound.
/// Its private fields bind the checked controls to the exact branch definition.
#[derive(Clone, Debug)]
pub struct SectionBranchApproximation {
    definition: PatchSectionBranch,
    candidate: SectionCandidate,
    bound_mm: Rational,
    rho: Rational,
    speed: Rational,
    slope_lower: Rational,
    cap_mm: Rational,
    work: u64,
    proof: GeometricProof<&'static str>,
}
impl SectionBranchApproximation {
    /// Exact branch represented by the controls.
    pub fn definition(&self) -> &PatchSectionBranch {
        &self.definition
    }
    /// Checked controls, unchanged from the supplied candidate.
    pub fn candidate(&self) -> &SectionCandidate {
        &self.candidate
    }
    /// Whole-component Euclidean distance bound in millimetres, both ways.
    pub fn bound_mm(&self) -> &Rational {
        &self.bound_mm
    }
    /// Proved parameter-space tube half-width.
    pub fn half_width(&self) -> &Rational {
        &self.rho
    }
    /// Proved Euclidean speed bound in the transverse patch parameter.
    pub fn speed_bound(&self) -> &Rational {
        &self.speed
    }
    /// Positive lower bound for the cut polynomial's transverse derivative.
    /// Its reciprocal bounds the graph equation's transverse inverse.
    pub fn transverse_slope_lower(&self) -> &Rational {
        &self.slope_lower
    }
    /// Kernel cap recorded when this representation was checked.
    pub fn cap_mm(&self) -> &Rational {
        &self.cap_mm
    }
    /// Deterministic work spent checking the candidate.
    pub fn work(&self) -> u64 {
        self.work
    }
    /// Producer-owned correspondence theorem with retained definition and data.
    pub fn proof(&self) -> &GeometricProof<&'static str> {
        &self.proof
    }
}

fn add(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a + b)
}
fn mul(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a * b)
}
fn product(a: &[Rational], b: &[Rational], dec: &mut Decider) -> Result<Vec<Rational>, Refusal> {
    let mut c = vec![Rational::zero(); a.len() + b.len() - 1];
    for (i, a) in a.iter().enumerate() {
        if a.is_zero() {
            continue;
        }
        for (j, b) in b.iter().enumerate() {
            if b.is_zero() {
                continue;
            }
            c[i + j] = add(&c[i + j], &mul(a, b, dec)?, dec)?;
        }
    }
    while c.len() > 1 && c.last().unwrap().is_zero() {
        c.pop();
    }
    Ok(c)
}
fn cubic(controls: &[Rational; 4], dec: &mut Decider) -> Result<Vec<Rational>, Refusal> {
    let [a, b, c, d] = controls;
    Ok(vec![
        a.clone(),
        mul(&Rational::from_i64(3), &add(b, &(-a), dec)?, dec)?,
        mul(
            &Rational::from_i64(3),
            &add(
                &add(a, &(-mul(&Rational::from_i64(2), b, dec)?), dec)?,
                c,
                dec,
            )?,
            dec,
        )?,
        add(
            &add(&(-a), &mul(&Rational::from_i64(3), b, dec)?, dec)?,
            &add(&(-mul(&Rational::from_i64(3), c, dec)?), d, dec)?,
            dec,
        )?,
    ])
}
fn binomial(n: usize, k: usize) -> Rational {
    let mut value = Integer::one();
    for i in 0..k {
        value = (value * Integer::from_u64((n - i) as u64))
            .div_rem_floor(&Integer::from_u64((i + 1) as u64))
            .expect("nonzero binomial divisor")
            .0;
    }
    Rational::from_integer(value)
}
fn bernstein(power: &[Rational], dec: &mut Decider) -> Result<Vec<Rational>, Refusal> {
    let n = power.len() - 1;
    let mut output = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let mut value = Rational::zero();
        for (i, c) in power.iter().enumerate().take(k + 1) {
            let a = binomial(k, i);
            let b = binomial(n, i);
            dec.reserve_exact_arithmetic(&a, &b, true)?;
            let ratio = &a / &b;
            value = add(&value, &mul(c, &ratio, dec)?, dec)?;
        }
        output.push(value);
    }
    Ok(output)
}
fn composed(
    equation: &RealPoly,
    u: &[Rational],
    v: &[Rational],
    dec: &mut Decider,
) -> Result<Vec<Rational>, Refusal> {
    let degree = degrees(equation, dec)?;
    let mut rows = vec![vec![Rational::zero(); degree[1] as usize + 1]; degree[0] as usize + 1];
    for (monomial, coefficient) in equation.terms() {
        let Some(c) = coefficient.rational_ref() else {
            return Err(failure(
                "the first tube checker requires rational cut-equation coefficients",
            ));
        };
        let e = monomial.exponents();
        rows[e[0] as usize][e[1] as usize] = c.clone();
    }
    let mut result = vec![Rational::zero()];
    for row in rows.iter().rev() {
        let mut p = vec![Rational::zero()];
        for c in row.iter().rev() {
            p = product(&p, v, dec)?;
            p[0] = add(&p[0], c, dec)?;
        }
        result = product(&result, u, dec)?;
        result.resize(result.len().max(p.len()), Rational::zero());
        for (value, c) in result.iter_mut().zip(p) {
            *value = add(value, &c, dec)?;
        }
    }
    while result.len() > 1 && result.last().unwrap().is_zero() {
        result.pop();
    }
    bernstein(&result, dec)
}
fn single_sign(values: &[Rational]) -> Option<Sign> {
    let sign = values.first()?.sign();
    (sign != Sign::Zero && values.iter().all(|c| c.sign() == sign)).then_some(sign)
}
fn speed_bound(
    branch: &PatchSectionBranch,
    solved: usize,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let span = &branch.patch.data.spans[branch.cells[0].span];
    let [du, dv] = branch.patch.carrier().degrees();
    let count = [du + 1, dv + 1];
    let controls = branch
        .patch
        .carrier()
        .span_bezier_controls(span.indices, dec)?;
    let mut h = Vec::new();
    for control in controls {
        let mut row: [Rational; 4] = std::array::from_fn(|_| Rational::zero());
        for (to, from) in row.iter_mut().zip(control) {
            *to = dec.rational_value(&from)?.ok_or_else(|| {
                failure("the first tube speed bound requires rational homogeneous controls")
            })?;
        }
        h.push(row);
    }
    let wmin = h.iter().map(|c| c[3].clone()).min().unwrap();
    if wmin <= Rational::zero() {
        return Err(failure(
            "the tube needs a positive homogeneous weight bound",
        ));
    }
    let mut derivative: [Rational; 4] = std::array::from_fn(|_| Rational::zero());
    for i in 0..count[0] - usize::from(solved == 0) {
        for j in 0..count[1] - usize::from(solved == 1) {
            let k = i * count[1] + j;
            let next = k + if solved == 0 { count[1] } else { 1 };
            let factor = Rational::from_i64((count[solved] - 1) as i64)
                / (&span.domain[solved][1] - &span.domain[solved][0]);
            for a in 0..4 {
                derivative[a] = derivative[a]
                    .clone()
                    .max(mul(&factor, &add(&h[next][a], &(-&h[k][a]), dec)?, dec)?.abs());
            }
        }
    }
    let mut squared = Rational::zero();
    for a in 0..3 {
        let coordinate = h.iter().map(|c| (&c[a] / &c[3]).abs()).max().unwrap();
        let numerator = add(&derivative[a], &mul(&coordinate, &derivative[3], dec)?, dec)?;
        dec.reserve_exact_arithmetic(&numerator, &wmin, true)?;
        let bound = &numerator / &wmin;
        squared = add(&squared, &mul(&bound, &bound, dec)?, dec)?;
    }
    Ok(math::sqrt_upper(&squared))
}
impl PatchSectionBranch {
    /// Checks opposite tube-wall signs, complete monotone coverage and exact
    /// border end caps. No residual-to-two-surfaces heuristic is used.
    /// This first representation checker covers single-span graph branches.
    pub fn certify_candidate(
        &self,
        candidate: &SectionCandidate,
        dec: &mut Decider,
    ) -> Result<SectionBranchApproximation, Refusal> {
        let start = dec.stats().work;
        let first = &self.cells[0];
        if self.closed || self.cells.iter().any(|c| c.span != first.span) {
            return Err(Refusal::new(
                self.surface.codes().candidate_domain,
                "the first tube representation requires one open single-span graph branch",
            ));
        }
        // The free parameter must run between opposite span borders. Cells may
        // have been classified along either axis; the tube below proves its own
        // transverse uniqueness over every candidate box.
        let endpoints = self.parameter_endpoints();
        let span = &self.patch.data.spans[first.span];
        let mut graph = None;
        for solved in [first.solved, 1 - first.solved] {
            let free = 1 - solved;
            let orientation = exact_sign(&(&endpoints[1][free] - &endpoints[0][free]), dec)?;
            let interval = &span.domain[free];
            let bounds = if orientation == Sign::Positive {
                [interval[0].clone(), interval[1].clone()]
            } else {
                [interval[1].clone(), interval[0].clone()]
            };
            let mut covers = true;
            for i in 0..2 {
                if exact_sign(
                    &(&endpoints[i][free] - &Real::from_rational(bounds[i].clone())),
                    dec,
                )? != Sign::Zero
                {
                    covers = false;
                    break;
                }
            }
            if covers {
                graph = Some((solved, orientation));
                break;
            }
        }
        let Some((solved, orientation)) = graph else {
            return Err(failure(
                "the exact branch does not cover a complete free-parameter interval",
            ));
        };
        let free = 1 - solved;
        let interval = &span.domain[free];
        for cell in self.cells.iter() {
            if exact_sign(&(&cell.ends[1][free] - &cell.ends[0][free]), dec)? != orientation {
                return Err(failure(
                    "the retained branch is not a single monotone graph",
                ));
            }
        }
        let Some(begin) = candidate.pieces.first() else {
            return Err(failure(
                "the candidate is empty; the exact branch remains available",
            ));
        };
        let end = candidate.pieces.last().unwrap();
        let forward = begin[0][free] == interval[0] && end[3][free] == interval[1];
        let backward = begin[0][free] == interval[1] && end[3][free] == interval[0];
        if !forward && !backward {
            return Err(Refusal::new(
                self.surface.codes().candidate_endpoint,
                "candidate end caps must reach both exact patch borders",
            ));
        }
        let direction = if forward {
            Sign::Positive
        } else {
            Sign::Negative
        };
        for (i, piece) in candidate.pieces.iter().enumerate() {
            if i > 0 && candidate.pieces[i - 1][3] != piece[0] {
                return Err(failure("the candidate pieces do not join exactly"));
            }
            for pair in piece.windows(2) {
                if (&pair[1][free] - &pair[0][free]).sign() != direction {
                    return Err(failure(
                        "the candidate is not strictly monotone in the free parameter",
                    ));
                }
            }
            for p in piece {
                for (a, value) in p.iter().enumerate() {
                    if *value < span.domain[a][0] || *value > span.domain[a][1] {
                        return Err(failure("candidate controls leave the retained patch span"));
                    }
                }
            }
        }
        let speed = speed_bound(self, solved, dec)?;
        if speed <= Rational::zero() {
            return Err(failure("the transverse speed bound is not positive"));
        }
        let settings = LocalToleranceSettings::default();
        let equation = self.surface.equation(span);
        for cap in [settings.kernel_target_mm(), settings.kernel_cap_mm()] {
            let mut rho = Rational::one();
            while mul(&rho, &speed, dec)? > *cap {
                dec.check_deadline()?;
                rho = rho.mul_pow2(-1);
            }
            let mut accepted = true;
            let mut slope_lower: Option<Rational> = None;
            for (index, piece) in candidate.pieces.iter().enumerate() {
                // The end caps bind this tube to the selected exact component.
                for (at, control) in [(0, 0), (1, 3)] {
                    if (at == 0 && index == 0) || (at == 1 && index + 1 == candidate.pieces.len()) {
                        let endpoint = if forward == (orientation == Sign::Positive) {
                            at
                        } else {
                            1 - at
                        };
                        let center = Real::from_rational(piece[control][solved].clone());
                        let r = Real::from_rational(rho.clone());
                        if exact_sign(&(&endpoints[endpoint][solved] - &(&center - &r)), dec)?
                            == Sign::Negative
                            || exact_sign(&(&(&center + &r) - &endpoints[endpoint][solved]), dec)?
                                == Sign::Negative
                        {
                            accepted = false;
                            break;
                        }
                    }
                }
                if !accepted {
                    break;
                }
                let mut box_domain = span.domain.clone();
                for a in 0..2 {
                    let lo = piece.iter().map(|p| p[a].clone()).min().unwrap();
                    let hi = piece.iter().map(|p| p[a].clone()).max().unwrap();
                    box_domain[a] = if a == solved {
                        [&lo - &rho, &hi + &rho]
                    } else {
                        [lo, hi]
                    };
                }
                if (0..2).any(|a| {
                    box_domain[a][0] < span.domain[a][0] || box_domain[a][1] > span.domain[a][1]
                }) {
                    accepted = false;
                    break;
                }
                let controls = self.patch.carrier().span_bezier_controls_on(
                    span.indices,
                    box_domain.each_ref().map(|b| [&b[0], &b[1]]),
                    dec,
                )?;
                let (values, count) =
                    self.surface
                        .net(&controls, self.patch.carrier().degrees(), dec)?;
                let net = Net {
                    count,
                    values,
                    domain: box_domain,
                    depth: 0,
                };
                let derivative = net.derivative(solved, dec)?;
                if uniform(&derivative, dec)?.is_none() {
                    accepted = false;
                    break;
                }
                for value in derivative {
                    let lower = dec
                        .rational_value(&value)?
                        .ok_or_else(|| {
                            failure("the tube needs an exact rational transverse slope bound")
                        })?
                        .abs()
                        * Rational::from_i64((net.count[solved] - 1) as i64)
                        / (&net.domain[solved][1] - &net.domain[solved][0]);
                    slope_lower = Some(slope_lower.map_or(lower.clone(), |old| old.min(lower)));
                }
                let free_controls = std::array::from_fn(|i| piece[i][free].clone());
                let solved_controls = std::array::from_fn(|i| piece[i][solved].clone());
                let free_poly = cubic(&free_controls, dec)?;
                let mut lower = cubic(&solved_controls, dec)?;
                let mut upper = lower.clone();
                lower[0] = add(&lower[0], &(-&rho), dec)?;
                upper[0] = add(&upper[0], &rho, dec)?;
                if solved_controls
                    .iter()
                    .any(|v| v - &rho < span.domain[solved][0] || v + &rho > span.domain[solved][1])
                {
                    accepted = false;
                    break;
                }
                let walls = if solved == 1 {
                    [
                        composed(&equation, &free_poly, &lower, dec)?,
                        composed(&equation, &free_poly, &upper, dec)?,
                    ]
                } else {
                    [
                        composed(&equation, &lower, &free_poly, dec)?,
                        composed(&equation, &upper, &free_poly, dec)?,
                    ]
                };
                let signs = [single_sign(&walls[0]), single_sign(&walls[1])];
                if !matches!(
                    signs,
                    [Some(Sign::Positive), Some(Sign::Negative)]
                        | [Some(Sign::Negative), Some(Sign::Positive)]
                ) {
                    accepted = false;
                    break;
                }
            }
            if accepted {
                let bound_mm = mul(&rho, &speed, dec)?;
                return Ok(SectionBranchApproximation {
                    definition: self.clone(),
                    candidate: candidate.clone(),
                    bound_mm,
                    rho,
                    speed,
                    slope_lower: slope_lower.expect("nonempty checked tube"),
                    cap_mm: settings.kernel_cap_mm().clone(),
                    work: dec.stats().work - start,
                    proof: GeometricProof::Approximation {
                        definition: "retained exact patch-section branch and cell cover",
                        theorem: "opposite Bernstein tube walls, monotone graph coverage, exact border caps and rational surface speed bound",
                    },
                });
            }
        }
        Err(Refusal::new(
            self.surface.codes().candidate_cap,
            "no complete two-sided candidate tube was certified within the unchanged 1 micrometre kernel cap",
        ))
    }
}
