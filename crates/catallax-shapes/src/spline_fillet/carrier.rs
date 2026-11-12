//! The spine's branch, the fit's nodes, the exact carrier and its rails, and
//! the carrier's regularity.

use super::jets::{Meter, ball};
use super::*;
use crate::rolling_envelope::binomial;

/// The M10 meeting of the two offsets: exactly one open branch, spanning the
/// requested range of `tau`.
pub(super) fn spine(
    def: &SplineFilletDefinition,
    dec: &mut Decider,
) -> Result<PatchPairIntersection, Refusal> {
    let offset = |k: usize| {
        let side = Rational::from_i64(if def.sides[k] == Sign::Positive {
            1
        } else {
            -1
        });
        def.supports[k].offset(&def.radius * &side)
    };
    let meeting = offset(0).offset_meeting(&offset(1), dec)?;
    let [branch] = meeting.branches() else {
        return Err(spine_refusal(format!(
            "the two offsets meet in {} branches; a fillet follows exactly one",
            meeting.branches().len()
        )));
    };
    if branch.is_closed() {
        return Err(spine_refusal("the offsets meet in a closed loop"));
    }
    let mut ends = branch.ends().map(|end| end.enclosure()[def.axis].clone());
    if ends[0][0] > ends[1][0] {
        ends.swap(0, 1);
    }
    if ends[0][1] > def.range[0] || ends[1][0] < def.range[1] {
        return Err(spine_refusal(format!(
            "the branch spans tau from about {} to {}, which does not cover [{}, {}]",
            ends[0][0].to_f64_nearest(),
            ends[1][1].to_f64_nearest(),
            def.range[0],
            def.range[1]
        )));
    }
    Ok(meeting)
}

fn spine_refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new(SPLINE_FILLET_SPINE, detail)
}

/// A binary64 seed on the branch: the arc end whose `tau` is nearest the
/// start of the range, as `(tau, x)`.
pub(super) fn seed(
    def: &SplineFilletDefinition,
    spine: &PatchPairIntersection,
) -> Result<(f64, [f64; 3]), Refusal> {
    let start = def.range[0].to_f64_nearest();
    let mut best: Option<(f64, f64, [f64; 3])> = None;
    for arc in spine.branches().iter().flat_map(|b| b.arcs()) {
        for end in &arc.ends {
            let e = end.enclosure();
            let mid = |axis: usize| {
                ((&e[axis][0] + &e[axis][1]) / Rational::from_i64(2)).to_f64_nearest()
            };
            let tau = mid(def.axis);
            let x = [mid(1 - def.axis), mid(2), mid(3)];
            let distance = (tau - start).abs();
            if best.is_none_or(|b| distance < b.0) {
                best = Some((distance, tau, x));
            }
        }
    }
    let (_, tau, x) = best.ok_or_else(|| spine_refusal("the branch has no arcs"))?;
    Ok((tau, x))
}

/// Binary64 spine estimates at every node, marching from the seed. They only
/// propose boxes; every node and span is proved exactly.
pub(super) fn estimates(
    exact: &Exact,
    range: &[Rational; 2],
    spans: usize,
    (tau, x): (f64, [f64; 3]),
) -> Result<Vec<[f64; 3]>, Refusal> {
    let [a, c] = range.each_ref().map(Rational::to_f64_nearest);
    let failed = || spine_refusal("binary64 marching along the spine failed");
    // Walk from the seed to the start in small steps.
    let mut current = x;
    let steps = 64;
    for k in 1..=steps {
        let t = tau + (a - tau) * f64::from(k) / f64::from(steps);
        current = exact.newton(t, current).ok_or_else(failed)?;
    }
    let mut out = vec![current];
    for i in 1..=spans {
        let t = a + (c - a) * i as f64 / spans as f64;
        let guess = if out.len() >= 2 {
            let (p, q) = (out[out.len() - 2], out[out.len() - 1]);
            std::array::from_fn(|k| 2.0 * q[k] - p[k])
        } else {
            out[out.len() - 1]
        };
        out.push(exact.newton(t, guess).ok_or_else(failed)?);
    }
    Ok(out)
}

pub(super) fn node_tau(range: &[Rational; 2], spans: usize, i: usize) -> Rational {
    let step = (&range[1] - &range[0]) / Rational::from_i64(spans as i64);
    &range[0] + &(step * Rational::from_i64(i as i64))
}

/// Rounds a ball's centre to a binary64 rational, keeping node data small.
fn short(value: &Ball) -> Rational {
    Rational::from_f64_exact(value.mid().to_f64_nearest()).unwrap_or_else(|| value.mid())
}

/// A rail border's node of `tau`, and the fitted row it pins with its value.
/// A border of `tau` itself is exact and pins nothing. Otherwise a binary64
/// march along the spine finds the one place in the range where the row
/// crosses the value; the node there is rounded to binary64, and the
/// two-sided bound measures the pinned value against the exact one.
pub(super) fn border(
    exact: &Exact,
    def: &SplineFilletDefinition,
    border: &RailBorder,
    seed: (f64, [f64; 3]),
) -> Result<(Rational, Option<(usize, Rational)>), Refusal> {
    if border.support == 0 && border.parameter == def.axis {
        return Ok((border.value.clone(), None));
    }
    let row = if border.support == 0 {
        0
    } else {
        1 + border.parameter
    };
    let steps = 256;
    let samples = estimates(exact, &def.range, steps, seed)?;
    let [a, c] = def.range.each_ref().map(Rational::to_f64_nearest);
    let at = |i: usize| a + (c - a) * i as f64 / steps as f64;
    let value = border.value.to_f64_nearest();
    let crossings: Vec<usize> = (0..steps)
        .filter(|&i| (samples[i][row] - value) * (samples[i + 1][row] - value) <= 0.0)
        .collect();
    let [i] = crossings[..] else {
        return Err(definition(format!(
            "support {}'s parameter {} meets {} at {} places in the range, not one",
            border.support,
            border.parameter,
            border.value,
            crossings.len()
        )));
    };
    let failed = || spine_refusal("binary64 Newton on the spine failed");
    let (mut lo, mut hi) = (at(i), at(i + 1));
    let mut guess = samples[i];
    let below = samples[i][row] < value;
    for _ in 0..80 {
        let middle = 0.5 * (lo + hi);
        guess = exact.newton(middle, guess).ok_or_else(failed)?;
        if (guess[row] < value) == below {
            lo = middle;
        } else {
            hi = middle;
        }
    }
    let tau = Rational::from_f64_exact(0.5 * (lo + hi)).ok_or_else(failed)?;
    if tau <= def.range[0] || tau >= def.range[1] {
        return Err(definition("a rail border must lie inside the range of tau"));
    }
    Ok((tau, Some((row, border.value.clone()))))
}

/// A binary64 spine estimate at `tau`, the given fraction of the way from
/// one node's estimate to the next.
pub(super) fn between(
    exact: &Exact,
    tau: &Rational,
    ends: [[f64; 3]; 2],
    fraction: f64,
) -> Result<[f64; 3], Refusal> {
    let guess = std::array::from_fn(|k| ends[0][k] + fraction * (ends[1][k] - ends[0][k]));
    exact
        .newton(tau.to_f64_nearest(), guess)
        .ok_or_else(|| spine_refusal("binary64 Newton on the spine failed"))
}

/// One node's proved point jet: its fit data, rounded, and the exact value
/// and slope enclosures the bound measures that data against.
pub(super) fn node(
    exact: &Exact,
    tau: &Rational,
    estimate: [f64; 3],
    dec: &mut Decider,
) -> Result<(Node, NodeData), Refusal> {
    let solution = exact.solve([tau, tau], [estimate, estimate], 1, dec)?;
    let mut meter = Meter::default();
    let rows = exact.rows(&solution, &mut meter)?;
    meter.charge(dec)?;
    let jets = fitted(&solution, &rows);
    let node = Node {
        tau: tau.clone(),
        values: std::array::from_fn(|j| short(&jets[j].0[0])),
        slopes: std::array::from_fn(|j| short(&jets[j].0[1])),
    };
    Ok((node, jets.map(|jet| jet.data())))
}

/// The fitted rows' exact jets, in node order.
pub(super) fn fitted(solution: &jets::Solution, rows: &jets::Rows) -> [Jet; ROWS] {
    [
        solution.x[0].clone(),
        solution.x[1].clone(),
        solution.x[2].clone(),
        rows.middle[0].clone(),
        rows.middle[1].clone(),
        rows.middle[2].clone(),
        rows.weight.clone(),
    ]
}

// Exact Bernstein arithmetic on [0, 1], each operation charged.
type Bern = Vec<Rational>;

fn mul(a: &Rational, c: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, c, true)?;
    Ok(a * c)
}
fn add(a: &Rational, c: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, c, false)?;
    Ok(a + c)
}

fn product(a: &Bern, c: &Bern, dec: &mut Decider) -> Result<Bern, Refusal> {
    let (m, n) = (a.len() - 1, c.len() - 1);
    let mut out = vec![Rational::zero(); m + n + 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in c.iter().enumerate() {
            if y.is_zero() {
                continue;
            }
            let factor = binomial(m, i) * binomial(n, j) / binomial(m + n, i + j);
            let term = mul(&mul(x, y, dec)?, &factor, dec)?;
            out[i + j] = add(&out[i + j], &term, dec)?;
        }
    }
    Ok(out)
}

fn elevate(a: &Bern, degree: usize, dec: &mut Decider) -> Result<Bern, Refusal> {
    let mut a = a.clone();
    while a.len() - 1 < degree {
        let n = a.len();
        let mut next = Vec::with_capacity(n + 1);
        next.push(a[0].clone());
        for i in 1..n {
            let left = Rational::from_ratio(i as i64, n as i64).expect("nonzero");
            let right = Rational::one() - left.clone();
            let value = add(&mul(&a[i - 1], &left, dec)?, &mul(&a[i], &right, dec)?, dec)?;
            next.push(value);
        }
        next.push(a[n - 1].clone());
        a = next;
    }
    Ok(a)
}

fn sum(a: &Bern, c: &Bern, dec: &mut Decider) -> Result<Bern, Refusal> {
    let degree = a.len().max(c.len()) - 1;
    let (a, c) = (elevate(a, degree, dec)?, elevate(c, degree, dec)?);
    a.iter().zip(&c).map(|(x, y)| add(x, y, dec)).collect()
}

/// A two-parameter polynomial at Bernstein polynomials in `tau`.
fn compose(poly: &MPoly, params: [&Bern; 2], dec: &mut Decider) -> Result<Bern, Refusal> {
    let mut powers: [Vec<Bern>; 2] = [vec![vec![Rational::one()]], vec![vec![Rational::one()]]];
    let mut total: Bern = vec![Rational::zero()];
    for (monomial, coefficient) in poly.terms() {
        let e = monomial.exponents();
        for axis in 0..2 {
            while powers[axis].len() <= e[axis] as usize {
                let next = product(powers[axis].last().expect("one"), params[axis], dec)?;
                powers[axis].push(next);
            }
        }
        let term = product(&powers[0][e[0] as usize], &powers[1][e[1] as usize], dec)?;
        let scaled: Bern = term
            .iter()
            .map(|x| mul(x, &coefficient, dec))
            .collect::<Result<_, _>>()?;
        total = sum(&total, &scaled, dec)?;
    }
    Ok(total)
}

/// The Hermite cubic of one row on `[t0, t1]`, in Bernstein form.
pub(super) fn cubic(span: [&Rational; 2], values: [&Rational; 2], slopes: [&Rational; 2]) -> Bern {
    let third = (span[1] - span[0]) / Rational::from_i64(3);
    vec![
        values[0].clone(),
        values[0] + &(&third * slopes[0]),
        values[1] - &(&third * slopes[1]),
        values[1].clone(),
    ]
}

/// Most halvings of one rail cubic while its domain containment is proved.
const RAIL_DEPTH: u32 = 12;

/// Why a rail cubic is not inside its domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Escape {
    /// A point of the cubic lies outside.
    Leaves,
    /// Neither inside nor outside was proved within the halvings.
    Unproved,
}

/// Proves one rail parameter's cubic inside `domain`. A cubic lies in the
/// hull of its Bernstein coefficients, so a piece whose coefficients all lie
/// inside is inside. An end coefficient is a point of the cubic, so one
/// outside proves the rail leaves. Other pieces are halved exactly.
pub(super) fn cubic_inside(
    coefficients: &[Rational; 4],
    domain: &[Rational; 2],
    dec: &mut Decider,
) -> Result<Option<Escape>, Refusal> {
    let [lo, hi] = domain;
    let inside = |value: &Rational| value >= lo && value <= hi;
    let mut pending = vec![(coefficients.to_vec(), 0u32)];
    let half = Rational::from_ratio(1, 2).expect("nonzero");
    while let Some((c, depth)) = pending.pop() {
        dec.check_deadline()?;
        if c.iter().all(inside) {
            continue;
        }
        if !inside(&c[0]) || !inside(&c[3]) {
            return Ok(Some(Escape::Leaves));
        }
        if depth >= RAIL_DEPTH {
            return Ok(Some(Escape::Unproved));
        }
        // De Casteljau at the middle.
        let mut rows = vec![c];
        for _ in 0..3 {
            let last = rows.last().expect("one row");
            let mut next = Vec::with_capacity(last.len() - 1);
            for pair in last.windows(2) {
                dec.reserve_exact_arithmetic(&pair[0], &pair[1], false)?;
                next.push((&pair[0] + &pair[1]) * half.clone());
            }
            rows.push(next);
        }
        let left: Vec<Rational> = rows.iter().map(|row| row[0].clone()).collect();
        let right: Vec<Rational> = rows
            .iter()
            .rev()
            .map(|row| row[row.len() - 1].clone())
            .collect();
        pending.push((right, depth + 1));
        pending.push((left, depth + 1));
    }
    Ok(None)
}

/// Proves every returned rail inside its support's domain, or refuses.
pub(super) fn rails_inside(
    rails: &[Vec<RailSpan>; 2],
    supports: &[RationalPatch; 2],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for (side, (rail, support)) in rails.iter().zip(supports).enumerate() {
        for span in rail {
            for (parameter, coefficients) in span.parameters.iter().enumerate() {
                let domain = &support.domain()[parameter];
                let Some(escape) = cubic_inside(coefficients, domain, dec)? else {
                    continue;
                };
                let how = match escape {
                    Escape::Leaves => "leaves",
                    Escape::Unproved => "was not proved inside",
                };
                return Err(Refusal::new(
                    SPLINE_FILLET_RAIL_DOMAIN,
                    format!(
                        "the fitted rail on support {side} {how} its parameter {parameter}'s \
                         domain [{}, {}] on tau in [{}, {}]",
                        domain[0], domain[1], span.tau[0], span.tau[1]
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// Homogeneous Bernstein rows of one span: first rail, middle, second rail.
pub(super) type SpanRows = [[Bern; 4]; 3];

/// The carrier, both rails and each span's homogeneous rows.
#[allow(clippy::type_complexity)]
pub(super) fn build(
    exact: &Exact,
    nodes: &[Node],
    dec: &mut Decider,
) -> Result<(BSplineSurface, [Vec<RailSpan>; 2], Vec<SpanRows>), Refusal> {
    let mut rows: Vec<SpanRows> = Vec::with_capacity(nodes.len() - 1);
    let mut rails: [Vec<RailSpan>; 2] = [Vec::new(), Vec::new()];
    for pair in nodes.windows(2) {
        dec.check_deadline()?;
        let span = [&pair[0].tau, &pair[1].tau];
        let row = |j: usize| {
            cubic(
                span,
                [&pair[0].values[j], &pair[1].values[j]],
                [&pair[0].slopes[j], &pair[1].slopes[j]],
            )
        };
        let tau = elevate(&vec![span[0].clone(), span[1].clone()], 3, dec)?;
        let first = exact.params(&[row(0), row(1), row(2)], &tau);
        let rail = |params: &[Bern; 2]| RailSpan {
            tau: [span[0].clone(), span[1].clone()],
            parameters: params
                .each_ref()
                .map(|p| p.clone().try_into().expect("cubic rail parameter")),
        };
        rails[0].push(rail(&first[0]));
        rails[1].push(rail(&first[1]));
        let mut span_rows: Vec<[Bern; 4]> = Vec::with_capacity(3);
        for (k, params) in first.iter().enumerate() {
            let support = &exact.supports[k];
            let mut coordinates = Vec::with_capacity(4);
            for poly in &support.exact {
                coordinates.push(compose(poly, [&params[0], &params[1]], dec)?);
            }
            span_rows.push(coordinates.try_into().expect("four coordinates"));
            if k == 0 {
                // The middle row: w1 times the middle point, and w1.
                let weight = row(WEIGHT);
                let mut middle = Vec::with_capacity(4);
                for i in 0..3 {
                    middle.push(product(&weight, &row(MIDDLE + i), dec)?);
                }
                middle.push(weight);
                span_rows.push(middle.try_into().expect("four coordinates"));
            }
        }
        rows.push(span_rows.try_into().expect("three rows"));
    }
    let degree = rows
        .iter()
        .flatten()
        .flatten()
        .map(|p: &Bern| p.len() - 1)
        .max()
        .expect("nonempty");
    for span in &mut rows {
        for row in span.iter_mut() {
            for coordinate in row.iter_mut() {
                *coordinate = elevate(coordinate, degree, dec)?;
            }
        }
    }
    let spans = rows.len();
    let count = spans * degree + 1;
    let mut controls = Vec::with_capacity(count * 3);
    let mut weights = Vec::with_capacity(count * 3);
    for index in 0..count {
        let (span, local) = if index == count - 1 {
            (spans - 1, degree)
        } else {
            (index / degree, index % degree)
        };
        for row in &rows[span] {
            let w = &row[3][local];
            if w.sign() != Sign::Positive {
                return Err(Refusal::new(
                    SPLINE_FILLET_IRREGULAR,
                    "a fitted carrier weight is not positive",
                ));
            }
            let mut coordinate = |i: usize| -> Result<Real, Refusal> {
                dec.reserve_exact_arithmetic(&row[i][local], w, true)?;
                Ok(Real::from_rational(&row[i][local] / w))
            };
            controls.push(Point3::new(coordinate(0)?, coordinate(1)?, coordinate(2)?));
            weights.push(Real::from_rational(w.clone()));
        }
    }
    let mut knots = vec![nodes[0].tau.clone(); degree + 1];
    for node in &nodes[1..nodes.len() - 1] {
        knots.extend(std::iter::repeat_n(node.tau.clone(), degree));
    }
    knots.extend(std::iter::repeat_n(
        nodes[nodes.len() - 1].tau.clone(),
        degree + 1,
    ));
    let side = vec![
        Rational::zero(),
        Rational::zero(),
        Rational::zero(),
        Rational::one(),
        Rational::one(),
        Rational::one(),
    ];
    let carrier = BSplineSurface::new(
        [degree, 2],
        [knots, side],
        [count, 3],
        controls,
        weights,
        dec,
    )?;
    Ok((carrier, rails, rows))
}

/// Power-basis coefficients in the span's local parameter `s` in `[0, 1]`,
/// as balls.
fn power(bern: &Bern, meter: &mut Meter) -> Vec<Ball> {
    let n = bern.len() - 1;
    let mut out = vec![ball(&Rational::zero()); n + 1];
    for (k, b) in bern.iter().enumerate() {
        if b.is_zero() {
            continue;
        }
        let b = ball(b);
        for (j, slot) in out.iter_mut().enumerate().skip(k) {
            let mut factor = binomial(n, k) * binomial(n - k, j - k);
            if (j - k) % 2 == 1 {
                factor = -factor;
            }
            *slot = &*slot + &(&b * &ball(&factor));
        }
    }
    let pairs = ((n + 1) * (n + 2) / 2) as u64;
    meter.balls(pairs, true);
    meter.balls(pairs, false);
    out
}

/// Coefficients of `p(c + sigma)` from those of `p(s)`, by repeated
/// synthetic division.
fn shifted(coefficients: &[Ball], c: &Ball, meter: &mut Meter) -> Vec<Ball> {
    let mut a = coefficients.to_vec();
    let n = a.len() - 1;
    for i in 0..n {
        for j in (i..n).rev() {
            a[j] = &a[j] + &(c * &a[j + 1]);
        }
    }
    let pairs = (n * (n + 1) / 2) as u64;
    meter.balls(pairs, true);
    meter.balls(pairs, false);
    a
}

/// One span's homogeneous rows in the power basis of `s`: row, coordinate.
type PowerRows = [[Vec<Ball>; 4]; 3];

/// The rows expanded about the middle of an interval of `s`, with the
/// interval's half width.
struct Slab {
    rows: PowerRows,
    half: Rational,
}

fn slab(rows: &PowerRows, s: &[Rational; 2], meter: &mut Meter) -> Slab {
    let two = Rational::from_i64(2);
    let centre = ball(&((&s[0] + &s[1]) / two.clone()));
    Slab {
        rows: rows.each_ref().map(|row| {
            row.each_ref()
                .map(|coordinate| shifted(coordinate, &centre, meter))
        }),
        half: (&s[1] - &s[0]) / two,
    }
}

/// The carrier's normal over a box, up to a positive factor: with
/// homogeneous coordinates `H` and weight `W`, the cross product
/// `(W H_s - W_s H) x (W H_t - W_t H)`. The rows are combined with the
/// section basis coefficient by coefficient, about the box's centre in both
/// directions, before anything is bounded over the box; the rows move
/// together along `s`, and their common motion cancels exactly there.
fn normal(slab: &Slab, theta: &[Rational; 2], meter: &mut Meter) -> [Ball; 3] {
    let two = Rational::from_i64(2);
    let t = (&theta[0] + &theta[1]) / two.clone();
    let half = (&theta[1] - &theta[0]) / two.clone();
    let one = Rational::one();
    let m = &one - &t;
    // (1-t)^2, 2t(1-t), t^2 at t + delta, by powers of delta.
    let basis = [
        [&m * &m, -(&two * &m), one.clone()],
        [
            &(&two * &t) * &m,
            &two * &(&one - &(&two * &t)),
            -two.clone(),
        ],
        [&t * &t, &two * &t, one],
    ]
    .map(|row| row.map(|value| ball(&value)));
    let delta = ball(&-half.clone()).hull(&ball(&half));
    let delta2 = delta.pow_u32(2);
    let sigma = ball(&-slab.half.clone()).hull(&ball(&slab.half));
    let degree = slab.rows[0][0].len();
    let powers: Vec<Ball> = (0..degree).map(|j| sigma.pow_u32(j as u32)).collect();
    let zero = ball(&Rational::zero());
    let twice = ball(&two);
    // Per coordinate: H, H_s, H_t over the box.
    let parts: [[Ball; 3]; 4] = std::array::from_fn(|i| {
        let mut out = [zero.clone(), zero.clone(), zero.clone()];
        for j in 0..degree {
            let c: [Ball; 3] = std::array::from_fn(|m| {
                (0..3).fold(zero.clone(), |sum, k| {
                    &sum + &(&basis[k][m] * &slab.rows[k][i][j])
                })
            });
            let value = &(&c[0] + &(&c[1] * &delta)) + &(&c[2] * &delta2);
            let slope = &c[1] + &(&(&c[2] * &twice) * &delta);
            out[0] = &out[0] + &(&value * &powers[j]);
            out[2] = &out[2] + &(&slope * &powers[j]);
            if j > 0 {
                let factor = ball(&Rational::from_i64(j as i64));
                out[1] = &out[1] + &(&(&value * &factor) * &powers[j - 1]);
            }
        }
        out
    });
    meter.balls(4 * degree as u64 * (9 + 9), true);
    meter.balls(4 * degree as u64 * (9 + 5), false);
    let w = &parts[3];
    let along: [[Ball; 3]; 2] = [
        std::array::from_fn(|i| &(&w[0] * &parts[i][1]) - &(&w[1] * &parts[i][0])),
        std::array::from_fn(|i| &(&w[0] * &parts[i][2]) - &(&w[2] * &parts[i][0])),
    ];
    meter.balls(12 + 6, true);
    meter.balls(6 + 3, false);
    std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &(&along[0][j] * &along[1][k]) - &(&along[0][k] * &along[1][j])
    })
}

/// Most halvings of one box of the regularity proof.
const REGULAR_DEPTH: u32 = 8;

/// Proves the carrier's normal nonzero: on each box of each span its
/// component along the middle of its own enclosure is positive.
pub(super) fn regular(rows: &[SpanRows], dec: &mut Decider) -> Result<(), Refusal> {
    order::in_order(
        rows.len(),
        |index, dec| {
            let mut meter = Meter::default();
            let mut power_rows: PowerRows =
                rows[index].each_ref().map(|row| row.each_ref().map(|coordinate| power(coordinate, &mut meter)));
            // The normal is unchanged by translating every point by the same
            // `p`, which takes `H` to `H - W p`; a point of the span keeps the
            // translated coordinates small, and with them the cancellation in
            // `W H_s - W_s H`.
            let middle = |coordinate: &Vec<Ball>| {
                coordinate.iter().rev().fold(0.0, |t, c| t * 0.5 + c.mid().to_f64_nearest())
            };
            let w = middle(&power_rows[1][3]);
            let p: [Ball; 3] = std::array::from_fn(|i| {
                ball(&Rational::from_f64_exact(middle(&power_rows[1][i]) / w).unwrap_or_else(Rational::zero))
            });
            for row in &mut power_rows {
                let weight = row[3].clone();
                for (i, offset) in p.iter().enumerate() {
                    for (value, w) in row[i].iter_mut().zip(&weight) {
                        *value = &*value - &(w * offset);
                    }
                }
                meter.balls(3 * weight.len() as u64, true);
                meter.balls(3 * weight.len() as u64, false);
            }
            let quarter = |k: i64| Rational::from_ratio(k, 4).expect("nonzero");
            let mut pending: Vec<([Rational; 2], Vec<[Rational; 2]>, u32)> =
                vec![([Rational::zero(), Rational::one()], (0..4).map(|k| [quarter(k), quarter(k + 1)]).collect(), 0)];
            let two = Rational::from_i64(2);
            let halves = |r: &[Rational; 2]| {
                let m = (&r[0] + &r[1]) / two.clone();
                [[r[0].clone(), m.clone()], [m, r[1].clone()]]
            };
            while let Some((s, thetas, depth)) = pending.pop() {
                dec.check_deadline()?;
                let slab = slab(&power_rows, &s, &mut meter);
                let mut failed = Vec::new();
                for theta in thetas {
                    let over = normal(&slab, &theta, &mut meter);
                    let direction = over.each_ref().map(|value| {
                        ball(&Rational::from_f64_exact(value.mid().to_f64_nearest()).unwrap_or_else(Rational::zero))
                    });
                    let along = (0..3).fold(ball(&Rational::zero()), |t, i| &t + &(&over[i] * &direction[i]));
                    meter.balls(3, true);
                    meter.balls(3, false);
                    if along.lo().sign() != Sign::Positive {
                        failed.push(theta);
                    }
                }
                meter.charge(dec)?;
                if failed.is_empty() {
                    continue;
                }
                if depth >= REGULAR_DEPTH {
                    let theta = &failed[0];
                    return Err(Refusal::new(
                        SPLINE_FILLET_IRREGULAR,
                        format!(
                            "the carrier's normal was not proved nonzero on span {index} over s [{}, {}] theta [{}, {}]",
                            s[0], s[1], theta[0], theta[1]
                        ),
                    ));
                }
                let split: Vec<[Rational; 2]> = failed.iter().flat_map(halves).collect();
                for half in halves(&s) {
                    pending.push((half, split.clone(), depth + 1));
                }
            }
            Ok(())
        },
        dec,
    )
    .map(|_| ())
}
