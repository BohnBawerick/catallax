//! Exact normal bounds and graph controls. Neither full paper algorithm is implemented.
use super::{Q, choose, q, r};
use catallax_numbers::{Decider, Refusal, with_operation_budget};
use std::cmp::Ordering;

#[path = "paper2026.rs"]
mod paper2026;

const ALLOWANCE: u64 = 20_000_000;
type Result<T> = std::result::Result<T, Refusal>;
fn refuse(code: &'static str) -> Refusal {
    Refusal::new(code, "the bounded R13 predicate did not certify this case")
}
type Net = Vec<Vec<[Q; 3]>>;
type Pair = [Net; 2];
type Range = [Q; 2];

// Match numbers::Work::exact: numerator plus denominator limbs, with one
// cancellation checkpoint per reservation. Parsing immutable input is excluded.
#[derive(Default)]
struct Meter(Decider);
impl Meter {
    fn charge(&mut self, a: &Q, b: &Q, product: bool) -> Result<()> {
        let limbs = |x: &Q| {
            x.numerator().bit_length().div_ceil(64).max(1)
                + x.denominator().bit_length().div_ceil(64).max(1)
        };
        let units = if product {
            limbs(a).saturating_mul(limbs(b))
        } else {
            limbs(a).saturating_add(limbs(b))
        };
        self.0.charge_work(units.saturating_add(1))
    }
    fn add(&mut self, a: &Q, b: &Q) -> Result<Q> {
        self.charge(a, b, false)?;
        Ok(a + b)
    }
    fn sub(&mut self, a: &Q, b: &Q) -> Result<Q> {
        self.charge(a, b, false)?;
        Ok(a - b)
    }
    fn mul(&mut self, a: &Q, b: &Q) -> Result<Q> {
        self.charge(a, b, true)?;
        Ok(a * b)
    }
    fn div(&mut self, a: &Q, b: &Q) -> Result<Q> {
        self.charge(a, b, true)?;
        a.checked_div(b)
            .ok_or_else(|| refuse("r13.zero-denominator"))
    }
    fn cmp(&mut self, a: &Q, b: &Q) -> Result<Ordering> {
        // Comparing fractions cross-multiplies; reserve two products.
        self.charge(a, b, true)?;
        self.charge(a, b, true)?;
        Ok(a.cmp(b))
    }
    fn equal(&mut self, a: &Q, b: &Q) -> Result<()> {
        if self.cmp(a, b)?.is_eq() {
            Ok(())
        } else {
            Err(refuse("r13.identity-not-proved"))
        }
    }
}

// Dense tensor power coefficients, in a unit knot span, lowest powers first.
#[derive(Clone)]
struct P(Vec<Vec<Q>>);
impl P {
    fn zero(u: usize, v: usize) -> Self {
        Self(vec![vec![q(0); v]; u])
    }
    fn line(a: Q, b: Q) -> Self {
        Self(vec![vec![a], vec![b]])
    }
    fn add(&self, b: &Self, m: &mut Meter) -> Result<Self> {
        let mut out = Self::zero(
            self.0.len().max(b.0.len()),
            self.0[0].len().max(b.0[0].len()),
        );
        for p in [self, b] {
            for (i, row) in p.0.iter().enumerate() {
                for (j, x) in row.iter().enumerate() {
                    out.0[i][j] = m.add(&out.0[i][j], x)?;
                }
            }
        }
        Ok(out)
    }
    fn scale(&self, a: &Q, m: &mut Meter) -> Result<Self> {
        let mut out = self.clone();
        for x in out.0.iter_mut().flatten() {
            *x = m.mul(x, a)?;
        }
        Ok(out)
    }
    fn sub(&self, b: &Self, m: &mut Meter) -> Result<Self> {
        self.add(&b.scale(&q(-1), m)?, m)
    }
    fn mul(&self, b: &Self, m: &mut Meter) -> Result<Self> {
        let mut out = Self::zero(
            self.0.len() + b.0.len() - 1,
            self.0[0].len() + b.0[0].len() - 1,
        );
        for (i, row) in self.0.iter().enumerate() {
            for (j, x) in row.iter().enumerate() {
                for (k, br) in b.0.iter().enumerate() {
                    for (l, y) in br.iter().enumerate() {
                        let xy = m.mul(x, y)?;
                        out.0[i + k][j + l] = m.add(&out.0[i + k][j + l], &xy)?;
                    }
                }
            }
        }
        Ok(out)
    }
    fn derivative(&self, axis: usize, m: &mut Meter) -> Result<Self> {
        let mut out = Self::zero(
            self.0.len().saturating_sub(usize::from(axis == 0)).max(1),
            self.0[0]
                .len()
                .saturating_sub(usize::from(axis == 1))
                .max(1),
        );
        for (i, row) in self.0.iter().enumerate() {
            for (j, x) in row.iter().enumerate() {
                let degree = if axis == 0 { i } else { j };
                if degree > 0 {
                    out.0[i - usize::from(axis == 0)][j - usize::from(axis == 1)] =
                        m.mul(x, &q(degree as i64))?;
                }
            }
        }
        Ok(out)
    }
    fn eval(&self, u: &Q, v: &Q, m: &mut Meter) -> Result<Q> {
        let mut value = q(0);
        for row in self.0.iter().rev() {
            let mut inner = q(0);
            for c in row.iter().rev() {
                inner = m.mul(&inner, v)?;
                inner = m.add(&inner, c)?;
            }
            value = m.mul(&value, u)?;
            value = m.add(&value, &inner)?;
        }
        Ok(value)
    }
    fn bernstein(&self, m: &mut Meter) -> Result<Self> {
        let (n, d) = (self.0.len() - 1, self.0[0].len() - 1);
        let mut out = Self::zero(n + 1, d + 1);
        for i in 0..=n {
            for j in 0..=d {
                for k in 0..=i {
                    for l in 0..=j {
                        let u = m.div(&q(choose(i, k)), &q(choose(n, k)))?;
                        let v = m.div(&q(choose(j, l)), &q(choose(d, l)))?;
                        let uv = m.mul(&u, &v)?;
                        let term = m.mul(&self.0[k][l], &uv)?;
                        out.0[i][j] = m.add(&out.0[i][j], &term)?;
                    }
                }
            }
        }
        Ok(out)
    }
    fn bounds(&self, m: &mut Meter) -> Result<Range> {
        let c = self.bernstein(m)?;
        let mut out = [c.0[0][0].clone(), c.0[0][0].clone()];
        for x in c.0.iter().flatten() {
            if m.cmp(x, &out[0])?.is_lt() {
                out[0] = x.clone();
            }
            if m.cmp(x, &out[1])?.is_gt() {
                out[1] = x.clone();
            }
        }
        Ok(out)
    }
    fn positive_inside(&self, m: &mut Meter) -> Result<()> {
        let c = self.bernstein(m)?;
        let mut nonzero = false;
        for x in c.0.iter().flatten() {
            match m.cmp(x, &q(0))? {
                Ordering::Less => return Err(refuse("r13.sign-not-proved")),
                Ordering::Greater => nonzero = true,
                Ordering::Equal => (),
            }
        }
        if nonzero {
            Ok(())
        } else {
            Err(refuse("r13.zero-polynomial"))
        }
    }
    // Only used on univariate polynomials, including oriented boundary edges.
    fn affine_u(&self, a: &Q, b: &Q, m: &mut Meter) -> Result<Self> {
        assert!(self.0.iter().all(|row| row.len() == 1));
        let mut out = Self::zero(1, 1);
        let mut power = Self(vec![vec![q(1)]]);
        let delta = m.sub(b, a)?;
        let linear = Self::line(a.clone(), delta);
        for row in &self.0 {
            out = out.add(&power.scale(&row[0], m)?, m)?;
            power = power.mul(&linear, m)?;
        }
        Ok(out)
    }
    fn edge(&self, varying: usize, fixed: &Q, m: &mut Meter) -> Result<Self> {
        let len = if varying == 0 {
            self.0.len()
        } else {
            self.0[0].len()
        };
        let mut out = Self::zero(len, 1);
        for (i, row) in self.0.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                let (target, degree) = if varying == 0 { (i, j) } else { (j, i) };
                let mut term = c.clone();
                for _ in 0..degree {
                    term = m.mul(&term, fixed)?;
                }
                out.0[target][0] = m.add(&out.0[target][0], &term)?;
            }
        }
        Ok(out)
    }
}

fn basis(i: usize, degree: usize, span: usize, m: &mut Meter) -> Result<P> {
    const K: [i64; 10] = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3];
    if degree == 0 {
        return Ok(P(vec![vec![q(i64::from(i == span))]]));
    }
    let mut out = P::zero(1, 1);
    for (child, a, b, denominator) in [
        (i, K[span] - K[i], 1, K[i + degree] - K[i]),
        (
            i + 1,
            K[i + degree + 1] - K[span],
            -1,
            K[i + degree + 1] - K[i + 1],
        ),
    ] {
        if denominator != 0 {
            let factor = P::line(
                m.div(&q(a), &q(denominator))?,
                m.div(&q(b), &q(denominator))?,
            );
            out = out.add(&factor.mul(&basis(child, degree - 1, span, m)?, m)?, m)?;
        }
    }
    Ok(out)
}
fn bases(m: &mut Meter) -> Result<Vec<Vec<P>>> {
    (3..6)
        .map(|s| (0..6).map(|i| basis(i, 3, s, m)).collect())
        .collect()
}
fn patches(net: &Net, b: &[Vec<P>], m: &mut Meter) -> Result<Vec<[P; 3]>> {
    let mut out = Vec::new();
    for u in b {
        for v in b {
            let mut xyz = std::array::from_fn(|_| P::zero(4, 4));
            for (i, row) in net.iter().enumerate() {
                for (j, point) in row.iter().enumerate() {
                    let mut product = P::zero(4, 4);
                    for k in 0..4 {
                        for l in 0..4 {
                            product.0[k][l] = m.mul(&u[i].0[k][0], &v[j].0[l][0])?;
                        }
                    }
                    for axis in 0..3 {
                        xyz[axis] = xyz[axis].add(&product.scale(&point[axis], m)?, m)?;
                    }
                }
            }
            out.push(xyz);
        }
    }
    Ok(out)
}
fn cross(a: &[P; 3], b: &[P; 3], m: &mut Meter) -> Result<[P; 3]> {
    let mut out = std::array::from_fn(|_| P::zero(1, 1));
    for (axis, p) in out.iter_mut().enumerate() {
        let (j, k) = ((axis + 1) % 3, (axis + 2) % 3);
        *p = a[j].mul(&b[k], m)?.sub(&a[k].mul(&b[j], m)?, m)?;
    }
    Ok(out)
}
fn normal_bounds(net: &Net, b: &[Vec<P>], m: &mut Meter) -> Result<Vec<[Range; 3]>> {
    let mut out = Vec::new();
    for p in patches(net, b, m)? {
        let mut du = std::array::from_fn(|_| P::zero(1, 1));
        let mut dv = du.clone();
        for axis in 0..3 {
            du[axis] = p[axis].derivative(0, m)?;
            dv[axis] = p[axis].derivative(1, m)?;
        }
        let n = cross(&du, &dv, m)?;
        out.push([n[0].bounds(m)?, n[1].bounds(m)?, n[2].bounds(m)?]);
    }
    Ok(out)
}
fn interval_product(a: &Range, b: &Range, m: &mut Meter) -> Result<Range> {
    let mut out = [m.mul(&a[0], &b[0])?, m.mul(&a[0], &b[0])?];
    for x in a {
        for y in b {
            let p = m.mul(x, y)?;
            if m.cmp(&p, &out[0])?.is_lt() {
                out[0] = p.clone();
            }
            if m.cmp(&p, &out[1])?.is_gt() {
                out[1] = p;
            }
        }
    }
    Ok(out)
}
fn normal_rule(pair: &Pair, m: &mut Meter) -> Result<String> {
    let b = bases(m)?;
    let a = normal_bounds(&pair[0], &b, m)?;
    let b = normal_bounds(&pair[1], &b, m)?;
    for axis in 0..3 {
        let mut positive = true;
        let mut negative = true;
        let (j, k) = ((axis + 1) % 3, (axis + 2) % 3);
        for x in &a {
            for y in &b {
                let p = interval_product(&x[j], &y[k], m)?;
                let n = interval_product(&x[k], &y[j], m)?;
                let lo = m.sub(&p[0], &n[1])?;
                let hi = m.sub(&p[1], &n[0])?;
                positive &= m.cmp(&lo, &q(0))?.is_gt();
                negative &= m.cmp(&hi, &q(0))?.is_lt();
            }
        }
        if positive || negative {
            return Ok(format!("NoClosedLoop(axis={axis})"));
        }
    }
    Err(refuse("r13.normal-cones-not-separated"))
}

fn certify_hill(b: &[Vec<P>], m: &mut Meter) -> Result<()> {
    let mut hs = Vec::new();
    let mut xs = Vec::new();
    for span in b {
        let h = span[2].add(&span[3], m)?;
        h.positive_inside(m)?;
        let mut x = P::zero(1, 1);
        for (i, p) in span.iter().enumerate() {
            x = x.add(&p.scale(&q(i as i64), m)?, m)?;
        }
        x.derivative(0, m)?.positive_inside(m)?;
        hs.push(h);
        xs.push(x);
    }
    for s in 1..3 {
        for ps in [&hs, &xs] {
            let left = ps[s - 1].eval(&q(1), &q(0), m)?;
            let right = ps[s].eval(&q(0), &q(0), m)?;
            m.equal(&left, &right)?;
        }
    }
    for (s, a, z, sign) in [
        (0, q(0), q(1), 1),
        (1, q(0), r("1/2"), 1),
        (1, r("1/2"), q(1), -1),
        (2, q(0), q(1), -1),
    ] {
        hs[s]
            .derivative(0, m)?
            .affine_u(&a, &z, m)?
            .scale(&q(sign), m)?
            .positive_inside(m)?;
    }
    for (s, t, expected) in [(0, q(0), q(0)), (2, q(1), q(0)), (1, r("1/2"), r("15/16"))] {
        let actual = hs[s].eval(&t, &q(0), m)?;
        m.equal(&actual, &expected)?;
    }
    Ok(())
}
fn shared_chart(pair: &Pair, m: &mut Meter) -> Result<()> {
    for net in pair {
        for (i, row) in net.iter().enumerate() {
            for (j, p) in row.iter().enumerate() {
                m.equal(&p[0], &q(i as i64))?;
                m.equal(&p[1], &q(j as i64))?;
            }
        }
    }
    Ok(())
}
fn graph_topology(pair: &Pair, m: &mut Meter) -> Result<String> {
    shared_chart(pair, m)?;
    let b = bases(m)?;
    certify_hill(&b, m)?;
    let difference =
        |i: usize, j: usize, m: &mut Meter| m.sub(&pair[0][i][j][2], &pair[1][i][j][2]);
    let beta = difference(0, 0, m)?;
    let peak_control = difference(2, 2, m)?;
    let alpha = m.sub(&peak_control, &beta)?;
    let mut hill = true;
    for i in 0..6 {
        for j in 0..6 {
            let expected = if (2..=3).contains(&i) && (2..=3).contains(&j) {
                m.add(&alpha, &beta)?
            } else {
                beta.clone()
            };
            let actual = difference(i, j, m)?;
            hill &= m.cmp(&actual, &expected)?.is_eq();
        }
    }
    if hill {
        if m.cmp(&alpha, &q(0))?.is_eq() {
            return Err(refuse("r13.constant-or-coincident-pair"));
        }
        let negative_beta = m.sub(&q(0), &beta)?;
        let c = m.div(&negative_beta, &alpha)?;
        return Ok(match (m.cmp(&c, &q(0))?, m.cmp(&c, &r("225/256"))?) {
            (Ordering::Less, _) | (_, Ordering::Greater) => "Disjoint",
            (_, Ordering::Equal) => "IsolatedContact",
            (Ordering::Greater, Ordering::Less) => "OneInteriorLoop",
            _ => return Err(refuse("r13.boundary-level")),
        }
        .into());
    }
    let next_control = difference(1, 0, m)?;
    let slope = m.sub(&next_control, &beta)?;
    if m.cmp(&slope, &q(0))?.is_eq() {
        return Err(refuse("r13.scalar-topology-not-implemented"));
    }
    for i in 0..6 {
        for j in 0..6 {
            let linear = m.mul(&slope, &q(i as i64))?;
            let expected = m.add(&beta, &linear)?;
            let actual = difference(i, j, m)?;
            if !m.cmp(&actual, &expected)?.is_eq() {
                return Err(refuse("r13.scalar-topology-not-implemented"));
            }
        }
    }
    let negative_beta = m.sub(&q(0), &beta)?;
    let root = m.div(&negative_beta, &slope)?;
    if m.cmp(&root, &q(0))?.is_gt() && m.cmp(&root, &q(5))?.is_lt() {
        Ok("OneOpenArc".into())
    } else {
        Err(refuse("r13.boundary-or-exterior-linear-level"))
    }
}

fn is_zero(p: &P, m: &mut Meter) -> Result<bool> {
    for c in p.0.iter().flatten() {
        if !m.cmp(c, &q(0))?.is_eq() {
            return Ok(false);
        }
    }
    Ok(true)
}
// At most a 5 by 5 Sylvester matrix here. Laplace expansion keeps every
// coefficient operation interruptible; the public MPoly resultant is unbounded.
fn determinant(matrix: &[Vec<P>], m: &mut Meter) -> Result<P> {
    if matrix.len() == 1 {
        return Ok(matrix[0][0].clone());
    }
    let mut out = P::zero(1, 1);
    for (column, p) in matrix[0].iter().enumerate() {
        if is_zero(p, m)? {
            continue;
        }
        let minor: Vec<Vec<_>> = matrix
            .iter()
            .skip(1)
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(j, _)| *j != column)
                    .map(|(_, x)| x.clone())
                    .collect()
            })
            .collect();
        let term = p.mul(&determinant(&minor, m)?, m)?;
        out = if column % 2 == 0 {
            out.add(&term, m)?
        } else {
            out.sub(&term, m)?
        };
    }
    Ok(out)
}
fn resultant_v(a: &P, b: &P, m: &mut Meter) -> Result<P> {
    let coefficients = |p: &P, m: &mut Meter| -> Result<Vec<P>> {
        let mut c: Vec<_> = (0..p.0[0].len())
            .map(|j| P(p.0.iter().map(|row| vec![row[j].clone()]).collect()))
            .collect();
        while c.last().is_some() {
            if !is_zero(c.last().unwrap(), m)? {
                break;
            }
            c.pop();
        }
        Ok(c)
    };
    let a = coefficients(a, m)?;
    let b = coefficients(b, m)?;
    if a.is_empty() || b.is_empty() {
        return Err(refuse("r13.zero-critical-derivative"));
    }
    let (n, d) = (a.len() - 1, b.len() - 1);
    let size = n + d;
    if size == 0 {
        return Err(refuse("r13.constant-resultant-input"));
    }
    assert!(size <= 5);
    let mut matrix = vec![vec![P::zero(1, 1); size]; size];
    for i in 0..d {
        for j in 0..=n {
            matrix[i][i + j] = a[n - j].clone();
        }
    }
    for i in 0..n {
        for j in 0..=d {
            matrix[d + i][i + j] = b[d - j].clone();
        }
    }
    determinant(&matrix, m)
}
fn critical_projection(pair: &Pair, m: &mut Meter) -> Result<String> {
    shared_chart(pair, m)?;
    let b = bases(m)?;
    certify_hill(&b, m)?;
    let a = patches(&pair[0], &b, m)?;
    let b = patches(&pair[1], &b, m)?;
    let mut max_degree = 0;
    for (a, b) in a.iter().zip(&b) {
        let g = a[2].sub(&b[2], m)?;
        let result = resultant_v(&g, &g.derivative(1, m)?, m)?;
        if is_zero(&result, m)? {
            return Err(refuse("r13.nonisolated-critical-projection"));
        }
        for (i, row) in result.0.iter().enumerate() {
            if !m.cmp(&row[0], &q(0))?.is_eq() {
                max_degree = max_degree.max(i);
            }
        }
    }
    Ok(format!(
        "CriticalProjection(cells=9,max-degree={max_degree});topology=unresolved"
    ))
}

// A curve and its endpoint chord in the same strict coordinate half-plane
// are homotopic away from zero. Subdivision plus this certificate makes the
// integer polygon winding exact, without atan2 or sampled-angle assumptions.
fn boundary_segments(
    field: &[P; 2],
    depth: usize,
    points: &mut Vec<[Q; 2]>,
    m: &mut Meter,
) -> Result<()> {
    let a = field[0].bounds(m)?;
    let b = field[1].bounds(m)?;
    let mut separated = false;
    for range in [a, b] {
        separated |= m.cmp(&range[0], &q(0))?.is_gt() || m.cmp(&range[1], &q(0))?.is_lt();
    }
    if separated {
        points.push([
            field[0].eval(&q(1), &q(0), m)?,
            field[1].eval(&q(1), &q(0), m)?,
        ]);
        return Ok(());
    }
    if depth == 12 {
        return Err(refuse("r13.boundary-field-not-separated"));
    }
    for (lo, hi) in [(q(0), r("1/2")), (r("1/2"), q(1))] {
        boundary_segments(
            &[
                field[0].affine_u(&lo, &hi, m)?,
                field[1].affine_u(&lo, &hi, m)?,
            ],
            depth + 1,
            points,
            m,
        )?;
    }
    Ok(())
}
fn winding(pair: &Pair, m: &mut Meter) -> Result<String> {
    shared_chart(pair, m)?;
    let b = bases(m)?;
    let a = patches(&pair[0], &b[1..2], m)?;
    let z = patches(&pair[1], &b[1..2], m)?;
    let g = a[0][2].sub(&z[0][2], m)?;
    let f = [g.derivative(0, m)?, g.derivative(1, m)?];
    let mut points = vec![[f[0].eval(&q(0), &q(0), m)?, f[1].eval(&q(0), &q(0), m)?]];
    for (axis, fixed, reverse) in [
        (0, q(0), false),
        (1, q(1), false),
        (0, q(1), true),
        (1, q(0), true),
    ] {
        let mut edge = [f[0].edge(axis, &fixed, m)?, f[1].edge(axis, &fixed, m)?];
        if reverse {
            edge = [
                edge[0].affine_u(&q(1), &q(0), m)?,
                edge[1].affine_u(&q(1), &q(0), m)?,
            ];
        }
        boundary_segments(&edge, 0, &mut points, m)?;
    }
    let mut degree = 0;
    for edge in points.windows(2) {
        let (a, b) = (&edge[0], &edge[1]);
        let lhs = m.mul(&a[0], &b[1])?;
        let rhs = m.mul(&a[1], &b[0])?;
        let orientation = m.cmp(&lhs, &rhs)?;
        let ay = m.cmp(&a[1], &q(0))?;
        let by = m.cmp(&b[1], &q(0))?;
        if !ay.is_gt() && by.is_gt() && orientation.is_gt() {
            degree += 1;
        }
        if ay.is_gt() && !by.is_gt() && orientation.is_lt() {
            degree -= 1;
        }
    }
    Ok(format!("Index({degree});intersection=unresolved"))
}

#[derive(Clone)]
struct Transform {
    rotation: [[Q; 3]; 3],
    scale: Q,
    translation: [Q; 3],
}
impl Transform {
    fn identity() -> Self {
        Self {
            rotation: std::array::from_fn(|i| std::array::from_fn(|j| q(i64::from(i == j)))),
            scale: q(1),
            translation: [q(0), q(0), q(0)],
        }
    }
    fn forward(&self, p: &[Q; 3]) -> [Q; 3] {
        std::array::from_fn(|i| {
            let rotated = self.rotation[i]
                .iter()
                .zip(p)
                .fold(q(0), |v, (a, b)| v + a * b);
            &(&rotated * &self.scale) + &self.translation[i]
        })
    }
    fn inverse_pair(&self, pair: &Pair, m: &mut Meter) -> Result<Pair> {
        if !m.cmp(&self.scale, &q(0))?.is_gt() {
            return Err(refuse("r13.nonpositive-scale"));
        }
        for i in 0..3 {
            for j in 0..3 {
                let mut dot = q(0);
                for k in 0..3 {
                    let term = m.mul(&self.rotation[k][i], &self.rotation[k][j])?;
                    dot = m.add(&dot, &term)?;
                }
                m.equal(&dot, &q(i64::from(i == j)))?;
            }
        }
        let mut out = pair.clone();
        for net in &mut out {
            for p in net.iter_mut().flatten() {
                let mut unscaled = [q(0), q(0), q(0)];
                for k in 0..3 {
                    let delta = m.sub(&p[k], &self.translation[k])?;
                    unscaled[k] = m.div(&delta, &self.scale)?;
                }
                for (i, x) in p.iter_mut().enumerate() {
                    *x = q(0);
                    for (k, y) in unscaled.iter().enumerate() {
                        let term = m.mul(&self.rotation[k][i], y)?;
                        *x = m.add(x, &term)?;
                    }
                }
            }
        }
        Ok(out)
    }
}
struct Case {
    name: String,
    pair: Pair,
    transform: Transform,
}
fn fixtures() -> Vec<Case> {
    let hill: Net = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| {
                    [
                        q(i),
                        q(j),
                        q(i64::from((2..=3).contains(&i) && (2..=3).contains(&j))),
                    ]
                })
                .collect()
        })
        .collect();
    let mut pairs = Vec::new();
    pairs.push(("m25-wavy", super::wavy_fixture::pair()));
    for (name, c) in [
        ("loop", r("3/10")),
        ("contact", r("225/256")),
        ("small-loop", r("225/256") - r("1/1000000")),
        ("disjoint", r("225/256") + r("1/1000000")),
    ] {
        let mut other = hill.clone();
        for p in other.iter_mut().flatten() {
            p[2] = (&p[2] + &c) / q(2);
        }
        pairs.push((name, [hill.clone(), other]));
    }
    let mut other = hill.clone();
    for p in other.iter_mut().flatten() {
        p[2] = &p[2] + &(&q(10) * &(&p[0] - &r("5/2")));
    }
    pairs.push(("open-control", [hill, other]));
    let mut moved = Transform::identity();
    moved.translation = [r("17/3"), q(-11), r("2/7")];
    let mut turned = Transform::identity();
    turned.rotation = [
        [r("3/5"), q(0), r("4/5")],
        [q(0), q(1), q(0)],
        [r("-4/5"), q(0), r("3/5")],
    ];
    let mut scaled = Transform::identity();
    scaled.scale = r("7/3");
    let mut cases = Vec::new();
    for (name, pair) in pairs {
        for (suffix, transform) in [
            ("base", Transform::identity()),
            ("moved", moved.clone()),
            ("turned", turned.clone()),
            ("scaled", scaled.clone()),
        ] {
            let mut physical = pair.clone();
            for net in &mut physical {
                for p in net.iter_mut().flatten() {
                    *p = transform.forward(p);
                }
            }
            cases.push(Case {
                name: format!("{name}/{suffix}"),
                pair: physical,
                transform,
            });
        }
    }
    cases
}
#[derive(Debug, PartialEq, Eq)]
struct Observation {
    answer: String,
    work: u64,
    detail: String,
}
fn render(result: Result<String>) -> String {
    match result {
        Ok(value) => value,
        Err(refusal) => format!("Refused({})", refusal.code),
    }
}
fn observe(
    case: &Case,
    limit: u64,
    method: fn(&Pair, &mut Meter) -> Result<String>,
) -> Observation {
    let (result, work) = with_operation_budget(limit, || {
        let mut m = Meter::default();
        let canonical = case.transform.inverse_pair(&case.pair, &mut m)?;
        method(&canonical, &mut m)
    });
    let detail = result
        .as_ref()
        .err()
        .map(|r| r.detail.clone())
        .unwrap_or_default();
    let answer = if work.exhausted {
        "Refused(numbers.operation-budget)".into()
    } else {
        render(result)
    };
    Observation {
        answer,
        work: work.used,
        detail,
    }
}
#[derive(Debug, PartialEq, Eq)]
struct Row {
    name: String,
    cone: Observation,
    graph: Observation,
    projection: Observation,
    winding: Observation,
    paper: Observation,
}
fn experiment(workers: usize) -> Vec<Row> {
    let cases = fixtures();
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|worker| {
                let cases = &cases;
                scope.spawn(move || {
                    (worker..cases.len())
                        .step_by(workers)
                        .map(|i| {
                            let c = &cases[i];
                            (
                                i,
                                Row {
                                    name: c.name.clone(),
                                    cone: observe(c, ALLOWANCE, normal_rule),
                                    graph: observe(c, ALLOWANCE, graph_topology),
                                    projection: observe(c, ALLOWANCE, critical_projection),
                                    winding: observe(c, ALLOWANCE, winding),
                                    paper: observe(c, ALLOWANCE, paper2026::experiment),
                                },
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        let mut rows: Vec<_> = handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect();
        rows.sort_by_key(|(i, _)| *i);
        rows.into_iter().map(|(_, row)| row).collect()
    })
}
pub(super) fn run() {
    let rows = experiment(1);
    for workers in [2, 4] {
        assert_eq!(rows, experiment(workers));
    }
    println!(
        "case|normal cone|work|scalar family|work|2023 critical projection|work|gradient winding control|work|2026 paper method|work"
    );
    for row in rows {
        println!(
            "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
            row.name,
            row.cone.answer,
            row.cone.work,
            row.graph.answer,
            row.graph.work,
            row.projection.answer,
            row.projection.work,
            row.winding.answer,
            row.winding.work,
            row.paper.answer,
            row.paper.work
        );
        println!("2026-evidence {} {}", row.name, row.paper.detail);
    }
    println!(
        "24 spline pairs; allowance=20000000 per method; workers=1/2/4 identical answers, work and evidence; both paper cores checked against full text; no complete SSI implementation"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn families_transforms_and_workers() {
        let rows = experiment(1);
        assert_eq!(rows.len(), 24);
        assert_eq!(rows, experiment(2));
        assert_eq!(rows, experiment(4));
        for row in rows {
            let family = row.name.split('/').next().unwrap();
            let expected = match family {
                "loop" | "small-loop" => "OneInteriorLoop",
                "contact" => "IsolatedContact",
                "disjoint" => "Disjoint",
                "open-control" => "OneOpenArc",
                _ => "Refused(r13.scalar-topology-not-implemented)",
            };
            assert_eq!(row.graph.answer, expected, "{}", row.name);
            if family == "open-control" {
                assert!(row.cone.answer.starts_with("NoClosedLoop"));
            } else {
                assert_eq!(row.cone.answer, "Refused(r13.normal-cones-not-separated)");
            }
            if matches!(family, "loop" | "small-loop" | "contact" | "disjoint") {
                assert_eq!(row.winding.answer, "Index(1);intersection=unresolved");
            }
            assert_eq!(
                row.paper.answer,
                "Refused(r13.distance-field-coverage-unproved)"
            );
            assert!(row.paper.work > 0);
            if matches!(family, "loop" | "small-loop" | "contact" | "disjoint") {
                assert!(row.paper.detail.contains("peak=exact;closest=unique"));
                let starts = match family {
                    "loop" | "small-loop" => 4,
                    "contact" => 1,
                    _ => 0,
                };
                assert!(row.paper.detail.contains(&format!("starts={starts}")));
            } else {
                assert!(
                    row.paper
                        .detail
                        .contains("peak=Refused(r13.identity-not-proved)")
                );
            }
            if family == "open-control" {
                assert_eq!(
                    row.projection.answer,
                    "Refused(r13.zero-critical-derivative)"
                );
            } else {
                let degree = if family == "m25-wavy" { 15 } else { 12 };
                assert_eq!(
                    row.projection.answer,
                    format!("CriticalProjection(cells=9,max-degree={degree});topology=unresolved")
                );
            }
            for value in [row.cone, row.graph, row.projection, row.winding, row.paper] {
                assert!(value.work < ALLOWANCE);
            }
        }
    }
    #[test]
    fn exhausted_and_invalid_inputs_refuse() {
        let mut case = fixtures().remove(4);
        let exhausted = observe(&case, 1, normal_rule);
        assert_eq!(exhausted.answer, "Refused(numbers.operation-budget)");
        assert!(exhausted.work <= 1);
        assert_eq!(
            observe(&case, 1, paper2026::experiment).answer,
            "Refused(numbers.operation-budget)"
        );
        case.pair[1][0][0][0] = q(-1);
        assert_eq!(
            observe(&case, ALLOWANCE, graph_topology).answer,
            "Refused(r13.identity-not-proved)"
        );
        case.transform.scale = q(0);
        assert_eq!(
            observe(&case, ALLOWANCE, normal_rule).answer,
            "Refused(r13.nonpositive-scale)"
        );
    }
    #[test]
    fn winding_refuses_boundary_zero() {
        let mut m = Meter::default();
        let field = [P::line(q(-1), q(2)), P::zero(1, 1)];
        assert_eq!(
            boundary_segments(&field, 0, &mut Vec::new(), &mut m)
                .unwrap_err()
                .code,
            "r13.boundary-field-not-separated"
        );
    }
    #[test]
    fn critical_resultant_matches_quadratic_identity() {
        let mut m = Meter::default();
        let g = P(vec![
            vec![q(-1), q(0), q(1)],
            vec![q(0), q(0), q(0)],
            vec![q(1), q(0), q(0)],
        ]);
        let actual = resultant_v(&g, &g.derivative(1, &mut m).unwrap(), &mut m).unwrap();
        let expected = P(vec![vec![q(-4)], vec![q(0)], vec![q(4)]]);
        assert!(is_zero(&actual.sub(&expected, &mut m).unwrap(), &mut m).unwrap());
    }
    #[test]
    fn wavy_resultants_match_kernel_elimination() {
        use catallax_numbers::MPoly;
        fn polynomial(p: &P) -> MPoly {
            MPoly::from_terms(
                2,
                p.0.iter()
                    .enumerate()
                    .flat_map(|(i, row)| {
                        row.iter()
                            .enumerate()
                            .map(move |(j, x)| (vec![i as u32, j as u32], x.clone()))
                    })
                    .collect(),
            )
        }
        let mut m = Meter::default();
        let b = bases(&mut m).unwrap();
        let pair = super::super::wavy_fixture::pair();
        let a = patches(&pair[0], &b, &mut m).unwrap();
        let b = patches(&pair[1], &b, &mut m).unwrap();
        for (a, b) in a.iter().zip(&b) {
            let g = a[2].sub(&b[2], &mut m).unwrap();
            let derivative = g.derivative(1, &mut m).unwrap();
            let actual = resultant_v(&g, &derivative, &mut m).unwrap();
            let expected = polynomial(&g).resultant(&polynomial(&derivative), 1);
            assert_eq!(polynomial(&actual), expected);
        }
    }
}
