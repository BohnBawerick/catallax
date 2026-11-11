//! Exact Bernstein nets of one patch over a parameter cell, with the bounds the
//! cell-pair certificates read. Every rational operation reserves its exact
//! work before it runs.

use super::*;
use std::sync::OnceLock;

/// Charged exact rational arithmetic.
pub(super) fn add(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a + b)
}
pub(super) fn sub(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a - b)
}
pub(super) fn mul(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a * b)
}
/// Comparing fractions cross-multiplies, so it reserves two products.
pub(super) fn less(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<bool, Refusal> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a < b)
}
fn lerp(a: &Rational, b: &Rational, t: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    let step = sub(b, a, dec)?;
    let step = mul(&step, t, dec)?;
    add(a, &step, dec)
}

/// Dense power-basis coefficients in local unit-square parameters,
/// `c[i * (deg[1] + 1) + j]` for `a^i b^j`.
#[derive(Clone)]
pub(super) struct Dense {
    deg: [usize; 2],
    c: Vec<Rational>,
}

impl Dense {
    fn zero(deg: [usize; 2]) -> Self {
        Self {
            deg,
            c: vec![Rational::zero(); (deg[0] + 1) * (deg[1] + 1)],
        }
    }
    fn at(&self, i: usize, j: usize) -> &Rational {
        &self.c[i * (self.deg[1] + 1) + j]
    }
    fn at_mut(&mut self, i: usize, j: usize) -> &mut Rational {
        let stride = self.deg[1] + 1;
        &mut self.c[i * stride + j]
    }

    /// `poly` in global parameters, rewritten for `u = u0 + (u1-u0) a`,
    /// `v = v0 + (v1-v0) b` on the given rectangle.
    pub(super) fn local(poly: &MPoly, domain: &Domain, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut deg = [0usize; 2];
        for (monomial, _) in poly.terms() {
            for (axis, d) in deg.iter_mut().enumerate() {
                *d = (*d).max(monomial.exponents()[axis] as usize);
            }
        }
        // powers[axis][i] holds the coefficients of (lo + h x)^i.
        let mut powers: [Vec<Vec<Rational>>; 2] = [Vec::new(), Vec::new()];
        for axis in 0..2 {
            let lo = &domain[axis][0];
            let h = sub(&domain[axis][1], lo, dec)?;
            let mut current = vec![Rational::one()];
            powers[axis].push(current.clone());
            for _ in 0..deg[axis] {
                let mut next = vec![Rational::zero(); current.len() + 1];
                for (k, value) in current.iter().enumerate() {
                    let low = mul(value, lo, dec)?;
                    next[k] = add(&next[k], &low, dec)?;
                    let high = mul(value, &h, dec)?;
                    next[k + 1] = add(&next[k + 1], &high, dec)?;
                }
                current = next;
                powers[axis].push(current.clone());
            }
        }
        let mut out = Self::zero(deg);
        for (monomial, coefficient) in poly.terms() {
            let [i, j] = [
                monomial.exponents()[0] as usize,
                monomial.exponents()[1] as usize,
            ];
            for (k, pu) in powers[0][i].iter().enumerate() {
                let scaled = mul(&coefficient, pu, dec)?;
                for (l, pv) in powers[1][j].iter().enumerate() {
                    let term = mul(&scaled, pv, dec)?;
                    let slot = out.at(k, l).clone();
                    *out.at_mut(k, l) = add(&slot, &term, dec)?;
                }
            }
        }
        Ok(out.trimmed())
    }

    /// The same polynomial with zero coefficients up to degree `deg`.
    fn padded(&self, deg: [usize; 2]) -> Self {
        let mut out = Self::zero(deg);
        for i in 0..=self.deg[0] {
            for j in 0..=self.deg[1] {
                *out.at_mut(i, j) = self.at(i, j).clone();
            }
        }
        out
    }

    /// Drops exactly zero highest-degree rows and columns.
    fn trimmed(mut self) -> Self {
        loop {
            let [m, n] = self.deg;
            if m > 0 && (0..=n).all(|j| self.at(m, j).is_zero()) {
                let mut next = Self::zero([m - 1, n]);
                for i in 0..m {
                    for j in 0..=n {
                        *next.at_mut(i, j) = self.at(i, j).clone();
                    }
                }
                self = next;
                continue;
            }
            if n > 0 && (0..=m).all(|i| self.at(i, n).is_zero()) {
                let mut next = Self::zero([m, n - 1]);
                for i in 0..=m {
                    for j in 0..n {
                        *next.at_mut(i, j) = self.at(i, j).clone();
                    }
                }
                self = next;
                continue;
            }
            return self;
        }
    }

    pub(super) fn mul(&self, other: &Self, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut out = Self::zero([self.deg[0] + other.deg[0], self.deg[1] + other.deg[1]]);
        for i in 0..=self.deg[0] {
            for j in 0..=self.deg[1] {
                let a = self.at(i, j);
                if a.is_zero() {
                    continue;
                }
                for k in 0..=other.deg[0] {
                    for l in 0..=other.deg[1] {
                        let b = other.at(k, l);
                        if b.is_zero() {
                            continue;
                        }
                        let term = mul(a, b, dec)?;
                        let slot = out.at(i + k, j + l).clone();
                        *out.at_mut(i + k, j + l) = add(&slot, &term, dec)?;
                    }
                }
            }
        }
        Ok(out.trimmed())
    }

    pub(super) fn sub(&self, other: &Self, dec: &mut Decider) -> Result<Self, Refusal> {
        let deg = [self.deg[0].max(other.deg[0]), self.deg[1].max(other.deg[1])];
        let mut out = Self::zero(deg);
        for i in 0..=deg[0] {
            for j in 0..=deg[1] {
                let a = (i <= self.deg[0] && j <= self.deg[1]).then(|| self.at(i, j));
                let b = (i <= other.deg[0] && j <= other.deg[1]).then(|| other.at(i, j));
                *out.at_mut(i, j) = match (a, b) {
                    (Some(a), Some(b)) if !b.is_zero() => sub(a, b, dec)?,
                    (Some(a), _) => a.clone(),
                    (None, Some(b)) => -b.clone(),
                    (None, None) => Rational::zero(),
                };
            }
        }
        Ok(out.trimmed())
    }

    pub(super) fn derivative(&self, axis: usize, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut deg = self.deg;
        if deg[axis] == 0 {
            return Ok(Self::zero([0, 0]));
        }
        deg[axis] -= 1;
        let mut out = Self::zero(deg);
        for i in 0..=deg[0] {
            for j in 0..=deg[1] {
                let (k, l, n) = if axis == 0 {
                    (i + 1, j, i + 1)
                } else {
                    (i, j + 1, j + 1)
                };
                *out.at_mut(i, j) = mul(self.at(k, l), &Rational::from_i64(n as i64), dec)?;
            }
        }
        Ok(out.trimmed())
    }

    /// Tensor Bernstein coefficients on the unit square.
    pub(super) fn bernstein(&self, dec: &mut Decider) -> Result<Grid, Refusal> {
        let [m, n] = self.deg;
        let mut values = self.c.clone();
        // b_i = sum_{k <= i} C(i,k) / C(m,k) a_k along each axis.
        for axis in 0..2 {
            let d = self.deg[axis];
            let lines = if axis == 0 { n + 1 } else { m + 1 };
            for line in 0..lines {
                let index = |k: usize| {
                    if axis == 0 {
                        k * (n + 1) + line
                    } else {
                        line * (n + 1) + k
                    }
                };
                let source: Vec<Rational> = (0..=d).map(|k| values[index(k)].clone()).collect();
                for i in 0..=d {
                    let mut total = Rational::zero();
                    for (k, a) in source.iter().enumerate().take(i + 1) {
                        if a.is_zero() {
                            continue;
                        }
                        let factor = Rational::from_ratio(choose(i, k), choose(d, k))
                            .expect("positive binomial");
                        let term = mul(a, &factor, dec)?;
                        total = add(&total, &term, dec)?;
                    }
                    values[index(i)] = total;
                }
            }
        }
        Ok(Grid {
            deg: self.deg,
            v: values,
        })
    }
}

fn choose(n: usize, k: usize) -> i64 {
    (0..k).fold(1i64, |a, i| a * (n - i) as i64 / (i + 1) as i64)
}

/// A Bernstein coefficient.
pub(super) trait Value: Clone + std::fmt::Debug {
    /// `a + t (b - a)`.
    fn lerp(a: &Self, b: &Self, t: &Rational, dec: &mut Decider) -> Result<Self, Refusal>;
    /// Whether two values are equal, at the cost of one comparison.
    fn same(a: &Self, b: &Self, dec: &mut Decider) -> Result<bool, Refusal>;
    /// Rational bounds of the value.
    fn bounds(&self) -> [Rational; 2];
}

impl Value for Rational {
    fn lerp(a: &Self, b: &Self, t: &Rational, dec: &mut Decider) -> Result<Self, Refusal> {
        lerp(a, b, t, dec)
    }
    fn same(a: &Self, b: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        dec.reserve_exact_arithmetic(a, b, false)?;
        Ok(a == b)
    }
    fn bounds(&self) -> [Rational; 2] {
        [self.clone(), self.clone()]
    }
}

/// Tensor Bernstein coefficients over one cell, `v[i * (deg[1] + 1) + j]`.
#[derive(Clone, Debug)]
pub(super) struct Grid<T = Rational> {
    pub(super) deg: [usize; 2],
    pub(super) v: Vec<T>,
}

impl Grid<Rational> {
    fn all_positive(&self) -> bool {
        self.v.iter().all(|value| value.sign() == Sign::Positive)
    }
}

impl<T: Value> Grid<T> {
    /// de Casteljau at local parameter `t` along `axis`: the nets on `[0,t]`
    /// and `[t,1]`.
    fn split(&self, axis: usize, t: &Rational, dec: &mut Decider) -> Result<[Self; 2], Refusal> {
        let [m, n] = self.deg;
        let d = self.deg[axis];
        let mut low = self.clone();
        let mut high = self.clone();
        let lines = if axis == 0 { n + 1 } else { m + 1 };
        for line in 0..lines {
            let index = |k: usize| {
                if axis == 0 {
                    k * (n + 1) + line
                } else {
                    line * (n + 1) + k
                }
            };
            let mut work: Vec<T> = (0..=d).map(|k| self.v[index(k)].clone()).collect();
            // A constant line, such as a constant weight or a coordinate that
            // does not change along `axis`, splits to itself.
            let mut constant = true;
            for k in 1..=d {
                if !T::same(&work[0], &work[k], dec)? {
                    constant = false;
                    break;
                }
            }
            if constant {
                continue;
            }
            low.v[index(0)] = work[0].clone();
            high.v[index(d)] = work[d].clone();
            for level in 1..=d {
                for k in 0..=d - level {
                    work[k] = T::lerp(&work[k], &work[k + 1], t, dec)?;
                }
                low.v[index(level)] = work[0].clone();
                high.v[index(d - level)] = work[d - level].clone();
            }
        }
        Ok([low, high])
    }

    fn range(&self, dec: &mut Decider) -> Result<[Rational; 2], Refusal> {
        let mut out = self.v[0].bounds();
        for value in &self.v[1..] {
            let [lo, hi] = value.bounds();
            if less(&lo, &out[0], dec)? {
                out[0] = lo;
            }
            if less(&out[1], &hi, dec)? {
                out[1] = hi;
            }
        }
        Ok(out)
    }
}

/// Splits each grid at `t`, keeping the low parts and leaving the high.
fn take_low<T: Value>(
    grids: &mut [Grid<T>],
    axis: usize,
    t: &Rational,
    dec: &mut Decider,
) -> Result<Vec<Grid<T>>, Refusal> {
    let mut taken = Vec::with_capacity(grids.len());
    for grid in grids.iter_mut() {
        let [low, high] = grid.split(axis, t, dec)?;
        taken.push(low);
        *grid = high;
    }
    Ok(taken)
}

fn values_range<'a>(
    mut values: impl Iterator<Item = &'a Rational>,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let first = values.next().expect("nonempty net").clone();
    let mut range = [first.clone(), first];
    for value in values {
        if less(value, &range[0], dec)? {
            range[0] = value.clone();
        } else if less(&range[1], value, dec)? {
            range[1] = value.clone();
        }
    }
    Ok(range)
}

/// Closed boxes, one interval per axis.
pub(super) type Bounds3 = [[Rational; 2]; 3];

/// The tangent and normal nets of one cell and their ranges.
#[derive(Clone, Debug)]
pub(super) struct Derived {
    tu: [Grid; 3],
    tv: [Grid; 3],
    normal: [Grid; 3],
    /// Ranges of the homogeneous tangent and normal components.
    pub(super) tu_range: Bounds3,
    pub(super) tv_range: Bounds3,
    pub(super) normal_range: Bounds3,
}

impl Derived {
    fn new(
        tu: [Grid; 3],
        tv: [Grid; 3],
        normal: [Grid; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let ranges = |grids: &[Grid; 3], dec: &mut Decider| -> Result<Bounds3, Refusal> {
            Ok([
                grids[0].range(dec)?,
                grids[1].range(dec)?,
                grids[2].range(dec)?,
            ])
        };
        Ok(Self {
            tu_range: ranges(&tu, dec)?,
            tv_range: ranges(&tv, dec)?,
            normal_range: ranges(&normal, dec)?,
            tu,
            tv,
            normal,
        })
    }

    /// The nets on a sub-rectangle given by fractions of this one's sides.
    fn restricted(&self, parts: &[[Rational; 2]; 2], dec: &mut Decider) -> Result<Self, Refusal> {
        let mut grids: Vec<Grid> = self
            .tu
            .iter()
            .chain(&self.tv)
            .chain(&self.normal)
            .cloned()
            .collect();
        for (axis, [lo, hi]) in parts.iter().enumerate() {
            if !lo.is_zero() {
                for grid in &mut grids {
                    *grid = grid.split(axis, lo, dec)?[1].clone();
                }
            }
            if *hi != Rational::one() {
                let width = sub(&Rational::one(), lo, dec)?;
                let offset = sub(hi, lo, dec)?;
                dec.reserve_exact_arithmetic(&offset, &width, true)?;
                let t = &offset / &width;
                for grid in &mut grids {
                    *grid = grid.split(axis, &t, dec)?[0].clone();
                }
            }
        }
        let mut it = grids.into_iter();
        let mut three = || -> [Grid; 3] { std::array::from_fn(|_| it.next().expect("nine grids")) };
        let tu = three();
        let tv = three();
        let normal = three();
        Self::new(tu, tv, normal, dec)
    }
}

/// The nets of one patch over one cell of one knot span, in that span's
/// local parameters, and the bounds that exclusion reads.
#[derive(Debug)]
pub(super) struct Cell {
    pub(super) span: usize,
    /// Global parameter rectangle.
    pub(super) rect: Domain,
    /// Local rectangle inside the span's unit square.
    pub(super) local: Domain,
    pub(super) depth: u32,
    unit_weight: bool,
    h: [Grid; 4],
    /// Box round the patch over the cell, from outward boxes round its
    /// controls.
    pub(super) bbox: Bounds3,
    /// Boxes round the four border curves: low u, high u, low v, high v.
    pub(super) edges: [Bounds3; 4],
    /// Projected controls, `points[i * (degree[1] + 1) + j]`.
    points: Vec<[Rational; 3]>,
    /// Outward dyadic boxes round them.
    coarse: Vec<[[Rational; 2]; 3]>,
    /// The exact box of the controls, for the touch proof only.
    exact_box: OnceLock<Bounds3>,
    degree: [usize; 2],
    /// A small integer direction near the cell's normal, from its corner
    /// controls. Slabs along any direction keep exclusion sound; this one is
    /// tight where the two patches nearly touch.
    pub(super) frame: [[Rational; 3]; 1],
    /// The cell's own slab along it.
    pub(super) own: [[Rational; 2]; 1],
    /// Binary64 copies of the boxes, the direction and the slab. They only
    /// skip exact slab tests that would clearly fail.
    rough: Rough,
    /// The largest side of its box.
    pub(super) extent: Rational,
    /// Tangent and normal nets, filled by the request's own thread when the
    /// cell is first proved against something.
    derived: OnceLock<Derived>,
    /// The cell this one was split from.
    parent: Option<Arc<Cell>>,
    /// The second patch's canonical children, filled once by the writer.
    children: OnceLock<Vec<Arc<Cell>>>,
}

impl Cell {
    /// The whole knot span, from its exact global homogeneous polynomials.
    pub(super) fn span(
        index: usize,
        domain: &Domain,
        homogeneous: &[MPoly; 4],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut local = Vec::with_capacity(4);
        for poly in homogeneous {
            local.push(Dense::local(poly, domain, dec)?);
        }
        // The four homogeneous coordinates share one net, so one degree.
        let deg = local
            .iter()
            .fold([0, 0], |d, p| [d[0].max(p.deg[0]), d[1].max(p.deg[1])]);
        let local: Vec<Dense> = local.iter().map(|p| p.padded(deg)).collect();
        let local: [Dense; 4] = local.try_into().ok().expect("four coordinates");
        let mut tangents: Vec<[Dense; 3]> = Vec::with_capacity(2);
        for axis in 0..2 {
            let dw = local[3].derivative(axis, dec)?;
            let mut parts = Vec::with_capacity(3);
            for coordinate in local.iter().take(3) {
                let dx = coordinate.derivative(axis, dec)?;
                let left = dx.mul(&local[3], dec)?;
                let right = coordinate.mul(&dw, dec)?;
                parts.push(left.sub(&right, dec)?);
            }
            tangents.push(parts.try_into().ok().expect("three components"));
        }
        let [tu, tv]: [[Dense; 3]; 2] = tangents.try_into().ok().expect("two tangents");
        let mut normal = Vec::with_capacity(3);
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let left = tu[j].mul(&tv[k], dec)?;
            let right = tu[k].mul(&tv[j], dec)?;
            normal.push(left.sub(&right, dec)?);
        }
        let mut grids = |polys: &[Dense]| -> Result<Vec<Grid>, Refusal> {
            polys.iter().map(|p| p.bernstein(dec)).collect()
        };
        let h: [Grid; 4] = grids(&local)?.try_into().expect("four");
        let tu: [Grid; 3] = grids(&tu)?.try_into().expect("three");
        let tv: [Grid; 3] = grids(&tv)?.try_into().expect("three");
        let normal: [Grid; 3] = grids(&normal)?.try_into().expect("three");
        if !h[3].all_positive() {
            return Err(Refusal::new(
                "intersect.spline-pair-weights",
                "a knot span has a homogeneous weight control that is not positive",
            ));
        }
        let unit_weight = h[3].v.iter().all(|w| *w == Rational::one());
        let unit = [Rational::zero(), Rational::one()];
        let cell = Self::finish(
            index,
            domain.clone(),
            [unit.clone(), unit],
            0,
            unit_weight,
            h,
            None,
            dec,
        )?;
        let derived = Derived::new(tu, tv, normal, dec)?;
        let _ = cell.derived.set(derived);
        Ok(cell)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        span: usize,
        rect: Domain,
        local: Domain,
        depth: u32,
        unit_weight: bool,
        h: [Grid; 4],
        parent: Option<Arc<Cell>>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let [nu, nv] = h[0].deg;
        let mut points = Vec::with_capacity((nu + 1) * (nv + 1));
        for k in 0..h[0].v.len() {
            let mut point: [Rational; 3] = std::array::from_fn(|i| h[i].v[k].clone());
            if !unit_weight {
                let w = &h[3].v[k];
                for value in &mut point {
                    dec.reserve_exact_arithmetic(value, w, true)?;
                    *value = &*value / w;
                }
            }
            points.push(point);
        }
        // Exclusion reads outward dyadic boxes round the controls, whose size
        // does not grow with depth; exact controls stay for the touch proof.
        let mut coarse: Vec<[[Rational; 2]; 3]> = Vec::with_capacity(points.len());
        for point in &points {
            let mut boxed: Vec<[Rational; 2]> = Vec::with_capacity(3);
            for value in point {
                boxed.push(round_out(value, dec)?);
            }
            coarse.push(boxed.try_into().expect("three axes"));
        }
        let box_of = |indices: &mut dyn Iterator<Item = usize>,
                      dec: &mut Decider|
         -> Result<Bounds3, Refusal> {
            let chosen: Vec<usize> = indices.collect();
            let mut out: Vec<[Rational; 2]> = Vec::with_capacity(3);
            for (axis, first) in coarse[chosen[0]].iter().enumerate() {
                let mut range = first.clone();
                for &k in &chosen[1..] {
                    let [lo, hi] = &coarse[k][axis];
                    if less(lo, &range[0], dec)? {
                        range[0] = lo.clone();
                    }
                    if less(&range[1], hi, dec)? {
                        range[1] = hi.clone();
                    }
                }
                out.push(range);
            }
            Ok(out.try_into().expect("three axes"))
        };
        let stride = nv + 1;
        let bbox = box_of(&mut (0..points.len()), dec)?;
        let edges = [
            box_of(&mut (0..=nv), dec)?,
            box_of(&mut (nu * stride..=nu * stride + nv), dec)?,
            box_of(&mut (0..=nu).map(|i| i * stride), dec)?,
            box_of(&mut (0..=nu).map(|i| i * stride + nv), dec)?,
        ];
        // The frame comes from the corner controls' lower bounds: any
        // directions are sound, so approximate corners do.
        let corner =
            |k: usize| -> [Rational; 3] { std::array::from_fn(|a| coarse[k][a][0].clone()) };
        let [p00, p10, p01, p11] = [
            corner(0),
            corner(nu * stride),
            corner(nv),
            corner(nu * stride + nv),
        ];
        let mut du: [Rational; 3] = std::array::from_fn(|_| Rational::zero());
        let mut dv: [Rational; 3] = std::array::from_fn(|_| Rational::zero());
        for a in 0..3 {
            let x = sub(&p10[a], &p00[a], dec)?;
            let y = sub(&p11[a], &p01[a], dec)?;
            du[a] = add(&x, &y, dec)?;
            let x = sub(&p01[a], &p00[a], dec)?;
            let y = sub(&p11[a], &p10[a], dec)?;
            dv[a] = add(&x, &y, dec)?;
        }
        let du = small_direction(&du, dec)?;
        let dv = small_direction(&dv, dec)?;
        let n = small_direction(&cross(&du, &dv, dec)?, dec)?;
        let frame = [n];
        let own = [slab_of(&coarse, &frame[0], None, [nu, nv], dec)?];
        let rough = Rough {
            boxes: coarse
                .iter()
                .map(|b| {
                    b.each_ref()
                        .map(|r| r.each_ref().map(Rational::to_f64_nearest))
                })
                .collect(),
            frame: frame[0].each_ref().map(Rational::to_f64_nearest),
            own: own[0].each_ref().map(Rational::to_f64_nearest),
        };
        let mut extent = Rational::zero();
        for [lo, hi] in &bbox {
            let side = sub(hi, lo, dec)?;
            if less(&extent, &side, dec)? {
                extent = side;
            }
        }
        Ok(Self {
            span,
            rect,
            local,
            depth,
            unit_weight,
            h,
            bbox,
            edges,
            points,
            coarse,
            exact_box: OnceLock::new(),
            degree: [nu, nv],
            frame,
            own,
            rough,
            extent,
            derived: OnceLock::new(),
            parent,
            children: OnceLock::new(),
        })
    }

    /// The tangent and normal nets; the writer filled them first.
    pub(super) fn derived(&self) -> &Derived {
        self.derived
            .get()
            .expect("the writer fills derived nets before proving a cell")
    }

    /// Fills the tangent and normal nets from the nearest ancestor that has
    /// them. Only the request's own thread calls this, in a fixed order.
    pub(super) fn ensure_derived(&self, dec: &mut Decider) -> Result<(), Refusal> {
        if self.derived.get().is_some() {
            return Ok(());
        }
        writer_only()?;
        let mut ancestor = self.parent.as_ref().expect("a span cell has derived nets");
        while ancestor.derived.get().is_none() {
            ancestor = ancestor
                .parent
                .as_ref()
                .expect("a span cell has derived nets");
        }
        let mut parts: Vec<[Rational; 2]> = Vec::with_capacity(2);
        for axis in 0..2 {
            let [a0, a1] = &ancestor.local[axis];
            let width = sub(a1, a0, dec)?;
            let mut part = Vec::with_capacity(2);
            for end in &self.local[axis] {
                let offset = sub(end, a0, dec)?;
                dec.reserve_exact_arithmetic(&offset, &width, true)?;
                part.push(&offset / &width);
            }
            parts.push(part.try_into().expect("two ends"));
        }
        let parts: [[Rational; 2]; 2] = parts.try_into().expect("two axes");
        let derived = ancestor.derived().restricted(&parts, dec)?;
        let _ = self.derived.set(derived);
        Ok(())
    }

    fn edge_indices(&self, edge: usize) -> Vec<usize> {
        edge_indices(self.degree, edge)
    }

    /// Whether, in binary64, this cell's slab (or one border's) along
    /// `other`'s direction overlaps `other`'s own slab by far more than any
    /// rounding, so the exact test could not separate them. A `false`
    /// decides nothing; the exact test then runs.
    pub(super) fn clearly_meets(
        &self,
        other: &Cell,
        edge: Option<usize>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let indices: Vec<usize> = match edge {
            Some(edge) => self.edge_indices(edge),
            None => (0..self.rough.boxes.len()).collect(),
        };
        dec.charge_work(indices.len() as u64 * 3)?;
        let d = other.rough.frame;
        let mut range = [f64::INFINITY, f64::NEG_INFINITY];
        let mut size = other.rough.own[0].abs().max(other.rough.own[1].abs());
        for k in indices {
            let mut total = [0.0f64; 2];
            let mut magnitude = 0.0f64;
            for (axis, [lo, hi]) in self.rough.boxes[k].iter().enumerate() {
                let (low, high) = if d[axis] >= 0.0 {
                    (*lo, *hi)
                } else {
                    (*hi, *lo)
                };
                total[0] += low * d[axis];
                total[1] += high * d[axis];
                magnitude += d[axis].abs() * lo.abs().max(hi.abs());
            }
            range = [range[0].min(total[0]), range[1].max(total[1])];
            size = size.max(magnitude);
        }
        let overlap = range[1].min(other.rough.own[1]) - range[0].max(other.rough.own[0]);
        Ok(overlap.is_finite() && overlap > size * (1.0 / (1u64 << 44) as f64))
    }

    /// The range of `direction . P` over the patch on the cell, or on one of
    /// its borders.
    pub(super) fn slab(
        &self,
        direction: &[Rational; 3],
        edge: Option<usize>,
        dec: &mut Decider,
    ) -> Result<[Rational; 2], Refusal> {
        slab_of(&self.coarse, direction, edge, self.degree, dec)
    }

    /// The exact box of the controls. Only the request's own thread asks,
    /// during cover refinement, so its one charge falls in a fixed order.
    pub(super) fn exact_box(&self, dec: &mut Decider) -> Result<&Bounds3, Refusal> {
        if let Some(found) = self.exact_box.get() {
            return Ok(found);
        }
        writer_only()?;
        let mut out: Vec<[Rational; 2]> = Vec::with_capacity(3);
        for axis in 0..3 {
            out.push(values_range(self.points.iter().map(|p| &p[axis]), dec)?);
        }
        let out: Bounds3 = out.try_into().expect("three axes");
        Ok(self.exact_box.get_or_init(|| out))
    }

    /// Whether the patch reaches its box's high (or low) side on `axis` only
    /// on one of the given borders of the cell: every control off that
    /// border lies strictly inside. The patch is a positive combination of
    /// its controls, so it then meets that side nowhere else.
    pub(super) fn extreme_on_border(
        &self,
        axis: usize,
        high: bool,
        edges: &[usize],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        if edges.is_empty() {
            return Ok(false);
        }
        let value = self.exact_box(dec)?[axis][usize::from(high)].clone();
        let value = &value;
        for &edge in edges {
            let on: Vec<usize> = self.edge_indices(edge);
            let mut strict = true;
            for k in 0..self.points.len() {
                if on.contains(&k) {
                    continue;
                }
                let inside = if high {
                    less(&self.points[k][axis], value, dec)?
                } else {
                    less(value, &self.points[k][axis], dec)?
                };
                if !inside {
                    strict = false;
                    break;
                }
            }
            if strict {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// The second patch's children of this cell, at two and four fifths of
    /// its longer side: computed once, by the request's own thread, and
    /// shared by every first-patch cell whose cover holds this cell. Only the
    /// computing call is charged; the writer fills them in a fixed order, so
    /// the work count does not depend on threads.
    pub(super) fn second_children(
        self: &Arc<Self>,
        dec: &mut Decider,
    ) -> Result<&[Arc<Cell>], Refusal> {
        if let Some(children) = self.children.get() {
            return Ok(children);
        }
        writer_only()?;
        let axis = self.split_axis(dec)?;
        let fifths = [
            Rational::from_ratio(2, 5).expect("nonzero"),
            Rational::from_ratio(4, 5).expect("nonzero"),
        ];
        let children = self
            .split(axis, &fifths, dec)?
            .into_iter()
            .map(Arc::new)
            .collect();
        Ok(self.children.get_or_init(|| children))
    }

    /// Longer local axis; ties take `u`.
    pub(super) fn split_axis(&self, dec: &mut Decider) -> Result<usize, Refusal> {
        let u = sub(&self.local[0][1], &self.local[0][0], dec)?;
        let v = sub(&self.local[1][1], &self.local[1][0], dec)?;
        Ok(usize::from(less(&u, &v, dec)?))
    }

    /// Children at the given fractions of the cell along `axis`, in order.
    /// The fractions are increasing and strictly inside `(0, 1)`. Only the
    /// control nets split here; tangent and normal nets follow on demand.
    pub(super) fn split(
        self: &Arc<Self>,
        axis: usize,
        fractions: &[Rational],
        dec: &mut Decider,
    ) -> Result<Vec<Self>, Refusal> {
        let mut children = Vec::with_capacity(fractions.len() + 1);
        let mut rest = self.h.clone();
        let mut done = Rational::zero();
        let mut cuts = Vec::with_capacity(fractions.len() + 2);
        cuts.push(Rational::zero());
        cuts.extend(fractions.iter().cloned());
        cuts.push(Rational::one());
        for (k, fraction) in fractions.iter().enumerate() {
            // The rest covers [done, 1] of the cell; cut at `fraction`.
            let width = sub(&Rational::one(), &done, dec)?;
            let offset = sub(fraction, &done, dec)?;
            dec.reserve_exact_arithmetic(&offset, &width, true)?;
            let t = &offset / &width;
            let h = take_low(&mut rest, axis, &t, dec)?;
            children.push(self.child(axis, [&cuts[k], &cuts[k + 1]], h, dec)?);
            done = fraction.clone();
        }
        let last = fractions.len();
        children.push(self.child(axis, [&cuts[last], &cuts[last + 1]], rest.to_vec(), dec)?);
        Ok(children)
    }

    fn child(
        self: &Arc<Self>,
        axis: usize,
        part: [&Rational; 2],
        h: Vec<Grid>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut rect = self.rect.clone();
        let mut local = self.local.clone();
        let global_width = sub(&self.rect[axis][1], &self.rect[axis][0], dec)?;
        let local_width = sub(&self.local[axis][1], &self.local[axis][0], dec)?;
        let mut ends = Vec::with_capacity(2);
        let mut local_ends = Vec::with_capacity(2);
        for fraction in part {
            let g = mul(&global_width, fraction, dec)?;
            ends.push(add(&self.rect[axis][0], &g, dec)?);
            let l = mul(&local_width, fraction, dec)?;
            local_ends.push(add(&self.local[axis][0], &l, dec)?);
        }
        rect[axis] = ends.try_into().expect("two ends");
        local[axis] = local_ends.try_into().expect("two ends");
        Self::finish(
            self.span,
            rect,
            local,
            self.depth + 1,
            self.unit_weight,
            h.try_into().expect("four"),
            Some(self.clone()),
            dec,
        )
    }
}

fn cross(
    a: &[Rational; 3],
    b: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    let mut out: [Rational; 3] = std::array::from_fn(|_| Rational::zero());
    for (i, value) in out.iter_mut().enumerate() {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        let left = mul(&a[j], &b[k], dec)?;
        let right = mul(&a[k], &b[j], dec)?;
        *value = sub(&left, &right, dec)?;
    }
    Ok(out)
}

/// Indices of the controls on one border: low u, high u, low v, high v.
fn edge_indices([nu, nv]: [usize; 2], edge: usize) -> Vec<usize> {
    let stride = nv + 1;
    match edge {
        0 => (0..=nv).collect(),
        1 => (nu * stride..=nu * stride + nv).collect(),
        2 => (0..=nu).map(|i| i * stride).collect(),
        _ => (0..=nu).map(|i| i * stride + nv).collect(),
    }
}

/// Bits of the dyadic grid that outward boxes round large controls to.
const COARSE_BITS: i64 = 62;

/// Limbs of a fraction's numerator and denominator, as exact arithmetic
/// charges them.
pub(super) fn limbs(value: &Rational) -> u64 {
    value.numerator().bit_length().div_ceil(64).max(1)
        + value.denominator().bit_length().div_ceil(64).max(1)
}

/// `value` itself while it is small; otherwise the closed dyadic interval
/// of step `2^-62` round it, so that deep cells keep small numbers.
fn round_out(value: &Rational, dec: &mut Decider) -> Result<[Rational; 2], Refusal> {
    if limbs(value) <= 4 {
        return Ok([value.clone(), value.clone()]);
    }
    dec.reserve_exact_arithmetic(value, value, false)?;
    let scaled = value.mul_pow2(COARSE_BITS);
    let lo = Rational::from_integer(scaled.floor()).mul_pow2(-COARSE_BITS);
    let hi = Rational::from_integer(scaled.ceil()).mul_pow2(-COARSE_BITS);
    Ok([lo, hi])
}

/// Binary64 copies of a cell's outward boxes, direction and own slab.
#[derive(Debug)]
struct Rough {
    boxes: Vec<[[f64; 2]; 3]>,
    frame: [f64; 3],
    own: [f64; 2],
}

/// An integer direction near `n`, with components below `2^30`. Slabs are
/// sound along any direction; a small one keeps their arithmetic small.
fn small_direction(n: &[Rational; 3], dec: &mut Decider) -> Result<[Rational; 3], Refusal> {
    let mut largest = Rational::zero();
    for value in n {
        let size = value.abs();
        if less(&largest, &size, dec)? {
            largest = size;
        }
    }
    if largest.is_zero() {
        return Ok(std::array::from_fn(|i| {
            Rational::from_i64(i64::from(i == 2))
        }));
    }
    let mut out: [Rational; 3] = std::array::from_fn(|_| Rational::zero());
    for (value, slot) in n.iter().zip(out.iter_mut()) {
        dec.reserve_exact_arithmetic(value, &largest, true)?;
        let scaled = (value / &largest).mul_pow2(30);
        *slot = Rational::from_integer(scaled.floor());
    }
    Ok(out)
}

/// The range of `direction . P` over boxes round the given controls: a
/// slab round the patch, since it is a positive combination of them.
fn slab_of(
    boxes: &[[[Rational; 2]; 3]],
    direction: &[Rational; 3],
    edge: Option<usize>,
    degree: [usize; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let indices: Vec<usize> = match edge {
        Some(edge) => edge_indices(degree, edge),
        None => (0..boxes.len()).collect(),
    };
    let mut out: Option<[Rational; 2]> = None;
    for k in indices {
        let mut total = [Rational::zero(), Rational::zero()];
        for axis in 0..3 {
            let d = &direction[axis];
            if d.is_zero() {
                continue;
            }
            let [lo, hi] = &boxes[k][axis];
            let (low, high) = if d.sign() == Sign::Positive {
                (lo, hi)
            } else {
                (hi, lo)
            };
            let low = mul(low, d, dec)?;
            let high = mul(high, d, dec)?;
            total = [add(&total[0], &low, dec)?, add(&total[1], &high, dec)?];
        }
        out = Some(match out {
            None => total,
            Some([lo, hi]) => [
                if less(&total[0], &lo, dec)? {
                    total[0].clone()
                } else {
                    lo
                },
                if less(&hi, &total[1], dec)? {
                    total[1].clone()
                } else {
                    hi
                },
            ],
        });
    }
    Ok(out.expect("nonempty controls"))
}

/// Whether two closed boxes are disjoint, with charged comparisons.
pub(super) fn disjoint(a: &Bounds3, b: &Bounds3, dec: &mut Decider) -> Result<bool, Refusal> {
    for axis in 0..3 {
        if less(&a[axis][1], &b[axis][0], dec)? || less(&b[axis][1], &a[axis][0], dec)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Exact interval product and sum, for the certificate bounds.
pub(super) fn interval_mul(
    a: &[Rational; 2],
    b: &[Rational; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let mut products = Vec::with_capacity(4);
    for x in a {
        for y in b {
            products.push(mul(x, y, dec)?);
        }
    }
    values_range(products.iter(), dec)
}

pub(super) fn interval_add(
    a: &[Rational; 2],
    b: &[Rational; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    Ok([add(&a[0], &b[0], dec)?, add(&a[1], &b[1], dec)?])
}

/// An interval round the dot product of two boxed vectors.
pub(super) fn interval_dot(
    a: &Bounds3,
    b: &Bounds3,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let mut total = interval_mul(&a[0], &b[0], dec)?;
    for axis in 1..3 {
        let term = interval_mul(&a[axis], &b[axis], dec)?;
        total = interval_add(&total, &term, dec)?;
    }
    Ok(total)
}

/// An interval round one component of the cross product of two boxed vectors.
pub(super) fn interval_cross(
    a: &Bounds3,
    b: &Bounds3,
    component: usize,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let (j, k) = ((component + 1) % 3, (component + 2) % 3);
    let left = interval_mul(&a[j], &b[k], dec)?;
    let right = interval_mul(&a[k], &b[j], dec)?;
    Ok([
        sub(&left[0], &right[1], dec)?,
        sub(&left[1], &right[0], dec)?,
    ])
}
