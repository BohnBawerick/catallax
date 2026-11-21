//! Certified stationary pairs of two strata whose cores have no closed form.
//!
//! Circles with non-parallel normals, a line in general position against a
//! circle, and any line, circle or cone against a cone have no closed form in
//! [`super::super::critical`]. OpenCASCADE answers them numerically:
//! `Extrema_GenExtCC` sweeps two curves with `math_GlobOptMin`, `Extrema_GenExtCS`
//! seeds `math_PSO` along a curve projected onto a cone or torus, and
//! `Extrema_GenExtSS` takes the best pair of two 20 x 20 grids; each refines a
//! candidate with Newton or BFGS steps. SMLib's `SmGlobalSolver` subdivides
//! both shapes, prunes a branch when its boxes are farther apart than the best
//! answer less a tolerance, and runs Newton at the leaves. Both then keep a
//! solution only when it lies on the trimmed face (`BRepClass_FaceClassifier`
//! in `BRepExtrema_DistanceSS`) and leave boundary minima to the edge and
//! vertex pairs, as the strata search does. Neither certifies a solution.
//!
//! This module follows SMLib's subdivision and pruning and OpenCASCADE's split
//! into strata pairs and trim classification, and proves each step:
//!
//! - Each core is parametrized: a line by arc length, a circle by the rational
//!   half-angle parameter of each quarter turn, a cone by that parameter and the
//!   height along its axis. Ball arithmetic encloses the core points, the
//!   stationary equations `F = 0` (half the gradient of the squared core
//!   distance, the cone's turning equation divided by its height so the apex is
//!   no singularity) and their Jacobian over a box of parameters.
//! - Boxes are visited nearest first. A box is dropped when its carrier points
//!   cannot lie in the strata's boxes, when no offset pair from it can come
//!   below the best distance found, or when some component of `F` keeps one
//!   sign on it, the monotonicity test of interval global optimization (Hansen
//!   and Walster, Global Optimization Using Interval Analysis, 2004, ch. 12).
//! - The Krawczyk operator (Krawczyk 1969) contracts every other box, and drops
//!   it when the image misses it. An image inside the box's interior proves
//!   exactly one solution there (Moore 1977). A box the operator shrank by half
//!   without proving one is retried at its own size around the shrunken box,
//!   which catches a solution on the box's face, as Rump's epsilon-inflation
//!   does (Rump 1983).
//! - The operator refines each proven solution. The mean value form encloses
//!   its offset distance. An exact pair of carrier points at the refined box's
//!   centre is classified on both trims: admitted, it gives the upper bound.
//!
//! A pair not admitted at the centre may still be admitted at the solution
//! when a trim boundary passes through the refined box. Then a point of the
//! stratum's boundary lies within the box's diagonal of the solution, on a
//! lower stratum the search also visits, which attains the solution's value
//! up to that diagonal. So the solution is dropped and its diagonals are added
//! to a slack subtracted from the final lower bound. Each drop reaches a lower
//! stratum, so the chain of drops behind any value is no longer than the
//! number of drops, and their summed diagonals bound its loss.

use super::super::critical::{Core, perpendicular};
use super::super::{root, unit};
use super::{Kind, Member, Scratch, Stratum, separation};
use crate::Accuracy;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::{Point3, Vector3};
use std::{cmp::Ordering, collections::BinaryHeap, rc::Rc};

/// One limb, the cheapest ball arithmetic. A proven solution moves to the
/// precisions `enclose` uses only when its value misses the requested accuracy.
const SEARCH_BITS: u32 = 64;
const LADDER: [u32; 2] = [128, 256];

type V3 = [Ball; 3];

/// Ball arithmetic charged by the decider's rule for logical work: a sum costs
/// `4L` units, a product, quotient or root `8L^2`, for `L` 64-bit limbs.
struct Meter {
    bits: u32,
    add: u64,
    mul: u64,
    units: u64,
}

impl Meter {
    fn new(bits: u32) -> Self {
        let l = u64::from(bits.div_ceil(64).max(1));
        Self {
            bits,
            add: 4 * l,
            mul: 8 * l * l,
            units: 0,
        }
    }
    fn ball(&self, r: &Rational) -> Ball {
        Ball::from_rational(r, self.bits)
    }
    fn span(&self, lo: &Rational, hi: &Rational) -> Ball {
        self.ball(lo).hull(&self.ball(hi))
    }
    fn int(&self, n: i64) -> Ball {
        self.ball(&Rational::from_i64(n))
    }
    fn zero3(&self) -> V3 {
        std::array::from_fn(|_| self.int(0))
    }
    fn add(&mut self, a: &Ball, b: &Ball) -> Ball {
        self.units += self.add;
        a + b
    }
    fn sub(&mut self, a: &Ball, b: &Ball) -> Ball {
        self.units += self.add;
        a - b
    }
    fn mul(&mut self, a: &Ball, b: &Ball) -> Ball {
        self.units += self.mul;
        a * b
    }
    fn div(&mut self, a: &Ball, b: &Ball) -> Option<Ball> {
        self.units += self.mul;
        a.checked_div(b)
    }
    fn sqrt(&mut self, a: &Ball) -> Option<Ball> {
        self.units += self.mul;
        a.sqrt()
    }
    fn square(&mut self, a: &Ball) -> Ball {
        self.units += self.mul;
        a.square()
    }
    fn plus(&mut self, a: &V3, b: &V3) -> V3 {
        std::array::from_fn(|i| self.add(&a[i], &b[i]))
    }
    fn minus(&mut self, a: &V3, b: &V3) -> V3 {
        std::array::from_fn(|i| self.sub(&a[i], &b[i]))
    }
    fn times(&mut self, s: &Ball, a: &V3) -> V3 {
        std::array::from_fn(|i| self.mul(s, &a[i]))
    }
    /// `a u + b v`.
    fn comb(&mut self, a: &Ball, u: &V3, b: &Ball, v: &V3) -> V3 {
        let (x, y) = (self.times(a, u), self.times(b, v));
        self.plus(&x, &y)
    }
    fn dot(&mut self, a: &V3, b: &V3) -> Ball {
        let mut sum = self.mul(&a[0], &b[0]);
        for i in 1..3 {
            let t = self.mul(&a[i], &b[i]);
            sum = self.add(&sum, &t);
        }
        sum
    }
    fn cross(&mut self, a: &V3, b: &V3) -> V3 {
        std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let (x, y) = (self.mul(&a[j], &b[k]), self.mul(&a[k], &b[j]));
            self.sub(&x, &y)
        })
    }
    fn norm2(&mut self, a: &V3) -> Ball {
        let mut sum = self.square(&a[0]);
        for x in &a[1..] {
            let t = self.square(x);
            sum = self.add(&sum, &t);
        }
        sum
    }
    fn charge(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        dec.charge_work(std::mem::take(&mut self.units))
    }
}

fn neg3(v: &V3) -> V3 {
    std::array::from_fn(|i| -&v[i])
}

/// The least `|z|` over a ball.
fn abs_lower(z: &Ball) -> Rational {
    z.lo().max(-z.hi()).max(Rational::zero())
}

fn abs_ball(z: &Ball, m: &Meter) -> Ball {
    let top = z.hi().max(-z.lo());
    m.span(&abs_lower(z), &top)
}

fn meet(a: &Ball, b: &Ball, m: &Meter) -> Option<Ball> {
    let (lo, hi) = (a.lo().max(b.lo()), a.hi().min(b.hi()));
    (lo <= hi).then(|| m.span(&lo, &hi))
}

fn overlaps(a: &V3, b: &V3) -> bool {
    (0..3).all(|i| a[i].lo() <= b[i].hi() && b[i].lo() <= a[i].hi())
}

/// An upper bound on the distance between two points of a box.
fn diagonal(a: &V3, m: &mut Meter) -> Rational {
    let mut sum = m.int(0);
    for x in a {
        let w = m.ball(&(x.hi() - x.lo()));
        let w = m.square(&w);
        sum = m.add(&sum, &w);
    }
    m.sqrt(&sum).map_or_else(Rational::zero, |r| r.hi())
}

type Cell = Vec<[Rational; 2]>;

fn size(cell: &Cell) -> Rational {
    cell.iter()
        .fold(Rational::zero(), |s, [lo, hi]| s + (hi - lo))
}

fn centre(cell: &Cell) -> Vec<Rational> {
    cell.iter().map(|[lo, hi]| (lo + hi).mul_pow2(-1)).collect()
}

fn inside(outer: &Cell, inner: &Cell) -> bool {
    outer
        .iter()
        .zip(inner)
        .all(|([a, b], [c, d])| a <= c && d <= b)
}

fn rational_turn(t: &Rational) -> (Rational, Rational) {
    let w = Rational::one() + t * t;
    (
        &(Rational::one() - t * t) / &w,
        &(t * &Rational::from_i64(2)) / &w,
    )
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Shape {
    Line,
    Circle,
    Cone,
}

/// A core's exact definition: `origin` is a line's point, a circle's centre or
/// a cone's apex; `axis` the line's direction, the circle's normal or the
/// cone's axis; `e1, e2` an orthonormal frame of the circle's plane, or
/// perpendicular to the cone's axis; `size` the radius or the slope; `offset`
/// the carrier's distance from the core.
struct Exact {
    origin: Point3,
    axis: Vector3,
    e1: Vector3,
    e2: Vector3,
    size: Real,
    offset: Option<Real>,
}

/// Ball enclosures of an [`Exact`] at one precision.
#[derive(Clone)]
struct Data {
    origin: V3,
    axis: V3,
    e1: V3,
    e2: V3,
    size: Ball,
    offset: Option<Ball>,
}

fn enclose3(x: [&Real; 3], bits: u32, dec: &mut Decider) -> Result<V3, Refusal> {
    Ok([
        dec.enclose_budgeted_at(x[0], bits)?,
        dec.enclose_budgeted_at(x[1], bits)?,
        dec.enclose_budgeted_at(x[2], bits)?,
    ])
}

impl Exact {
    fn enclose(&self, bits: u32, dec: &mut Decider) -> Result<Data, Refusal> {
        fn v(v: &Vector3) -> [&Real; 3] {
            [&v.x, &v.y, &v.z]
        }
        let p = &self.origin;
        Ok(Data {
            origin: enclose3([&p.x, &p.y, &p.z], bits, dec)?,
            axis: enclose3(v(&self.axis), bits, dec)?,
            e1: enclose3(v(&self.e1), bits, dec)?,
            e2: enclose3(v(&self.e2), bits, dec)?,
            size: dec.enclose_budgeted_at(&self.size, bits)?,
            offset: match &self.offset {
                Some(r) => Some(dec.enclose_budgeted_at(r, bits)?),
                None => None,
            },
        })
    }

    /// The frame of quarter turn `q`.
    fn frame(&self, q: usize) -> (Vector3, Vector3) {
        let neg = |v: &Vector3| v.scaled(&Real::from_i64(-1));
        match q {
            0 => (self.e1.clone(), self.e2.clone()),
            1 => (self.e2.clone(), neg(&self.e1)),
            2 => (neg(&self.e1), neg(&self.e2)),
            _ => (neg(&self.e2), self.e1.clone()),
        }
    }
}

fn frame(d: &Data, q: usize) -> (V3, V3) {
    match q {
        0 => (d.e1.clone(), d.e2.clone()),
        1 => (d.e2.clone(), neg3(&d.e1)),
        2 => (neg3(&d.e1), neg3(&d.e2)),
        _ => (neg3(&d.e2), d.e1.clone()),
    }
}

/// The quarter-turn chart `t -> ((1 - t^2), 2t) / (1 + t^2)` and its first two
/// derivatives, as `[[c, s], [c', s'], [c'', s'']]`.
fn turn(t: &Ball, m: &mut Meter) -> [[Ball; 2]; 3] {
    let one = m.int(1);
    let t2 = m.square(t);
    let w = m.add(&one, &t2);
    let w2 = m.square(&w);
    let w3 = m.mul(&w2, &w);
    let u = m.sub(&one, &t2);
    let two_t = m.mul(&m.int(2), t);
    let q = |a: &Ball, b: &Ball, m: &mut Meter| m.div(a, b).expect("1 + t^2 >= 1");
    let c = q(&u, &w, m);
    let s = q(&two_t, &w, m);
    let n = m.mul(&m.int(-4), t);
    let c1 = q(&n, &w2, m);
    let n = m.mul(&m.int(2), &u);
    let s1 = q(&n, &w2, m);
    let n = m.mul(&m.int(12), &t2);
    let n = m.sub(&n, &m.int(4));
    let c2 = q(&n, &w3, m);
    let n = m.sub(&t2, &m.int(3));
    let n = m.mul(&n, t);
    let n = m.mul(&m.int(4), &n);
    let s2 = q(&n, &w3, m);
    [[c, s], [c1, s1], [c2, s2]]
}

/// The normal space of a core at a point, as in `critical::meet`.
enum Normal {
    /// Every vector perpendicular to this one.
    Perp(V3),
    /// Multiples of this vector.
    Span(V3),
}

/// A core's point, its partial derivatives, the tangents `g` of its stationary
/// equations (the partials, with the cone's turning partial divided by its
/// height) and their partials `dg[i][k]`, over a box of parameters.
struct Jet {
    p: V3,
    dp: Vec<V3>,
    g: Vec<V3>,
    dg: Vec<Vec<V3>>,
    normal: Normal,
    /// The unit vector whose perpendicular plane holds the carrier offsets:
    /// a line's direction or a circle's tangent.
    along: V3,
}

struct Patch {
    shape: Shape,
    exact: Exact,
    data: Data,
    /// Range of the linear parameter: a line's arc length, a cone's height.
    range: Option<[Rational; 2]>,
}

impl Patch {
    fn new(s: &Stratum, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        let zero = Vector3::zero;
        let (shape, origin, axis, e1, e2, size) = match &s.core {
            Some(Core::Line { origin, dir }) => (
                Shape::Line,
                origin.clone(),
                dir.clone(),
                zero(),
                zero(),
                Real::zero(),
            ),
            Some(Core::Circle {
                centre,
                normal,
                radius,
                reference,
            }) => (
                Shape::Circle,
                centre.clone(),
                normal.clone(),
                reference.clone(),
                normal.cross(reference),
                radius.clone(),
            ),
            Some(Core::Cone { apex, axis, slope }) => {
                let e1 = perpendicular(axis, dec)?;
                let e2 = axis.cross(&e1);
                (
                    Shape::Cone,
                    apex.clone(),
                    axis.clone(),
                    e1,
                    e2,
                    slope.clone(),
                )
            }
            _ => return Ok(None),
        };
        let exact = Exact {
            origin,
            axis,
            e1,
            e2,
            size,
            offset: s.offset.clone(),
        };
        let data = exact.enclose(SEARCH_BITS, dec)?;
        let range = if shape == Shape::Circle {
            None
        } else {
            // A carrier point's offset is perpendicular to the axis, so its
            // projection onto the axis is the parameter.
            let Some(extent) = &s.extent else {
                return Ok(None);
            };
            let mut m = Meter::new(SEARCH_BITS);
            let rel = m.minus(extent, &data.origin);
            let t = m.dot(&rel, &data.axis);
            m.charge(dec)?;
            let lo = if shape == Shape::Cone {
                t.lo().max(Rational::zero())
            } else {
                t.lo()
            };
            if lo > t.hi() {
                return Ok(None);
            }
            Some([lo, t.hi()])
        };
        Ok(Some(Self {
            shape,
            exact,
            data,
            range,
        }))
    }

    fn dims(&self) -> usize {
        if self.shape == Shape::Cone { 2 } else { 1 }
    }

    fn charts(&self) -> usize {
        if self.shape == Shape::Line { 1 } else { 4 }
    }

    fn start(&self) -> Cell {
        let turn = [Rational::zero(), Rational::one()];
        match self.shape {
            Shape::Line => vec![self.range.clone().expect("a line has a range")],
            Shape::Circle => vec![turn],
            Shape::Cone => vec![turn, self.range.clone().expect("a cone has a range")],
        }
    }

    fn jet(&self, d: &Data, q: usize, x: &[Ball], m: &mut Meter) -> Jet {
        match self.shape {
            Shape::Line => {
                let step = m.times(&x[0], &d.axis);
                Jet {
                    p: m.plus(&d.origin, &step),
                    dp: vec![d.axis.clone()],
                    g: vec![d.axis.clone()],
                    dg: vec![vec![m.zero3()]],
                    normal: Normal::Perp(d.axis.clone()),
                    along: d.axis.clone(),
                }
            }
            Shape::Circle => {
                let (u, v) = frame(d, q);
                let [[c, s], [c1, s1], [c2, s2]] = turn(&x[0], m);
                let e = m.comb(&c, &u, &s, &v);
                let e1 = m.comb(&c1, &u, &s1, &v);
                let e2 = m.comb(&c2, &u, &s2, &v);
                let re = m.times(&d.size, &e);
                let tangent = m.times(&d.size, &e1);
                let bend = m.times(&d.size, &e2);
                let ns = -&s;
                Jet {
                    p: m.plus(&d.origin, &re),
                    dp: vec![tangent.clone()],
                    g: vec![tangent.clone()],
                    dg: vec![vec![bend]],
                    normal: Normal::Perp(tangent),
                    along: m.comb(&ns, &u, &c, &v),
                }
            }
            Shape::Cone => {
                let (u, v) = frame(d, q);
                let [[c, s], [c1, s1], [c2, s2]] = turn(&x[0], m);
                let e = m.comb(&c, &u, &s, &v);
                let e1 = m.comb(&c1, &u, &s1, &v);
                let e2 = m.comb(&c2, &u, &s2, &v);
                let se = m.times(&d.size, &e);
                let generator = m.plus(&d.axis, &se);
                let lift = m.times(&x[1], &generator);
                let swing = m.times(&d.size, &e1);
                let bend = m.times(&d.size, &e2);
                let sa = m.times(&d.size, &d.axis);
                Jet {
                    p: m.plus(&d.origin, &lift),
                    dp: vec![m.times(&x[1], &swing), generator.clone()],
                    g: vec![swing.clone(), generator],
                    dg: vec![vec![bend, m.zero3()], vec![swing, m.zero3()]],
                    normal: Normal::Span(m.minus(&e, &sa)),
                    along: m.zero3(),
                }
            }
        }
    }

    /// The exact core point at `at` in chart `q`, and the unit vector whose
    /// perpendicular plane holds the carrier offsets there.
    fn point(&self, q: usize, at: &[Rational]) -> (Point3, Vector3) {
        let e = &self.exact;
        match self.shape {
            Shape::Line => (
                e.origin.plus(&e.axis.scaled(&Real::from(at[0].clone()))),
                e.axis.clone(),
            ),
            Shape::Circle | Shape::Cone => {
                let (u, v) = e.frame(q);
                let (c, s) = rational_turn(&at[0]);
                let (c, s) = (Real::from(c), Real::from(s));
                let dir = u.scaled(&c).plus(&v.scaled(&s));
                if self.shape == Shape::Circle {
                    let tangent = u.scaled(&-&s).plus(&v.scaled(&c));
                    (e.origin.plus(&dir.scaled(&e.size)), tangent)
                } else {
                    let generator = e.axis.plus(&dir.scaled(&e.size));
                    (
                        e.origin.plus(&generator.scaled(&Real::from(at[1].clone()))),
                        Vector3::zero(),
                    )
                }
            }
        }
    }
}

/// The part of `w` perpendicular to the unit vector `t`.
fn normal_part(w: &V3, t: &V3, m: &mut Meter) -> V3 {
    let k = m.dot(w, t);
    let along = m.times(&k, t);
    m.minus(w, &along)
}

/// Every pair of points the solution's offsets can reach must lie in the
/// stratum's box; the box of the core points inflated by the offset holds them.
fn cover(p: &V3, offset: &Option<Ball>, extent: Option<&V3>, m: &mut Meter) -> Option<V3> {
    let mut out = p.clone();
    if let Some(r) = offset {
        let reach = m.span(&-r.hi(), &r.hi());
        out = std::array::from_fn(|i| m.add(&p[i], &reach));
    }
    match extent {
        None => Some(out),
        Some(e) => {
            let mut clipped = Vec::with_capacity(3);
            for i in 0..3 {
                clipped.push(meet(&out[i], &e[i], m)?);
            }
            Some(std::array::from_fn(|i| clipped[i].clone()))
        }
    }
}

fn signs(offset: &Option<Ball>) -> &'static [i64] {
    if offset.is_some() { &[1, -1] } else { &[0] }
}

/// The offset combination `sx rx + sy ry` of each pair of carrier sides.
fn combos(dx: &Data, dy: &Data, m: &mut Meter) -> Vec<(i64, i64, Ball)> {
    let mut out = Vec::new();
    for &sx in signs(&dx.offset) {
        for &sy in signs(&dy.offset) {
            let part = |o: &Option<Ball>, s: i64, m: &mut Meter| match o {
                Some(r) => m.mul(&m.int(s), r),
                None => m.int(0),
            };
            let (a, b) = (part(&dx.offset, sx, m), part(&dy.offset, sy, m));
            out.push((sx, sy, m.add(&a, &b)));
        }
    }
    out
}

/// Half the gradient of `|P(x) - Q(y)|^2`, with the cone's turning component
/// divided by its height, and its Jacobian when asked.
fn equations(jx: &Jet, jy: &Jet, jacobian: bool, m: &mut Meter) -> (Vec<Ball>, Vec<Vec<Ball>>) {
    let d = m.minus(&jx.p, &jy.p);
    let (kx, ky) = (jx.g.len(), jy.g.len());
    let mut f = Vec::with_capacity(kx + ky);
    for g in &jx.g {
        f.push(m.dot(&d, g));
    }
    for h in &jy.g {
        f.push(-&m.dot(&d, h));
    }
    if !jacobian {
        return (f, Vec::new());
    }
    let n = kx + ky;
    let mut jac = vec![vec![m.int(0); n]; n];
    for i in 0..kx {
        for k in 0..kx {
            let a = m.dot(&jx.dp[k], &jx.g[i]);
            let b = m.dot(&d, &jx.dg[i][k]);
            jac[i][k] = m.add(&a, &b);
        }
        for l in 0..ky {
            jac[i][kx + l] = -&m.dot(&jy.dp[l], &jx.g[i]);
        }
    }
    for j in 0..ky {
        for k in 0..kx {
            jac[kx + j][k] = -&m.dot(&jx.dp[k], &jy.g[j]);
        }
        for l in 0..ky {
            let a = m.dot(&jy.dp[l], &jy.g[j]);
            let b = m.dot(&d, &jy.dg[j][l]);
            jac[kx + j][kx + l] = m.sub(&a, &b);
        }
    }
    (f, jac)
}

/// A point approximate inverse of the midpoint matrix, by Gauss-Jordan
/// elimination with partial pivoting, or `None` when a pivot cannot be told
/// from zero. Any point matrix keeps the Krawczyk operator valid; this one
/// makes it contract.
fn preconditioner(a: &[Vec<Ball>], m: &mut Meter) -> Option<Vec<Vec<Ball>>> {
    let n = a.len();
    let mut rows: Vec<Vec<Ball>> = (0..n)
        .map(|i| {
            let mut row: Vec<Ball> = a[i].iter().map(|b| m.ball(&b.mid())).collect();
            row.extend((0..n).map(|j| m.int(i64::from(i == j))));
            row
        })
        .collect();
    for c in 0..n {
        let p = (c..n).max_by(|&x, &y| {
            rows[x][c]
                .mid()
                .abs()
                .cmp(&rows[y][c].mid().abs())
                .then(y.cmp(&x))
        })?;
        rows.swap(c, p);
        let pivot = rows[c][c].clone();
        if pivot.contains_zero() {
            return None;
        }
        for j in 0..2 * n {
            rows[c][j] = m.div(&rows[c][j], &pivot)?;
        }
        for r in 0..n {
            if r == c {
                continue;
            }
            let f = rows[r][c].clone();
            for j in 0..2 * n {
                let t = m.mul(&f, &rows[c][j]);
                rows[r][j] = m.sub(&rows[r][j], &t);
            }
        }
    }
    Some(
        rows.into_iter()
            .map(|row| row[n..].iter().map(|b| m.ball(&b.mid())).collect())
            .collect(),
    )
}

enum Step {
    /// The image misses the box: no solution in it.
    Empty,
    /// The image lies in the box's interior: exactly one solution in it.
    Unique,
    /// Every solution in the box lies in this smaller box.
    Contracted(Cell),
}

fn classify(k: &[Ball], cell: &Cell) -> Step {
    let mut unique = true;
    let mut out = Vec::with_capacity(cell.len());
    for (ki, [lo, hi]) in k.iter().zip(cell) {
        let (klo, khi) = (ki.lo(), ki.hi());
        if khi < *lo || klo > *hi {
            return Step::Empty;
        }
        unique &= klo > *lo && khi < *hi;
        out.push([klo.max(lo.clone()), khi.min(hi.clone())]);
    }
    if unique {
        Step::Unique
    } else {
        Step::Contracted(out)
    }
}

/// Two strata with parametrized cores.
struct Pair<'s> {
    x: Rc<Patch>,
    y: Rc<Patch>,
    extents: [Option<&'s V3>; 2],
    members: [&'s Member; 2],
}

/// What a parameter box gives before its Krawczyk step.
struct Assessed {
    bound: Rational,
    weights: Vec<Rational>,
    jac: Vec<Vec<Ball>>,
}

/// One side combination at a solution: its signs, its offset distance and
/// the boxes of its two carrier points, when both offset planes are proven.
type Sided = (i64, i64, Ball, Option<(V3, V3)>);

/// The offset distance and carrier boxes at a refined solution.
struct Solution {
    dist: Ball,
    heights: Vec<Ball>,
    normals: (Normal, Normal),
    combos: Vec<Sided>,
}

impl Pair<'_> {
    fn split(&self) -> usize {
        self.x.dims()
    }

    fn jets(
        &self,
        charts: [usize; 2],
        at: &[Ball],
        d: (&Data, &Data),
        m: &mut Meter,
    ) -> (Jet, Jet) {
        let k = self.split();
        (
            self.x.jet(d.0, charts[0], &at[..k], m),
            self.y.jet(d.1, charts[1], &at[k..], m),
        )
    }

    fn balls(cell: &Cell, m: &Meter) -> Vec<Ball> {
        cell.iter().map(|[lo, hi]| m.span(lo, hi)).collect()
    }

    /// The lower bound of every admitted offset pair from the box, or `None`
    /// when the box holds no solution worth visiting below `upper`.
    fn assess(
        &self,
        charts: [usize; 2],
        cell: &Cell,
        upper: Option<&Rational>,
        m: &mut Meter,
    ) -> Option<Assessed> {
        let d = (&self.x.data, &self.y.data);
        let at = Self::balls(cell, m);
        let (jx, jy) = self.jets(charts, &at, d, m);
        let cx = cover(&jx.p, &d.0.offset, self.extents[0], m)?;
        let cy = cover(&jy.p, &d.1.offset, self.extents[1], m)?;
        let mut bound = separation(&cx, &cy);
        let gap = m.minus(&jx.p, &jy.p);
        let g = m.norm2(&gap);
        let dist = m.sqrt(&g)?;
        let least = combos(d.0, d.1, m)
            .iter()
            .map(|(_, _, c)| abs_lower(&m.sub(&dist, c)))
            .min()
            .expect("at least one combination");
        bound = bound.max(least);
        if upper.is_some_and(|u| bound >= *u) {
            return None;
        }
        let (f, jac) = equations(&jx, &jy, true, m);
        if f.iter().any(|v| !v.contains_zero()) {
            return None;
        }
        let weights = jx
            .dp
            .iter()
            .chain(&jy.dp)
            .map(|v| {
                v.iter()
                    .map(|c| c.hi().max(-c.lo()))
                    .max()
                    .expect("three components")
            })
            .collect();
        Some(Assessed {
            bound,
            weights,
            jac,
        })
    }

    /// The Krawczyk image `m - Y F(m) + (I - Y J(X)) (X - m)` of a box.
    fn krawczyk(
        &self,
        charts: [usize; 2],
        cell: &Cell,
        jac: Option<&Vec<Vec<Ball>>>,
        d: (&Data, &Data),
        m: &mut Meter,
    ) -> Option<Vec<Ball>> {
        let at = Self::balls(cell, m);
        let owned;
        let jac = match jac {
            Some(j) => j,
            None => {
                let (jx, jy) = self.jets(charts, &at, d, m);
                owned = equations(&jx, &jy, true, m).1;
                &owned
            }
        };
        let mid: Vec<Ball> = centre(cell).iter().map(|c| m.ball(c)).collect();
        let (jx, jy) = self.jets(charts, &mid, d, m);
        let (fm, _) = equations(&jx, &jy, false, m);
        let y = preconditioner(jac, m)?;
        let n = cell.len();
        let delta: Vec<Ball> = (0..n).map(|k| m.sub(&at[k], &mid[k])).collect();
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let mut ki = mid[i].clone();
            for k in 0..n {
                let t = m.mul(&y[i][k], &fm[k]);
                ki = m.sub(&ki, &t);
            }
            for k in 0..n {
                let mut entry = m.int(i64::from(i == k));
                for l in 0..n {
                    let t = m.mul(&y[i][l], &jac[l][k]);
                    entry = m.sub(&entry, &t);
                }
                let t = m.mul(&entry, &delta[k]);
                ki = m.add(&ki, &t);
            }
            out.push(ki);
        }
        Some(out)
    }

    /// Shrinks a box holding exactly one solution with the same operator until
    /// a step no longer halves it, or every side is as narrow as the working
    /// precision resolves on its chart's range. Each step halves it, so the
    /// precision's bit count bounds the steps.
    fn refine(&self, charts: [usize; 2], mut cell: Cell, d: (&Data, &Data), m: &mut Meter) -> Cell {
        let mut start = self.x.start();
        start.extend(self.y.start());
        let floor: Vec<Rational> = start
            .iter()
            .map(|[lo, hi]| (hi - lo).mul_pow2(-i64::from(m.bits)))
            .collect();
        for _ in 0..m.bits {
            if cell.iter().zip(&floor).all(|([lo, hi], f)| &(hi - lo) <= f) {
                break;
            }
            let Some(k) = self.krawczyk(charts, &cell, None, d, m) else {
                break;
            };
            let next = match classify(&k, &cell) {
                Step::Empty => break,
                Step::Unique => k.iter().map(|b| [b.lo(), b.hi()]).collect(),
                Step::Contracted(c) => c,
            };
            let halved = size(&next).mul_pow2(1) < size(&cell);
            cell = next;
            if !halved {
                break;
            }
        }
        cell
    }

    /// The offset distance at a refined solution, by the mean value form, and
    /// each side combination's value and carrier boxes.
    fn solution(
        &self,
        charts: [usize; 2],
        cell: &Cell,
        d: (&Data, &Data),
        m: &mut Meter,
    ) -> Solution {
        let at = Self::balls(cell, m);
        let mid: Vec<Ball> = centre(cell).iter().map(|c| m.ball(c)).collect();
        let (jx, jy) = self.jets(charts, &at, d, m);
        let (mx, my) = self.jets(charts, &mid, d, m);
        let gap = m.minus(&jx.p, &jy.p);
        let natural = m.norm2(&gap);
        let gm = m.minus(&mx.p, &my.p);
        let mut form = m.norm2(&gm);
        let k = self.split();
        for (i, dp) in jx.dp.iter().chain(&jy.dp).enumerate() {
            let slope = m.dot(&gap, dp);
            let sign = if i < k { 2 } else { -2 };
            let slope = m.mul(&m.int(sign), &slope);
            let step = m.sub(&at[i], &mid[i]);
            let t = m.mul(&slope, &step);
            form = m.add(&form, &t);
        }
        let g = meet(&natural, &form, m).unwrap_or(natural);
        let g = m.span(&g.lo().max(Rational::zero()), &g.hi().max(Rational::zero()));
        let dist = m.sqrt(&g).expect("nonnegative");
        let mut heights = Vec::new();
        if self.x.shape == Shape::Cone {
            heights.push(at[1].clone());
        }
        if self.y.shape == Shape::Cone {
            heights.push(at[k + 1].clone());
        }
        // The common direction from the first core point to the second,
        // projected onto each side's offset plane.
        let w = m.minus(&jy.p, &jx.p);
        let side = |data: &Data, jet: &Jet, m: &mut Meter| -> Option<Option<V3>> {
            if data.offset.is_none() {
                return Some(None);
            }
            let part = normal_part(&w, &jet.along, m);
            let n2 = m.norm2(&part);
            let n = m.sqrt(&n2)?;
            if n.contains_zero() {
                return None;
            }
            let mut u = Vec::with_capacity(3);
            for x in &part {
                u.push(m.div(x, &n)?);
            }
            Some(Some(std::array::from_fn(|i| u[i].clone())))
        };
        let ux = side(d.0, &jx, m);
        let uy = side(d.1, &jy, m);
        let mut out = Vec::new();
        for (sx, sy, c) in combos(d.0, d.1, m) {
            let z = m.sub(&dist, &c);
            let value = abs_ball(&z, m);
            let carriers =
                match (&ux, &uy) {
                    (Some(ux), Some(uy)) => {
                        let reach =
                            |base: &V3, u: &Option<V3>, r: &Option<Ball>, s: i64, m: &mut Meter| {
                                match (u, r) {
                                    (Some(u), Some(r)) => {
                                        let k = m.mul(&m.int(s), r);
                                        let step = m.times(&k, u);
                                        m.plus(base, &step)
                                    }
                                    _ => base.clone(),
                                }
                            };
                        let pc = reach(&jx.p, ux, &d.0.offset, sx, m);
                        let qc = reach(&jy.p, uy, &d.1.offset, -sy, m);
                        Some((pc, qc))
                    }
                    _ => None,
                };
            out.push((sx, sy, value, carriers));
        }
        Solution {
            dist,
            heights,
            normals: (jx.normal, jy.normal),
            combos: out,
        }
    }

    /// Exact carrier points of combination `(sx, sy)` at `at`: the core points
    /// moved along the part of their difference in each offset plane.
    fn witness(&self, charts: [usize; 2], at: &[Rational], sx: i64, sy: i64) -> (Point3, Point3) {
        let k = self.split();
        let (a, ta) = self.x.point(charts[0], &at[..k]);
        let (b, tb) = self.y.point(charts[1], &at[k..]);
        let w = b.minus(&a);
        let moved = |base: &Point3, t: &Vector3, offset: &Option<Real>, s: i64| match offset {
            Some(r) => {
                let part = w.minus(&t.scaled(&w.dot(t)));
                base.plus(&unit(&part).scaled(&(r * &Real::from_i64(s))))
            }
            None => base.clone(),
        };
        (
            moved(&a, &ta, &self.x.exact.offset, sx),
            moved(&b, &tb, &self.y.exact.offset, -sy),
        )
    }
}

/// Whether no common normal of the two cores can exist on the boxes, as
/// `critical::meet` decides it for exact points.
fn no_common_normal(n: &(Normal, Normal), m: &mut Meter) -> bool {
    match n {
        (Normal::Perp(_), Normal::Perp(_)) => false,
        (Normal::Perp(v), Normal::Span(w)) | (Normal::Span(w), Normal::Perp(v)) => {
            !m.dot(v, w).contains_zero()
        }
        (Normal::Span(v), Normal::Span(w)) => m.cross(v, w).iter().any(|c| !c.contains_zero()),
    }
}

struct Node {
    bound: Rational,
    sequence: usize,
    charts: [usize; 2],
    cell: Cell,
    weights: Vec<Rational>,
    jac: Vec<Vec<Ball>>,
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.bound == other.bound && self.sequence == other.sequence
    }
}
impl Eq for Node {}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .bound
            .cmp(&self.bound)
            .then_with(|| other.sequence.cmp(&self.sequence))
    }
}

/// What one strata pair gives.
pub(super) struct Settled {
    /// A lower bound of every admitted stationary value of the pair below the
    /// bound in force, or `None` when no such value exists.
    pub(super) low: Option<Rational>,
    /// Whether the cell budget ran out before every box was settled.
    pub(super) exhausted: bool,
}

/// The certified search over strata pairs without closed forms. One budget of
/// boxes serves every pair, in the order the caller gives them.
pub(super) struct Pass<'s, 'a> {
    scratch: &'s mut Scratch<'a>,
    /// An upper bound attained before the pass.
    before: Option<Rational>,
    /// The least upper bound attained by an admitted pair of this pass.
    pub(super) witness: Option<Rational>,
    /// The summed diagonals of dropped solutions, see the module notes.
    pub(super) slack: Rational,
    cells: usize,
    /// Each stratum's parametrized core, built once for all its pairs.
    patches: Vec<(Kind, Option<Rc<Patch>>)>,
}

impl<'s, 'a> Pass<'s, 'a> {
    pub(super) fn new(
        scratch: &'s mut Scratch<'a>,
        before: Option<Rational>,
        cells: usize,
    ) -> Self {
        Self {
            scratch,
            before,
            witness: None,
            slack: Rational::zero(),
            cells,
            patches: Vec::new(),
        }
    }

    fn patch(&mut self, s: &Stratum, dec: &mut Decider) -> Result<Option<Rc<Patch>>, Refusal> {
        if let Some((_, p)) = self.patches.iter().find(|(k, _)| *k == s.kind) {
            return Ok(p.clone());
        }
        let p = Patch::new(s, dec)?.map(Rc::new);
        self.patches.push((s.kind, p.clone()));
        Ok(p)
    }

    fn bound(&self) -> Option<Rational> {
        self.before.iter().chain(&self.witness).min().cloned()
    }

    /// Settles one strata pair, or `None` when a core is not a line, circle
    /// or cone, or a range cannot be bounded.
    pub(super) fn settle(
        &mut self,
        x: &Stratum,
        y: &Stratum,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Option<Settled>, Refusal> {
        let (Some(px), Some(py)) = (self.patch(x, dec)?, self.patch(y, dec)?) else {
            return Ok(None);
        };
        let pair = Pair {
            x: px,
            y: py,
            extents: [x.extent.as_ref(), y.extent.as_ref()],
            members: [&x.member, &y.member],
        };
        let mut m = Meter::new(SEARCH_BITS);
        let mut heap = BinaryHeap::new();
        let mut sequence = 0usize;
        let mut push = |heap: &mut BinaryHeap<Node>,
                        charts: [usize; 2],
                        cell: Cell,
                        upper: Option<&Rational>,
                        m: &mut Meter| {
            if let Some(a) = pair.assess(charts, &cell, upper, m) {
                heap.push(Node {
                    bound: a.bound,
                    sequence,
                    charts,
                    cell,
                    weights: a.weights,
                    jac: a.jac,
                });
                sequence += 1;
            }
        };
        for qx in 0..pair.x.charts() {
            for qy in 0..pair.y.charts() {
                let mut cell = pair.x.start();
                cell.extend(pair.y.start());
                push(&mut heap, [qx, qy], cell, self.bound().as_ref(), &mut m);
            }
        }
        m.charge(dec)?;
        let mut unique: Vec<([usize; 2], Cell)> = Vec::new();
        let mut low: Option<Rational> = None;
        let mut exhausted = false;
        let data = (&pair.x.data, &pair.y.data);
        while let Some(node) = heap.pop() {
            dec.check_deadline()?;
            if self.bound().is_some_and(|u| node.bound >= u) {
                continue;
            }
            if self.cells == 0 {
                exhausted = true;
                low = Some(low.map_or(node.bound.clone(), |l| l.min(node.bound.clone())));
                continue;
            }
            self.cells -= 1;
            let step = pair.krawczyk(node.charts, &node.cell, Some(&node.jac), data, &mut m);
            let contracted = match step.as_ref().map(|k| classify(k, &node.cell)) {
                Some(Step::Empty) => None,
                Some(Step::Unique) => {
                    m.charge(dec)?;
                    self.found(
                        &pair,
                        node.charts,
                        node.cell,
                        &mut unique,
                        &mut low,
                        accuracy,
                        dec,
                    )?;
                    continue;
                }
                Some(Step::Contracted(c)) => Some(c),
                None => Some(node.cell.clone()),
            };
            let Some(cell) = contracted else {
                m.charge(dec)?;
                continue;
            };
            if size(&cell).mul_pow2(1) <= size(&node.cell) {
                // Retry at the visited box's size around the shrunken one.
                let mid = centre(&cell);
                let around: Cell = mid
                    .iter()
                    .zip(&node.cell)
                    .map(|(c, [lo, hi])| {
                        let r = (hi - lo).mul_pow2(-1);
                        [c - &r, c + &r]
                    })
                    .collect();
                let proven = pair
                    .krawczyk(node.charts, &around, None, data, &mut m)
                    .is_some_and(|k| matches!(classify(&k, &around), Step::Unique));
                m.charge(dec)?;
                if proven {
                    let covers = inside(&around, &cell);
                    self.found(
                        &pair,
                        node.charts,
                        around,
                        &mut unique,
                        &mut low,
                        accuracy,
                        dec,
                    )?;
                    if covers {
                        continue;
                    }
                }
            }
            let axis = (0..cell.len())
                .max_by(|&a, &b| {
                    let wa = (&cell[a][1] - &cell[a][0]) * node.weights[a].clone();
                    let wb = (&cell[b][1] - &cell[b][0]) * node.weights[b].clone();
                    wa.cmp(&wb).then(b.cmp(&a))
                })
                .expect("a box has an axis");
            let half = (&cell[axis][0] + &cell[axis][1]).mul_pow2(-1);
            let (mut left, mut right) = (cell.clone(), cell);
            left[axis][1] = half.clone();
            right[axis][0] = half;
            for part in [left, right] {
                if unique
                    .iter()
                    .any(|(c, b)| *c == node.charts && inside(b, &part))
                {
                    continue;
                }
                push(&mut heap, node.charts, part, self.bound().as_ref(), &mut m);
            }
            m.charge(dec)?;
        }
        Ok(Some(Settled { low, exhausted }))
    }

    /// Records a box proven to hold exactly one solution: refines it, unless
    /// an earlier proof's box already holds it, and resolves it.
    #[allow(clippy::too_many_arguments)]
    fn found(
        &mut self,
        pair: &Pair,
        charts: [usize; 2],
        cell: Cell,
        unique: &mut Vec<([usize; 2], Cell)>,
        low: &mut Option<Rational>,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let mut m = Meter::new(SEARCH_BITS);
        let star = pair.refine(charts, cell.clone(), (&pair.x.data, &pair.y.data), &mut m);
        m.charge(dec)?;
        let seen = unique.iter().any(|(c, b)| *c == charts && inside(b, &star));
        unique.push((charts, cell));
        if seen {
            return Ok(());
        }
        self.resolve(pair, charts, star, low, accuracy, dec)
    }

    /// Bounds one solution's admitted offset pairs: an exact witness for the
    /// upper bound when both trims admit it, the enclosure for the lower.
    #[allow(clippy::too_many_arguments)]
    fn resolve(
        &mut self,
        pair: &Pair,
        charts: [usize; 2],
        mut star: Cell,
        low: &mut Option<Rational>,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let lower = |low: &mut Option<Rational>, v: Rational| {
            *low = Some(low.take().map_or(v.clone(), |l| l.min(v)));
        };
        let mut bits = SEARCH_BITS;
        let mut data = (pair.x.data.clone(), pair.y.data.clone());
        let mut m = Meter::new(bits);
        let mut found = pair.solution(charts, &star, (&data.0, &data.1), &mut m);
        m.charge(dec)?;
        for next in LADDER {
            let upper = self.bound();
            let wide = found.combos.iter().any(|(_, _, v, _)| {
                upper.as_ref().is_none_or(|u| v.lo() < *u)
                    && !accuracy.reached(&crate::Enclosure {
                        lo: v.lo(),
                        hi: v.hi(),
                    })
            });
            if !wide || next > dec.policy().max_precision_bits {
                break;
            }
            bits = next;
            data = (
                pair.x.exact.enclose(bits, dec)?,
                pair.y.exact.enclose(bits, dec)?,
            );
            let mut m = Meter::new(bits);
            star = pair.refine(charts, star, (&data.0, &data.1), &mut m);
            found = pair.solution(charts, &star, (&data.0, &data.1), &mut m);
            m.charge(dec)?;
        }
        let mut m = Meter::new(bits);
        // A cone's turning equation was divided by its height: a solution at or
        // below the apex is not on the open stratum.
        if found.heights.iter().any(|h| h.hi() <= Rational::zero()) {
            return Ok(());
        }
        let straddles = found.heights.iter().any(|h| h.lo() <= Rational::zero());
        if straddles || found.dist.lo() <= Rational::zero() {
            // The cores may meet. Offset pairs need a common normal there.
            if !straddles && no_common_normal(&found.normals, &mut m) {
                m.charge(dec)?;
                return Ok(());
            }
            m.charge(dec)?;
            if let Some(v) = found.combos.iter().map(|(_, _, v, _)| v.lo()).min() {
                lower(low, v);
            }
            return Ok(());
        }
        let at = centre(&star);
        for (sx, sy, value, carriers) in &found.combos {
            if self.bound().is_some_and(|u| value.lo() >= u) {
                continue;
            }
            let Some((pc, qc)) = carriers else {
                lower(low, value.lo());
                continue;
            };
            // A carrier box outside its stratum's box holds no admitted point.
            let outside = [(pc, pair.extents[0]), (qc, pair.extents[1])]
                .iter()
                .any(|(b, e)| e.is_some_and(|e| !overlaps(b, e)));
            if outside {
                continue;
            }
            let (p, q) = pair.witness(charts, &at, *sx, *sy);
            let mut admitted = Some(true);
            for (member, point) in [(pair.members[0], &p), (pair.members[1], &q)] {
                match member.contains(point, self.scratch, dec) {
                    Ok(true) => {}
                    Ok(false) => {
                        admitted = Some(false);
                        break;
                    }
                    Err(_) => {
                        dec.check_deadline()?;
                        admitted = None;
                    }
                }
            }
            match admitted {
                Some(true) => {
                    let distance = root(&p.minus(&q).norm_squared());
                    let hi = dec
                        .enclose_budgeted_at(&distance, bits.max(LADDER[0]))?
                        .hi();
                    self.witness = Some(self.witness.take().map_or(hi.clone(), |w| w.min(hi)));
                    lower(low, value.lo());
                }
                Some(false) => {
                    let reach = &diagonal(pc, &mut m) + &diagonal(qc, &mut m);
                    self.slack = &self.slack + &reach;
                }
                None => lower(low, value.lo()),
            }
        }
        m.charge(dec)
    }
}
