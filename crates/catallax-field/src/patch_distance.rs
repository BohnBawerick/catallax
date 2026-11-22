//! Certified distance from a point to a whole rational spline patch: step 7.1
//! of 4.2a ([Checkpoint 2 decisions, section 4.6], and the milestone note
//! `docs/M42A-FIELD-FILLETS.md`).
//!
//! A patch is the union of its Bézier spans. Each span is a tensor net of
//! positive homogeneous controls `(X, Y, Z, W)` on the unit square. For a point
//! `p`, the squared distance to the span's point `X / W` is `N / D` with
//! `N = |X - W p|^2` and `D = W^2`. Both are polynomials of twice the span's
//! degrees. On any parameter cell their Bernstein coefficients follow from the
//! cell's net by the product rule for Bernstein polynomials (Farouki and Rajan,
//! "Algorithms for polynomials in Bernstein form", CAGD 5, 1988). When every
//! coefficient `D_k` is positive, `N / D = sum_k (D_k B_k / D) (N_k / D_k)` is a
//! convex combination of the ratios `N_k / D_k`. So the smallest ratio bounds
//! `N / D` from below on the cell, and the four corner ratios are values of
//! `N / D`. Both bounds approach the cell's minimum as the square of the cell
//! size (Garloff, "Convergent bounds for the range of multivariate
//! polynomials", Interval Mathematics 1985, LNCS 212).
//!
//! The search is a best-first branch and bound. It always halves the cell with
//! the smallest lower bound, across the parameter side whose controls spread
//! wider, and drops a cell whose lower bound exceeds the smallest corner value
//! found. It stops when the square roots of the smallest lower bound and the
//! smallest corner value differ by at most the requested width. All arithmetic
//! is outward-rounded binary64 intervals, so the enclosure holds whatever the
//! rounding. A query is one thread's sequential work, so its answer and its
//! work are the same at every worker count.
//!
//! OpenCASCADE (`Extrema_GenExtPS`) and SMLib (`SmSurfPtGlobalSolver`, in
//! `SmSurfacePointSolve.cpp`) take the nearest sample of a grid, or of the
//! control polygon, and polish it with Newton steps. Neither proves that no
//! other point of the surface is nearer. This routine replaces the guess with
//! the subdivision bound, which is a proof. It needs no Newton step, since the
//! corner values reach the minimum as the cells shrink.
//!
//! [Checkpoint 2 decisions, section 4.6]: ../../../docs/checkpoints/checkpoint-2-decisions.md

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use catallax_numbers::{Decider, Interval, Rational, Real, Refusal, StrictSign};
use catallax_shapes::{RPoint3, RationalPatch};

/// The default width of a distance enclosure: `2^-24` mm, about 0.06 nm. It is
/// the largest power of two below a thousandth of the decisions' 100 nm spine
/// tube (section 4.6, step 7.3), so a support's distance uses a negligible part
/// of the tube's error budget.
pub const WIDTH: f64 = 1.0 / 16_777_216.0;

/// Most cells one query evaluates, the panel cap of the rolling region's cover
/// (`rolling_region`, 4,096 panels). A query that reaches it returns the
/// enclosure it has, marked as not converged.
pub const MAX_CELLS: usize = 4096;

/// Deepest halving of a span, the binary64 significand: a cell halved more
/// often has a side below one unit in the last place of the span's parameter.
const MAX_DEPTH: u32 = 52;

fn refusal(detail: &str) -> Refusal {
    Refusal::new("field.patch-distance", detail)
}

/// A rational patch prepared for distance queries: the outward enclosures of
/// its spans' homogeneous Bézier controls.
#[derive(Clone, Debug)]
pub struct PatchDistance {
    degrees: [usize; 2],
    spans: Vec<Vec<[Interval; 4]>>,
    /// `products[axis][i][j]` = `C(n, i) C(n, j) / C(2n, i + j)` for the axis
    /// degree `n`: the coefficient of `B_{i+j}^{2n}` in `B_i^n B_j^n`.
    products: [Vec<Vec<Interval>>; 2],
}

/// Bounds on the distance from one point to a patch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Enclosure {
    /// Bounds on the squared distance, in square millimetres.
    pub squared: Interval,
    /// Bounds on the distance, in millimetres: the square roots of `squared`.
    pub distance: Interval,
    /// Cells whose bounds the query computed, spans included.
    pub cells: usize,
    /// Whether `distance` is no wider than the requested width. A query that
    /// reached `MAX_CELLS` or `MAX_DEPTH` first returns wider, still valid,
    /// bounds.
    pub converged: bool,
}

/// One parameter cell of a span on the search heap.
struct Cell {
    lower: f64,
    order: u64,
    depth: u32,
    net: Vec<[Interval; 4]>,
}

impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Cell {}

impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cell {
    /// The heap pops the smallest lower bound first, the earliest cell on a tie.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .lower
            .total_cmp(&self.lower)
            .then_with(|| other.order.cmp(&self.order))
    }
}

fn binomial(n: usize, k: usize) -> Rational {
    (0..k).fold(Rational::one(), |out, i| {
        &(&out * &Rational::from_i64((n - i) as i64)) / &Rational::from_i64((i + 1) as i64)
    })
}

fn products(n: usize) -> Vec<Vec<Interval>> {
    (0..=n)
        .map(|i| {
            (0..=n)
                .map(|j| {
                    let ratio = &(&binomial(n, i) * &binomial(n, j)) / &binomial(2 * n, i + j);
                    Interval::from_rational(&ratio)
                })
                .collect()
        })
        .collect()
}

fn enclose(value: &Real, dec: &mut Decider) -> Result<Interval, Refusal> {
    Ok(match value.rational_ref() {
        Some(rational) => Interval::from_rational(rational),
        None => dec.enclose_budgeted_at(value, 128)?.to_interval(),
    })
}

fn half() -> Interval {
    Interval::point(0.5).expect("a half")
}

impl PatchDistance {
    /// Encloses the exact homogeneous Bézier controls of every span of `patch`.
    pub fn new(patch: &RationalPatch, dec: &mut Decider) -> Result<Self, Refusal> {
        let degrees = patch.carrier().degrees();
        if degrees.contains(&0) {
            return Err(refusal("a patch of degree zero along a side has no area"));
        }
        let mut spans = Vec::new();
        for (_, controls) in patch.bezier_spans(dec)? {
            let mut net = Vec::with_capacity(controls.len());
            for control in controls {
                let mut point = [Interval::point(0.0).expect("zero"); 4];
                for (out, value) in point.iter_mut().zip(control) {
                    *out = enclose(value, dec)?;
                }
                if point[3].lo() <= 0.0 {
                    return Err(refusal("a span weight is not proved positive"));
                }
                net.push(point);
            }
            spans.push(net);
        }
        if spans.is_empty() {
            return Err(refusal("the patch has no span"));
        }
        Ok(Self {
            degrees,
            spans,
            products: degrees.map(products),
        })
    }

    /// Bounds on the patch: the hull of its Euclidean controls, which holds
    /// the patch since every point is a convex combination of them with
    /// positive weights.
    pub fn hull(&self) -> [Interval; 3] {
        let mut out: Option<[Interval; 3]> = None;
        for c in self.spans.iter().flatten() {
            let p: [Interval; 3] = std::array::from_fn(|k| c[k] / c[3]);
            out = Some(match out {
                None => p,
                Some(o) => std::array::from_fn(|k| o[k].hull(&p[k])),
            });
        }
        out.expect("a span")
    }

    /// Boxes whose union holds the patch: the control hulls of each span's
    /// 64 cells after three halvings along each side.
    pub fn boxes(&self) -> Vec<[Interval; 3]> {
        let mut out = Vec::new();
        for span in &self.spans {
            let mut cells = vec![span.clone()];
            for depth in 0..6 {
                cells = cells
                    .iter()
                    .flat_map(|net| halve(self.degrees, net, depth % 2 == 1))
                    .collect();
            }
            for net in cells {
                let mut hull: Option<[Interval; 3]> = None;
                for c in &net {
                    let p: [Interval; 3] = std::array::from_fn(|k| c[k] / c[3]);
                    hull = Some(match hull {
                        None => p,
                        Some(h) => std::array::from_fn(|k| h[k].hull(&p[k])),
                    });
                }
                out.push(hull.expect("a control"));
            }
        }
        out
    }

    /// Tentative squared distance from `q` to the patch and the side of the
    /// natural normal `S_u x S_v` at the nearest sample, from binary64
    /// arithmetic only: a 9 by 9 grid on each span, then eight 5 by 5 grids,
    /// each half as wide, round the best sample, as SMLib's
    /// `FindRoughGuess` zooms its grids.
    pub fn fast_nearest(&self, q: [f64; 3]) -> (f64, f64) {
        let (distance, side, _) = self.fast_sample(q);
        (distance, side)
    }

    pub(crate) fn fast_sample(&self, q: [f64; 3]) -> (f64, f64, [f64; 3]) {
        let [m, n] = self.degrees;
        let rows = n + 1;
        let mut best = (f64::INFINITY, 1.0, q);
        for span in &self.spans {
            let net: Vec<[f64; 4]> = span.iter().map(|c| c.map(|x| x.mid())).collect();
            let at = |u: f64, v: f64| -> [f64; 3] {
                let column: Vec<[f64; 4]> = (0..=m)
                    .map(|i| casteljau(&net[i * rows..(i + 1) * rows], v))
                    .collect();
                let h = casteljau(&column, u);
                std::array::from_fn(|k| h[k] / h[3])
            };
            let gap = |u: f64, v: f64| {
                let s = at(u, v);
                (0..3).map(|k| (q[k] - s[k]).powi(2)).sum::<f64>()
            };
            let mut here = (f64::INFINITY, 0.0, 0.0);
            for i in 0..=8 {
                for j in 0..=8 {
                    let (u, v) = (f64::from(i) / 8.0, f64::from(j) / 8.0);
                    let value = gap(u, v);
                    if value < here.0 {
                        here = (value, u, v);
                    }
                }
            }
            let mut step = 1.0 / 16.0;
            for _ in 0..8 {
                let (_, cu, cv) = here;
                for i in -2..=2 {
                    for j in -2..=2 {
                        let u = (cu + f64::from(i) * step).clamp(0.0, 1.0);
                        let v = (cv + f64::from(j) * step).clamp(0.0, 1.0);
                        let value = gap(u, v);
                        if value < here.0 {
                            here = (value, u, v);
                        }
                    }
                }
                step /= 2.0;
            }
            if here.0 < best.0 {
                let (_, u, v) = here;
                let h = 1.0 / 1024.0;
                let s = at(u, v);
                let du = at((u + h).min(1.0), v);
                let du0 = at((u + h).min(1.0) - h, v);
                let dv = at(u, (v + h).min(1.0));
                let dv0 = at(u, (v + h).min(1.0) - h);
                let a: [f64; 3] = std::array::from_fn(|k| du[k] - du0[k]);
                let b: [f64; 3] = std::array::from_fn(|k| dv[k] - dv0[k]);
                let normal = [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ];
                let side = (0..3).map(|k| (q[k] - s[k]) * normal[k]).sum::<f64>();
                best = (here.0, side, s);
            }
        }
        best
    }

    /// The distance from `p` to the patch, within `WIDTH`.
    pub fn distance(&self, p: &RPoint3, dec: &mut Decider) -> Result<Enclosure, Refusal> {
        self.distance_within(p, WIDTH, dec)
    }

    /// The distance from `p` to the patch, within `width` millimetres. Each
    /// cell evaluated charges one unit of work to `dec`.
    pub fn distance_within(
        &self,
        p: &RPoint3,
        width: f64,
        dec: &mut Decider,
    ) -> Result<Enclosure, Refusal> {
        if width.is_nan() || width <= 0.0 {
            return Err(refusal("the requested width must be positive"));
        }
        let p = [&p.x, &p.y, &p.z].map(Interval::from_rational);
        self.distance_enclosed_within(&p, width, dec)
    }

    pub(crate) fn distance_enclosed(
        &self,
        p: &[Interval; 3],
        dec: &mut Decider,
    ) -> Result<Enclosure, Refusal> {
        self.distance_enclosed_within(p, WIDTH, dec)
    }

    fn distance_enclosed_within(
        &self,
        p: &[Interval; 3],
        width: f64,
        dec: &mut Decider,
    ) -> Result<Enclosure, Refusal> {
        let mut heap = BinaryHeap::new();
        let mut upper = f64::INFINITY;
        let mut order = 0u64;
        let mut cells = 0usize;
        for net in &self.spans {
            dec.charge_work(1)?;
            cells += 1;
            let (lower, corner) = self.bounds(net, p)?;
            upper = upper.min(corner);
            heap.push(Cell {
                lower,
                order,
                depth: 0,
                net: net.clone(),
            });
            order += 1;
        }
        // The smallest lower bound of a cell that may not be halved again.
        let mut settled = f64::INFINITY;
        // Every point of the patch lies in a cell on the heap, a settled cell
        // or a dropped one, and a dropped cell's points are all farther than
        // a corner. So the smallest lower bound on the heap or settled is at
        // most the distance, which is at most the smallest corner value.
        let enclosure = |lower: f64, upper: f64, cells: usize| {
            let squared = Interval::new(lower, upper)
                .ok_or_else(|| refusal("the lower bound passed a corner value"))?;
            let distance = squared.sqrt().expect("non-negative bounds");
            let spread = Interval::point(distance.hi()).expect("finite")
                - Interval::point(distance.lo()).expect("finite");
            Ok::<_, Refusal>(Enclosure {
                squared,
                distance,
                cells,
                converged: spread.hi() <= width,
            })
        };
        loop {
            let Some(cell) = heap.pop() else {
                return enclosure(settled, upper, cells);
            };
            let lower = cell.lower.min(settled);
            let found = enclosure(lower, upper, cells)?;
            if found.converged {
                return Ok(found);
            }
            if cell.depth >= MAX_DEPTH {
                settled = settled.min(cell.lower);
                continue;
            }
            if cells + 2 > MAX_CELLS {
                return Ok(found);
            }
            for net in self.split(&cell.net) {
                dec.charge_work(1)?;
                cells += 1;
                let (lower, corner) = self.bounds(&net, p)?;
                upper = upper.min(corner);
                if lower <= upper {
                    heap.push(Cell {
                        lower,
                        order,
                        depth: cell.depth + 1,
                        net,
                    });
                    order += 1;
                }
            }
        }
    }

    /// The smallest Bernstein ratio of `N / D` on a cell's net, a lower bound
    /// on the squared distance there, and the smallest upper end of its four
    /// corner values.
    fn bounds(&self, net: &[[Interval; 4]], p: &[Interval; 3]) -> Result<(f64, f64), Refusal> {
        let [m, n] = self.degrees;
        let rows = n + 1;
        // The Bernstein net of `X - W p`, one vector per control.
        let offsets: Vec<[Interval; 3]> = net
            .iter()
            .map(|c| std::array::from_fn(|axis| c[axis] - c[3] * p[axis]))
            .collect();
        let zero = Interval::point(0.0).expect("zero");
        let mut lower = f64::INFINITY;
        let mut corner = f64::INFINITY;
        for k in 0..=2 * m {
            for l in 0..=2 * n {
                let mut numerator = zero;
                let mut denominator = zero;
                for i in k.saturating_sub(m)..=k.min(m) {
                    for j in l.saturating_sub(n)..=l.min(n) {
                        let a = i * rows + j;
                        let b = (k - i) * rows + (l - j);
                        let weight = self.products[0][i][k - i] * self.products[1][j][l - j];
                        let dot = offsets[a][0] * offsets[b][0]
                            + offsets[a][1] * offsets[b][1]
                            + offsets[a][2] * offsets[b][2];
                        numerator = numerator + weight * dot;
                        denominator = denominator + weight * (net[a][3] * net[b][3]);
                    }
                }
                if denominator.lo() <= 0.0 {
                    return Err(refusal("a weight product is not proved positive"));
                }
                let ratio = numerator / denominator;
                lower = lower.min(ratio.lo());
                if (k == 0 || k == 2 * m) && (l == 0 || l == 2 * n) {
                    corner = corner.min(ratio.hi());
                }
            }
        }
        // `N / D` is a sum of squares over a square.
        Ok((lower.max(0.0), corner))
    }

    /// Halves a cell's net at the middle of the side whose controls spread
    /// wider.
    fn split(&self, net: &[[Interval; 4]]) -> [Vec<[Interval; 4]>; 2] {
        let [m, n] = self.degrees;
        let rows = n + 1;
        let point = |c: &[Interval; 4]| -> [f64; 3] {
            let w = c[3].mid();
            std::array::from_fn(|axis| c[axis].mid() / w)
        };
        let spread = |a: &[Interval; 4], b: &[Interval; 4]| -> f64 {
            let (a, b) = (point(a), point(b));
            (0..3).map(|axis| (a[axis] - b[axis]).abs()).sum()
        };
        let across_u = (0..=n)
            .map(|j| spread(&net[j], &net[m * rows + j]))
            .fold(0.0, f64::max);
        let across_v = (0..=m)
            .map(|i| spread(&net[i * rows], &net[i * rows + n]))
            .fold(0.0, f64::max);
        halve(self.degrees, net, across_v > across_u)
    }

    /// The signed number of points where the open ray `p + t d`, `t > 0`,
    /// crosses the patch: one for each crossing where `(S_u x S_v) . d` is
    /// positive, minus one for each where it is negative. `None` where the
    /// search cannot prove every crossing transverse and inside one span, as
    /// where the ray meets the patch's border, touches the patch, or starts on
    /// it; another ray decides those. Each cell charges one unit of work.
    ///
    /// With `e1 = d x a` for the coordinate axis `a` along which `d` is
    /// shortest and `e2 = d x e1`, the patch point `X / W` lies on the ray's
    /// line where `F_k = (X - W p) . e_k` vanish, `k = 1, 2`, and on the ray
    /// itself where also `G = (X - W p) . d` is positive. All three are
    /// polynomials with Bernstein nets on each cell. A cell where `F_1` or
    /// `F_2` has coefficients of one strict sign holds no crossing. Otherwise
    /// the Krawczyk operator over the cell (Krawczyk 1969; Moore, Kearfott and
    /// Cloud, "Introduction to Interval Analysis", 2009, section 8.2) either
    /// excludes a root, proves exactly one root strictly inside the cell, or
    /// asks for halving. At a root `X - W p` is parallel to `d`, so the
    /// Jacobian of `(F_1, F_2)` is `W^2 |e1|^2 (S_u x S_v) . d` in
    /// determinant, since `e1 x e2 = |e1|^2 d`: its sign over the cell is the
    /// crossing's sign, and `G`'s sign says whether the root lies ahead.
    ///
    /// OpenCASCADE (`IntCurveSurface_Inter`, `IntCurveSurface_TheExactHInter`)
    /// and SMLib (`SmRayTracer.cpp`) find line-surface hits by sampling and
    /// Newton steps, which proves no count. The field needs the exact count for
    /// the winding number, so it counts proved roots only.
    pub fn crossings(
        &self,
        p: &[Rational; 3],
        d: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Option<i64>, Refusal> {
        let shortest = (0..3)
            .min_by(|&i, &j| d[i].abs().cmp(&d[j].abs()))
            .expect("three axes");
        let mut a: [Rational; 3] = std::array::from_fn(|_| Rational::zero());
        a[shortest] = Rational::one();
        let e1 = cross(d, &a);
        let e2 = cross(d, &e1);
        if e1.iter().all(Rational::is_zero) {
            return Err(refusal("a ray direction must be nonzero"));
        }
        let [p, d, e1, e2] = [p, d, &e1, &e2].map(|v| v.each_ref().map(Interval::from_rational));
        let dot = |a: &[Interval; 3], b: &[Interval; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let mut count = 0i64;
        let mut cells = 0usize;
        for span in &self.spans {
            let start: Vec<[Interval; 3]> = span
                .iter()
                .map(|c| {
                    let offset: [Interval; 3] = std::array::from_fn(|k| c[k] - c[3] * p[k]);
                    [dot(&offset, &e1), dot(&offset, &e2), dot(&offset, &d)]
                })
                .collect();
            // Depth first, low half before high half: a fixed order.
            let mut stack = vec![(start, 0u32)];
            while let Some((net, depth)) = stack.pop() {
                dec.charge_work(1)?;
                cells += 1;
                if cells > MAX_CELLS {
                    return Ok(None);
                }
                match self.root(&net) {
                    Root::None => (),
                    Root::One { sign, ahead } => {
                        if ahead {
                            count += sign;
                        }
                    }
                    Root::Unknown => {
                        if depth >= CROSSING_DEPTH {
                            return Ok(None);
                        }
                        let [low, high] = halve(self.degrees, &net, depth % 2 == 1);
                        stack.push((high, depth + 1));
                        stack.push((low, depth + 1));
                    }
                }
            }
        }
        Ok(Some(count))
    }

    /// What one cell of `(F_1, F_2, G)` nets holds, in the cell's own unit
    /// square.
    fn root(&self, net: &[[Interval; 3]]) -> Root {
        let one_sign =
            |k: usize| net.iter().all(|c| c[k].lo() > 0.0) || net.iter().all(|c| c[k].hi() < 0.0);
        if one_sign(0) || one_sign(1) {
            return Root::None;
        }
        let [m, n] = self.degrees;
        let rows = n + 1;
        // Bernstein nets of the partials: `m` times the forward differences
        // along `u`, `n` times those along `v`; their hulls enclose the
        // Jacobian over the cell.
        let hull = |values: &mut dyn Iterator<Item = Interval>| {
            values.reduce(|a, b| a.hull(&b)).expect("a nonempty net")
        };
        let scale = |k: usize| Interval::point(k as f64).expect("a degree");
        let jacobian: [[Interval; 2]; 2] = std::array::from_fn(|k| {
            [
                hull(
                    &mut (0..m)
                        .flat_map(|i| (0..=n).map(move |j| (i, j)))
                        .map(|(i, j)| {
                            scale(m) * (net[(i + 1) * rows + j][k] - net[i * rows + j][k])
                        }),
                ),
                hull(
                    &mut (0..=m)
                        .flat_map(|i| (0..n).map(move |j| (i, j)))
                        .map(|(i, j)| scale(n) * (net[i * rows + j + 1][k] - net[i * rows + j][k])),
                ),
            ]
        });
        let [[a, b], [c, d]] = jacobian.map(|row| row.map(|x| x.mid()));
        let det = a * d - b * c;
        if !det.is_finite() || det == 0.0 {
            return Root::Unknown;
        }
        // Any matrix serves as Krawczyk's preconditioner; the inverse of the
        // Jacobian's midpoint makes the operator contract.
        let y = [[d / det, -b / det], [-c / det, a / det]].map(|row| {
            row.map(|x| Interval::point(x).unwrap_or_else(|| Interval::point(0.0).expect("zero")))
        });
        let centre = centre(self.degrees, net);
        let radius = Interval::new(-0.5, 0.5).expect("a radius");
        let one = |row: usize, column: usize| {
            Interval::point(if row == column { 1.0 } else { 0.0 }).expect("a unit")
        };
        let image: [Interval; 2] = std::array::from_fn(|row| {
            let step = y[row][0] * centre[0] + y[row][1] * centre[1];
            let spread = |column: usize| {
                one(row, column)
                    - (y[row][0] * jacobian[0][column] + y[row][1] * jacobian[1][column])
            };
            half() - step + spread(0) * radius + spread(1) * radius
        });
        if image.iter().any(|k| k.hi() < 0.0 || k.lo() > 1.0) {
            return Root::None;
        }
        if image.iter().any(|k| k.lo() <= 0.0 || k.hi() >= 1.0) {
            return Root::Unknown;
        }
        let determinant = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
        let ahead = if net.iter().all(|c| c[2].lo() > 0.0) {
            true
        } else if net.iter().all(|c| c[2].hi() < 0.0) {
            false
        } else {
            return Root::Unknown;
        };
        match determinant.sign() {
            Some(StrictSign::Positive) => Root::One { sign: 1, ahead },
            Some(StrictSign::Negative) => Root::One { sign: -1, ahead },
            None => Root::Unknown,
        }
    }
}

/// De Casteljau's value of a Bernstein row at `t`, in binary64.
fn casteljau(row: &[[f64; 4]], t: f64) -> [f64; 4] {
    let mut row = row.to_vec();
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|w| std::array::from_fn(|c| w[0][c] + t * (w[1][c] - w[0][c])))
            .collect();
    }
    row[0]
}

/// What a crossing cell holds.
enum Root {
    None,
    One { sign: i64, ahead: bool },
    Unknown,
}

/// Most halvings of a span in the crossing search, alternating the sides:
/// twenty per side. A crossing not proved by then lies on or within a few
/// units in the last place of a dyadic line of the span, and a ray in
/// another direction decides.
const CROSSING_DEPTH: u32 = 40;

fn cross(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

/// Halves a tensor Bernstein net of the given degrees, u-major, along `v`
/// or along `u`, by de Casteljau's algorithm at one half.
fn halve<const K: usize>(
    degrees: [usize; 2],
    net: &[[Interval; K]],
    along_v: bool,
) -> [Vec<[Interval; K]>; 2] {
    let [m, n] = degrees;
    let rows = n + 1;
    let (count, lines) = if along_v {
        (n + 1, m + 1)
    } else {
        (m + 1, n + 1)
    };
    let index = |line: usize, at: usize| {
        if along_v {
            line * rows + at
        } else {
            at * rows + line
        }
    };
    let mut low = net.to_vec();
    let mut high = net.to_vec();
    let mut row = Vec::with_capacity(count);
    for line in 0..lines {
        row.clear();
        row.extend((0..count).map(|at| net[index(line, at)]));
        // After pass `r`, `row[0]` is the low half's control `r` and
        // `row[count - 1 - r]` the high half's control `count - 1 - r`.
        low[index(line, 0)] = row[0];
        high[index(line, count - 1)] = row[count - 1];
        for r in 1..count {
            for at in 0..count - r {
                row[at] = std::array::from_fn(|c| (row[at][c] + row[at + 1][c]) * half());
            }
            low[index(line, r)] = row[0];
            high[index(line, count - 1 - r)] = row[count - 1 - r];
        }
    }
    [low, high]
}

/// The net's value at the middle of its square, by de Casteljau's algorithm
/// at one half along `v`, then along `u`.
fn centre<const K: usize>(degrees: [usize; 2], net: &[[Interval; K]]) -> [Interval; K] {
    let [m, n] = degrees;
    let rows = n + 1;
    let middle = |mut row: Vec<[Interval; K]>| {
        while row.len() > 1 {
            row = row
                .windows(2)
                .map(|w| std::array::from_fn(|c| (w[0][c] + w[1][c]) * half()))
                .collect();
        }
        row[0]
    };
    middle(
        (0..=m)
            .map(|i| middle(net[i * rows..(i + 1) * rows].to_vec()))
            .collect(),
    )
}
