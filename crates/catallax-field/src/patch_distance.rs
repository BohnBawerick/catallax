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

use catallax_numbers::{Decider, Interval, Rational, Real, Refusal};
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
        let mut heap = BinaryHeap::new();
        let mut upper = f64::INFINITY;
        let mut order = 0u64;
        let mut cells = 0usize;
        for net in &self.spans {
            dec.charge_work(1)?;
            cells += 1;
            let (lower, corner) = self.bounds(net, &p)?;
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
                let (lower, corner) = self.bounds(&net, &p)?;
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
    /// wider, by de Casteljau's algorithm at one half.
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
        let along_v = across_v > across_u;
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
}
