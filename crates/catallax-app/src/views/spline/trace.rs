//! The zero set of a carrier's contour function, traced in binary64.
//!
//! The method is SMLib's `SmSurfaceSilhouette`, a `SmSurfaceTracer`: start
//! points on boundary curves, interior loops found through a test on the
//! surface normals, and a march along `(F_v, -F_u)` with a local point solve
//! after each step. OpenCASCADE's `Contap_Contour` marches the same function
//! from start points on the face restrictions and on an interior sample grid
//! (`Contap_TheSearchInside`). Here both kinds of start points come from one
//! grid whose cells provably hold no closed loop:
//!
//! 1. Every nonzero span of the carrier is a Bezier patch with the Bernstein
//!    form of `K . d` (`bernstein.rs`). A cell is split into quarters until
//!    `K . d` has one sign on it (no generator), or one of its partial
//!    derivatives does (no closed loop: a closed level curve encloses a
//!    critical point). This is the convex hull exclusion SMLib's normal cone
//!    test performs.
//! 2. The grid lines are the span boundaries and, in each span, the finest
//!    split any patch in that row or column needed. Every grid cell then lies
//!    inside one tested cell, so every generator, open or closed, crosses a
//!    grid line.
//! 3. Each grid segment's roots are isolated by Descartes' rule of signs in
//!    Bernstein form (Lane and Riesenfeld, "Bounds on a polynomial", BIT 21,
//!    1981; Mourrain, Rouillier and Roy, "The Bernstein basis and real root
//!    isolation", 2005) and refined by bisection.
//! 4. Tracing starts from roots on the domain boundary, inward, then from
//!    roots still unused, which lie on closed loops. Each step predicts along
//!    the tangent, corrects by Newton's method toward `F = 0` along the
//!    gradient, and marks every root whose grid segment it crosses as used. A
//!    trace ends where it crosses the domain boundary at a root, or at its
//!    own start.
//!
//! The traced points only guide the fit (`spline.rs`); the fitted carrier is
//! measured against the true generator, solved afresh at each parameter.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Refusal};

use super::super::float::{above, dot};
use super::bernstein::{Carrier, Poly1, Poly2};

/// Source: SMLib's surface tracer bounds the turn of each step by its
/// angular tolerance, "usually 10 or 20 degrees" (`SmSurfaceSilhouette.cpp`,
/// `ComputePointValues`); the smaller value.
const TURN_DEGREES: f64 = 10.0;

/// Source: SMLib `SmSurfaceSilhouette::LocalPointSolve`: the first step of a
/// trace is at most one hundredth of the domain's size.
const FIRST_STEP_SHARE: f64 = 100.0;

/// Source: judgement. Newton's method converges quadratically from a
/// prediction within [`TURN_DEGREES`] of the curve; eight steps only end a
/// solve that does not converge, and the step is then halved.
const NEWTON_STEPS: usize = 8;

/// Source: the 53-bit significand of IEEE 754 binary64 plus margin: a
/// bisection halves its bracket at most this often before the bracket is
/// below binary64 resolution of any parameter in it.
const BISECTIONS: usize = 64;

fn refuse(code: &'static str, detail: &str) -> Refusal {
    Refusal::new(code, detail.to_owned())
}

/// A grid line: a parameter value on one axis, with the span holding it and
/// the local parameter there.
#[derive(Clone, Copy, Debug)]
struct Line {
    value: f64,
    span: usize,
    local: f64,
}

/// A root of the contour function on one grid segment.
#[derive(Clone, Debug)]
pub(super) struct Seed {
    /// The axis of the line's own coordinate: 0 for a line `u = const`.
    axis: usize,
    line: usize,
    /// Where it lies along the line, and the interval isolating it there.
    at: f64,
    interval: [f64; 2],
    uv: [f64; 2],
    /// The trace that used it.
    owner: Option<usize>,
}

/// A traced point and the tangent `(F_v, -F_u)` oriented along the trace.
#[derive(Clone, Copy, Debug)]
pub(super) struct Vertex {
    pub uv: [f64; 2],
    pub tangent: [f64; 2],
}

/// One traced generator: from boundary to boundary, or closed on itself with
/// its last vertex equal to its first.
#[derive(Clone, Debug)]
pub(super) struct Branch {
    pub vertices: Vec<Vertex>,
    pub closed: bool,
    /// For an open branch, the domain sides its two ends lie on, as
    /// `(axis, high)`: axis 0 for a side `u = const`, `high` for the upper
    /// end of that parameter.
    pub sides: Option<[(usize, bool); 2]>,
}

/// The binary64 resolution of each parameter axis: the unit roundoff of the
/// domain's largest end, doubled as in `fit.rs`.
pub(super) fn resolution(carrier: &Carrier, axis: usize) -> f64 {
    let b = &carrier.breaks[axis];
    b[0].abs().max(b[b.len() - 1].abs()) * f64::EPSILON
}

struct Meter<'a> {
    work: u64,
    cap: u64,
    dec: &'a mut Decider,
}

impl Meter<'_> {
    fn charge(&mut self, units: u64) -> Result<(), Refusal> {
        self.work += units;
        self.dec.charge_work(units)?;
        if self.work > self.cap {
            return Err(refuse(
                "app.contour-trace-work",
                "a traced contour generator exceeds its work allowance",
            ));
        }
        Ok(())
    }
}

/// The deepest split a patch needs before every cell holds no closed loop.
fn depth(poly: &Poly2, level: u32, limit: u32, meter: &mut Meter<'_>) -> Result<u32, Refusal> {
    meter.charge(poly.c.len() as u64)?;
    if poly.one_sign() || poly.derivative(0).one_sign() || poly.derivative(1).one_sign() {
        return Ok(level);
    }
    if level >= limit {
        return Err(refuse(
            "app.contour-singular",
            "a contour generator crosses itself or closes in a loop smaller than binary64 \
             resolves",
        ));
    }
    let mut deepest = level;
    let (a, b) = poly.halves(0);
    for half in [a, b] {
        let (c, d) = half.halves(1);
        for quarter in [c, d] {
            deepest = deepest.max(depth(&quarter, level + 1, limit, meter)?);
        }
    }
    Ok(deepest)
}

/// Roots of `poly` over `[a, b]` (local to the segment), ascending, each with
/// its isolating interval. A root at `a` itself is kept; one at `b` belongs
/// to the next interval.
fn isolate(poly: &Poly1, a: f64, b: f64, width: f64, out: &mut Vec<(f64, [f64; 2])>) {
    let c = &poly.0;
    if c[0] == 0.0 {
        out.push((a, [a, a]));
    }
    let sign = |x: f64| if x > 0.0 { 1_i8 } else { -1 };
    let signs: Vec<i8> = c.iter().filter(|x| **x != 0.0).map(|&x| sign(x)).collect();
    let changes = signs.windows(2).filter(|w| w[0] != w[1]).count();
    if changes == 0 {
        return;
    }
    let (first, last) = (signs[0], signs[signs.len() - 1]);
    if changes == 1 {
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        for _ in 0..BISECTIONS {
            let middle = lo + (hi - lo) / 2.0;
            if (hi - lo) * (b - a) <= width || !(lo < middle && middle < hi) {
                break;
            }
            let value = poly.eval(middle);
            if value == 0.0 {
                (lo, hi) = (middle, middle);
                break;
            }
            if sign(value) == first {
                lo = middle;
            } else {
                hi = middle;
            }
        }
        out.push((a + (b - a) * (lo + (hi - lo) / 2.0), [a, b]));
        return;
    }
    if b - a <= width {
        if first != last {
            out.push((a + (b - a) / 2.0, [a, b]));
        }
        return;
    }
    let middle = a + (b - a) / 2.0;
    let (left, right) = poly.split(0.5);
    isolate(&left, a, middle, width, out);
    isolate(&right, middle, b, width, out);
}

/// The grid, its roots and the traced branches of one carrier.
pub(super) struct Tracer<'a> {
    carrier: &'a Carrier,
    lines: [Vec<Line>; 2],
    levels: [Vec<u32>; 2],
    seeds: Vec<Seed>,
    /// Seeds by axis, line and segment.
    index: BTreeMap<(usize, usize, usize), Vec<usize>>,
    floor: f64,
    first_step: f64,
    traces: usize,
}

impl<'a> Tracer<'a> {
    /// Splits, grids and roots the carrier. `diagonal` is the size of its
    /// control box and `floor` rule T1's floor, both in millimetres.
    pub fn new(
        carrier: &'a Carrier,
        diagonal: f64,
        floor: f64,
        cap: u64,
        dec: &mut Decider,
    ) -> Result<(Self, u64), Refusal> {
        let mut meter = Meter { work: 0, cap, dec };
        let spans = [carrier.spans(0), carrier.spans(1)];
        let resolutions = [resolution(carrier, 0), resolution(carrier, 1)];
        let mut levels: [Vec<u32>; 2] = [vec![0; spans[0]], vec![0; spans[1]]];
        for k in 0..spans[0] {
            for l in 0..spans[1] {
                meter.dec.check_deadline()?;
                // Halve no further than twice the binary64 resolution.
                let mut limit = 0;
                let (mut wu, mut wv) = (carrier.width(0, k), carrier.width(1, l));
                while wu > 4.0 * resolutions[0] && wv > 4.0 * resolutions[1] {
                    (wu, wv) = (wu / 2.0, wv / 2.0);
                    limit += 1;
                }
                let level = depth(&carrier.patch(k, l).contour, 0, limit, &mut meter)?;
                levels[0][k] = levels[0][k].max(level);
                levels[1][l] = levels[1][l].max(level);
            }
        }
        let lines: [Vec<Line>; 2] = std::array::from_fn(|axis| {
            let mut out = Vec::new();
            for k in 0..spans[axis] {
                let n = 1_u64 << levels[axis][k];
                let (start, width) = (carrier.breaks[axis][k], carrier.width(axis, k));
                for i in 0..n {
                    let local = i as f64 / n as f64;
                    out.push(Line {
                        value: start + local * width,
                        span: k,
                        local,
                    });
                }
            }
            out.push(Line {
                value: carrier.breaks[axis][spans[axis]],
                span: spans[axis] - 1,
                local: 1.0,
            });
            out
        });
        let mut tracer = Self {
            carrier,
            lines,
            levels,
            seeds: Vec::new(),
            index: BTreeMap::new(),
            floor,
            first_step: diagonal / FIRST_STEP_SHARE,
            traces: 0,
        };
        tracer.root(resolutions, &mut meter)?;
        Ok((tracer, meter.work))
    }

    /// Every root on every grid segment, by axis, line and segment.
    fn root(&mut self, resolutions: [f64; 2], meter: &mut Meter<'_>) -> Result<(), Refusal> {
        for axis in 0..2 {
            let other = 1 - axis;
            for (i, line) in self.lines[axis].iter().enumerate() {
                meter.dec.check_deadline()?;
                let segments = self.lines[other].len() - 1;
                for j in 0..segments {
                    let (from, to) = (self.lines[other][j], self.lines[other][j + 1]);
                    let span = from.span;
                    let end = if to.span == span { to.local } else { 1.0 };
                    let (k, l) = if axis == 0 {
                        (line.span, span)
                    } else {
                        (span, line.span)
                    };
                    let contour = &self.carrier.patch(k, l).contour;
                    let along = if axis == 0 {
                        contour.at_s(line.local)
                    } else {
                        contour.at_t(line.local)
                    };
                    let piece = along.restrict(from.local, end);
                    meter.charge((contour.c.len() + 4 * piece.0.len()) as u64)?;
                    let mut roots = Vec::new();
                    let width = resolutions[other] / (to.value - from.value);
                    isolate(&piece, 0.0, 1.0, width, &mut roots);
                    if j + 1 == segments && piece.0[piece.0.len() - 1] == 0.0 {
                        roots.push((1.0, [1.0, 1.0]));
                    }
                    let scale = to.value - from.value;
                    for (r, [lo, hi]) in roots {
                        let at = from.value + r * scale;
                        let mut uv = [0.0; 2];
                        uv[axis] = line.value;
                        uv[other] = at;
                        self.index
                            .entry((axis, i, j))
                            .or_default()
                            .push(self.seeds.len());
                        self.seeds.push(Seed {
                            axis,
                            line: i,
                            at,
                            interval: [from.value + lo * scale, from.value + hi * scale],
                            uv,
                            owner: None,
                        });
                    }
                }
            }
        }
        Ok(())
    }

    fn is_boundary(&self, axis: usize, line: usize) -> bool {
        line == 0 || line + 1 == self.lines[axis].len()
    }

    /// `(F_v, -F_u)` at `uv`, times `sense`.
    fn tangent(&self, uv: [f64; 2], sense: f64) -> [f64; 2] {
        let level = self.carrier.level(uv[0], uv[1]);
        [sense * level.grad[1], -sense * level.grad[0]]
    }

    fn velocity(&self, uv: [f64; 2], tangent: [f64; 2]) -> [f64; 3] {
        let jet = self.carrier.jet(uv[0], uv[1]);
        std::array::from_fn(|c| jet.du[c] * tangent[0] + jet.dv[c] * tangent[1])
    }

    fn cell_width(&self, axis: usize, x: f64) -> f64 {
        let (k, _) = self.carrier.locate(axis, x);
        self.carrier.width(axis, k) / (1_u64 << self.levels[axis][k]) as f64
    }

    /// The seed a step from `p` to `q` meets where it crosses line `line` of
    /// `axis`, if any: the one whose isolating interval holds the crossing,
    /// else the nearest on that segment.
    fn seed_at(&self, axis: usize, line: usize, p: [f64; 2], q: [f64; 2]) -> Option<usize> {
        let other = 1 - axis;
        let x = self.lines[axis][line].value;
        let share = (x - p[axis]) / (q[axis] - p[axis]);
        let y = p[other] + share * (q[other] - p[other]);
        let lines = &self.lines[other];
        let j = lines
            .partition_point(|l| l.value <= y)
            .saturating_sub(1)
            .min(lines.len() - 2);
        let candidates = self.index.get(&(axis, line, j))?;
        candidates
            .iter()
            .copied()
            .find(|&s| self.seeds[s].interval[0] <= y && y <= self.seeds[s].interval[1])
            .or_else(|| {
                candidates.iter().copied().min_by(|&a, &b| {
                    (self.seeds[a].at - y)
                        .abs()
                        .total_cmp(&(self.seeds[b].at - y).abs())
                })
            })
    }

    /// Lines crossed by a step from `p` to `q`, in order along it. A line
    /// through `p` itself is not crossed; one through `q` is.
    fn crossed(&self, p: [f64; 2], q: [f64; 2]) -> Vec<(f64, usize, usize)> {
        let mut out = Vec::new();
        for axis in 0..2 {
            // No line passes both tests when `a` equals `b`.
            let (a, b) = (p[axis], q[axis]);
            for (i, line) in self.lines[axis].iter().enumerate() {
                let x = line.value;
                if (a < x && x <= b) || (b <= x && x < a) {
                    out.push(((x - a) / (b - a), axis, i));
                }
            }
        }
        out.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        out
    }

    /// Traces from seed `start` in the sense `sense`. `None` for a trace from
    /// an interior root that runs into roots already used: a generator traced
    /// before.
    fn march(
        &mut self,
        start: usize,
        sense: f64,
        meter: &mut Meter<'_>,
    ) -> Result<Option<Branch>, Refusal> {
        let cos_turn = libm::cos(TURN_DEGREES.to_radians());
        let interior = !self.is_boundary(self.seeds[start].axis, self.seeds[start].line);
        let trace = self.traces;
        self.traces += 1;
        self.seeds[start].owner = Some(trace);
        let mut p = self.seeds[start].uv;
        let mut vertices = vec![Vertex {
            uv: p,
            tangent: self.tangent(p, sense),
        }];
        let mut length = self.first_step;
        loop {
            meter.dec.check_deadline()?;
            meter.charge(1)?;
            let t = vertices[vertices.len() - 1].tangent;
            let velocity = self.velocity(p, t);
            let speed = dot(velocity, velocity).sqrt();
            if !above(speed, 0.0) || !above(length, self.floor) {
                return Err(refuse(
                    "app.contour-singular",
                    "a traced contour generator reached a singular point",
                ));
            }
            let mut step = [t[0] * length / speed, t[1] * length / speed];
            // At most one grid cell per step on each axis.
            for axis in 0..2 {
                let width = self.cell_width(axis, p[axis]);
                if step[axis].abs() > width {
                    let shrink = width / step[axis].abs();
                    step = [step[0] * shrink, step[1] * shrink];
                }
            }
            let predicted = [p[0] + step[0], p[1] + step[1]];
            let mut q = predicted;
            let mut converged = false;
            for _ in 0..NEWTON_STEPS {
                meter.charge(1)?;
                let level = self.carrier.level(q[0], q[1]);
                let norm = level.grad[0] * level.grad[0] + level.grad[1] * level.grad[1];
                if !above(norm, 0.0) {
                    break;
                }
                let delta = [
                    -level.value * level.grad[0] / norm,
                    -level.value * level.grad[1] / norm,
                ];
                q = [q[0] + delta[0], q[1] + delta[1]];
                let moved = self.velocity(q, delta);
                if dot(moved, moved).sqrt() <= self.floor {
                    converged = true;
                    break;
                }
            }
            let drift = self.velocity(p, [q[0] - predicted[0], q[1] - predicted[1]]);
            let tangent = self.tangent(q, sense);
            let next = self.velocity(q, tangent);
            let turn = dot(velocity, next) / (speed * dot(next, next).sqrt());
            if !converged
                || dot(drift, drift).sqrt() > length / 2.0
                || turn.is_nan()
                || turn < cos_turn
            {
                length /= 2.0;
                continue;
            }
            for (_, axis, line) in self.crossed(p, q) {
                let boundary = self.is_boundary(axis, line);
                let Some(seed) = self.seed_at(axis, line, p, q) else {
                    if boundary {
                        return Err(refuse(
                            "app.contour-trace",
                            "a traced contour generator left the face where no boundary root \
                             was found",
                        ));
                    }
                    continue;
                };
                if seed == start {
                    if vertices.len() > 1 {
                        let uv = self.seeds[start].uv;
                        vertices.push(Vertex {
                            uv,
                            tangent: self.tangent(uv, sense),
                        });
                        return Ok(Some(Branch {
                            vertices,
                            closed: true,
                            sides: None,
                        }));
                    }
                    continue;
                }
                match self.seeds[seed].owner {
                    Some(owner) if owner == trace => {
                        return Err(refuse(
                            "app.contour-trace",
                            "a traced contour generator passed its own start unseen",
                        ));
                    }
                    Some(_) if interior => return Ok(None),
                    Some(_) => {
                        return Err(refuse(
                            "app.contour-trace",
                            "two traced contour generators met one root",
                        ));
                    }
                    None => self.seeds[seed].owner = Some(trace),
                }
                if boundary {
                    let uv = self.seeds[seed].uv;
                    vertices.push(Vertex {
                        uv,
                        tangent: self.tangent(uv, sense),
                    });
                    if interior {
                        // An open generator traced from inside: it was traced
                        // from its boundary roots already.
                        return Ok(None);
                    }
                    let side = |s: &Seed| (s.axis, s.line != 0);
                    return Ok(Some(Branch {
                        vertices,
                        closed: false,
                        sides: Some([side(&self.seeds[start]), side(&self.seeds[seed])]),
                    }));
                }
            }
            vertices.push(Vertex { uv: q, tangent });
            p = q;
            if turn >= libm::cos((TURN_DEGREES / 2.0).to_radians()) {
                length *= 2.0;
            }
        }
    }

    /// Every generator: from boundary roots first, inward, then closed loops
    /// from the roots left.
    pub fn branches(&mut self, cap: u64, dec: &mut Decider) -> Result<(Vec<Branch>, u64), Refusal> {
        let mut meter = Meter { work: 0, cap, dec };
        let mut out = Vec::new();
        let boundary: Vec<usize> = (0..self.seeds.len())
            .filter(|&s| self.is_boundary(self.seeds[s].axis, self.seeds[s].line))
            .collect();
        for s in boundary {
            if self.seeds[s].owner.is_some() {
                continue;
            }
            let Seed { axis, line, uv, .. } = self.seeds[s];
            // Inward: increasing along `axis` from its first line.
            let inward = if line == 0 { 1.0 } else { -1.0 };
            let along = self.tangent(uv, 1.0)[axis];
            if along == 0.0 {
                // The generator only touches the boundary here.
                self.seeds[s].owner = Some(usize::MAX);
                continue;
            }
            let sense = if along * inward > 0.0 { 1.0 } else { -1.0 };
            if let Some(branch) = self.march(s, sense, &mut meter)? {
                out.push(branch);
            }
        }
        for s in 0..self.seeds.len() {
            if self.seeds[s].owner.is_some() {
                continue;
            }
            if let Some(branch) = self.march(s, 1.0, &mut meter)? {
                out.push(branch);
            }
        }
        Ok((out, meter.work))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Rational, Real};
    use catallax_shapes::{BSplineSurface, Point3};

    /// A bicubic bump over `[0, 3]^2`: `z = B(x) B(y)` with `B` the sum of the
    /// two middle basis functions. Height and slope vanish on the boundary,
    /// and the largest slope `B'(x) B(y)` is about 0.94.
    fn bump() -> (Carrier, f64) {
        let q = |s: &str| Rational::parse(s).unwrap();
        let knots: Vec<Rational> = ["0", "0", "0", "0", "1", "2", "3", "3", "3", "3"]
            .map(q)
            .to_vec();
        // Greville abscissae, so that x = u and y = v.
        let grid = ["0", "1/3", "1", "2", "8/3", "3"];
        let mut controls = Vec::new();
        for (i, x) in grid.iter().enumerate() {
            for (j, y) in grid.iter().enumerate() {
                let z = i64::from((2..4).contains(&i) && (2..4).contains(&j));
                controls.push(Point3::new(
                    Real::from(q(x)),
                    Real::from(q(y)),
                    Real::from_i64(z),
                ));
            }
        }
        let mut dec = Decider::default();
        let surface = BSplineSurface::new(
            [3, 3],
            [knots.clone(), knots],
            [6, 6],
            controls,
            vec![Real::one(); 36],
            &mut dec,
        )
        .unwrap();
        // The view (2, 0, 1): N . d = 0 where the slope in x is one half.
        let d = [2.0 / 5.0_f64.sqrt(), 0.0, 1.0 / 5.0_f64.sqrt()];
        let zero = Rational::zero();
        let floor = Rational::parse("1e-7").unwrap();
        let (carrier, _) =
            Carrier::new(&surface, [&zero, &zero, &zero], d, &floor, &mut dec).unwrap();
        (carrier, 1e-7)
    }

    #[test]
    fn an_interior_loop_is_found_and_traced_closed() {
        let (carrier, floor) = bump();
        let mut dec = Decider::default();
        let (mut tracer, _) = Tracer::new(&carrier, 3.0, floor, u64::MAX, &mut dec).unwrap();
        // No root lies on the domain boundary, where the slope is zero.
        assert!(
            tracer
                .seeds
                .iter()
                .all(|s| !tracer.is_boundary(s.axis, s.line))
        );
        let (branches, _) = tracer.branches(u64::MAX, &mut dec).unwrap();
        assert_eq!(branches.len(), 1);
        let branch = &branches[0];
        assert!(branch.closed);
        assert!(branch.vertices.len() > 16);
        let ends = [
            branch.vertices[0].uv,
            branch.vertices[branch.vertices.len() - 1].uv,
        ];
        assert_eq!(
            ends[0].partial_cmp(&ends[1]),
            Some(std::cmp::Ordering::Equal)
        );
        for vertex in &branch.vertices {
            let level = carrier.level(vertex.uv[0], vertex.uv[1]);
            let gradient = level.grad[0].hypot(level.grad[1]);
            assert!(level.value.abs() <= 1e-9 * gradient, "{vertex:?}");
            // On the rising side of the bump.
            assert!(vertex.uv[0] > 0.0 && vertex.uv[0] < 1.5);
        }
    }
}
