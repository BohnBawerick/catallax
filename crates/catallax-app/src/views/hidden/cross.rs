//! Binary64 crossings of two strokes on the sheet, for the pairs the exact
//! cuts do not cover, and the folds and cusps of spline strokes.
//!
//! Crossings follow OpenCASCADE's `IntCurve_IntPolyPolyGen`: localise by
//! polygons, then refine each candidate with Newton's method on the two
//! curves (`IntCurve_ExactIntersectionPoint`). The polygons here are the
//! Bezier control polygons of the strokes' rational pieces, halved until each
//! is flat within rule T1's floor (Lane and Riesenfeld, "A theoretical
//! development for the computer generation and display of piecewise
//! polynomial surfaces", IEEE PAMI 2, 1980; the convex hull property bounds
//! each half, Farin, section 15.3). SMLib's curve intersector subdivides the
//! same way.

use super::super::float::above;
use super::geo::{Geo, Sheet, distance, split};

/// Source: judgement, as in `spline.rs`. Newton's method converges
/// quadratically from a flat leaf; sixteen steps only end one that stalls.
pub(super) const SOLVE_STEPS: usize = 16;

/// Source: the binary64 significand. After 52 halvings a piece's parameter
/// interval is below the resolution of its parameter, so halving further
/// changes nothing.
const DEPTH: usize = 52;

/// A crossing of two strokes on the sheet: the parameters on each and the
/// sheet distance left between the two points.
#[derive(Clone, Copy, Debug)]
pub(super) struct Crossing {
    pub t: f64,
    pub u: f64,
    pub residual: f64,
}

type H2 = [f64; 3];

fn project(sheet: &Sheet, ctrl: &[[f64; 4]]) -> Vec<H2> {
    ctrl.iter()
        .map(|h| {
            let p = [h[0], h[1], h[2]];
            [
                super::super::float::dot(p, sheet.e1),
                super::super::float::dot(p, sheet.e2),
                h[3],
            ]
        })
        .collect()
}

fn cartesian(h: &H2) -> [f64; 2] {
    [h[0] / h[2], h[1] / h[2]]
}

fn bounds(ctrl: &[H2]) -> [f64; 4] {
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for h in ctrl {
        let [x, y] = cartesian(h);
        b = [b[0].min(x), b[1].min(y), b[2].max(x), b[3].max(y)];
    }
    b
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// Whether every control point lies within `tol` of the chord, so that the
/// piece does by the convex hull property.
fn flat(ctrl: &[H2], tol: f64) -> bool {
    let (p0, p1) = (cartesian(&ctrl[0]), cartesian(&ctrl[ctrl.len() - 1]));
    let chord = sub2(p1, p0);
    let length = dot2(chord, chord).sqrt();
    ctrl[1..ctrl.len() - 1].iter().all(|h| {
        let q = sub2(cartesian(h), p0);
        let off = if length > 0.0 {
            cross2(chord, q).abs() / length
        } else {
            dot2(q, q).sqrt()
        };
        off <= tol
    })
}

/// How two flat pieces' chords meet.
enum Meet {
    /// They cross at these local parameters.
    Cross(f64, f64),
    /// They run together over this local stretch of the first.
    Run(f64, f64),
}

/// Where two flat pieces' chords meet, as local parameters.
fn chords(a: &[H2], b: &[H2], tol: f64) -> Option<Meet> {
    let (p0, p1) = (cartesian(&a[0]), cartesian(&a[a.len() - 1]));
    let (q0, q1) = (cartesian(&b[0]), cartesian(&b[b.len() - 1]));
    let (dp, dq) = (sub2(p1, p0), sub2(q1, q0));
    let (lp, lq) = (dot2(dp, dp).sqrt(), dot2(dq, dq).sqrt());
    let den = cross2(dp, dq);
    let r = sub2(q0, p0);
    if den.abs() <= tol * (lp + lq) {
        // Parallel within the floor: a crossing only where they run
        // together, and then at both ends of the common stretch.
        if lp == 0.0 || lq == 0.0 {
            return (dot2(r, r).sqrt() <= 2.0 * tol).then_some(Meet::Cross(0.0, 0.0));
        }
        if cross2(dp, r).abs() / lp > 2.0 * tol {
            return None;
        }
        let s0 = dot2(r, dp) / (lp * lp);
        let s1 = dot2(sub2(q1, p0), dp) / (lp * lp);
        let (lo, hi) = (s0.min(s1).max(0.0), s0.max(s1).min(1.0));
        return (lo <= hi).then_some(Meet::Run(lo, hi));
    }
    let s = cross2(r, dq) / den;
    let t = cross2(r, dp) / den;
    let (ea, eb) = (
        if lp > 0.0 { tol / lp } else { 1.0 },
        if lq > 0.0 { tol / lq } else { 1.0 },
    );
    if s < -ea || s > 1.0 + ea || t < -eb || t > 1.0 + eb {
        return None;
    }
    Some(Meet::Cross(s.clamp(0.0, 1.0), t.clamp(0.0, 1.0)))
}

/// What the search found on one pair of pieces, in local parameters.
#[derive(Default)]
struct Found {
    crossings: Vec<(f64, f64)>,
    runs: Vec<[f64; 2]>,
}

#[allow(clippy::too_many_arguments)]
fn search(
    a: &[H2],
    ra: [f64; 2],
    b: &[H2],
    rb: [f64; 2],
    depth: usize,
    tol: f64,
    out: &mut Found,
    work: &mut u64,
) {
    *work += 1;
    let (ba, bb) = (bounds(a), bounds(b));
    if ba[0] > bb[2] + tol || bb[0] > ba[2] + tol || ba[1] > bb[3] + tol || bb[1] > ba[3] + tol {
        return;
    }
    let fa = depth >= DEPTH || flat(a, tol);
    let fb = depth >= DEPTH || flat(b, tol);
    if fa && fb {
        match chords(a, b, tol) {
            Some(Meet::Cross(s, t)) => out
                .crossings
                .push((ra[0] + s * (ra[1] - ra[0]), rb[0] + t * (rb[1] - rb[0]))),
            Some(Meet::Run(lo, hi)) => out
                .runs
                .push([ra[0] + lo * (ra[1] - ra[0]), ra[0] + hi * (ra[1] - ra[0])]),
            None => {}
        }
        return;
    }
    let size = |x: [f64; 4]| (x[2] - x[0]).max(x[3] - x[1]);
    if !fa && (fb || size(ba) >= size(bb)) {
        let (left, right) = split(a, 0.5);
        let mid = (ra[0] + ra[1]) / 2.0;
        search(&left, [ra[0], mid], b, rb, depth + 1, tol, out, work);
        search(&right, [mid, ra[1]], b, rb, depth + 1, tol, out, work);
    } else {
        let (left, right) = split(b, 0.5);
        let mid = (rb[0] + rb[1]) / 2.0;
        search(a, ra, &left, [rb[0], mid], depth + 1, tol, out, work);
        search(a, ra, &right, [mid, rb[1]], depth + 1, tol, out, work);
    }
}

/// The unit roundoff of a parameter range's largest end, doubled, as in
/// `trace.rs`.
fn resolution(range: [f64; 2]) -> f64 {
    2.0 * f64::EPSILON * range[0].abs().max(range[1].abs()).max(1.0)
}

/// Newton's method on `P_a(t) - P_b(u) = 0` on the sheet, kept inside both
/// ranges. Returns the better of the start and the result.
fn polish(
    sheet: &Sheet,
    a: &Geo,
    ra: [f64; 2],
    b: &Geo,
    rb: [f64; 2],
    start: (f64, f64),
    work: &mut u64,
) -> Crossing {
    let gap = |t: f64, u: f64| {
        let (p, q) = (sheet.project(a.point(t)), sheet.project(b.point(u)));
        let d = sub2(p, q);
        dot2(d, d).sqrt()
    };
    let (mut t, mut u) = start;
    let first = Crossing {
        t,
        u,
        residual: gap(t, u),
    };
    let (res_a, res_b) = (resolution(ra), resolution(rb));
    for _ in 0..SOLVE_STEPS {
        *work += 1;
        let ((pa, da), (pb, db)) = (a.jet(t), b.jet(u));
        let f = sub2(sheet.project(pa), sheet.project(pb));
        let (ja, jb) = (sheet.project(da), sheet.project(db));
        // Columns J = [P_a', -P_b'].
        let det = ja[0] * -jb[1] - -jb[0] * ja[1];
        if det == 0.0 || !det.is_finite() {
            break;
        }
        let dt = (f[0] * -jb[1] - -jb[0] * f[1]) / det;
        let du = (ja[0] * f[1] - ja[1] * f[0]) / det;
        let next = ((t - dt).clamp(ra[0], ra[1]), (u - du).clamp(rb[0], rb[1]));
        let done = (next.0 - t).abs() <= res_a && (next.1 - u).abs() <= res_b;
        (t, u) = next;
        if done {
            break;
        }
    }
    let last = Crossing {
        t,
        u,
        residual: gap(t, u),
    };
    if last.residual <= first.residual {
        last
    } else {
        first
    }
}

/// Crossings of `a` over `ra` with `b` over `rb` on the sheet, within twice
/// rule T1's floor `tol`. For `same`, the stroke against itself: distinct
/// pieces only, and both parameters are cuts. Where the two run together on
/// the sheet, only the ends of each common stretch count, with their own
/// parameter as `u` too.
#[allow(clippy::too_many_arguments)]
pub(super) fn crossings(
    sheet: &Sheet,
    a: &Geo,
    ra: [f64; 2],
    b: &Geo,
    rb: [f64; 2],
    same: bool,
    tol: f64,
    work: &mut u64,
) -> Vec<Crossing> {
    let pa = a.pieces(ra[0], ra[1]);
    let pb = b.pieces(rb[0], rb[1]);
    let mut starts = Vec::new();
    let mut runs: Vec<[f64; 2]> = Vec::new();
    for (i, x) in pa.iter().enumerate() {
        let xa = project(sheet, &x.ctrl);
        for (j, y) in pb.iter().enumerate() {
            if same && j <= i {
                continue;
            }
            let yb = project(sheet, &y.ctrl);
            let mut local = Found::default();
            search(&xa, [0.0, 1.0], &yb, [0.0, 1.0], 0, tol, &mut local, work);
            starts.extend(
                local
                    .crossings
                    .into_iter()
                    .map(|(s, t)| (x.map.native(s), y.map.native(t))),
            );
            if !same {
                runs.extend(
                    local
                        .runs
                        .into_iter()
                        .map(|[lo, hi]| [x.map.native(lo), x.map.native(hi)]),
                );
            }
        }
    }
    // Leaves of one common stretch touch end to end; merge them.
    runs.sort_by(|x, y| x[0].total_cmp(&y[0]));
    let slack = resolution(ra);
    let mut merged: Vec<[f64; 2]> = Vec::new();
    for run in runs {
        match merged.last_mut() {
            Some(last) if run[0] <= last[1] + slack => last[1] = last[1].max(run[1]),
            _ => merged.push(run),
        }
    }
    let within = |t: f64| merged.iter().any(|r| t > r[0] + slack && t < r[1] - slack);
    let mut found: Vec<Crossing> = Vec::new();
    for run in &merged {
        for t in *run {
            found.push(Crossing {
                t,
                u: t,
                residual: tol,
            });
        }
    }
    for start in starts {
        if within(start.0) {
            continue;
        }
        let c = polish(sheet, a, ra, b, rb, start, work);
        if c.residual > 2.0 * tol {
            continue;
        }
        // Against itself, neighbouring pieces share their joint: one point
        // in space, not a crossing of two points at different depths.
        if same && distance(a.point(c.t), a.point(c.u)) <= tol {
            continue;
        }
        let p = a.point(c.t);
        match found.iter_mut().find(|f| distance(a.point(f.t), p) <= tol) {
            Some(f) if f.residual > c.residual => *f = c,
            Some(_) => {}
            None => found.push(c),
        }
    }
    found.sort_by(|x, y| x.t.total_cmp(&y.t));
    found
}

fn binomial(n: usize, k: usize) -> f64 {
    let mut value = 1.0;
    for i in 0..k {
        value = value * (n - i) as f64 / (i + 1) as f64;
    }
    value
}

/// The product of two Bernstein polynomials (Farouki and Rajan, "Algorithms
/// for polynomials in Bernstein form", CAGD 5, 1988, section 2.4).
fn times(a: &[f64], b: &[f64]) -> Vec<f64> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    (0..=m + n)
        .map(|k| {
            let lo = k.saturating_sub(n);
            (lo..=m.min(k))
                .map(|i| binomial(m, i) * binomial(n, k - i) / binomial(m + n, k) * a[i] * b[k - i])
                .sum()
        })
        .collect()
}

/// The derivative of a Bernstein polynomial (Farouki and Rajan, section 2.5).
fn derivative(a: &[f64]) -> Vec<f64> {
    let n = a.len() - 1;
    if n == 0 {
        return vec![0.0];
    }
    (0..n).map(|i| n as f64 * (a[i + 1] - a[i])).collect()
}

/// Local minima of a polynomial on `[0, 1]`: where its derivative `g`
/// crosses from negative to positive, isolated by halving while the
/// Bernstein coefficients change sign (Descartes' rule of signs in Bernstein
/// form, Farouki and Rajan, section 3.2). Near a root the halves'
/// coefficients shrink with the interval, so only the caller judges binary64
/// noise, once, on the whole piece.
fn minima(g: &[f64], range: [f64; 2], depth: usize, out: &mut Vec<f64>, work: &mut u64) {
    *work += 1;
    if g.iter().all(|&c| c >= 0.0) || g.iter().all(|&c| c <= 0.0) {
        return;
    }
    if depth >= DEPTH {
        if g[0] < 0.0 && g[g.len() - 1] > 0.0 {
            out.push((range[0] + range[1]) / 2.0);
        }
        return;
    }
    let as_rows: Vec<[f64; 1]> = g.iter().map(|&c| [c]).collect();
    let (left, right) = split(&as_rows, 0.5);
    let mid = (range[0] + range[1]) / 2.0;
    let flat = |v: Vec<[f64; 1]>| v.into_iter().map(|c| c[0]).collect::<Vec<_>>();
    let (left, right) = (flat(left), flat(right));
    // A root on the halving point itself is in neither half's interior.
    let n = left.len();
    if n >= 2 && !above(left[n - 1].abs(), 0.0) && left[n - 2] < 0.0 && right[1] > 0.0 {
        out.push(mid);
    }
    minima(&left, [range[0], mid], depth + 1, out, work);
    minima(&right, [mid, range[1]], depth + 1, out, work);
}

/// Where a spline stroke's sheet speed has a local minimum at which the
/// sheet direction reverses: the cusps of its projection, where the sheet
/// derivative vanishes (Koenderink, "What does the occluding contour tell
/// us about solid shape?", Perception 13, 1984: the apparent contour ends
/// or cusps only there). A fitted carrier passes near such a point rather
/// than through it, so its speed only dips; the reversal ([`cusp`]) tells a
/// cusp from an ordinary slow stretch.
pub(super) fn slow_points(
    sheet: &Sheet,
    geo: &Geo,
    range: [f64; 2],
    closed: bool,
    work: &mut u64,
) -> Vec<f64> {
    if !matches!(geo, Geo::Spline { .. }) {
        return Vec::new();
    }
    let pieces = geo.pieces(range[0], range[1]);
    let bounds: Vec<[f64; 2]> = pieces
        .iter()
        .map(|p| [p.map.native(0.0), p.map.native(1.0)])
        .collect();
    let mut out = Vec::new();
    // The slope of the squared sheet speed at the end of the previous piece.
    let mut falling = false;
    for (k, piece) in pieces.iter().enumerate() {
        let ctrl = project(sheet, &piece.ctrl);
        let column = |k: usize| ctrl.iter().map(|h| h[k]).collect::<Vec<_>>();
        let (x, y, w) = (column(0), column(1), column(2));
        let (dx, dy, dw) = (derivative(&x), derivative(&y), derivative(&w));
        // The numerator of the sheet derivative of (x / w, y / w).
        let minus =
            |a: Vec<f64>, b: Vec<f64>| a.iter().zip(&b).map(|(p, q)| p - q).collect::<Vec<_>>();
        let nx = minus(times(&dx, &w), times(&x, &dw));
        let ny = minus(times(&dy, &w), times(&y, &dw));
        let speed: Vec<f64> = times(&nx, &nx)
            .iter()
            .zip(times(&ny, &ny))
            .map(|(p, q)| p + q)
            .collect();
        let scale = speed.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
        let g = derivative(&speed);
        let noise = 4.0 * g.len() as f64 * f64::EPSILON * scale;
        let mut local = Vec::new();
        // A slowest point on the joint with the previous piece, where the
        // speed's slope may jump, is in neither piece's interior.
        if k > 0 && falling && g[0] > noise {
            local.push(0.0);
        }
        falling = g[g.len() - 1] < -noise;
        // Where every coefficient is binary64 noise the speed is flat: no
        // isolated minimum.
        if g.iter().any(|c| c.abs() > noise) {
            minima(&g, [0.0, 1.0], 0, &mut local, work);
        }
        if local.is_empty() {
            continue;
        }
        // The window reaches into both neighbouring pieces, so a cusp near a
        // joint is judged on its whole turn.
        let window = [
            bounds[k.saturating_sub(1)][0],
            bounds[(k + 1).min(bounds.len() - 1)][1],
        ];
        for s in local {
            let t = piece.map.native(s);
            if cusp(sheet, geo, t, window, work) {
                out.push(t);
            }
        }
    }
    // A fitted carrier may also turn back at a joint where only its point is
    // shared (a triple knot, or a closed carrier's own ends): the one-sided
    // sheet directions there point apart.
    let ends: Vec<[[f64; 2]; 2]> = pieces
        .iter()
        .map(|piece| {
            let ctrl = project(sheet, &piece.ctrl);
            [ends_direction(&ctrl, 0.0), ends_direction(&ctrl, 1.0)]
        })
        .collect();
    let joints = if closed {
        pieces.len()
    } else {
        pieces.len().saturating_sub(1)
    };
    for k in 0..joints {
        *work += 1;
        let next = (k + 1) % pieces.len();
        if dot2(ends[k][1], ends[next][0]) < 0.0 {
            out.push(if next == 0 { range[0] } else { bounds[next][0] });
        }
    }
    out
}

/// The sheet direction of homogeneous sheet controls at local `s`.
fn ends_direction(ctrl: &[H2], s: f64) -> [f64; 2] {
    let rows: Vec<[f64; 4]> = ctrl.iter().map(|h| [h[0], h[1], 0.0, h[2]]).collect();
    let (h, dh) = super::geo::casteljau(&rows, s);
    [
        (dh[0] * h[3] - h[0] * dh[3]) / (h[3] * h[3]),
        (dh[1] * h[3] - h[1] * dh[3]) / (h[3] * h[3]),
    ]
}

/// Whether the slow point `t` is a cusp of the projection: whether the
/// sheet direction has reversed by where the sheet speed has doubled on
/// either side. Near a cusp the sheet derivative is `m e + c (s - t) f`
/// (Whitney's cusp, with `m` zero on the true curve and small on a fitted
/// one); the speed doubles at `|s - t| = sqrt(3) m / c`, where each side has
/// turned by sixty degrees from the slowest direction, so the two
/// directions there are opposite to within sixty degrees. An ordinary slow
/// stretch keeps its direction.
fn cusp(sheet: &Sheet, geo: &Geo, t: f64, window: [f64; 2], work: &mut u64) -> bool {
    let speed = |x: f64| {
        let v = sheet.project(geo.jet(x).1);
        dot2(v, v).sqrt()
    };
    let slowest = speed(t);
    // The nearest parameter on one side where the speed reaches twice the
    // slowest, by halving the distance to the window's end.
    let mut reach = |end: f64| {
        let mut far = end;
        for _ in 0..DEPTH {
            *work += 1;
            let near = t + (far - t) / 2.0;
            if speed(near) >= 2.0 * slowest && above((near - t).abs(), 0.0) {
                far = near;
            } else {
                break;
            }
        }
        far
    };
    let (before, after) = (reach(window[0]), reach(window[1]));
    let (u, v) = (
        sheet.project(geo.jet(before).1),
        sheet.project(geo.jet(after).1),
    );
    dot2(u, v) < 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sheet() -> Sheet {
        Sheet {
            e1: [1.0, 0.0, 0.0],
            e2: [0.0, 1.0, 0.0],
        }
    }

    #[test]
    fn a_line_crosses_a_circle_twice_where_the_closed_form_says() {
        let line = Geo::Line {
            o: [-3.0, 0.5, 0.0],
            v: [1.0, 0.0, 0.0],
        };
        let circle = Geo::Conic {
            c: [0.0, 0.0, 5.0],
            u: [1.0, 0.0, 0.0],
            w: [0.0, 1.0, 0.0],
        };
        let mut work = 0;
        let found = crossings(
            &sheet(),
            &line,
            [0.0, 6.0],
            &circle,
            [0.0, std::f64::consts::TAU],
            false,
            1e-7,
            &mut work,
        );
        assert_eq!(found.len(), 2);
        let x = 0.75_f64.sqrt();
        assert!((found[0].t - (3.0 - x)).abs() < 1e-14 && (found[1].t - (3.0 + x)).abs() < 1e-14);
        for c in &found {
            assert!(c.residual < 1e-14);
        }
    }

    #[test]
    fn a_cusped_cubic_has_its_cusp_among_its_slow_points() {
        // (t^2, t^3) about t = 0, over [-1, 1]: Bezier controls of x = t^2
        // and y = t^3 on s in [0, 1] with t = 2 s - 1.
        let ctrl = vec![
            [1.0, -1.0, 0.0, 1.0],
            [-1.0 / 3.0, 1.0, 0.0, 1.0],
            [-1.0 / 3.0, -1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
        ];
        let geo = Geo::Spline {
            spans: vec![super::super::geo::Span {
                lo: 0.0,
                hi: 1.0,
                ctrl,
            }],
        };
        let mut work = 0;
        let slow = slow_points(&sheet(), &geo, [0.0, 1.0], false, &mut work);
        assert!(slow.iter().any(|&s| (s - 0.5).abs() < 1e-9), "{slow:?}");
    }
}
