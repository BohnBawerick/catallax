//! Fitted carriers for contour generators that have no exact curve.
//!
//! The caller knows the true generator pointwise, with its first derivative,
//! in the generator's own parameter. The fit is a chain of cubic Hermite spans
//! in that parameter, the form SMLib's surface tracer gives a traced curve
//! (`SmSurfaceTracer`, spans through `SmHermiteCurve`, a span split where its
//! accuracy test fails). Written as a B-spline, the chain has a double knot at
//! every joint, so it is C1 (Piegl and Tiller, "The NURBS Book", 2nd ed.,
//! section 5.2: a knot of multiplicity `k` leaves a degree-`p` curve
//! `C^(p-k)`).
//!
//! Each span's deviation from the true generator is measured as decision Q6
//! (Checkpoint 2 decisions, section 4.2, setup A+) measures an edge gap: the
//! checker's own sample count for a cubic piece at Chebyshev points of the
//! first kind, both ends, the coefficient tail as the estimate of what the
//! samples leave unresolved, and Brent's method on every local maximum. A
//! span passes when the polished maximum, the unresolved estimate and rule
//! T1's binary64 floor together stay within rule T2's fit target; otherwise
//! it is halved. The tail can only make a span fail, so it never ends a
//! search early. The result keeps that sum, the reported deviation.
//!
//! All arithmetic is IEEE binary64 basic operations and `libm`, in one
//! thread, so a fit is the same bits on every machine and at any thread
//! count.

use catallax_numbers::{Decider, Refusal};
use catallax_topology::numerical::gap::PieceShape;

/// A point of the true generator and its first derivative in the fit
/// parameter, in the caller's local frame.
pub(super) type Jet = ([f64; 3], [f64; 3]);

/// A fitted chain of cubic Hermite spans.
#[derive(Clone, Debug)]
pub(super) struct Fitted {
    /// Span joints in the fit parameter, ascending, both ends included.
    pub breaks: Vec<f64>,
    /// B-spline controls in the local frame: the first point, the two inner
    /// Bezier controls of every span, and the last point.
    pub controls: Vec<[f64; 3]>,
    /// The largest polished deviation plus its unresolved estimate, over all
    /// spans, plus rule T1's floor.
    pub deviation_mm: f64,
}

/// Source: Mason and Handscomb, "Chebyshev Polynomials", 2003, section 6.3:
/// the `count` first-kind points in `(0, 1)`, ascending,
/// `(1 - cos((2 j + 1) pi / (2 count))) / 2`. The checker of 3.0A samples a
/// piece at the same points.
fn chebyshev_shares(count: usize) -> Vec<f64> {
    (0..count)
        .map(|j| {
            let angle = std::f64::consts::PI * (2 * j + 1) as f64 / (2 * count) as f64;
            (1.0 - libm::cos(angle)) / 2.0
        })
        .collect()
}

/// Source: Orszag, "On the elimination of aliasing in finite-difference
/// schemes by filtering high-wavenumber components", J. Atmos. Sci. 28,
/// 1971: aliased content lands in the top third of the retained modes, so
/// that third is read as the tail. The 3.0A checker reads the same third.
const TAIL_SHARE_DENOMINATOR: usize = 3;

/// Source: Trefethen, "Approximation Theory and Approximation Practice",
/// 2013, chapter 4: a Chebyshev interpolant is off by at most twice the sum
/// of the magnitudes of the coefficients it does not carry; the last
/// retained coefficients estimate those (Boyd, "Chebyshev and Fourier
/// Spectral Methods", 2nd ed., section 2.12).
const INTERPOLATION_FACTOR: f64 = 2.0;

/// Source: the triangle inequality. If the interpolant is within `e` of the
/// residual, the true maximum is at most `e` above the interpolant's and the
/// polished maximum at most `e` below it, so the search can be short by two
/// interpolation errors.
const SEARCH_FACTOR: f64 = 2.0;

/// The sum of the Euclidean lengths of the top third of the Chebyshev
/// coefficients of a vector function sampled at first-kind points, from the
/// discrete orthogonality `c_k = (2 / n) sum_j f(x_j) T_k(x_j)` (Mason and
/// Handscomb, section 6.3). Samples ascend in the share, which only flips
/// the sign of odd terms.
fn chebyshev_tail(samples: &[[f64; 3]]) -> f64 {
    let n = samples.len();
    let table: Vec<f64> = (0..4 * n)
        .map(|m| libm::cos(std::f64::consts::PI * m as f64 / (2 * n) as f64))
        .collect();
    let first = n - n / TAIL_SHARE_DENOMINATOR;
    (first..n)
        .map(|k| {
            let mut coefficient = [0.0; 3];
            for (j, sample) in samples.iter().enumerate() {
                let weight = table[(k * (2 * j + 1)) % (4 * n)];
                for axis in 0..3 {
                    coefficient[axis] += sample[axis] * weight;
                }
            }
            length(coefficient) * 2.0 / n as f64
        })
        .sum()
}

/// Source: Brent, "Algorithms for Minimization without Derivatives", 1973,
/// chapter 5: the golden section ratio `(3 - sqrt(5)) / 2`.
const GOLDEN: f64 = 0.381_966_011_250_105_1;

/// Source: judgement, as in the 3.0A checker. Brent's method needs about ten
/// to forty evaluations on a smooth bracket; one hundred steps only ends a
/// search that stalls.
const BRENT_STEPS: u32 = 100;

/// The maximum of `f` on `[0, 1]` by Brent's procedure `localmin` applied to
/// `-f` (Brent 1973, chapter 5), with the relative tolerance the square root
/// of the machine epsilon, the smallest Brent recommends. Returns the best
/// value and the evaluations spent.
// Brent's procedure tests whether two stored abscissae are the same point,
// so these comparisons are exact by design.
#[allow(clippy::float_cmp)]
fn brent_maximum(mut f: impl FnMut(f64) -> f64) -> (f64, u64) {
    let tolerance = f64::EPSILON.sqrt();
    let (mut a, mut b) = (0.0_f64, 1.0_f64);
    let mut x = a + GOLDEN * (b - a);
    let (mut w, mut v) = (x, x);
    let mut fx = -f(x);
    let (mut fw, mut fv) = (fx, fx);
    let (mut d, mut e) = (0.0_f64, 0.0_f64);
    let mut evaluations = 1;
    for _ in 0..BRENT_STEPS {
        let middle = (a + b) / 2.0;
        let tol = tolerance * x.abs() + tolerance / 2.0;
        if (x - middle).abs() <= 2.0 * tol - (b - a) / 2.0 {
            break;
        }
        let mut parabolic = false;
        if e.abs() > tol {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            if p.abs() < (q * e / 2.0).abs() && p > q * (a - x) && p < q * (b - x) {
                e = d;
                d = p / q;
                let u = x + d;
                if u - a < 2.0 * tol || b - u < 2.0 * tol {
                    d = if x < middle { tol } else { -tol };
                }
                parabolic = true;
            }
        }
        if !parabolic {
            e = if x < middle { b - x } else { a - x };
            d = GOLDEN * e;
        }
        let u = if d.abs() >= tol {
            x + d
        } else if d > 0.0 {
            x + tol
        } else {
            x - tol
        };
        let fu = -f(u);
        evaluations += 1;
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            (v, fv) = (w, fw);
            (w, fw) = (x, fx);
            (x, fx) = (u, fu);
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                (v, fv) = (w, fw);
                (w, fw) = (u, fu);
            } else if fu <= fv || v == x || v == w {
                (v, fv) = (u, fu);
            }
        }
    }
    (-fx, evaluations)
}

fn length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn axpy(p: [f64; 3], s: f64, d: [f64; 3]) -> [f64; 3] {
    [p[0] + s * d[0], p[1] + s * d[1], p[2] + s * d[2]]
}

/// The cubic Bezier form of a Hermite span over `[a, b]`: `b1 = p0 + h/3 d0`
/// and `b2 = p1 - h/3 d1` with `h = b - a` (Piegl and Tiller, section 1.4,
/// the derivative at the ends of a Bezier curve).
fn bezier(from: &(f64, Jet), to: &(f64, Jet)) -> [[f64; 3]; 4] {
    let third = (to.0 - from.0) / 3.0;
    [
        from.1.0,
        axpy(from.1.0, third, from.1.1),
        axpy(to.1.0, -third, to.1.1),
        to.1.0,
    ]
}

/// De Casteljau's algorithm (Piegl and Tiller, section 1.4).
fn bezier_point(b: &[[f64; 3]; 4], t: f64) -> [f64; 3] {
    let mut level = *b;
    for round in 1..4 {
        for i in 0..4 - round {
            for axis in 0..3 {
                level[i][axis] = level[i][axis] + t * (level[i + 1][axis] - level[i][axis]);
            }
        }
    }
    level[0]
}

/// Measures one span: the polished largest deviation, the unresolved
/// estimate, and the work spent.
fn measure(
    generator: &impl Fn(f64) -> Jet,
    from: &(f64, Jet),
    to: &(f64, Jet),
    controls: &[[f64; 3]; 4],
) -> (f64, f64, u64) {
    let span = to.0 - from.0;
    let residual = |t: f64| -> [f64; 3] {
        let fitted = bezier_point(controls, t);
        let (exact, _) = generator(from.0 + t * span);
        [
            fitted[0] - exact[0],
            fitted[1] - exact[1],
            fitted[2] - exact[2],
        ]
    };
    // The 3.0A checker's count for a piece whose 3D copy is a cubic.
    let shares = chebyshev_shares(
        PieceShape {
            curve: 3,
            face: 0,
            surface: 0,
        }
        .points(),
    );
    let interior: Vec<[f64; 3]> = shares.iter().map(|&t| residual(t)).collect();
    let pairs = (shares.len() * (shares.len() / TAIL_SHARE_DENOMINATOR)) as u64;
    let unresolved = SEARCH_FACTOR * INTERPOLATION_FACTOR * chebyshev_tail(&interior);
    let mut ts = vec![0.0];
    ts.extend(&shares);
    ts.push(1.0);
    let mut values = vec![length(residual(0.0))];
    values.extend(interior.iter().map(|&r| length(r)));
    values.push(length(residual(1.0)));
    let mut work = values.len() as u64 + pairs;
    let mut best = values.iter().copied().fold(0.0_f64, f64::max);
    // Polish every local maximum; a run of equal values is one candidate.
    let n = values.len();
    for j in 0..n {
        let rises = j == 0 || values[j] > values[j - 1];
        let holds = j + 1 == n || values[j] >= values[j + 1];
        if !(rises && holds) {
            continue;
        }
        let lo = ts[j.saturating_sub(1)];
        let hi = ts[(j + 1).min(n - 1)];
        if hi <= lo {
            continue;
        }
        let (peak, spent) = brent_maximum(|x| length(residual(lo + x * (hi - lo))));
        work += spent;
        best = best.max(peak);
    }
    (best, unresolved, work)
}

/// Fits the generator from `from` to `to`, both given with their jets, to
/// `target_mm`, adding `floor_mm` (rule T1) to every measured deviation.
/// Refuses with `app.contour-fit` when a failing span cannot be halved in
/// binary64, and with `app.contour-fit-work` when the fit would spend more
/// than `cap` units.
pub(super) fn fit(
    generator: impl Fn(f64) -> Jet,
    from: (f64, Jet),
    to: (f64, Jet),
    target_mm: f64,
    floor_mm: f64,
    cap: u64,
    dec: &mut Decider,
) -> Result<Fitted, Refusal> {
    let mut breaks = vec![from.0];
    let mut controls = vec![from.1.0];
    let mut deviation = 0.0_f64;
    let mut work = 0_u64;
    // Spans narrower than binary64 resolves at the domain's magnitude are not
    // split again: the unit roundoff of IEEE 754 binary64 is half of
    // `f64::EPSILON`, relative to the largest parameter.
    let resolution = from.0.abs().max(to.0.abs()) * f64::EPSILON;
    let mut left = from;
    // Right ends still to reach, nearest last: the leftmost failing span is
    // always split first, so the order of work is fixed.
    let mut pending = vec![to];
    while let Some(right) = pending.last().copied() {
        dec.check_deadline()?;
        let bezier = bezier(&left, &right);
        let (peak, unresolved, spent) = measure(&generator, &left, &right, &bezier);
        work += spent;
        dec.charge_work(spent)?;
        if work > cap {
            return Err(Refusal::new(
                "app.contour-fit-work",
                "a fitted contour generator exceeds its work allowance",
            ));
        }
        if peak + unresolved + floor_mm <= target_mm {
            deviation = deviation.max(peak + unresolved);
            controls.push(bezier[1]);
            controls.push(bezier[2]);
            breaks.push(right.0);
            left = right;
            pending.pop();
            continue;
        }
        let middle = left.0 + (right.0 - left.0) / 2.0;
        if right.0 - left.0 <= 2.0 * resolution || !(left.0 < middle && middle < right.0) {
            return Err(Refusal::new(
                "app.contour-fit",
                "a fitted contour generator cannot meet rule T2's fit target in binary64",
            ));
        }
        work += 1;
        dec.charge_work(1)?;
        pending.push((middle, generator(middle)));
    }
    controls.push(left.1.0);
    Ok(Fitted {
        breaks,
        controls,
        deviation_mm: (deviation + floor_mm).next_up(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A helix `(cos u, sin u, u / 4)`: smooth, so the fit converges, and not
    /// a polynomial, so the first span fails.
    fn helix(u: f64) -> Jet {
        (
            [libm::cos(u), libm::sin(u), u / 4.0],
            [-libm::sin(u), libm::cos(u), 0.25],
        )
    }

    #[test]
    fn a_helix_fits_within_the_target_and_reports_a_bounded_deviation() {
        let (target, floor) = (1e-6, 1e-7);
        let fitted = fit(
            helix,
            (0.0, helix(0.0)),
            (3.0, helix(3.0)),
            target,
            floor,
            1_000_000,
            &mut Decider::default(),
        )
        .unwrap();
        let spans = fitted.breaks.len() - 1;
        assert!(spans > 1);
        assert_eq!(fitted.controls.len(), 2 * spans + 2);
        assert!(fitted.deviation_mm <= target && fitted.deviation_mm >= floor);
        // A dense independent check: no point of the chain is further from
        // the helix than the reported deviation.
        for i in 0..spans {
            let b = [
                if i == 0 {
                    fitted.controls[0]
                } else {
                    helix(fitted.breaks[i]).0
                },
                fitted.controls[1 + 2 * i],
                fitted.controls[2 + 2 * i],
                helix(fitted.breaks[i + 1]).0,
            ];
            for k in 0..=1000 {
                let t = k as f64 / 1000.0;
                let u = fitted.breaks[i] + t * (fitted.breaks[i + 1] - fitted.breaks[i]);
                let p = bezier_point(&b, t);
                let q = helix(u).0;
                let gap = length([p[0] - q[0], p[1] - q[1], p[2] - q[2]]);
                assert!(gap <= fitted.deviation_mm, "{gap} at {u}");
            }
        }
    }

    #[test]
    fn an_unreachable_target_and_a_starved_cap_refuse_by_name() {
        let unreachable = fit(
            helix,
            (0.0, helix(0.0)),
            (3.0, helix(3.0)),
            1e-30,
            0.0,
            u64::MAX,
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(unreachable.code, "app.contour-fit");
        let starved = fit(
            helix,
            (0.0, helix(0.0)),
            (3.0, helix(3.0)),
            1e-6,
            1e-7,
            100,
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(starved.code, "app.contour-fit-work");
    }

    #[test]
    fn brent_finds_an_interior_maximum() {
        let (peak, evaluations) = brent_maximum(|x| 1.0 - (x - 0.3) * (x - 0.3));
        assert!((peak - 1.0).abs() < 1e-15);
        assert!(evaluations < 40);
    }
}
