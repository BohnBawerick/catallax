//! Binary64 building blocks of the gap checker. Each names its source.
//!
//! Only IEEE basic operations and square roots are used, so every machine
//! and every thread count returns the same bits. No library sine or cosine
//! is called.

#![allow(clippy::float_cmp)] // Identities of stored values, as Brent's method states them.

/// Source: IEEE 754 binary64 unit roundoff, 2^-53 (Higham, "Accuracy and
/// Stability of Numerical Algorithms", 2nd ed., section 2.1).
pub(super) const UNIT_ROUNDOFF: f64 = f64::EPSILON / 2.0;

/// Source: fdlibm `e_rem_pio2.c`, `pio2_1` (1.57079632673412561417, bits
/// 0x3FF921FB 0x54400000): the first 33 bits of pi / 2, so
/// its product with a quadrant count below 2^20 is exact (Cody and Waite,
/// "Software Manual for the Elementary Functions", 1980, chapter 8).
const HALF_PI_HEAD: f64 = f64::from_bits(0x3FF9_21FB_5440_0000);
/// Source: fdlibm `e_rem_pio2.c`, `pio2_1t` (6.07710050650619224932e-11,
/// bits 0x3DD0B461 0x1A626331): pi / 2 minus `HALF_PI_HEAD`.
const HALF_PI_TAIL: f64 = f64::from_bits(0x3DD0_B461_1A62_6331);
/// Source: the Taylor remainder of sine and cosine on [-pi/4, pi/4]. After
/// this many terms it is below (pi/4)^19 / 19! = 8e-20, under the unit
/// roundoff.
const TAYLOR_TERMS: u32 = 9;
/// Source: an operation count of [`sin_cos`]. The reduction rounds twice and
/// each Taylor sum accumulates under two units, so four units bound a value.
pub(super) const TRIG_ROUNDINGS: f64 = 4.0;

/// Source: a turn has four quadrants of pi / 2.
const QUADRANTS: f64 = 4.0;

/// Sine and cosine from basic operations only: Cody-Waite reduction to
/// [-pi/4, pi/4], then Taylor sums.
pub(super) fn sin_cos(x: f64) -> (f64, f64) {
    let quadrant = (x / (HALF_PI_HEAD + HALF_PI_TAIL)).round();
    let r = (x - quadrant * HALF_PI_HEAD) - quadrant * HALF_PI_TAIL;
    let z = r * r;
    let (mut sine, mut cosine) = (r, 1.0);
    let (mut odd, mut even) = (r, 1.0);
    for n in 1..=TAYLOR_TERMS {
        let k = f64::from(2 * n);
        even *= -z / (k * (k - 1.0));
        odd *= -z / (k * (k + 1.0));
        cosine += even;
        sine += odd;
    }
    match (quadrant.rem_euclid(QUADRANTS)) as u8 {
        0 => (sine, cosine),
        1 => (cosine, -sine),
        2 => (-sine, -cosine),
        _ => (-cosine, sine),
    }
}

/// The knot span holding `t`: the last `i` in `degree..count` with
/// `knots[i] <= t`. A parameter outside the domain gets the end span, whose
/// polynomial then continues past the border.
/// Source: Piegl and Tiller, "The NURBS Book", 2nd ed., algorithm A2.1.
pub(super) fn find_span(knots: &[f64], degree: usize, count: usize, t: f64) -> usize {
    let (mut lo, mut hi) = (degree, count);
    while hi - lo > 1 {
        let middle = lo + (hi - lo) / 2;
        if knots[middle] <= t {
            lo = middle;
        } else {
            hi = middle;
        }
    }
    // A clamped end repeats its knot; step back to the last nonempty span.
    while lo > degree && knots[lo] >= knots[count] {
        lo -= 1;
    }
    lo
}

/// The `degree + 1` basis functions that are nonzero on a span.
/// Source: Piegl and Tiller, algorithm A2.2 (Cox-de Boor recurrence).
///
/// Each knot is a binary64 value in `knots` plus its remainder in `low`, the
/// exact knot minus that value. A difference `t - knot` is taken as
/// `(t - value) - remainder`. The first step is exact for a knot within a
/// factor of two of `t` (Sterbenz' lemma; Higham, theorem 2.5), so a
/// difference keeps its relative accuracy however narrow the span is. With
/// rounded knots alone, a span 5e-7 wide at parameter 2 is known to one part
/// in 1e9 only.
pub(super) fn basis_functions(
    knots: &[f64],
    low: &[f64],
    degree: usize,
    span: usize,
    t: f64,
    values: &mut Vec<f64>,
    scratch: &mut Vec<f64>,
) {
    values.clear();
    values.resize(degree + 1, 0.0);
    scratch.clear();
    scratch.resize(2 * (degree + 1), 0.0);
    let (left, right) = scratch.split_at_mut(degree + 1);
    values[0] = 1.0;
    for j in 1..=degree {
        left[j] = (t - knots[span + 1 - j]) - low[span + 1 - j];
        right[j] = (knots[span + j] - t) + low[span + j];
        let mut saved = 0.0;
        for r in 0..j {
            let share = values[r] / (right[r + 1] + left[j - r]);
            values[r] = saved + right[r + 1] * share;
            saved = left[j - r] * share;
        }
        values[j] = saved;
    }
}

/// Rounding errors that one basis-weighted sum of the given degree can
/// accumulate, in units of the roundoff times the sum of magnitudes.
/// Source: an operation count of [`basis_functions`] and the sum over its
/// values. Each of `degree` levels rounds a difference, a quotient and a
/// product-sum (three), the weighted sum adds one product and one addition
/// per term (two), and the stored control and the final division each round
/// once. Sums of nonnegative terms bound the result as in Higham, lemma 3.1,
/// and corner-cutting evaluation has this form of bound (Mainar and Pena,
/// "Error analysis of corner cutting algorithms", Numerical Algorithms 22,
/// 1999).
pub(super) fn roundings(degree: usize) -> f64 {
    (3 * degree + 2 * (degree + 1) + 2) as f64
}

/// The blossom of one polynomial span at the given arguments: de Boor's
/// algorithm with a different argument at each level.
/// Source: Ramshaw, "Blossoming: a connect-the-dots approach to splines",
/// DEC SRC report 19, 1987.
pub(super) fn blossom<const D: usize>(
    knots: &[f64],
    degree: usize,
    controls: &[[f64; D]],
    span: usize,
    arguments: &[f64],
) -> [f64; D] {
    let mut d: Vec<[f64; D]> = controls[span - degree..=span].to_vec();
    for (level, t) in arguments.iter().enumerate() {
        let r = level + 1;
        for j in (r..=degree).rev() {
            let i = span - degree + j;
            let share = (t - knots[i]) / (knots[i + degree - r + 1] - knots[i]);
            for axis in 0..D {
                d[j][axis] = (1.0 - share) * d[j - 1][axis] + share * d[j][axis];
            }
        }
    }
    d[degree]
}

/// Bezier coefficients of one span restricted to `[from, to]`: the blossom at
/// `degree - i` copies of `from` and `i` copies of `to`.
pub(super) fn bezier_piece<const D: usize>(
    knots: &[f64],
    degree: usize,
    controls: &[[f64; D]],
    span: usize,
    from: f64,
    to: f64,
) -> Vec<[f64; D]> {
    (0..=degree)
        .map(|i| {
            let arguments: Vec<f64> = (0..degree)
                .map(|k| if k < degree - i { from } else { to })
                .collect();
            blossom(knots, degree, controls, span, &arguments)
        })
        .collect()
}

/// Source: judgement. Forty halvings locate a crossing to 2^-40 of its
/// piece, about 1e-12, far inside the spacing of any sample plan.
const ROOT_DEPTH: u32 = 40;

/// Shares in `(0, 1)` where a Bernstein polynomial is zero, by recursive de
/// Casteljau halving. A piece whose coefficients share one strict sign has
/// no root (variation diminishing property).
/// Source: Lane and Riesenfeld, "Bounds on a polynomial", BIT 21, 1981.
/// `zero` is the magnitude below which a value counts as zero. A root at an
/// end of the piece is the piece's own border and is not reported, and a
/// polynomial that is zero throughout has no crossings.
///
/// `allowance` is the number of halvings the caller can still pay for. One
/// is taken from it before each halving. Returns the halvings done, or
/// `None` when the allowance ran out first; the search stops there.
pub(super) fn sign_changes(
    coefficients: &[f64],
    zero: f64,
    roots: &mut Vec<f64>,
    allowance: &mut u64,
) -> Option<u64> {
    let before = *allowance;
    // A polynomial of degree n has at most n roots.
    let most = roots.len() + coefficients.len() - 1;
    halve(coefficients, zero, [0.0, 1.0], 0, most, roots, allowance).then(|| before - *allowance)
}

fn halve(
    coefficients: &[f64],
    zero: f64,
    cell: [f64; 2],
    depth: u32,
    most: usize,
    roots: &mut Vec<f64>,
    allowance: &mut u64,
) -> bool {
    // A zero at an end of the cell was reported by the caller, or is the
    // border of the whole piece. Divide it out: b(s) = s q(s) or (1 - s) q(s).
    let mut c = coefficients.to_vec();
    while c.len() > 1 && c[0].abs() <= zero {
        let n = (c.len() - 1) as f64;
        c = (1..c.len()).map(|i| c[i] * n / i as f64).collect();
    }
    while c.len() > 1 && c[c.len() - 1].abs() <= zero {
        let n = (c.len() - 1) as f64;
        c = (0..c.len() - 1)
            .map(|i| c[i] * n / (n - i as f64))
            .collect();
    }
    if c.len() <= 1 || roots.len() >= most {
        return true;
    }
    if c.iter().all(|x| *x > zero) || c.iter().all(|x| *x < -zero) {
        return true;
    }
    let middle = cell[0] + (cell[1] - cell[0]) / 2.0;
    if depth == ROOT_DEPTH {
        roots.push(middle);
        return true;
    }
    // Paid for before it is done.
    if *allowance == 0 {
        return false;
    }
    *allowance -= 1;
    let n = c.len();
    let (mut left, mut right, mut row) = (c.clone(), c.clone(), c);
    for level in 1..n {
        for i in 0..n - level {
            row[i] = (row[i] + row[i + 1]) / 2.0;
        }
        left[level] = row[0];
        right[n - 1 - level] = row[n - 1 - level];
    }
    // The value at the middle is the shared end of the two halves.
    if right[0].abs() <= zero {
        roots.push(middle);
    }
    halve(
        &left,
        zero,
        [cell[0], middle],
        depth + 1,
        most,
        roots,
        allowance,
    ) && halve(
        &right,
        zero,
        [middle, cell[1]],
        depth + 1,
        most,
        roots,
        allowance,
    )
}

/// The `count` Chebyshev points of the first kind in `(0, 1)`, ascending:
/// `(1 - cos((2 j + 1) pi / (2 count))) / 2`.
/// Source: Mason and Handscomb, "Chebyshev Polynomials", 2003, section 6.3.
/// With an even count no point is a rational share of the piece: the cosine
/// of a rational multiple of pi is rational only at multiples of pi / 3 and
/// pi / 2 (Niven's theorem), and `(2 j + 1) / (2 count)` is neither when the
/// count is even. So no point can coincide with a builder that checks at
/// rational shares such as eighths.
pub(super) fn chebyshev_shares(count: usize) -> Vec<f64> {
    (0..count)
        .map(|j| {
            let angle = std::f64::consts::PI * (2 * j + 1) as f64 / (2 * count) as f64;
            (1.0 - sin_cos(angle).1) / 2.0
        })
        .collect()
}

/// Source: Orszag, "On the elimination of aliasing in finite-difference
/// schemes by filtering high-wavenumber components", J. Atmos. Sci. 28,
/// 1971: aliased content lands in the top third of the retained modes, so
/// that third is the part of the spectrum read as the tail.
const TAIL_SHARE_DENOMINATOR: usize = 3;

/// The size of the high-order tail of the Chebyshev series through the
/// samples of a vector function at first-kind points: the sum over the top
/// third of the coefficients of their Euclidean length.
///
/// Coefficients come from the discrete orthogonality of the Chebyshev
/// polynomials on these points, `c_k = (2 / n) sum_j f(x_j) T_k(x_j)`
/// (Mason and Handscomb, section 6.3). The interpolant's error is at most
/// twice the sum of the magnitudes of the neglected coefficients (Trefethen,
/// "Approximation Theory and Approximation Practice", 2013, chapter 4,
/// the aliasing formula and the bound that follows from it), and the last
/// retained coefficients estimate the neglected ones (Boyd, "Chebyshev and
/// Fourier Spectral Methods", 2nd ed., section 2.12).
/// The sample and mode pairs one [`chebyshev_tail`] of `count` samples
/// visits: every sample for each coefficient of the top third. It is the
/// work of the transform, counted before the transform runs.
pub(super) fn chebyshev_tail_pairs(count: usize) -> u64 {
    (count * (count / TAIL_SHARE_DENOMINATOR)) as u64
}

pub(super) fn chebyshev_tail(samples: &[[f64; 3]]) -> f64 {
    let n = samples.len();
    // cos(pi m / (2 n)) for every m that k (2 j + 1) can reach modulo 4 n.
    let table: Vec<f64> = (0..4 * n)
        .map(|m| sin_cos(std::f64::consts::PI * m as f64 / (2 * n) as f64).1)
        .collect();
    let first = n - n / TAIL_SHARE_DENOMINATOR;
    (first..n)
        .map(|k| {
            let mut coefficient = [0.0; 3];
            // Samples are ascending in the share, which is descending in the
            // Chebyshev variable; the order only flips the sign of odd terms.
            for (j, sample) in samples.iter().enumerate() {
                let weight = table[(k * (2 * j + 1)) % (4 * n)];
                for axis in 0..3 {
                    coefficient[axis] += sample[axis] * weight;
                }
            }
            let scale = 2.0 / n as f64;
            (coefficient[0] * coefficient[0]
                + coefficient[1] * coefficient[1]
                + coefficient[2] * coefficient[2])
                .sqrt()
                * scale
        })
        .sum()
}

/// Source: Brent, "Algorithms for Minimization without Derivatives", 1973,
/// chapter 5: the golden section ratio `(3 - sqrt(5)) / 2`.
const GOLDEN: f64 = 0.381_966_011_250_105_1;

/// Source: judgement. Brent's method needs about ten to forty evaluations on
/// a smooth bracket; one hundred steps only ends a search that stalls.
const BRENT_STEPS: u32 = 100;

/// The maximum of `f` on `[0, 1]` by Brent's `localmin` applied to `-f`:
/// golden section steps with parabolic interpolation where it is safe.
/// Source: Brent 1973, chapter 5, procedure `localmin`. The relative
/// tolerance is the square root of the machine epsilon, the smallest Brent
/// recommends, because a smooth function is flat to rounding within that
/// distance of its extremum. The caller maps its bracket onto `[0, 1]`, so
/// the tolerance is a share of the bracket and not of the parameter.
/// Returns the best share, its value and the number of evaluations.
pub(super) fn brent_maximum(f: &mut dyn FnMut(f64) -> f64) -> (f64, f64, u64) {
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
    (x, -fx, evaluations)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sine_and_cosine_match_the_library_to_rounding() {
        for i in -4000..=4000 {
            let x = f64::from(i) * 0.013_7;
            let (sine, cosine) = sin_cos(x);
            assert!(
                (sine - x.sin()).abs() <= TRIG_ROUNDINGS * UNIT_ROUNDOFF,
                "sin {x}"
            );
            assert!(
                (cosine - x.cos()).abs() <= TRIG_ROUNDINGS * UNIT_ROUNDOFF,
                "cos {x}"
            );
        }
    }

    #[test]
    fn basis_functions_sum_to_one_and_continue_past_the_border() {
        let knots = [0.0, 0.0, 0.0, 0.0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0];
        let (mut values, mut scratch) = (Vec::new(), Vec::new());
        for t in [0.0, 0.1, 0.3, 0.5, 0.7, 0.99, 1.0, 1.2, -0.2] {
            let span = find_span(&knots, 3, 6, t);
            basis_functions(&knots, &[0.0; 10], 3, span, t, &mut values, &mut scratch);
            assert!((values.iter().sum::<f64>() - 1.0).abs() < 1e-14, "{t}");
        }
        assert_eq!(find_span(&knots, 3, 6, 1.0), 5);
        assert_eq!(find_span(&knots, 3, 6, 1.2), 5);
        assert_eq!(find_span(&knots, 3, 6, -0.2), 3);
    }

    #[test]
    fn knot_remainders_keep_narrow_spans_accurate() {
        // A hat 2^-40 wide whose three knots are not binary64 numbers. With
        // the remainders the hat's value at a parameter is right to rounding;
        // without them the knots are known to a few percent of the span.
        let centre = 2.2;
        let step = 2.0_f64.powi(-41);
        let third = 2.0_f64.powi(-43) / 3.0;
        let knots = [0.0, 0.0, centre - step, centre, centre + step, 4.0, 4.0];
        let low = [0.0, 0.0, third, -third, third, 0.0, 0.0];
        let (mut values, mut scratch) = (Vec::new(), Vec::new());
        let t = centre - step / 4.0;
        // Span 2 is [centre - step, centre]; the hat is its second function.
        basis_functions(&knots, &low, 1, 2, t, &mut values, &mut scratch);
        let expected = (step * 0.75 - third) / (step - 2.0 * third);
        assert!((values[1] - expected).abs() <= 4.0 * f64::EPSILON);
        basis_functions(&knots, &[0.0; 7], 1, 2, t, &mut values, &mut scratch);
        assert!((values[1] - expected).abs() > 1e-3);
    }

    #[test]
    fn a_bezier_piece_reproduces_its_span() {
        let knots = [0.0, 0.0, 0.0, 0.0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0];
        let controls = [[1.0], [-2.0], [0.5], [4.0], [-1.0], [3.0]];
        let piece = bezier_piece(&knots, 3, &controls, 4, 0.4, 0.6);
        let (mut values, mut scratch) = (Vec::new(), Vec::new());
        for (t, expected) in [(0.4, piece[0][0]), (0.6, piece[3][0])] {
            basis_functions(&knots, &[0.0; 10], 3, 4, t, &mut values, &mut scratch);
            let value: f64 = (0..4).map(|i| values[i] * controls[1 + i][0]).sum();
            assert!((value - expected).abs() < 1e-14);
        }
    }

    #[test]
    fn sign_changes_finds_every_crossing_and_none_on_a_zero_polynomial() {
        // (s - 1/4)(s - 3/5) in Bernstein form on [0, 1].
        let power = [0.15, -0.85, 1.0];
        let bernstein = [
            power[0],
            power[0] + power[1] / 2.0,
            power[0] + power[1] + power[2],
        ];
        let mut roots = Vec::new();
        let mut allowance = u64::MAX;
        let halvings = sign_changes(&bernstein, 1e-15, &mut roots, &mut allowance).expect("paid");
        assert_eq!(u64::MAX - allowance, halvings);
        roots.sort_by(f64::total_cmp);
        assert_eq!(roots.len(), 2);
        assert!((roots[0] - 0.25).abs() < 1e-10 && (roots[1] - 0.6).abs() < 1e-10);
        let mut none = Vec::new();
        let free = |coefficients: &[f64], roots: &mut Vec<f64>| {
            sign_changes(coefficients, 1e-15, roots, &mut u64::MAX.clone())
        };
        assert_eq!(free(&[0.0, 0.0, 0.0], &mut none), Some(0));
        // A root at the piece's own end is not a crossing inside it.
        assert_eq!(free(&[0.0, 1.05], &mut none), Some(0));
        assert_eq!(free(&[-0.3, -0.1, 0.0], &mut none), Some(0));
        free(&[1.0, 2.0, 0.5], &mut none);
        assert!(none.is_empty());
        // The allowance is taken before each halving. With exactly enough
        // the search finishes on its last unit; with one less it stops, and
        // with none it does not halve at all.
        let mut exact = halvings;
        assert_eq!(
            sign_changes(&bernstein, 1e-15, &mut Vec::new(), &mut exact),
            Some(halvings)
        );
        assert_eq!(exact, 0);
        let mut short = halvings - 1;
        assert_eq!(
            sign_changes(&bernstein, 1e-15, &mut Vec::new(), &mut short),
            None
        );
        assert_eq!(short, 0);
        let (mut empty, mut found) = (0, Vec::new());
        assert_eq!(
            sign_changes(&bernstein, 1e-15, &mut found, &mut empty),
            None
        );
        assert!(found.is_empty());
    }

    #[test]
    fn chebyshev_tail_is_quiet_on_a_resolved_polynomial_and_loud_on_a_spike() {
        let shares = chebyshev_shares(18);
        assert!(shares.windows(2).all(|pair| pair[0] < pair[1]));
        let smooth: Vec<[f64; 3]> = shares
            .iter()
            .map(|s| [s * s * s - s, 2.0 * s, 1.0])
            .collect();
        assert!(chebyshev_tail(&smooth) < 1e-14);
        let spike: Vec<[f64; 3]> = shares
            .iter()
            .map(|s| [(-((s - 0.4137) / 0.05).powi(2)).exp(), 0.0, 0.0])
            .collect();
        assert!(chebyshev_tail(&spike) > 1e-3);
    }

    #[test]
    fn no_sample_share_is_a_builder_fraction() {
        // The fitter checks at eighths today. Thirds, sixteenths and
        // sixty-fourths are covered as well.
        for count in (18..=162).step_by(2) {
            for share in chebyshev_shares(count) {
                for denominator in [3.0, 64.0] {
                    let scaled = share * denominator;
                    assert!(
                        (scaled - scaled.round()).abs() > 1e-9 * denominator,
                        "{count}: {share}"
                    );
                }
            }
        }
    }

    #[test]
    fn brent_finds_a_smooth_and_a_kinked_maximum() {
        let (x, value, evaluations) = brent_maximum(&mut |s| 1.0 - (s - 0.3217) * (s - 0.3217));
        assert!((x - 0.3217).abs() < 1e-7 && (value - 1.0).abs() < 1e-14);
        assert!(evaluations < 40);
        let (x, _, _) = brent_maximum(&mut |s| -(s - 0.7).abs());
        assert!((x - 0.7).abs() < 1e-6);
    }
}
