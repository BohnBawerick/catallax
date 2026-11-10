//! Exact orientation predicates on file coordinates, with a binary64 filter.
//!
//! Every coordinate is an exact rational. Its binary64 copy is the nearest
//! double, so each copy has relative error at most 2^-53. The filter accepts a
//! binary64 determinant only when it exceeds a bound that covers that input
//! error and the rounding of the evaluation; otherwise the determinant is
//! recomputed in exact rational arithmetic. No answer depends on the filter.

use std::cell::Cell;

use catallax_numbers::{Rational, Sign};

/// Relative bound on the filtered determinant error, in units of the
/// permanent. The analysis gives about 1.2e-15 for orient3d and less for
/// orient2d; the margin is deliberate.
const BOUND: f64 = 1e-14;
/// Every nonzero difference and magnitude the filter multiplies must lie in
/// `[2^-300, 2^300]`. A product of three such numbers is then a normal
/// binary64, so no product underflows or overflows and the relative error
/// analysis behind [`BOUND`] holds. Anything outside goes to exact arithmetic.
const LOW: f64 = 4.909093465297727e-91; // 2^-300
const HIGH: f64 = 2.037035976334486e90; // 2^300

thread_local! {
    /// Exact rational work done since the last [`take_exact_work`], in the
    /// units of the numbers crate: arithmetic operations times limb products.
    static EXACT_WORK: Cell<u64> = const { Cell::new(0) };
}

/// Takes the exact work recorded on this thread since the last call.
pub(crate) fn take_exact_work() -> u64 {
    EXACT_WORK.with(|w| w.replace(0))
}

/// Exact work recorded on this thread and not yet taken.
pub(crate) fn pending_exact_work() -> u64 {
    EXACT_WORK.with(Cell::get)
}

fn limbs(r: &Rational) -> u64 {
    (r.numerator().bit_length() + r.denominator().bit_length()) / 64 + 1
}

/// Records `operations` exact operations on numbers as large as the largest of `inputs`.
pub(crate) fn record_exact<'a>(operations: u64, inputs: impl IntoIterator<Item = &'a Rational>) {
    let size = inputs.into_iter().map(limbs).max().unwrap_or(1);
    let cost = operations.saturating_mul(size.saturating_mul(size));
    EXACT_WORK.with(|w| w.set(w.get().saturating_add(cost)));
}

pub(crate) fn in_range<'a>(values: impl IntoIterator<Item = &'a f64>) -> bool {
    values.into_iter().all(|&x| {
        let x = x.abs();
        x == 0.0 || (LOW..=HIGH).contains(&x)
    })
}

/// One point: its exact coordinates and their nearest binary64 copies.
#[derive(Clone, Copy)]
pub(crate) struct P<'a> {
    pub e: &'a [Rational; 3],
    pub f: [f64; 3],
}

impl<'a> P<'a> {
    pub(crate) fn new(e: &'a [Rational; 3]) -> Self {
        Self {
            e,
            f: e.each_ref().map(Rational::to_f64_nearest),
        }
    }
}

/// Sign of `det[b - a, c - a, d - a]`: positive when `d` lies on the side of
/// the plane `abc` that the normal `(b - a) x (c - a)` points to.
pub(crate) fn orient3d(a: P, b: P, c: P, d: P) -> Sign {
    let diff = |p: P| std::array::from_fn::<f64, 3, _>(|i| p.f[i] - a.f[i]);
    let size = |p: P| std::array::from_fn::<f64, 3, _>(|i| p.f[i].abs() + a.f[i].abs());
    let (u, v, w) = (diff(b), diff(c), diff(d));
    let (su, sv, sw) = (size(b), size(c), size(d));
    if in_range([u, v, w, su, sv, sw].iter().flatten()) {
        let det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0]);
        let permanent = su[0] * (sv[1] * sw[2] + sv[2] * sw[1])
            + su[1] * (sv[0] * sw[2] + sv[2] * sw[0])
            + su[2] * (sv[0] * sw[1] + sv[1] * sw[0]);
        let bound = BOUND * permanent;
        if det > bound {
            return Sign::Positive;
        }
        if det < -bound {
            return Sign::Negative;
        }
    }
    record_exact(23, a.e.iter().chain(b.e).chain(c.e).chain(d.e));
    let diff = |p: P| std::array::from_fn::<Rational, 3, _>(|i| &p.e[i] - &a.e[i]);
    let (u, v, w) = (diff(b), diff(c), diff(d));
    let det = &(&u[0] * &(&(&v[1] * &w[2]) - &(&v[2] * &w[1])))
        - &(&u[1] * &(&(&v[0] * &w[2]) - &(&v[2] * &w[0])));
    let det = &det + &(&u[2] * &(&(&v[0] * &w[1]) - &(&v[1] * &w[0])));
    det.sign()
}

/// The two kept coordinates when axis `k` is dropped, in cyclic order, so the
/// projected orientation of a triangle equals the sign of its normal's `k`
/// component.
pub(crate) fn kept(k: usize) -> (usize, usize) {
    ((k + 1) % 3, (k + 2) % 3)
}

/// Sign of the projected `det[b - a, c - a]` with axis `k` dropped.
pub(crate) fn orient2d(k: usize, a: P, b: P, c: P) -> Sign {
    let (i, j) = kept(k);
    let entries = [
        b.f[i] - a.f[i],
        c.f[j] - a.f[j],
        b.f[j] - a.f[j],
        c.f[i] - a.f[i],
    ];
    let sizes = [
        b.f[i].abs() + a.f[i].abs(),
        c.f[j].abs() + a.f[j].abs(),
        b.f[j].abs() + a.f[j].abs(),
        c.f[i].abs() + a.f[i].abs(),
    ];
    if in_range(&entries) && in_range(&sizes) {
        let det = entries[0] * entries[1] - entries[2] * entries[3];
        let bound = BOUND * (sizes[0] * sizes[1] + sizes[2] * sizes[3]);
        if det > bound {
            return Sign::Positive;
        }
        if det < -bound {
            return Sign::Negative;
        }
    }
    orient2d_exact(k, a, b, c).sign()
}

/// The exact projected determinant.
pub(crate) fn orient2d_exact(k: usize, a: P, b: P, c: P) -> Rational {
    let (i, j) = kept(k);
    record_exact(7, [&a.e[i], &a.e[j], &b.e[i], &b.e[j], &c.e[i], &c.e[j]]);
    &(&(&b.e[i] - &a.e[i]) * &(&c.e[j] - &a.e[j])) - &(&(&b.e[j] - &a.e[j]) * &(&c.e[i] - &a.e[i]))
}

/// A projection in which the triangle `abc` keeps a nonzero area, and the sign
/// of that area. `None` exactly when the triangle is degenerate.
pub(crate) fn projection(a: P, b: P, c: P) -> Option<(usize, Sign)> {
    // Try the axis with the largest approximate normal component first; the
    // exact sign decides, so the order only saves work.
    let u: [f64; 3] = std::array::from_fn(|i| b.f[i] - a.f[i]);
    let v: [f64; 3] = std::array::from_fn(|i| c.f[i] - a.f[i]);
    let n = [
        (u[1] * v[2] - u[2] * v[1]).abs(),
        (u[2] * v[0] - u[0] * v[2]).abs(),
        (u[0] * v[1] - u[1] * v[0]).abs(),
    ];
    let mut order = [0, 1, 2];
    order.sort_by(|&x, &y| n[y].total_cmp(&n[x]));
    order.into_iter().find_map(|k| match orient2d(k, a, b, c) {
        Sign::Zero => None,
        s => Some((k, s)),
    })
}

fn agrees(s: Sign, reference: Sign) -> bool {
    s == Sign::Zero || s == reference
}

/// Whether `p`, coplanar with the triangle, lies in the closed triangle.
/// `k` and `s` come from [`projection`] for that triangle.
pub(crate) fn in_triangle_2d(k: usize, s: Sign, t: [P; 3], p: P) -> bool {
    (0..3).all(|i| agrees(orient2d(k, t[i], t[(i + 1) % 3], p), s))
}

/// Whether coplanar point `p` lies on the closed segment `ab` under projection `k`.
fn on_segment_2d(k: usize, a: P, b: P, p: P) -> bool {
    if orient2d(k, a, b, p) != Sign::Zero {
        return false;
    }
    let (i, j) = kept(k);
    [i, j].into_iter().all(|axis| {
        let (lo, hi) = if a.e[axis] <= b.e[axis] {
            (&a.e[axis], &b.e[axis])
        } else {
            (&b.e[axis], &a.e[axis])
        };
        lo <= &p.e[axis] && &p.e[axis] <= hi
    })
}

/// Whether the closed coplanar segments `ab` and `cd` meet, under projection `k`.
fn segments_meet_2d(k: usize, a: P, b: P, c: P, d: P) -> bool {
    let o1 = orient2d(k, a, b, c);
    let o2 = orient2d(k, a, b, d);
    let o3 = orient2d(k, c, d, a);
    let o4 = orient2d(k, c, d, b);
    if o1 != Sign::Zero
        && o2 != Sign::Zero
        && o1 != o2
        && o3 != Sign::Zero
        && o4 != Sign::Zero
        && o3 != o4
    {
        return true;
    }
    on_segment_2d(k, a, b, c)
        || on_segment_2d(k, a, b, d)
        || on_segment_2d(k, c, d, a)
        || on_segment_2d(k, c, d, b)
}

/// Whether the closed segment `ab` meets the closed nondegenerate triangle `t`.
pub(crate) fn segment_meets_triangle(a: P, b: P, t: [P; 3]) -> bool {
    let sa = orient3d(t[0], t[1], t[2], a);
    let sb = orient3d(t[0], t[1], t[2], b);
    if sa == sb && sa != Sign::Zero {
        return false;
    }
    if sa == Sign::Zero && sb == Sign::Zero {
        let (k, s) = projection(t[0], t[1], t[2]).expect("triangles are nondegenerate");
        return in_triangle_2d(k, s, t, a)
            || in_triangle_2d(k, s, t, b)
            || (0..3).any(|i| segments_meet_2d(k, a, b, t[i], t[(i + 1) % 3]));
    }
    // The segment reaches the plane at one point. The line through it meets the
    // closed triangle exactly when no two edge orientations have strict
    // opposite signs.
    let o = [
        orient3d(a, b, t[0], t[1]),
        orient3d(a, b, t[1], t[2]),
        orient3d(a, b, t[2], t[0]),
    ];
    !(o.contains(&Sign::Positive) && o.contains(&Sign::Negative))
}

/// Whether coplanar point `a` lies in the closed angle of triangle `t` at its
/// vertex `t[0]`, so the segment from `t[0]` towards `a` enters the triangle.
pub(crate) fn in_corner_2d(k: usize, s: Sign, t: [P; 3], a: P) -> bool {
    agrees(orient2d(k, t[0], t[1], a), s) && agrees(orient2d(k, t[0], a, t[2]), s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exact(x: &str) -> Rational {
        Rational::parse(x).unwrap()
    }

    #[test]
    fn an_underflowing_product_falls_back_to_exact_arithmetic() {
        // Exact det = 1e-100 - 1e-200 > 0. In binary64, v[1] * w[2] = 1e-400
        // underflows to zero and the filtered determinant is -1e-200.
        let points = [
            ["0", "0", "0"],
            ["1e300", "1e-100", "0"],
            ["1e100", "1e-200", "0"],
            ["0", "0", "1e-200"],
        ]
        .map(|p| p.map(exact));
        let [a, b, c, d] = points.each_ref().map(P::new);
        assert_eq!(orient3d(a, b, c, d), Sign::Positive);
        assert_eq!(orient3d(a, c, b, d), Sign::Negative);
        assert!(
            take_exact_work() > 0,
            "the exact path was taken and metered"
        );
    }

    #[test]
    fn an_overflowing_product_falls_back_to_exact_arithmetic() {
        let points = [
            ["0", "0", "0"],
            ["1e200", "0", "0"],
            ["0", "1e200", "0"],
            ["0", "0", "1e-300"],
        ]
        .map(|p| p.map(exact));
        let [a, b, c, d] = points.each_ref().map(P::new);
        assert_eq!(orient3d(a, b, c, d), Sign::Positive);
        let [a, b, c] = [&points[0], &points[1], &points[2]].map(P::new);
        assert_eq!(orient2d(2, a, b, c), Sign::Positive);
    }
}
