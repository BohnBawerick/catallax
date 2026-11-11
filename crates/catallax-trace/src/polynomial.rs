// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: results are returned values instead of output
// references; the solvers always succeed, as in SMLib.
// Translated from SMLib SmPolynomial.cpp (SmPolynomial::SolveCubicEqn,
// sm_AddSolution, SmPolynomial::SolveQuarticEqn); rust-brep has no port.
//! The cubic and quartic root finders of SMLib (Graphics Gems I), which the
//! conic intersections of the Dupin indicatrices use. The quadratic one is
//! [`crate::geometry::solve_quadratic`].

use crate::geometry::solve_quadratic;
use crate::scalar::{PI, acos, cbrt, cos, sqrt};

/// `SmPolynomial::SolveCubicEqn`: the real roots of
/// `c[0] + c[1] x + c[2] x^2 + c[3] x^3`, a double root once. A negligible
/// leading coefficient falls back to the quadratic.
pub fn solve_cubic(c: [f64; 4], zero: f64) -> (usize, [f64; 3]) {
    let mut roots = [0.0; 3];
    if c[3].abs() < zero {
        let (count, quadratic) = solve_quadratic([c[0], c[1], c[2]], zero);
        roots[..count].copy_from_slice(&quadratic[..count]);
        return (count, roots);
    }
    let a = c[2] / c[3];
    let b = c[1] / c[3];
    let c0 = c[0] / c[3];
    let sq_a = a * a;
    let p = 1.0 / 3.0 * (-1.0 / 3.0 * sq_a + b);
    let q = 1.0 / 2.0 * (2.0 / 27.0 * a * sq_a - 1.0 / 3.0 * a * b + c0);
    let cb_p = p * p * p;
    let d = q * q + cb_p;
    let count;
    if d.abs() < zero {
        if q.abs() < zero {
            roots[0] = 0.0;
            count = 1;
        } else {
            let u = cbrt(-q);
            roots[0] = 2.0 * u;
            roots[1] = -u;
            count = 2;
        }
    } else if d < 0.0 {
        // SMLib takes the library square root here, not `sm_sqrt`.
        let phi = 1.0 / 3.0 * acos(-q / (-cb_p).sqrt());
        let t = 2.0 * sqrt(-p);
        roots[0] = t * cos(phi);
        roots[1] = -t * cos(phi + PI / 3.0);
        roots[2] = -t * cos(phi - PI / 3.0);
        count = 3;
    } else {
        let sqrt_d = sqrt(d);
        let u = cbrt(sqrt_d - q);
        let v = -cbrt(sqrt_d + q);
        roots[0] = u + v;
        count = 1;
    }
    let sub = 1.0 / 3.0 * a;
    for root in &mut roots[..count] {
        *root -= sub;
    }
    (count, roots)
}

/// `sm_AddSolution`: appends `value` unless an earlier root lies within
/// `zero` of it.
fn add_solution(value: f64, zero: f64, count: &mut usize, roots: &mut [f64; 4]) {
    if roots[..*count]
        .iter()
        .any(|root| (value - root).abs() < zero)
    {
        return;
    }
    roots[*count] = value;
    *count += 1;
}

/// `SmPolynomial::SolveQuarticEqn`: the distinct real roots of
/// `c[0] + c[1] x + ... + c[4] x^4` by Ferrari's resolvent cubic. A
/// negligible leading coefficient falls back to the cubic.
pub fn solve_quartic(c: [f64; 5], zero: f64) -> (usize, [f64; 4]) {
    let mut roots = [0.0; 4];
    let mut count = 0;
    if c[4].abs() < zero {
        let (cubic_count, cubic) = solve_cubic([c[0], c[1], c[2], c[3]], zero);
        for &root in &cubic[..cubic_count] {
            add_solution(root, zero, &mut count, &mut roots);
        }
        return (count, roots);
    }
    let a = c[3] / c[4];
    let b = c[2] / c[4];
    let c1 = c[1] / c[4];
    let d = c[0] / c[4];
    let sq_a = a * a;
    let p = -3.0 / 8.0 * sq_a + b;
    let q = 1.0 / 8.0 * sq_a * a - 1.0 / 2.0 * a * b + c1;
    let r = (-3.0 / 256.0 * sq_a * sq_a) + (1.0 / 16.0 * sq_a * b) - (1.0 / 4.0 * a * c1) + d;
    if r.abs() < zero {
        let (cubic_count, cubic) = solve_cubic([q, p, 0.0, 1.0], zero);
        for &root in &cubic[..cubic_count] {
            add_solution(root, zero, &mut count, &mut roots);
        }
        add_solution(0.0, zero, &mut count, &mut roots);
    } else {
        let (_, cubic) = solve_cubic(
            [
                1.0 / 2.0 * r * p - 1.0 / 8.0 * q * q,
                -r,
                -1.0 / 2.0 * p,
                1.0,
            ],
            zero,
        );
        let z = cubic[0];
        let mut u = z * z - r;
        let mut v = 2.0 * z - p;
        if u.abs() < zero {
            u = 0.0;
        } else if u > 0.0 {
            u = sqrt(u);
        } else {
            return (0, roots);
        }
        if v.abs() < zero {
            v = 0.0;
        } else if v > 0.0 {
            v = sqrt(v);
        } else {
            return (0, roots);
        }
        let (n, first) = solve_quadratic([z - u, if q < 0.0 { -v } else { v }, 1.0], zero);
        for &root in &first[..n] {
            add_solution(root, zero, &mut count, &mut roots);
        }
        let (n, second) = solve_quadratic([z + u, if q < 0.0 { v } else { -v }, 1.0], zero);
        for &root in &second[..n] {
            add_solution(root, zero, &mut count, &mut roots);
        }
    }
    let sub = 1.0 / 4.0 * a;
    for root in &mut roots[..count] {
        *root -= sub;
    }
    (count, roots)
}
