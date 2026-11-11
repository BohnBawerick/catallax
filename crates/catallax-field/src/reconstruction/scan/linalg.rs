//! Small binary64 vector and matrix helpers for candidate fitting. Every loop
//! runs in a fixed order, so equal input gives equal bits.

pub(super) type V3 = [f64; 3];

pub(super) fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(super) fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub(super) fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub(super) fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(super) fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// `a` scaled to unit length; `None` for a zero or nonfinite vector.
pub(super) fn unit(a: V3) -> Option<V3> {
    let n = norm(a);
    (n > 0.0 && n.is_finite()).then(|| scale(a, 1.0 / n))
}

/// The sign that makes the largest component of a direction positive, so a
/// line's direction has one spelling. Ties go to the lowest axis.
pub(super) fn canonical(a: V3) -> V3 {
    let mut k = 0;
    for i in 1..3 {
        if a[i].abs() > a[k].abs() {
            k = i;
        }
    }
    if a[k] < 0.0 { scale(a, -1.0) } else { a }
}

/// Two unit vectors completing `a` to a right-handed orthonormal frame.
pub(super) fn frame(a: V3) -> (V3, V3) {
    let mut k = 0;
    for i in 1..3 {
        if a[i].abs() < a[k].abs() {
            k = i;
        }
    }
    let mut e = [0.0; 3];
    e[k] = 1.0;
    let u = unit(cross(a, e)).expect("a unit axis is not parallel to its smallest axis");
    (u, cross(a, u))
}

/// Mean and covariance of points, summed in index order.
pub(super) fn moments(points: impl Iterator<Item = V3> + Clone) -> (V3, [[f64; 3]; 3], usize) {
    let mut n = 0usize;
    let mut sum = [0.0; 3];
    for p in points.clone() {
        sum = add(sum, p);
        n += 1;
    }
    let mean = scale(sum, 1.0 / n.max(1) as f64);
    let mut c = [[0.0; 3]; 3];
    for p in points {
        let d = sub(p, mean);
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] += d[i] * d[j];
            }
        }
    }
    let inv = 1.0 / n.max(1) as f64;
    (mean, c.map(|row| row.map(|x| x * inv)), n)
}

/// Eigenvalues in ascending order and their unit eigenvectors of a symmetric
/// 3x3 matrix, by cyclic Jacobi rotations.
pub(super) fn eigen(mut a: [[f64; 3]; 3]) -> ([f64; 3], [V3; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..64 {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        let diag = a[0][0].abs() + a[1][1].abs() + a[2][2].abs();
        if off <= f64::EPSILON * diag * 1e-3 || off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for row in &mut a {
                let (akp, akq) = (row[p], row[q]);
                row[p] = c * akp - s * akq;
                row[q] = s * akp + c * akq;
            }
            let (rp, rq) = (a[p], a[q]);
            for k in 0..3 {
                a[p][k] = c * rp[k] - s * rq[k];
                a[q][k] = s * rp[k] + c * rq[k];
            }
            for row in &mut v {
                let (vp, vq) = (row[p], row[q]);
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
    }
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]).then(i.cmp(&j)));
    let values = order.map(|i| a[i][i]);
    let vectors = order.map(|i| [v[0][i], v[1][i], v[2][i]]);
    (values, vectors)
}

/// Solves the `n`-by-`n` system `m x = b` by Gaussian elimination with
/// partial pivoting; `None` when it is singular to working precision.
pub(super) fn solve<const N: usize>(mut m: [[f64; N]; N], mut b: [f64; N]) -> Option<[f64; N]> {
    for col in 0..N {
        let mut pivot = col;
        for row in col + 1..N {
            if m[row][col].abs() > m[pivot][col].abs() {
                pivot = row;
            }
        }
        if m[pivot][col].abs() < 1e-300 || !m[pivot][col].is_finite() {
            return None;
        }
        m.swap(col, pivot);
        b.swap(col, pivot);
        let pivot_row = m[col];
        for row in col + 1..N {
            let f = m[row][col] / pivot_row[col];
            for (target, source) in m[row][col..].iter_mut().zip(&pivot_row[col..]) {
                *target -= f * source;
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = [0.0; N];
    for row in (0..N).rev() {
        let mut s = b[row];
        for k in row + 1..N {
            s -= m[row][k] * x[k];
        }
        x[row] = s / m[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}
