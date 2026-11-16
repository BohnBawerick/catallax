//! Binary64 evaluation and inversion that only propose sample locations.
//! Every proposed value is rechecked by numerical admission at 128 bits.

use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::{BSplineCurve, BSplineSurface};

fn real(value: &Real) -> f64 {
    value.to_f64()
}

fn span(knots: &[f64], degree: usize, count: usize, t: f64) -> usize {
    if t >= knots[count] {
        let mut index = count - 1;
        while index > degree && knots[index] >= knots[count] {
            index -= 1;
        }
        return index;
    }
    let mut index = degree;
    while index + 1 < count && knots[index + 1] <= t {
        index += 1;
    }
    index
}

// Values and first derivatives of the nonzero basis functions at `t`
// (The NURBS Book, A2.3, first derivatives only).
fn basis(knots: &[f64], degree: usize, index: usize, t: f64) -> (Vec<f64>, Vec<f64>) {
    let p = degree;
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots[index + 1 - j];
        right[j] = knots[index + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let values: Vec<f64> = (0..=p).map(|j| ndu[j][p]).collect();
    let mut first = vec![0.0; p + 1];
    if p > 0 {
        for (r, slot) in first.iter_mut().enumerate() {
            let mut d = 0.0;
            if r >= 1 {
                d += ndu[r - 1][p - 1] / ndu[p][r - 1];
            }
            if r < p {
                d -= ndu[r][p - 1] / ndu[p][r];
            }
            *slot = d * p as f64;
        }
    }
    (values, first)
}

/// A binary64 copy of a finite rational B-spline surface.
pub(super) struct Surface64 {
    degrees: [usize; 2],
    knots: [Vec<f64>; 2],
    count: [usize; 2],
    controls: Vec<[f64; 4]>,
}

impl Surface64 {
    pub fn new(source: &BSplineSurface) -> Self {
        Self {
            degrees: source.degrees(),
            knots: source
                .knots()
                .map(|row| row.iter().map(|q| q.to_f64_nearest()).collect()),
            count: source.control_count(),
            controls: source
                .controls()
                .iter()
                .zip(source.weights())
                .map(|(p, w)| {
                    let w = real(w);
                    [real(&p.x) * w, real(&p.y) * w, real(&p.z) * w, w]
                })
                .collect(),
        }
    }

    /// Position and first partial derivatives at a clamped location.
    pub fn jet(&self, uv: [f64; 2]) -> [[f64; 3]; 3] {
        let mut basis_rows = Vec::with_capacity(2);
        let mut spans = [0; 2];
        for axis in 0..2 {
            let knots = &self.knots[axis];
            let degree = self.degrees[axis];
            let count = self.count[axis];
            let t = uv[axis].clamp(knots[degree], knots[count]);
            spans[axis] = span(knots, degree, count, t);
            basis_rows.push(basis(knots, degree, spans[axis], t));
        }
        let mut a = [[0.0; 4]; 3];
        let [pu, pv] = self.degrees;
        for i in 0..=pu {
            let row = spans[0] - pu + i;
            for j in 0..=pv {
                let column = spans[1] - pv + j;
                let control = &self.controls[row * self.count[1] + column];
                let weights = [
                    basis_rows[0].0[i] * basis_rows[1].0[j],
                    basis_rows[0].1[i] * basis_rows[1].0[j],
                    basis_rows[0].0[i] * basis_rows[1].1[j],
                ];
                for (k, weight) in weights.into_iter().enumerate() {
                    for c in 0..4 {
                        a[k][c] += weight * control[c];
                    }
                }
            }
        }
        let w = a[0][3];
        let point: [f64; 3] = std::array::from_fn(|c| a[0][c] / w);
        let du: [f64; 3] = std::array::from_fn(|c| (a[1][c] - a[1][3] * point[c]) / w);
        let dv: [f64; 3] = std::array::from_fn(|c| (a[2][c] - a[2][3] * point[c]) / w);
        [point, du, dv]
    }

    pub fn point(&self, uv: [f64; 2]) -> [f64; 3] {
        self.jet(uv)[0]
    }

    /// Closest location in `region`, refined by Gauss-Newton from `seed`, or
    /// from the best location of a coarse grid when no seed is given.
    pub fn invert(
        &self,
        target: [f64; 3],
        region: [[f64; 2]; 2],
        seed: Option<[f64; 2]>,
        dec: &mut Decider,
    ) -> Result<[f64; 2], Refusal> {
        let mut uv = match seed {
            Some(seed) => seed,
            None => {
                const GRID: usize = 24;
                let mut best = (f64::INFINITY, [region[0][0], region[1][0]]);
                dec.charge_work(((GRID + 1) * (GRID + 1)) as u64)?;
                for i in 0..=GRID {
                    for j in 0..=GRID {
                        let uv = [
                            lerp(region[0], i as f64 / GRID as f64),
                            lerp(region[1], j as f64 / GRID as f64),
                        ];
                        let d = distance(self.point(uv), target);
                        if d < best.0 {
                            best = (d, uv);
                        }
                    }
                }
                best.1
            }
        };
        for _ in 0..64 {
            dec.charge_work(1)?;
            let [p, du, dv] = self.jet(uv);
            let r = sub(target, p);
            let a = [[dot(du, du), dot(du, dv)], [dot(du, dv), dot(dv, dv)]];
            let b = [dot(du, r), dot(dv, r)];
            let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
            if det.is_nan() || det.abs() <= 0.0 {
                break;
            }
            let step = [
                (b[0] * a[1][1] - b[1] * a[0][1]) / det,
                (a[0][0] * b[1] - a[1][0] * b[0]) / det,
            ];
            let next = [
                (uv[0] + step[0]).clamp(region[0][0], region[0][1]),
                (uv[1] + step[1]).clamp(region[1][0], region[1][1]),
            ];
            let moved = (next[0] - uv[0]).abs() + (next[1] - uv[1]).abs();
            uv = next;
            if moved <= 1e-15 * (1.0 + uv[0].abs() + uv[1].abs()) {
                break;
            }
        }
        // A location on a region side within binary64 noise is that side.
        for axis in 0..2 {
            let width = region[axis][1] - region[axis][0];
            for side in region[axis] {
                if (uv[axis] - side).abs() <= 1e-12 * width {
                    uv[axis] = side;
                }
            }
        }
        Ok(uv)
    }
}

/// A binary64 copy of a finite rational B-spline curve.
pub(super) struct Curve64 {
    degree: usize,
    knots: Vec<f64>,
    controls: Vec<[f64; 4]>,
}

impl Curve64 {
    pub fn new(source: &BSplineCurve) -> Self {
        Self {
            degree: source.degree(),
            knots: source.knots().iter().map(|q| q.to_f64_nearest()).collect(),
            controls: source
                .controls()
                .zip(source.weights())
                .map(|(p, w)| {
                    let w = real(w);
                    [real(&p.x) * w, real(&p.y) * w, real(&p.z) * w, w]
                })
                .collect(),
        }
    }

    pub fn point(&self, t: f64) -> [f64; 3] {
        let count = self.controls.len();
        let t = t.clamp(self.knots[self.degree], self.knots[count]);
        let index = span(&self.knots, self.degree, count, t);
        let (values, _) = basis(&self.knots, self.degree, index, t);
        let mut a = [0.0; 4];
        for (i, value) in values.iter().enumerate() {
            let control = &self.controls[index - self.degree + i];
            for c in 0..4 {
                a[c] += value * control[c];
            }
        }
        std::array::from_fn(|c| a[c] / a[3])
    }
}

pub(super) fn lerp(range: [f64; 2], share: f64) -> f64 {
    range[0] + share * (range[1] - range[0])
}

pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    dot(d, d).sqrt()
}

/// Cubic not-a-knot interpolation at strictly increasing parameters
/// (de Boor, A Practical Guide to Splines, ch. XIII). The collocation
/// matrix is totally positive, so banded elimination needs no pivoting.
/// The chord through the end values is removed before the solve and its
/// Greville controls restored after it, so a constant coordinate, such as an
/// iso-parameter trim, keeps exactly its value.
/// Returns the expanded knots and one control row per data row.
pub(super) fn interpolate<const D: usize>(
    parameters: &[f64],
    values: &[[f64; D]],
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<[f64; D]>), Refusal> {
    let n = parameters.len();
    assert!(n >= 4 && values.len() == n);
    let mut knots = vec![parameters[0]; 4];
    knots.extend_from_slice(&parameters[2..n - 2]);
    knots.extend(std::iter::repeat_n(parameters[n - 1], 4));
    let width = parameters[n - 1] - parameters[0];
    let chord = |t: f64| -> [f64; D] {
        let share = (t - parameters[0]) / width;
        std::array::from_fn(|c| values[0][c] + share * (values[n - 1][c] - values[0][c]))
    };
    // Band storage: column j of row i lives at slot j + 3 - i.
    let mut band = vec![[0.0f64; 7]; n];
    let mut rhs: Vec<[f64; D]> = parameters
        .iter()
        .zip(values)
        .map(|(t, value)| {
            let line = chord(*t);
            std::array::from_fn(|c| value[c] - line[c])
        })
        .collect();
    for (i, t) in parameters.iter().enumerate() {
        dec.charge_work(1)?;
        let index = span(&knots, 3, n, *t);
        let (basis, _) = basis(&knots, 3, index, *t);
        for (k, value) in basis.into_iter().enumerate() {
            let column = index - 3 + k;
            let slot = column + 3 - i;
            band[i][slot] = value;
        }
    }
    for k in 0..n {
        dec.charge_work(1)?;
        let pivot = band[k][3];
        if pivot.is_nan() || pivot.abs() <= 1e-300 {
            return Err(Refusal::new(
                "step.numerical-fit",
                "interpolation matrix lost a pivot",
            ));
        }
        for i in k + 1..(k + 4).min(n) {
            let factor = band[i][k + 3 - i] / pivot;
            for j in k..(k + 4).min(n) {
                band[i][j + 3 - i] -= factor * band[k][j + 3 - k];
            }
            for c in 0..D {
                rhs[i][c] -= factor * rhs[k][c];
            }
        }
    }
    let mut controls = vec![[0.0; D]; n];
    for k in (0..n).rev() {
        let mut value = rhs[k];
        for j in k + 1..(k + 4).min(n) {
            for c in 0..D {
                value[c] -= band[k][j + 3 - k] * controls[j][c];
            }
        }
        for c in 0..D {
            controls[k][c] = value[c] / band[k][3];
        }
    }
    for (i, control) in controls.iter_mut().enumerate() {
        let line = chord((knots[i + 1] + knots[i + 2] + knots[i + 3]) / 3.0);
        for c in 0..D {
            control[c] += line[c];
        }
    }
    controls[0] = values[0];
    controls[n - 1] = values[n - 1];
    Ok((knots, controls))
}

/// Joins clamped cubic pieces, each starting where the previous one ends, into
/// one curve that is C0 at the joins (knot multiplicity three).
pub(super) fn concatenate<const D: usize>(
    pieces: Vec<(Vec<f64>, Vec<[f64; D]>)>,
) -> (Vec<f64>, Vec<[f64; D]>) {
    let count = pieces.len();
    let mut knots = Vec::new();
    let mut controls = Vec::new();
    for (index, (piece_knots, piece_controls)) in pieces.into_iter().enumerate() {
        let first = if index == 0 { 0 } else { 4 };
        let last = if index + 1 == count {
            piece_knots.len()
        } else {
            piece_knots.len() - 1
        };
        knots.extend_from_slice(&piece_knots[first..last]);
        controls.extend_from_slice(&piece_controls[usize::from(index != 0)..]);
    }
    (knots, controls)
}

/// Binary64 evaluation of an [`interpolate`] result.
pub(super) fn evaluate<const D: usize>(knots: &[f64], controls: &[[f64; D]], t: f64) -> [f64; D] {
    let n = controls.len();
    let t = t.clamp(knots[3], knots[n]);
    let index = span(knots, 3, n, t);
    let (basis, _) = basis(knots, 3, index, t);
    let mut result = [0.0; D];
    for (k, value) in basis.into_iter().enumerate() {
        for c in 0..D {
            result[c] += value * controls[index - 3 + k][c];
        }
    }
    result
}
