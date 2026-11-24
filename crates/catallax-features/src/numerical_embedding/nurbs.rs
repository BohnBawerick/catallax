//! Binary64 evaluation and inversion that only propose sample locations.
//! Every proposed value is rechecked by numerical admission at 128 bits.

use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::{BSplineCurve, BSplineSurface};

fn real(value: &Real, dec: &Decider) -> Result<f64, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.to_f64_nearest());
    }
    Ok(super::scalar(value, dec)?.to_f64_nearest())
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

// One unit per converted scalar, charged before the conversion starts.
fn charge(scalars: usize, dec: &mut Decider) -> Result<(), Refusal> {
    dec.charge_work(scalars as u64)
}

impl Surface64 {
    pub fn new(source: &BSplineSurface, dec: &mut Decider) -> Result<Self, Refusal> {
        let [u, v] = source.knots();
        charge(u.len() + v.len() + 4 * source.controls().len(), dec)?;
        let columns = source.control_count()[1];
        let mut controls = Vec::with_capacity(source.controls().len());
        for (row, (points, weights)) in source
            .controls()
            .chunks(columns)
            .zip(source.weights().chunks(columns))
            .enumerate()
        {
            if row % 16 == 0 {
                dec.check_deadline()?;
            }
            for (p, w) in points.iter().zip(weights) {
                let w = real(w, dec)?;
                controls.push([
                    real(&p.x, dec)? * w,
                    real(&p.y, dec)? * w,
                    real(&p.z, dec)? * w,
                    w,
                ]);
            }
        }
        Ok(Self {
            degrees: source.degrees(),
            knots: [u, v].map(|row| row.iter().map(|q| q.to_f64_nearest()).collect()),
            count: source.control_count(),
            controls,
        })
    }
    pub fn internal_knots(&self, axis: usize) -> &[f64] {
        &self.knots[axis][self.degrees[axis] + 1..self.count[axis]]
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
    pub fn new(source: &BSplineCurve, dec: &mut Decider) -> Result<Self, Refusal> {
        charge(source.knots().len() + 4 * source.controls().len(), dec)?;
        let mut controls = Vec::with_capacity(source.controls().len());
        for (i, (p, w)) in source.controls().zip(source.weights()).enumerate() {
            if i % 64 == 0 {
                dec.check_deadline()?;
            }
            let w = real(w, dec)?;
            controls.push([
                real(&p.x, dec)? * w,
                real(&p.y, dec)? * w,
                real(&p.z, dec)? * w,
                w,
            ]);
        }
        Ok(Self {
            degree: source.degree(),
            knots: source.knots().iter().map(|q| q.to_f64_nearest()).collect(),
            controls,
        })
    }
    pub fn internal_knots(&self) -> &[f64] {
        &self.knots[self.degree + 1..self.controls.len()]
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
pub(crate) fn interpolate<const D: usize>(
    parameters: &[f64],
    values: &[[f64; D]],
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<[f64; D]>), Refusal> {
    interpolate_degree::<D, 3, 7>(parameters, values, dec)
}

pub(super) fn interpolate_quintic<const D: usize>(
    parameters: &[f64],
    values: &[[f64; D]],
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<[f64; D]>), Refusal> {
    interpolate_degree::<D, 5, 11>(parameters, values, dec)
}

fn interpolate_degree<const D: usize, const P: usize, const W: usize>(
    parameters: &[f64],
    values: &[[f64; D]],
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<[f64; D]>), Refusal> {
    let n = parameters.len();
    assert!(n > P && values.len() == n);
    let skip = P.div_ceil(2);
    let mut knots = vec![parameters[0]; P + 1];
    knots.extend_from_slice(&parameters[skip..n - skip]);
    knots.extend(std::iter::repeat_n(parameters[n - 1], P + 1));
    let width = parameters[n - 1] - parameters[0];
    let chord = |t: f64| -> [f64; D] {
        let share = (t - parameters[0]) / width;
        std::array::from_fn(|c| values[0][c] + share * (values[n - 1][c] - values[0][c]))
    };
    // Band storage: column j of row i lives at slot j + P - i.
    let mut band = vec![[0.0f64; W]; n];
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
        let index = span(&knots, P, n, *t);
        let (basis, _) = basis(&knots, P, index, *t);
        for (k, value) in basis.into_iter().enumerate() {
            let column = index - P + k;
            let slot = column + P - i;
            band[i][slot] = value;
        }
    }
    for k in 0..n {
        dec.charge_work(1)?;
        let pivot = band[k][P];
        if pivot.is_nan() || pivot.abs() <= 1e-300 {
            return Err(Refusal::new(
                "step.numerical-fit",
                "interpolation matrix lost a pivot",
            ));
        }
        for i in k + 1..(k + P + 1).min(n) {
            let factor = band[i][k + P - i] / pivot;
            for j in k..(k + P + 1).min(n) {
                band[i][j + P - i] -= factor * band[k][j + P - k];
            }
            for c in 0..D {
                rhs[i][c] -= factor * rhs[k][c];
            }
        }
    }
    let mut controls = vec![[0.0; D]; n];
    for k in (0..n).rev() {
        let mut value = rhs[k];
        for j in k + 1..(k + P + 1).min(n) {
            for c in 0..D {
                value[c] -= band[k][j + P - k] * controls[j][c];
            }
        }
        for c in 0..D {
            controls[k][c] = value[c] / band[k][P];
        }
    }
    for (i, control) in controls.iter_mut().enumerate() {
        let line = chord(knots[i + 1..=i + P].iter().sum::<f64>() / P as f64);
        for c in 0..D {
            control[c] += line[c];
        }
    }
    controls[0] = values[0];
    controls[n - 1] = values[n - 1];
    Ok((knots, controls))
}

/// A quartic Bezier interpolant through five parameter locations.
pub(super) fn interpolate_quartic<const D: usize>(
    parameters: &[f64; 5],
    values: &[[f64; D]; 5],
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<[f64; D]>), Refusal> {
    let width = parameters[4] - parameters[0];
    let chord = |share: f64| -> [f64; D] {
        std::array::from_fn(|axis| values[0][axis] + share * (values[4][axis] - values[0][axis]))
    };
    let mut matrix = [[0.0; 3]; 3];
    let mut rhs = [[0.0; D]; 3];
    for i in 0..3 {
        dec.charge_work(1)?;
        let t = (parameters[i + 1] - parameters[0]) / width;
        let s = 1.0 - t;
        matrix[i] = [
            4.0 * t * s * s * s,
            6.0 * t * t * s * s,
            4.0 * t * t * t * s,
        ];
        let line = chord(t);
        rhs[i] = std::array::from_fn(|axis| values[i + 1][axis] - line[axis]);
    }
    for k in 0..3 {
        dec.charge_work(1)?;
        let pivot = matrix[k][k];
        if !pivot.is_finite() || pivot.abs() <= 1e-300 {
            return Err(super::fit_failure(
                "quartic interpolation lost a separated pivot",
            ));
        }
        for i in k + 1..3 {
            let factor = matrix[i][k] / pivot;
            for j in k..3 {
                matrix[i][j] -= factor * matrix[k][j];
            }
            for axis in 0..D {
                rhs[i][axis] -= factor * rhs[k][axis];
            }
        }
    }
    let mut controls = [[0.0; D]; 5];
    for k in (0..3).rev() {
        for axis in 0..D {
            let mut value = rhs[k][axis];
            for j in k + 1..3 {
                value -= matrix[k][j] * controls[j + 1][axis];
            }
            controls[k + 1][axis] = value / matrix[k][k];
        }
    }
    for (i, control) in controls.iter_mut().enumerate() {
        let line = chord(i as f64 / 4.0);
        for axis in 0..D {
            control[axis] += line[axis];
        }
    }
    controls[0] = values[0];
    controls[4] = values[4];
    let mut knots = vec![parameters[0]; 5];
    knots.extend(std::iter::repeat_n(parameters[4], 5));
    Ok((knots, controls.to_vec()))
}

/// Joins equal-degree clamped pieces with C0 continuity at their shared ends.
pub(super) fn concatenate<const D: usize>(
    pieces: Vec<(Vec<f64>, Vec<[f64; D]>)>,
) -> (Vec<f64>, Vec<[f64; D]>) {
    let count = pieces.len();
    let mut knots = Vec::new();
    let mut controls = Vec::new();
    for (index, (piece_knots, piece_controls)) in pieces.into_iter().enumerate() {
        let degree = piece_knots.len() - piece_controls.len() - 1;
        let first = if index == 0 { 0 } else { degree + 1 };
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
pub(crate) fn evaluate<const D: usize>(knots: &[f64], controls: &[[f64; D]], t: f64) -> [f64; D] {
    let n = controls.len();
    let degree = knots.len() - n - 1;
    let t = t.clamp(knots[degree], knots[n]);
    let index = span(knots, degree, n, t);
    let (basis, _) = basis(knots, degree, index, t);
    let mut result = [0.0; D];
    for (k, value) in basis.into_iter().enumerate() {
        for c in 0..D {
            result[c] += value * controls[index - degree + k][c];
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Rational, with_operation_budget};
    use catallax_shapes::{Point3, SplineBoundary};

    fn grid(side: usize) -> BSplineSurface {
        let knots: Vec<_> = std::iter::once(0)
            .chain(0..side as i64)
            .chain(std::iter::once(side as i64 - 1))
            .map(Rational::from_i64)
            .collect();
        let controls = (0..side * side)
            .map(|i| {
                Point3::new(
                    Real::from((i / side) as i64),
                    Real::from((i % side) as i64),
                    Real::from(1),
                )
            })
            .collect();
        BSplineSurface::new(
            [1, 1],
            [knots.clone(), knots],
            [side, side],
            controls,
            vec![Real::from(2); side * side],
            &mut Decider::default(),
        )
        .unwrap()
    }
    #[test]
    fn retained_curve_and_surface_proposals_refine_algebraic_controls_and_weights() {
        let offset = Real::from(Rational::parse("1000000000000333333333333").unwrap()).sqrt()
            - Real::from_i64(1_000_000_000_000);
        let point = |x, y| {
            Point3::new(
                &offset + &Real::from_i64(x),
                Real::from_i64(y),
                Real::zero(),
            )
        };
        let mut dec = Decider::default();
        let source = BSplineSurface::new(
            [1, 1],
            std::array::from_fn(|_| {
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                ]
            }),
            [2, 2],
            vec![point(0, 0), point(0, 1), point(1, 0), point(1, 1)],
            (0..4)
                .map(|i| Real::one() + &offset * &Real::from_i64(i))
                .collect(),
            &mut dec,
        )
        .unwrap();
        let surface = Surface64::new(&source, &mut dec).unwrap();
        let boundary = BSplineCurve::boundary(&source, SplineBoundary::UMin);
        let curve = Curve64::new(&boundary, &mut dec).unwrap();
        let mut error = 0.0_f64;
        for [u, v] in [[1, 2], [4, 6], [7, 5]] {
            let u = Rational::from_ratio(u, 8).unwrap();
            let v = Rational::from_ratio(v, 8).unwrap();
            let expected = source
                .evaluate_box([&u, &u], [&v, &v], 128, &dec)
                .unwrap()
                .each_ref()
                .map(|value| value.mid().to_f64_nearest());
            error = error.max(distance(
                surface.point([u.to_f64_nearest(), v.to_f64_nearest()]),
                expected,
            ));
            let expected = boundary
                .evaluate_box([&v, &v], 128, &dec)
                .unwrap()
                .each_ref()
                .map(|value| value.mid().to_f64_nearest());
            error = error.max(distance(curve.point(v.to_f64_nearest()), expected));
        }
        assert!(
            error <= 2.5e-8,
            "retained geometry proposal error={error:e} mm"
        );
        println!(
            "refined proposal error={error:e} mm; coarse offset error={:e} mm",
            (offset.to_f64()
                - super::super::scalar(&offset, &dec)
                    .unwrap()
                    .to_f64_nearest())
            .abs()
        );
    }

    #[test]
    fn binary64_copies_charge_their_payload_before_converting() {
        let side = 40;
        let source = grid(side);
        let payload = (2 * (side + 2) + 4 * side * side) as u64;
        let ((), work) = with_operation_budget(2 * payload, || {
            let copy = Surface64::new(&source, &mut Decider::default()).unwrap();
            assert_eq!(
                copy.point([1.5, 2.5]).partial_cmp(&[1.5, 2.5, 1.0]),
                Some(std::cmp::Ordering::Equal)
            );
        });
        assert!(!work.exhausted && work.used >= payload, "{work:?}");
        let (result, work) = with_operation_budget(payload - 1, || {
            Surface64::new(&source, &mut Decider::default()).map(|_| ())
        });
        assert!(result.is_err() && work.exhausted, "{work:?}");
        let curve = BSplineCurve::boundary(&source, SplineBoundary::UMin);
        let curve_payload = (curve.knots().len() + 4 * curve.controls().len()) as u64;
        let ((), work) = with_operation_budget(2 * curve_payload, || {
            let copy = Curve64::new(&curve, &mut Decider::default()).unwrap();
            assert_eq!(
                copy.point(3.5).partial_cmp(&[0.0, 3.5, 1.0]),
                Some(std::cmp::Ordering::Equal)
            );
        });
        assert!(!work.exhausted && work.used >= curve_payload, "{work:?}");
        let (result, work) = with_operation_budget(curve_payload - 1, || {
            Curve64::new(&curve, &mut Decider::default()).map(|_| ())
        });
        assert!(result.is_err() && work.exhausted, "{work:?}");
    }
}
