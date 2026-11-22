use std::ops::{Add, Mul, Sub};

use catallax_numbers::{Decider, Interval, Rational, Refusal};

use crate::BSplineSurface;

type Net<T> = Vec<[T; 4]>;

/// A knot rectangle whose homogeneous Bézier net uses local unit coordinates.
pub struct RationalSpan {
    domain: [[Rational; 2]; 2],
    degree: [usize; 2],
    centre: [f64; 3],
    control: Net<Interval>,
    first: [Net<Interval>; 2],
    position_first: [Interval; 2],
    hessian: [Interval; 3],
    affine_hessian: [Interval; 3],
    radius: Interval,
    min_weight: Interval,
    coordinate_bounds: [Interval; 3],
}

/// A floating-point chart parameter and an outward bound on its point residual.
pub struct Projection {
    /// Local coordinates on the selected knot span.
    pub uv: [f64; 2],
    /// Distance in world coordinates from the chart evaluation to the input.
    pub residual: Interval,
}

fn refuse(detail: &'static str) -> Refusal {
    Refusal::new("mesh.spline-chart", detail)
}

fn point(value: f64) -> Interval {
    Interval::point(value).expect("validated finite coordinate")
}

fn magnitude(value: Interval) -> Interval {
    point(value.lo().abs().max(value.hi().abs()))
}

fn vector_length(value: [Interval; 3]) -> Interval {
    (value[0].square() + value[1].square() + value[2].square())
        .sqrt()
        .expect("square sum is nonnegative")
}

// The Bernstein derivative net is the degree times the forward differences.
fn difference(net: &[[Interval; 4]], degree: [usize; 2], axis: usize) -> Net<Interval> {
    if net.is_empty() || degree[axis] == 0 {
        return Vec::new();
    }
    let [u, v] = degree;
    let scale = point(degree[axis] as f64);
    let stride = if axis == 0 { v + 1 } else { 1 };
    let mut result = Vec::with_capacity((u + 1) * (v + 1));
    for i in 0..=u - usize::from(axis == 0) {
        for j in 0..=v - usize::from(axis == 1) {
            let index = i * (v + 1) + j;
            let a = net[index];
            let b = net[index + stride];
            result.push(std::array::from_fn(|k| (b[k] - a[k]) * scale));
        }
    }
    result
}

fn lower_degree(mut degree: [usize; 2], axis: usize) -> [usize; 2] {
    degree[axis] = degree[axis].saturating_sub(1);
    degree
}

fn derivative_bound(net: &[[Interval; 4]], radius: Interval, weight: Interval) -> Interval {
    let upper = net
        .iter()
        .map(|&[x, y, z, w]| ((vector_length([x, y, z]) + radius * magnitude(w)) / weight).hi())
        .fold(0.0_f64, f64::max);
    point(upper)
}

/// Prepares the carrier's complete nonzero knot spans in u-major order.
pub fn rational_spans(
    surface: &BSplineSurface,
    dec: &mut Decider,
) -> Result<Vec<RationalSpan>, Refusal> {
    let degree = surface.degrees();
    surface
        .span_indices()
        .map(|span| {
            let domain = std::array::from_fn(|axis| {
                [span[axis], span[axis] + 1].map(|i| surface.knots()[axis][i].clone())
            });
            let exact = surface.span_bezier_controls(span, dec)?;
            let mut control = Vec::with_capacity(exact.len());
            for homogeneous in exact {
                let mut values = [point(0.0); 4];
                for (slot, value) in values.iter_mut().zip(&homogeneous) {
                    *slot = if let Some(rational) = value.rational_ref() {
                        Interval::from_rational(rational)
                    } else {
                        let bound = dec.enclose_at(value, f64::MANTISSA_DIGITS)?;
                        Interval::from_rational(&bound.lo())
                            .hull(&Interval::from_rational(&bound.hi()))
                    };
                }
                if values.iter().any(|value| !value.is_finite()) || values[3].lo() <= 0.0 {
                    return Err(refuse(
                        "homogeneous controls need finite coordinates and positive weights",
                    ));
                }
                control.push(values);
            }
            RationalSpan::prepare(domain, degree, control, dec)
        })
        .collect()
}

impl RationalSpan {
    fn prepare(
        domain: [[Rational; 2]; 2],
        degree: [usize; 2],
        mut control: Net<Interval>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let [x, y, z, w] = control[0];
        let centre = [x, y, z].map(|value| (value / w).mid());
        if centre.iter().any(|value| !value.is_finite()) {
            return Err(refuse("span translation exceeds binary64"));
        }
        let mut radius = 0.0_f64;
        let mut min_weight = f64::INFINITY;
        let mut coordinate_bounds = [x, y, z].map(|value| value / w);
        // Translation (3 products) and the squared translated norm (3).
        dec.charge_work(6 * control.len() as u64)?;
        for homogeneous in &mut control {
            let weight = homogeneous[3];
            min_weight = min_weight.min(weight.lo());
            for i in 0..3 {
                homogeneous[i] = homogeneous[i] - point(centre[i]) * weight;
            }
            let relative = [homogeneous[0], homogeneous[1], homogeneous[2]].map(|h| h / weight);
            radius = radius.max(vector_length(relative).hi());
            for i in 0..3 {
                coordinate_bounds[i] = coordinate_bounds[i].hull(&(relative[i] + point(centre[i])));
            }
        }

        let first = [
            difference(&control, degree, 0),
            difference(&control, degree, 1),
        ];
        let axes = [[0, 0], [0, 1], [1, 1]];
        let second = axes.map(|[a, b]| difference(&first[a], lower_degree(degree, a), b));
        // Four scalings per homogeneous difference. A bound adds three
        // squared components and one product with the radius.
        dec.charge_work(
            4 * first.iter().map(Vec::len).sum::<usize>() as u64
                + 8 * second.iter().map(Vec::len).sum::<usize>() as u64,
        )?;
        let radius = point(radius);
        let min_weight = point(min_weight);
        let hessian = second
            .each_ref()
            .map(|net| derivative_bound(net, radius, min_weight));
        dec.charge_work(4 * first.iter().map(Vec::len).sum::<usize>() as u64 + 2 * 3)?;
        let position_first = first
            .each_ref()
            .map(|net| derivative_bound(net, radius, min_weight));
        let weight_first = first.each_ref().map(|net| {
            point(
                net.iter()
                    .map(|v| magnitude(v[3]).hi())
                    .fold(0.0_f64, f64::max),
            )
        });
        // Apply the quotient rule to N = S * W, using the translated net.
        let affine_hessian = std::array::from_fn(|i| {
            let [a, b] = axes[i];
            hessian[i]
                + (position_first[a] * weight_first[b] + position_first[b] * weight_first[a])
                    / min_weight
        });
        if hessian
            .iter()
            .chain(&affine_hessian)
            .any(|value| !value.is_finite())
        {
            return Err(refuse("span derivative bound exceeds binary64"));
        }
        Ok(Self {
            domain,
            degree,
            centre,
            control,
            first,
            position_first,
            hessian,
            affine_hessian,
            radius,
            min_weight,
            coordinate_bounds,
        })
    }

    /// The rectangle in the original surface parameters.
    pub fn domain(&self) -> &[[Rational; 2]; 2] {
        &self.domain
    }

    /// Homogeneous derivative bounds for use with the chord formula.
    ///
    /// Zheng and Sederberg, Theorem 4. The positive denominator is bounded
    /// below by the smallest weight. Translating the net tightens the radius.
    /// Using the full radius in place of `r - eps` gives a bound independent
    /// of the requested tolerance. Positive vertex weights reparameterize the
    /// flat triangle without changing the set of points it covers.
    pub fn chord_coefficients(&self) -> [Interval; 3] {
        self.hessian
    }

    /// Cartesian control hull enclosing the whole span.
    pub fn coordinate_bounds(&self) -> [Interval; 3] {
        self.coordinate_bounds
    }

    /// Bounds on the speed along the two local span coordinates.
    pub fn first_derivative_bounds(&self) -> [Interval; 2] {
        self.position_first
    }

    /// Cartesian second derivative bounds for affine parameter interpolation.
    pub fn affine_hessian(&self) -> [Interval; 3] {
        self.affine_hessian
    }

    /// Bounds the jump in the transverse first derivative across a common knot.
    ///
    /// The two spans must share homogeneous boundary controls and the same
    /// transverse parameter interval. Widths map unit-span derivatives into
    /// the caller's coordinates. Bernstein differences bound the quotient
    /// derivative on the entire common seam.
    pub fn first_jump(
        &self,
        next: &Self,
        axis: usize,
        widths: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        if axis > 1 {
            return Err(refuse("a seam axis must be 0 or 1"));
        }
        let secondary = 1 - axis;
        if self.degree != next.degree
            || self.domain[secondary] != next.domain[secondary]
            || self.domain[axis][1] != next.domain[axis][0]
            || widths.iter().any(|width| width.lo() <= 0.0)
        {
            return Err(refuse(
                "derivative jump requires adjacent spans and positive widths",
            ));
        }
        if self.first[axis].is_empty() {
            return Ok(point(0.0));
        }
        let degree = lower_degree(self.degree, axis);
        let v_count = degree[1] + 1;
        let delta_centre: [Interval; 3] =
            std::array::from_fn(|i| point(next.centre[i]) - point(self.centre[i]));
        // Translate the next span's first derivative net to this span's origin.
        // Per sample there are three translation, three squared norm, and
        // one radius product.
        dec.charge_work(7 * (degree[secondary] as u64 + 1))?;
        let mut bound = 0.0_f64;
        for i in 0..=degree[secondary] {
            let [end, start] = if axis == 0 {
                [degree[0] * v_count + i, i]
            } else {
                [i * v_count + degree[1], i * v_count]
            };
            let before = self.first[axis][end].map(|value| value / widths[0]);
            let after = next.first[axis][start].map(|value| value / widths[1]);
            let xyz = std::array::from_fn(|k| after[k] + delta_centre[k] * after[3] - before[k]);
            let weight = magnitude(after[3] - before[3]);
            let jump = (vector_length(xyz) + self.radius * weight) / self.min_weight;
            bound = bound.max(jump.hi());
        }
        Ok(point(bound))
    }

    /// Outward-enclosed point at local coordinates on the unit knot rectangle.
    pub fn evaluate(&self, uv: [f64; 2], dec: &mut Decider) -> Result<[Interval; 3], Refusal> {
        if uv.iter().any(|x| !(0.0..=1.0).contains(x)) {
            return Err(refuse("span evaluation leaves its unit rectangle"));
        }
        self.evaluate_box(uv.map(point), dec)
    }

    /// Outward-enclosed positions over a cell in the unit knot rectangle.
    pub fn evaluate_box(
        &self,
        uv: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        if uv.iter().any(|axis| axis.lo() < 0.0 || axis.hi() > 1.0) {
            return Err(refuse("span evaluation leaves its unit rectangle"));
        }
        dec.charge_work(eval_products(self.degree))?;
        let [x, y, z, weight] = evaluate_range(&self.control, self.degree, uv);
        let weight = Interval::new(weight.lo().max(self.min_weight.lo()), weight.hi())
            .ok_or_else(|| refuse("evaluated weight and control bounds do not overlap"))?;
        Ok(std::array::from_fn(|i| {
            [x, y, z][i] / weight + point(self.centre[i])
        }))
    }

    /// Encloses a positive multiple of the parameter normal over a cell.
    ///
    /// The rational hodograph multiplies each first derivative by the square
    /// of the positive weight. Its cross product retains the normal direction
    /// without a division. See Saito, Wang and Sederberg (1995).
    pub fn normal_box(
        &self,
        uv: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        if uv.iter().any(|axis| axis.lo() < 0.0 || axis.hi() > 1.0) {
            return Err(refuse("normal enclosure leaves its knot span"));
        }
        dec.charge_work(eval_products(self.degree))?;
        let homogeneous = evaluate_range(&self.control, self.degree, uv);
        let mut tangents = [[point(0.0); 3]; 2];
        for (axis, tangent) in tangents.iter_mut().enumerate() {
            if self.first[axis].is_empty() {
                continue;
            }
            let degree = lower_degree(self.degree, axis);
            dec.charge_work(eval_products(degree) + 2 * 3)?;
            let first = evaluate_range(&self.first[axis], degree, uv);
            *tangent =
                std::array::from_fn(|i| homogeneous[3] * first[i] - homogeneous[i] * first[3]);
        }
        dec.charge_work(2 * 3)?;
        let [u, v] = tangents;
        Ok([
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ])
    }

    /// Projects from a neighbouring sample's local parameter, or the centre
    /// of this span when starting a new edge. The final point residual is
    /// enclosed independently of the solver. It can exceed the requested goal.
    #[allow(clippy::float_cmp)] // Exact equality detects a step smaller than one representable increment.
    pub fn project(
        &self,
        target: [f64; 3],
        seed: Option<[f64; 2]>,
        goal: f64,
        dec: &mut Decider,
    ) -> Result<Projection, Refusal> {
        if !goal.is_finite() || goal < 0.0 || target.iter().any(|value| !value.is_finite()) {
            return Err(refuse(
                "projection requires a finite point and residual goal",
            ));
        }
        let mut uv = seed.unwrap_or([0.5; 2]);
        if uv.iter().any(|value| !value.is_finite()) {
            return Err(refuse("projection seed exceeds binary64"));
        }
        uv = uv.map(|value| value.clamp(0.0, 1.0));
        let relative = std::array::from_fn(|i| target[i] - self.centre[i]);
        let mut difference = self.offset(uv, relative, dec)?;

        // Damped Gauss-Newton (Nocedal and Wright, Numerical Optimization,
        // chapter 10). A bounded span may put the minimum on an edge. Clamp
        // each trial and stop if there is no representable improving step.
        for _ in 0..f64::MANTISSA_DIGITS {
            dec.charge_work(3)?;
            let distance_squared = dot(difference, difference);
            if distance_squared.sqrt() <= goal {
                break;
            }
            let [du, dv] = self.jacobian(uv, dec)?;
            // Five dot products, the determinant and the two numerators.
            dec.charge_work(5 * 3 + 2 + 2 * 2)?;
            let (a, b, c) = (dot(du, du), dot(du, dv), dot(dv, dv));
            let (d, e) = (dot(du, difference), dot(dv, difference));
            let determinant = a * c - b * b;
            if !determinant.is_finite() || determinant <= 0.0 {
                break;
            }
            let step = [(c * d - b * e) / determinant, (a * e - b * d) / determinant];
            let mut improved = false;
            let mut damping = 1.0;
            for _ in 0..f64::MANTISSA_DIGITS {
                dec.charge_work(2)?;
                let trial = std::array::from_fn(|i| (uv[i] - damping * step[i]).clamp(0.0, 1.0));
                if trial == uv || trial.iter().any(|value| !value.is_finite()) {
                    break;
                }
                let next = self.offset(trial, relative, dec)?;
                dec.charge_work(4)?;
                if dot(next, next) < distance_squared {
                    uv = trial;
                    difference = next;
                    improved = true;
                    break;
                }
                damping *= 0.5;
            }
            if !improved {
                break;
            }
        }

        let sample = self.evaluate(uv, dec)?;
        dec.charge_work(3)?;
        let residual = vector_length(std::array::from_fn(|i| sample[i] - point(target[i])));
        Ok(Projection { uv, residual })
    }

    fn offset(
        &self,
        uv: [f64; 2],
        target: [f64; 3],
        dec: &mut Decider,
    ) -> Result<[f64; 3], Refusal> {
        dec.charge_work(eval_products(self.degree))?;
        let homogeneous = evaluate(&self.control, self.degree, uv, |x| x, |x| x.mid());
        Ok(std::array::from_fn(|i| {
            homogeneous[i] / homogeneous[3] - target[i]
        }))
    }

    fn jacobian(&self, uv: [f64; 2], dec: &mut Decider) -> Result<[[f64; 3]; 2], Refusal> {
        dec.charge_work(eval_products(self.degree))?;
        let [x, y, z, weight] = evaluate(&self.control, self.degree, uv, |x| x, |x| x.mid());
        let position = [x / weight, y / weight, z / weight];
        let mut derivatives = [[0.0; 3]; 2];
        for (axis, row) in derivatives.iter_mut().enumerate() {
            if self.first[axis].is_empty() {
                continue;
            }
            let degree = lower_degree(self.degree, axis);
            dec.charge_work(eval_products(degree) + 3)?;
            let value = evaluate(&self.first[axis], degree, uv, |x| x, |x| x.mid());
            *row = std::array::from_fn(|i| (value[i] - position[i] * value[3]) / weight);
        }
        Ok(derivatives)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// Bernstein basis recurrence. Each order contributes two products per basis
// entry, then the tensor evaluation pays one product for the basis and four
// for the homogeneous coordinates of each control.
fn eval_products([u, v]: [usize; 2]) -> u64 {
    (u * (u + 1) + v * (v + 1) + 5 * (u + 1) * (v + 1)) as u64
}

fn basis<T>(degree: usize, t: T, one: T, zero: T) -> Vec<T>
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T>,
{
    let mut values = vec![zero; degree + 1];
    values[0] = one;
    let complement = one - t;
    for order in 1..=degree {
        values[order] = values[order - 1] * t;
        for i in (1..order).rev() {
            values[i] = values[i] * complement + values[i - 1] * t;
        }
        values[0] = values[0] * complement;
    }
    values
}

fn evaluate<T>(
    net: &[[Interval; 4]],
    [u, v]: [usize; 2],
    uv: [T; 2],
    scalar: impl Fn(f64) -> T,
    value: impl Fn(Interval) -> T,
) -> [T; 4]
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T>,
{
    let zero = scalar(0.0);
    let one = scalar(1.0);
    let rows = basis(u, uv[0], one, zero);
    let columns = basis(v, uv[1], one, zero);
    let mut sum = [zero; 4];
    for (row, &row_weight) in net.chunks_exact(v + 1).zip(&rows) {
        for (control, &column_weight) in row.iter().zip(&columns) {
            let weight = row_weight * column_weight;
            for (slot, &control) in sum.iter_mut().zip(control) {
                *slot = *slot + value(control) * weight;
            }
        }
    }
    sum
}

fn evaluate_range(net: &[[Interval; 4]], degree: [usize; 2], uv: [Interval; 2]) -> [Interval; 4] {
    let evaluated = evaluate(net, degree, uv, point, |i| i);
    let hull = net.iter().skip(1).fold(net[0], |bounds, control| {
        std::array::from_fn(|i| bounds[i].hull(&control[i]))
    });
    std::array::from_fn(|i| {
        // The whole patch is also inside its Bernstein control hull. Preserve
        // exact constants when wide parameter intervals obscure partition of unity.
        let lo = evaluated[i].lo().max(hull[i].lo());
        let hi = evaluated[i].hi().min(hull[i].hi());
        Interval::new(lo, hi).expect("point and control hull enclose the same patch")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bilinear_chart_bounds_triangle_interior_and_projects_an_exact_point() {
        // S(u, v) = (u, v, uv / 4). The triangle over (0,0), (1,0), (0,1)
        // is on z = 0. Its largest error is 1/16, halfway along the diagonal.
        let control = [
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.25, 1.0],
        ];
        let mut dec = Decider::default();
        let span = RationalSpan::prepare(
            std::array::from_fn(|_| [Rational::zero(), Rational::one()]),
            [1, 1],
            control.map(|h| h.map(point)).to_vec(),
            &mut dec,
        )
        .unwrap();
        let bound =
            crate::tessellation::chord_bound_interval([point(1.0); 2], &span.chord_coefficients());
        assert!(bound.hi() >= 1.0 / 16.0);
        assert!(bound.hi() < 1.0 / 16.0 + f64::EPSILON.sqrt());

        // A residual goal smaller than the distance between any of these
        // sample points prevents early acceptance of the seed.
        let uv = [0.25, 0.75];
        let target = [uv[0], uv[1], 3.0 / 64.0];
        let goal = f64::EPSILON.sqrt();
        let projection = span
            .project(target, Some([0.125, 0.875]), goal, &mut dec)
            .unwrap();
        assert!(projection.residual.hi() <= goal);
        for (actual, expected) in projection.uv.into_iter().zip(uv) {
            assert!((actual - expected).abs() <= goal);
        }
    }
}
