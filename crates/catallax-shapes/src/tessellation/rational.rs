use std::{
    cell::RefCell,
    ops::{Add, Div, Mul, Sub},
};

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
    control_hull: [Interval; 4],
    first_hulls: [[Interval; 4]; 2],
    second_hulls: [[Interval; 4]; 3],
    position_first: [Interval; 2],
    hessian: [Interval; 3],
    affine_hessian: [Interval; 3],
    radius: Interval,
    min_weight: Interval,
    coordinate_bounds: [Interval; 3],
    position_row: RefCell<Option<PositionRow>>,
}

struct PositionRow {
    u: u64,
    controls: Net<Interval>,
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
    rational_span_range(surface, 0..=usize::MAX, dec)
}

pub(crate) fn rational_span_range(
    surface: &BSplineSurface,
    range: std::ops::RangeInclusive<usize>,
    dec: &mut Decider,
) -> Result<Vec<RationalSpan>, Refusal> {
    let degree = surface.degrees();
    let count = range.end().saturating_sub(*range.start()).saturating_add(1);
    surface
        .span_indices()
        .skip(*range.start())
        .take(count)
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
    /// Prepares a normalized rational Bézier curve as a one-parameter chart.
    pub fn from_curve(controls: Vec<[Interval; 4]>, dec: &mut Decider) -> Result<Self, Refusal> {
        let degree = controls
            .len()
            .checked_sub(1)
            .ok_or_else(|| refuse("empty curve net"))?;
        if controls
            .iter()
            .any(|point| point.iter().any(|x| !x.is_finite()) || point[3].lo() <= 0.0)
        {
            return Err(refuse(
                "homogeneous controls need finite coordinates and positive weights",
            ));
        }
        Self::prepare(
            std::array::from_fn(|_| [Rational::zero(), Rational::one()]),
            [degree, 0],
            controls,
            dec,
        )
    }

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
            control_hull: control_hull(&control),
            first_hulls: first.each_ref().map(|net| control_hull(net)),
            second_hulls: second.each_ref().map(|net| control_hull(net)),
            control,
            first,
            position_first,
            hessian,
            affine_hessian,
            radius,
            min_weight,
            coordinate_bounds,
            position_row: RefCell::default(),
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
        let center = uv.map(|value| point(value.mid()));
        let delta = std::array::from_fn(|i| uv[i] - center[i]);
        let at = self.homogeneous_position(center.map(|axis| axis.mid()), dec)?;
        let [x, y, z, weight] = expand_center(at, delta, self.control_hull, self.first_hulls, dec)?;
        let weight = Interval::new(weight.lo().max(self.min_weight.lo()), weight.hi())
            .ok_or_else(|| refuse("evaluated weight and control bounds do not overlap"))?;
        Ok(std::array::from_fn(|i| {
            [x, y, z][i] / weight + point(self.centre[i])
        }))
    }

    fn homogeneous_position(
        &self,
        [u, v]: [f64; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 4], Refusal> {
        // On a tensor grid, consecutive points have the same u. Collapse
        // that axis once when the saved tensor products cover the contraction.
        if self.degree[0] < 3 || self.degree[1] == 0 || !(0.0 < u && u < 1.0) {
            return evaluate(&self.control, self.degree, [u, v], point, |i| i, dec);
        }
        let key = u.to_bits();
        let mut row = self.position_row.borrow_mut();
        if row.as_ref().is_none_or(|entry| entry.u != key) {
            if v <= 0.0 || v >= 1.0 {
                return evaluate(&self.control, self.degree, [u, v], point, |i| i, dec);
            }
            let (weights, products) = basis(self.degree[0], u, &point);
            dec.charge_work((products + 4 * self.control.len()) as u64)?;
            let mut controls = vec![[point(0.0); 4]; self.degree[1] + 1];
            for (values, weight) in self.control.chunks_exact(controls.len()).zip(weights) {
                for (column, value) in controls.iter_mut().zip(values) {
                    for (out, component) in column.iter_mut().zip(value) {
                        *out = *out + weight * *component;
                    }
                }
            }
            *row = Some(PositionRow { u: key, controls });
        }
        evaluate(
            &row.as_ref().unwrap().controls,
            [0, self.degree[1]],
            [0.0, v],
            point,
            |i| i,
            dec,
        )
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
        let center = uv.map(|value| point(value.mid()));
        let delta = std::array::from_fn(|i| uv[i] - center[i]);
        let homogeneous = centered_range(
            &self.control,
            self.degree,
            [center, delta],
            self.control_hull,
            self.first_hulls,
            dec,
        )?;
        let mut tangents = [[point(0.0); 3]; 2];
        for (axis, tangent) in tangents.iter_mut().enumerate() {
            if self.first[axis].is_empty() {
                continue;
            }
            let degree = lower_degree(self.degree, axis);
            let first = centered_range(
                &self.first[axis],
                degree,
                [center, delta],
                self.first_hulls[axis],
                [self.second_hulls[axis], self.second_hulls[axis + 1]],
                dec,
            )?;
            dec.charge_work(2 * 3)?;
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
        // A merely acceptable residual can leave an exact rail on alternating
        // sides of its chart line. Continue toward the rounding scale before
        // using the residual to decide whether the projection is admissible.
        let rounding_goal = goal.min(self.radius.hi() * f64::EPSILON);

        // Damped Gauss-Newton (Nocedal and Wright, Numerical Optimization,
        // chapter 10). A bounded span may put the minimum on an edge. Clamp
        // each trial and stop if there is no representable improving step.
        for _ in 0..f64::MANTISSA_DIGITS {
            dec.charge_work(3)?;
            let distance_squared = dot(difference, difference);
            if distance_squared.sqrt() <= rounding_goal {
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
        let homogeneous = evaluate(&self.control, self.degree, uv, |x| x, |x| x.mid(), dec)?;
        Ok(std::array::from_fn(|i| {
            homogeneous[i] / homogeneous[3] - target[i]
        }))
    }

    fn jacobian(&self, uv: [f64; 2], dec: &mut Decider) -> Result<[[f64; 3]; 2], Refusal> {
        let [x, y, z, weight] = evaluate(&self.control, self.degree, uv, |x| x, |x| x.mid(), dec)?;
        let position = [x / weight, y / weight, z / weight];
        let mut derivatives = [[0.0; 3]; 2];
        for (axis, row) in derivatives.iter_mut().enumerate() {
            if self.first[axis].is_empty() {
                continue;
            }
            let degree = lower_degree(self.degree, axis);
            let value = evaluate(&self.first[axis], degree, uv, |x| x, |x| x.mid(), dec)?;
            dec.charge_work(3)?;
            *row = std::array::from_fn(|i| (value[i] - position[i] * value[3]) / weight);
        }
        Ok(derivatives)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn basis<T>(degree: usize, at: f64, scalar: &impl Fn(f64) -> T) -> (Vec<T>, usize)
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T>,
{
    let one = scalar(1.0);
    let t = scalar(at);
    let complement = one - t;
    if degree > 2 {
        // Consecutive Bernstein basis functions have the ratio in Woźny and
        // Chudy (2020), section 2. Reverse the order to divide by max(t, 1-t).
        let reversed = at > 0.5;
        let (major, minor) = if reversed {
            (t, complement)
        } else {
            (complement, t)
        };
        let odds = minor / major;
        let mut power = major;
        let mut exponent = degree;
        let mut weight = one;
        let mut products = 0;
        while exponent != 0 {
            if exponent & 1 == 1 {
                weight = weight * power;
                products += 1;
            }
            exponent >>= 1;
            if exponent != 0 {
                power = power * power;
                products += 1;
            }
        }
        let mut values = Vec::with_capacity(degree + 1);
        values.push(weight);
        for i in 1..=degree {
            let ratio = scalar((degree - i + 1) as f64) / scalar(i as f64);
            weight = weight * (ratio * odds);
            values.push(weight);
        }
        if reversed {
            values.reverse();
        }
        return (values, products + 2 * degree);
    }
    let zero = scalar(0.0);
    let mut values = vec![zero; degree + 1];
    values[0] = one;
    for order in 1..=degree {
        values[order] = values[order - 1] * t;
        for i in (1..order).rev() {
            values[i] = values[i] * complement + values[i - 1] * t;
        }
        values[0] = values[0] * complement;
    }
    (values, degree * (degree + 1))
}

fn evaluate<T>(
    net: &[[Interval; 4]],
    [u, v]: [usize; 2],
    uv: [f64; 2],
    scalar: impl Fn(f64) -> T,
    value: impl Fn(Interval) -> T,
    dec: &mut Decider,
) -> Result<[T; 4], Refusal>
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T>,
{
    let zero = scalar(0.0);
    let one = scalar(1.0);
    let active_basis = |degree, parameter: f64| {
        if degree == 0 || parameter <= 0.0 {
            (0, vec![one], 0)
        } else if parameter >= 1.0 {
            (degree, vec![one], 0)
        } else {
            let (weights, products) = basis(degree, parameter, &scalar);
            (0, weights, products)
        }
    };
    let (first_row, rows, u_products) = active_basis(u, uv[0]);
    let (first_column, columns, v_products) = active_basis(v, uv[1]);
    let products = match (rows.len(), columns.len()) {
        (1, 1) => return Ok(net[first_row * (v + 1) + first_column].map(value)),
        (1, columns) => 4 * columns,
        (rows, 1) => 4 * rows,
        (rows, columns) => 5 * rows * columns,
    };
    dec.charge_work((u_products + v_products + products) as u64)?;
    let mut sum = [zero; 4];
    let row_iter = net.chunks_exact(v + 1).skip(first_row).zip(&rows);
    for (row, &row_weight) in row_iter {
        for (control, &column_weight) in row[first_column..].iter().zip(&columns) {
            let weight = if rows.len() == 1 {
                column_weight
            } else if columns.len() == 1 {
                row_weight
            } else {
                row_weight * column_weight
            };
            for (slot, &control) in sum.iter_mut().zip(control) {
                *slot = *slot + value(control) * weight;
            }
        }
    }
    Ok(sum)
}

fn control_hull(net: &[[Interval; 4]]) -> [Interval; 4] {
    let Some((first, rest)) = net.split_first() else {
        return [point(0.0); 4];
    };
    rest.iter().fold(*first, |bounds, control| {
        std::array::from_fn(|i| bounds[i].hull(&control[i]))
    })
}

fn centered_range(
    net: &[[Interval; 4]],
    degree: [usize; 2],
    [center, [du, dv]]: [[Interval; 2]; 2],
    hull: [Interval; 4],
    [u, v]: [[Interval; 4]; 2],
    dec: &mut Decider,
) -> Result<[Interval; 4], Refusal> {
    // Mean value form on the box. The Bernstein hulls bound the two
    // derivatives along both coordinate displacements from the centre.
    let at = evaluate(
        net,
        degree,
        center.map(|axis| axis.mid()),
        point,
        |i| i,
        dec,
    )?;
    expand_center(at, [du, dv], hull, [u, v], dec)
}

fn expand_center(
    at: [Interval; 4],
    [du, dv]: [Interval; 2],
    hull: [Interval; 4],
    [u, v]: [[Interval; 4]; 2],
    dec: &mut Decider,
) -> Result<[Interval; 4], Refusal> {
    let zero = point(0.0);
    let u_delta = if du == zero {
        [zero; 4]
    } else {
        dec.charge_work(4)?;
        u.map(|d| du * d)
    };
    let v_delta = if dv == zero {
        [zero; 4]
    } else {
        dec.charge_work(4)?;
        v.map(|d| dv * d)
    };
    Ok(intersect_hull(
        std::array::from_fn(|i| at[i] + u_delta[i] + v_delta[i]),
        hull,
    ))
}

fn intersect_hull(evaluated: [Interval; 4], hull: [Interval; 4]) -> [Interval; 4] {
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
    fn bernstein_basis_encloses_cubic_weights_from_both_ends() {
        let exact = [27.0, 27.0, 9.0, 1.0].map(|entry| entry / 64.0);
        for (at, expected) in [
            (0.25, exact),
            (0.75, [exact[3], exact[2], exact[1], exact[0]]),
        ] {
            let (weights, _) = basis(3, at, &point);
            assert!(
                weights
                    .into_iter()
                    .zip(expected)
                    .all(|(weight, value)| weight.contains(value))
            );
        }
    }

    #[test]
    fn contracted_rows_enclose_positions_before_and_after_row_changes() {
        // The controls represent S(u, v) = (u³, v, u³v).
        let controls = (0..4)
            .flat_map(|i| {
                let x = f64::from(u8::from(i == 3));
                [0.0, 1.0].map(move |v| [x, v, x * v, 1.0].map(point))
            })
            .collect();
        let mut dec = Decider::default();
        let span = RationalSpan::prepare(
            std::array::from_fn(|_| [Rational::zero(), Rational::one()]),
            [3, 1],
            controls,
            &mut dec,
        )
        .unwrap();
        for u in [0.25_f64, 0.75, 0.25] {
            let x = u.powi(3);
            for v in [0.25, 0.5, 1.0] {
                let xyz = span.evaluate([u, v], &mut dec).unwrap();
                assert!(
                    xyz.into_iter()
                        .zip([x, v, x * v])
                        .all(|(a, b)| a.contains(b))
                );
            }
        }
    }

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
