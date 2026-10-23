//! Finite rational tensor-product B-spline definitions and certified chart jets.
//!
//! Controls remain exact definitions. A represented face additionally needs
//! regularity, joins, finite trims and its correspondence certificate.

use std::sync::Arc;

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};

use crate::{Point3, Pose, math};

mod curve;
pub use curve::{BSplineCurve, BSplineCurveJet, SplineBoundary, SplineCurveParameter};

fn validate_axis(
    degree: usize,
    knots: &[Rational],
    n: usize,
    dec: &Decider,
) -> Result<(), Refusal> {
    if degree == 0
        || degree >= n
        || n.checked_add(degree).and_then(|x| x.checked_add(1)) != Some(knots.len())
    {
        return Err(invalid("invalid degree or expanded knot count"));
    }
    for pair in knots.windows(2) {
        dec.check_deadline()?;
        if pair[0] > pair[1] {
            return Err(invalid("knot vector is not nondecreasing"));
        }
    }
    if knots[degree] >= knots[n] {
        return Err(invalid("zero-width active knot domain"));
    }
    if !knots[..=degree].iter().all(|x| x == &knots[degree])
        || !knots[n..].iter().all(|x| x == &knots[n])
    {
        return Err(invalid("finite spline endpoints must be clamped"));
    }
    if knots[degree + 1] == knots[degree] || knots[n - 1] == knots[n] {
        return Err(invalid(
            "endpoint knot multiplicity exceeds degree plus one",
        ));
    }
    let mut multiplicity = 0;
    let mut previous = None;
    for knot in &knots[degree + 1..n] {
        multiplicity = if previous == Some(knot) {
            multiplicity + 1
        } else {
            1
        };
        if multiplicity > degree {
            return Err(invalid("discontinuous internal knot"));
        }
        previous = Some(knot);
    }
    Ok(())
}

fn validate_controls(
    controls: &[Point3],
    weights: &[Real],
    dec: &mut Decider,
) -> Result<(Rational, bool), Refusal> {
    let mut minimum: Option<Rational> = None;
    for (point, weight) in controls.iter().zip(weights) {
        dec.check_deadline()?;
        math::validate_point(point, dec)?;
        math::validate(weight, dec)?;
        math::positive(weight, "shapes.bspline-weight", dec)?;
        let lower = positive_weight_bound(weight, dec)?;
        minimum = Some(match minimum {
            Some(old) => old.min(lower),
            None => lower,
        });
    }
    let constant = weights[1..].iter().all(|w| w.same_definition(&weights[0]));
    Ok((minimum.expect("validated nonempty control net"), constant))
}

fn invalid(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.bspline-data", detail)
}
fn bounds(lo: &Rational, hi: &Rational, precision: u32) -> Ball {
    Ball::from_rational(lo, precision).hull(&Ball::from_rational(hi, precision))
}
fn zero(precision: u32) -> Ball {
    Ball::from_rational(&Rational::zero(), precision)
}

fn positive_weight_bound(weight: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    if let Some(value) = weight.rational_ref() {
        return Ok(value.clone());
    }
    let maximum = dec.policy().max_precision_bits;
    if maximum == 0 {
        return Err(Refusal::new(
            "numbers.precision-limit",
            "no spline weight precision available",
        ));
    }
    let mut precision = maximum.min(64);
    loop {
        let lower = dec.enclose_budgeted_at(weight, precision)?.lo();
        if lower.sign() == Sign::Positive {
            return Ok(lower);
        }
        if precision == maximum {
            return Err(Refusal::new(
                "numbers.geometric-gap-unproven",
                "positive spline weight lower bound unresolved within caller policy",
            ));
        }
        precision = precision.saturating_mul(2).min(maximum);
    }
}

#[derive(Clone, Debug)]
struct Net {
    degrees: [usize; 2],
    knots: [Arc<[Rational]>; 2],
    count: [usize; 2],
    controls: Vec<[Real; 4]>,
}

impl Net {
    fn derivative(&self, axis: usize) -> Option<Self> {
        let degree = self.degrees[axis];
        if degree == 0 {
            return None;
        }
        let mut count = self.count;
        count[axis] -= 1;
        let scale = Real::from_rational(Rational::from_integer((degree as u64).into()));
        let mut controls = Vec::with_capacity(count[0] * count[1]);
        for u in 0..count[0] {
            for v in 0..count[1] {
                let i = if axis == 0 { u } else { v };
                let denominator = &self.knots[axis][i + degree + 1] - &self.knots[axis][i + 1];
                let a = &self.controls[u * self.count[1] + v];
                let b = &self.controls
                    [(u + usize::from(axis == 0)) * self.count[1] + v + usize::from(axis == 1)];
                controls.push(if denominator.is_zero() {
                    // A repeated-knot zero-support basis contributes nothing.
                    std::array::from_fn(|_| Real::zero())
                } else {
                    let factor = &scale / &Real::from_rational(denominator);
                    std::array::from_fn(|k| &(&b[k] - &a[k]) * &factor)
                });
            }
        }
        let mut degrees = self.degrees;
        degrees[axis] -= 1;
        let mut knots = self.knots.clone();
        knots[axis] = Arc::from(&self.knots[axis][1..self.knots[axis].len() - 1]);
        Some(Self {
            degrees,
            knots,
            count,
            controls,
        })
    }

    fn spans(&self, axis: usize) -> impl Iterator<Item = usize> + '_ {
        (self.degrees[axis]..self.count[axis])
            .filter(move |&i| self.knots[axis][i] < self.knots[axis][i + 1])
    }

    fn evaluate_span(
        &self,
        span: [usize; 2],
        parameter: [&Ball; 2],
        precision: u32,
        dec: &Decider,
        workspace: &mut [Vec<[Ball; 4]>; 2],
    ) -> Result<[Ball; 4], Refusal> {
        let [row, columns] = workspace;
        columns.clear();
        for j in span[1] - self.degrees[1]..=span[1] {
            row.clear();
            for i in span[0] - self.degrees[0]..=span[0] {
                dec.check_deadline()?;
                let point = &self.controls[i * self.count[1] + j];
                row.push([
                    dec.enclose_at(&point[0], precision)?,
                    dec.enclose_at(&point[1], precision)?,
                    dec.enclose_at(&point[2], precision)?,
                    dec.enclose_at(&point[3], precision)?,
                ]);
            }
            columns.push(de_boor(
                row,
                &self.knots[0],
                span[0],
                parameter[0],
                precision,
                dec,
            )?);
        }
        de_boor(
            columns,
            &self.knots[1],
            span[1],
            parameter[1],
            precision,
            dec,
        )
    }

    fn enclose(
        &self,
        parameter: [&Ball; 2],
        precision: u32,
        dec: &Decider,
        workspace: &mut [Vec<[Ball; 4]>; 2],
    ) -> Result<[Ball; 4], Refusal> {
        let input = parameter.map(|t| (t.lo(), t.hi()));
        let mut result: Option<[Ball; 4]> = None;
        for u in self.spans(0) {
            let ulo = input[0].0.clone().max(self.knots[0][u].clone());
            let uhi = input[0].1.clone().min(self.knots[0][u + 1].clone());
            if ulo > uhi {
                continue;
            }
            let u_parameter = bounds(&ulo, &uhi, precision);
            for v in self.spans(1) {
                let vlo = input[1].0.clone().max(self.knots[1][v].clone());
                let vhi = input[1].1.clone().min(self.knots[1][v + 1].clone());
                if vlo > vhi {
                    continue;
                }
                let v_parameter = bounds(&vlo, &vhi, precision);
                let value = self.evaluate_span(
                    [u, v],
                    [&u_parameter, &v_parameter],
                    precision,
                    dec,
                    workspace,
                )?;
                match &mut result {
                    None => result = Some(value),
                    Some(previous) => {
                        for (old, new) in previous.iter_mut().zip(value) {
                            *old = old.hull(&new);
                        }
                    }
                }
            }
        }
        result.ok_or_else(|| invalid("parameter box meets no nonzero knot span"))
    }
}

fn de_boor(
    values: &mut Vec<[Ball; 4]>,
    knots: &[Rational],
    span: usize,
    parameter: &Ball,
    precision: u32,
    dec: &Decider,
) -> Result<[Ball; 4], Refusal> {
    let degree = values.len() - 1;
    let one = Ball::from_rational(&Rational::one(), precision);
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            dec.check_deadline()?;
            let i = span - degree + j;
            let denominator = &knots[i + degree + 1 - r] - &knots[i];
            if denominator.sign() != Sign::Positive {
                return Err(invalid("zero-support interval in de Boor evaluation"));
            }
            let alpha = (parameter - &Ball::from_rational(&knots[i], precision))
                .checked_div(&Ball::from_rational(&denominator, precision))
                .ok_or_else(|| invalid("knot denominator was not separated from zero"))?;
            let complement = &one - &alpha;
            let (before, current) = values.split_at_mut(j);
            for (value, previous) in current[0].iter_mut().zip(&before[j - 1]) {
                *value = &complement * previous + &alpha * &*value;
            }
        }
    }
    Ok(values.pop().expect("nonempty de Boor control row"))
}

#[derive(Debug)]
struct Data {
    degrees: [usize; 2],
    knots: [Arc<[Rational]>; 2],
    count: [usize; 2],
    controls: Vec<Point3>,
    weights: Arc<[Real]>,
    weight_lower: Rational,
    constant_weight: bool,
    // Value, u, v, uu, uv and vv homogeneous derivative nets.
    jets: [Option<Net>; 6],
}

/// Exact definition of a finite rational tensor-product B-spline surface.
///
/// Positive weights prove a positive rational denominator over every span.
/// Construction validates the representation, not regularity or injectivity.
/// Topology must also prove joins, orientation and finite-domain correspondence.
#[derive(Clone, Debug)]
pub struct BSplineSurface {
    data: Arc<Data>,
}

/// Outward chart enclosures, with one-sided hulls at internal knots.
#[derive(Clone, Debug)]
pub struct BSplineJet {
    /// Position in millimetres.
    pub value: [Ball; 3],
    /// First derivatives in u, v order.
    pub first: [[Ball; 3]; 2],
    /// Second derivatives in uu, uv, vv order.
    pub second: [[Ball; 3]; 3],
}

impl BSplineSurface {
    /// Validates degrees, expanded knots, rectangular controls and positive weights.
    /// Controls and weights are in u-major, v-minor order.
    pub fn new(
        degrees: [usize; 2],
        knots: [Vec<Rational>; 2],
        count: [usize; 2],
        controls: Vec<Point3>,
        weights: Vec<Real>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let size = count[0]
            .checked_mul(count[1])
            .ok_or_else(|| invalid("control-net size overflow"))?;
        if controls.len() != size || weights.len() != size {
            return Err(invalid(
                "control and weight counts do not match the rectangular net",
            ));
        }
        for axis in 0..2 {
            validate_axis(degrees[axis], &knots[axis], count[axis], dec)?;
        }
        let (weight_lower, constant_weight) = validate_controls(&controls, &weights, dec)?;
        Ok(Self::from_validated(
            degrees,
            knots.map(Arc::from),
            count,
            controls,
            weights.into(),
            weight_lower,
            constant_weight,
        ))
    }

    fn from_validated(
        degrees: [usize; 2],
        knots: [Arc<[Rational]>; 2],
        count: [usize; 2],
        controls: Vec<Point3>,
        weights: Arc<[Real]>,
        weight_lower: Rational,
        constant_weight: bool,
    ) -> Self {
        let net = Net {
            degrees,
            knots: knots.clone(),
            count,
            controls: controls
                .iter()
                .zip(weights.iter())
                .map(|(p, w)| [&p.x * w, &p.y * w, &p.z * w, w.clone()])
                .collect(),
        };
        let u = net.derivative(0);
        let v = net.derivative(1);
        let uu = u.as_ref().and_then(|n| n.derivative(0));
        let uv = u.as_ref().and_then(|n| n.derivative(1));
        let vv = v.as_ref().and_then(|n| n.derivative(1));
        Self {
            data: Arc::new(Data {
                degrees,
                knots,
                count,
                controls,
                weights,
                weight_lower,
                constant_weight,
                jets: [Some(net), u, v, uu, uv, vv],
            }),
        }
    }

    /// Tensor degrees, with no inference from control-point count.
    pub fn degrees(&self) -> [usize; 2] {
        self.data.degrees
    }
    /// Rectangular control counts in u, v order.
    pub fn control_count(&self) -> [usize; 2] {
        self.data.count
    }
    /// Exact expanded knot vectors.
    pub fn knots(&self) -> [&[Rational]; 2] {
        self.data.knots.each_ref().map(AsRef::as_ref)
    }
    /// Exact controls in u-major, v-minor order.
    pub fn controls(&self) -> &[Point3] {
        &self.data.controls
    }
    /// Exact positive weights, one for a polynomial carrier.
    pub fn weights(&self) -> &[Real] {
        &self.data.weights
    }
    /// Cheap positive definition identity, not a negative geometry comparison.
    pub fn same_definition(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
            || (self.data.degrees == other.data.degrees
                && self.data.knots == other.data.knots
                && self.data.count == other.data.count
                && self.controls().iter().zip(other.controls()).all(|(a, b)| {
                    [&a.x, &a.y, &a.z]
                        .into_iter()
                        .zip([&b.x, &b.y, &b.z])
                        .all(|(a, b)| a.same_definition(b))
                })
                && self
                    .weights()
                    .iter()
                    .zip(other.weights())
                    .all(|(a, b)| a.same_definition(b)))
    }

    /// Whether every control and weight has a cached exact rational value.
    pub fn has_rational_definition(&self) -> bool {
        self.controls()
            .iter()
            .all(|p| p.rational_components().is_some())
            && self.weights().iter().all(|w| w.rational_ref().is_some())
    }
    /// Closed active parameter intervals.
    pub fn domain(&self) -> [[&Rational; 2]; 2] {
        std::array::from_fn(|axis| {
            [
                &self.data.knots[axis][self.data.degrees[axis]],
                &self.data.knots[axis][self.data.count[axis]],
            ]
        })
    }
    /// Ordered nonzero span intervals in the selected parameter direction.
    pub fn spans(&self, axis: usize) -> impl Iterator<Item = [&Rational; 2]> {
        self.data.jets[0]
            .as_ref()
            .unwrap()
            .spans(axis)
            .map(move |i| [&self.data.knots[axis][i], &self.data.knots[axis][i + 1]])
    }

    /// Applies an exact rigid motion, sharing unchanged knots and scalar weights.
    pub fn transformed(&self, pose: &Pose) -> Self {
        Self::from_validated(
            self.data.degrees,
            self.data.knots.clone(),
            self.data.count,
            self.data
                .controls
                .iter()
                .map(|p| pose.apply_point(p))
                .collect(),
            self.data.weights.clone(),
            self.data.weight_lower.clone(),
            self.data.constant_weight,
        )
    }

    fn check_binary_box(&self, u: &Ball, v: &Ball, dec: &Decider) -> Result<u32, Refusal> {
        dec.check_deadline()?;
        let precision = u.precision().max(v.precision());
        if precision > dec.policy().max_precision_bits {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "spline chart precision exceeds caller policy",
            ));
        }
        for (parameter, domain) in [u, v].into_iter().zip(self.domain()) {
            if parameter.lo() < *domain[0] || parameter.hi() > *domain[1] {
                return Err(Refusal::new(
                    "shapes.bspline-domain",
                    "chart box leaves the finite spline domain",
                ));
            }
        }
        Ok(precision)
    }

    /// Encloses position without evaluating unrequested derivatives.
    pub fn evaluate(&self, u: &Ball, v: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        let precision = self.check_binary_box(u, v, dec)?;
        let ([value], denominator) = self.homogeneous_cover(u, v, precision, dec)?;
        cartesian(&value, &denominator)
    }

    /// Encloses both first derivatives, without evaluating second derivatives.
    pub fn derivatives(
        &self,
        u: &Ball,
        v: &Ball,
        dec: &Decider,
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        let precision = self.check_binary_box(u, v, dec)?;
        let ([position, du, dv], denominator) = self.homogeneous_cover(u, v, precision, dec)?;
        if self.data.constant_weight {
            return Ok([cartesian(&du, &denominator)?, cartesian(&dv, &denominator)?]);
        }
        let value = cartesian(&position, &denominator)?;
        Ok([
            quotient_first(&value, &du, &denominator)?,
            quotient_first(&value, &dv, &denominator)?,
        ])
    }

    /// Encloses position and all first/second chart derivatives on a finite box.
    /// Boxes crossing knots are covered by every intersected nonzero span.
    pub fn jet(&self, u: &Ball, v: &Ball, dec: &Decider) -> Result<BSplineJet, Refusal> {
        let precision = self.check_binary_box(u, v, dec)?;
        self.jet_cover(u, v, precision, dec)
    }

    /// Encloses an exact rational parameter box, including non-dyadic finite
    /// endpoints. Outward rounding does not turn a valid endpoint into an
    /// out-of-domain request; only the declared rational box is being covered.
    pub fn jet_box(
        &self,
        u: [&Rational; 2],
        v: [&Rational; 2],
        precision: u32,
        dec: &Decider,
    ) -> Result<BSplineJet, Refusal> {
        dec.check_deadline()?;
        if precision == 0 || precision > dec.policy().max_precision_bits {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "spline chart precision exceeds caller policy",
            ));
        }
        for (parameter, domain) in [u, v].into_iter().zip(self.domain()) {
            if parameter[0] > parameter[1] || parameter[0] < domain[0] || parameter[1] > domain[1] {
                return Err(Refusal::new(
                    "shapes.bspline-domain",
                    "rational chart box leaves the finite spline domain",
                ));
            }
        }
        self.jet_cover(
            &bounds(u[0], u[1], precision),
            &bounds(v[0], v[1], precision),
            precision,
            dec,
        )
    }

    fn jet_cover(
        &self,
        u: &Ball,
        v: &Ball,
        precision: u32,
        dec: &Decider,
    ) -> Result<BSplineJet, Refusal> {
        let (homogeneous, denominator) = self.homogeneous_cover::<6>(u, v, precision, dec)?;
        let divide = |value: Ball| {
            value
                .checked_div(&denominator)
                .ok_or_else(|| invalid("positive rational chart denominator was not separated"))
        };
        let value = [
            divide(homogeneous[0][0].clone())?,
            divide(homogeneous[0][1].clone())?,
            divide(homogeneous[0][2].clone())?,
        ];
        let mut first: [[Ball; 3]; 2] =
            std::array::from_fn(|_| std::array::from_fn(|_| zero(precision)));
        for axis in 0..2 {
            for k in 0..3 {
                let correction = &value[k] * &homogeneous[axis + 1][3];
                first[axis][k] = divide(&homogeneous[axis + 1][k] - &correction)?;
            }
        }
        let mut second: [[Ball; 3]; 3] =
            std::array::from_fn(|_| std::array::from_fn(|_| zero(precision)));
        for (output, (a, b)) in [(0, 0), (0, 1), (1, 1)].into_iter().enumerate() {
            for k in 0..3 {
                let correction = &first[a][k] * &homogeneous[b + 1][3]
                    + &first[b][k] * &homogeneous[a + 1][3]
                    + &value[k] * &homogeneous[output + 3][3];
                second[output][k] = divide(&homogeneous[output + 3][k] - &correction)?;
            }
        }
        Ok(BSplineJet {
            value,
            first,
            second,
        })
    }

    fn homogeneous_cover<const N: usize>(
        &self,
        u: &Ball,
        v: &Ball,
        precision: u32,
        dec: &Decider,
    ) -> Result<([[Ball; 4]; N], Ball), Refusal> {
        let mut values: [Option<[Ball; 4]>; N] = [const { None }; N];
        let mut workspace = self
            .data
            .degrees
            .map(|degree| Vec::with_capacity(degree + 1));
        for (output, net) in values.iter_mut().zip(&self.data.jets) {
            *output = Some(match net {
                Some(net) => net.enclose([u, v], precision, dec, &mut workspace)?,
                None => std::array::from_fn(|_| zero(precision)),
            });
        }
        let homogeneous =
            values.map(|value| value.expect("every requested homogeneous derivative assigned"));
        // Partition of unity bounds the denominator independently of interval
        // dependencies in de Boor's expression.
        let denominator = if self.data.constant_weight {
            dec.enclose_at(&self.data.weights[0], precision)?
        } else {
            bounds(
                &homogeneous[0][3].lo().max(self.data.weight_lower.clone()),
                &homogeneous[0][3].hi(),
                precision,
            )
        };
        Ok((homogeneous, denominator))
    }
}

fn cartesian(value: &[Ball; 4], denominator: &Ball) -> Result<[Ball; 3], Refusal> {
    let divide = |i: usize| {
        value[i]
            .checked_div(denominator)
            .ok_or_else(|| invalid("positive chart denominator was not separated"))
    };
    Ok([divide(0)?, divide(1)?, divide(2)?])
}

fn quotient_first(
    value: &[Ball; 3],
    derivative: &[Ball; 4],
    denominator: &Ball,
) -> Result<[Ball; 3], Refusal> {
    let component = |i: usize| {
        let correction = &value[i] * &derivative[3];
        (&derivative[i] - &correction)
            .checked_div(denominator)
            .ok_or_else(|| invalid("positive chart denominator was not separated"))
    };
    Ok([component(0)?, component(1)?, component(2)?])
}
