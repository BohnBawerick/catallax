//! One-dimensional carriers, including exact views of clamped surface boundaries.

use super::{
    BSplineSurface, Net, bounds, de_boor, invalid, validate_axis, validate_controls, zero,
};
use crate::{Point3, Pose, math};
use catallax_numbers::{Ball, Decider, Rational, Real, RealPoly, Refusal, Sign};
use std::sync::Arc;

/// One clamped boundary of a tensor-product surface. A view shares its controls.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplineBoundary {
    /// First u row; the curve parameter is v.
    UMin,
    /// Last u row; the curve parameter is v.
    UMax,
    /// First v column; the curve parameter is u.
    VMin,
    /// Last v column; the curve parameter is u.
    VMax,
}
impl SplineBoundary {
    fn axis(self) -> usize {
        match self {
            Self::UMin | Self::UMax => 1,
            _ => 0,
        }
    }
    fn fixed(self, counts: [usize; 2]) -> usize {
        match self {
            Self::UMin | Self::VMin => 0,
            Self::UMax => counts[0] - 1,
            Self::VMax => counts[1] - 1,
        }
    }
}

#[derive(Debug)]
struct Owned {
    degree: usize,
    knots: Arc<[Rational]>,
    controls: Vec<Point3>,
    weights: Arc<[Real]>,
    lower: Rational,
    constant: bool,
    nets: [Option<Net>; 3],
    transport: Option<CurveTransport>,
}
#[derive(Debug)]
struct CurveTransport {
    source: BSplineCurve,
    pose: Pose,
    inverse: Pose,
}
#[derive(Clone, Debug)]
enum Definition {
    Boundary {
        surface: BSplineSurface,
        side: SplineBoundary,
    },
    Owned(Arc<Owned>),
}

/// Finite positive-weight B-spline curve. A surface-boundary view does not fit,
/// copy or independently move an edge. Regularity and closure are separate proofs.
#[derive(Clone, Debug)]
pub struct BSplineCurve {
    definition: Definition,
}

/// Position and first/second parameter derivatives, including knot-side hulls.
#[derive(Clone, Debug)]
pub struct BSplineCurveJet {
    /// Curve position in millimetres.
    pub value: [Ball; 3],
    /// First parameter derivative, not a unit tangent.
    pub first: [Ball; 3],
    /// Second parameter derivative, not curvature.
    pub second: [Ball; 3],
}

impl BSplineCurve {
    /// Validates a standalone finite clamped curve with exact controls/weights.
    pub fn new(
        degree: usize,
        knots: Vec<Rational>,
        controls: Vec<Point3>,
        weights: Vec<Real>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if controls.len() != weights.len() {
            return Err(invalid("curve control and weight counts differ"));
        }
        validate_axis(degree, &knots, controls.len(), dec)?;
        let (lower, constant) = validate_controls(&controls, &weights, dec)?;
        Ok(Self::from_validated(
            degree,
            knots.into(),
            controls,
            weights.into(),
            lower,
            constant,
            None,
        ))
    }
    fn from_validated(
        degree: usize,
        knots: Arc<[Rational]>,
        controls: Vec<Point3>,
        weights: Arc<[Real]>,
        lower: Rational,
        constant: bool,
        transport: Option<CurveTransport>,
    ) -> Self {
        let net = Net {
            degrees: [degree, 0],
            knots: [
                knots.clone(),
                Arc::from([Rational::zero(), Rational::one()]),
            ],
            count: [controls.len(), 1],
            controls: controls
                .iter()
                .zip(weights.iter())
                .map(|(p, w)| [&p.x * w, &p.y * w, &p.z * w, w.clone()])
                .collect(),
        };
        let first = net.derivative(0);
        let second = first.as_ref().and_then(|n| n.derivative(0));
        Self {
            definition: Definition::Owned(Arc::new(Owned {
                degree,
                knots,
                controls,
                weights,
                lower,
                constant,
                nets: [Some(net), first, second],
                transport,
            })),
        }
    }
    /// Exact boundary view of an existing clamped surface. No controls are copied.
    pub fn boundary(surface: &BSplineSurface, side: SplineBoundary) -> Self {
        Self {
            definition: Definition::Boundary {
                surface: surface.clone(),
                side,
            },
        }
    }
    /// Whether this is the named shared data view; false is not geometric inequality.
    pub fn is_boundary_of(&self, surface: &BSplineSurface, side: SplineBoundary) -> bool {
        matches!(&self.definition, Definition::Boundary { surface: original, side: boundary }
            if *boundary == side && Arc::ptr_eq(&original.data, &surface.data))
    }
    fn transport(&self) -> Option<&CurveTransport> {
        match &self.definition {
            Definition::Owned(data) => data.transport.as_ref(),
            Definition::Boundary { .. } => None,
        }
    }
    fn axis(&self) -> usize {
        match &self.definition {
            Definition::Boundary { side, .. } => side.axis(),
            _ => 0,
        }
    }
    fn fixed(&self, net: &Net) -> usize {
        match &self.definition {
            Definition::Boundary { side, .. } => side.fixed(net.count),
            _ => 0,
        }
    }
    fn net(&self, order: usize) -> Option<&Net> {
        match &self.definition {
            Definition::Owned(data) => data.nets[order].as_ref(),
            Definition::Boundary { surface, side } => {
                let index = if side.axis() == 0 {
                    [0, 1, 3]
                } else {
                    [0, 2, 5]
                };
                surface.data.jets[index[order]].as_ref()
            }
        }
    }
    /// Polynomial degree, independent of control count.
    pub fn degree(&self) -> usize {
        match &self.definition {
            Definition::Owned(d) => d.degree,
            Definition::Boundary { surface, side } => surface.data.degrees[side.axis()],
        }
    }
    /// Number of controls along the curve.
    pub fn control_count(&self) -> usize {
        match &self.definition {
            Definition::Owned(d) => d.controls.len(),
            Definition::Boundary { surface, side } => surface.data.count[side.axis()],
        }
    }
    /// Exact expanded knot vector.
    pub fn knots(&self) -> &[Rational] {
        match &self.definition {
            Definition::Owned(d) => &d.knots,
            Definition::Boundary { surface, side } => &surface.data.knots[side.axis()],
        }
    }
    fn index(i: usize, side: SplineBoundary, counts: [usize; 2]) -> usize {
        assert!(i < counts[side.axis()]);
        if side.axis() == 0 {
            i * counts[1] + side.fixed(counts)
        } else {
            side.fixed(counts) * counts[1] + i
        }
    }
    /// One exact control. Out-of-range indices are programming errors.
    pub fn control(&self, i: usize) -> &Point3 {
        match &self.definition {
            Definition::Owned(d) => &d.controls[i],
            Definition::Boundary { surface, side } => {
                &surface.data.controls[Self::index(i, *side, surface.data.count)]
            }
        }
    }
    /// One exact positive weight.
    pub fn weight(&self, i: usize) -> &Real {
        match &self.definition {
            Definition::Owned(d) => &d.weights[i],
            Definition::Boundary { surface, side } => {
                &surface.data.weights[Self::index(i, *side, surface.data.count)]
            }
        }
    }
    /// Shared exact controls in parameter order, without a temporary array.
    pub fn controls(&self) -> impl ExactSizeIterator<Item = &Point3> {
        (0..self.control_count()).map(move |i| self.control(i))
    }
    /// Exact positive weights in parameter order.
    pub fn weights(&self) -> impl ExactSizeIterator<Item = &Real> {
        (0..self.control_count()).map(move |i| self.weight(i))
    }
    /// Closed active finite parameter interval.
    pub fn domain(&self) -> [&Rational; 2] {
        [
            &self.knots()[self.degree()],
            &self.knots()[self.control_count()],
        ]
    }
    /// Ordered nonzero spans, with their expanded-knot indices.
    pub fn spans(&self) -> impl Iterator<Item = (usize, [&Rational; 2])> {
        (self.degree()..self.control_count())
            .filter(move |&i| self.knots()[i] < self.knots()[i + 1])
            .map(move |i| (i, [&self.knots()[i], &self.knots()[i + 1]]))
    }
    fn weight_bound(&self) -> (&Rational, Option<&Real>) {
        match &self.definition {
            Definition::Owned(d) => (&d.lower, d.constant.then(|| &d.weights[0])),
            Definition::Boundary { surface, .. } => (
                &surface.data.weight_lower,
                surface
                    .data
                    .constant_weight
                    .then(|| &surface.data.weights[0]),
            ),
        }
    }
    fn denominator(&self, raw: &Ball, precision: u32, dec: &Decider) -> Result<Ball, Refusal> {
        let (lower, constant) = self.weight_bound();
        match constant {
            Some(w) => dec.enclose_at(w, precision),
            None => Ok(bounds(&raw.lo().max(lower.clone()), &raw.hi(), precision)),
        }
    }
    fn homogeneous<const N: usize>(
        &self,
        parameter: [&Rational; 2],
        precision: u32,
        dec: &Decider,
        row: &mut Vec<[Ball; 4]>,
    ) -> Result<[[Ball; 4]; N], Refusal> {
        dec.check_deadline()?;
        if precision == 0 || precision > dec.policy().max_precision_bits {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "curve chart precision exceeds caller policy",
            ));
        }
        let domain = self.domain();
        if parameter[0] > parameter[1] || parameter[0] < domain[0] || parameter[1] > domain[1] {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "curve box leaves the finite knot domain",
            ));
        }
        let mut values: [Option<[Ball; 4]>; N] = [const { None }; N];
        for (order, value) in values.iter_mut().enumerate() {
            *value = Some(match self.net(order) {
                Some(net) => enclose_row(
                    net,
                    self.axis(),
                    self.fixed(net),
                    parameter,
                    precision,
                    dec,
                    row,
                )?,
                None => std::array::from_fn(|_| zero(precision)),
            });
        }
        Ok(values.map(|value| value.expect("every requested derivative was assigned")))
    }
    /// Certified position enclosure at the supplied binary-box precision.
    pub fn evaluate(&self, parameter: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        let interval = [parameter.lo(), parameter.hi()];
        self.evaluate_box(interval.each_ref(), parameter.precision(), dec)
    }
    /// Certified position over an exact rational interval, including finite endpoints.
    pub fn evaluate_box(
        &self,
        parameter: [&Rational; 2],
        precision: u32,
        dec: &Decider,
    ) -> Result<[Ball; 3], Refusal> {
        let mut row = Vec::with_capacity(self.degree() + 1);
        let [value] = self.homogeneous(parameter, precision, dec, &mut row)?;
        cartesian(&value, &self.denominator(&value[3], precision, dec)?)
    }
    /// Position and first/second derivatives over a finite exact parameter interval.
    pub fn jet_box(
        &self,
        parameter: [&Rational; 2],
        precision: u32,
        dec: &Decider,
    ) -> Result<BSplineCurveJet, Refusal> {
        let mut row = Vec::with_capacity(self.degree() + 1);
        let [p, d, dd] = self.homogeneous(parameter, precision, dec, &mut row)?;
        let denominator = self.denominator(&p[3], precision, dec)?;
        let value = cartesian(&p, &denominator)?;
        let derivative = |i: usize| {
            let correction = &value[i] * &d[3];
            (&d[i] - &correction)
                .checked_div(&denominator)
                .ok_or_else(|| invalid("unresolved positive curve denominator"))
        };
        let first = [derivative(0)?, derivative(1)?, derivative(2)?];
        let twice = Ball::from_rational(&Rational::from_i64(2), precision);
        let derivative = |i: usize| {
            let first_term = &(&first[i] * &d[3]) * &twice;
            let second_term = &value[i] * &dd[3];
            let correction = &first_term + &second_term;
            (&dd[i] - &correction)
                .checked_div(&denominator)
                .ok_or_else(|| invalid("unresolved positive curve denominator"))
        };
        let second = [derivative(0)?, derivative(1)?, derivative(2)?];
        Ok(BSplineCurveJet {
            value,
            first,
            second,
        })
    }
    /// Reuses the span workspace and evaluates no unrequested second derivative.
    pub(crate) fn first_box_reusing(
        &self,
        parameter: [&Rational; 2],
        precision: u32,
        dec: &Decider,
        row: &mut Vec<[Ball; 4]>,
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        let [p, d] = self.homogeneous(parameter, precision, dec, row)?;
        let denominator = self.denominator(&p[3], precision, dec)?;
        let value = cartesian(&p, &denominator)?;
        let first = if self.weight_bound().1.is_some() {
            cartesian(&d, &denominator)?
        } else {
            let derivative = |i: usize| {
                let correction = &value[i] * &d[3];
                (&d[i] - &correction)
                    .checked_div(&denominator)
                    .ok_or_else(|| invalid("unresolved positive curve denominator"))
            };
            [derivative(0)?, derivative(1)?, derivative(2)?]
        };
        Ok([value, first])
    }
    /// Exact rigid transport, retaining the source frame for inverse queries.
    /// Repeated motions compose against the original controls, not the last fit.
    pub fn transformed(&self, pose: &Pose) -> Self {
        let (source, pose) = match self.transport() {
            None => (self, pose.clone()),
            Some(transport) => (&transport.source, transport.pose.then(pose)),
        };
        let (lower, constant) = source.weight_bound();
        let knots = match &source.definition {
            Definition::Owned(data) => data.knots.clone(),
            Definition::Boundary { surface, side } => surface.data.knots[side.axis()].clone(),
        };
        let weights = match &source.definition {
            Definition::Owned(data) => data.weights.clone(),
            Definition::Boundary { .. } => source.weights().cloned().collect(),
        };
        let inverse = pose.inverse();
        Self::from_validated(
            source.degree(),
            knots,
            source.controls().map(|p| pose.apply_point(p)).collect(),
            weights,
            lower.clone(),
            constant.is_some(),
            Some(CurveTransport {
                source: source.clone(),
                pose,
                inverse,
            }),
        )
    }
    /// Safe positive definition identity; false does not prove different point sets.
    pub fn same_definition(&self, other: &Self) -> bool {
        self.degree() == other.degree()
            && self.knots() == other.knots()
            && self.control_count() == other.control_count()
            && self.controls().zip(other.controls()).all(|(a, b)| {
                [&a.x, &a.y, &a.z]
                    .into_iter()
                    .zip([&b.x, &b.y, &b.z])
                    .all(|(a, b)| a.same_definition(b))
            })
            && self
                .weights()
                .zip(other.weights())
                .all(|(a, b)| a.same_definition(b))
    }
    /// Whether all defining fields have cheap exact rational values.
    pub fn has_rational_definition(&self) -> bool {
        self.controls().all(|p| p.rational_components().is_some())
            && self.weights().all(|w| w.rational_ref().is_some())
    }
    /// Positive construction proof of closure; otherwise a geometric query is needed.
    pub fn is_closed(&self) -> Option<bool> {
        let a = self.control(0);
        let b = self.control(self.control_count() - 1);
        [&a.x, &a.y, &a.z]
            .into_iter()
            .zip([&b.x, &b.y, &b.z])
            .all(|(a, b)| a.same_definition(b))
            .then_some(true)
    }
}

fn cartesian(value: &[Ball; 4], denominator: &Ball) -> Result<[Ball; 3], Refusal> {
    let divide = |i: usize| {
        value[i]
            .checked_div(denominator)
            .ok_or_else(|| invalid("unresolved positive curve denominator"))
    };
    Ok([divide(0)?, divide(1)?, divide(2)?])
}
fn enclose_row(
    net: &Net,
    axis: usize,
    fixed: usize,
    parameter: [&Rational; 2],
    precision: u32,
    dec: &Decider,
    row: &mut Vec<[Ball; 4]>,
) -> Result<[Ball; 4], Refusal> {
    let mut result: Option<[Ball; 4]> = None;
    for span in net.spans(axis) {
        let lo = parameter[0].clone().max(net.knots[axis][span].clone());
        let hi = parameter[1].clone().min(net.knots[axis][span + 1].clone());
        if lo > hi {
            continue;
        }
        row.clear();
        for i in span - net.degrees[axis]..=span {
            dec.check_deadline()?;
            let index = if axis == 0 {
                i * net.count[1] + fixed
            } else {
                fixed * net.count[1] + i
            };
            let p = &net.controls[index];
            row.push([
                dec.enclose_at(&p[0], precision)?,
                dec.enclose_at(&p[1], precision)?,
                dec.enclose_at(&p[2], precision)?,
                dec.enclose_at(&p[3], precision)?,
            ]);
        }
        let next = de_boor(
            row,
            &net.knots[axis],
            span,
            &bounds(&lo, &hi, precision),
            precision,
            dec,
        )?;
        match &mut result {
            None => result = Some(next),
            Some(old) => {
                for (a, b) in old.iter_mut().zip(next) {
                    *a = a.hull(&b);
                }
            }
        }
    }
    result.ok_or_else(|| invalid("curve box meets no nonzero knot span"))
}

/// One exact preimage of a point, bound to its immutable curve and knot span.
/// The value is recovered by a common polynomial factor, not a rounded inverse.
#[derive(Debug)]
pub struct SplineCurveParameter<'a> {
    curve: &'a BSplineCurve,
    span: usize,
    value: Real,
}
impl SplineCurveParameter<'_> {
    /// Exact global knot parameter.
    pub fn value(&self) -> &Real {
        &self.value
    }
    /// Expanded-knot index selecting the certified branch.
    pub fn span(&self) -> usize {
        self.span
    }
    /// Exact point on the selected branch.
    pub fn point(&self, dec: &Decider) -> Result<Point3, Refusal> {
        self.curve.point_on_span(self.span, &self.value, dec)
    }
    /// Exact first derivative of the selected branch.
    pub fn first(&self, dec: &Decider) -> Result<crate::Vector3, Refusal> {
        self.curve.first_on_span(self.span, &self.value, dec)
    }
    /// Arriving and leaving derivatives at this parameter, including a closed
    /// seam. A stationary branch refuses instead of inventing a tangent.
    pub fn tangents(&self, dec: &mut Decider) -> Result<[crate::Vector3; 2], Refusal> {
        let mut arriving = None;
        let mut leaving = None;
        for (span, interval) in self.curve.spans() {
            let from = math::sign(
                &(&self.value - &Real::from_rational(interval[0].clone())),
                dec,
            )?;
            let to = math::sign(
                &(&self.value - &Real::from_rational(interval[1].clone())),
                dec,
            )?;
            if from == Sign::Positive && to == Sign::Negative {
                let tangent = self.curve.first_on_span(span, &self.value, dec)?;
                if math::sign(&tangent.norm_squared(), dec)? == Sign::Zero {
                    return Err(invalid("stationary spline branch has no regular tangent"));
                }
                return Ok([tangent.clone(), tangent]);
            }
            if to == Sign::Zero {
                arriving = Some(self.curve.first_on_span(span, &self.value, dec)?);
            }
            if from == Sign::Zero {
                leaving = Some(self.curve.first_on_span(span, &self.value, dec)?);
            }
        }
        let domain = self.curve.domain();
        if self.curve.is_closed() == Some(true) {
            if math::sign(
                &(&self.value - &Real::from_rational(domain[0].clone())),
                dec,
            )? == Sign::Zero
            {
                arriving = Some(self.curve.first_on_span(
                    self.curve.control_count() - 1,
                    &Real::from_rational(domain[1].clone()),
                    dec,
                )?);
            }
            if math::sign(
                &(&self.value - &Real::from_rational(domain[1].clone())),
                dec,
            )? == Sign::Zero
            {
                leaving = Some(self.curve.first_on_span(
                    self.curve.degree(),
                    &Real::from_rational(domain[0].clone()),
                    dec,
                )?);
            }
        }
        let tangents = match (arriving, leaving) {
            (Some(a), Some(b)) => [a, b],
            (Some(tangent), None) | (None, Some(tangent)) => [tangent.clone(), tangent],
            (None, None) => return Err(invalid("spline preimage leaves its finite domain")),
        };
        for tangent in &tangents {
            if math::sign(&tangent.norm_squared(), dec)? == Sign::Zero {
                return Err(invalid("stationary spline branch has no regular tangent"));
            }
        }
        Ok(tangents)
    }
}

impl BSplineCurve {
    /// Homogeneous coordinate polynomials on one nonzero knot span.
    /// The variable is the original global parameter; no span is extrapolated.
    pub fn span_polynomials(&self, span: usize, dec: &Decider) -> Result<[RealPoly; 4], Refusal> {
        dec.check_deadline()?;
        if span < self.degree()
            || span >= self.control_count()
            || self.knots()[span] >= self.knots()[span + 1]
        {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "not a nonzero curve knot span",
            ));
        }
        let net = self.net(0).expect("curve value net");
        let axis = self.axis();
        let fixed = self.fixed(net);
        let degree = self.degree();
        let mut values = Vec::with_capacity(degree + 1);
        for i in span - degree..=span {
            let index = if axis == 0 {
                i * net.count[1] + fixed
            } else {
                fixed * net.count[1] + i
            };
            values.push(
                net.controls[index]
                    .each_ref()
                    .map(|v| RealPoly::constant(1, v.clone())),
            );
        }
        let variable = RealPoly::var(1, 0);
        let one = RealPoly::constant(1, Real::one());
        let constant = self.weight_bound().1;
        let components = if constant.is_some() { 3 } else { 4 };
        for r in 1..=degree {
            for j in (r..=degree).rev() {
                dec.check_deadline()?;
                let i = span - degree + j;
                let width = &self.knots()[i + degree + 1 - r] - &self.knots()[i];
                if width <= Rational::zero() {
                    return Err(invalid("zero-support symbolic curve basis"));
                }
                let alpha = (&variable - &RealPoly::constant(1, self.knots()[i].clone().into()))
                    .scale(&Real::from_rational(width.recip().unwrap()));
                let complement = &one - &alpha;
                let (before, current) = values.split_at_mut(j);
                for (next, previous) in current[0][..components]
                    .iter_mut()
                    .zip(&before[j - 1][..components])
                {
                    *next = &complement * previous + &alpha * &*next;
                }
            }
        }
        let mut value = values.pop().expect("nonempty curve basis");
        if let Some(weight) = constant {
            value[3] = RealPoly::constant(1, weight.clone());
        }
        Ok(value)
    }

    fn span_for(&self, parameter: &Real, dec: &mut Decider) -> Result<usize, Refusal> {
        dec.check_deadline()?;
        math::validate(parameter, dec)?;
        for (span, interval) in self.spans() {
            dec.check_deadline()?;
            let contained = if let Some(value) = parameter.rational_ref() {
                value >= interval[0] && value <= interval[1]
            } else {
                math::sign(
                    &(parameter - &Real::from_rational(interval[0].clone())),
                    dec,
                )? != Sign::Negative
                    && math::sign(
                        &(parameter - &Real::from_rational(interval[1].clone())),
                        dec,
                    )? != Sign::Positive
            };
            if contained {
                return Ok(span);
            }
        }
        Err(Refusal::new(
            "shapes.bspline-domain",
            "exact parameter leaves the finite curve",
        ))
    }
    fn point_on_span(
        &self,
        span: usize,
        parameter: &Real,
        dec: &Decider,
    ) -> Result<Point3, Refusal> {
        if let Some(transport) = self.transport() {
            return Ok(transport
                .pose
                .apply_point(&transport.source.point_on_span(span, parameter, dec)?));
        }
        let p = self.span_polynomials(span, dec)?;
        let denominator = p[3].eval_real(std::slice::from_ref(parameter));
        Ok(Point3::new(
            &p[0].eval_real(std::slice::from_ref(parameter)) / &denominator,
            &p[1].eval_real(std::slice::from_ref(parameter)) / &denominator,
            &p[2].eval_real(std::slice::from_ref(parameter)) / &denominator,
        ))
    }
    fn first_on_span(
        &self,
        span: usize,
        parameter: &Real,
        dec: &Decider,
    ) -> Result<crate::Vector3, Refusal> {
        if let Some(transport) = self.transport() {
            return Ok(transport
                .pose
                .apply_vector(&transport.source.first_on_span(span, parameter, dec)?));
        }
        let p = self.span_polynomials(span, dec)?;
        let denominator = p[3].eval_real(std::slice::from_ref(parameter));
        let weight_derivative = p[3].partial(0).eval_real(std::slice::from_ref(parameter));
        let derivative = |i: usize| {
            let value = p[i].eval_real(std::slice::from_ref(parameter));
            (&p[i].partial(0).eval_real(std::slice::from_ref(parameter)) * &denominator
                - &value * &weight_derivative)
                / denominator.square()
        };
        Ok(crate::Vector3::new(
            derivative(0),
            derivative(1),
            derivative(2),
        ))
    }
    /// Exact position at a proven finite parameter, without fitting or rounding.
    pub fn point_at(&self, parameter: &Real, dec: &mut Decider) -> Result<Point3, Refusal> {
        let span = self.span_for(parameter, dec)?;
        self.point_on_span(span, parameter, dec)
    }
    /// Exact first derivative on the selected knot side.
    pub fn first_at(&self, parameter: &Real, dec: &mut Decider) -> Result<crate::Vector3, Refusal> {
        let span = self.span_for(parameter, dec)?;
        self.first_on_span(span, parameter, dec)
    }

    /// Finds an exact point preimage in deterministic knot-span order.
    /// Coordinate-polynomial gcd recovers a common parameter directly in the
    /// existing coefficient field, avoiding unrelated one-coordinate roots.
    /// A closed seam selects its first preimage; this is not a uniqueness proof
    /// for an arbitrary unchecked, possibly self-intersecting curve.
    pub fn parameter_of<'a>(
        &'a self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<SplineCurveParameter<'a>>, Refusal> {
        self.preimage(point, dec, false)
    }
    /// Finds a unique preimage, identifying structurally shared closed endpoints.
    /// Other multiple preimages refuse; ordering never selects an arbitrary branch.
    pub fn unique_parameter_of<'a>(
        &'a self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<SplineCurveParameter<'a>>, Refusal> {
        self.preimage(point, dec, true)
    }
    fn preimage<'a>(
        &'a self,
        point: &Point3,
        dec: &mut Decider,
        unique: bool,
    ) -> Result<Option<SplineCurveParameter<'a>>, Refusal> {
        dec.check_deadline()?;
        if let Some(transport) = self.transport() {
            return Ok(transport
                .source
                .preimage(&transport.inverse.apply_point(point), dec, unique)?
                .map(|parameter| SplineCurveParameter {
                    curve: self,
                    span: parameter.span,
                    value: parameter.value,
                }));
        }
        math::validate_point(point, dec)?;
        let bits = 128.min(dec.policy().max_precision_bits);
        if bits == 0 {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "no curve-inverse working precision",
            ));
        }
        let query = [
            dec.enclose_at(&point.x, bits)?,
            dec.enclose_at(&point.y, bits)?,
            dec.enclose_at(&point.z, bits)?,
        ];
        let mut found: Option<SplineCurveParameter<'a>> = None;
        let canonical_seam = unique && self.is_closed() == Some(true);
        for (span, interval) in self.spans() {
            dec.check_deadline()?;
            let mut possible = true;
            for axis in 0..3 {
                let mut hull: Option<Ball> = None;
                for i in span - self.degree()..=span {
                    let control = self.control(i);
                    let coordinate = [&control.x, &control.y, &control.z][axis];
                    let value = dec.enclose_at(coordinate, bits)?;
                    hull = Some(match hull {
                        None => value,
                        Some(old) => old.hull(&value),
                    });
                }
                let hull = hull.expect("nonempty positive-weight span hull");
                if query[axis].hi() < hull.lo() || query[axis].lo() > hull.hi() {
                    possible = false;
                    break;
                }
            }
            if !possible {
                continue;
            }
            let coordinates = self.span_polynomials(span, dec)?;
            let mut common: Option<RealPoly> = None;
            let mut excluded = false;
            for (coordinate, target) in coordinates[..3].iter().zip([&point.x, &point.y, &point.z])
            {
                let equation = coordinate - &coordinates[3].scale(target);
                if equation.is_zero(dec)? {
                    continue;
                }
                let factor = match &common {
                    None => equation,
                    Some(old) => old.gcd(&equation, dec)?,
                };
                if factor.total_degree(dec)? == Some(0) {
                    excluded = true;
                    break;
                }
                common = Some(factor);
            }
            if excluded {
                continue;
            }
            let Some(common) = common else {
                return Err(Refusal::new(
                    "shapes.spline-inverse",
                    "constant curve has no isolated point preimage",
                ));
            };
            for parameter in small_roots(&common, dec)? {
                if math::sign(
                    &(&parameter - &Real::from_rational(interval[0].clone())),
                    dec,
                )? != Sign::Negative
                    && math::sign(
                        &(&parameter - &Real::from_rational(interval[1].clone())),
                        dec,
                    )? != Sign::Positive
                {
                    let mut candidate = SplineCurveParameter {
                        curve: self,
                        span,
                        value: parameter,
                    };
                    if !unique {
                        return Ok(Some(candidate));
                    }
                    if canonical_seam
                        && math::sign(
                            &(&candidate.value - &Real::from_rational(self.domain()[1].clone())),
                            dec,
                        )? == Sign::Zero
                    {
                        candidate.span = self.degree();
                        candidate.value = Real::from_rational(self.domain()[0].clone());
                    }
                    if let Some(old) = &found {
                        if math::sign(&(&candidate.value - &old.value), dec)? != Sign::Zero {
                            return Err(Refusal::new(
                                "shapes.spline-ambiguous-parameter",
                                "point has distinct finite spline preimages",
                            ));
                        }
                    } else {
                        found = Some(candidate);
                    }
                }
            }
        }
        Ok(found)
    }
}

fn small_roots(polynomial: &RealPoly, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
    let degree = polynomial.total_degree(dec)?;
    let terms = polynomial.terms();
    let coefficient = |power| {
        terms
            .iter()
            .find(|(m, _)| m.exponents()[0] == power)
            .map_or_else(Real::zero, |(_, c)| c.clone())
    };
    match degree {
        None | Some(0) => Ok(Vec::new()),
        Some(1) => Ok(vec![-coefficient(0) / coefficient(1)]),
        Some(2) => {
            let a = coefficient(2);
            let b = coefficient(1);
            let c = coefficient(0);
            let discriminant = b.square() - &Real::from_i64(4) * &(&a * &c);
            let denominator = &Real::from_i64(2) * &a;
            match math::sign(&discriminant, dec)? {
                Sign::Negative => Ok(Vec::new()),
                Sign::Zero => Ok(vec![-b / denominator]),
                Sign::Positive => {
                    let root = discriminant.sqrt();
                    let negative = -b;
                    let mut values = vec![
                        (&negative - &root) / denominator.clone(),
                        (&negative + &root) / denominator,
                    ];
                    if math::sign(&a, dec)? == Sign::Negative {
                        values.swap(0, 1);
                    }
                    Ok(values)
                }
            }
        }
        _ => polynomial.real_roots(dec),
    }
}
