//! Exact homogeneous equations on finite tensor-product knot spans.

use super::*;
use catallax_numbers::RealPoly;

fn symbolic_de_boor(
    mut values: Vec<[RealPoly; 4]>,
    knots: &[Rational],
    span: usize,
    axis: usize,
    dec: &Decider,
) -> Result<[RealPoly; 4], Refusal> {
    let degree = values.len() - 1;
    let variable = RealPoly::var(2, axis);
    let one = RealPoly::constant(2, Real::one());
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            dec.check_deadline()?;
            let index = span - degree + j;
            let width = &knots[index + degree + 1 - r] - &knots[index];
            if width <= Rational::zero() {
                return Err(invalid("zero-support symbolic surface basis"));
            }
            let alpha = (&variable - &RealPoly::constant(2, knots[index].clone().into())).scale(
                &Real::from_rational(width.recip().expect("positive knot width")),
            );
            let complement = &one - &alpha;
            let (before, current) = values.split_at_mut(j);
            for (value, previous) in current[0].iter_mut().zip(&before[j - 1]) {
                *value = &complement * previous + &alpha * &*value;
            }
        }
    }
    Ok(values.pop().expect("nonempty symbolic basis"))
}

pub(super) fn blossom(
    values: &mut [[Real; 4]],
    knots: &[Rational],
    span: usize,
    bounds: [&Rational; 2],
    high_arguments: usize,
    dec: &Decider,
) -> Result<(), Refusal> {
    let degree = values.len() - 1;
    let one = Rational::one();
    for r in 1..=degree {
        let parameter = bounds[usize::from(r > degree - high_arguments)];
        for j in (r..=degree).rev() {
            dec.check_deadline()?;
            let index = span - degree + j;
            let width = &knots[index + degree + 1 - r] - &knots[index];
            let alpha = &(parameter - &knots[index]) / &width;
            if alpha == one {
                continue;
            }
            let (before, current) = values.split_at_mut(j);
            if alpha.is_zero() {
                current[0] = before[j - 1].clone();
                continue;
            }
            let beta = Real::from_rational(&one - &alpha);
            let alpha = Real::from_rational(alpha);
            for (value, previous) in current[0].iter_mut().zip(&before[j - 1]) {
                *value = &beta * previous + &alpha * &*value;
            }
        }
    }
    Ok(())
}

impl BSplineSurface {
    /// Expanded-knot indices of every nonzero tensor span, in u-major order.
    pub fn span_indices(&self) -> impl Iterator<Item = [usize; 2]> + '_ {
        let net = self.data.jets[0].as_ref().expect("surface value net");
        net.spans(0)
            .flat_map(move |u| net.spans(1).map(move |v| [u, v]))
    }

    /// Homogeneous [X,Y,Z,W] polynomials in the original global parameters.
    /// The caller retains the selected finite knot rectangle; these equations
    /// do not authorise extrapolation beyond it.
    pub fn span_polynomials(
        &self,
        span: [usize; 2],
        dec: &Decider,
    ) -> Result<[RealPoly; 4], Refusal> {
        dec.check_deadline()?;
        for (axis, &index) in span.iter().enumerate() {
            if index < self.data.degrees[axis]
                || index >= self.data.count[axis]
                || self.data.knots[axis][index] >= self.data.knots[axis][index + 1]
            {
                return Err(Refusal::new(
                    "shapes.bspline-domain",
                    "not a nonzero surface knot span",
                ));
            }
        }
        let net = self.data.jets[0].as_ref().expect("surface value net");
        let mut columns = Vec::with_capacity(self.data.degrees[1] + 1);
        for v in span[1] - self.data.degrees[1]..=span[1] {
            let controls = (span[0] - self.data.degrees[0]..=span[0])
                .map(|u| {
                    net.controls[u * self.data.count[1] + v]
                        .each_ref()
                        .map(|value| RealPoly::constant(2, value.clone()))
                })
                .collect();
            columns.push(symbolic_de_boor(
                controls,
                &self.data.knots[0],
                span[0],
                0,
                dec,
            )?);
        }
        let mut value = symbolic_de_boor(columns, &self.data.knots[1], span[1], 1, dec)?;
        if self.data.constant_weight {
            value[3] = RealPoly::constant(2, self.data.weights[0].clone());
        }
        Ok(value)
    }

    /// Exact positive homogeneous Bézier net on one complete knot rectangle.
    /// Both parameter axes are normalised to [0,1], with u-major controls.
    /// Blossom extraction uses convex rational blends, not fitted samples or
    /// global-power coefficients with large translated knot coordinates.
    pub fn span_bezier_controls(
        &self,
        span: [usize; 2],
        dec: &Decider,
    ) -> Result<Vec<[Real; 4]>, Refusal> {
        let bounds = self.span_domain(span)?;
        self.bezier_controls_on(span, bounds, dec)
    }

    /// Exact positive homogeneous controls on a nonempty subrectangle of one
    /// knot span. The returned net uses local [0,1] parameters; the source knot
    /// coordinates and defining surface are unchanged.
    pub fn span_bezier_controls_on(
        &self,
        span: [usize; 2],
        bounds: [[&Rational; 2]; 2],
        dec: &Decider,
    ) -> Result<Vec<[Real; 4]>, Refusal> {
        let original = self.span_domain(span)?;
        for (part, whole) in bounds.iter().zip(original) {
            if part[0] < whole[0] || part[1] > whole[1] || part[0] >= part[1] {
                return Err(Refusal::new(
                    "shapes.bspline-domain",
                    "a Bezier subrectangle must lie within its nonzero knot span",
                ));
            }
        }
        self.bezier_controls_on(span, bounds, dec)
    }

    /// Validates an exact positive homogeneous tensor net without multiplying
    /// its numerators by the weights a second time.
    pub fn from_homogeneous(
        degrees: [usize; 2],
        knots: [Vec<Rational>; 2],
        count: [usize; 2],
        homogeneous: Vec<[Real; 4]>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if count[0].checked_mul(count[1]) != Some(homogeneous.len()) {
            return Err(invalid("homogeneous tensor control counts differ"));
        }
        for axis in 0..2 {
            validate_axis(degrees[axis], &knots[axis], count[axis], dec)?;
        }
        let data = validate_homogeneous(&homogeneous, dec)?;
        let net = Net {
            degrees,
            knots: knots.map(Arc::from),
            count,
            controls: homogeneous,
        };
        Ok(Self::from_net(
            net,
            data.controls,
            data.weights,
            data.lower,
            data.constant,
        ))
    }

    /// Retains precisely a subrectangle of one knot span as a clamped exact
    /// surface, in the original global parameters. Positivity and constant
    /// weights are inherited through the convex blossom construction.
    pub fn bezier_patch(
        &self,
        span: [usize; 2],
        bounds: [[&Rational; 2]; 2],
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let homogeneous = self.span_bezier_controls_on(span, bounds, dec)?;
        let count = self.data.degrees.map(|degree| degree + 1);
        let knots = std::array::from_fn(|axis| {
            std::iter::repeat_n(bounds[axis][0].clone(), count[axis])
                .chain(std::iter::repeat_n(bounds[axis][1].clone(), count[axis]))
                .collect::<Vec<_>>()
                .into()
        });
        let mut controls = Vec::with_capacity(homogeneous.len());
        let mut weights = Vec::with_capacity(homogeneous.len());
        for point in &homogeneous {
            let weight = &point[3];
            controls.push(Point3::new(
                &point[0] / weight,
                &point[1] / weight,
                &point[2] / weight,
            ));
            weights.push(weight.clone());
        }
        let net = Net {
            degrees: self.data.degrees,
            knots,
            count,
            controls: homogeneous,
        };
        Ok(Self::from_net(
            net,
            controls,
            weights.into(),
            self.data.weight_lower.clone(),
            self.data.constant_weight,
        ))
    }

    fn span_domain(&self, span: [usize; 2]) -> Result<[[&Rational; 2]; 2], Refusal> {
        for (axis, &index) in span.iter().enumerate() {
            if index < self.data.degrees[axis]
                || index >= self.data.count[axis]
                || self.data.knots[axis][index] >= self.data.knots[axis][index + 1]
            {
                return Err(Refusal::new(
                    "shapes.bspline-domain",
                    "not a nonzero surface knot span",
                ));
            }
        }
        Ok(std::array::from_fn(|axis| {
            [
                &self.data.knots[axis][span[axis]],
                &self.data.knots[axis][span[axis] + 1],
            ]
        }))
    }

    fn bezier_controls_on(
        &self,
        span: [usize; 2],
        bounds: [[&Rational; 2]; 2],
        dec: &Decider,
    ) -> Result<Vec<[Real; 4]>, Refusal> {
        dec.check_deadline()?;
        let count = self.data.degrees.map(|degree| degree + 1);
        let source = self.data.jets[0].as_ref().expect("surface value net");
        let mut net = Vec::with_capacity(count[0] * count[1]);
        for u in span[0] - self.data.degrees[0]..=span[0] {
            for v in span[1] - self.data.degrees[1]..=span[1] {
                net.push(source.controls[u * self.data.count[1] + v].clone());
            }
        }
        let maximum = count[0].max(count[1]);
        let mut row = Vec::new();
        let mut work = Vec::new();
        for axis in 0..2 {
            let degree = self.data.degrees[axis];
            let knots = &self.data.knots[axis];
            if bounds[axis] == [&knots[span[axis]], &knots[span[axis] + 1]]
                && knots[span[axis] - degree] == knots[span[axis]]
                && knots[span[axis] + 1] == knots[span[axis] + degree + 1]
            {
                continue;
            }
            if row.capacity() == 0 {
                row.reserve_exact(maximum);
                work.reserve_exact(maximum);
            }
            for fixed in 0..count[1 - axis] {
                let at = |i| {
                    if axis == 0 {
                        i * count[1] + fixed
                    } else {
                        fixed * count[1] + i
                    }
                };
                row.clear();
                row.extend((0..count[axis]).map(|i| net[at(i)].clone()));
                for i in 0..count[axis] {
                    work.clone_from(&row);
                    blossom(&mut work, knots, span[axis], bounds[axis], i, dec)?;
                    net[at(i)] = work.last().expect("positive degree").clone();
                }
            }
        }
        Ok(net)
    }

    pub(crate) fn span_workspace(&self) -> [Vec<[Ball; 4]>; 2] {
        self.data
            .degrees
            .map(|degree| Vec::with_capacity(degree + 1))
    }

    fn homogeneous_span<const N: usize>(
        &self,
        span: [usize; 2],
        domain: &[[Rational; 2]; 2],
        precision: u32,
        dec: &Decider,
        workspace: &mut [Vec<[Ball; 4]>; 2],
    ) -> Result<([[Ball; 4]; N], Ball), Refusal> {
        for axis in 0..2 {
            if span[axis] < self.data.degrees[axis]
                || span[axis] >= self.data.count[axis]
                || domain[axis][0] > domain[axis][1]
                || domain[axis][0] < self.data.knots[axis][span[axis]]
                || domain[axis][1] > self.data.knots[axis][span[axis] + 1]
            {
                return Err(invalid("parameter box leaves its owned knot rectangle"));
            }
        }
        let orders = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]];
        let mut values: [Option<[Ball; 4]>; N] = [const { None }; N];
        for (i, output) in values.iter_mut().enumerate() {
            *output = Some(match &self.data.jets[i] {
                Some(net) => net.evaluate_span(
                    [span[0] - orders[i][0], span[1] - orders[i][1]],
                    domain.each_ref().map(|bounds| bounds.each_ref()),
                    precision,
                    dec,
                    workspace,
                )?,
                None => std::array::from_fn(|_| zero(precision)),
            });
        }
        let values = values.map(|value| value.expect("requested owned derivative assigned"));
        let denominator = self.denominator(&values[0][3], precision, dec)?;
        Ok((values, denominator))
    }

    pub(crate) fn value_on_span(
        &self,
        span: [usize; 2],
        domain: &[[Rational; 2]; 2],
        precision: u32,
        dec: &Decider,
        workspace: &mut [Vec<[Ball; 4]>; 2],
    ) -> Result<[Ball; 3], Refusal> {
        let ([value], denominator) =
            self.homogeneous_span(span, domain, precision, dec, workspace)?;
        cartesian(&value, &denominator)
    }

    pub(crate) fn derivatives_on_span(
        &self,
        span: [usize; 2],
        domain: &[[Rational; 2]; 2],
        precision: u32,
        dec: &Decider,
        workspace: &mut [Vec<[Ball; 4]>; 2],
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        let ([position, du, dv], denominator) =
            self.homogeneous_span(span, domain, precision, dec, workspace)?;
        if self.data.constant_weight {
            return Ok([cartesian(&du, &denominator)?, cartesian(&dv, &denominator)?]);
        }
        let value = cartesian(&position, &denominator)?;
        Ok([
            quotient_first(&value, &du, &denominator)?,
            quotient_first(&value, &dv, &denominator)?,
        ])
    }

    pub(crate) fn jet_on_span(
        &self,
        span: [usize; 2],
        domain: &[[Rational; 2]; 2],
        precision: u32,
        dec: &Decider,
        workspace: &mut [Vec<[Ball; 4]>; 2],
    ) -> Result<BSplineJet, Refusal> {
        let (values, denominator) =
            self.homogeneous_span(span, domain, precision, dec, workspace)?;
        Self::cartesian_jet(values, denominator, precision)
    }
}
