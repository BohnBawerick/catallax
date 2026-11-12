//! Complete finite analytic sections with exact rational chart reductions.

use super::*;
use crate::Plane;
use catallax_numbers::Integer;

mod branch;
mod surface;
pub(crate) use branch::{Cover, search_cover};
pub use branch::{
    PatchSectionBranch, SectionBranchApproximation, SectionCandidate, SectionHits, SectionSearch,
    SectionTracePiece,
};
pub use surface::SectionSurface;

/// A complete checked finite section of one retained rational patch by an
/// admitted analytic surface. Constructors are private so raw curves cannot
/// claim this correspondence.
#[derive(Clone, Debug)]
pub struct RationalPatchSection {
    patch: RationalPatch,
    surface: SectionSurface,
    curves: Vec<(BSplineCurve, bool)>,
    points: Vec<Point3>,
    branches: Vec<PatchSectionBranch>,
    search: Option<SectionSearch>,
}

impl RationalPatchSection {
    /// The unchanged generating patch.
    pub fn patch(&self) -> &RationalPatch {
        &self.patch
    }
    /// The exact cutting surface.
    pub fn surface(&self) -> &SectionSurface {
        &self.surface
    }
    /// Exact rational spline components and whether they cross the surface.
    /// Algebraic components are retained separately in `branches`. `false`
    /// is a proved even-order tangency.
    pub fn curves(&self) -> impl Iterator<Item = (&BSplineCurve, bool)> {
        self.curves
            .iter()
            .map(|(curve, crossing)| (curve, *crossing))
    }
    /// Isolated contacts not covered by the returned curve pieces.
    pub fn points(&self) -> &[Point3] {
        &self.points
    }
    /// Exact algebraic components with a complete parameter-cell cover.
    pub fn branches(&self) -> &[PatchSectionBranch] {
        &self.branches
    }
    /// Work and cell counts for a bivariate section search.
    pub fn search(&self) -> Option<&SectionSearch> {
        self.search.as_ref()
    }
}

fn crossing(
    polynomial: &RealPoly,
    root: &Real,
    span: &[Rational; 2],
    domain: &[Rational; 2],
    regular_contact_code: Option<&'static str>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let degree = polynomial
        .total_degree(dec)?
        .ok_or_else(|| failure("a contact equation is identically zero"))?;
    let mut derivative = polynomial.clone();
    for order in 1..=degree {
        derivative = derivative.partial(0);
        if exact_sign(&derivative.eval_real(std::slice::from_ref(root)), dec)? == Sign::Zero {
            continue;
        }
        if order > 1 {
            if let Some(code) = regular_contact_code {
                return Err(Refusal::new(
                    code,
                    "a multiple section root needs a separate contact certificate",
                ));
            }
            let at = |bound: &Rational, dec: &mut Decider| {
                exact_sign(&(root - &Real::from_rational(bound.clone())), dec)
                    .map(|sign| sign == Sign::Zero)
            };
            if (at(&span[0], dec)? || at(&span[1], dec)?)
                && !(at(&domain[0], dec)? || at(&domain[1], dec)?)
            {
                return Err(failure(
                    "a higher-order knot contact needs a two-sided branch certificate",
                ));
            }
        }
        return Ok(order % 2 == 1);
    }
    Err(failure("the root has no certified finite contact order"))
}

fn failure(detail: &str) -> Refusal {
    Refusal::new("intersect.rational-patch-section", detail)
}

fn axis_polynomial(poly: &RealPoly, axis: usize, dec: &mut Decider) -> Result<RealPoly, Refusal> {
    let mut terms = Vec::new();
    for (monomial, coefficient) in poly.terms() {
        if monomial.exponents()[1 - axis] != 0 {
            if exact_sign(&coefficient, dec)? == Sign::Zero {
                continue;
            }
            return Err(failure(
                "section reduction retained a second free parameter",
            ));
        }
        terms.push((vec![monomial.exponents()[axis]], coefficient));
    }
    Ok(RealPoly::from_terms(1, terms))
}

fn degrees(poly: &RealPoly, dec: &mut Decider) -> Result<[u32; 2], Refusal> {
    let mut result = [0, 0];
    for (monomial, coefficient) in poly.terms() {
        if exact_sign(&coefficient, dec)? != Sign::Zero {
            for (degree, &exponent) in result.iter_mut().zip(monomial.exponents()) {
                *degree = (*degree).max(exponent);
            }
        }
    }
    Ok(result)
}

fn in_range(value: &Real, domain: &[Rational; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        exact_sign(&(value - &Real::from_rational(domain[0].clone())), dec)? != Sign::Negative
            && exact_sign(&(&Real::from_rational(domain[1].clone()) - value), dec)?
                != Sign::Negative,
    )
}

fn parameter_curve(
    coordinates: [RealPoly; 4],
    interval: [Real; 2],
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let scale = &interval[1] - &interval[0];
    let variable = RealPoly::var(1, 0).scale(&scale) + RealPoly::constant(1, interval[0].clone());
    let mut coordinates = coordinates.map(|p| p.substitute(0, &variable));
    let middle = coordinates[3].eval_real(&[Real::from_rational(
        Rational::from_ratio(1, 2).expect("nonzero"),
    )]);
    match exact_sign(&middle, dec)? {
        Sign::Zero => return Err(failure("a section has no positive rational denominator")),
        Sign::Negative => coordinates = coordinates.map(|p| p.scale(&Real::from_i64(-1))),
        Sign::Positive => (),
    }
    let mut degree = 1usize;
    for polynomial in &coordinates {
        degree = degree.max(polynomial.total_degree(dec)?.unwrap_or(0) as usize);
    }
    let mut powers = vec![std::array::from_fn::<_, 4, _>(|_| Real::zero()); degree + 1];
    for (axis, polynomial) in coordinates.iter().enumerate() {
        for (monomial, coefficient) in polynomial.terms() {
            powers[monomial.exponents()[0] as usize][axis] = coefficient;
        }
    }
    // A strictly positive polynomial has positive Bernstein coefficients after
    // enough degree elevation. Failure to prove them never changes the curve.
    for elevated in degree..=degree.saturating_add(64) {
        dec.check_deadline()?;
        let mut homogeneous: Vec<[Real; 4]> = Vec::with_capacity(elevated + 1);
        for k in 0..=elevated {
            let mut value = std::array::from_fn(|_| Real::zero());
            let mut factor = Rational::one();
            for (power, coefficients) in powers.iter().enumerate().take(k.min(degree) + 1) {
                if power > 0 {
                    factor = factor
                        * Rational::from_integer(Integer::from_u64((k - power + 1) as u64))
                        / Rational::from_integer(Integer::from_u64((elevated - power + 1) as u64));
                }
                let factor = Real::from_rational(factor.clone());
                for (value, coefficient) in value.iter_mut().zip(coefficients) {
                    *value = &*value + &(coefficient * &factor);
                }
            }
            homogeneous.push(value);
        }
        let mut positive = true;
        for control in &homogeneous {
            if exact_sign(&control[3], dec)? != Sign::Positive {
                positive = false;
                break;
            }
        }
        if !positive {
            continue;
        }
        let mut controls = Vec::with_capacity(homogeneous.len());
        let mut weights = Vec::with_capacity(homogeneous.len());
        for [x, y, z, w] in homogeneous {
            controls.push(Point3::new(&x / &w, &y / &w, &z / &w));
            weights.push(w);
        }
        let knots = std::iter::repeat_n(Rational::zero(), elevated + 1)
            .chain(std::iter::repeat_n(Rational::one(), elevated + 1))
            .collect();
        return BSplineCurve::new(elevated, knots, controls, weights, dec);
    }
    Err(failure(
        "the exact section denominator has no proved positive Bernstein net within the degree allowance",
    ))
}

fn iso(span: &Span, axis: usize, value: &Real, dec: &mut Decider) -> Result<BSplineCurve, Refusal> {
    let coordinates = span
        .coordinates
        .each_ref()
        .map(|p| p.substitute(axis, &RealPoly::constant(2, value.clone())));
    let mut reduced = Vec::with_capacity(4);
    for p in coordinates {
        reduced.push(axis_polynomial(&p, 1 - axis, dec)?);
    }
    let coordinates: [RealPoly; 4] = reduced.try_into().expect("four homogeneous coordinates");
    parameter_curve(
        coordinates,
        span.domain[1 - axis].clone().map(Real::from_rational),
        dec,
    )
}

fn graph_coordinates(
    span: &Span,
    axis: usize,
    numerator: &RealPoly,
    denominator: &RealPoly,
) -> [RealPoly; 4] {
    let power = span
        .coordinates
        .iter()
        .flat_map(RealPoly::terms)
        .map(|(monomial, _)| monomial.exponents()[axis])
        .max()
        .unwrap_or(0);
    let variable = RealPoly::var(1, 0);
    span.coordinates.each_ref().map(|p| {
        p.terms()
            .into_iter()
            .fold(RealPoly::zero(1), |sum, (monomial, coefficient)| {
                let exponents = monomial.exponents();
                sum + (&(&variable.pow_u32(exponents[1 - axis])
                    * &numerator.pow_u32(exponents[axis]))
                    * &denominator.pow_u32(power - exponents[axis]))
                    .scale(&coefficient)
            })
    })
}

fn graph(
    span: &Span,
    equation: &RealPoly,
    solved: usize,
    dec: &mut Decider,
) -> Result<(Vec<BSplineCurve>, Vec<Point3>), Refusal> {
    let free = 1 - solved;
    let mut coefficients = [RealPoly::zero(1), RealPoly::zero(1)];
    for (monomial, coefficient) in equation.terms() {
        let exponents = monomial.exponents();
        if exponents[solved] > 1 && exact_sign(&coefficient, dec)? != Sign::Zero {
            return Err(failure("a section is not affine in the selected parameter"));
        }
        if exponents[solved] <= 1 {
            coefficients[exponents[solved] as usize] = &coefficients[exponents[solved] as usize]
                + &RealPoly::from_terms(1, vec![(vec![exponents[free]], coefficient)]);
        }
    }
    let [a, b] = coefficients;
    if b.is_zero(dec)? {
        return Err(failure("a graph section has no solving coefficient"));
    }
    let numerator = a.scale(&Real::from_i64(-1));
    let mut breaks = span.domain[free].clone().map(Real::from_rational).to_vec();
    for root in b.real_roots(dec)? {
        if in_range(&root, &span.domain[free], dec)? {
            if exact_sign(&a.eval_real(std::slice::from_ref(&root)), dec)? == Sign::Zero {
                return Err(failure(
                    "the section contains an unresolved whole parameter row",
                ));
            }
            insert(&mut breaks, root, dec)?;
        }
    }
    for bound in &span.domain[solved] {
        let crossing = &a + &b.scale(&Real::from_rational(bound.clone()));
        if crossing.is_zero(dec)? {
            continue;
        }
        for root in crossing.real_roots(dec)? {
            if in_range(&root, &span.domain[free], dec)? {
                insert(&mut breaks, root, dec)?;
            }
        }
    }
    let graph_value = |value: &Real| -> Real {
        &numerator.eval_real(std::slice::from_ref(value))
            / &b.eval_real(std::slice::from_ref(value))
    };
    let mut curves = Vec::new();
    let mut covered = vec![false; breaks.len()];
    for i in 0..breaks.len() - 1 {
        let middle = (&breaks[i] + &breaks[i + 1]) / Real::from_i64(2);
        if !in_range(&graph_value(&middle), &span.domain[solved], dec)? {
            continue;
        }
        for endpoint in [&breaks[i], &breaks[i + 1]] {
            if exact_sign(&b.eval_real(std::slice::from_ref(endpoint)), dec)? == Sign::Zero {
                return Err(failure(
                    "a finite section branch ends at an unresolved rational pole",
                ));
            }
            if !in_range(&graph_value(endpoint), &span.domain[solved], dec)? {
                return Err(failure(
                    "a section endpoint leaves the retained parameter rectangle",
                ));
            }
        }
        curves.push(parameter_curve(
            graph_coordinates(span, solved, &numerator, &b),
            [breaks[i].clone(), breaks[i + 1].clone()],
            dec,
        )?);
        covered[i] = true;
        covered[i + 1] = true;
    }
    let mut points = Vec::new();
    for (value, covered) in breaks.iter().zip(covered) {
        if covered || exact_sign(&b.eval_real(std::slice::from_ref(value)), dec)? == Sign::Zero {
            continue;
        }
        let other = graph_value(value);
        if in_range(&other, &span.domain[solved], dec)? {
            let mut at = [value.clone(), value.clone()];
            at[solved] = other;
            points.push(span.point(&at));
        }
    }
    Ok((curves, points))
}

fn same_curve(a: &BSplineCurve, b: &BSplineCurve, dec: &mut Decider) -> Result<bool, Refusal> {
    if a.same_definition(b) {
        return Ok(true);
    }
    let a = a.span_polynomials(a.degree(), dec)?;
    let b = b.span_polynomials(b.degree(), dec)?;
    for axis in 0..3 {
        if !(&a[axis] * &b[3] - &b[axis] * &a[3]).is_zero(dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether the cut equation vanishes along a whole side of some span. Such a
/// row is a section component lying on a cell edge, which the cell search
/// cannot classify; callers refuse it by name before searching.
pub(crate) fn vanishing_side_row(
    patch: &RationalPatch,
    surface: &SectionSurface,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for span in &patch.data.spans {
        let equation = surface.equation(span);
        for side in 0..4 {
            let bound = Real::from_rational(span.domain[side / 2][side % 2].clone());
            if equation
                .substitute(side / 2, &RealPoly::constant(2, bound))
                .is_zero(dec)?
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

impl RationalPatch {
    /// Complete finite sections by an admitted analytic surface, with rational
    /// splines for reducible equations and exact parameter-cell branches for
    /// regular bivariate cuts. Unresolved contacts and exhausted proof
    /// allowances return typed refusals.
    pub fn sections(
        &self,
        surface: &SectionSurface,
        dec: &mut Decider,
    ) -> Result<RationalPatchSection, Refusal> {
        surface.validate(dec)?;
        let mut section = RationalPatchSection {
            patch: self.clone(),
            surface: surface.clone(),
            curves: Vec::new(),
            points: Vec::new(),
            branches: Vec::new(),
            search: None,
        };
        for span in &self.data.spans {
            if matches!(surface, SectionSurface::Cone(_)) {
                let h = self.carrier().span_bezier_controls(span.indices, dec)?;
                if !surface.retains_span(&h, dec)? {
                    continue;
                }
            }
            let equation = surface.equation(span);
            if equation.is_zero(dec)? {
                return Err(
                    surface.unresolved("the cutting surface contains an entire finite patch span")
                );
            }
            let bits = dec.policy().max_precision_bits.min(128);
            if bits == 0 {
                return Err(failure("no section enclosure precision is available"));
            }
            if !equation
                .eval_ball(&box_balls(&span.domain, bits), dec)?
                .contains_zero()
            {
                continue;
            }
            let degree = degrees(&equation, dec)?;
            let mut curves = Vec::new();
            let mut points = Vec::new();
            if degree[0] == 0 || degree[1] == 0 {
                let axis = usize::from(degree[0] == 0);
                let polynomial = axis_polynomial(&equation, axis, dec)?;
                for root in polynomial.real_roots(dec)? {
                    if in_range(&root, &span.domain[axis], dec)? {
                        let crossing = crossing(
                            &polynomial,
                            &root,
                            &span.domain[axis],
                            &self.data.domain[axis],
                            surface.regular_contact_code(),
                            dec,
                        )?;
                        curves.push((iso(span, axis, &root, dec)?, crossing));
                    }
                }
            } else {
                let axis = if degree[1] == 1 {
                    1
                } else if degree[0] == 1 {
                    0
                } else {
                    return branch::search(self, surface, dec);
                };
                let (pieces, contacts) = graph(span, &equation, axis, dec)?;
                if !contacts.is_empty()
                    && let Some(code) = surface.regular_contact_code()
                {
                    return Err(Refusal::new(
                        code,
                        "an isolated boundary contact needs a separate contact certificate",
                    ));
                }
                curves = pieces.into_iter().map(|curve| (curve, true)).collect();
                points = contacts;
            }
            for (curve, crossing) in curves {
                let mut duplicate = false;
                for (old, old_crossing) in &section.curves {
                    if same_curve(old, &curve, dec)? {
                        if *old_crossing != crossing {
                            return Err(failure(
                                "one-sided section contact classifications disagree",
                            ));
                        }
                        duplicate = true;
                        break;
                    }
                }
                if !duplicate {
                    section.curves.push((curve, crossing));
                }
            }
            for point in points {
                let mut duplicate = false;
                for old in &section.points {
                    if exact_sign(&old.minus(&point).norm_squared(), dec)? == Sign::Zero {
                        duplicate = true;
                        break;
                    }
                }
                if !duplicate {
                    section.points.push(point);
                }
            }
        }
        let mut points = Vec::new();
        for point in section.points {
            let mut covered = false;
            for (curve, _) in &section.curves {
                if curve.unique_parameter_of(&point, dec)?.is_some() {
                    covered = true;
                    break;
                }
            }
            if !covered {
                points.push(point);
            }
        }
        section.points = points;
        Ok(section)
    }
}
