//! Whole-curve boundary identities, independent of sampled endpoint agreement.

use super::*;
use crate::Curve;

fn normalized_map(curve: &BSplineCurve, reverse: bool) -> RealPoly {
    let [lo, hi] = curve.domain();
    let width = hi - lo;
    let (offset, scale) = if reverse {
        (hi.clone(), -width)
    } else {
        (lo.clone(), width)
    };
    RealPoly::var(1, 0).scale(&Real::from_rational(scale))
        + RealPoly::constant(1, Real::from_rational(offset))
}

fn polynomial_at(
    curve: &BSplineCurve,
    at: &Rational,
    map: &RealPoly,
    dec: &mut Decider,
) -> Result<[RealPoly; 4], Refusal> {
    let parameter = map.eval_real(&[Real::from_rational(at.clone())]);
    let parameter = parameter
        .as_rational()
        .expect("rational normalized knot map");
    let span = curve
        .spans()
        .find_map(|(span, [lo, hi])| (&parameter >= lo && &parameter <= hi).then_some(span))
        .ok_or_else(|| domain_error("a normalized boundary parameter leaves its source span"))?;
    Ok(curve
        .span_polynomials(span, dec)?
        .map(|p| p.substitute(0, map)))
}

fn normalized_equal(
    a: &BSplineCurve,
    b: &BSplineCurve,
    reverse: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if !reverse && a.same_definition(b) {
        return Ok(true);
    }
    let mut cuts = vec![Rational::zero(), Rational::one()];
    for (curve, backwards) in [(a, false), (b, reverse)] {
        let [lo, hi] = curve.domain();
        let width = hi - lo;
        for knot in curve.knots() {
            let mut at = &(knot - lo) / &width;
            if backwards {
                at = Rational::one() - at;
            }
            cuts.push(at);
        }
    }
    cuts.sort();
    cuts.dedup();
    let map_a = normalized_map(a, false);
    let map_b = normalized_map(b, reverse);
    for interval in cuts.windows(2) {
        let at = (&interval[0] + &interval[1]).mul_pow2(-1);
        let pa = polynomial_at(a, &at, &map_a, dec)?;
        let pb = polynomial_at(b, &at, &map_b, dec)?;
        for axis in 0..3 {
            if !(&pa[axis] * &pb[3] - &pb[axis] * &pa[3]).is_zero(dec)? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn vector_polynomial(coordinates: &[RealPoly; 4], origin: &Point3) -> [RealPoly; 3] {
    let origin = [&origin.x, &origin.y, &origin.z];
    std::array::from_fn(|i| &coordinates[i] - &coordinates[3].scale(origin[i]))
}

fn polynomial_dot(vector: &[RealPoly; 3], direction: &Vector3) -> RealPoly {
    let direction = [&direction.x, &direction.y, &direction.z];
    (0..3).fold(RealPoly::zero(1), |sum, i| {
        sum + vector[i].scale(direction[i])
    })
}

fn conic_orientation(
    boundary: &BSplineCurve,
    center: &Point3,
    normal: &Vector3,
    ellipse: Option<(&Vector3, &Real, &Real)>,
    radius: Option<&Real>,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    if exact_sign(&normal.norm_squared(), dec)? != Sign::Positive {
        return Ok(None);
    }
    let mut orientation = None;
    for (span, interval) in boundary.spans() {
        let p = boundary.span_polynomials(span, dec)?;
        let radial = vector_polynomial(&p, center);
        if !polynomial_dot(&radial, normal).is_zero(dec)? {
            return Ok(None);
        }
        let equation = if let Some((major, major_sq, minor_sq)) = ellipse {
            if exact_sign(major_sq, dec)? != Sign::Positive
                || exact_sign(minor_sq, dec)? != Sign::Positive
                || exact_sign(&major.norm_squared(), dec)? != Sign::Positive
                || exact_sign(&major.dot(normal), dec)? != Sign::Zero
            {
                return Ok(None);
            }
            let across = normal.cross(major);
            let a = polynomial_dot(&radial, major);
            let b = polynomial_dot(&radial, &across);
            a.pow_u32(2).scale(&(&normal.norm_squared() * minor_sq)) + b.pow_u32(2).scale(major_sq)
                - p[3].pow_u32(2).scale(
                    &(&(&normal.norm_squared() * &major.norm_squared()) * &(major_sq * minor_sq)),
                )
        } else {
            let radius = radius.expect("circle radius supplied");
            if exact_sign(radius, dec)? != Sign::Positive {
                return Ok(None);
            }
            (0..3).fold(RealPoly::zero(1), |sum, i| sum + radial[i].pow_u32(2))
                - p[3].pow_u32(2).scale(radius)
        };
        if !equation.is_zero(dec)? {
            return Ok(None);
        }
        let derivative = radial.each_ref().map(|p| p.partial(0));
        let turning = [
            &radial[1] * &derivative[2] - &radial[2] * &derivative[1],
            &radial[2] * &derivative[0] - &radial[0] * &derivative[2],
            &radial[0] * &derivative[1] - &radial[1] * &derivative[0],
        ];
        let turning = polynomial_dot(&turning, normal);
        if turning.is_zero(dec)? {
            return Ok(None);
        }
        for root in turning.real_roots(dec)? {
            if exact_sign(&(&root - &Real::from_rational(interval[0].clone())), dec)?
                != Sign::Negative
                && exact_sign(&(&Real::from_rational(interval[1].clone()) - &root), dec)?
                    != Sign::Negative
            {
                return Ok(None);
            }
        }
        let middle = Real::from_rational((interval[0] + interval[1]).mul_pow2(-1));
        let sign = exact_sign(&turning.eval_real(&[middle]), dec)?;
        if sign == Sign::Zero {
            return Ok(None);
        }
        let forward = sign == Sign::Positive;
        if orientation.is_some_and(|old| old != forward) {
            return Ok(None);
        }
        orientation = Some(forward);
    }
    Ok(orientation)
}

impl RationalPatch {
    /// Proves that a finite edge carrier follows one complete boundary row.
    /// The Boolean says whether the carrier's natural direction agrees with
    /// increasing row parameter. Physical endpoints and loop use are checked
    /// separately by topology. Matching endpoints alone never proves a row.
    pub fn boundary_orientation(
        &self,
        side: SplineBoundary,
        curve: &Curve,
        dec: &mut Decider,
    ) -> Result<Option<bool>, Refusal> {
        self.carrier().boundary_orientation(side, curve, dec)
    }
}

impl crate::BSplineSurface {
    /// Checks whole-row carrier identity and orientation without certifying
    /// the surface interior or its embedding.
    pub fn boundary_orientation(
        &self,
        side: SplineBoundary,
        curve: &Curve,
        dec: &mut Decider,
    ) -> Result<Option<bool>, Refusal> {
        let boundary = BSplineCurve::boundary(self, side);
        match curve {
            Curve::BSpline(curve) => {
                if curve.is_boundary_of(self, side)
                    || normalized_equal(curve, &boundary, false, dec)?
                {
                    return Ok(Some(true));
                }
                if normalized_equal(curve, &boundary, true, dec)? {
                    return Ok(Some(false));
                }
                Ok(None)
            }
            Curve::Line(line) => {
                if exact_sign(&line.dir.norm_squared(), dec)? != Sign::Positive {
                    return Ok(None);
                }
                let d = [&line.dir.x, &line.dir.y, &line.dir.z];
                for (span, _) in boundary.spans() {
                    let p = boundary.span_polynomials(span, dec)?;
                    let radial = vector_polynomial(&p, &line.origin);
                    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
                        if !(&radial[a].scale(d[b]) - &radial[b].scale(d[a])).is_zero(dec)? {
                            return Ok(None);
                        }
                    }
                }
                let [lo, hi] = boundary.domain().map(|p| Real::from_rational(p.clone()));
                let a = boundary.point_at(&lo, dec)?;
                let b = boundary.point_at(&hi, dec)?;
                let direction = exact_sign(&b.minus(&a).dot(&line.dir), dec)?;
                Ok((direction != Sign::Zero).then_some(direction == Sign::Positive))
            }
            Curve::Circle(circle) => conic_orientation(
                &boundary,
                &circle.center,
                &circle.normal,
                None,
                Some(&circle.radius_sq),
                dec,
            ),
            Curve::Ellipse(ellipse) => conic_orientation(
                &boundary,
                &ellipse.center,
                &ellipse.normal,
                Some((&ellipse.major_dir, &ellipse.major_sq, &ellipse.minor_sq)),
                None,
                dec,
            ),
            Curve::Intersection(_) | Curve::PatchSection(_) => Err(unproven(
                "an intersection-defined edge needs a finite boundary correspondence proof",
            )),
        }
    }
}
