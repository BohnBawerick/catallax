//! Interior witnesses along exact finite rational-chart paths.

use super::*;
use catallax_shapes::{BSplineCurve, RationalPatch, SplineSupport};

fn compose(polynomial: &RealPoly, parameters: &[P; 2], dec: &Decider) -> Result<P, Refusal> {
    let mut value = P::constant(Real::zero());
    for (monomial, coefficient) in polynomial.terms() {
        dec.check_deadline()?;
        let mut term = P::constant(coefficient);
        for (parameter, &power) in parameters.iter().zip(monomial.exponents()) {
            term = term.mul(&parameter.pow(power));
        }
        value = value.add(&term);
    }
    Ok(value)
}

fn inside(
    at: &[Real; 2],
    domain: &[[Rational; 2]; 2],
    strict: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for (value, [lo, hi]) in at.iter().zip(domain) {
        for margin in [
            value - &Real::from(lo.clone()),
            &Real::from(hi.clone()) - value,
        ] {
            let side = coefficient_sign(&margin, dec)?;
            if side == Sign::Negative || (strict && side == Sign::Zero) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn leaves_left(path: &Path, left: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    let w = path.denominator.coefficient(0);
    let dw = path.denominator.coefficient(1);
    let [x, y, z] = path
        .numerators
        .each_ref()
        .map(|coordinate| &coordinate.coefficient(1) * &w - &coordinate.coefficient(0) * &dw);
    Ok(coefficient_sign(&left.dot(&Vector3::new(x, y, z)), dec)? == Sign::Positive)
}

/// Positive spline weights put the complete curve in its control hull.
/// Separation from one hull coordinate interval excludes the whole curve,
/// including an allowed contact at the path's initial point.
fn clears_native(
    path: &Path,
    q: &Point3,
    curve: &BSplineCurve,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for axis in 0..3 {
        let mut controls = curve.controls();
        let first = controls.next().expect("a native spline has controls");
        let mut lo = point(first)[axis].clone();
        let mut hi = lo.clone();
        for control in controls {
            let value = point(control)[axis];
            if coefficient_sign(&(value - &lo), dec)? == Sign::Negative {
                lo = value.clone();
            }
            if coefficient_sign(&(value - &hi), dec)? == Sign::Positive {
                hi = value.clone();
            }
        }
        for (bound, outside) in [(lo, Sign::Negative), (hi, Sign::Positive)] {
            if coefficient_sign(&(point(q)[axis] - &bound), dec)? != outside {
                continue;
            }
            let separation = path.numerators[axis].sub(&path.denominator.scale(bound));
            if separation.excludes_zero(dec)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn clears(
    model: &Model,
    path: &Path,
    q: &Point3,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let curve_id = model.edge(edge).curve;
    match model.geom().curve(curve_id) {
        Curve::BSpline(curve) => clears_native(path, q, curve, dec),
        Curve::PatchSection(branch) => {
            let equation = branch.surface().surface().implicit().ok_or_else(|| {
                crate::geometry::refusal(
                    "boolean.spline-face-probe",
                    "a retained section needs an exact cutting equation for path clearance",
                )
            })?;
            numerator(&equation, path, dec)?.excludes_zero(dec)
        }
        _ => {
            let Some(equations) = model.geom().curve_equations(curve_id) else {
                return Ok(false);
            };
            for equation in equations {
                if numerator(equation, path, dec)?.excludes_zero(dec)? {
                    return Ok(true);
                }
            }
            Ok(false)
        }
    }
}

fn on_patch(
    model: &Model,
    patch: &RationalPatch,
    p: &Point3,
    left: &Vector3,
    step: &Real,
    barriers: &mut [EdgeId],
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    let Some(at) = patch.parameters_at(p, dec)? else {
        return Ok(None);
    };
    let [du, dv] = patch.derivatives_at(&at, dec)?;
    let a = du.norm_squared();
    let b = du.dot(&dv);
    let c = dv.norm_squared();
    let determinant = &a * &c - b.square();
    if coefficient_sign(&determinant, dec)? != Sign::Positive {
        return Ok(None);
    }
    let u = du.dot(left);
    let v = dv.dot(left);
    let velocity = [
        &(&c * &u - &b * &v) / &determinant,
        &(&a * &v - &b * &u) / &determinant,
    ];
    let to = std::array::from_fn(|axis| &at[axis] + &(&velocity[axis] * step));
    // The open segment stays inside this one checked chart and one knot span.
    // It needs no seam choice, extrapolation or globally injective atlas map.
    if !inside(&to, patch.domain(), true, dec)? {
        return Ok(None);
    }
    let carrier = patch.carrier();
    let knots = carrier.knots();
    for span in carrier.span_indices() {
        let domain = std::array::from_fn(|axis| {
            [
                knots[axis][span[axis]].clone(),
                knots[axis][span[axis] + 1].clone(),
            ]
        });
        if !inside(&at, &domain, false, dec)? || !inside(&to, &domain, true, dec)? {
            continue;
        }
        let parameters = std::array::from_fn(|axis| P::line(at[axis].clone(), to[axis].clone()));
        let homogeneous = carrier.span_polynomials(span, dec)?;
        let mut coordinates = Vec::with_capacity(4);
        for coordinate in &homogeneous {
            coordinates.push(compose(coordinate, &parameters, dec)?);
        }
        let [x, y, z, denominator]: [P; 4] = coordinates
            .try_into()
            .expect("four homogeneous coordinates");
        let path = Path {
            numerators: [x, y, z],
            denominator,
            axial: None,
            cylindrical: None,
        };
        if !leaves_left(&path, left, dec)? {
            continue;
        }
        let q = patch.point_at(&to, dec)?;
        let mut blocked = None;
        for (index, &edge) in barriers.iter().enumerate() {
            if !clears(model, &path, &q, edge, dec)? {
                blocked = Some(index);
                break;
            }
        }
        match blocked {
            None => return Ok(Some(q)),
            Some(index) => barriers[..=index].rotate_right(1),
        }
    }
    Ok(None)
}

/// A point on the finite carrier, reached strictly left of the source arc by
/// a complete exact path avoiding every boundary edge for 0 < t <= 1.
pub(crate) fn probe(
    model: &Model,
    surface: &Surface,
    p: &Point3,
    left: &Vector3,
    step: &Real,
    barriers: &mut [EdgeId],
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    let before_separation = dec.stats().coincident_separation;
    let result = match surface {
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
            on_patch(model, patch, p, left, step, barriers, dec)?
        }
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => {
            let mut found = None;
            for chart in atlas.charts() {
                if let Some(candidate) =
                    on_patch(model, chart.patch(), p, left, step, barriers, dec)?
                {
                    found = Some(candidate);
                    break;
                }
            }
            found
        }
        _ => None,
    };
    if result.is_some() && dec.stats().coincident_separation != before_separation {
        return Err(crate::geometry::refusal(
            "boolean.spline-face-probe",
            "a spline interior path needs exact identities, not separation-bound coincidences",
        ));
    }
    Ok(result)
}
