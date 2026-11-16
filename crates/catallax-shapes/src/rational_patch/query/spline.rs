//! Finite spline points on graphs with a proved affine projection inverse.

use super::*;
use crate::{BSplineCurve, SectionHits};

struct Pullback {
    parameters: [RealPoly; 2],
    denominator: RealPoly,
    equation: RealPoly,
    domain: Domain,
    contained: bool,
}

impl Pullback {
    fn new(span: &Span, curve: &[RealPoly; 4], dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        let Some(frame) = affine::frame(span, dec)? else {
            return Ok(None);
        };
        let rhs: [RealPoly; 2] = std::array::from_fn(|row| {
            curve[frame.axes[row]].scale(&frame.weight) - curve[3].scale(&frame.coordinates[row][0])
        });
        let [a, b] = &frame.coordinates;
        let parameters = [
            rhs[0].scale(&b[2]) - rhs[1].scale(&a[2]),
            rhs[1].scale(&a[1]) - rhs[0].scale(&b[1]),
        ];
        let denominator = curve[3].scale(&frame.determinant);
        let remaining = (0..3).find(|axis| !frame.axes.contains(axis)).unwrap();
        let terms = span.coordinates[remaining].terms();
        let degree = terms
            .iter()
            .map(|(m, _)| m.exponents().iter().sum::<u32>())
            .max()
            .unwrap_or(0);
        let mut value = RealPoly::zero(1);
        for (m, c) in terms {
            let [u, v] = [m.exponents()[0], m.exponents()[1]];
            let term = &parameters[0].pow_u32(u) * &parameters[1].pow_u32(v);
            value = value + (&term * &denominator.pow_u32(degree - u - v)).scale(&c);
        }
        let equation = &value * &curve[3]
            - (&curve[remaining] * &denominator.pow_u32(degree)).scale(&frame.weight);
        let contained = equation.is_zero(dec)?;
        Ok(Some(Self {
            parameters,
            denominator,
            equation,
            domain: span.domain.clone(),
            contained,
        }))
    }

    fn inside(&self, parameter: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
        let denominator = self.denominator.eval_real(std::slice::from_ref(parameter));
        for (value, [lo, hi]) in self.parameters.iter().zip(&self.domain) {
            let at = &value.eval_real(std::slice::from_ref(parameter)) / &denominator;
            if exact_sign(&(&at - &Real::from(lo.clone())), dec)? == Sign::Negative
                || exact_sign(&(&Real::from(hi.clone()) - &at), dec)? == Sign::Negative
            {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

fn insert(parameters: &mut Vec<Real>, value: Real, dec: &mut Decider) -> Result<(), Refusal> {
    for index in 0..parameters.len() {
        match exact_sign(&(&value - &parameters[index]), dec)? {
            Sign::Zero => return Ok(()),
            Sign::Negative => {
                parameters.insert(index, value);
                return Ok(());
            }
            Sign::Positive => (),
        }
    }
    parameters.push(value);
    Ok(())
}

fn roots(
    poly: &RealPoly,
    domain: [&Rational; 2],
    parameters: &mut Vec<Real>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if poly.is_zero(dec)? {
        return Ok(());
    }
    let precision = dec.policy().max_precision_bits.min(128);
    let local = if precision == 0 {
        CellRoots::Unknown
    } else {
        monotone_root(&mut CellPolynomial::new(poly, precision, dec)?, domain, dec)?
    };
    let roots = match local {
        CellRoots::Empty => return Ok(()),
        CellRoots::One(root) => vec![root],
        CellRoots::Unknown | CellRoots::Split => poly.real_roots(dec)?,
    };
    for parameter in roots {
        if exact_sign(&(&parameter - &Real::from(domain[0].clone())), dec)? != Sign::Negative
            && exact_sign(&(&Real::from(domain[1].clone()) - &parameter), dec)? != Sign::Negative
        {
            insert(parameters, parameter, dec)?;
        }
    }
    Ok(())
}

pub(crate) enum SplineParameterHits {
    Contained,
    Points(Vec<Real>),
}

pub(crate) fn spline_parameters(
    patches: &[&RationalPatch],
    spline: &BSplineCurve,
    dec: &mut Decider,
) -> Result<Option<SplineParameterHits>, Refusal> {
    let mut unique: Vec<&RationalPatch> = Vec::new();
    for &patch in patches {
        let mut duplicate = false;
        for &previous in &unique {
            dec.charge_work(1)?;
            if std::ptr::eq(patch, previous) {
                duplicate = true;
                break;
            }
        }
        if !duplicate {
            unique.push(patch);
        }
    }
    let mut parameters = Vec::new();
    let (mut any_contained, mut all_contained) = (false, true);
    for (index, domain) in spline.spans() {
        let curve = spline.query_span_polynomials(index, dec)?;
        let mut pullbacks = Vec::new();
        for patch in &unique {
            for span in &patch.data.spans {
                let Some(pullback) = Pullback::new(span, &curve, dec)? else {
                    return Ok(None);
                };
                pullbacks.push(pullback);
            }
        }
        let mut partition = vec![Real::from(domain[0].clone()), Real::from(domain[1].clone())];
        for pullback in &pullbacks {
            if pullback.contained {
                for (value, bounds) in pullback.parameters.iter().zip(&pullback.domain) {
                    for bound in bounds {
                        roots(
                            &(value - &pullback.denominator.scale(&Real::from(bound.clone()))),
                            domain,
                            &mut partition,
                            dec,
                        )?;
                    }
                }
            } else {
                let mut candidates = Vec::new();
                roots(&pullback.equation, domain, &mut candidates, dec)?;
                // Each root already belongs to this equation. Keep that
                // ownership while checking its finite source domain.
                for parameter in candidates {
                    if pullback.inside(&parameter, dec)? {
                        insert(&mut parameters, parameter, dec)?;
                    }
                }
            }
        }
        // Complete domain-boundary roots make membership constant inside
        // each open interval; the midpoint selects that proved interval.
        for interval in partition.windows(2) {
            let middle = (&interval[0] + &interval[1]) / Real::from_i64(2);
            let mut contained = false;
            for pullback in pullbacks.iter().filter(|p| p.contained) {
                contained |= pullback.inside(&middle, dec)?;
            }
            any_contained |= contained;
            all_contained &= contained;
        }
        for parameter in partition {
            for pullback in pullbacks.iter().filter(|p| p.contained) {
                if pullback.inside(&parameter, dec)? {
                    insert(&mut parameters, parameter, dec)?;
                    break;
                }
            }
        }
    }
    if all_contained {
        return Ok(Some(SplineParameterHits::Contained));
    }
    if any_contained {
        return Err(unproven(
            "partial spline containment needs finite interval ownership",
        ));
    }
    Ok(Some(SplineParameterHits::Points(parameters)))
}

pub(crate) fn spline_points(
    patches: &[&RationalPatch],
    spline: &BSplineCurve,
    dec: &mut Decider,
) -> Result<Option<SectionHits>, Refusal> {
    let Some(answer) = spline_parameters(patches, spline, dec)? else {
        return Ok(None);
    };
    let SplineParameterHits::Points(parameters) = answer else {
        return Ok(Some(SectionHits::Contained));
    };
    let mut points: Vec<Point3> = Vec::new();
    for parameter in parameters {
        let point = spline.point_at(&parameter, dec)?;
        let mut duplicate = false;
        for previous in &points {
            duplicate |= exact_sign(&point.minus(previous).norm_squared(), dec)? == Sign::Zero;
        }
        if !duplicate {
            points.push(point)
        }
    }
    Ok(Some(SectionHits::Points(points)))
}

impl RationalPatch {
    /// Complete finite points through a proved affine projection. No contact
    /// label is inferred. Unsupported projections return `None`.
    pub fn spline_hits(
        &self,
        spline: &BSplineCurve,
        dec: &mut Decider,
    ) -> Result<Option<SectionHits>, Refusal> {
        spline_points(&[self], spline, dec)
    }
}
