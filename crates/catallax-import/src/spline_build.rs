//! Exact admission checks for finite STEP spline geometry.

use catallax_numbers::{
    Ball, Decider, Outcome, Rational, Real, RealPoly, Refusal, StrictSign, ZeroProof,
    current_operation_work,
};
use catallax_shapes::{BSplineCurve, Curve, Line, RPoint3, RVector3, RationalPatch, Surface};

use crate::{ImportRefusal, OrientedEdge, StepSplineSurface};

pub(crate) fn refused(entity: &str, kind: &str, step: &str, error: Refusal) -> ImportRefusal {
    if current_operation_work().exhausted {
        return ImportRefusal::ResourceLimit;
    }
    ImportRefusal::Unsupported {
        entity: entity.into(),
        type_name: format!("{kind} (spline {step}: {}: {})", error.code, error.detail),
    }
}

pub(crate) fn surface(
    entity: &str,
    source: &StepSplineSurface,
    dec: &mut Decider,
) -> Result<Surface, ImportRefusal> {
    let support = source
        .geometry
        .to_surface(dec)
        .and_then(|carrier| RationalPatch::new(carrier, dec))
        .map_err(|error| refused(entity, "ADVANCED_FACE", "support admission", error))?;
    Ok(Surface::BSpline(support.into()))
}

fn point(p: &[Rational; 3]) -> catallax_shapes::Point3 {
    RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()).to_real()
}

pub(crate) fn curve(edge: &OrientedEdge, dec: &mut Decider) -> Result<Curve, ImportRefusal> {
    let crate::Curve::BSpline(source) = &edge.curve else {
        unreachable!()
    };
    let curve = source
        .geometry
        .to_curve(dec)
        .map_err(|error| refused(&edge.edge, "EDGE_CURVE", "curve admission", error))?;
    let mut parameters = Vec::new();
    for endpoint in [&edge.start.point, &edge.end.point] {
        let parameter = curve
            .unique_parameter_of(&point(endpoint), dec)
            .map_err(|error| refused(&edge.edge, "EDGE_CURVE", "endpoint incidence", error))?;
        if let Some(parameter) = parameter {
            parameters.push(parameter.value().clone());
        } else {
            return Err(ImportRefusal::NotSolid {
                detail: format!(
                    "spline endpoint incidence: a vertex of edge {} is outside its finite curve",
                    edge.edge
                ),
            });
        }
    }
    let order = dec.sign(&(&parameters[1] - &parameters[0]));
    let expected = if edge.curve_same_sense {
        StrictSign::Positive
    } else {
        StrictSign::Negative
    };
    match order {
        Outcome::Proven(actual) if actual == expected => (),
        Outcome::ProvenCoincident(proof) if proof.proof == ZeroProof::SeparationBound => {
            return Err(refused(
                &edge.edge,
                "EDGE_CURVE",
                "endpoint order",
                Refusal::new(
                    "numbers.geometric-gap-unproven",
                    "finite spline endpoint order needs an exact sign",
                ),
            ));
        }
        Outcome::Proven(_) | Outcome::ProvenCoincident(_) => {
            return Err(ImportRefusal::NotSolid {
                detail: format!(
                    "spline endpoint order: edge {} disagrees with its finite curve direction",
                    edge.edge
                ),
            });
        }
        Outcome::NotExact(error) => {
            return Err(refused(
                &edge.edge,
                "EDGE_CURVE",
                "endpoint order",
                Refusal::from(error),
            ));
        }
        Outcome::Refused(error) => {
            return Err(refused(&edge.edge, "EDGE_CURVE", "endpoint order", error));
        }
    }
    embedding(&curve, source.geometry.controls(), dec)
        .map_err(|error| refused(&edge.edge, "EDGE_CURVE", "curve embedding", error))?;
    let controls = source.geometry.controls();
    let origin = &controls[0];
    let direction: [Rational; 3] =
        std::array::from_fn(|i| &controls[controls.len() - 1][i] - &origin[i]);
    let straight = controls.iter().all(|p| {
        let d: [Rational; 3] = std::array::from_fn(|i| &p[i] - &origin[i]);
        crate::exact::is_zero(&crate::exact::cross(&d, &direction))
    });
    if straight {
        Ok(Curve::Line(Line {
            origin: point(origin),
            dir: RVector3::new(
                direction[0].clone(),
                direction[1].clone(),
                direction[2].clone(),
            )
            .to_real(),
        }))
    } else {
        Ok(Curve::BSpline(curve))
    }
}

// One strictly monotone projection over every span proves a regular embedded arc.
fn embedding(
    curve: &BSplineCurve,
    controls: &[[Rational; 3]],
    dec: &Decider,
) -> Result<(), Refusal> {
    let bits = dec.policy().max_precision_bits.min(128);
    if bits == 0 {
        return Err(Refusal::new(
            "numbers.precision-limit",
            "no spline embedding working precision",
        ));
    }
    let direction: [Rational; 3] =
        std::array::from_fn(|i| &controls[controls.len() - 1][i] - &controls[0][i]);
    for (span, interval) in curve.spans() {
        let p = curve.span_polynomials(span, dec)?;
        let mut numerator = RealPoly::zero(1);
        for (coordinate, direction) in p[..3].iter().zip(&direction) {
            numerator = numerator
                + (&coordinate.partial(0) * &p[3] - coordinate * &p[3].partial(0))
                    .scale(&Real::from_rational(direction.clone()));
        }
        let parameter =
            Ball::from_rational(interval[0], bits).hull(&Ball::from_rational(interval[1], bits));
        if numerator.eval_ball(&[parameter], dec)?.lo() <= Rational::zero() {
            return Err(Refusal::new(
                "import.spline-embedding-unproven",
                "no strictly monotone endpoint-chord projection over every knot span",
            ));
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Incidence {
    OffSurface,
    FiniteEdge,
    WholeCurve,
}

pub(crate) fn on_surface(
    curve: &Curve,
    surface: &Surface,
    edge: &OrientedEdge,
    dec: &mut Decider,
) -> Result<Incidence, Refusal> {
    if matches!(surface, Surface::Cone(_)) {
        return Err(Refusal::new(
            "import.spline-cone-nappe-unproven",
            "implicit incidence alone does not prove the retained cone nappe",
        ));
    }
    let equation = surface.implicit().ok_or_else(|| {
        Refusal::new(
            "import.spline-incidence-unproven",
            "the finite spline support needs a whole-boundary proof",
        )
    })?;
    match curve {
        Curve::BSpline(curve) => {
            let parameters = [&edge.start.point, &edge.end.point].map(|endpoint| {
                curve
                    .unique_parameter_of(&point(endpoint), dec)
                    .and_then(|parameter| {
                        parameter.map(|p| p.value().clone()).ok_or_else(|| {
                            Refusal::new(
                                "import.spline-endpoint-unproven",
                                "a finite spline endpoint has no preimage",
                            )
                        })
                    })
            });
            let [a, b] = parameters;
            let [a, b] = [a?, b?];
            let (lo, hi) = if edge.curve_same_sense {
                (a, b)
            } else {
                (b, a)
            };
            let mut whole_curve = true;
            for (span, [start, end]) in curve.spans() {
                if !positive(&(&hi - &Real::from_rational(start.clone())), dec)?
                    || !positive(&(&Real::from_rational(end.clone()) - &lo), dec)?
                {
                    whole_curve = false;
                    continue;
                }
                if !substitute(&equation, &curve.span_polynomials(span, dec)?).is_zero(dec)? {
                    return Ok(Incidence::OffSurface);
                }
            }
            Ok(if whole_curve {
                Incidence::WholeCurve
            } else {
                Incidence::FiniteEdge
            })
        }
        Curve::Line(line) => {
            let origin = [&line.origin.x, &line.origin.y, &line.origin.z];
            let direction = [&line.dir.x, &line.dir.y, &line.dir.z];
            let p = std::array::from_fn(|i| {
                if i == 3 {
                    RealPoly::constant(1, Real::one())
                } else {
                    RealPoly::constant(1, origin[i].clone())
                        + RealPoly::var(1, 0).scale(direction[i])
                }
            });
            Ok(if substitute(&equation, &p).is_zero(dec)? {
                Incidence::WholeCurve
            } else {
                Incidence::OffSurface
            })
        }
        _ => Err(Refusal::new(
            "import.spline-incidence-unproven",
            "unsupported spline boundary carrier",
        )),
    }
}

fn substitute(equation: &RealPoly, coordinates: &[RealPoly; 4]) -> RealPoly {
    let terms = equation.terms();
    let degree = terms
        .iter()
        .map(|(m, _)| m.exponents().iter().sum::<u32>())
        .max()
        .unwrap_or(0);
    terms.into_iter().fold(RealPoly::zero(1), |sum, (m, c)| {
        let mut term = RealPoly::constant(1, c);
        let mut weight_power = degree;
        for (i, &power) in m.exponents().iter().enumerate() {
            term = term * coordinates[i].pow_u32(power);
            weight_power -= power;
        }
        sum + term * coordinates[3].pow_u32(weight_power)
    })
}

fn positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(true),
        Outcome::Proven(StrictSign::Negative) => Ok(false),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(false),
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            "numbers.geometric-gap-unproven",
            "a finite spline span needs exact trim ordering",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Plane, spline_tools::RationalSplineCurve};

    #[test]
    fn incidence_checks_only_spans_crossed_by_the_finite_edge() {
        let mut dec = Decider::default();
        let controls = [[0, 0, 0], [1, 0, 0], [2, 1, 0]].map(|p| p.map(Rational::from_i64));
        let source = RationalSplineCurve::new(
            1,
            [0, 0, 1, 2, 2].map(Rational::from_i64).to_vec(),
            controls.to_vec(),
            vec![Rational::one(); 3],
            &dec,
        )
        .unwrap();
        let mut edge = OrientedEdge {
            edge: "synthetic".into(),
            start: crate::Vertex {
                name: "start".into(),
                point: controls[0].clone(),
            },
            end: crate::Vertex {
                name: "end".into(),
                point: controls[1].clone(),
            },
            curve_same_sense: true,
            forward: true,
            curve: crate::Curve::BSpline(crate::StepSplineCurve {
                geometry: source,
                closed: Some(false),
                self_intersect: Some(false),
            }),
        };
        let carrier = curve(&edge, &mut dec).unwrap();
        let plane = Surface::Plane(Plane {
            origin: point(&controls[0]),
            normal: RVector3::new(Rational::zero(), Rational::one(), Rational::zero()).to_real(),
            ref_dir: RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
        });
        assert_eq!(
            on_surface(&carrier, &plane, &edge, &mut dec).unwrap(),
            Incidence::FiniteEdge
        );
        edge.end.point = controls[2].clone();
        assert_eq!(
            on_surface(&carrier, &plane, &edge, &mut dec).unwrap(),
            Incidence::OffSurface
        );
        let global_plane = Surface::Plane(Plane {
            origin: point(&controls[0]),
            normal: RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
            ref_dir: RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
        });
        assert_eq!(
            on_surface(&carrier, &global_plane, &edge, &mut dec).unwrap(),
            Incidence::WholeCurve
        );
    }
}
