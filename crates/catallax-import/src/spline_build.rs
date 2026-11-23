//! Exact admission checks for finite STEP spline geometry.

use catallax_numbers::{
    Ball, Decider, Outcome, Rational, Real, RealPoly, Refusal, StrictSign, ZeroProof,
    current_operation_work,
};
use catallax_shapes::{BSplineCurve, Curve, Line, RPoint3, RVector3, RationalPatch, Surface};

use crate::{Face, ImportRefusal, OrientedEdge, SourceFormat, StepSplineSurface};

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
    face: &Face,
    source: &StepSplineSurface,
    dec: &mut Decider,
) -> Result<Surface, ImportRefusal> {
    let carrier = &source.geometry;
    let result = (|| {
        let trimmed = rectangular_trim(face, source, dec)?
            .map(|domain| carrier.extracted(domain.each_ref().map(|axis| axis.each_ref()), dec))
            .transpose()?;
        RationalPatch::new(trimmed.as_ref().unwrap_or(carrier).to_surface(dec)?, dec)
    })();
    let support =
        result.map_err(|error| refused(&face.name, "ADVANCED_FACE", "support admission", error))?;
    Ok(Surface::BSpline(support.into()))
}

// A nonplanar bilinear carrier has unique coordinates in the basis
// (P10-P00, P01-P00, P11-P10-P01+P00). A straight physical edge
// can follow an iso-row exactly. The four full rows are checked later
// against the clipped carrier by RationalPatchFace.
fn rectangular_trim(
    face: &Face,
    source: &StepSplineSurface,
    dec: &mut Decider,
) -> Result<Option<[[Rational; 2]; 2]>, Refusal> {
    use crate::exact::{cross, dot, sub};
    if source.source_format != SourceFormat::Step {
        return Ok(None);
    }
    let [ring] = face.loops.as_slice() else {
        return Ok(None);
    };
    let carrier = &source.geometry;
    let [p00, p01, p10, p11] = carrier.controls() else {
        return Ok(None);
    };
    if carrier.degrees() != [1, 1]
        || carrier.weights().iter().any(|w| w != &carrier.weights()[0])
        || ring.edges.len() != 4
        || !ring.edges.iter().all(|edge| match &edge.curve {
            crate::Curve::Line(..) => true,
            crate::Curve::BSpline(spline) => {
                spline.geometry.degree() == 1 && spline.geometry.controls().len() == 2
            }
            _ => false,
        })
        || ring
            .edges
            .iter()
            .all(|edge| carrier.controls().contains(&edge.start.point))
    {
        return Ok(None);
    }

    let limbs: u64 = carrier
        .controls()
        .iter()
        .chain(
            ring.edges
                .iter()
                .flat_map(|edge| [&edge.start.point, &edge.end.point]),
        )
        .flatten()
        .map(Rational::limbs)
        .max()
        .unwrap_or(1);
    // The only unbounded arithmetic here is on the input rational limbs.
    // There are four edges and a fixed number of dot and cross products.
    dec.charge_work(limbs.saturating_mul(limbs).saturating_mul(128))?;

    Ok((|| {
        let u = sub(p10, p00);
        let v = sub(p01, p00);
        let uv = sub(&sub(p11, p10), &v);
        let dual_u = cross(&v, &uv);
        let dual_v = cross(&uv, &u);
        let dual_uv = cross(&u, &v);
        let volume = dot(&u, &dual_u);
        if volume.is_zero() {
            return None;
        }
        let coordinate = |p: &crate::Vec3| {
            let q = sub(p, p00);
            let a = &dot(&q, &dual_u) / &volume;
            let b = &dot(&q, &dual_v) / &volume;
            let in_range = |x: &Rational| Rational::zero() <= *x && *x <= Rational::one();
            (in_range(&a) && in_range(&b) && dot(&q, &dual_uv) == &(&a * &b) * &volume)
                .then_some([a, b])
        };
        let mut vertices = std::collections::BTreeSet::new();
        for edge in &ring.edges {
            let a = coordinate(&edge.start.point)?;
            let b = coordinate(&edge.end.point)?;
            if (a[0] == b[0]) == (a[1] == b[1]) {
                return None;
            }
            vertices.extend([a, b]);
        }
        if vertices.len() != 4 {
            return None;
        }
        let unit_bounds: [[Rational; 2]; 2] = std::array::from_fn(|axis| {
            let values = || vertices.iter().map(|uv| &uv[axis]);
            [
                values().min().unwrap().clone(),
                values().max().unwrap().clone(),
            ]
        });
        if unit_bounds.iter().any(|[lo, hi]| lo >= hi)
            || !vertices
                .iter()
                .all(|uv| (0..2).all(|axis| unit_bounds[axis].contains(&uv[axis])))
        {
            return None;
        }
        let domain = carrier.domain();
        let trimmed = std::array::from_fn(|axis| {
            let [lo, hi] = domain[axis];
            unit_bounds[axis].each_ref().map(|t| lo + &(t * &(hi - lo)))
        });
        (trimmed.each_ref().map(|axis| axis.each_ref()) != domain).then_some(trimmed)
    })())
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
    embedding(&curve, &source.geometry, dec)
        .map_err(|error| refused(&edge.edge, "EDGE_CURVE", "curve embedding", error))?;
    let controls = source.geometry.controls();
    let origin = &controls[0];
    let direction: [Rational; 3] =
        std::array::from_fn(|i| &controls[controls.len() - 1][i] - &origin[i]);
    let straight = controls.iter().all(|p| {
        let d: [Rational; 3] = std::array::from_fn(|i| &p[i] - &origin[i]);
        crate::exact::is_zero(&crate::exact::cross(&d, &direction))
    });
    if straight && !crate::exact::is_zero(&direction) {
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
    source: &catallax_shapes::spline_tools::RationalSplineCurve,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    match embedding_fast(curve, source.controls(), dec) {
        Err(error) if error.code == "import.spline-embedding-unproven" => {
            crate::spline_embedding::prove(source, dec).map_err(|_| error)
        }
        result => result,
    }
}

fn embedding_fast(
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
    fn source_edge(controls: Vec<[Rational; 3]>) -> OrientedEdge {
        let geometry = RationalSplineCurve::new(
            3,
            [0, 0, 0, 0, 1, 1, 1, 1].map(Rational::from_i64).to_vec(),
            controls.clone(),
            vec![Rational::one(); 4],
            &Decider::default(),
        )
        .unwrap();
        OrientedEdge {
            edge: "synthetic".into(),
            start: crate::Vertex {
                name: "start".into(),
                point: controls[0].clone(),
            },
            end: crate::Vertex {
                name: "end".into(),
                point: controls[3].clone(),
            },
            curve_same_sense: true,
            forward: true,
            curve: crate::Curve::BSpline(crate::StepSplineCurve {
                geometry,
                closed: Some(false),
                self_intersect: Some(false),
            }),
        }
    }
    #[test]
    fn quarter_circle_uses_the_stronger_certificate_after_the_fast_check_refuses() {
        let edge = source_edge(
            [
                ["1", "0", "0"],
                ["1", "0.55228475", "0"],
                ["0.55228475", "1", "0"],
                ["0", "1", "0"],
            ]
            .map(|p| p.map(|v| Rational::parse(v).unwrap()))
            .to_vec(),
        );
        let crate::Curve::BSpline(s) = &edge.curve else {
            unreachable!()
        };
        let carrier = s.geometry.to_curve(&mut Decider::default()).unwrap();
        assert_eq!(
            embedding_fast(&carrier, s.geometry.controls(), &Decider::default())
                .unwrap_err()
                .code,
            "import.spline-embedding-unproven"
        );
        assert!(matches!(
            curve(&edge, &mut Decider::default()).unwrap(),
            Curve::BSpline(_)
        ));
    }
    #[test]
    fn a_true_self_intersection_keeps_the_typed_embedding_refusal() {
        let edge = source_edge(
            [
                ["2/9", "-1/9", "0"],
                ["-1/9", "7/54", "0"],
                ["-1/9", "-7/54", "0"],
                ["2/9", "1/9", "0"],
            ]
            .map(|p| p.map(|v| Rational::parse(v).unwrap()))
            .to_vec(),
        );
        let crate::Curve::BSpline(source) = &edge.curve else {
            unreachable!()
        };
        for t in [
            Rational::from_ratio(1, 3).unwrap(),
            Rational::from_ratio(2, 3).unwrap(),
        ] {
            let mut points = source.geometry.controls().to_vec();
            let complement = &Rational::one() - &t;
            while points.len() > 1 {
                points = points
                    .windows(2)
                    .map(|p| std::array::from_fn(|i| &(&p[0][i] * &complement) + &(&p[1][i] * &t)))
                    .collect();
            }
            assert!(points[0].iter().all(Rational::is_zero));
        }
        let error = curve(&edge, &mut Decider::default()).unwrap_err();
        assert!(
            matches!(error,ImportRefusal::Unsupported{ref type_name,..} if type_name.contains("import.spline-embedding-unproven")),
            "{error:?}"
        );
    }
    #[test]
    fn a_successful_fast_embedding_check_keeps_its_work_count() {
        let edge = source_edge(
            [
                ["0", "0", "0"],
                ["1", "0", "0"],
                ["2", "0", "0"],
                ["3", "0", "0"],
            ]
            .map(|p| p.map(|v| Rational::parse(v).unwrap()))
            .to_vec(),
        );
        let crate::Curve::BSpline(s) = &edge.curve else {
            unreachable!()
        };
        let carrier = s.geometry.to_curve(&mut Decider::default()).unwrap();
        let (old, old_work) = catallax_numbers::with_operation_budget(20_000_000, || {
            embedding_fast(&carrier, s.geometry.controls(), &Decider::default())
        });
        let (new, new_work) = catallax_numbers::with_operation_budget(20_000_000, || {
            embedding(&carrier, &s.geometry, &mut Decider::default())
        });
        assert!(old.is_ok());
        assert!(new.is_ok());
        assert_eq!(old_work.used, new_work.used);
        assert_eq!(old_work.exhausted, new_work.exhausted);
    }
    #[test]
    fn an_unclamped_arc_with_repeated_outer_controls_keeps_its_spline_carrier() {
        let geometry = RationalSplineCurve::new(
            3,
            (0..8).map(Rational::from_i64).collect(),
            [[0, 0, 0], [-1, 0, 0], [1, 1, 0], [0, 0, 0]]
                .map(|p| p.map(Rational::from_i64))
                .to_vec(),
            vec![Rational::one(); 4],
            &Decider::default(),
        )
        .unwrap();
        let start = ["-1/2", "1/6", "0"].map(|v| Rational::parse(v).unwrap());
        let end = ["1/2", "2/3", "0"].map(|v| Rational::parse(v).unwrap());
        let edge = OrientedEdge {
            edge: "unclamped".into(),
            start: crate::Vertex {
                name: "start".into(),
                point: start.clone(),
            },
            end: crate::Vertex {
                name: "end".into(),
                point: end.clone(),
            },
            curve_same_sense: true,
            forward: true,
            curve: crate::Curve::BSpline(crate::StepSplineCurve {
                geometry,
                closed: Some(false),
                self_intersect: Some(false),
            }),
        };
        let mut dec = Decider::default();
        let Curve::BSpline(carrier) = curve(&edge, &mut dec).unwrap() else {
            panic!("curved source collapsed to a degenerate line")
        };
        for (endpoint, parameter) in [(start, 3), (end, 4)] {
            let actual = carrier
                .unique_parameter_of(&point(&endpoint), &mut dec)
                .unwrap()
                .unwrap();
            assert_eq!(
                actual.value().as_rational(),
                Some(Rational::from_i64(parameter))
            );
        }
    }
}
