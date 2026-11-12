//! Interior witnesses along exact finite rational-chart paths.

use super::*;
use catallax_shapes::{BSplineCurve, RationalPatch, SplineBoundary, SplineSupport};

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
    strict: [bool; 4],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for (axis, (value, [lo, hi])) in at.iter().zip(domain).enumerate() {
        for (boundary, margin) in [
            value - &Real::from(lo.clone()),
            &Real::from(hi.clone()) - value,
        ]
        .into_iter()
        .enumerate()
        {
            let side = coefficient_sign(&margin, dec)?;
            if side == Sign::Negative || (strict[2 * axis + boundary] && side == Sign::Zero) {
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
    (patch, physical): (&RationalPatch, [bool; 4]),
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
    // The segment stays in one closed chart and knot span. Only a certified
    // atlas join may hold an interior witness on a chart's boundary.
    if !inside(&to, patch.domain(), physical, dec)? {
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
        if !inside(&at, &domain, [false; 4], dec)? || !inside(&to, &domain, [false; 4], dec)? {
            continue;
        }
        let parameters = std::array::from_fn(|axis| P::line(at[axis].clone(), to[axis].clone()));
        let homogeneous = carrier.span_polynomials(span, dec)?;
        let path = Path {
            numerators: [
                compose(&homogeneous[0], &parameters, dec)?,
                compose(&homogeneous[1], &parameters, dec)?,
                compose(&homogeneous[2], &parameters, dec)?,
            ],
            denominator: compose(&homogeneous[3], &parameters, dec)?,
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
            on_patch(model, (patch, [true; 4]), p, left, step, barriers, dec)?
        }
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => {
            let mut found = None;
            for chart in atlas.charts() {
                let physical = [
                    SplineBoundary::UMin,
                    SplineBoundary::UMax,
                    SplineBoundary::VMin,
                    SplineBoundary::VMax,
                ]
                .map(|side| chart.neighbor(side).is_none());
                if let Some(candidate) = on_patch(
                    model,
                    (chart.patch(), physical),
                    p,
                    left,
                    step,
                    barriers,
                    dec,
                )? {
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{BSplineSurface, SplineBoundary};
    use catallax_topology::{Edge, Sense};

    fn fixture(dec: &mut Decider) -> (Model, Surface, Point3, Vec<EdgeId>) {
        let half = Rational::from_ratio(1, 2).unwrap();
        let mut controls = Vec::new();
        for (x, z) in [(r(0), r(0)), (half.clone(), r(0)), (r(1), r(1))] {
            for y in [0, 1] {
                controls.push(Point3::new(
                    x.clone().into(),
                    Real::from_i64(y),
                    z.clone().into(),
                ));
            }
        }
        let patch = RationalPatch::new(
            BSplineSurface::new(
                [2, 1],
                [
                    vec![r(0), r(0), r(0), r(1), r(1), r(1)],
                    vec![r(0), r(0), r(1), r(1)],
                ],
                [3, 2],
                controls,
                vec![Real::one(); 6],
                dec,
            )
            .unwrap(),
            dec,
        )
        .unwrap();
        let p = patch.point_at(&[Real::zero(), half.into()], dec).unwrap();
        let mut model = Model::new();
        let edges = [
            SplineBoundary::UMin,
            SplineBoundary::UMax,
            SplineBoundary::VMin,
            SplineBoundary::VMax,
        ]
        .map(|side| {
            let curve = model
                .geom_mut()
                .add_curve(Curve::BSpline(patch.boundary(side)));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            })
        })
        .to_vec();
        (
            model,
            Surface::BSpline(SplineSupport::RationalPatch(patch)),
            p,
            edges,
        )
    }

    #[test]
    fn spline_probe_stays_in_the_finite_chart_and_leaves_its_native_boundary() {
        let mut dec = Decider::default();
        let (model, surface, p, mut edges) = fixture(&mut dec);
        let left = unit(0);
        let quarter = Real::from(Rational::from_ratio(1, 4).unwrap());
        let q = probe(&model, &surface, &p, &left, &quarter, &mut edges, &mut dec)
            .unwrap()
            .unwrap();
        assert_eq!(
            coefficient_sign(&(&q.x - &quarter), &mut dec).unwrap(),
            Sign::Zero
        );
        assert_eq!(
            coefficient_sign(&(&q.z - &quarter.square()), &mut dec).unwrap(),
            Sign::Zero
        );
        for step in [Real::one(), Real::from_i64(2), -quarter] {
            assert!(
                probe(&model, &surface, &p, &left, &step, &mut edges, &mut dec)
                    .unwrap()
                    .is_none()
            );
        }
    }

    #[test]
    fn spline_probe_rejects_an_intermediate_native_boundary_hit() {
        let mut dec = Decider::default();
        let (mut model, surface, p, mut edges) = fixture(&mut dec);
        let x = Real::from(Rational::from_ratio(1, 8).unwrap());
        let curve = model.geom_mut().add_curve(Curve::BSpline(
            BSplineCurve::new(
                1,
                vec![r(0), r(0), r(1), r(1)],
                vec![
                    Point3::new(x.clone(), Real::zero(), x.square()),
                    Point3::new(x.clone(), Real::one(), x.square()),
                ],
                vec![Real::one(); 2],
                &mut dec,
            )
            .unwrap(),
        ));
        edges.push(model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        }));
        for (denominator, clear) in [(16, true), (8, false), (4, false)] {
            let step = Real::from(Rational::from_ratio(1, denominator).unwrap());
            assert_eq!(
                probe(&model, &surface, &p, &unit(0), &step, &mut edges, &mut dec)
                    .unwrap()
                    .is_some(),
                clear
            );
        }
    }

    #[test]
    fn spline_probe_can_follow_a_certified_internal_chart_join() {
        use catallax_shapes::RationalAtlas;
        let mut dec = Decider::default();
        let (model, _, p, mut edges) = fixture(&mut dec);
        let half = Rational::from_ratio(1, 2).unwrap();
        let mut controls = Vec::new();
        for (x, z) in [(r(0), r(0)), (half.clone(), r(0)), (r(1), r(1))] {
            for y in [r(0), half.clone(), r(1)] {
                controls.push(Point3::new(x.clone().into(), y.into(), z.clone().into()));
            }
        }
        let source = catallax_shapes::BSplineSurface::new(
            [2, 1],
            [
                vec![r(0), r(0), r(0), r(1), r(1), r(1)],
                vec![r(0), r(0), half.clone(), r(1), r(1)],
            ],
            [3, 3],
            controls,
            vec![Real::one(); 9],
            &mut dec,
        )
        .unwrap();
        let patch = RationalPatch::new(source.clone(), &mut dec).unwrap();
        let atlas = RationalAtlas::partition(source, false, &mut dec)
            .unwrap()
            .remove(0);
        let quarter = Real::from(Rational::from_ratio(1, 4).unwrap());
        for support in [
            SplineSupport::RationalPatch(patch),
            SplineSupport::RationalAtlas(atlas),
        ] {
            let q = probe(
                &model,
                &Surface::BSpline(support),
                &p,
                &unit(0),
                &quarter,
                &mut edges,
                &mut dec,
            )
            .unwrap()
            .expect("an internal chart join is not a face boundary");
            for difference in [
                &q.x - &quarter,
                &q.y - &Real::from(half.clone()),
                &q.z - &quarter.square(),
            ] {
                assert_eq!(coefficient_sign(&difference, &mut dec).unwrap(), Sign::Zero);
            }
        }
    }
}
