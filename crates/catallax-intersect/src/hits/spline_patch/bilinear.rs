//! A constant-degree pullback for finite polynomial bilinear patches.

use super::*;
use crate::spline_arithmetic::{arithmetic, cut, derivative_value, difference, quotient, value};
use catallax_shapes::{BSplineCurve, RationalPatch};
use std::collections::{BTreeMap, btree_map::Entry};

const CODE: &str = "intersect.spline-patch";

fn exact_sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    crate::spline_arithmetic::sign(value, dec, CODE)
}

fn coordinates(point: &Point3) -> [&Real; 3] {
    [&point.x, &point.y, &point.z]
}

fn subtract(a: &[Real; 3], b: &[Real; 3], dec: &mut Decider) -> Result<[Real; 3], Refusal> {
    a.iter()
        .zip(b)
        .map(|(a, b)| difference(a, b, dec))
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| unreachable!())
}

fn cross(a: &[Real; 3], b: &[Real; 3], dec: &mut Decider) -> Result<[Real; 3], Refusal> {
    let mut result = Vec::new();
    for (i, j) in [(1, 2), (2, 0), (0, 1)] {
        let left = arithmetic(&a[i], &b[j], true, dec)?;
        let right = arithmetic(&a[j], &b[i], true, dec)?;
        result.push(difference(&left, &right, dec)?);
    }
    Ok(result.try_into().unwrap())
}

fn dot(a: &[Real; 3], b: &[Real; 3], dec: &mut Decider) -> Result<Real, Refusal> {
    let mut sum = Real::zero();
    for (a, b) in a.iter().zip(b) {
        let product = arithmetic(a, b, true, dec)?;
        sum = arithmetic(&sum, &product, false, dec)?;
    }
    Ok(sum)
}

struct Frame {
    factors: [[Real; 4]; 3],
    determinant: Real,
}

impl Frame {
    fn new(patch: &RationalPatch, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        dec.charge_work(1)?;
        let source = patch.carrier();
        if source.degrees() != [1, 1] || source.control_count() != [2, 2] {
            return Ok(None);
        }
        for (retained, original) in patch.domain().iter().zip(source.domain()) {
            dec.charge_work(1)?;
            if retained[0] != *original[0] || retained[1] != *original[1] {
                return Ok(None);
            }
        }
        let Some(weight) = source.weights()[0].rational_ref() else {
            return Ok(None);
        };
        for (point, w) in source.controls().iter().zip(source.weights()) {
            dec.charge_work(1)?;
            if w.rational_ref() != Some(weight)
                || coordinates(point)
                    .iter()
                    .any(|v| v.rational_ref().is_none())
            {
                return Ok(None);
            }
        }
        let poles: [[Real; 3]; 4] =
            std::array::from_fn(|i| coordinates(&source.controls()[i]).map(Clone::clone));
        let b = subtract(&poles[2], &poles[0], dec)?;
        let c = subtract(&poles[1], &poles[0], dec)?;
        let diagonal = subtract(&poles[3], &poles[2], dec)?;
        let d = subtract(&diagonal, &c, dec)?;
        let dual = [
            cross(&c, &d, dec)?,
            cross(&d, &b, dec)?,
            cross(&b, &c, dec)?,
        ];
        let determinant = dot(&b, &dual[0], dec)?;
        if exact_sign(&determinant, dec)? == Sign::Zero {
            return Ok(None);
        }
        let mut factors = Vec::new();
        for row in dual {
            let offset = dot(&row, &poles[0], dec)?;
            factors.push([row[0].clone(), row[1].clone(), row[2].clone(), -offset]);
        }
        Ok(Some(Self {
            factors: factors.try_into().unwrap(),
            determinant,
        }))
    }

    fn inside(
        &self,
        transformed: &[RealPoly; 3],
        weight: &RealPoly,
        at: &Real,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let denominator = arithmetic(
            &self.determinant,
            &value(weight, at, false, dec)?,
            true,
            dec,
        )?;
        for component in &transformed[..2] {
            let coordinate = quotient(&value(component, at, false, dec)?, &denominator, dec)?;
            if exact_sign(&coordinate, dec)? == Sign::Negative
                || exact_sign(&difference(&Real::one(), &coordinate, dec)?, dec)? == Sign::Negative
            {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn contains_control_hull(
        &self,
        curve: &BSplineCurve,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for point in curve.controls() {
            dec.charge_work(1)?;
            for factors in &self.factors[..2] {
                let mut coordinate = factors[3].clone();
                for (value, factor) in coordinates(point).into_iter().zip(factors) {
                    let term = arithmetic(value, factor, true, dec)?;
                    coordinate = arithmetic(&coordinate, &term, false, dec)?;
                }
                let coordinate = quotient(&coordinate, &self.determinant, dec)?;
                if exact_sign(&coordinate, dec)? == Sign::Negative
                    || exact_sign(&difference(&Real::one(), &coordinate, dec)?, dec)?
                        == Sign::Negative
                {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
}

fn product(a: &RealPoly, b: &RealPoly, dec: &mut Decider) -> Result<RealPoly, Refusal> {
    let mut terms = BTreeMap::new();
    for (ma, ca) in a.terms() {
        for (mb, cb) in b.terms() {
            dec.charge_work(1)?;
            let power = ma.exponents()[0] + mb.exponents()[0];
            let term = arithmetic(&ca, &cb, true, dec)?;
            match terms.entry(vec![power]) {
                Entry::Vacant(entry) => {
                    entry.insert(term);
                }
                Entry::Occupied(mut entry) => {
                    let sum = arithmetic(entry.get(), &term, false, dec)?;
                    entry.insert(sum);
                }
            }
        }
    }
    Ok(RealPoly::from_terms(1, terms.into_iter().collect()))
}

struct Root {
    parameter: Real,
    position: Point3,
    left: Option<Sign>,
    right: Option<Sign>,
    extension: Contact,
}

fn opposite(sign: Sign) -> Sign {
    match sign {
        Sign::Negative => Sign::Positive,
        Sign::Positive => Sign::Negative,
        Sign::Zero => unreachable!(),
    }
}

fn line_window(
    line: &Line,
    patch: &RationalPatch,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let direction = [line.dir.x.clone(), line.dir.y.clone(), line.dir.z.clone()];
    let origin = coordinates(&line.origin).map(Clone::clone);
    let squared = dot(&direction, &direction, dec)?;
    if exact_sign(&squared, dec)? == Sign::Zero {
        return Err(Refusal::new(
            "shapes.zero-vector",
            "a line direction must be nonzero",
        ));
    }
    let mut bounds: Option<[Rational; 2]> = None;
    for point in patch.carrier().controls() {
        dec.charge_work(1)?;
        let relative = subtract(&coordinates(point).map(Clone::clone), &origin, dec)?;
        let parameter = quotient(&dot(&relative, &direction, dec)?, &squared, dec)?;
        let interval = if let Some(value) = parameter.rational_ref() {
            [value.clone(), value.clone()]
        } else {
            let bits = dec.policy().max_precision_bits.min(128);
            if bits == 0 {
                return Err(Refusal::new(
                    CODE,
                    "no precision is available for a finite line window",
                ));
            }
            let value = dec.enclose_budgeted_at(&parameter, bits)?;
            [value.lo(), value.hi()]
        };
        bounds = Some(match bounds {
            None => interval,
            Some(old) => [
                old[0].clone().min(interval[0].clone()),
                old[1].clone().max(interval[1].clone()),
            ],
        });
    }
    let [lo, hi] = bounds.expect("a bilinear patch has four controls");
    if lo >= hi {
        return Err(Refusal::new(
            CODE,
            "the finite line window needs increasing bounds",
        ));
    }
    let mut controls = Vec::new();
    for parameter in [&lo, &hi] {
        let mut point = Vec::new();
        for (origin, direction) in origin.iter().zip(&direction) {
            let shift = arithmetic(direction, &Real::from(parameter.clone()), true, dec)?;
            point.push(arithmetic(origin, &shift, false, dec)?);
        }
        let [x, y, z]: [Real; 3] = point.try_into().unwrap();
        controls.push(Point3::new(x, y, z));
    }
    // Projection of a positive-weight surface hull bounds every meeting's
    // original line parameter. The exact line definition is retained.
    BSplineCurve::new(
        1,
        vec![lo.clone(), lo, hi.clone(), hi],
        controls,
        vec![Real::one(); 2],
        dec,
    )
}

pub(super) fn meetings(
    geom: &mut GeomStore,
    curve_id: CurveId,
    surface: SurfaceId,
    dec: &mut Decider,
) -> Result<Option<CurveSurface>, Refusal> {
    let Surface::BSpline(SplineSupport::RationalPatch(patch)) = geom.surface(surface) else {
        return Ok(None);
    };
    match geom.curve(curve_id) {
        Curve::BSpline(curve) if curve.degree() <= 2 => {
            for (point, weight) in curve.controls().zip(curve.weights()) {
                dec.charge_work(1)?;
                if weight.rational_ref().is_none()
                    || coordinates(point)
                        .iter()
                        .any(|v| v.rational_ref().is_none())
                {
                    return Ok(None);
                }
            }
        }
        Curve::Line(_) => {}
        _ => return Ok(None),
    }
    let Some(frame) = Frame::new(patch, dec)? else {
        return Ok(None);
    };
    let window = match geom.curve(curve_id) {
        Curve::Line(line) => Some(line_window(line, patch, dec)?),
        _ => None,
    };
    let is_line = window.is_some();
    let curve = if let Some(curve) = &window {
        curve
    } else {
        let Curve::BSpline(curve) = geom.curve(curve_id) else {
            unreachable!()
        };
        curve
    };
    let mut roots: Vec<Root> = Vec::new();
    let (mut spans, mut contained) = (0, 0);
    for (index, domain) in curve.spans() {
        dec.charge_work(1)?;
        spans += 1;
        let source = curve.query_span_polynomials(index, dec)?;
        let transformed = [
            cut(&source, &frame.factors[0], dec)?,
            cut(&source, &frame.factors[1], dec)?,
            cut(&source, &frame.factors[2], dec)?,
        ];
        // In the invertible basis (b,c,d), the patch is (u,v,u*v).
        // Positive finite weights introduce no denominator roots.
        let zw = product(&transformed[2], &source[3], dec)?;
        let xy = product(&transformed[0], &transformed[1], dec)?;
        let equation = cut(
            &[zw, xy, RealPoly::zero(1), RealPoly::zero(1)],
            &[
                frame.determinant.clone(),
                Real::from_i64(-1),
                Real::zero(),
                Real::zero(),
            ],
            dec,
        )?;
        if equation.is_zero(dec)? {
            contained += 1;
            continue;
        }
        let degree = equation.total_degree(dec)?.unwrap();
        for parameter in equation.real_roots(dec)? {
            dec.charge_work(1)?;
            let low = exact_sign(
                &difference(&parameter, &Real::from(domain[0].clone()), dec)?,
                dec,
            )?;
            let high = exact_sign(
                &difference(&Real::from(domain[1].clone()), &parameter, dec)?,
                dec,
            )?;
            if low == Sign::Negative
                || high == Sign::Negative
                || !frame.inside(&transformed, &source[3], &parameter, dec)?
            {
                continue;
            }
            let mut order = 1;
            let coefficient = loop {
                dec.charge_work(1)?;
                let coefficient =
                    exact_sign(&derivative_value(&equation, &parameter, order, dec)?, dec)?;
                if coefficient != Sign::Zero {
                    break coefficient;
                }
                order += 1;
                assert!(
                    order <= degree as usize,
                    "a nonzero polynomial has a nonzero leading derivative"
                );
            };
            let left = (low != Sign::Zero).then_some(if order % 2 == 0 {
                coefficient
            } else {
                opposite(coefficient)
            });
            let right = (high != Sign::Zero).then_some(coefficient);
            if let Some(previous) = roots.last_mut() {
                dec.charge_work(1)?;
                if exact_sign(&difference(&parameter, &previous.parameter, dec)?, dec)?
                    == Sign::Zero
                {
                    previous.left = previous.left.or(left);
                    previous.right = previous.right.or(right);
                    continue;
                }
            }
            roots.push(Root {
                position: curve.point_at(&parameter, dec)?,
                parameter,
                left,
                right,
                extension: if order % 2 == 0 {
                    Contact::Tangent
                } else {
                    Contact::Crossing
                },
            });
        }
    }
    if contained > 0 {
        if is_line {
            return Err(Refusal::new(
                "intersect.line-patch-overlap",
                "a line on a finite patch needs interval ownership",
            ));
        }
        if contained == spans && frame.contains_control_hull(curve, dec)? {
            return Ok(Some(CurveSurface::Contained));
        }
        return Err(Refusal::new(
            "intersect.spline-patch-overlap",
            "partial spline/patch containment needs finite interval ownership",
        ));
    }
    let mut hits = Vec::new();
    for root in roots {
        dec.charge_work(1)?;
        let contact = match (root.left, root.right) {
            (Some(a), Some(b)) if a == b => Contact::Tangent,
            (Some(_), Some(_)) => Contact::Crossing,
            _ => root.extension,
        };
        let point = geom.add_point(PointDef::Meet(Meet {
            of: MeetOf::CurveSurface(curve_id, surface),
            position: Position::Closed(root.position),
        }));
        hits.push(CurveHit { point, contact });
    }
    Ok(Some(if hits.is_empty() {
        CurveSurface::Disjoint
    } else {
        CurveSurface::Meets(hits)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    use catallax_shapes::BSplineSurface;

    fn scalar(value: &str) -> Real {
        Rational::parse(value).unwrap().into()
    }

    fn spline(points: [[&str; 3]; 3], dec: &mut Decider) -> BSplineCurve {
        BSplineCurve::new(
            2,
            [2, 2, 2, 5, 5, 5].map(Rational::from_i64).to_vec(),
            points
                .into_iter()
                .map(|p| {
                    let [x, y, z] = p.map(scalar);
                    Point3::new(x, y, z)
                })
                .collect(),
            vec![Real::one(); 3],
            dec,
        )
        .unwrap()
    }

    #[test]
    fn finite_bilinear_meetings_keep_touches_bounds_and_partial_overlap() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let source = BSplineSurface::new(
                [1, 1],
                std::array::from_fn(|_| [0, 0, 1, 1].map(Rational::from_i64).to_vec()),
                [2, 2],
                vec![
                    Point3::new(0.into(), 0.into(), 0.into()),
                    Point3::new(0.into(), 1.into(), 0.into()),
                    Point3::new(1.into(), 0.into(), 0.into()),
                    Point3::new(1.into(), 1.into(), 1.into()),
                ],
                vec![Real::one(); 4],
                &mut dec,
            )
            .unwrap();
            let patch = RationalPatch::new(source, &mut dec).unwrap();
            let mut geom = GeomStore::new();
            let surface = geom.add_surface(Surface::BSpline(patch.into()));
            for (points, contact) in [
                (
                    [
                        ["0", "1/2", "1/4"],
                        ["1/2", "1/2", "1/4"],
                        ["1", "1/2", "1/4"],
                    ],
                    Contact::Crossing,
                ),
                (
                    [
                        ["0", "1/2", "1/4"],
                        ["1/2", "1/2", "0"],
                        ["1", "1/2", "3/4"],
                    ],
                    Contact::Tangent,
                ),
            ] {
                let curve = geom.add_curve(Curve::BSpline(spline(points, &mut dec)));
                assert!(!super::super::disjoint(&geom, curve, surface, &mut dec).unwrap());
                let Some(CurveSurface::Meets(hits)) =
                    meetings(&mut geom, curve, surface, &mut dec).unwrap()
                else {
                    panic!("one finite meeting")
                };
                assert_eq!(hits.len(), 1);
                assert_eq!(hits[0].contact, contact);
                let point = geom.position(hits[0].point);
                assert_eq!(
                    exact_sign(
                        &point
                            .minus(&Point3::new(scalar("1/2"), scalar("1/2"), scalar("1/4")))
                            .norm_squared(),
                        &mut dec
                    )
                    .unwrap(),
                    Sign::Zero
                );
            }
            let finite = geom.add_curve(Curve::BSpline(spline(
                [
                    ["0", "1/2", "1/4"],
                    ["1/8", "1/2", "1/4"],
                    ["1/4", "1/2", "1/4"],
                ],
                &mut dec,
            )));
            assert!(matches!(
                meetings(&mut geom, finite, surface, &mut dec).unwrap(),
                Some(CurveSurface::Disjoint)
            ));
            let far = geom.add_curve(Curve::BSpline(spline(
                [
                    ["2", "1/2", "1/4"],
                    ["5/2", "1/2", "1/4"],
                    ["3", "1/2", "1/4"],
                ],
                &mut dec,
            )));
            assert!(super::super::disjoint(&geom, far, surface, &mut dec).unwrap());
            for (last, contained) in [("1", true), ("2", false)] {
                let curve = geom.add_curve(Curve::BSpline(spline(
                    [["0", "0", "0"], ["1/2", "0", "0"], [last, "0", "0"]],
                    &mut dec,
                )));
                assert!(!super::super::disjoint(&geom, curve, surface, &mut dec).unwrap());
                let answer = meetings(&mut geom, curve, surface, &mut dec);
                if contained {
                    assert!(matches!(answer.unwrap(), Some(CurveSurface::Contained)));
                } else {
                    assert_eq!(answer.unwrap_err().code, "intersect.spline-patch-overlap");
                }
            }
            for (z, direction, expected) in [
                ("0", [0, 0, 1], Contact::Crossing),
                ("1/4", [1, -1, 0], Contact::Tangent),
            ] {
                let [x, y, z_direction] = direction.map(Real::from_i64);
                let curve = geom.add_curve(Curve::Line(Line {
                    origin: Point3::new(scalar("1/2"), scalar("1/2"), scalar(z)),
                    dir: Vector3::new(x, y, z_direction),
                }));
                let Some(CurveSurface::Meets(hits)) =
                    meetings(&mut geom, curve, surface, &mut dec).unwrap()
                else {
                    panic!("one finite line meeting")
                };
                assert_eq!(hits.len(), 1);
                assert_eq!(hits[0].contact, expected);
                let point = geom.position(hits[0].point);
                assert_eq!(
                    exact_sign(
                        &point
                            .minus(&Point3::new(scalar("1/2"), scalar("1/2"), scalar("1/4")))
                            .norm_squared(),
                        &mut dec
                    )
                    .unwrap(),
                    Sign::Zero
                );
            }
            let ruling = geom.add_curve(Curve::Line(Line {
                origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
                dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
            }));
            assert_eq!(
                meetings(&mut geom, ruling, surface, &mut dec)
                    .unwrap_err()
                    .code,
                "intersect.line-patch-overlap"
            );
        });
        assert!(!work.exhausted, "{work:?}");
    }
}
