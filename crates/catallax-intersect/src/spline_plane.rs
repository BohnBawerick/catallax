//! Complete transverse interior meetings of a finite spline with a plane.

use crate::spline_arithmetic::{arithmetic, cut, difference, quotient, value};
use crate::{Contact, CurveHit, CurveSurface};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{
    Curve, CurveId, GeomStore, Meet, MeetOf, Point3, PointDef, Position, Surface, SurfaceId,
};

fn missing(detail: &str) -> Refusal {
    Refusal::new("intersect.spline-plane", detail)
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    crate::spline_arithmetic::sign(value, dec, "intersect.spline-plane")
}

pub(super) fn hits(
    geom: &mut GeomStore,
    curve: CurveId,
    surface: SurfaceId,
    dec: &mut Decider,
) -> Result<Option<CurveSurface>, Refusal> {
    let (Curve::BSpline(spline), Surface::Plane(plane)) =
        (geom.curve(curve), geom.surface(surface))
    else {
        return Ok(None);
    };
    let (spline, plane) = (spline.clone(), plane.clone());
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
    if normal
        .iter()
        .chain(origin.iter())
        .any(|v| v.rational_ref().is_none())
    {
        return Err(missing(
            "finite spline/plane meetings need rational plane coefficients",
        ));
    }
    if normal.iter().all(|v| v.rational_ref().unwrap().is_zero()) {
        return Err(missing("a plane normal must be nonzero"));
    }
    let mut offset = Real::zero();
    for axis in 0..3 {
        let term = arithmetic(normal[axis], origin[axis], true, dec)?;
        offset = arithmetic(&offset, &term, false, dec)?;
    }
    let factors = [
        normal[0].clone(),
        normal[1].clone(),
        normal[2].clone(),
        -offset,
    ];
    let mut positions: Vec<Point3> = Vec::new();
    let (mut spans, mut contained) = (0, 0);
    for (span, [lo, hi]) in spline.spans() {
        spans += 1;
        let coordinates = spline.query_span_polynomials(span, dec)?;
        let equation = cut(&coordinates, &factors, dec)?;
        let rational = equation.rational().ok_or_else(|| {
            missing("finite spline/plane meetings need rational homogeneous cut coefficients")
        })?;
        if rational.is_zero() {
            contained += 1;
            continue;
        }
        for parameter in equation.real_roots(dec)? {
            let lower = sign(
                &difference(&parameter, &Real::from_rational(lo.clone()), dec)?,
                dec,
            )?;
            let upper = sign(
                &difference(&Real::from_rational(hi.clone()), &parameter, dec)?,
                dec,
            )?;
            if lower == Sign::Negative || upper == Sign::Negative {
                continue;
            }
            if lower == Sign::Zero || upper == Sign::Zero {
                return Err(missing(
                    "knot or endpoint meetings need a one-sided contact proof",
                ));
            }
            if sign(&value(&equation, &parameter, true, dec)?, dec)? == Sign::Zero {
                return Err(missing(
                    "finite spline/plane tangency needs a contact proof",
                ));
            }
            let mut homogeneous = Vec::new();
            for coordinate in &coordinates {
                homogeneous.push(value(coordinate, &parameter, false, dec)?);
            }
            // Positive spline weights prove W nonzero on the retained span.
            let point = Point3::new(
                quotient(&homogeneous[0], &homogeneous[3], dec)?,
                quotient(&homogeneous[1], &homogeneous[3], dec)?,
                quotient(&homogeneous[2], &homogeneous[3], dec)?,
            );
            for previous in &positions {
                let mut same = true;
                for (a, b) in [
                    (&point.x, &previous.x),
                    (&point.y, &previous.y),
                    (&point.z, &previous.z),
                ] {
                    if sign(&difference(a, b, dec)?, dec)? != Sign::Zero {
                        same = false;
                        break;
                    }
                }
                if same {
                    return Err(missing(
                        "multiple finite parameters at one point need a branch contact proof",
                    ));
                }
            }
            positions.push(point);
        }
    }
    if contained == spans {
        return Ok(Some(CurveSurface::Contained));
    }
    if contained != 0 {
        return Err(missing(
            "partial-span containment needs a finite overlap solver",
        ));
    }
    let hits = positions
        .into_iter()
        .map(|position| CurveHit {
            point: geom.add_point(PointDef::Meet(Meet {
                of: MeetOf::CurveSurface(curve, surface),
                position: Position::Closed(position),
            })),
            contact: Contact::Crossing,
        })
        .collect::<Vec<_>>();
    Ok(Some(if hits.is_empty() {
        CurveSurface::Disjoint
    } else {
        CurveSurface::Meets(hits)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Rational, with_operation_budget};
    use catallax_shapes::{BSplineCurve, Plane, Vector3};

    fn fixture(z: [i64; 3], weight: i64, dec: &mut Decider) -> (GeomStore, CurveId, SurfaceId) {
        let mut geom = GeomStore::new();
        let spline = BSplineCurve::new(
            2,
            [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
            z.into_iter()
                .enumerate()
                .map(|(i, z)| {
                    Point3::new(
                        Rational::from_ratio(i as i64, 2).unwrap().into(),
                        Real::zero(),
                        Real::from_i64(z),
                    )
                })
                .collect(),
            vec![Real::one(), Real::from_i64(weight), Real::one()],
            dec,
        )
        .unwrap();
        let curve = geom.add_curve(Curve::BSpline(spline));
        let surface = geom.add_surface(Surface::Plane(Plane {
            origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
            normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
            ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        }));
        (geom, curve, surface)
    }

    #[test]
    fn transverse_finite_spline_plane_roots_are_complete_with_positive_weights() {
        for weight in [1, 2] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let (mut geom, curve, surface) = fixture([1, -2, 1], weight, &mut dec);
                let answer =
                    crate::intersect_curve_surface(&mut geom, curve, surface, &mut dec).unwrap();
                assert!(answer.marks.is_empty());
                let CurveSurface::Meets(hits) = answer.value else {
                    panic!("expected two finite crossings")
                };
                assert_eq!(hits.len(), 2);
                for hit in hits {
                    assert!(matches!(hit.contact, Contact::Crossing));
                    assert_eq!(
                        sign(&geom.position(hit.point).z, &mut dec).unwrap(),
                        Sign::Zero
                    );
                }
            });
            assert!(!work.exhausted, "{} units", work.used);
        }
    }

    #[test]
    fn containment_pays_and_unproved_contacts_refuse_by_name() {
        let mut dec = Decider::default();
        let (mut geom, curve, surface) = fixture([0, 0, 0], 1, &mut dec);
        let (answer, work) = with_operation_budget(0, || hits(&mut geom, curve, surface, &mut dec));
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(work.exhausted);
        assert!(matches!(
            hits(&mut geom, curve, surface, &mut dec).unwrap(),
            Some(CurveSurface::Contained)
        ));
        for (z, reason) in [([1, -1, 1], "tangency"), ([0, 1, 1], "one-sided")] {
            let (mut geom, curve, surface) = fixture(z, 1, &mut dec);
            let error = hits(&mut geom, curve, surface, &mut dec).unwrap_err();
            assert_eq!(error.code, "intersect.spline-plane");
            assert!(error.detail.contains(reason), "{}", error.detail);
        }
    }
}
