//! Curve meetings needed when a non-maximal operand has coplanar faces.
//! Rational elimination remains the fast path; algebraic coefficients use
//! shared certified operations and every candidate is checked before incidence.

use catallax_intersect::{
    Component, CurveSurface, SurfacePair, intersect_curve_surface, intersect_surfaces, on_curve,
    same_curve,
};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, MPoly, Outcome, Real, RealPoly, Refusal, Sign,
};
use catallax_shapes::{
    Curve, CurveId, GeomKey, Meet, MeetOf, Point3, PointDef, PointId, Position, Vector3,
};
use catallax_topology::Model;

use crate::geometry::{canonical_point, coincident, refusal, sign};

fn retain_incidence(
    model: &mut Model,
    curve: CurveId,
    point: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    match on_curve(model.geom(), curve, point, dec) {
        Outcome::ProvenCoincident(proof) => {
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Point(point),
                b: GeomKey::Curve(curve),
                at: None,
                origin: FactOrigin::Proof(proof.proof),
            });
            Ok(true)
        }
        other => coincident(other, "boolean.imprint-ambiguous"),
    }
}

fn substitute(poly: &RealPoly, origin: &Point3, u: &Vector3, v: Option<&Vector3>) -> RealPoly {
    let n = if v.is_some() { 2 } else { 1 };
    let xs = [&origin.x, &origin.y, &origin.z];
    let us = [&u.x, &u.y, &u.z];
    let vs = v.map(|v| [&v.x, &v.y, &v.z]);
    let coordinates: Vec<_> = (0..3)
        .map(|i| {
            let mut p = RealPoly::constant(n, xs[i].clone()) + RealPoly::var(n, 0).scale(us[i]);
            if let Some(vs) = vs {
                p = p + RealPoly::var(n, 1).scale(vs[i]);
            }
            p
        })
        .collect();
    let mut result = RealPoly::zero(n);
    for (monomial, coefficient) in poly.terms() {
        let mut term = RealPoly::constant(n, coefficient);
        for (axis, &power) in monomial.exponents().iter().enumerate() {
            for _ in 0..power {
                term = &term * &coordinates[axis];
            }
        }
        result = result + term;
    }
    result
}

fn univariate(poly: &RealPoly, axis: usize) -> RealPoly {
    let mut result = Vec::new();
    for (monomial, coefficient) in poly.terms() {
        assert!(
            monomial
                .exponents()
                .iter()
                .enumerate()
                .all(|(i, &p)| i == axis || p == 0)
        );
        let power = monomial.exponents()[axis];
        result.push((vec![power], coefficient));
    }
    RealPoly::from_terms(1, result)
}

fn candidate(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    p: Point3,
    dec: &mut Decider,
) -> Result<Option<PointId>, Refusal> {
    for curve in [a, b] {
        for equation in model.geom().curve(curve).equations(model.geom()).ok_or_else(|| catallax_numbers::Refusal::new("boolean.parametric-query", "curve equation unavailable; this boolean path needs a finite parametric solver"))? {
            if sign(
                &equation.eval_real(&[p.x.clone(), p.y.clone(), p.z.clone()]),
                dec,
                "boolean.imprint-ambiguous",
            )? != Sign::Zero
            {
                return Ok(None);
            }
        }
    }
    // Only explicit analytic carriers reach this path. Their common zero
    // set is the curve itself, so the preceding exact equation proofs also
    // certify membership. Intersection branches are handled separately.
    Ok(Some(model.geom_mut().add_point(PointDef::Meet(Meet {
        of: MeetOf::Curves(a, b),
        position: Position::Closed(p),
    }))))
}

fn roots(poly: &RealPoly, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
    if poly.is_zero(dec)? {
        return Ok(Vec::new());
    }
    poly.real_roots(dec)
}

fn affine(origin: &Point3, u: &Vector3, v: Option<&Vector3>, x: &Real, y: Option<&Real>) -> Point3 {
    let p = origin.plus(&u.scaled(x));
    match (v, y) {
        (Some(v), Some(y)) => p.plus(&v.scaled(y)),
        _ => p,
    }
}

fn quadratic(poly: &RealPoly, x: &Real, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
    let mut c = [Real::zero(), Real::zero(), Real::zero()];
    for (m, k) in poly.terms() {
        let e = m.exponents();
        assert!(e[1] <= 2);
        c[e[1] as usize] = &c[e[1] as usize] + &(k * x.pow_u32(e[0]));
    }
    if sign(&c[2], dec, "boolean.imprint-ambiguous")? == Sign::Zero {
        if sign(&c[1], dec, "boolean.imprint-ambiguous")? == Sign::Zero {
            return Ok(Vec::new());
        }
        return Ok(vec![-(&c[0] / &c[1])]);
    }
    let discriminant = c[1].square() - (&c[2] * &c[0]) * Real::from_i64(4);
    let denominator = &Real::from_i64(2) * &c[2];
    match sign(&discriminant, dec, "boolean.imprint-ambiguous")? {
        Sign::Negative => Ok(Vec::new()),
        Sign::Zero => Ok(vec![-&c[1] / denominator]),
        Sign::Positive => Ok(vec![
            &(-&c[1] - discriminant.sqrt()) / &denominator,
            (-&c[1] + discriminant.sqrt()) / denominator,
        ]),
    }
}

pub(crate) enum Meetings {
    Coincident,
    Points(Vec<PointId>),
}

pub(crate) fn meetings(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Meetings, Refusal> {
    if coincident(
        same_curve(model.geom(), a, b, dec),
        "boolean.imprint-ambiguous",
    )? {
        return Ok(Meetings::Coincident);
    }
    let mut result = Vec::new();
    if let Curve::Intersection(branch) = model.geom().curve(b).clone() {
        for surface in branch.surfaces {
            let answer = intersect_curve_surface(model.geom_mut(), a, surface, dec)?;
            if !answer.marks.is_empty() {
                return Err(refusal(
                    "boolean.imprint-ambiguous",
                    "curve pair elimination required unresolved incidence",
                ));
            }
            let finite = !matches!(&answer.value, CurveSurface::Contained);
            if let CurveSurface::Meets(hits) = answer.value {
                for hit in hits {
                    if retain_incidence(model, b, hit.point, dec)? {
                        canonical_point(
                            model,
                            &mut result,
                            hit.point,
                            dec,
                            "boolean.imprint-ambiguous",
                        )?;
                    }
                }
            }
            // b lies on both supports. Any finite meeting with one support
            // exhausts a's possible meetings with b; the branch predicate
            // above still checks the other support and the component.
            if finite {
                break;
            }
        }
        let pair = intersect_surfaces(
            model.geom_mut(),
            branch.surfaces[0],
            branch.surfaces[1],
            dec,
        )?;
        if let SurfacePair::Meets(parts) = pair.value {
            for part in parts {
                let points = match part {
                    Component::Curve { singular, .. } => singular,
                    Component::Point { point } => vec![point],
                };
                for point in points {
                    if retain_incidence(model, a, point, dec)?
                        && retain_incidence(model, b, point, dec)?
                    {
                        canonical_point(
                            model,
                            &mut result,
                            point,
                            dec,
                            "boolean.imprint-ambiguous",
                        )?;
                    }
                }
            }
        }
        return Ok(Meetings::Points(result));
    }
    if matches!(model.geom().curve(a), Curve::Intersection(_)) {
        return meetings(model, b, a, dec);
    }
    if let Curve::Line(line) = model.geom().curve(a).clone() {
        let mut common: Option<RealPoly> = None;
        let mut common_proven = true;
        for equation in model.geom().curve(b).equations(model.geom()).ok_or_else(|| catallax_numbers::Refusal::new("boolean.parametric-query", "curve equation unavailable; this boolean path needs a finite parametric solver"))? {
            let p = univariate(&substitute(&equation, &line.origin, &line.dir, None), 0);
            if !p.is_zero(dec)? {
                common = Some(match common {
                    None => p,
                    Some(q) => {
                        if let (Some(a), Some(b)) = (q.rational(), p.rational()) {
                            let gcd = a.to_upoly(0).unwrap().gcd(&b.to_upoly(0).unwrap());
                            RealPoly::from(MPoly::from_terms(
                                1,
                                gcd.coefficients()
                                    .iter()
                                    .enumerate()
                                    .map(|(i, c)| (vec![i as u32], c.clone()))
                                    .collect(),
                            ))
                        } else {
                            common_proven = false;
                            q
                        }
                    }
                });
            }
        }
        if let Some(poly) = common {
            for t in roots(&poly, dec)? {
                let position = affine(&line.origin, &line.dir, None, &t, None);
                if common_proven {
                    // Every substituted equation is divisible by the retained
                    // GCD (or is the sole nonzero equation). Its certified root
                    // and the affine line definition prove both incidences.
                    result.push(model.geom_mut().add_point(PointDef::Meet(Meet {
                        of: MeetOf::Curves(a, b),
                        position: Position::Closed(position),
                    })));
                } else if let Some(point) = candidate(model, a, b, position, dec)? {
                    // A non-rational equation was not reduced into the GCD.
                    // Its membership must still be proved by shared predicates.
                    result.push(point);
                }
            }
        }
        return Ok(Meetings::Points(result));
    }
    if matches!(model.geom().curve(b), Curve::Line(_)) {
        return meetings(model, b, a, dec);
    }
    let (origin, u, normal) = match model.geom().curve(a) {
        Curve::Circle(c) => (c.center.clone(), c.ref_dir.clone(), c.normal.clone()),
        Curve::Ellipse(c) => (c.center.clone(), c.major_dir.clone(), c.normal.clone()),
        _ => {
            return Err(refusal(
                "boolean.imprint-ambiguous",
                "curve does not provide a supported elimination chart",
            ));
        }
    };
    let v = normal.cross(&u);
    let mut first = None;
    for equation in model
        .geom()
        .curve(a)
        .equations(model.geom())
        .ok_or_else(|| {
            catallax_numbers::Refusal::new(
                "boolean.parametric-query",
                "curve equation unavailable; this boolean path needs a finite parametric solver",
            )
        })?
    {
        let polynomial = substitute(&equation, &origin, &u, Some(&v));
        if !polynomial.is_zero(dec)? {
            first = Some(polynomial);
            break;
        }
    }
    let first = first.ok_or_else(|| {
        refusal(
            "boolean.imprint-ambiguous",
            "analytic conic equation was not certified",
        )
    })?;
    for equation in model
        .geom()
        .curve(b)
        .equations(model.geom())
        .ok_or_else(|| {
            catallax_numbers::Refusal::new(
                "boolean.parametric-query",
                "curve equation unavailable; this boolean path needs a finite parametric solver",
            )
        })?
    {
        let second = substitute(&equation, &origin, &u, Some(&v));
        if second.is_zero(dec)? {
            continue;
        }
        let eliminated = first.resultant(&second, 1, dec)?;
        if eliminated.is_zero(dec)? {
            continue;
        }
        for x in roots(&univariate(&eliminated, 0), dec)? {
            let mut ys = quadratic(&first, &x, dec)?;
            if ys.is_empty() {
                ys = quadratic(&second, &x, dec)?;
            }
            for y in ys {
                if let Some(point) = candidate(
                    model,
                    a,
                    b,
                    affine(&origin, &u, Some(&v), &x, Some(&y)),
                    dec,
                )? {
                    canonical_point(model, &mut result, point, dec, "boolean.imprint-ambiguous")?;
                }
            }
        }
        return Ok(Meetings::Points(result));
    }
    // Distinct coplanar analytic conics cannot share a one-dimensional
    // component. A missing eliminant therefore needs a proof, not a guess.
    Err(refusal(
        "boolean.imprint-ambiguous",
        "curve equations did not separate distinct carriers",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Rational, StrictSign};
    use catallax_shapes::{
        Circle, Cylinder, Ellipse, IntersectionCurve, Line, Plane, RPoint3, RVector3, Surface,
    };

    fn r(n: i64) -> Rational {
        Rational::from_i64(n)
    }

    fn points(model: &mut Model, a: CurveId, b: CurveId) -> Vec<Point3> {
        let Meetings::Points(points) = meetings(model, a, b, &mut Decider::default()).unwrap()
        else {
            panic!("distinct carriers cannot coincide");
        };
        points
            .into_iter()
            .map(|p| model.geom().position(p))
            .collect()
    }

    #[test]
    fn line_conic_gcd_keeps_every_carrier_constraint() {
        for reverse in [false, true] {
            let mut model = Model::new();
            let mut dec = Decider::default();
            let ellipse = model.geom_mut().add_curve(Curve::Ellipse(
                Ellipse::new(
                    RPoint3::new(r(1), r(-2), r(3)).into(),
                    RVector3::new(r(0), r(0), r(1)).into(),
                    RVector3::new(r(1), r(0), r(0)).into(),
                    r(5).into(),
                    r(2).into(),
                    &mut dec,
                )
                .unwrap(),
            ));
            for z in [3, 4] {
                let line = model.geom_mut().add_curve(Curve::Line(
                    Line::new(
                        RPoint3::new(r(0), r(-1), r(z)).into(),
                        RVector3::new(r(1), r(0), r(0)).into(),
                        &mut dec,
                    )
                    .unwrap(),
                ));
                let found = if reverse {
                    points(&mut model, ellipse, line)
                } else {
                    points(&mut model, line, ellipse)
                };
                if z == 4 {
                    assert!(found.is_empty(), "radial roots miss the ellipse plane");
                    continue;
                }
                assert_eq!(found.len(), 2);
                let mut dec = Decider::default();
                let mut signs = Vec::new();
                for p in found {
                    let x = p.x - Real::one();
                    assert!(
                        dec.sign(&(x.square() * Real::from_i64(2) - Real::from_i64(5)))
                            .is_coincident()
                    );
                    assert!(dec.sign(&(p.y + Real::one())).is_coincident());
                    assert!(dec.sign(&(p.z - Real::from_i64(3))).is_coincident());
                    signs.push(*dec.sign(&x).proven().unwrap());
                }
                signs.sort();
                assert_eq!(signs, [StrictSign::Negative, StrictSign::Positive]);
            }
        }
    }

    #[test]
    fn algebraic_line_candidates_keep_every_carrier_equation() {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let z = Vector3::new(Real::zero(), Real::zero(), Real::one());
        let circle = model.geom_mut().add_curve(Curve::Circle(
            Circle::new(
                Point3::new(Real::zero(), Real::zero(), Real::from_i64(2).sqrt()),
                z.clone(),
                Real::one(),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap(),
        ));
        let line = model.geom_mut().add_curve(Curve::Line(
            Line::new(
                Point3::new(Real::zero(), Real::zero(), Real::zero()),
                z,
                &mut dec,
            )
            .unwrap(),
        ));
        assert!(points(&mut model, line, circle).is_empty());
    }

    #[test]
    fn conic_candidates_reject_extraneous_resultant_roots() {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let mut circles = Vec::new();
        for x in [0, 1] {
            circles.push(
                model.geom_mut().add_curve(Curve::Circle(
                    Circle::new(
                        RPoint3::new(r(x), r(x), r(0)).into(),
                        RVector3::new(r(0), r(0), r(1)).into(),
                        r(1).into(),
                        RVector3::new(r(1), r(0), r(0)).into(),
                        &mut dec,
                    )
                    .unwrap(),
                )),
            );
        }
        let found = points(&mut model, circles[0], circles[1]);
        assert_eq!(found.len(), 2);
        let mut x_zeroes = 0;
        for p in found {
            assert!(dec.sign(&(&p.x + &p.y - Real::one())).is_coincident());
            assert!(dec.sign(&(&p.x * &p.y)).is_coincident());
            assert!(dec.sign(&p.z).is_coincident());
            x_zeroes += usize::from(dec.sign(&p.x).is_coincident());
        }
        assert_eq!(x_zeroes, 1);
    }

    #[test]
    fn finite_support_meetings_still_check_the_other_support_and_branch() {
        for reverse in [false, true] {
            let mut model = Model::new();
            let mut dec = Decider::default();
            let origin: Point3 = RPoint3::new(r(0), r(0), r(0)).into();
            let x: Vector3 = RVector3::new(r(1), r(0), r(0)).into();
            let z: Vector3 = RVector3::new(r(0), r(0), r(1)).into();
            let cylinder = model.geom_mut().add_surface(Surface::Cylinder(
                Cylinder::new(origin.clone(), z.clone(), r(1).into(), x.clone(), &mut dec).unwrap(),
            ));
            let plane = model.geom_mut().add_surface(Surface::Plane(
                Plane::new(origin.clone(), x.clone(), z.clone(), &mut dec).unwrap(),
            ));
            let witness = model
                .geom_mut()
                .add_point(PointDef::Explicit(RPoint3::new(r(0), r(1), r(0)).into()));
            let surfaces = if reverse {
                [plane, cylinder]
            } else {
                [cylinder, plane]
            };
            let branch = model
                .geom_mut()
                .add_curve(Curve::Intersection(IntersectionCurve { surfaces, witness }));
            for (center, radius_sq, expected) in [
                (origin, r(1), 1),
                (RPoint3::new(r(1), r(0), r(0)).into(), r(4), 0),
            ] {
                let circle = model.geom_mut().add_curve(Curve::Circle(
                    Circle::new(center, z.clone(), radius_sq.into(), x.clone(), &mut dec).unwrap(),
                ));
                let found = points(&mut model, circle, branch);
                assert_eq!(found.len(), expected);
                for p in found {
                    let mut dec = Decider::default();
                    assert!(dec.sign(&p.x).is_coincident());
                    assert!(dec.sign(&(p.y - Real::one())).is_coincident());
                    assert!(dec.sign(&p.z).is_coincident());
                }
            }
        }
    }
}
