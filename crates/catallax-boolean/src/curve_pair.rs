//! Curve meetings needed when a non-maximal operand has coplanar faces.
//! Elimination is over rational polynomials; candidate coordinates remain
//! exact real expressions, and every candidate is checked before incidence.

use catallax_intersect::{
    Component, CurveSurface, SurfacePair, intersect_curve_surface, intersect_surfaces, on_curve,
    same_curve,
};
use catallax_numbers::{Decider, MPoly, Outcome, Rational, Real, Refusal, Sign, UPoly};
use catallax_shapes::{
    Curve, CurveId, Meet, MeetOf, Point3, PointDef, PointId, Position, RPoint3, RVector3,
};
use catallax_topology::Model;

use crate::geometry::{canonical_point, coincident, refusal, sign};

fn substitute(poly: &MPoly, origin: &RPoint3, u: &RVector3, v: Option<&RVector3>) -> MPoly {
    let n = if v.is_some() { 2 } else { 1 };
    let xs = [&origin.x, &origin.y, &origin.z];
    let us = [&u.x, &u.y, &u.z];
    let vs = v.map(|v| [&v.x, &v.y, &v.z]);
    let coordinates: Vec<_> = (0..3)
        .map(|i| {
            let mut p = MPoly::constant(n, xs[i].clone()) + MPoly::var(n, 0).scale(us[i]);
            if let Some(vs) = vs {
                p = p + MPoly::var(n, 1).scale(vs[i]);
            }
            p
        })
        .collect();
    let mut result = MPoly::zero(n);
    for (monomial, coefficient) in poly.terms() {
        let mut term = MPoly::constant(n, coefficient);
        for (axis, &power) in monomial.exponents().iter().enumerate() {
            for _ in 0..power {
                term = &term * &coordinates[axis];
            }
        }
        result = result + term;
    }
    result
}

fn univariate(poly: &MPoly, axis: usize) -> UPoly {
    let mut result = Vec::new();
    for (monomial, coefficient) in poly.terms() {
        assert!(
            monomial
                .exponents()
                .iter()
                .enumerate()
                .all(|(i, &p)| i == axis || p == 0)
        );
        let power = monomial.exponents()[axis] as usize;
        result.resize(result.len().max(power + 1), Rational::zero());
        result[power] = &result[power] + &coefficient;
    }
    UPoly::new(result)
}

fn candidate(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    p: Point3,
    dec: &mut Decider,
) -> Result<Option<PointId>, Refusal> {
    for curve in [a, b] {
        for equation in model.geom().curve(curve).equations(model.geom()) {
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
    let of = match model.geom().curve(a) {
        Curve::Intersection(c) => MeetOf::Surfaces(c.surfaces.to_vec()),
        _ => MeetOf::Curves(a, a),
    };
    let point = model.geom_mut().add_point(PointDef::Meet(Meet {
        of,
        position: Position::Closed(p),
    }));
    if !coincident(
        on_curve(model.geom(), a, point, dec),
        "boolean.imprint-ambiguous",
    )? || !coincident(
        on_curve(model.geom(), b, point, dec),
        "boolean.imprint-ambiguous",
    )? {
        return Ok(None);
    }
    Ok(Some(point))
}

fn roots(poly: &UPoly, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
    if poly.is_zero() {
        return Ok(Vec::new());
    }
    match dec.real_roots(poly) {
        Outcome::Proven(roots) => Ok(roots
            .into_iter()
            .map(|root| Real::from_algebraic(root.value))
            .collect()),
        Outcome::NotExact(n) => Err(refusal("boolean.imprint-ambiguous", n.detail)),
        Outcome::Refused(e) => Err(e),
        Outcome::ProvenCoincident(_) => unreachable!("nonzero polynomial roots are strict answers"),
    }
}

fn affine(
    origin: &RPoint3,
    u: &RVector3,
    v: Option<&RVector3>,
    x: &Real,
    y: Option<&Real>,
) -> Point3 {
    let p = origin.to_real().plus(&u.to_real().scaled(x));
    match (v, y) {
        (Some(v), Some(y)) => p.plus(&v.to_real().scaled(y)),
        _ => p,
    }
}

fn quadratic(poly: &MPoly, x: &Real, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
    let mut c = [Real::zero(), Real::zero(), Real::zero()];
    for (m, k) in poly.terms() {
        let e = m.exponents();
        assert!(e[1] <= 2);
        c[e[1] as usize] = &c[e[1] as usize] + &(Real::from(k) * x.pow_u32(e[0]));
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
            if let CurveSurface::Meets(hits) = answer.value {
                for hit in hits {
                    if coincident(
                        on_curve(model.geom(), b, hit.point, dec),
                        "boolean.imprint-ambiguous",
                    )? {
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
                    if coincident(
                        on_curve(model.geom(), a, point, dec),
                        "boolean.imprint-ambiguous",
                    )? && coincident(
                        on_curve(model.geom(), b, point, dec),
                        "boolean.imprint-ambiguous",
                    )? {
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
        let mut common: Option<UPoly> = None;
        for equation in model.geom().curve(b).equations(model.geom()) {
            let p = univariate(&substitute(&equation, &line.origin, &line.dir, None), 0);
            if !p.is_zero() {
                common = Some(common.map_or(p.clone(), |q| q.gcd(&p)));
            }
        }
        if let Some(poly) = common {
            for t in roots(&poly, dec)? {
                if let Some(point) = candidate(
                    model,
                    a,
                    b,
                    affine(&line.origin, &line.dir, None, &t, None),
                    dec,
                )? {
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
    let first = model
        .geom()
        .curve(a)
        .equations(model.geom())
        .into_iter()
        .map(|p| substitute(&p, &origin, &u, Some(&v)))
        .find(|p| !p.is_zero())
        .expect("analytic conic equation");
    for equation in model.geom().curve(b).equations(model.geom()) {
        let second = substitute(&equation, &origin, &u, Some(&v));
        if second.is_zero() {
            continue;
        }
        let eliminated = first.resultant(&second, 1);
        if eliminated.is_zero() {
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
