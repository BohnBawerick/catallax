use crate::{
    CurveForm, Cyclic,
    chart::{self, Chart},
    math::*,
    pair::{self, Frame, Plan},
};
use catallax_numbers::{
    Ball, Decider, FactKind, Outcome, Rational, Real, Refusal, Sign, StrictOrder, StrictSign,
};
use catallax_shapes::{
    Apart, BSplineCurve, Curve, CurveId, GeomKey, GeomStore, Point3, PointId, SplineCurveParameter,
    Surface, Vector3,
};
use std::sync::Arc;

pub(super) fn validate(
    geom: &GeomStore,
    c: CurveId,
    dec: &mut Decider,
) -> Result<Option<Arc<Plan>>, Refusal> {
    match geom.curve(c) {
        Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) | Curve::BSpline(_) => Ok(None),
        Curve::Intersection(i) => {
            if !crate::supports(
                geom.surface(i.surfaces[0]).kind(),
                geom.surface(i.surfaces[1]).kind(),
            ) {
                return Err(refused_curve());
            }
            let p = pair::plan(geom, i.surfaces, dec)?;
            if let (Surface::Cylinder(a), Surface::Cylinder(b)) =
                (geom.surface(i.surfaces[0]), geom.surface(i.surfaces[1]))
                && p.frame.is_none()
                && sign(&a.axis.cross(&b.axis).norm_squared(), dec)? != Sign::Zero
                && sign(&(&a.radius_sq - &b.radius_sq), dec)? == Sign::Zero
            {
                return Err(refused_curve());
            }
            Ok(Some(p))
        }
        _ => Err(refused_curve()),
    }
}
pub(super) fn form(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<CurveForm, Refusal> {
    let p = validate(geom, c, dec)?;
    Ok(match geom.curve(c) {
        Curve::Line(l) => CurveForm::Line {
            origin: l.origin.clone(),
            dir: l.dir.clone(),
        },
        Curve::Circle(l) => CurveForm::Circle {
            center: l.center.clone(),
            normal: l.normal.clone(),
            radius_sq: real(&l.radius_sq),
            ref_dir: l.ref_dir.clone(),
        },
        Curve::Ellipse(l) => CurveForm::Ellipse {
            center: l.center.clone(),
            normal: l.normal.clone(),
            major_dir: l.major_dir.clone(),
            major_sq: real(&l.major_sq),
            minor_sq: real(&l.minor_sq),
        },
        Curve::BSpline(curve) => CurveForm::BSpline(curve.clone()),
        Curve::Intersection(i) => {
            let plan = p.unwrap();
            if plan.frame.is_some() {
                return Ok(CurveForm::Branch);
            }
            let point = plan.divider.as_ref().map(|_| geom.position(i.witness));
            let index = if let Some(d) = &plan.divider {
                if sign(&d.value(point.as_ref().unwrap()), dec)? == Sign::Negative {
                    1
                } else {
                    0
                }
            } else {
                0
            };
            let mut form = plan.curves[index].form.clone();
            if let CurveForm::Line { origin, .. } = &mut form {
                *origin = point.unwrap_or_else(|| geom.position(i.witness));
            }
            form
        }
        _ => return Err(refused_curve()),
    })
}
// Comparing point sets does not require choosing a natural orientation.
fn carrier(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<CurveForm, Refusal> {
    let Some(plan) = validate(geom, c, dec)? else {
        return form(geom, c, dec);
    };
    if plan.frame.is_some() {
        return Ok(CurveForm::Branch);
    }
    if plan
        .divider
        .as_ref()
        .is_some_and(|d| !d.canonical_direction)
    {
        return form(geom, c, dec);
    }
    let Curve::Intersection(i) = geom.curve(c) else {
        unreachable!()
    };
    let witness = geom.position(i.witness);
    let mut shape = plan.curves[0].form.clone();
    match &mut shape {
        CurveForm::Line { origin, .. } => *origin = witness,
        CurveForm::Circle { center, normal, .. } if plan.divider.is_some() => {
            let sphere = i
                .surfaces
                .iter()
                .find_map(|s| {
                    if let Surface::Sphere(s) = geom.surface(*s) {
                        Some(s)
                    } else {
                        None
                    }
                })
                .expect("CS1 sphere");
            let base = sphere.center.clone();
            let z = witness.minus(&base).dot(&normal.clone()) / real(&normal.norm_squared());
            *center = base.plus(&normal.clone().scaled(&z));
        }
        _ => {}
    }
    Ok(shape)
}
fn orientation(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<Real, Refusal> {
    let Some(plan) = validate(geom, c, dec)? else {
        return Ok(Real::one());
    };
    let Some(divider) = &plan.divider else {
        return Ok(Real::one());
    };
    if !divider.canonical_direction {
        return Ok(Real::one());
    }
    let Curve::Intersection(i) = geom.curve(c) else {
        unreachable!()
    };
    let value = divider.value(&geom.position(i.witness));
    Ok(&value / &value.square().sqrt())
}
pub(super) fn surface_domain(
    surface: &Surface,
    point: &Point3,
    dec: &mut Decider,
) -> Result<bool, Failure> {
    let Surface::Cone(cone) = surface else {
        return Ok(true);
    };
    let normalized;
    let (axis, lower) = if let Some([x, y, z]) = cone.axis.rational_components() {
        (&cone.axis, x.abs().max(y.abs()).max(z.abs()))
    } else {
        normalized = unit(&cone.axis);
        (&normalized, r(1))
    };
    match sign(&point.minus(&cone.apex).dot(axis), dec) {
        Ok(value) => Ok(value != Sign::Negative),
        Err(Failure::Uncertain(mut n)) => {
            n.bound = n.bound / lower;
            Err(Failure::Uncertain(n))
        }
        Err(e) => Err(e),
    }
}
pub(super) fn curve_domain(
    geom: &GeomStore,
    curve: CurveId,
    point: &Point3,
    dec: &mut Decider,
) -> Result<bool, Failure> {
    if let Curve::Intersection(i) = geom.curve(curve) {
        for &surface in &i.surfaces {
            if !surface_domain(geom.surface(surface), point, dec)? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
pub(super) fn interval_component_contains(
    geom: &GeomStore,
    curve: CurveId,
    point: &Point3,
    plan: Option<&Plan>,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let Some(Frame::ConeTorus(frame)) = plan.and_then(|plan| plan.frame.as_ref()) else {
        return Ok(None);
    };
    let Curve::Intersection(branch) = geom.curve(curve) else {
        unreachable!("interval-set branch");
    };
    let component = frame
        .component(&geom.position(branch.witness), dec)
        .map_err(Refusal::from)?;
    frame
        .contains(component, point, dec)
        .map(Some)
        .map_err(Refusal::from)
}

pub(super) fn on(geom: &GeomStore, c: CurveId, p: PointId, dec: &mut Decider) -> Outcome<Apart> {
    let plan = match validate(geom, c, dec) {
        Ok(p) => p,
        Err(e) => return Outcome::Refused(e),
    };
    if geom
        .facts()
        .holds(FactKind::Incident, GeomKey::Point(p), GeomKey::Curve(c))
    {
        return fact(dec);
    }
    let carrier = geom.point_on_curve(c, p, dec);
    let Outcome::ProvenCoincident(proof) = carrier else {
        return carrier;
    };
    if cone_support(geom, c).is_some() {
        match curve_domain(geom, c, &geom.position(p), dec) {
            Ok(true) => {}
            Ok(false) => return Outcome::Proven(Apart),
            Err(e) => return e.outcome(),
        }
    }
    match interval_component_contains(geom, c, &geom.position(p), plan.as_deref(), dec) {
        Ok(Some(true)) => return Outcome::ProvenCoincident(proof),
        Ok(Some(false)) => return Outcome::Proven(Apart),
        Ok(None) => {}
        Err(refusal) => return Outcome::Refused(refusal),
    }
    let Some(divider) = plan.as_ref().and_then(|p| p.divider.as_ref()) else {
        return Outcome::ProvenCoincident(proof);
    };
    let Curve::Intersection(i) = geom.curve(c) else {
        unreachable!()
    };
    let point_sign = match dec.sign(&divider.value(&geom.position(p))) {
        Outcome::Proven(s) => s,
        Outcome::ProvenCoincident(c) => return Outcome::ProvenCoincident(c),
        Outcome::Refused(e) => return Outcome::Refused(e),
        Outcome::NotExact(mut n) => {
            n.bound = &n.bound / &divider.length_lower;
            return Outcome::NotExact(n);
        }
    };
    match dec.sign(&divider.value(&geom.position(i.witness))) {
        Outcome::Proven(s) if s == point_sign => Outcome::ProvenCoincident(proof),
        Outcome::Proven(_) => Outcome::Proven(Apart),
        Outcome::ProvenCoincident(_) => panic!("regular witness is not on divider"),
        Outcome::Refused(e) => Outcome::Refused(e),
        Outcome::NotExact(mut n) => match dec.enclose(&divider.value(&geom.position(p)), &r(1)) {
            Ok(b) => {
                n.bound = &b.lo().abs().max(b.hi().abs()) / &divider.length_lower;
                Outcome::NotExact(n)
            }
            Err(e) => Outcome::Refused(e),
        },
    }
}
fn sqrt_upper(r: &Rational) -> Rational {
    Ball::from_rational(r, 64).sqrt().unwrap().hi()
}
fn equal_dimensions(a: &Real, b: &Real, dec: &mut Decider) -> Outcome<Apart> {
    if let (Some(a), Some(b)) = (a.as_rational(), b.as_rational()) {
        return if a == b {
            zero()
        } else {
            Outcome::Proven(Apart)
        };
    }
    match dec.sign(&(a - b)).map(|_| Apart) {
        Outcome::NotExact(mut n) => {
            // |sqrt(a) - sqrt(b)| <= sqrt(|a - b|) for nonnegative a, b.
            n.bound = sqrt_upper(&n.bound);
            Outcome::NotExact(n)
        }
        other => other,
    }
}
fn equal_points(a: &Point3, b: &Point3, dec: &mut Decider) -> Outcome<Apart> {
    if let Some(value) = a.minus(b).norm_squared().as_rational() {
        return if value.is_zero() {
            zero()
        } else {
            Outcome::Proven(Apart)
        };
    }
    match dec.sign(&a.minus(b).norm_squared()).map(|_| Apart) {
        Outcome::NotExact(mut n) => {
            n.bound = sqrt_upper(&n.bound);
            Outcome::NotExact(n)
        }
        other => other,
    }
}
fn equal_lines(a: &Point3, b: &Point3, d: &Vector3, dec: &mut Decider) -> Outcome<Apart> {
    match dec.sign(&a.minus(b).cross(d).norm_squared()).map(|_| Apart) {
        Outcome::NotExact(mut n) => {
            if let Some(d) = d.as_rational() {
                let lower = d.x.abs().max(d.y.abs()).max(d.z.abs());
                n.bound = sqrt_upper(&n.bound) / lower;
                Outcome::NotExact(n)
            } else {
                let value = (&Real::from(n.bound.clone()) / &d.norm_squared()).sqrt();
                match dec.enclose_at(&value, 64) {
                    Ok(bound) => {
                        n.bound = bound.hi();
                        Outcome::NotExact(n)
                    }
                    Err(refusal) => Outcome::Refused(refusal),
                }
            }
        }
        other => other,
    }
}
fn cone_support(geom: &GeomStore, curve: CurveId) -> Option<&catallax_shapes::Cone> {
    let Curve::Intersection(i) = geom.curve(curve) else {
        return None;
    };
    i.surfaces.iter().find_map(|s| match geom.surface(*s) {
        Surface::Cone(cone) => Some(cone),
        _ => None,
    })
}
pub(super) fn same(geom: &GeomStore, a: CurveId, b: CurveId, dec: &mut Decider) -> Outcome<Apart> {
    if let Err(e) = validate(geom, a, dec).and_then(|_| validate(geom, b, dec)) {
        return Outcome::Refused(e);
    }
    if geom
        .facts()
        .holds(FactKind::Coincident, GeomKey::Curve(a), GeomKey::Curve(b))
    {
        return fact(dec);
    }
    if geom.curve(a).same_definition(geom.curve(b)) {
        if geom.curve(a).has_rational_definition() {
            return zero();
        }
        return Outcome::ProvenCoincident(catallax_numbers::Coincidence {
            proof: catallax_numbers::ZeroProof::Algebraic,
        });
    }
    if matches!(geom.curve(a), Curve::BSpline(_)) || matches!(geom.curve(b), Curve::BSpline(_)) {
        return Outcome::Refused(Refusal::new(
            "intersect.parametric-coincidence",
            "different spline definitions require a finite curve-equivalence proof",
        ));
    }
    let fa = match carrier(geom, a, dec) {
        Ok(f) => f,
        Err(e) => return Outcome::Refused(e),
    };
    let fb = match carrier(geom, b, dec) {
        Ok(f) => f,
        Err(e) => return Outcome::Refused(e),
    };
    let known = |kind, dec: &mut Decider| {
        let holds = geom
            .facts()
            .holds(kind, GeomKey::Curve(a), GeomKey::Curve(b));
        if holds {
            dec.note_fact_decision();
        }
        holds
    };
    match (fa, fb) {
        (
            CurveForm::Line {
                origin: pa,
                dir: da,
            },
            CurveForm::Line {
                origin: pb,
                dir: db,
            },
        ) => {
            match sign(&da.cross(&db).norm_squared(), dec) {
                Ok(Sign::Zero) => {}
                Ok(_) => return Outcome::Proven(Apart),
                Err(error) => return error.outcome(),
            }
            match (cone_support(geom, a), cone_support(geom, b)) {
                (None, None) => equal_lines(&pa, &pb, &da, dec),
                (Some(ca), Some(cb)) => {
                    let first = match sign(&da.dot(&ca.axis), dec) {
                        Ok(value) => value,
                        Err(error) => return error.outcome(),
                    };
                    let second = match sign(&da.dot(&cb.axis), dec) {
                        Ok(value) => value,
                        Err(error) => return error.outcome(),
                    };
                    if first != second {
                        return Outcome::Proven(Apart);
                    }
                    match equal_points(&ca.apex, &cb.apex, dec) {
                        Outcome::ProvenCoincident(_) => equal_lines(&pa, &pb, &da, dec),
                        other => other,
                    }
                }
                _ => Outcome::Proven(Apart),
            }
        }
        (
            CurveForm::Circle {
                center: a,
                normal: na,
                radius_sq: ra,
                ..
            },
            CurveForm::Circle {
                center: b,
                normal: nb,
                radius_sq: rb,
                ..
            },
        ) => {
            if !known(FactKind::Coaxial, dec) {
                match sign(&na.cross(&nb).norm_squared(), dec) {
                    Ok(Sign::Zero) => {}
                    Ok(_) => return Outcome::Proven(Apart),
                    Err(error) => return error.outcome(),
                }
            }
            match equal_dimensions(&ra, &rb, dec) {
                Outcome::ProvenCoincident(_) if known(FactKind::Concentric, dec) => zero(),
                Outcome::ProvenCoincident(_) => equal_points(&a, &b, dec),
                other => other,
            }
        }
        (
            CurveForm::Ellipse {
                center: a,
                normal: na,
                major_dir: da,
                major_sq: ma,
                minor_sq: sa,
            },
            CurveForm::Ellipse {
                center: b,
                normal: nb,
                major_dir: db,
                major_sq: mb,
                minor_sq: sb,
            },
        ) => {
            for cross in [na.cross(&nb), da.cross(&db)] {
                match sign(&cross.norm_squared(), dec) {
                    Ok(Sign::Zero) => {}
                    Ok(_) => return Outcome::Proven(Apart),
                    Err(error) => return error.outcome(),
                }
            }
            for (a, b) in [(&ma, &mb), (&sa, &sb)] {
                match equal_dimensions(a, b, dec) {
                    Outcome::ProvenCoincident(_) => {}
                    other => return other,
                }
            }
            equal_points(&a, &b, dec)
        }
        (
            CurveForm::Parabola {
                vertex: a,
                normal: na,
                axis_dir: da,
                focal: fa,
            },
            CurveForm::Parabola {
                vertex: b,
                normal: nb,
                axis_dir: db,
                focal: fb,
            },
        ) => {
            for cross in [na.cross(&nb), da.cross(&db)] {
                match sign(&cross.norm_squared(), dec) {
                    Ok(Sign::Zero) => {}
                    Ok(_) => return Outcome::Proven(Apart),
                    Err(error) => return error.outcome(),
                }
            }
            match sign(&da.dot(&db), dec) {
                Ok(Sign::Positive) => {}
                Ok(_) => return Outcome::Proven(Apart),
                Err(error) => return error.outcome(),
            }
            match dec.sign(&(&fa - &fb)).map(|_| Apart) {
                Outcome::ProvenCoincident(_) => equal_points(&a, &b, dec),
                other => other,
            }
        }
        (
            CurveForm::Hyperbola {
                center: a,
                normal: na,
                axis_dir: da,
                transverse_sq: ta,
                conjugate_sq: ca,
                branch: ba,
            },
            CurveForm::Hyperbola {
                center: b,
                normal: nb,
                axis_dir: db,
                transverse_sq: tb,
                conjugate_sq: cb,
                branch: bb,
            },
        ) => {
            for cross in [na.cross(&nb), da.cross(&db)] {
                match sign(&cross.norm_squared(), dec) {
                    Ok(Sign::Zero) => {}
                    Ok(_) => return Outcome::Proven(Apart),
                    Err(error) => return error.outcome(),
                }
            }
            match sign(&da.dot(&db), dec) {
                Ok(value) if (value == Sign::Positive) == (ba == bb) => {}
                Ok(_) => return Outcome::Proven(Apart),
                Err(error) => return error.outcome(),
            }
            for (a, b) in [(&ta, &tb), (&ca, &cb)] {
                match equal_dimensions(a, b, dec) {
                    Outcome::ProvenCoincident(_) => {}
                    other => return other,
                }
            }
            equal_points(&a, &b, dec)
        }
        (CurveForm::Branch, CurveForm::Branch) => {
            let (Curve::Intersection(ia), Curve::Intersection(ib)) = (geom.curve(a), geom.curve(b))
            else {
                unreachable!()
            };
            for permutation in [[0, 1], [1, 0]] {
                let mut matched = true;
                for (j, &other) in permutation.iter().enumerate() {
                    match geom.same_surface(ia.surfaces[j], ib.surfaces[other], dec) {
                        Outcome::Proven(_) => {
                            matched = false;
                            break;
                        }
                        Outcome::ProvenCoincident(_) => {}
                        Outcome::NotExact(n) => return Outcome::NotExact(n),
                        Outcome::Refused(e) => return Outcome::Refused(e),
                    }
                }
                if matched {
                    return on(geom, a, ib.witness, dec);
                }
            }
            Outcome::Proven(Apart)
        }
        _ => Outcome::Proven(Apart),
    }
}
pub(super) fn spline_parameter<'a>(
    curve: &'a BSplineCurve,
    point: &Point3,
    dec: &mut Decider,
) -> Result<SplineCurveParameter<'a>, Refusal> {
    curve.unique_parameter_of(point, dec)?.ok_or_else(|| {
        Refusal::new(
            "intersect.spline-parameter",
            "point is not on the finite spline curve",
        )
    })
}
pub(super) fn order(
    geom: &GeomStore,
    c: CurveId,
    p: PointId,
    q: PointId,
    dec: &mut Decider,
) -> Outcome<StrictOrder> {
    if let Curve::BSpline(curve) = geom.curve(c) {
        let parameters = (|| -> Result<(Real, Real), Refusal> {
            if crate::curve_is_closed(geom, c, dec)? {
                return Err(Refusal::new(
                    "intersect.spline-parameter",
                    "linear order requires an open spline curve",
                ));
            }
            Ok((
                spline_parameter(curve, &geom.position(p), dec)?
                    .value()
                    .clone(),
                spline_parameter(curve, &geom.position(q), dec)?
                    .value()
                    .clone(),
            ))
        })();
        return match parameters {
            Ok((a, b)) => dec.sign(&(b - a)).map(|sign| match sign {
                StrictSign::Positive => StrictOrder::Less,
                StrictSign::Negative => StrictOrder::Greater,
            }),
            Err(error) => Outcome::Refused(error),
        };
    }
    let dir = match carrier(geom, c, dec) {
        Ok(CurveForm::Line { dir, .. }) => dir,
        Ok(
            CurveForm::Parabola {
                normal, axis_dir, ..
            }
            | CurveForm::Hyperbola {
                normal, axis_dir, ..
            },
        ) => normal.cross(&axis_dir),
        Ok(_) => panic!("order requires an open curve"),
        Err(e) => return Outcome::Refused(e),
    };
    let sense = match orientation(geom, c, dec) {
        Ok(s) => s,
        Err(e) => return Outcome::Refused(e),
    };
    dec.sign(&(geom.position(q).minus(&geom.position(p)).dot(&dir.clone()) * sense))
        .map(|s| {
            if s == StrictSign::Positive {
                StrictOrder::Less
            } else {
                StrictOrder::Greater
            }
        })
}
pub(super) fn charts(
    geom: &GeomStore,
    c: CurveId,
    dec: &mut Decider,
) -> Result<Vec<Chart>, Failure> {
    let plan = validate(geom, c, dec)?;
    let f = form(geom, c, dec)?;
    if matches!(f, CurveForm::Branch) {
        chart::branch_charts(geom, c, &plan.unwrap(), dec)
    } else if matches!(f, CurveForm::BSpline(_)) {
        Err(Refusal::new(
            "intersect.spline-chart",
            "finite spline tracing requires explicit knot-span ranges",
        )
        .into())
    } else {
        Ok(chart::form_charts(&f))
    }
}
pub(super) fn cyclic(
    geom: &GeomStore,
    c: CurveId,
    ids: [PointId; 3],
    dec: &mut Decider,
) -> Outcome<Cyclic> {
    let f = match carrier(geom, c, dec) {
        Ok(f) => f,
        Err(e) => return Outcome::Refused(e),
    };
    let [p, q, rp] = ids.map(|p| geom.position(p));
    match f {
        CurveForm::Circle { normal, .. } | CurveForm::Ellipse { normal, .. } => {
            let sense = match orientation(geom, c, dec) {
                Ok(s) => s,
                Err(e) => return Outcome::Refused(e),
            };
            dec.sign(&(q.minus(&p).cross(&rp.minus(&p)).dot(&normal.clone()) * sense))
                .map(|s| {
                    if s == StrictSign::Positive {
                        Cyclic::Forward
                    } else {
                        Cyclic::Backward
                    }
                })
        }
        CurveForm::Line { .. } | CurveForm::Parabola { .. } | CurveForm::Hyperbola { .. } => {
            panic!("cyclic order requires closed curve")
        }
        CurveForm::BSpline(curve) => {
            let result = (|| -> Result<Real, Refusal> {
                if !crate::curve_is_closed(geom, c, dec)? {
                    return Err(Refusal::new(
                        "intersect.spline-parameter",
                        "cyclic order requires a closed spline curve",
                    ));
                }
                let a = spline_parameter(&curve, &p, dec)?.value().clone();
                let b = spline_parameter(&curve, &q, dec)?.value().clone();
                let d = spline_parameter(&curve, &rp, dec)?.value().clone();
                Ok((&b - &a) * ((&d - &a) * (&d - &b)))
            })();
            match result {
                Ok(value) => dec.sign(&value).map(|sign| match sign {
                    StrictSign::Positive => Cyclic::Forward,
                    StrictSign::Negative => Cyclic::Backward,
                }),
                Err(error) => Outcome::Refused(error),
            }
        }
        CurveForm::Branch => {
            let result = (|| -> Result<Outcome<Cyclic>, Failure> {
                let charts = charts(geom, c, dec)?;
                let positions = [p, q, rp];
                let mut locations = vec![];
                for p in &positions {
                    locations.push(chart::locate(&charts, p, dec)?);
                }
                let mut comparisons = [[std::cmp::Ordering::Equal; 3]; 3];
                for i in 0..3 {
                    for j in i + 1..3 {
                        let cmp = if locations[i].0 == locations[j].0 {
                            match sign(&(&locations[j].1 - &locations[i].1), dec)? {
                                Sign::Positive => std::cmp::Ordering::Less,
                                Sign::Negative => std::cmp::Ordering::Greater,
                                Sign::Zero => return Ok(zero()),
                            }
                        } else {
                            locations[i].0.cmp(&locations[j].0)
                        };
                        comparisons[i][j] = cmp;
                        comparisons[j][i] = cmp.reverse();
                    }
                }
                let less = |i: usize, j: usize| comparisons[i][j] == std::cmp::Ordering::Less;
                let forward = (less(0, 1) && less(1, 2))
                    || (less(1, 2) && less(2, 0))
                    || (less(2, 0) && less(0, 1));
                Ok(Outcome::Proven(if forward {
                    Cyclic::Forward
                } else {
                    Cyclic::Backward
                }))
            })();
            match result {
                Ok(o) => o,
                Err(e) => e.outcome(),
            }
        }
    }
}
pub(super) fn tangent(
    geom: &GeomStore,
    c: CurveId,
    p: PointId,
    dec: &mut Decider,
) -> Result<[Vector3; 2], Refusal> {
    validate(geom, c, dec)?;
    let position = geom.position(p);
    if let Curve::BSpline(curve) = geom.curve(c) {
        return spline_parameter(curve, &position, dec)?.tangents(dec);
    }
    if let Curve::Intersection(i) = geom.curve(c) {
        // Tangency is local to its recorded component, not the entire pair.
        let tangent = geom.facts().iter().any(|f| {
            f.kind == FactKind::Tangent
                && ((f.a == GeomKey::Surface(i.surfaces[0])
                    && f.b == GeomKey::Surface(i.surfaces[1]))
                    || (f.b == GeomKey::Surface(i.surfaces[0])
                        && f.a == GeomKey::Surface(i.surfaces[1])))
                && match f.at {
                    Some(GeomKey::Point(at)) => at == p,
                    Some(GeomKey::Curve(at)) => {
                        at == c
                            || geom.facts().holds(
                                FactKind::Incident,
                                GeomKey::Point(p),
                                GeomKey::Curve(at),
                            )
                    }
                    _ => false,
                }
        });
        let mut regular = false;
        let cross = if tangent {
            dec.note_fact_decision();
            None
        } else {
            let cross = geom
                .surface(i.surfaces[0])
                .normal_at(&position, dec)?
                .cross(&geom.surface(i.surfaces[1]).normal_at(&position, dec)?);
            for component in [&cross.x, &cross.y, &cross.z] {
                if sign(component, dec)? != Sign::Zero {
                    regular = true;
                    break;
                }
            }
            Some(cross)
        };
        if regular {
            let cross = cross.expect("regular gradient cross product");
            return Ok([cross.clone(), cross]);
        }
        let closed_form = form(geom, c, dec)?;
        if !matches!(closed_form, CurveForm::Branch) {
            let t = analytic_tangent(closed_form, &position);
            return Ok([t.clone(), t]);
        }
        let charts = charts(geom, c, dec)?;
        let mut arriving = None;
        let mut leaving = None;
        for ch in charts {
            for (t, is_end) in [(Real::zero(), false), (Real::one(), true)] {
                let here = ch.position(&t);
                if sign(&here.minus(&position).norm_squared(), dec)? == Sign::Zero {
                    let v = ch.endpoint_tangent(is_end);
                    if is_end {
                        arriving = Some(v);
                    } else {
                        leaving = Some(v);
                    }
                }
            }
        }
        return Ok([
            arriving.expect("singular end of branch"),
            leaving.expect("singular start of branch"),
        ]);
    }
    let t = analytic_tangent(form(geom, c, dec)?, &position);
    Ok([t.clone(), t])
}
fn analytic_tangent(form: CurveForm, position: &Point3) -> Vector3 {
    match form {
        CurveForm::Line { dir, .. } => dir,
        CurveForm::Circle { center, normal, .. } => normal.cross(&position.minus(&center)),
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            let d = position.minus(&center);
            let v = normal.cross(&major_dir);
            let u = major_dir.clone().scaled(
                &(d.dot(&major_dir.clone()) / (real(&major_dir.norm_squared()) * major_sq)),
            );
            let v = v
                .clone()
                .scaled(&(d.dot(&v.clone()) / (real(&v.norm_squared()) * minor_sq)));
            normal.clone().cross(&u.plus(&v))
        }
        CurveForm::Parabola {
            vertex,
            normal,
            axis_dir,
            focal,
        } => {
            let u = unit(&axis_dir);
            let v = unit(&normal.cross(&axis_dir));
            let y = position.minus(&vertex).dot(&v);
            u.scaled(&(y / (Real::from_i64(2) * focal))).plus(&v)
        }
        CurveForm::Hyperbola {
            center,
            normal,
            axis_dir,
            transverse_sq,
            conjugate_sq,
            ..
        } => {
            let u = unit(&axis_dir);
            let v = unit(&normal.cross(&axis_dir));
            let d = position.minus(&center);
            let derivative = transverse_sq * d.dot(&v) / (conjugate_sq * d.dot(&u));
            u.scaled(&derivative).plus(&v)
        }
        CurveForm::Branch => unreachable!("traced branch tangent"),
        CurveForm::BSpline(_) => unreachable!("finite spline tangent dispatched by preimage"),
    }
}
