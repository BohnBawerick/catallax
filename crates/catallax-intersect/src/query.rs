use crate::{
    CurveForm, Cyclic,
    chart::{self, Chart},
    math::*,
    pair::{self, Plan},
};
use catallax_numbers::{
    Ball, Decider, FactKind, Outcome, Rational, Real, Refusal, Sign, StrictOrder, StrictSign,
};
use catallax_shapes::{
    Apart, Curve, CurveId, GeomKey, GeomStore, Point3, PointId, Surface, Vector3,
};
use std::sync::Arc;

pub(super) fn validate(geom: &GeomStore, c: CurveId) -> Result<Option<Arc<Plan>>, Refusal> {
    match geom.curve(c) {
        Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) => Ok(None),
        Curve::Intersection(i) => {
            if !crate::supports(
                geom.surface(i.surfaces[0]).kind(),
                geom.surface(i.surfaces[1]).kind(),
            ) {
                return Err(refused_curve());
            }
            let p = pair::plan(geom.surface(i.surfaces[0]), geom.surface(i.surfaces[1]))
                .map_err(|_| refused_curve())?;
            if matches!((geom.surface(i.surfaces[0]),geom.surface(i.surfaces[1])),(Surface::Cylinder(a),Surface::Cylinder(b)) if !a.axis.cross(&b.axis).is_zero() && a.radius==b.radius && p.frame.is_none())
            {
                return Err(refused_curve());
            }
            Ok(Some(p))
        }
        _ => Err(refused_curve()),
    }
}
pub(super) fn form(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<CurveForm, Refusal> {
    let p = validate(geom, c)?;
    Ok(match geom.curve(c) {
        Curve::Line(l) => CurveForm::Line {
            origin: l.origin.to_real(),
            dir: l.dir.clone(),
        },
        Curve::Circle(l) => CurveForm::Circle {
            center: l.center.to_real(),
            normal: l.normal.clone(),
            radius_sq: real(&l.radius_sq),
            ref_dir: l.ref_dir.clone(),
        },
        Curve::Ellipse(l) => CurveForm::Ellipse {
            center: l.center.to_real(),
            normal: l.normal.clone(),
            major_dir: l.major_dir.clone(),
            major_sq: real(&l.major_sq),
            minor_sq: real(&l.minor_sq),
        },
        Curve::Intersection(i) => {
            let plan = p.unwrap();
            if plan.frame.is_some() {
                return Ok(CurveForm::Branch);
            }
            let point = geom.position(i.witness);
            let index = if let Some(d) = &plan.divider {
                if sign(&d.value(&point), dec)? == Sign::Negative {
                    1
                } else {
                    0
                }
            } else {
                0
            };
            let mut form = plan.curves[index].form.clone();
            if let CurveForm::Line { origin, .. } = &mut form {
                *origin = point;
            }
            form
        }
        _ => return Err(refused_curve()),
    })
}
// Comparing point sets does not require choosing a natural orientation.
fn carrier(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<CurveForm, Refusal> {
    let Some(plan) = validate(geom, c)? else {
        return form(geom, c, dec);
    };
    if plan.frame.is_some() {
        return Ok(CurveForm::Branch);
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
            let base = sphere.center.to_real();
            let z = witness.minus(&base).dot(&normal.to_real()) / real(&normal.norm_squared());
            *center = base.plus(&normal.to_real().scaled(&z));
        }
        _ => {}
    }
    Ok(shape)
}
fn orientation(geom: &GeomStore, c: CurveId) -> Result<Real, Refusal> {
    let Some(plan) = validate(geom, c)? else {
        return Ok(Real::one());
    };
    let Some(divider) = &plan.divider else {
        return Ok(Real::one());
    };
    let Curve::Intersection(i) = geom.curve(c) else {
        unreachable!()
    };
    let value = divider.value(&geom.position(i.witness));
    Ok(&value / &value.square().sqrt())
}
pub(super) fn on(geom: &GeomStore, c: CurveId, p: PointId, dec: &mut Decider) -> Outcome<Apart> {
    let plan = match validate(geom, c) {
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
    match dec.sign(&a.minus(b).norm_squared()).map(|_| Apart) {
        Outcome::NotExact(mut n) => {
            n.bound = sqrt_upper(&n.bound);
            Outcome::NotExact(n)
        }
        other => other,
    }
}
fn equal_lines(
    a: &Point3,
    b: &Point3,
    d: &catallax_shapes::RVector3,
    dec: &mut Decider,
) -> Outcome<Apart> {
    match dec
        .sign(&a.minus(b).cross(&d.to_real()).norm_squared())
        .map(|_| Apart)
    {
        Outcome::NotExact(mut n) => {
            let lower = d.x.abs().max(d.y.abs()).max(d.z.abs());
            n.bound = sqrt_upper(&n.bound) / lower;
            Outcome::NotExact(n)
        }
        other => other,
    }
}
pub(super) fn same(geom: &GeomStore, a: CurveId, b: CurveId, dec: &mut Decider) -> Outcome<Apart> {
    if let Err(e) = validate(geom, a).and_then(|_| validate(geom, b)) {
        return Outcome::Refused(e);
    }
    if geom
        .facts()
        .holds(FactKind::Coincident, GeomKey::Curve(a), GeomKey::Curve(b))
    {
        return fact(dec);
    }
    if !matches!(geom.curve(a), Curve::Intersection(_))
        && !matches!(geom.curve(b), Curve::Intersection(_))
    {
        let same = match (geom.curve(a), geom.curve(b)) {
            (Curve::Line(a), Curve::Line(b)) => {
                a.dir.cross(&b.dir).is_zero() && a.origin.minus(&b.origin).cross(&a.dir).is_zero()
            }
            (Curve::Circle(a), Curve::Circle(b)) => {
                a.center == b.center
                    && a.radius_sq == b.radius_sq
                    && a.normal.cross(&b.normal).is_zero()
            }
            (Curve::Ellipse(a), Curve::Ellipse(b)) => {
                a.center == b.center
                    && a.major_sq == b.major_sq
                    && a.minor_sq == b.minor_sq
                    && a.normal.cross(&b.normal).is_zero()
                    && a.major_dir.cross(&b.major_dir).is_zero()
            }
            _ => false,
        };
        return if same { zero() } else { Outcome::Proven(Apart) };
    }
    let fa = match carrier(geom, a, dec) {
        Ok(f) => f,
        Err(e) => return Outcome::Refused(e),
    };
    let fb = match carrier(geom, b, dec) {
        Ok(f) => f,
        Err(e) => return Outcome::Refused(e),
    };
    match (fa, fb) {
        (CurveForm::Line { origin: a, dir: da }, CurveForm::Line { origin: b, dir: db }) => {
            if da.cross(&db).is_zero() {
                equal_lines(&a, &b, &da, dec)
            } else {
                Outcome::Proven(Apart)
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
            if !na.cross(&nb).is_zero() {
                return Outcome::Proven(Apart);
            }
            match equal_dimensions(&ra, &rb, dec) {
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
            if !na.cross(&nb).is_zero() || !da.cross(&db).is_zero() {
                return Outcome::Proven(Apart);
            }
            for (a, b) in [(&ma, &mb), (&sa, &sb)] {
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
pub(super) fn order(
    geom: &GeomStore,
    c: CurveId,
    p: PointId,
    q: PointId,
    dec: &mut Decider,
) -> Outcome<StrictOrder> {
    let dir = match carrier(geom, c, dec) {
        Ok(CurveForm::Line { dir, .. }) => dir,
        Ok(_) => panic!("order requires line"),
        Err(e) => return Outcome::Refused(e),
    };
    let sense = match orientation(geom, c) {
        Ok(s) => s,
        Err(e) => return Outcome::Refused(e),
    };
    dec.sign(
        &(geom
            .position(q)
            .minus(&geom.position(p))
            .dot(&dir.to_real())
            * sense),
    )
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
    let plan = validate(geom, c)?;
    let f = form(geom, c, dec)?;
    if matches!(f, CurveForm::Branch) {
        chart::branch_charts(geom, c, &plan.unwrap(), dec)
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
            let sense = match orientation(geom, c) {
                Ok(s) => s,
                Err(e) => return Outcome::Refused(e),
            };
            dec.sign(&(q.minus(&p).cross(&rp.minus(&p)).dot(&normal.to_real()) * sense))
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
    validate(geom, c)?;
    let position = geom.position(p);
    if let Curve::Intersection(i) = geom.curve(c) {
        let cross = geom
            .surface(i.surfaces[0])
            .gradient(&position)
            .cross(&geom.surface(i.surfaces[1]).gradient(&position));
        let mut regular = false;
        for component in [&cross.x, &cross.y, &cross.z] {
            if sign(component, dec)? != Sign::Zero {
                regular = true;
                break;
            }
        }
        if regular {
            return Ok([cross.clone(), cross]);
        }
        let charts = charts(geom, c, dec)?;
        let mut arriving = None;
        let mut leaving = None;
        for ch in charts {
            for (t, is_end) in [(Real::zero(), false), (Real::one(), true)] {
                let here = ch.position(&t);
                if sign(&here.minus(&position).norm_squared(), dec)? == Sign::Zero {
                    let j = ch.exact_jets(&t, 1);
                    let v = Vector3::new(j[1][0].clone(), j[1][1].clone(), j[1][2].clone());
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
    let t = match form(geom, c, dec)? {
        CurveForm::Line { dir, .. } => dir.to_real(),
        CurveForm::Circle { center, normal, .. } => {
            normal.to_real().cross(&position.minus(&center))
        }
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            let d = position.minus(&center);
            let v = normal.cross(&major_dir);
            let u = major_dir.to_real().scaled(
                &(d.dot(&major_dir.to_real()) / (real(&major_dir.norm_squared()) * major_sq)),
            );
            let v = v
                .to_real()
                .scaled(&(d.dot(&v.to_real()) / (real(&v.norm_squared()) * minor_sq)));
            normal.to_real().cross(&u.plus(&v))
        }
        CurveForm::Branch | CurveForm::Parabola { .. } | CurveForm::Hyperbola { .. } => {
            unreachable!("non-analytic forms are stored as intersection curves")
        }
    };
    Ok([t.clone(), t])
}
