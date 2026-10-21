//! Rational probes obtained from rational surface parametrisations.

use catallax_intersect::{arc_point, tangent_at};
use catallax_numbers::{
    Answer, Answered, Decider, Fact, FactKind, FactOrigin, Integer, Outcome, Rational, Real,
    Refusal, Sign,
};
use catallax_shapes::{
    GeomKey, Meet, MeetOf, Point3, PointDef, PointId, Position, RPoint3, Rotation, Surface,
};
use catallax_topology::{FaceId, Loop, Model, Sense};

use crate::geometry::{natural_ends, oriented, product, r, refusal, unit};

fn integer_sqrt(n: &Integer) -> Option<Integer> {
    if n.sign() == Sign::Negative {
        return None;
    }
    if n.is_zero() {
        return Some(Integer::zero());
    }
    let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
    loop {
        let q = n.div_rem_floor(&x)?.0;
        let y = (&x + &q).div_rem_floor(&Integer::from_i64(2))?.0;
        if y >= x {
            return (&x * &x == *n).then_some(x);
        }
        x = y;
    }
}

pub(crate) fn sqrt_rational(n: &Rational) -> Option<Rational> {
    Rational::new(integer_sqrt(n.numerator())?, integer_sqrt(n.denominator())?)
}

pub(crate) fn curve_sample(
    model: &mut Model,
    curve: catallax_shapes::CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Result<PointId, Refusal> {
    if from.is_none()
        && to.is_none()
        && let catallax_shapes::Curve::Circle(circle) = model.geom().curve(curve)
        && let Some(scale) = sqrt_rational(&(&circle.radius_sq / &circle.ref_dir.norm_squared()))
    {
        let point = circle.center.plus(&circle.ref_dir.scaled(&scale));
        return Ok(model.geom_mut().add_point(PointDef::Explicit(point)));
    }
    let answer = arc_point(model.geom(), curve, from, to, dec)?;
    if !answer.marks.is_empty() {
        return Err(refusal(
            "boolean.imprint-ambiguous",
            "a marked arc probe cannot establish incidence",
        ));
    }
    Ok(model.geom_mut().add_point(PointDef::Meet(Meet {
        of: MeetOf::Curves(curve, curve),
        position: Position::Closed(answer.value),
    })))
}

fn rationalise(
    surface: &Surface,
    target: &Point3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<RPoint3>, Refusal> {
    let approximate = |p: &Point3, dec: &mut Decider| -> Result<RPoint3, Refusal> {
        Ok(RPoint3::new(
            dec.enclose(&p.x, width)?.mid(),
            dec.enclose(&p.y, width)?.mid(),
            dec.enclose(&p.z, width)?.mid(),
        ))
    };
    match surface {
        Surface::Plane(s) => {
            let p = approximate(target, dec)?;
            let displacement = p.minus(&s.origin);
            Ok(Some(p.plus(&s.normal.scaled(
                &(-displacement.dot(&s.normal) / s.normal.norm_squared()),
            ))))
        }
        Surface::Cylinder(s) => {
            let Some(length) = sqrt_rational(&s.ref_dir.norm_squared()) else {
                return Ok(None);
            };
            let base = s.ref_dir.scaled(&(&s.radius / &length));
            let tangent = s.axis.cross(&base);
            let squared = base.norm_squared();
            let tangent_squared = tangent.norm_squared();
            let d = target.minus(&s.origin.to_real());
            let h = d.dot(&s.axis.to_real()) / Real::from(s.axis.norm_squared());
            for sense in [Sense::Same, Sense::Reversed] {
                let b = if sense == Sense::Same {
                    base.clone()
                } else {
                    base.scaled(&r(-1))
                };
                let v = s.axis.cross(&b);
                let denominator = Real::from(squared.clone()) + d.dot(&b.to_real());
                match dec.sign(&denominator) {
                    Outcome::ProvenCoincident(_) | Outcome::NotExact(_) => continue,
                    Outcome::Refused(e) => return Err(e),
                    Outcome::Proven(_) => {}
                }
                let t =
                    (d.dot(&v.to_real()) / denominator) * Real::from(&squared / &tangent_squared);
                let t = dec.enclose(&t, width)?.mid();
                let h = dec.enclose(&h, width)?.mid();
                let z = &tangent_squared * &(&t * &t);
                let denominator = &squared + &z;
                let radial = b
                    .scaled(&(&(&squared - &z) / &denominator))
                    .plus(&v.scaled(&((&squared * &t) * r(2) / denominator)));
                return Ok(Some(s.origin.plus(&radial).plus(&s.axis.scaled(&h))));
            }
            Ok(None)
        }
        Surface::Sphere(s) => {
            let d = target.minus(&s.center.to_real());
            for choice in 0..3 {
                let rotation = Rotation::from_quaternion(
                    r(1),
                    r(choice),
                    r(choice * choice),
                    r(choice * choice * choice),
                )?;
                let u = rotation.apply(&unit(0));
                let v = rotation.apply(&unit(1));
                let pole = rotation.apply(&unit(2));
                let denominator = Real::from(s.radius.clone()) - d.dot(&pole.to_real());
                match dec.sign(&denominator) {
                    Outcome::ProvenCoincident(_) => continue,
                    Outcome::Refused(e) => return Err(e),
                    Outcome::NotExact(_) => continue,
                    Outcome::Proven(_) => {}
                }
                let x = dec
                    .enclose(&(&d.dot(&u.to_real()) / &denominator), width)?
                    .mid();
                let y = dec
                    .enclose(&(&d.dot(&v.to_real()) / &denominator), width)?
                    .mid();
                let norm = &x * &x + &y * &y;
                let factor = &s.radius / &(Rational::one() + norm.clone());
                let radial = u
                    .scaled(&(r(2) * x))
                    .plus(&v.scaled(&(r(2) * y)))
                    .plus(&pole.scaled(&(norm - Rational::one())));
                return Ok(Some(s.center.plus(&radial.scaled(&factor))));
            }
            Ok(None)
        }
        _ => Ok(None),
    }
}

fn search<T>(
    model: &mut Model,
    face: FaceId,
    dec: &mut Decider,
    barriers: &[catallax_topology::EdgeId],
    mut accept: impl FnMut(&mut Model, PointId, &RPoint3, &mut Decider) -> Result<Option<T>, Refusal>,
) -> Result<T, Refusal> {
    let f = model.face(face).clone();
    let surface = model.geom().surface(f.surface).clone();
    if f.loops.is_empty() {
        let Surface::Sphere(s) = surface else {
            return Err(refusal(
                "boolean.no-rational-point",
                "unbounded or unsupported whole face",
            ));
        };
        for k in 0..crate::geometry::DIRECTION_ATTEMPTS {
            let rotation = Rotation::from_quaternion(r(1), r(i64::from(k)), r(1), r(0))?;
            let p = s.center.plus(&rotation.apply(&unit(0)).scaled(&s.radius));
            let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
            if let Some(value) = accept(model, point, &p, dec)? {
                return Ok(value);
            }
        }
        return Err(refusal(
            "boolean.no-rational-point",
            "no whole-face probe settled the query",
        ));
    }
    for &loop_id in &f.loops {
        let Loop::Coedges(cs) = model.loop_(loop_id) else {
            continue;
        };
        for coedge_id in cs.clone() {
            let coedge = model.coedge(coedge_id).clone();
            let edge = model.edge(coedge.edge).clone();
            let (from, to) = natural_ends(model, coedge.edge);
            let point = curve_sample(model, edge.curve, from, to, dec)?;
            let p = model.geom().position(point);
            let tangent = tangent_at(model.geom(), edge.curve, point, dec)?[1].clone();
            let tangent = oriented(&tangent, product(edge.sense, coedge.sense));
            let inward = oriented(&surface.gradient(&p), f.sense).cross(&tangent);
            let scale = match &surface {
                Surface::Cylinder(s) => s.radius.clone(),
                Surface::Sphere(s) => s.radius.clone(),
                _ => {
                    if let (Some(a), Some(b)) = (from, to) {
                        let delta = model.geom().position(a).minus(&model.geom().position(b));
                        let bounds = [&delta.x, &delta.y, &delta.z]
                            .map(|x| dec.enclose(x, &Rational::one()));
                        let mut extent = Rational::zero();
                        for bound in bounds {
                            let b = bound?;
                            extent = extent.max(b.lo().abs().max(b.hi().abs()));
                        }
                        if extent.is_zero() {
                            Rational::one()
                        } else {
                            extent
                        }
                    } else {
                        Rational::one()
                    }
                }
            };
            let mut upper = Rational::zero();
            for x in [&inward.x, &inward.y, &inward.z] {
                let b = dec.enclose(x, &scale)?;
                upper = upper.max(b.lo().abs().max(b.hi().abs()));
            }
            if upper.is_zero() {
                continue;
            }
            for depth in 0..crate::geometry::PROBE_DEPTHS {
                let step = scale.mul_pow2(-i64::from(depth));
                let target = p.plus(&inward.scaled(&Real::from(&step / &upper)));
                let width = step.mul_pow2(-8);
                let Some(candidate) = rationalise(&surface, &target, &width, dec)? else {
                    continue;
                };
                let query = model
                    .geom_mut()
                    .add_point(PointDef::Explicit(candidate.clone()));
                match model.geom().side(f.surface, query, dec) {
                    Outcome::ProvenCoincident(proof) => {
                        model.geom_mut().record_fact(Fact {
                            kind: FactKind::Incident,
                            a: GeomKey::Point(query),
                            b: GeomKey::Surface(f.surface),
                            at: None,
                            origin: FactOrigin::Proof(proof.proof),
                        });
                    }
                    Outcome::Refused(e) => return Err(e),
                    _ => continue,
                }
                match crate::probe_path::adjacent(
                    model,
                    &surface,
                    (point, &p, edge.curve),
                    &candidate,
                    &inward,
                    barriers,
                    dec,
                ) {
                    Ok(true) => {}
                    Ok(false) => continue,
                    Err(e)
                        if e.code == "boolean.imprint-ambiguous"
                            || e.code == "boolean.classify-ambiguous" =>
                    {
                        continue;
                    }
                    Err(e) => return Err(e),
                }
                if let Some(result) = accept(model, query, &candidate, dec)? {
                    return Ok(result);
                }
            }
        }
    }
    Err(refusal(
        "boolean.no-rational-point",
        "no rational probe was certified strictly inside the face",
    ))
}

pub(crate) fn find<T>(
    model: &mut Model,
    face: FaceId,
    dec: &mut Decider,
    accept: impl FnMut(&mut Model, PointId, &RPoint3, &mut Decider) -> Result<Option<T>, Refusal>,
) -> Result<T, Refusal> {
    let mut barriers = Vec::new();
    for &loop_id in &model.face(face).loops {
        if let Loop::Coedges(coedges) = model.loop_(loop_id) {
            barriers.extend(coedges.iter().map(|c| model.coedge(*c).edge));
        }
    }
    barriers.sort_unstable();
    barriers.dedup();
    // A boundary coedge has the face on its left. A path starting strictly
    // left and meeting no boundary certifies the witness without re-solving
    // global face membership at that same witness.
    search(model, face, dec, &barriers, accept)
}

pub(crate) fn left_probe<T>(
    model: &mut Model,
    face: FaceId,
    barriers: &[catallax_topology::EdgeId],
    dec: &mut Decider,
    accept: impl FnMut(&mut Model, PointId, &RPoint3, &mut Decider) -> Result<Option<T>, Refusal>,
) -> Result<T, Refusal> {
    search(model, face, dec, barriers, accept)
}

pub(crate) fn face_point(model: &mut Model, face: FaceId, dec: &mut Decider) -> Answered<RPoint3> {
    find(model, face, dec, |_, _, p, _| Ok(Some(p.clone()))).map(Answer::exact)
}
