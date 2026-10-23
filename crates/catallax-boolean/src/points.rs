//! Rational probes obtained from rational surface parametrisations.

use catallax_intersect::{arc_point, tangent_at};
use catallax_numbers::{
    Answer, Answered, Decider, Fact, FactKind, FactOrigin, Integer, Outcome, Rational, Real,
    Refusal, Sign,
};
use catallax_shapes::{
    GeomKey, Meet, MeetOf, Point3, PointDef, PointId, Position, RPoint3, Rotation, Surface, Torus,
    Vector3,
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
pub(crate) fn analytic_frame(axis: &Vector3, reference: &Vector3) -> Option<[Vector3; 3]> {
    let axis_sq = axis.norm_squared();
    let reference_sq = reference.norm_squared();
    let (axis_length, reference_length) =
        if axis.rational_components().is_some() && reference.rational_components().is_some() {
            // Preserve the native frame boundary; algebraic definitions use exact roots.
            (
                Real::from(sqrt_rational(axis_sq.rational_ref()?)?),
                Real::from(sqrt_rational(reference_sq.rational_ref()?)?),
            )
        } else {
            (root(&axis_sq), root(&reference_sq))
        };
    let axis = axis.scaled(&(Real::one() / axis_length));
    let u = reference.scaled(&(Real::one() / reference_length));
    let v = axis.cross(&u);
    Some([u, v, axis])
}
fn rational_circle(t: &Rational) -> [Rational; 2] {
    let squared = t * t;
    let denominator = &r(1) + &squared;
    [
        &(&r(1) - &squared) / &denominator,
        &(&r(2) * t) / &denominator,
    ]
}
pub(crate) fn stereographic(
    x: &Real,
    y: &Real,
    radius: &Real,
    side: i64,
    dec: &mut Decider,
) -> Result<Option<Real>, Refusal> {
    let denominator = radius + &(&Real::from_i64(side) * x);
    match dec.sign(&denominator) {
        Outcome::ProvenCoincident(_) | Outcome::NotExact(_) => Ok(None),
        Outcome::Refused(e) => Err(e),
        Outcome::Proven(_) => Ok(Some((&Real::from_i64(side) * y) / denominator)),
    }
}
fn rational_phase(
    x: &Real,
    y: &Real,
    radius: &Real,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<[Rational; 2]>, Refusal> {
    for side in [1, -1] {
        let Some(t) = stereographic(x, y, radius, side, dec)? else {
            continue;
        };
        let t = dec.enclose(&t, width)?.mid();
        let [cos, sin] = rational_circle(&t);
        return Ok(Some(if side > 0 { [cos, sin] } else { [-cos, -sin] }));
    }
    Ok(None)
}
fn torus_point(torus: &Torus, frame: &[Vector3; 3], u: [Rational; 2], v: [Rational; 2]) -> Point3 {
    let u = u.map(Real::from);
    let v = v.map(Real::from);
    let radial = frame[0].scaled(&u[0]).plus(&frame[1].scaled(&u[1]));
    let radius = &torus.major + &(&torus.minor * &v[0]);
    torus
        .center
        .plus(&radial.scaled(&radius))
        .plus(&frame[2].scaled(&(&torus.minor * &v[1])))
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
    {
        let square = &circle.radius_sq / &circle.ref_dir.norm_squared();
        if let Some(scale) = square.as_rational().and_then(|q| sqrt_rational(&q)) {
            let point = circle
                .center
                .plus(&circle.ref_dir.scaled(&Real::from(scale)));
            return Ok(model.geom_mut().add_point(PointDef::Explicit(point)));
        }
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

pub(crate) fn root(value: &Real) -> Real {
    value
        .as_rational()
        .and_then(|q| sqrt_rational(&q))
        .map_or_else(|| value.sqrt(), Real::from)
}

fn exact_probe(
    surface: &Surface,
    target: &Point3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    match surface {
        Surface::Plane(s) => {
            let p = RPoint3::new(
                dec.enclose(&target.x, width)?.mid(),
                dec.enclose(&target.y, width)?.mid(),
                dec.enclose(&target.z, width)?.mid(),
            )
            .to_real();
            let displacement = p.minus(&s.origin);
            Ok(Some(p.plus(&s.normal.scaled(
                &(-displacement.dot(&s.normal) / s.normal.norm_squared()),
            ))))
        }
        Surface::Cylinder(s) => {
            let base = s.radial_reference();
            let tangent = s.axis.cross(&base);
            let squared = base.norm_squared();
            let tangent_squared = tangent.norm_squared();
            let d = target.minus(&s.origin);
            let h = d.dot(&s.axis) / s.axis.norm_squared();
            let along = d.dot(&base);
            // Keep the stereographic denominator away from the opposite pole,
            // including when the target is a rounded point just off the surface.
            let senses = if crate::geometry::sign(&along, dec, "boolean.imprint-ambiguous")?
                == Sign::Negative
            {
                [Sense::Reversed, Sense::Same]
            } else {
                [Sense::Same, Sense::Reversed]
            };
            for sense in senses {
                let b = if sense == Sense::Same {
                    base.clone()
                } else {
                    base.scaled(&Real::from_i64(-1))
                };
                let v = s.axis.cross(&b);
                let denominator = if sense == Sense::Same {
                    &squared + &along
                } else {
                    &squared - &along
                };
                match dec.sign(&denominator) {
                    Outcome::ProvenCoincident(_) | Outcome::NotExact(_) => continue,
                    Outcome::Refused(e) => return Err(e),
                    Outcome::Proven(_) => {}
                }
                let t = (d.dot(&v) / denominator) * (&squared / &tangent_squared);
                let t = Real::from(dec.enclose(&t, width)?.mid());
                let h = Real::from(dec.enclose(&h, width)?.mid());
                let z = &tangent_squared * &t.square();
                let denominator = &squared + &z;
                let radial = b
                    .scaled(&(&(&squared - &z) / &denominator))
                    .plus(&v.scaled(&((&squared * &t) * Real::from_i64(2) / denominator)));
                return Ok(Some(s.origin.plus(&radial).plus(&s.axis.scaled(&h))));
            }
            Ok(None)
        }
        Surface::Cone(s) => {
            let Some([u, v, axis]) = analytic_frame(&s.axis, &s.ref_dir) else {
                return Ok(None);
            };
            let d = target.minus(&s.apex);
            let height = dec.enclose(&d.dot(&axis), width)?.mid();
            if height.sign() != Sign::Positive {
                return Ok(None);
            }
            let height = Real::from(height);
            let radius = &s.slope * &height;
            let Some([cos, sin]) = rational_phase(&d.dot(&u), &d.dot(&v), &radius, width, dec)?
            else {
                return Ok(None);
            };
            let radial = u
                .scaled(&Real::from(cos))
                .plus(&v.scaled(&Real::from(sin)))
                .scaled(&radius);
            Ok(Some(s.apex.plus(&axis.scaled(&height)).plus(&radial)))
        }
        Surface::Torus(s) => {
            let Some(frame) = analytic_frame(&s.axis, &s.ref_dir) else {
                return Ok(None);
            };
            let d = target.minus(&s.center);
            let x = d.dot(&frame[0]);
            let y = d.dot(&frame[1]);
            let z = d.dot(&frame[2]);
            let radius = (x.square() + y.square()).sqrt();
            let Some(u) = rational_phase(&x, &y, &radius, width, dec)? else {
                return Ok(None);
            };
            let Some(v) = rational_phase(&(&radius - &s.major), &z, &s.minor, width, dec)? else {
                return Ok(None);
            };
            Ok(Some(torus_point(s, &frame, u, v)))
        }
        Surface::Sphere(s) => {
            let d = target.minus(&s.center);
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
                let denominator = &s.radius - &d.dot(&pole);
                match dec.sign(&denominator) {
                    Outcome::ProvenCoincident(_) | Outcome::NotExact(_) => continue,
                    Outcome::Refused(e) => return Err(e),
                    Outcome::Proven(_) => {}
                }
                let x = dec.enclose(&(&d.dot(&u) / &denominator), width)?.mid();
                let y = dec.enclose(&(&d.dot(&v) / &denominator), width)?.mid();
                let norm = &x * &x + &y * &y;
                let factor = &s.radius / &Real::from(Rational::one() + norm.clone());
                let radial = u
                    .scaled(&Real::from(r(2) * x))
                    .plus(&v.scaled(&Real::from(r(2) * y)))
                    .plus(&pole.scaled(&Real::from(norm - Rational::one())));
                return Ok(Some(s.center.plus(&radial.scaled(&factor))));
            }
            Ok(None)
        }
        Surface::EllipticCylinder(s) => {
            let first = s
                .major_dir
                .scaled(&(&Real::one() / &root(&s.major_dir.norm_squared())));
            let axis = s
                .axis
                .scaled(&(&Real::one() / &root(&s.axis.norm_squared())));
            let second = axis.cross(&first);
            let major = root(&s.major_sq);
            let minor = root(&s.minor_sq);
            let d = target.minus(&s.origin);
            let height = d.dot(&s.axis) / s.axis.norm_squared();
            for sense in [1, -1] {
                let cosine = &(&d.dot(&first) / &major) * &Real::from_i64(sense);
                let sine = &(&d.dot(&second) / &minor) * &Real::from_i64(sense);
                let denominator = &Real::one() + &cosine;
                if crate::geometry::sign(&denominator, dec, "boolean.imprint-ambiguous")?
                    == Sign::Zero
                {
                    continue;
                }
                let t = Real::from(dec.enclose(&(&sine / &denominator), width)?.mid());
                let height = Real::from(dec.enclose(&height, width)?.mid());
                let denominator = &Real::one() + &t.square();
                let cosine = (&Real::one() - &t.square()) / denominator.clone();
                let sine = (Real::from_i64(2) * t) / denominator;
                let radial = first
                    .scaled(&(&major * &cosine))
                    .plus(&second.scaled(&(&minor * &sine)))
                    .scaled(&Real::from_i64(sense));
                return Ok(Some(s.origin.plus(&radial).plus(&s.axis.scaled(&height))));
            }
            Ok(None)
        }
        _ => Ok(None),
    }
}

pub(crate) fn comparison_witness(
    surface: &Surface,
    target: &Point3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<RPoint3>, Refusal> {
    match surface {
        Surface::Cone(s) => {
            let (Some(apex), Some(axis), Some(reference), Some(slope)) = (
                s.apex.as_rational(),
                s.axis.as_rational(),
                s.ref_dir.as_rational(),
                s.slope.as_rational(),
            ) else {
                return Ok(None);
            };
            let Some(axis_length) = sqrt_rational(&axis.norm_squared()) else {
                return Ok(None);
            };
            let axis = axis.scaled(&axis_length.recip().unwrap());
            let Some(radial) = circle_direction(
                &axis,
                &reference,
                &target.minus(&apex.to_real()),
                width,
                dec,
            )?
            else {
                return Ok(None);
            };
            let height = target.minus(&apex.to_real()).dot(&axis.to_real());
            let height = dec.enclose(&height, width)?.mid().max(Rational::zero());
            Ok(Some(
                apex.plus(&axis.scaled(&height))
                    .plus(&radial.scaled(&(&slope * &height))),
            ))
        }
        Surface::Torus(s) => {
            let (Some(center), Some(axis), Some(reference), Some(major), Some(minor)) = (
                s.center.as_rational(),
                s.axis.as_rational(),
                s.ref_dir.as_rational(),
                s.major.as_rational(),
                s.minor.as_rational(),
            ) else {
                return Ok(None);
            };
            let Some(axis_length) = sqrt_rational(&axis.norm_squared()) else {
                return Ok(None);
            };
            let axis = axis.scaled(&axis_length.recip().unwrap());
            let d = target.minus(&center.to_real());
            let Some(radial) = circle_direction(&axis, &reference, &d, width, dec)? else {
                return Ok(None);
            };
            let height = d.dot(&axis.to_real());
            let planar = d.minus(&axis.to_real().scaled(&height));
            let x = planar.norm_squared().sqrt() - Real::from(major.clone());
            let tube = catallax_shapes::Vector3::new(x, height, Real::zero());
            let z = catallax_shapes::RVector3::new(r(0), r(0), r(1));
            let x = catallax_shapes::RVector3::new(r(1), r(0), r(0));
            let Some(tube) = circle_direction(&z, &x, &tube, width, dec)? else {
                return Ok(None);
            };
            Ok(Some(
                center
                    .plus(&radial.scaled(&(&major + &(&minor * &tube.x))))
                    .plus(&axis.scaled(&(&minor * &tube.y))),
            ))
        }
        _ => Ok(exact_probe(surface, target, width, dec)?.and_then(|point| point.as_rational())),
    }
}

fn circle_direction(
    axis: &catallax_shapes::RVector3,
    reference: &catallax_shapes::RVector3,
    target: &catallax_shapes::Vector3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<catallax_shapes::RVector3>, Refusal> {
    let Some(a) = sqrt_rational(&axis.norm_squared()) else {
        return Ok(None);
    };
    let Some(b) = sqrt_rational(&reference.norm_squared()) else {
        return Ok(None);
    };
    let axis = axis.scaled(&a.recip().unwrap());
    let u = reference.scaled(&b.recip().unwrap());
    let v = axis.cross(&u);
    let x = target.dot(&u.to_real());
    let y = target.dot(&v.to_real());
    let length = (&x.square() + &y.square()).sqrt();
    for sign in [1, -1] {
        let denominator = &length + &(&x * &Real::from_i64(sign));
        if !matches!(dec.sign(&denominator), Outcome::Proven(_)) {
            continue;
        }
        let t = dec.enclose(&(&y / &denominator), width)?.mid();
        let squared = &t * &t;
        let denominator = &Rational::one() + &squared;
        let numerator = Rational::one() - squared;
        let cosine = (&numerator / &denominator) * r(sign);
        let sine = r(2) * t / denominator;
        return Ok(Some(u.scaled(&cosine).plus(&v.scaled(&sine))));
    }
    Ok(None)
}

fn size_upper(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.clone());
    }
    Ok(dec.enclose(value, &Rational::one())?.hi())
}

fn search<T>(
    model: &mut Model,
    face: FaceId,
    dec: &mut Decider,
    barriers: &[catallax_topology::EdgeId],
    mut accept: impl FnMut(&mut Model, PointId, &Point3, &mut Decider) -> Result<Option<T>, Refusal>,
) -> Result<T, Refusal> {
    let f = model.face(face).clone();
    let surface = model.geom().surface(f.surface).clone();
    if f.loops.is_empty() {
        if let Surface::Torus(torus) = &surface {
            let frame = analytic_frame(&torus.axis, &torus.ref_dir).ok_or_else(|| {
                refusal(
                    "boolean.no-rational-point",
                    "torus frame requires algebraic face probes from the shared foundation",
                )
            })?;
            for k in 0..crate::geometry::DIRECTION_ATTEMPTS {
                let k = i64::from(k);
                let p = torus_point(
                    torus,
                    &frame,
                    rational_circle(&r(k)),
                    rational_circle(&r(k * k + 1)),
                );
                let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
                if let Some(value) = accept(model, point, &p, dec)? {
                    return Ok(value);
                }
            }
            return Err(refusal(
                "boolean.no-rational-point",
                "no whole-torus probe settled the query",
            ));
        }
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
                Surface::Cylinder(s) => size_upper(&s.radius(), dec)?,
                Surface::Sphere(s) => size_upper(&s.radius, dec)?,
                Surface::EllipticCylinder(s) => size_upper(&root(&s.minor_sq), dec)?,
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
                let Some(candidate) = exact_probe(&surface, &target, &width, dec)? else {
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
    accept: impl FnMut(&mut Model, PointId, &Point3, &mut Decider) -> Result<Option<T>, Refusal>,
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
    accept: impl FnMut(&mut Model, PointId, &Point3, &mut Decider) -> Result<Option<T>, Refusal>,
) -> Result<T, Refusal> {
    search(model, face, dec, barriers, accept)
}

pub(crate) fn face_point(model: &mut Model, face: FaceId, dec: &mut Decider) -> Answered<Point3> {
    find(model, face, dec, |_, _, p, _| Ok(Some(p.clone()))).map(Answer::exact)
}
