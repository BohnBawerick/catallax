//! Deterministic per-property refinement and orientation handling.

use catallax_intersect::{CurveTrace, trace};
use catallax_numbers::{Answer, Answered, Ball, Decider, Outcome, Rational, Real, Refusal};
use catallax_shapes::{Circle, Curve, Ellipse, Line, Point3, RPoint3, RVector3, Surface, Vector3};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};
use std::collections::{BTreeSet, HashMap, hash_map::Entry};

use crate::integrate::{Stop, Work};
use crate::planar;
use crate::surface::{self, Compiled, Property};
use crate::{Accuracy, Enclosure};

#[derive(Default)]
pub(crate) struct TraceCache {
    curves: HashMap<EdgeId, std::sync::Arc<CurveTrace>>,
    failed_gauges: BTreeSet<GaugeAttempt>,
    reuse_gauges: bool,
    torus_rectangles: HashMap<FaceId, Option<std::sync::Arc<crate::torus::rectangle::Domain>>>,
}

impl TraceCache {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn grouped() -> Self {
        Self {
            reuse_gauges: true,
            ..Self::default()
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct GaugeAttempt {
    face: FaceId,
    depth: u32,
    rectangle: bool,
    precision: u32,
    candidate: usize,
}

struct Arc<'a> {
    trace: Option<std::sync::Arc<CurveTrace>>,
    reversed: bool,
    analytic: Option<crate::integrate::Program>,
    integral: Option<crate::integrate::Analytic>,
    exact: Option<[Rational; 5]>,
    meridian: Option<(&'a Point3, i64)>,
    parallel: Option<Parallel>,
    /// A meridian arc of an outer torus sheet: `du` vanishes along it, so it
    /// only advances the lifted tube angle to its exact endpoint value.
    still: bool,
}
struct Parallel {
    cosine: Real,
    sine: Real,
    orientation: i64,
}

struct Curved<'a> {
    id: FaceId,
    surface: &'a Surface,
    sense: Sense,
    boundaries: Vec<Vec<Arc<'a>>>,
    compiled: Vec<Compiled>,
    torus: Option<crate::torus::Programs>,
    rectangle: Option<std::sync::Arc<crate::torus::rectangle::Domain>>,
    analytic_only: bool,
}

enum Prepared<'a> {
    Curved(Curved<'a>),
    PlanarArea {
        area_numerator: Rational,
        normal_squared: Rational,
    },
}

pub(crate) fn anchor(model: &Model, faces: &[FaceId]) -> RPoint3 {
    let point = match model.geom().surface(model.face(faces[0]).surface) {
        Surface::Plane(p) => &p.origin,
        Surface::Cylinder(c) => &c.origin,
        Surface::Cone(c) => &c.apex,
        Surface::Sphere(s) => &s.center,
        Surface::Torus(t) => &t.center,
        Surface::TorusSheet(t) => &t.center,
        Surface::EllipticCylinder(c) => &c.origin,
        _ => return RPoint3::origin(),
    };
    point.as_rational().unwrap_or_else(RPoint3::origin)
}

fn budget() -> Refusal {
    Refusal::new(
        "measure.budget-exhausted",
        "no enclosure reached within the measurement budget",
    )
}

fn finite_distance_face(model: &Model, face: FaceId) -> bool {
    let boundary = model.face(face);
    if matches!(model.geom().surface(boundary.surface), Surface::Plane(_)) {
        return false;
    }
    distance_generated(model, face)
        && boundary.loops.iter().any(|lp| {
            matches!(model.loop_(*lp), Loop::Coedges(cs) if cs.iter().any(|c| {
                let edge = model.edge(model.coedge(*c).edge);
                edge.start.is_some() || edge.end.is_some()
            }))
        })
}

fn distance_generated(model: &Model, face: FaceId) -> bool {
    model.face_origins(face).iter().any(|origin| {
        let mut role = &origin.role;
        while let catallax_recipe::Role::Copy { of, .. } = role {
            role = &of.role;
        }
        matches!(
            role,
            catallax_recipe::Role::Round { .. }
                | catallax_recipe::Role::RoundCorner { .. }
                | catallax_recipe::Role::Inner { .. }
                | catallax_recipe::Role::Offset { .. }
        )
    })
}

fn analytic_carrier(
    model: &Model,
    curve: catallax_shapes::CurveId,
    dec: &mut Decider,
) -> Result<Option<Curve>, Refusal> {
    if !matches!(model.geom().curve(curve), Curve::Intersection(_)) {
        return Ok(None);
    }
    Ok(Some(
        match catallax_intersect::curve_form(model.geom(), curve, dec)? {
            catallax_intersect::CurveForm::Line { origin, dir } => {
                Curve::Line(Line { origin, dir })
            }
            catallax_intersect::CurveForm::Circle {
                center,
                normal,
                radius_sq,
                ref_dir,
            } => Curve::Circle(Circle {
                center,
                normal,
                radius_sq,
                ref_dir,
            }),
            catallax_intersect::CurveForm::Ellipse {
                center,
                normal,
                major_dir,
                major_sq,
                minor_sq,
            } => Curve::Ellipse(Ellipse {
                center,
                normal,
                major_dir,
                major_sq,
                minor_sq,
            }),
            _ => return Ok(None),
        },
    ))
}

pub(crate) fn measure(
    model: &Model,
    faces: &[FaceId],
    anchor: &RPoint3,
    accuracy: &Accuracy,
    property: Property,
    traces: &mut TraceCache,
    dec: &mut Decider,
) -> Answered<Vec<Enclosure>> {
    let indices = match property {
        Property::Area => vec![0],
        Property::Volume => vec![1],
        Property::Centroid | Property::FaceCentroid => vec![1, 2, 3, 4],
        Property::SecondMoment(_, _) => vec![2],
    };
    if let Some(values) = planar::body(model, faces, anchor, property) {
        let result = indices
            .iter()
            .map(|i| Enclosure::exact(values[*i].clone()))
            .collect();
        return finish(result, anchor, property);
    }
    let mut exact: [Rational; 5] = std::array::from_fn(|_| Rational::zero());
    let mut pending = Vec::new();
    let reduce_distance = faces.iter().any(|face| distance_generated(model, *face));
    for id in faces {
        if let Some(planar) = planar::face(model, *id, anchor, property) {
            match property {
                Property::Area => {
                    if let Some(length) = planar::sqrt(&planar.normal_squared) {
                        exact[0] = &exact[0] + &(planar.area_numerator / length);
                    } else {
                        pending.push((*id, Some(planar)));
                    }
                }
                Property::Volume | Property::Centroid => {
                    exact[1] = &exact[1] + &planar.volume;
                    if property == Property::Centroid {
                        for (sum, moment) in exact[2..].iter_mut().zip(planar.moments) {
                            *sum = &*sum + &moment;
                        }
                    }
                }
                Property::SecondMoment(_, _) => {
                    exact[2] = &exact[2] + &planar.moments[0];
                }
                Property::FaceCentroid => {
                    exact[1] = &exact[1] + &planar.area_numerator;
                    let sense = planar::sign(model.face(*id).sense);
                    for (sum, moment) in exact[2..].iter_mut().zip(planar.moments) {
                        *sum = &*sum + &(&moment * &sense);
                    }
                }
            }
        } else {
            pending.push((*id, None));
        }
    }
    if pending.is_empty() {
        let result = indices
            .iter()
            .map(|i| Enclosure::exact(exact[*i].clone()))
            .collect::<Vec<_>>();
        return finish(result, anchor, property);
    }
    let policy = *dec.policy();
    if policy.max_precision_bits == 0 || policy.max_work == 0 {
        return Err(budget());
    }
    // An outer torus sheet uses the ring torus's (u, v) chart and forms. Its
    // trimmed face stays where `major + minor cos v > 0`, so boundary
    // coordinates `cos v = (rho - major)/minor` are those of its own sheet.
    let sheets: HashMap<FaceId, Surface> = pending
        .iter()
        .filter_map(
            |(id, _)| match model.geom().surface(model.face(*id).surface) {
                Surface::TorusSheet(s) if !model.face(*id).loops.is_empty() => Some((
                    *id,
                    Surface::Torus(catallax_shapes::Torus {
                        center: s.center.clone(),
                        axis: s.axis.clone(),
                        major: s.major.clone(),
                        minor: s.minor.clone(),
                        ref_dir: s.ref_dir.clone(),
                    }),
                )),
                _ => None,
            },
        )
        .collect();
    let mut prepared = Vec::new();
    for (id, planar) in pending {
        if let Some(planar) = planar {
            prepared.push(Prepared::PlanarArea {
                area_numerator: planar.area_numerator,
                normal_squared: planar.normal_squared,
            });
            continue;
        }
        let face = model.face(id);
        let surface = model.geom().surface(face.surface);
        let closed = face.loops.is_empty();
        let surface = sheets.get(&id).unwrap_or(surface);
        let analytic_only = finite_distance_face(model, id);
        let rectangle = match surface {
            Surface::Torus(torus) if analytic_only => {
                let domain = match traces.torus_rectangles.entry(id) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => entry.insert(
                        crate::torus::rectangle::prepare(model, id, torus, dec)?
                            .map(std::sync::Arc::new),
                    ),
                };
                domain.clone()
            }
            _ => None,
        };
        let torus = match surface {
            Surface::Torus(t) if !closed => Some(crate::torus::compile(
                t,
                anchor,
                property,
                rectangle.is_some(),
            )),
            _ => None,
        };
        let compiled = surface::compile(
            surface,
            anchor,
            property,
            closed || torus.is_some(),
            &gauge(surface, 0),
        )
        .ok_or_else(|| {
            Refusal::new(
                "measure.unsupported-surface",
                "surface kind has no measurement definition",
            )
        })?;
        let sheet = sheets.contains_key(&id);
        let mut boundaries = Vec::new();
        let loops = if rectangle.is_some() {
            &[][..]
        } else {
            &face.loops
        };
        for id in loops {
            let mut arcs = Vec::new();
            match model.loop_(*id) {
                Loop::Vertex(v) => {
                    let width = if accuracy.absolute > Rational::zero() {
                        accuracy.absolute.clone()
                    } else {
                        Rational::one().mul_pow2(-64)
                    };
                    model
                        .geom()
                        .enclose_point(model.vertex(*v).point, &width, dec)?;
                }
                Loop::Coedges(coedges) => {
                    for c in coedges {
                        let c = model.coedge(*c);
                        let edge = model.edge(c.edge);
                        let mut from = edge.start.map(|v| model.vertex(v).point);
                        let mut to = edge.end.map(|v| model.vertex(v).point);
                        if edge.sense == Sense::Reversed {
                            std::mem::swap(&mut from, &mut to);
                        }
                        let exact = if matches!(
                            model.geom().curve(edge.curve),
                            catallax_shapes::Curve::Line(_)
                        ) {
                            planar::line(surface, model.geom(), from, to, anchor, property)
                        } else {
                            None
                        };
                        let exact = if exact.is_none() {
                            cone_generator(surface, model, edge.curve, dec)?
                        } else {
                            exact
                        };
                        let circle = match model.geom().curve(edge.curve) {
                            catallax_shapes::Curve::Circle(circle) if edge.start == edge.end => {
                                Some(circle)
                            }
                            _ => None,
                        };
                        // A whole analytic circle has no endpoint-location
                        // question. Its carrier integral is evaluated directly.
                        let analytic = if let Some(circle) = circle {
                            crate::ring::compile(
                                surface,
                                circle,
                                anchor,
                                &compiled.indices,
                                property,
                                dec,
                            )?
                        } else {
                            None
                        };
                        let meridian = if coedges.len() == 1 {
                            circle.and_then(|circle| torus_meridian(surface, circle))
                        } else {
                            None
                        };
                        let parallel =
                            if torus.is_some() && coedges.len() == 1 && edge.start == edge.end {
                                torus_parallel(surface, model, edge.curve, dec)?
                            } else {
                                None
                            };
                        let mut integral =
                            if torus.is_none() && exact.is_none() && analytic.is_none() {
                                from.zip(to).and_then(|(from, to)| {
                                    match (model.geom().point(from), model.geom().point(to)) {
                                        (
                                            catallax_shapes::PointDef::Explicit(from),
                                            catallax_shapes::PointDef::Explicit(to),
                                        ) => from.as_rational().zip(to.as_rational()).and_then(
                                            |(from, to)| {
                                                compiled.boundary.analytic(
                                                    model.geom().curve(edge.curve),
                                                    &from,
                                                    &to,
                                                )
                                            },
                                        ),
                                        _ => None,
                                    }
                                })
                            } else {
                                None
                            };
                        if reduce_distance
                            && torus.is_none()
                            && exact.is_none()
                            && analytic.is_none()
                            && integral.is_none()
                            && let Some((from, to)) = from.zip(to)
                        {
                            let carrier = analytic_carrier(model, edge.curve, dec)?;
                            integral = compiled.boundary.analytic_reduced(
                                carrier.as_ref().unwrap_or(model.geom().curve(edge.curve)),
                                &model.geom().position(from),
                                &model.geom().position(to),
                                dec,
                            )?;
                        }
                        // Admit a finite distance trim only with a complete
                        // antiderivative. Incidence alone proves no integral.
                        if analytic_only
                            && exact.is_none()
                            && analytic.is_none()
                            && integral.is_none()
                            && meridian.is_none()
                            && parallel.is_none()
                        {
                            return Err(Refusal::new(
                                "measure.budget-exhausted",
                                "finite curved distance-feature trim has no proven analytic property integral",
                            ));
                        }
                        let trace = if exact.is_none()
                            && analytic.is_none()
                            && meridian.is_none()
                            && parallel.is_none()
                            && (integral.is_none()
                                || (matches!(surface, Surface::Sphere(_)) && !analytic_only))
                        {
                            let traced = match traces.curves.entry(c.edge) {
                                Entry::Occupied(entry) => entry.into_mut(),
                                Entry::Vacant(entry) => entry.insert(std::sync::Arc::new(
                                    trace(model.geom(), edge.curve, from, to, dec)?.value,
                                )),
                            };
                            Some(traced.clone())
                        } else {
                            None
                        };
                        let still =
                            sheet && sheet_meridian(surface, model.geom().curve(edge.curve), dec)?;
                        arcs.push(Arc {
                            trace,
                            reversed: (edge.sense == Sense::Reversed)
                                != (c.sense == Sense::Reversed),
                            analytic,
                            integral,
                            exact,
                            meridian,
                            parallel,
                            still,
                        });
                    }
                }
            }
            if !arcs.is_empty() {
                boundaries.push(arcs);
            }
        }
        prepared.push(Prepared::Curved(Curved {
            id,
            surface,
            sense: face.sense,
            boundaries,
            compiled: vec![compiled],
            torus,
            rectangle,
            analytic_only,
        }));
    }
    let mut best = None;
    // Six independent scalar refinements share the tensor's arithmetic allowance.
    let mut remaining = if matches!(property, Property::SecondMoment(_, _)) {
        policy.max_work / 6
    } else {
        policy.max_work
    };
    let mut depth = 0u32;
    let mut rectangle = true;
    let mut working_precision = 128u32.min(policy.max_precision_bits);
    let mut split_precision = false;
    let mut volume_enclosure_reached = false;
    loop {
        dec.check_deadline()?;
        let precision = if rectangle {
            64u32.min(policy.max_precision_bits)
        } else {
            working_precision
        };
        let mut work = match Work::new(&policy, precision, remaining, dec) {
            Ok(work) => work,
            Err(Stop::Refusal(refusal)) => return Err(refusal),
            Err(_) => break,
        };
        let attempt = evaluate(
            &mut prepared,
            anchor,
            property,
            &exact,
            &indices,
            depth,
            rectangle,
            traces,
            &mut work,
        );
        remaining = work.remaining();
        dec.check_deadline()?;
        match attempt {
            Ok(result) => {
                volume_enclosure_reached = true;
                let result = finish(result, anchor, property);
                match result {
                    Ok(answer) => {
                        // Standalone volume must retain the separation needed
                        // by mass_properties, even below the absolute floor.
                        let separated = property != Property::Volume
                            || !answer.value[0].contains(&Rational::zero());
                        let reached = separated && answer.value.iter().all(|v| accuracy.reached(v));
                        best = Some(answer);
                        if reached {
                            return Ok(best.unwrap());
                        }
                    }
                    Err(r) if r.code == "measure.centroid-not-separated" => {}
                    Err(r) => return Err(r),
                }
                if rectangle {
                    rectangle = false;
                    depth = 0;
                } else if work.quadrature_dominates {
                    depth += 1;
                } else if working_precision < policy.max_precision_bits {
                    working_precision = working_precision
                        .saturating_mul(2)
                        .min(policy.max_precision_bits);
                } else {
                    break;
                }
            }
            Err(Stop::Budget) => break,
            Err(Stop::Refusal(refusal)) => return Err(refusal),
            Err(Stop::Split | Stop::GaugeSplit) => {
                if rectangle && depth < 2 {
                    depth += 1;
                } else if rectangle {
                    rectangle = false;
                    depth = 0;
                } else if !split_precision && working_precision < policy.max_precision_bits {
                    working_precision = working_precision
                        .saturating_mul(2)
                        .min(policy.max_precision_bits);
                    split_precision = true;
                } else {
                    depth += 1;
                    split_precision = false;
                }
            }
        }
        if remaining == 0 {
            break;
        }
    }
    best.ok_or_else(|| {
        if matches!(property, Property::Centroid | Property::FaceCentroid)
            && volume_enclosure_reached
        {
            Refusal::new(
                "measure.centroid-not-separated",
                if property == Property::FaceCentroid {
                    "area could not be separated from zero"
                } else {
                    "volume could not be separated from zero"
                },
            )
        } else {
            budget()
        }
    })
}

fn finish(
    result: Vec<Enclosure>,
    anchor: &RPoint3,
    property: Property,
) -> Answered<Vec<Enclosure>> {
    if !matches!(property, Property::Centroid | Property::FaceCentroid) {
        return Ok(Answer::exact(result));
    }
    let volume = &result[0];
    if volume.lo <= Rational::zero() && volume.hi >= Rational::zero() {
        return Err(Refusal::new(
            "measure.centroid-not-separated",
            if property == Property::FaceCentroid {
                "area enclosure includes zero"
            } else {
                "volume enclosure includes zero"
            },
        ));
    }
    let coordinates = planar::components(anchor);
    let result = result[1..]
        .iter()
        .zip(coordinates)
        .map(|(moment, offset)| {
            let quotients = [
                &moment.lo / &volume.lo,
                &moment.lo / &volume.hi,
                &moment.hi / &volume.lo,
                &moment.hi / &volume.hi,
            ];
            Enclosure {
                lo: quotients.iter().min().unwrap() + offset,
                hi: quotients.iter().max().unwrap() + offset,
            }
        })
        .collect();
    Ok(Answer::exact(result))
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    prepared: &mut [Prepared<'_>],
    anchor: &RPoint3,
    property: Property,
    exact: &[Rational; 5],
    indices: &[usize],
    depth: u32,
    rectangle: bool,
    cache: &mut TraceCache,
    work: &mut Work,
) -> Result<Vec<Enclosure>, Stop> {
    let mut sums: Vec<Ball> = indices
        .iter()
        .map(|i| work.rational(&exact[*i]))
        .collect::<Result<_, _>>()?;
    for face in prepared {
        let values = match face {
            Prepared::PlanarArea {
                area_numerator,
                normal_squared,
            } => {
                let numerator = work.rational(area_numerator)?;
                let square = work.rational(normal_squared)?;
                let norm = work.sqrt(&square)?;
                vec![work.div(&numerator, &norm)?]
            }
            Prepared::Curved(face) => {
                let mut answer = None;
                let candidates = if matches!(face.surface, Surface::Sphere(_))
                    && !face.boundaries.is_empty()
                    && !face.analytic_only
                {
                    6
                } else {
                    1
                };
                let cacheable = cache.reuse_gauges
                    && candidates == 6
                    && !matches!(property, Property::SecondMoment(_, _))
                    && face.boundaries.iter().flatten().all(|arc| {
                        arc.integral.is_none() && arc.exact.is_none() && arc.analytic.is_none()
                    });
                let face_id = face.id;
                let precision = work.precision;
                let key = |candidate| GaugeAttempt {
                    face: face_id,
                    depth,
                    rectangle,
                    precision,
                    candidate,
                };
                // The area-gauge denominator and its trace zeroth coefficient
                // are identical across these properties. Reuse only a complete
                // set of failed gauges, never a geometric nonmembership claim.
                // Skipping individual gauges in a successful attempt could
                // change its quadrature-dominance flag.
                if cacheable
                    && (0..candidates)
                        .all(|candidate| cache.failed_gauges.contains(&key(candidate)))
                {
                    return Err(Stop::Split);
                }
                for candidate in 0..candidates {
                    if face.compiled.len() <= candidate {
                        face.compiled.push(
                            surface::compile(
                                face.surface,
                                anchor,
                                property,
                                false,
                                &gauge(face.surface, candidate),
                            )
                            .unwrap(),
                        );
                    }
                    match integrate_face(
                        face,
                        &face.compiled[candidate],
                        indices,
                        depth,
                        rectangle,
                        work,
                    ) {
                        Ok(values) => {
                            answer = Some(values);
                            break;
                        }
                        Err(Stop::GaugeSplit) => {
                            if cacheable && property != Property::Centroid {
                                cache.failed_gauges.insert(key(candidate));
                            }
                        }
                        Err(Stop::Split) => {}
                        Err(e) => return Err(e),
                    }
                }
                answer.ok_or(Stop::Split)?
            }
        };
        for (sum, value) in sums.iter_mut().zip(values) {
            *sum = work.add(sum, &value)?;
        }
    }
    Ok(sums.iter().map(Enclosure::from_ball).collect())
}

fn integrate_face(
    face: &Curved<'_>,
    compiled: &Compiled,
    indices: &[usize],
    depth: u32,
    rectangle: bool,
    work: &mut Work,
) -> Result<Vec<Ball>, Stop> {
    if let Some(torus) = &face.torus {
        return match torus {
            crate::torus::Programs::Boundary(programs) => {
                integrate_torus(face, programs, indices, depth, rectangle, work)
            }
            crate::torus::Programs::Rectangle(programs) => {
                let domain = face.rectangle.as_ref().unwrap();
                let extras = domain.extras(work)?;
                let zero = work.rational(&Rational::zero())?;
                let mut sums =
                    programs
                        .integral
                        .point(&[zero.clone(), zero.clone(), zero], &extras, work)?;
                let orientation = work.rational(&Rational::from_i64(domain.orientation))?;
                for sum in &mut sums {
                    *sum = work.mul(sum, &orientation)?;
                }
                finish_torus(face, sums, &programs.full, &programs.indices, indices, work)
            }
        };
    }
    let zero = work.rational(&Rational::zero())?;
    let sense = work.rational(&planar::sign(face.sense))?;
    let mut sums = vec![zero.clone(); compiled.indices.len()];
    let pi = if compiled.full.is_some()
        || face
            .boundaries
            .iter()
            .flatten()
            .any(|arc| arc.analytic.is_some())
    {
        Some(work.pi()?)
    } else {
        None
    };
    let full = if let Some(full) = &compiled.full {
        Some(full.point(
            &[zero.clone(), zero.clone(), zero.clone()],
            &[pi.as_ref().unwrap().clone()],
            work,
        )?)
    } else {
        None
    };
    if compiled.closed {
        sums = full.unwrap();
    } else {
        let count = 1u64.checked_shl(depth).ok_or(Stop::Budget)?;
        for boundary in &face.boundaries {
            for arc in boundary {
                if let Some(exact) = &arc.exact {
                    for (sum, index) in sums.iter_mut().zip(&compiled.indices) {
                        let value = work.rational(&exact[*index])?;
                        *sum = if arc.reversed {
                            work.sub(sum, &value)?
                        } else {
                            work.add(sum, &value)?
                        };
                    }
                    continue;
                }
                // Each alternate sphere gauge has a different boundary form.
                if std::ptr::eq(compiled, &face.compiled[0])
                    && let Some(integral) = &arc.integral
                {
                    for (sum, value) in sums.iter_mut().zip(integral.evaluate(work)?) {
                        *sum = if arc.reversed {
                            work.sub(sum, &value)?
                        } else {
                            work.add(sum, &value)?
                        };
                    }
                    continue;
                }
                if let Some(analytic) = &arc.analytic {
                    let values = analytic.point(
                        &[zero.clone(), zero.clone(), zero.clone()],
                        &[pi.as_ref().unwrap().clone()],
                        work,
                    )?;
                    for (sum, value) in sums.iter_mut().zip(values) {
                        *sum = if arc.reversed {
                            work.sub(sum, &value)?
                        } else {
                            work.add(sum, &value)?
                        };
                    }
                    continue;
                }
                for piece in arc.trace.as_ref().unwrap().pieces() {
                    let (lo, hi) = piece.range();
                    let step = (&hi - &lo)
                        / Rational::from_integer(catallax_numbers::Integer::from_u64(count));
                    for i in 0..count {
                        let a = &lo
                            + &(&step
                                * &Rational::from_integer(catallax_numbers::Integer::from_u64(i)));
                        let b = &a + &step;
                        let values = if rectangle {
                            compiled.boundary.rectangle(piece, &a, &b, work)?
                        } else {
                            compiled.boundary.panel(piece, &a, &b, work)?
                        };
                        for (sum, value) in sums.iter_mut().zip(values) {
                            *sum = if arc.reversed {
                                work.sub(sum, &value)?
                            } else {
                                work.add(sum, &value)?
                            };
                        }
                    }
                }
            }
        }
        // Integrals above followed the coedges. Convert to the surface's
        // orientation before selecting the bounded region of a closed surface.
        for sum in &mut sums {
            *sum = work.mul(sum, &sense)?;
        }
        if compiled.sphere {
            let full = full.as_ref().unwrap();
            let area_index = compiled.indices.iter().position(|i| *i == 0).unwrap();
            correct_whole_surface(&mut sums, full, area_index, work)?;
        }
    }
    indices
        .iter()
        .map(|i| {
            let value = &sums[compiled.indices.iter().position(|j| i == j).unwrap()];
            if *i == 0 {
                Ok(value.clone())
            } else {
                work.mul(value, &sense)
            }
        })
        .collect()
}

fn gauge(surface: &Surface, index: usize) -> Vector3 {
    if let Surface::Sphere(sphere) = surface {
        let vector = match index % 3 {
            0 => sphere.axis.clone(),
            1 => sphere.ref_dir.clone(),
            _ => sphere.axis.cross(&sphere.ref_dir),
        };
        return if index < 3 {
            vector
        } else {
            Vector3::new(-vector.x, -vector.y, -vector.z)
        };
    }
    let axes = [
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, -1],
        [-1, 0, 0],
        [0, -1, 0],
    ];
    let [x, y, z] = axes[index];
    RVector3::new(
        Rational::from_i64(x),
        Rational::from_i64(y),
        Rational::from_i64(z),
    )
    .to_real()
}

fn integrate_torus(
    face: &Curved<'_>,
    programs: &crate::torus::Boundary,
    indices: &[usize],
    depth: u32,
    rectangle: bool,
    work: &mut Work,
) -> Result<Vec<Ball>, Stop> {
    let zero = work.rational(&Rational::zero())?;
    let mut sums = vec![zero.clone(); programs.indices.len()];
    let count = 1u64.checked_shl(depth).ok_or(Stop::Budget)?;
    for boundary in &face.boundaries {
        let first = &boundary[0];
        if let Some(parallel) = &first.parallel {
            // A full parallel has constant tube angle and makes one exact
            // azimuthal turn. Integrate its Fourier mean, not a traced curve.
            let cosine = work.real(&parallel.cosine)?;
            let sine = work.real(&parallel.sine)?;
            let angle = crate::integrate::angle(&cosine, &sine, work)?;
            let pi = work.pi()?;
            let values = programs.parallel.point(
                &[zero.clone(), zero.clone(), zero.clone()],
                &[angle, cosine, sine, pi],
                work,
            )?;
            let orientation = if first.reversed {
                -parallel.orientation
            } else {
                parallel.orientation
            };
            let orientation = work.rational(&Rational::from_i64(orientation))?;
            for (sum, value) in sums.iter_mut().zip(values) {
                let value = work.mul(&value, &orientation)?;
                *sum = work.add(sum, &value)?;
            }
            continue;
        }
        if let Some((core, orientation)) = first.meridian {
            // Along a meridian du is identically zero. Its lifted tube angle
            // makes exactly one turn, so only the holonomy term remains.
            let position = [
                work.real(&core.x)?,
                work.real(&core.y)?,
                work.real(&core.z)?,
            ];
            let coordinates = programs.coordinates.point(&position, &[], work)?;
            let u = crate::integrate::angle(&coordinates[0], &coordinates[1], work)?;
            let correction = programs.correction.point(&position, &[u], work)?;
            let pi = work.pi()?;
            let turns = work.rational(&Rational::from_i64(if first.reversed {
                -2 * orientation
            } else {
                2 * orientation
            }))?;
            let winding = work.mul(&pi, &turns)?;
            for (sum, correction) in sums.iter_mut().zip(correction) {
                let correction = work.mul(&winding, &correction)?;
                *sum = work.add(sum, &correction)?;
            }
            continue;
        }
        let first_trace = first.trace.as_ref().unwrap();
        let piece = if first.reversed {
            first_trace.pieces().last().unwrap()
        } else {
            &first_trace.pieces()[0]
        };
        let (a, b) = piece.range();
        let argument = work.rational(&if first.reversed { b } else { a })?;
        let position = piece.evaluate(&argument, work.dec)?;
        let coordinates = programs.coordinates.point(&position, &[], work)?;
        let u = crate::integrate::angle(&coordinates[0], &coordinates[1], work)?;
        let initial_v = crate::integrate::angle(&coordinates[2], &coordinates[3], work)?;
        let correction = programs.correction.point(&position, &[u], work)?;
        let mut v = initial_v.clone();
        for arc in boundary {
            let pieces = arc.trace.as_ref().unwrap().pieces();
            if arc.still
                && let Some(turn) = sheet_turn(arc, pieces, programs, work)?
            {
                v = work.add(&v, &turn)?;
                continue;
            }
            for j in 0..pieces.len() {
                let piece = &pieces[if arc.reversed {
                    pieces.len() - 1 - j
                } else {
                    j
                }];
                let (lo, hi) = piece.range();
                let (lo, hi) = if arc.reversed { (hi, lo) } else { (lo, hi) };
                let step = (&hi - &lo)
                    / Rational::from_integer(catallax_numbers::Integer::from_u64(count));
                for i in 0..count {
                    let a = &lo
                        + &(&step
                            * &Rational::from_integer(catallax_numbers::Integer::from_u64(i)));
                    let b = &a + &step;
                    let (values, end) = programs.boundary.lifted_panel(
                        &programs.angular_derivative,
                        piece,
                        &a,
                        &b,
                        &v,
                        rectangle,
                        work,
                    )?;
                    for (sum, value) in sums.iter_mut().zip(values) {
                        *sum = work.add(sum, &value)?;
                    }
                    v = end;
                }
            }
        }
        let winding = work.sub(&v, &initial_v)?;
        for (sum, correction) in sums.iter_mut().zip(correction) {
            let correction = work.mul(&winding, &correction)?;
            *sum = work.add(sum, &correction)?;
        }
    }
    finish_torus(face, sums, &programs.full, &programs.indices, indices, work)
}

/// The tube-angle change along an outer-sheet meridian arc. Every point of
/// the sheet has `major + minor cos v > 0`, so its tube angle has a unique
/// representative in `(-pi, pi)`, and the arc stays inside that interval.
/// An enclosure not separated from the cut at `pi` falls back to quadrature.
fn sheet_turn(
    arc: &Arc<'_>,
    pieces: &[catallax_intersect::TracePiece],
    programs: &crate::torus::Boundary,
    work: &mut Work,
) -> Result<Option<Ball>, Stop> {
    let pi = work.pi()?;
    let two = work.rational(&Rational::from_i64(2))?;
    let turn = work.mul(&pi, &two)?;
    let mut ends = Vec::with_capacity(2);
    for last in [false, true] {
        let piece = if last != arc.reversed {
            pieces.last().unwrap()
        } else {
            &pieces[0]
        };
        let (lo, hi) = piece.range();
        let argument = work.rational(&if last != arc.reversed { hi } else { lo })?;
        let position = piece.evaluate(&argument, work.dec)?;
        let coordinates = programs.coordinates.point(&position, &[], work)?;
        // `angle` returns `(-pi/2, 3pi/2)`; fold `(pi, 3pi/2)` below the cut.
        let mut angle = crate::integrate::angle(&coordinates[2], &coordinates[3], work)?;
        if angle.lo() > pi.hi() {
            angle = work.sub(&angle, &turn)?;
        } else if angle.hi() >= pi.lo() {
            return Ok(None);
        }
        ends.push(angle);
    }
    Ok(Some(work.sub(&ends[1], &ends[0])?))
}

/// Whether a circle is a tube circle of an outer sheet's chart torus: radius
/// `minor`, centred on the core circle, in the meridian plane through it.
fn sheet_meridian(surface: &Surface, curve: &Curve, dec: &mut Decider) -> Result<bool, Refusal> {
    let (Surface::Torus(torus), Curve::Circle(circle)) = (surface, curve) else {
        return Ok(false);
    };
    let offset = circle.center.minus(&torus.center);
    let zero = |value: Real, dec: &mut Decider| -> Result<bool, Refusal> {
        match dec.sign(&value) {
            Outcome::ProvenCoincident(_) => Ok(true),
            Outcome::Proven(_) => Ok(false),
            Outcome::Refused(e) => Err(e),
            Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        }
    };
    Ok(zero(&circle.radius_sq - &torus.minor.square(), dec)?
        && zero(offset.dot(&torus.axis), dec)?
        && zero(offset.norm_squared() - torus.major.square(), dec)?
        && zero(
            circle
                .normal
                .cross(&offset.cross(&torus.axis))
                .norm_squared(),
            dec,
        )?)
}

fn finish_torus(
    face: &Curved<'_>,
    mut sums: Vec<Ball>,
    full: &crate::integrate::Program,
    source_indices: &[usize],
    indices: &[usize],
    work: &mut Work,
) -> Result<Vec<Ball>, Stop> {
    let sense = work.rational(&planar::sign(face.sense))?;
    for sum in &mut sums {
        *sum = work.mul(sum, &sense)?;
    }
    let zero = work.rational(&Rational::zero())?;
    let pi = work.pi()?;
    let full = full.point(&[zero.clone(), zero.clone(), zero], &[pi], work)?;
    correct_whole_surface(&mut sums, &full, 0, work)?;
    indices
        .iter()
        .map(|i| {
            let value = &sums[source_indices.iter().position(|j| i == j).unwrap()];
            if *i == 0 {
                Ok(value.clone())
            } else {
                work.mul(value, &sense)
            }
        })
        .collect()
}

fn correct_whole_surface(
    sums: &mut [Ball],
    full: &[Ball],
    area_index: usize,
    work: &mut Work,
) -> Result<(), Stop> {
    let ratio = work.div(&sums[area_index], &full[area_index])?;
    // 0 <= raw_area + k * full_area <= full_area. Both ends are closed:
    // an integer lift can describe either an empty region or a whole surface
    // cut along a seam. Never silently discard the whole-surface possibility.
    let shift_lo = -ratio.hi().floor();
    let shift_hi = catallax_numbers::Integer::one() - ratio.lo().ceil();
    let shift = work.interval(
        &Rational::from_integer(shift_lo),
        &Rational::from_integer(shift_hi),
    )?;
    for (sum, whole) in sums.iter_mut().zip(full) {
        let correction = work.mul(&shift, whole)?;
        *sum = work.add(sum, &correction)?;
    }
    Ok(())
}

fn torus_meridian<'a>(
    surface: &Surface,
    circle: &'a catallax_shapes::Circle,
) -> Option<(&'a Point3, i64)> {
    let Surface::Torus(torus) = surface else {
        return None;
    };
    let offset = circle.center.minus(&torus.center);
    let normal = offset.cross(&torus.axis);
    if !(&circle.radius_sq - &torus.minor.square())
        .as_rational()?
        .is_zero()
        || !(offset.norm_squared() - torus.major.square())
            .as_rational()?
            .is_zero()
        || !offset.dot(&torus.axis).as_rational()?.is_zero()
        || !circle
            .normal
            .cross(&normal)
            .norm_squared()
            .as_rational()?
            .is_zero()
    {
        return None;
    }
    let orientation =
        if circle.normal.dot(&normal).as_rational()?.sign() == catallax_numbers::Sign::Positive {
            1
        } else {
            -1
        };
    Some((&circle.center, orientation))
}

fn torus_parallel(
    surface: &Surface,
    model: &Model,
    curve: catallax_shapes::CurveId,
    dec: &mut Decider,
) -> Result<Option<Parallel>, Refusal> {
    let Surface::Torus(torus) = surface else {
        return Ok(None);
    };
    let catallax_intersect::CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ..
    } = catallax_intersect::curve_form(model.geom(), curve, dec)?
    else {
        return Ok(None);
    };
    match dec.sign(&normal.cross(&torus.axis).norm_squared()) {
        Outcome::ProvenCoincident(_) => {}
        Outcome::Proven(_) => return Ok(None),
        Outcome::Refused(e) => return Err(e),
        Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
    }
    let axis = &torus.axis;
    let offset = center.minus(&torus.center);
    match dec.sign(&offset.cross(axis).norm_squared()) {
        Outcome::ProvenCoincident(_) => {}
        Outcome::Proven(_) => return Ok(None),
        Outcome::Refused(e) => return Err(e),
        Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
    }
    let minor = &torus.minor;
    let cosine = &(&radius_sq.sqrt() - &torus.major) / minor;
    let sine = offset.dot(axis) / (minor * &torus.axis.norm_squared().sqrt());
    let orientation = match dec.sign(&normal.dot(axis)) {
        Outcome::Proven(catallax_numbers::StrictSign::Positive) => 1,
        Outcome::Proven(_) | Outcome::ProvenCoincident(_) => -1,
        Outcome::Refused(e) => return Err(e),
        Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
    };
    Ok(Some(Parallel {
        cosine,
        sine,
        orientation,
    }))
}

fn cone_generator(
    surface: &Surface,
    model: &Model,
    curve: catallax_shapes::CurveId,
    dec: &mut Decider,
) -> Result<Option<[Rational; 5]>, Refusal> {
    let Surface::Cone(cone) = surface else {
        return Ok(None);
    };
    if !matches!(
        model.geom().curve(curve),
        catallax_shapes::Curve::Line(_) | catallax_shapes::Curve::Intersection(_)
    ) {
        return Ok(None);
    }
    let catallax_intersect::CurveForm::Line { origin, dir } =
        catallax_intersect::curve_form(model.geom(), curve, dec)?
    else {
        return Ok(None);
    };
    let dot = dir.dot(&cone.axis);
    let equation =
        dir.cross(&cone.axis).norm_squared() - (&cone.slope * &cone.slope) * (&dot * &dot);
    match dec.sign(&equation) {
        Outcome::ProvenCoincident(_) => {}
        Outcome::Proven(_) => return Ok(None),
        Outcome::Refused(e) => return Err(e),
        Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
    }
    let cross = origin.minus(&cone.apex).cross(&dir);
    match dec.sign(&cross.norm_squared()) {
        Outcome::ProvenCoincident(_) => {}
        Outcome::Proven(_) => return Ok(None),
        Outcome::Refused(e) => return Err(e),
        Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
    }
    // Along a generator q x dq = 0. Every cone boundary form has that
    // factor, so its integral is exactly zero, including an apex endpoint.
    // Evaluating the removable 1/z coordinate singularity loses this proof.
    Ok(Some(std::array::from_fn(|_| Rational::zero())))
}
