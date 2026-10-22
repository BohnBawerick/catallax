//! Deterministic per-property refinement and orientation handling.

use catallax_intersect::{CurveTrace, trace};
use catallax_numbers::{Answer, Answered, Ball, Decider, Outcome, Rational, Real, Refusal};
use catallax_shapes::{Point3, RPoint3, RVector3, Surface, Vector3};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};
use std::collections::{HashMap, hash_map::Entry};

use crate::integrate::{Stop, Work};
use crate::planar::{self, Planar};
use crate::surface::{self, Compiled, Property};
use crate::{Accuracy, Enclosure};

pub(crate) type TraceCache = HashMap<EdgeId, std::sync::Arc<CurveTrace>>;

struct Arc<'a> {
    trace: Option<std::sync::Arc<CurveTrace>>,
    reversed: bool,
    analytic: Option<crate::integrate::Program>,
    integral: Option<crate::integrate::Analytic>,
    exact: Option<[Rational; 5]>,
    meridian: Option<(&'a Point3, i64)>,
    parallel: Option<Parallel>,
}
struct Parallel {
    cosine: Real,
    sine: Real,
    orientation: i64,
}

struct Curved<'a> {
    surface: &'a Surface,
    sense: Sense,
    boundaries: Vec<Vec<Arc<'a>>>,
    compiled: Vec<Compiled>,
    torus: Option<crate::torus::Programs>,
}

enum Prepared<'a> {
    Curved(Curved<'a>),
    PlanarArea(Planar),
}

pub(crate) fn anchor(model: &Model, faces: &[FaceId]) -> RPoint3 {
    let point = match model.geom().surface(model.face(faces[0]).surface) {
        Surface::Plane(p) => &p.origin,
        Surface::Cylinder(c) => &c.origin,
        Surface::Cone(c) => &c.apex,
        Surface::Sphere(s) => &s.center,
        Surface::Torus(t) => &t.center,
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
        Property::Centroid => vec![1, 2, 3, 4],
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
    let mut prepared = Vec::new();
    for (id, planar) in pending {
        if let Some(planar) = planar {
            prepared.push(Prepared::PlanarArea(planar));
            continue;
        }
        let face = model.face(id);
        let surface = model.geom().surface(face.surface);
        let closed = face.loops.is_empty();
        let torus = match surface {
            Surface::Torus(t) if !closed => Some(crate::torus::compile(t, anchor, property)),
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
        let mut boundaries = Vec::new();
        for id in &face.loops {
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
                        let analytic = circle.and_then(|circle| {
                            crate::ring::compile(surface, circle, anchor, &compiled.indices)
                        });
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
                        let integral = if torus.is_none() && exact.is_none() && analytic.is_none() {
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
                        let trace = if exact.is_none()
                            && analytic.is_none()
                            && meridian.is_none()
                            && parallel.is_none()
                            && (integral.is_none() || matches!(surface, Surface::Sphere(_)))
                        {
                            let traced = match traces.entry(c.edge) {
                                Entry::Occupied(entry) => entry.into_mut(),
                                Entry::Vacant(entry) => entry.insert(std::sync::Arc::new(
                                    trace(model.geom(), edge.curve, from, to, dec)?.value,
                                )),
                            };
                            Some(traced.clone())
                        } else {
                            None
                        };
                        arcs.push(Arc {
                            trace,
                            reversed: (edge.sense == Sense::Reversed)
                                != (c.sense == Sense::Reversed),
                            analytic,
                            integral,
                            exact,
                            meridian,
                            parallel,
                        });
                    }
                }
            }
            if !arcs.is_empty() {
                boundaries.push(arcs);
            }
        }
        prepared.push(Prepared::Curved(Curved {
            surface,
            sense: face.sense,
            boundaries,
            compiled: vec![compiled],
            torus,
        }));
    }
    let mut best = None;
    let mut remaining = policy.max_work;
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
            Err(Stop::Split) => {
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
        if property == Property::Centroid && volume_enclosure_reached {
            Refusal::new(
                "measure.centroid-not-separated",
                "volume could not be separated from zero",
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
    if property != Property::Centroid {
        return Ok(Answer::exact(result));
    }
    let volume = &result[0];
    if volume.lo <= Rational::zero() && volume.hi >= Rational::zero() {
        return Err(Refusal::new(
            "measure.centroid-not-separated",
            "volume enclosure includes zero",
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
    work: &mut Work,
) -> Result<Vec<Enclosure>, Stop> {
    let mut sums: Vec<Ball> = indices
        .iter()
        .map(|i| work.rational(&exact[*i]))
        .collect::<Result<_, _>>()?;
    for face in prepared {
        let values = match face {
            Prepared::PlanarArea(p) => {
                let numerator = work.rational(&p.area_numerator)?;
                let square = work.rational(&p.normal_squared)?;
                let norm = work.sqrt(&square)?;
                vec![work.div(&numerator, &norm)?]
            }
            Prepared::Curved(face) => {
                let mut answer = None;
                let candidates =
                    if matches!(face.surface, Surface::Sphere(_)) && !face.boundaries.is_empty() {
                        6
                    } else {
                        1
                    };
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
        return integrate_torus(face, torus, indices, depth, rectangle, work);
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
    programs: &crate::torus::Programs,
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
    let sense = work.rational(&planar::sign(face.sense))?;
    for sum in &mut sums {
        *sum = work.mul(sum, &sense)?;
    }
    let pi = work.pi()?;
    let full = programs
        .full
        .point(&[zero.clone(), zero.clone(), zero], &[pi], work)?;
    correct_whole_surface(&mut sums, &full, 0, work)?;
    indices
        .iter()
        .map(|i| {
            let value = &sums[programs.indices.iter().position(|j| i == j).unwrap()];
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
