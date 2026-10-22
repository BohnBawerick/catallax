//! Closed-form winding for polygon loops and complete circular boundaries.
//! Other loops use certified traces, with the same orientation convention.

use crate::geometry::{product, refusal, sign};
use crate::region::Location;
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Refusal, Sign};
use catallax_shapes::{PointId, Surface};
use catallax_topology::{FaceId, Loop, LoopId, Model, Sense};

fn ring(model: &Model, l: LoopId) -> Option<(catallax_topology::EdgeId, Sense)> {
    let Loop::Coedges(cs) = model.loop_(l) else {
        return None;
    };
    if cs.len() != 1 {
        return None;
    }
    let c = model.coedge(cs[0]);
    let e = model.edge(c.edge);
    if e.start.is_some() || e.end.is_some() {
        return None;
    }
    Some((c.edge, product(e.sense, c.sense)))
}

fn ring_side(
    model: &Model,
    l: LoopId,
    p: PointId,
    dec: &mut Decider,
) -> Result<Option<(Sign, Sense)>, Refusal> {
    let Some((edge, travel)) = ring(model, l) else {
        return Ok(None);
    };
    let CurveForm::Circle { center, normal, .. } =
        curve_form(model.geom(), model.edge(edge).curve, dec)?
    else {
        return Ok(None);
    };
    let value = model
        .geom()
        .position(p)
        .minus(&center)
        .dot(&normal.to_real());
    Ok(Some((
        sign(&value, dec, "boolean.classify-ambiguous")?,
        travel,
    )))
}

fn plane_winding(
    model: &Model,
    surface: &catallax_shapes::Plane,
    l: LoopId,
    p: PointId,
    dec: &mut Decider,
) -> Result<Option<i64>, Refusal> {
    let query = model.geom().position(p);
    if let Some((edge, travel)) = ring(model, l) {
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = curve_form(model.geom(), model.edge(edge).curve, dec)?
        else {
            return Ok(None);
        };
        let side = sign(
            &(query.minus(&center).norm_squared() - radius_sq),
            dec,
            "boolean.classify-ambiguous",
        )?;
        if side == Sign::Zero {
            return Ok(None);
        }
        if side == Sign::Positive {
            return Ok(Some(0));
        }
        let agree = sign(
            &surface.normal.to_real().dot(&normal.to_real()),
            dec,
            "boolean.classify-ambiguous",
        )? == Sign::Positive;
        return Ok(Some(if agree == (travel == Sense::Same) {
            1
        } else {
            -1
        }));
    }
    let Loop::Coedges(cs) = model.loop_(l) else {
        return Ok(None);
    };
    let u = surface.ref_dir.to_real();
    let v = surface.normal.cross(&surface.ref_dir).to_real();
    let origin = surface.origin.to_real();
    let d = query.minus(&origin);
    let qx = d.dot(&u);
    let qy = d.dot(&v);
    let project = |vertex| {
        let d = model
            .geom()
            .position(model.vertex(vertex).point)
            .minus(&origin);
        (d.dot(&u), d.dot(&v))
    };
    let mut previous = None;
    let mut result = 0;
    for &coedge in cs {
        let coedge = model.coedge(coedge);
        let edge = model.edge(coedge.edge);
        if !matches!(
            curve_form(model.geom(), edge.curve, dec)?,
            CurveForm::Line { .. }
        ) {
            return Ok(None);
        }
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Ok(None);
        };
        let (start, end) = if coedge.sense == Sense::Same {
            (start, end)
        } else {
            (end, start)
        };
        let (ax, ay) = match previous.take() {
            Some((vertex, coordinates)) if vertex == start => coordinates,
            _ => project(start),
        };
        let (bx, by) = project(end);
        let dx = &qx - &ax;
        let dy = &qy - &ay;
        let cross = &(&bx - &ax) * &dy - &(&by - &ay) * &dx;
        let cross = sign(&cross, dec, "boolean.classify-ambiguous")?;
        let ex = &qx - &bx;
        if cross == Sign::Zero {
            let between = &dx * &ex + &dy * &(&qy - &by);
            if sign(&between, dec, "boolean.classify-ambiguous")? != Sign::Positive {
                return Ok(None);
            }
        }
        let left = sign(&dx, dec, "boolean.classify-ambiguous")?;
        let right = sign(&ex, dec, "boolean.classify-ambiguous")?;
        if left != Sign::Negative && right == Sign::Negative && cross == Sign::Negative {
            result -= 1;
        }
        if right != Sign::Negative && left == Sign::Negative && cross == Sign::Positive {
            result += 1;
        }
        previous = Some((end, (bx, by)));
    }
    Ok(Some(result))
}

pub(crate) fn signature(
    model: &Model,
    surface: &Surface,
    loops: &[LoopId],
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<i64>, Refusal> {
    match surface {
        Surface::Plane(plane) => {
            let mut value = 0;
            for &l in loops {
                let Some(w) = plane_winding(model, plane, l, point, dec)? else {
                    return Ok(None);
                };
                value += w;
            }
            Ok(Some(value))
        }
        Surface::Sphere(_) | Surface::Cylinder(_) | Surface::Cone(_) => {
            if loops.len() != 1 {
                return Ok(None);
            }
            if matches!(model.loop_(loops[0]), Loop::Vertex(_)) {
                return Ok(Some(0));
            }
            let Some((side, _)) = ring_side(model, loops[0], point, dec)? else {
                return Ok(None);
            };
            match side {
                Sign::Positive => Ok(Some(1)),
                Sign::Negative => Ok(Some(0)),
                Sign::Zero => Err(refusal(
                    "boolean.classify-ambiguous",
                    "cell probe is on a circular boundary",
                )),
            }
        }
        _ => Ok(None),
    }
}

pub(crate) fn locate(
    model: &Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let f = model.face(face);
    let surface = model.geom().surface(f.surface);
    match surface {
        Surface::Plane(plane) => {
            let mut value = 0;
            for &l in &f.loops {
                let Some(w) = plane_winding(model, plane, l, point, dec)? else {
                    return Ok(None);
                };
                value += w;
            }
            if f.sense == Sense::Reversed {
                value = -value;
            }
            match value {
                0 => Ok(Some(Location::Outside)),
                1 => Ok(Some(Location::Inside)),
                _ => Err(refusal(
                    "boolean.classify-ambiguous",
                    "polygon winding is not a unit face region",
                )),
            }
        }
        Surface::Sphere(_) | Surface::Cylinder(_) | Surface::Cone(_) => {
            if let Surface::Cone(cone) = surface
                && sign(
                    &model
                        .geom()
                        .position(point)
                        .minus(&cone.apex.to_real())
                        .dot(&cone.axis.to_real()),
                    dec,
                    "boolean.classify-ambiguous",
                )? == Sign::Negative
            {
                return Ok(Some(Location::Outside));
            }
            let mut boundary = false;
            for &l in &f.loops {
                if let Loop::Vertex(vertex) = model.loop_(l) {
                    if crate::geometry::coincident(
                        model
                            .geom()
                            .same_point(model.vertex(*vertex).point, point, dec),
                        "boolean.classify-ambiguous",
                    )? {
                        return Ok(Some(Location::Boundary));
                    }
                    continue;
                }
                let Some((side, travel)) = ring_side(model, l, point, dec)? else {
                    return Ok(None);
                };
                let front = travel == f.sense;
                match side {
                    Sign::Zero => boundary = true,
                    Sign::Positive if !front => return Ok(Some(Location::Outside)),
                    Sign::Negative if front => return Ok(Some(Location::Outside)),
                    _ => {}
                }
            }
            Ok(Some(if boundary {
                Location::Boundary
            } else {
                Location::Inside
            }))
        }
        _ => Ok(None),
    }
}
