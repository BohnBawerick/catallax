//! Product trims bounded by two circular arcs and two axial generators.

use super::{Location, product, sign};
use crate::geometry::natural_ends;
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, FactKind, Refusal, Sign};
use catallax_shapes::{Circle, Curve, Cylinder, GeomKey, Point3, PointId, Surface, Vector3};
use catallax_topology::{Coedge, Edge, FaceId, Loop, Model, Sense};

const CODE: &str = "boolean.classify-ambiguous";

struct Cap<'a> {
    center: &'a Point3,
    normal: &'a Vector3,
    from: PointId,
    to: PointId,
    above: bool,
}

fn cap<'a>(
    model: &'a Model,
    coedge: &Coedge,
    edge: &Edge,
    cylinder: &Cylinder,
    sense: Sense,
    dec: &mut Decider,
) -> Result<Option<Cap<'a>>, Refusal> {
    let Curve::Circle(Circle {
        center,
        normal,
        radius_sq,
        ..
    }) = model.geom().curve(edge.curve)
    else {
        return Ok(None);
    };
    if sign(&normal.cross(&cylinder.axis).norm_squared(), dec, CODE)? != Sign::Zero
        || sign(
            &center
                .minus(&cylinder.origin)
                .cross(&cylinder.axis)
                .norm_squared(),
            dec,
            CODE,
        )? != Sign::Zero
        || sign(&(radius_sq - &cylinder.radius_sq), dec, CODE)? != Sign::Zero
    {
        return Ok(None);
    }
    let align = sign(&normal.dot(&cylinder.axis), dec, CODE)?;
    if align == Sign::Zero {
        return Ok(None);
    }
    let (Some(from), Some(to)) = natural_ends(model, coedge.edge) else {
        return Ok(None);
    };
    let front = product(edge.sense, coedge.sense) == sense;
    Ok(Some(Cap {
        center,
        normal,
        from,
        to,
        above: (align == Sign::Positive) == front,
    }))
}

/// Like the other native face queries, the query is on the face's support.
/// Incidence only verifies the carriers here. Closure, opposite cap travel,
/// ordered heights and distinct oriented generators establish the finite trim.
pub(super) fn locate(
    model: &Model,
    face: FaceId,
    cylinder: &Cylinder,
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let face = model.face(face);
    let [loop_id] = face.loops.as_slice() else {
        return Ok(None);
    };
    let Loop::Coedges(coedges) = model.loop_(*loop_id) else {
        return Ok(None);
    };
    if coedges.len() != 4 {
        return Ok(None);
    }
    let native_circle = |index: usize| {
        matches!(
            model
                .geom()
                .curve(model.edge(model.coedge(coedges[index]).edge).curve),
            Curve::Circle(_)
        )
    };
    let offset = usize::from(!native_circle(0));
    let edges: [(&Coedge, &Edge); 4] = std::array::from_fn(|index| {
        let coedge = model.coedge(coedges[(offset + index) % 4]);
        (coedge, model.edge(coedge.edge))
    });
    for (index, (coedge, edge)) in edges.iter().enumerate() {
        let supported = matches!(
            (index % 2, model.geom().curve(edge.curve)),
            (0, Curve::Circle(_)) | (1, Curve::Line(_) | Curve::Intersection(_))
        );
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Ok(None);
        };
        if !supported {
            return Ok(None);
        }
        if let Curve::Intersection(branch) = model.geom().curve(edge.curve) {
            let (Surface::Cylinder(a), Surface::Cylinder(b)) = (
                model.geom().surface(branch.surfaces[0]),
                model.geom().surface(branch.surfaces[1]),
            ) else {
                return Ok(None);
            };
            let parallel = a.axis.cross(&b.axis).norm_squared();
            if parallel.rational_ref().is_none() || sign(&parallel, dec, CODE)? != Sign::Zero {
                return Ok(None);
            }
        }
        let end = if coedge.sense == Sense::Same {
            end
        } else {
            start
        };
        let (next_coedge, next_edge) = edges[(index + 1) % 4];
        let next_start = if next_coedge.sense == Sense::Same {
            next_edge.start
        } else {
            next_edge.end
        };
        if Some(end) != next_start {
            return Ok(None);
        }
        // Missing retained incidence is inconclusive, not permission to infer
        // that topological connectivity establishes geometric connectivity.
        for vertex in [edge.start.unwrap(), edge.end.unwrap()] {
            if !model.geom().facts().holds(
                FactKind::Incident,
                GeomKey::Point(model.vertex(vertex).point),
                GeomKey::Curve(edge.curve),
            ) {
                return Ok(None);
            }
            dec.note_fact_decision();
        }
    }
    let Some(first) = cap(model, edges[0].0, edges[0].1, cylinder, face.sense, dec)? else {
        return Ok(None);
    };
    let Some(second) = cap(model, edges[2].0, edges[2].1, cylinder, face.sense, dec)? else {
        return Ok(None);
    };
    if first.above == second.above {
        return Ok(None);
    }
    let (lower, upper) = if first.above {
        (&first, &second)
    } else {
        (&second, &first)
    };
    if sign(
        &upper.center.minus(lower.center).dot(&cylinder.axis),
        dec,
        CODE,
    )? != Sign::Positive
    {
        return Ok(None);
    }
    for (index, expected_up) in [(1, first.above), (3, !first.above)] {
        let (coedge, edge) = edges[index];
        let form;
        let dir = match model.geom().curve(edge.curve) {
            Curve::Line(line) => &line.dir,
            _ => {
                form = curve_form(model.geom(), edge.curve, dec)?;
                let CurveForm::Line { dir, .. } = &form else {
                    return Ok(None);
                };
                dir
            }
        };
        if sign(&dir.cross(&cylinder.axis).norm_squared(), dec, CODE)? != Sign::Zero
            || sign(&dir.norm_squared(), dec, CODE)? != Sign::Positive
        {
            return Ok(None);
        }
        // An intersection form's direction need not be its natural direction.
        // Finite line travel comes from the coedge endpoints, not that direction.
        let start = model
            .geom()
            .position(model.vertex(edge.start.unwrap()).point);
        let end = model.geom().position(model.vertex(edge.end.unwrap()).point);
        let direction = sign(&end.minus(&start).dot(&cylinder.axis), dec, CODE)?;
        if direction == Sign::Zero
            || ((direction == Sign::Positive) == (coedge.sense == Sense::Same)) != expected_up
        {
            return Ok(None);
        }
    }
    let from = model.geom().position(first.from);
    let to = model.geom().position(first.to);
    // Distinct points on one perpendicular cap cannot share an axial line.
    if sign(&to.minus(&from).norm_squared(), dec, CODE)? != Sign::Positive {
        return Ok(None);
    }
    // The generator lines connect corresponding cap endpoints. Opposite cap
    // travel then selects the same angular arc at both heights, even for a
    // major arc or a reversed normal. The region is their Cartesian product.
    let query = model.geom().position(point);
    let bottom = sign(&query.minus(lower.center).dot(&cylinder.axis), dec, CODE)?;
    let top = sign(&query.minus(upper.center).dot(&cylinder.axis), dec, CODE)?;
    if bottom == Sign::Negative || top == Sign::Positive {
        return Ok(Some(Location::Outside));
    }
    // On a circle, its forward arc lies to the right of the oriented chord.
    // The axial component of query-from vanishes in this scalar triple product.
    let angular = sign(
        &to.minus(&from).cross(&query.minus(&from)).dot(first.normal),
        dec,
        CODE,
    )?;
    Ok(Some(match angular {
        Sign::Positive => Location::Outside,
        Sign::Zero => Location::Boundary,
        Sign::Negative if bottom == Sign::Zero || top == Sign::Zero => Location::Boundary,
        Sign::Negative => Location::Inside,
    }))
}

#[cfg(test)]
mod tests;
