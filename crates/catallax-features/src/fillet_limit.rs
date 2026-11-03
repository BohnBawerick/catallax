//! Proofs that a fillet radius is too large for its edge.
//!
//! A proof shows that a rolling ball of the requested radius cannot rest on
//! both faces of the edge: either no ball centre touches both supports on
//! their material sides, or the ball's contact with one face lies wholly
//! outside that face. Each proof also gives the radius past which it holds.
//! Proofs only relabel a refusal the construction already returned, so they
//! never change what builds.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Curve, Cylinder, Line, Plane, Point3, Surface, Vector3};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};

use crate::{
    fillet::{FilletDefinition, RADIUS_TOO_LARGE, too_large},
    profile::sign,
};

type Uses = BTreeMap<EdgeId, Vec<(FaceId, Sense)>>;

const NO_CENTRE: &str = "no ball of this radius rolls along the edge against both faces";
const NO_CONTACT: &str = "the ball of this radius cannot touch one face of the edge";

fn root(value: &Real) -> Real {
    value
        .as_rational()
        .and_then(|r| crate::extrude::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| value.sqrt())
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &root(&v.norm_squared())))
}
fn toward(same: bool) -> Real {
    Real::from_i64(if same { 1 } else { -1 })
}
fn larger(a: Real, b: Real, dec: &mut Decider) -> Result<Real, Refusal> {
    Ok(if sign(&(&b - &a), dec)? == Sign::Positive {
        b
    } else {
        a
    })
}
fn smaller(a: Real, b: Real, dec: &mut Decider) -> Result<Real, Refusal> {
    Ok(if sign(&(&b - &a), dec)? == Sign::Negative {
        b
    } else {
        a
    })
}

/// A condition `start + slope r > 0` that every ball centre of radius `r` needs.
struct Affine {
    start: Real,
    slope: Real,
}
impl Affine {
    fn new(start: Real, slope: Real) -> Self {
        Self { start, slope }
    }
    /// The offset `radius - s r` of a cylinder, where `s` is 1 when the ball
    /// sits inside it.
    fn offset(cylinder: &Cylinder, inside: &Real) -> Self {
        Self::new(cylinder.radius(), -inside)
    }
}

/// The radius at which the first condition fails, when every condition holds
/// strictly at radius zero, so on the whole interval below it.
fn first_failure(conditions: &[Affine], dec: &mut Decider) -> Result<Option<Real>, Refusal> {
    let mut limit: Option<Real> = None;
    for condition in conditions {
        if sign(&condition.start, dec)? != Sign::Positive {
            return Ok(None);
        }
        if sign(&condition.slope, dec)? != Sign::Negative {
            continue;
        }
        let at = &condition.start / &(-&condition.slope);
        limit = Some(match limit {
            Some(limit) => smaller(limit, at, dec)?,
            None => at,
        });
    }
    Ok(limit)
}

/// The foot on a cylinder's axis of the normal section through a point.
fn axis_foot(cylinder: &Cylinder, at: &Point3) -> Point3 {
    cylinder.origin.plus(
        &cylinder.axis.scaled(
            &(at.minus(&cylinder.origin).dot(&cylinder.axis) / cylinder.axis.norm_squared()),
        ),
    )
}

/// The conditions on a ball centre for a straight generator between a plane
/// and a parallel cylinder: the centre lies on the plane's offset line and on
/// the cylinder's offset circle, which must cross.
fn plane_cylinder(
    plane: &Plane,
    cylinder: &Cylinder,
    edge: &Line,
    same: [bool; 2],
    k: &Real,
) -> Vec<Affine> {
    let inside = k * &toward(same[1]);
    let outward = unit(&plane.normal).scaled(&toward(same[0]));
    let height = edge
        .origin
        .minus(&axis_foot(cylinder, &edge.origin))
        .dot(&outward);
    let radius = cylinder.radius();
    vec![
        Affine::offset(cylinder, &inside),
        // offset - across, with across = height - k r.
        Affine::new(&radius - &height, k - &inside),
        // offset + across.
        Affine::new(&radius + &height, -&(k + &inside)),
    ]
}

/// The conditions on a ball centre for a straight generator between two
/// parallel cylinders: their offset circles must cross.
fn cylinders(a: &Cylinder, b: &Cylinder, edge: &Line, same: [bool; 2], k: &Real) -> Vec<Affine> {
    let inside = same.map(|same| k * &toward(same));
    let apart = root(
        &axis_foot(b, &edge.origin)
            .minus(&axis_foot(a, &edge.origin))
            .norm_squared(),
    );
    let [ra, rb] = [a.radius(), b.radius()];
    vec![
        Affine::offset(a, &inside[0]),
        Affine::offset(b, &inside[1]),
        Affine::new(&(&ra + &rb) - &apart, -&(&inside[0] + &inside[1])),
        Affine::new(&(&apart - &ra) + &rb, &inside[0] - &inside[1]),
        Affine::new(&(&apart + &ra) - &rb, &inside[1] - &inside[0]),
    ]
}

/// The conditions on a ball centre for the ring where two nonparallel
/// cylinders cross: the offsets must cross the way the cylinders themselves
/// do, as a ring through the wide one or as one loop.
fn saddle(
    first: &Cylinder,
    second: &Cylinder,
    same: [bool; 2],
    k: &Real,
    dec: &mut Decider,
) -> Result<Vec<Affine>, Refusal> {
    let first_wide = sign(&(&first.radius_sq - &second.radius_sq), dec)? != Sign::Negative;
    let ([wide, narrow], sides) = if first_wide {
        ([first, second], same)
    } else {
        ([second, first], [same[1], same[0]])
    };
    let [sw, sn] = sides.map(|same| k * &toward(same));
    let normal = wide.axis.cross(&narrow.axis);
    let squared = normal.norm_squared();
    if sign(&squared, dec)? == Sign::Zero {
        return Ok(Vec::new());
    }
    let reach = root(&(narrow.origin.minus(&wide.origin).dot(&normal).square() / squared));
    let [rw, rn] = [wide.radius(), narrow.radius()];
    let mut conditions = vec![Affine::offset(wide, &sw), Affine::offset(narrow, &sn)];
    let clear = &(&rw - &rn) - &reach;
    if sign(&clear, dec)? == Sign::Positive {
        conditions.push(Affine::new(clear, &sn - &sw));
    } else {
        conditions.extend([
            Affine::new(&(&reach - &rw) + &rn, &sw - &sn),
            Affine::new(&(&reach + &rw) - &rn, &sn - &sw),
            Affine::new(&(&rw + &rn) - &reach, -&(&sw + &sn)),
        ]);
    }
    Ok(conditions)
}

/// The radius past which no ball centre touches both supports on their
/// material sides, when one exists. A cylinder the ball sits inside caps the
/// radius at its own; the other conditions depend on the family.
fn contact_limit(
    supports: &[Surface; 2],
    same: [bool; 2],
    edge: &Curve,
    convex: bool,
    dec: &mut Decider,
) -> Result<Option<Real>, Refusal> {
    let k = Real::from_i64(if convex { 1 } else { -1 });
    let conditions = match (&supports[0], &supports[1], edge) {
        (Surface::Plane(p), Surface::Cylinder(c), Curve::Line(e)) => {
            plane_cylinder(p, c, e, same, &k)
        }
        (Surface::Cylinder(c), Surface::Plane(p), Curve::Line(e)) => {
            plane_cylinder(p, c, e, [same[1], same[0]], &k)
        }
        (Surface::Cylinder(a), Surface::Cylinder(b), Curve::Line(e)) => {
            cylinders(a, b, e, same, &k)
        }
        (Surface::Cylinder(a), Surface::Cylinder(b), Curve::Intersection(_)) => {
            saddle(a, b, same, &k, dec)?
        }
        _ => supports
            .iter()
            .zip(same)
            .filter_map(|(support, same)| match support {
                Surface::Cylinder(c) => Some(Affine::offset(c, &(&k * &toward(same)))),
                _ => None,
            })
            .collect(),
    };
    first_failure(&conditions, dec)
}

/// Whether a circular edge's convexity, found from both faces' normals at one
/// point, agrees with `convex`. The construction reads a circle's convexity
/// from the plane face alone, which a hole's mouth or floor can contradict; no
/// proof is made then. Other edges already use both normals.
pub(crate) fn convexity_agrees(
    model: &Model,
    edge: EdgeId,
    incident: &[(FaceId, Sense)],
    convex: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let e = model.edge(edge);
    let Curve::Circle(c) = model.geom().curve(e.curve) else {
        return Ok(true);
    };
    let [(a, sense), (b, _)] = incident else {
        return Ok(false);
    };
    let outward = |face: FaceId, at: &Point3, dec: &mut Decider| -> Result<Vector3, Refusal> {
        let f = model.face(face);
        let n = model.geom().surface(f.surface).normal_at(at, dec)?;
        Ok(n.scaled(&toward(f.sense == Sense::Same)))
    };
    let at = c.center.plus(&unit(&c.ref_dir).scaled(&root(&c.radius_sq)));
    let along = c.normal.cross(&c.ref_dir);
    let orientation = toward(*sense == e.sense);
    let value = outward(*b, &at, dec)?.dot(&outward(*a, &at, dec)?.cross(&along)) * orientation;
    Ok(match sign(&value, dec)? {
        Sign::Negative => convex,
        Sign::Positive => !convex,
        Sign::Zero => false,
    })
}

/// Relabels a refusal from building a strip's definition when no ball centre
/// of this radius exists. Edge and radius validation keep their own codes.
pub(crate) fn contact(
    refusal: Refusal,
    supports: &[Surface; 2],
    same: [bool; 2],
    edge: &Curve,
    radius: &Real,
    convex: bool,
    dec: &mut Decider,
) -> Refusal {
    if matches!(
        refusal.code,
        "fillet.invalid-edge" | "fillet.non-positive-radius" | RADIUS_TOO_LARGE
    ) {
        return refusal;
    }
    let proof = (|| -> Result<Option<Real>, Refusal> {
        let Some(limit) = contact_limit(supports, same, edge, convex, dec)? else {
            return Ok(None);
        };
        Ok((sign(&(radius - &limit), dec)? == Sign::Positive).then_some(limit))
    })();
    match proof {
        Ok(Some(limit)) => too_large(NO_CENTRE, Some(&limit), dec),
        _ => refusal,
    }
}

/// A point of an edge: its start, or a point of a closed conic.
fn point_on(model: &Model, edge: EdgeId) -> Option<Point3> {
    let e = model.edge(edge);
    if let Some(start) = e.start {
        return Some(model.geom().position(model.vertex(start).point));
    }
    match model.geom().curve(e.curve) {
        Curve::Circle(c) => Some(c.center.plus(&unit(&c.ref_dir).scaled(&root(&c.radius_sq)))),
        Curve::Ellipse(c) => Some(
            c.center
                .plus(&unit(&c.major_dir).scaled(&root(&c.major_sq))),
        ),
        _ => None,
    }
}

/// Whether a face meets its neighbour across an edge at a fixed nonzero
/// angle, so a contact leaving the face cannot continue smoothly onto the
/// neighbour. Only pairs whose angle cannot reach zero part way along the
/// edge are decided from one point; others prove nothing.
fn sharp(
    model: &Model,
    face: FaceId,
    other: FaceId,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [a, b] = [face, other].map(|f| model.geom().surface(model.face(f).surface));
    let curve = model.geom().curve(model.edge(edge).curve);
    let fixed = match (a, b) {
        (Surface::Plane(_), Surface::Plane(_)) => true,
        (Surface::Plane(_), Surface::Cylinder(_)) | (Surface::Cylinder(_), Surface::Plane(_)) => {
            matches!(curve, Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_))
        }
        (Surface::Cylinder(_), Surface::Cylinder(_)) => matches!(curve, Curve::Line(_)),
        _ => false,
    };
    let Some(at) = point_on(model, edge).filter(|_| fixed) else {
        return Ok(false);
    };
    let [na, nb] = [a.normal_at(&at, dec)?, b.normal_at(&at, dec)?];
    Ok(sign(&na.cross(&nb).norm_squared(), dec)? != Sign::Zero)
}

/// The edges bounding a face, or `None` when the face has no edge loop or a
/// neighbour may continue it smoothly.
fn boundary(
    model: &Model,
    uses: &Uses,
    face: FaceId,
    selected: EdgeId,
    dec: &mut Decider,
) -> Result<Option<Vec<EdgeId>>, Refusal> {
    let loops = &model.face(face).loops;
    if loops.is_empty() {
        return Ok(None);
    }
    let mut edges = Vec::new();
    for lp in loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Ok(None);
        };
        for coedge in coedges {
            let edge = model.coedge(*coedge).edge;
            edges.push(edge);
            if edge == selected {
                continue;
            }
            let Some(incident) = uses.get(&edge) else {
                return Ok(None);
            };
            for (other, _) in incident {
                if *other != face && !sharp(model, face, *other, edge, dec)? {
                    return Ok(None);
                }
            }
        }
    }
    Ok(Some(edges))
}

/// An upper bound of `(p - origin) . w` over a face's boundary edges, which
/// bounds it over the face: a linear function has no maximum inside a plane
/// region, nor inside a cylinder region when `w` runs along the axis.
fn highest(
    model: &Model,
    edges: &[EdgeId],
    origin: &Point3,
    w: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Real>, Refusal> {
    let f = |p: &Point3| p.minus(origin).dot(w);
    let mut best: Option<Real> = None;
    for edge in edges {
        let e = model.edge(*edge);
        let value = match model.geom().curve(e.curve) {
            Curve::Line(_) => {
                let (Some(start), Some(end)) = (e.start, e.end) else {
                    return Ok(None);
                };
                let [a, b] = [start, end].map(|v| f(&model.geom().position(model.vertex(v).point)));
                larger(a, b, dec)?
            }
            Curve::Circle(c) => {
                let n = c.normal.dot(w);
                let across = &w.norm_squared() - &(n.square() / c.normal.norm_squared());
                let centre = f(&c.center);
                if sign(&across, dec)? == Sign::Zero {
                    centre
                } else {
                    &centre + &(&root(&c.radius_sq) * &root(&across))
                }
            }
            Curve::Ellipse(c) => &f(&c.center) + &(&root(&c.major_sq) * &root(&w.norm_squared())),
            _ => return Ok(None),
        };
        best = Some(match best {
            Some(best) => larger(best, value, dec)?,
            None => value,
        });
    }
    Ok(best)
}

/// Lower and upper bounds of the distance from `centre` to a face's boundary
/// edges.
fn distances(
    model: &Model,
    edges: &[EdgeId],
    centre: &Point3,
    dec: &mut Decider,
) -> Result<Option<(Real, Real)>, Refusal> {
    let mut bounds: Option<(Real, Real)> = None;
    for edge in edges {
        let e = model.edge(*edge);
        let (low, high) = match model.geom().curve(e.curve) {
            Curve::Line(_) => {
                let (Some(start), Some(end)) = (e.start, e.end) else {
                    return Ok(None);
                };
                let [a, b] = [start, end].map(|v| model.geom().position(model.vertex(v).point));
                let along = b.minus(&a);
                let offset = centre.minus(&a);
                let foot = offset.dot(&along);
                let length_sq = along.norm_squared();
                let nearest_sq = if sign(&foot, dec)? != Sign::Positive {
                    offset.norm_squared()
                } else if sign(&(&foot - &length_sq), dec)? != Sign::Negative {
                    centre.minus(&b).norm_squared()
                } else {
                    offset.norm_squared() - foot.square() / length_sq
                };
                let far_sq = larger(offset.norm_squared(), centre.minus(&b).norm_squared(), dec)?;
                (root(&nearest_sq), root(&far_sq))
            }
            Curve::Circle(c) => {
                let d = root(&c.center.minus(centre).norm_squared());
                let r = root(&c.radius_sq);
                let gap = &d - &r;
                let low = if sign(&gap, dec)? == Sign::Negative {
                    -&gap
                } else {
                    gap
                };
                (low, &d + &r)
            }
            Curve::Ellipse(c) => {
                let d = root(&c.center.minus(centre).norm_squared());
                let a = root(&c.major_sq);
                (&d - &a, &d + &a)
            }
            _ => return Ok(None),
        };
        bounds = Some(match bounds {
            Some((lo, hi)) => (smaller(lo, low, dec)?, larger(hi, high, dec)?),
            None => (low, high),
        });
    }
    Ok(bounds)
}

/// A proof that one contact misses its face, with the radius past which it
/// does when the contact moves linearly with the radius.
fn misses(
    model: &mut Model,
    uses: &Uses,
    selected: EdgeId,
    face: FaceId,
    index: usize,
    definition: &FilletDefinition,
    dec: &mut Decider,
) -> Result<Option<Option<Real>>, Refusal> {
    let support = &definition.supports[index];
    let contact = &definition.contacts[index];
    let radius = &definition.radius;
    let Some(edges) = boundary(model, uses, face, selected, dec)? else {
        return Ok(None);
    };
    match (support, contact, &definition.edge) {
        // A straight contact on a plane, at a fixed distance from the edge.
        (Surface::Plane(_), Curve::Line(rail), Curve::Line(edge)) => {
            let d = rail.origin.minus(&edge.origin);
            let w = d.minus(
                &edge
                    .dir
                    .scaled(&(d.dot(&edge.dir) / edge.dir.norm_squared())),
            );
            let reach = w.norm_squared();
            if sign(&reach, dec)? == Sign::Zero {
                return Ok(None);
            }
            let Some(high) = highest(model, &edges, &edge.origin, &w, dec)? else {
                return Ok(None);
            };
            if sign(&(&reach - &high), dec)? != Sign::Positive {
                return Ok(None);
            }
            // Between two planes every contact scales with the radius.
            let planar = matches!(definition.supports, [Surface::Plane(_), Surface::Plane(_)]);
            Ok(Some(planar.then(|| &(radius * &high) / &reach)))
        }
        // The wall contact of a torus strip, a radius along the axis from the edge.
        (Surface::Cylinder(_), Curve::Circle(rail), Curve::Circle(edge)) => {
            let w = rail.center.minus(&edge.center);
            let reach = w.norm_squared();
            if sign(&reach, dec)? == Sign::Zero {
                return Ok(None);
            }
            let Some(high) = highest(model, &edges, &edge.center, &w, dec)? else {
                return Ok(None);
            };
            if sign(&(&reach - &high), dec)? != Sign::Positive {
                return Ok(None);
            }
            Ok(Some(Some(&(radius * &high) / &reach)))
        }
        // The plane contact of a torus strip, a circle about the edge's centre.
        (Surface::Plane(_), Curve::Circle(rail), Curve::Circle(edge)) => {
            let Some((low, high)) = distances(model, &edges, &edge.center, dec)? else {
                return Ok(None);
            };
            let [contact, own] = [&rail.radius_sq, &edge.radius_sq].map(root);
            let moved = &contact - &own;
            let bound = if sign(&(&contact - &high), dec)? == Sign::Positive {
                high
            } else if sign(&(&low - &contact), dec)? == Sign::Positive
                && !catallax_boolean::face_contains(model, face, &edge.center, dec)?
            {
                low
            } else {
                return Ok(None);
            };
            // The contact's radius is affine in the ball's.
            let limit = &(radius * &(&bound - &own)) / &moved;
            Ok(Some(
                (sign(&limit, dec)? == Sign::Positive).then_some(limit),
            ))
        }
        _ => Ok(None),
    }
}

/// Relabels a refusal from building strips when a contact of one of them
/// lies wholly outside its face, so the ball cannot touch that face at all.
pub(crate) fn width(
    model: &mut Model,
    uses: &Uses,
    strips: &[(EdgeId, [FaceId; 2], FilletDefinition)],
    refusal: Refusal,
    dec: &mut Decider,
) -> Refusal {
    if refusal.code == RADIUS_TOO_LARGE {
        return refusal;
    }
    let proof = (|| -> Result<Option<Option<Real>>, Refusal> {
        let mut found: Option<Option<Real>> = None;
        for (edge, faces, definition) in strips {
            let Some(incident) = uses.get(edge) else {
                continue;
            };
            if !convexity_agrees(model, *edge, incident, definition.convex, dec)? {
                continue;
            }
            for (index, face) in faces.iter().enumerate() {
                let Some(limit) = misses(model, uses, *edge, *face, index, definition, dec)? else {
                    continue;
                };
                found = Some(match (found, limit) {
                    (Some(Some(a)), Some(b)) => Some(smaller(a, b, dec)?),
                    (Some(a), None) => a,
                    (None | Some(None), b) => b,
                });
            }
        }
        Ok(found)
    })();
    match proof {
        Ok(Some(limit)) => too_large(NO_CONTACT, limit.as_ref(), dec),
        _ => refusal,
    }
}
