//! A closed plane/cylinder fillet ring whose plane rail touches countersink
//! mouths at isolated points.
//!
//! The rail and each mouth are coplanar circles outside each other at the sum
//! of their radii, so they share one point. The torus keeps its whole strip:
//! it lies above the plane, the cones below it, and the two meet only at those
//! contacts, where the trim circles are tangent and the supports are not. No
//! rolling envelope and no fitted surface takes part.
//!
//! Left whole, the plane face would pinch at every contact: its two cusps meet
//! there, and one boundary walk would pass through the vertex twice. A short
//! straight seam beside each contact cuts one cusp off as a plane cap, so every
//! loop stays a simple closed curve. The cap and the main plane patch cover the
//! original face exactly: the seam runs through the open face, and both keep
//! its origins. Around a contact the vertex link is one cycle: torus, cap,
//! cone, main plane patch, torus.
//!
//! The contact vertices are shared by the rail arcs and the split mouth circle,
//! so tangency is topology, not coincidence of coordinates. Every other point
//! of the rail lies in the open plane face, because `strip_clear` proved each
//! other boundary clear of the rail disc.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{Circle, Curve, Line, Point3, Surface, Vector3};
use catallax_topology::{BodyId, Edge, EdgeId, FaceId, Loop, Model, Sense, VertexId};

use super::fillet_envelope::{
    Ring, Turning, Uses, bad, ccw_sense, chain, checked, narrowing, other_edges, positive,
    real_sqrt, rebuild, ring, turn, unit, vertex, with, with_faces,
};
use crate::{extrude::reverse, fillet::FilletDefinition};

// Halvings of the seam before the contact is too crowded for this proof.
const SEAM_ATTEMPTS: usize = 6;

// A mouth touching the plane rail, with its cone below the plane.
struct Mouth {
    face: FaceId,
    edge: EdgeId,
    circle: Circle,
    contact: Point3,
    // Seam ends on the rail and on the mouth, counter-clockwise of the contact.
    seam: [Point3; 2],
}

// One contact's vertices.
struct Corner {
    contact: VertexId,
    rail: VertexId,
    mouth: VertexId,
}

// The cosine and sine of the angle whose half has this tangent.
fn half_angle(tangent: &Real) -> (Real, Real) {
    let square = tangent.square();
    let over = &Real::one() + &square;
    (
        &(&Real::one() - &square) / &over,
        &(&Real::from_i64(2) * tangent) / &over,
    )
}

// The seam ends at a size, when the straight seam leaves both circles
// transversally: a ray leaving a disc outward never re-enters it, so the open
// seam lies outside the rail disc and the mouth.
fn seam_ends(
    ring: &Ring,
    centre: &Point3,
    mouth: &Real,
    [toward, side]: [&Vector3; 2],
    size: &Real,
    dec: &mut Decider,
) -> Result<Option<[Point3; 2]>, Refusal> {
    let rail = &ring.centre;
    let (cos, sin) = half_angle(&(size / rail));
    let on_rail = ring.origin.plus(
        &toward
            .scaled(&(rail * &cos))
            .plus(&side.scaled(&(rail * &sin))),
    );
    let (cos, sin) = half_angle(&(size / mouth));
    let on_mouth = centre.plus(
        &toward
            .scaled(&-&(mouth * &cos))
            .plus(&side.scaled(&(mouth * &sin))),
    );
    let across = on_mouth.minus(&on_rail);
    let leaves_rail = across.dot(&on_rail.minus(&ring.origin));
    let leaves_mouth = -&across.dot(&on_mouth.minus(centre));
    Ok(
        (positive(&leaves_rail, dec)? && positive(&leaves_mouth, dec)?)
            .then_some([on_rail, on_mouth]),
    )
}

// Whether a disc about a point misses every circle of these edges. The seam
// ends lie in the disc, so the whole seam does too, and the disc touches the
// plane face's boundary only along the rail and its own mouth.
fn disc_clear(
    model: &Model,
    edges: &[EdgeId],
    at: &Point3,
    reach: &Real,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for edge in edges {
        let catallax_intersect::CurveForm::Circle {
            center, radius_sq, ..
        } = catallax_intersect::curve_form(model.geom(), model.edge(*edge).curve, dec)?
        else {
            return Ok(false);
        };
        let gap = &real_sqrt(&center.minus(at).norm_squared()) - &real_sqrt(&radius_sq);
        if !positive(&(gap.square() - reach.square()), dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

// The largest seam near half the smaller radius that stays clear of the
// plane's other boundaries.
fn seam(
    model: &Model,
    ring: &Ring,
    edge: EdgeId,
    mouth: &Circle,
    others: &[EdgeId],
    dec: &mut Decider,
) -> Result<Option<(Point3, [Point3; 2])>, Refusal> {
    let offset = mouth.center.minus(&ring.origin);
    let toward = unit(&offset);
    let side = ring.axis.cross(&toward);
    let contact = ring.origin.plus(&toward.scaled(&ring.centre));
    let radius = real_sqrt(&mouth.radius_sq);
    let smaller = if positive(&(&radius - &ring.centre), dec)? {
        ring.centre.clone()
    } else {
        radius.clone()
    };
    let mut size = &smaller / &Real::from_i64(2);
    let half = Real::from_rational(Rational::from_ratio(1, 2).expect("nonzero denominator"));
    for _ in 0..SEAM_ATTEMPTS {
        if let Some(ends) = seam_ends(ring, &mouth.center, &radius, [&toward, &side], &size, dec)? {
            let distances = ends.clone().map(|p| p.minus(&contact).norm_squared());
            let reach = real_sqrt(if positive(&(&distances[0] - &distances[1]), dec)? {
                &distances[0]
            } else {
                &distances[1]
            });
            let clear: Vec<_> = others.iter().copied().filter(|e| *e != edge).collect();
            if disc_clear(model, &clear, &contact, &reach, dec)? {
                return Ok(Some((contact, ends)));
            }
        }
        size = &size * &half;
    }
    Ok(None)
}

// The ring edge of a plane hole, when its other face is a cone narrowing below
// the plane.
fn touched_mouth(
    model: &Model,
    ring: &Ring,
    edge: EdgeId,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Option<(FaceId, Circle)>, Refusal> {
    let single = model.face(ring.plane_face).loops.iter().any(|lp| {
        matches!(model.loop_(*lp), Loop::Coedges(cs)
            if matches!(cs.as_slice(), [c] if model.coedge(*c).edge == edge))
    });
    let e = model.edge(edge);
    let (true, None, None, Curve::Circle(circle)) =
        (single, e.start, e.end, model.geom().curve(e.curve))
    else {
        return Ok(None);
    };
    let Some(cone_face) = uses
        .get(&edge)
        .and_then(|incident| incident.iter().find(|(f, _)| *f != ring.plane_face))
        .map(|(f, _)| *f)
    else {
        return Ok(None);
    };
    let Surface::Cone(cone) = model.geom().surface(model.face(cone_face).surface) else {
        return Ok(None);
    };
    Ok(narrowing(cone, circle, &ring.axis, dec)?.map(|_| (cone_face, circle.clone())))
}

/// Builds the ring with its contacts on `touched`, or returns `None` when the
/// contact is outside this stratum so the caller keeps its refusal.
#[allow(clippy::too_many_arguments)]
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    edge: EdgeId,
    faces: [FaceId; 2],
    definition: &FilletDefinition,
    touched: &[EdgeId],
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Option<BodyId>, Refusal> {
    let Some(ring) = ring(model, edge, faces, definition, dec)? else {
        return Ok(None);
    };
    let plane_edges = other_edges(model, ring.plane_face, edge)
        .ok_or_else(|| bad("a plane loop has no edges"))?;
    let mut mouths = Vec::with_capacity(touched.len());
    for mouth in touched {
        let Some((face, circle)) = touched_mouth(model, &ring, *mouth, uses, dec)? else {
            return Ok(None);
        };
        let Some((contact, seam)) = seam(model, &ring, *mouth, &circle, &plane_edges, dec)? else {
            return Ok(None);
        };
        mouths.push(Mouth {
            face,
            edge: *mouth,
            circle,
            contact,
            seam,
        });
    }

    // Disjoint mouths leave the plane connected between neighbouring contacts.
    for (i, a) in mouths.iter().enumerate() {
        for c in &mouths[i + 1..] {
            let reach = real_sqrt(&a.circle.radius_sq) + real_sqrt(&c.circle.radius_sq);
            let apart = a.circle.center.minus(&c.circle.center).norm_squared();
            if !positive(&(apart - reach.square()), dec)? {
                return Ok(None);
            }
        }
    }

    // Counter-clockwise order of the contacts about the axis.
    let reference = unit(&mouths[0].contact.minus(&ring.origin));
    let turning = Turning {
        origin: &ring.origin,
        reference: &reference,
        second: ring.axis.cross(&reference),
    };
    let mut order: Vec<usize> = (0..mouths.len()).collect();
    for i in 1..order.len() {
        let mut j = i;
        while j > 0
            && turning.before(
                &mouths[order[j]].contact,
                &mouths[order[j - 1]].contact,
                dec,
            )?
        {
            order.swap(j, j - 1);
            j -= 1;
        }
    }

    let corners: Vec<Corner> = mouths
        .iter()
        .map(|mouth| {
            let [rail, on_mouth] = mouth.seam.clone();
            Corner {
                contact: vertex(model, mouth.contact.clone()),
                rail: vertex(model, rail),
                mouth: vertex(model, on_mouth),
            }
        })
        .collect();
    let axis = &ring.axis;
    let [cylinder_contact, plane_contact] = if ring.plane_index == 0 {
        [
            definition.contacts[1].clone(),
            definition.contacts[0].clone(),
        ]
    } else {
        [
            definition.contacts[0].clone(),
            definition.contacts[1].clone(),
        ]
    };
    let rail_sense = ccw_sense(&plane_contact, axis, dec)?;
    let rail = model.geom_mut().add_curve(plane_contact);
    let cylinder_rail = model.geom_mut().add_curve(cylinder_contact);
    let plane_turn = turn(model, uses, axis, ring.plane_face, edge, dec)?;

    let arc = |model: &mut Model, curve, sense, start, end| {
        model.add_edge(Edge {
            curve,
            sense,
            start: Some(start),
            end: Some(end),
        })
    };
    // Each cap holds the rail arc from its contact to the seam; the main patch
    // holds the arc from the seam to the next contact.
    let mut cap_arcs = vec![None; mouths.len()];
    let mut main_arcs = Vec::with_capacity(mouths.len());
    for (k, &i) in order.iter().enumerate() {
        let next = order[(k + 1) % order.len()];
        let (contact, seam_foot) = (corners[i].contact, corners[i].rail);
        cap_arcs[i] = Some(arc(model, rail, rail_sense, contact, seam_foot));
        main_arcs.push(arc(
            model,
            rail,
            rail_sense,
            seam_foot,
            corners[next].contact,
        ));
    }

    let mut plane_main = Vec::new();
    let mut new_faces = Vec::new();
    let mut replaced = BTreeMap::new();
    let surface = model
        .geom()
        .surface(model.face(ring.plane_face).surface)
        .clone();
    let plane_sense = model.face(ring.plane_face).sense;
    let origins = model.face_origins(ring.plane_face).to_vec();
    for main in &main_arcs {
        plane_main.push((*main, Some(with(plane_turn))));
    }
    let mut rail_uses: Vec<(EdgeId, Option<Sense>)> = main_arcs
        .iter()
        .map(|main| (*main, Some(reverse(with(plane_turn)))))
        .collect();
    for (i, mouth) in mouths.iter().enumerate() {
        let corner = &corners[i];
        let curve = model.edge(mouth.edge).curve;
        let sense = ccw_sense(&Curve::Circle(mouth.circle.clone()), axis, dec)?;
        // Counter-clockwise the mouth reaches its seam foot just before the
        // contact, since the foot lies on the rail's counter-clockwise side.
        let short = arc(model, curve, sense, corner.mouth, corner.contact);
        let long = arc(model, curve, sense, corner.contact, corner.mouth);
        let line = model.geom_mut().add_curve(Curve::Line(Line {
            origin: mouth.seam[0].clone(),
            dir: mouth.seam[1].minus(&mouth.seam[0]),
        }));
        let cut = model.add_edge(Edge {
            curve: line,
            sense: Sense::Same,
            start: Some(corner.rail),
            end: Some(corner.mouth),
        });
        let hole_turn = turn(model, uses, axis, ring.plane_face, mouth.edge, dec)?;
        let cap_arc = cap_arcs[i].expect("every mouth has a place in the order");
        let cap = chain(
            model,
            vec![
                (cap_arc, Some(with(plane_turn))),
                (short, Some(with(hole_turn))),
                (cut, None),
            ],
        )?;
        let across = cap
            .iter()
            .find(|(e, _)| *e == cut)
            .map(|(_, s)| reverse(*s))
            .expect("the cap holds its seam");
        plane_main.push((long, Some(with(hole_turn))));
        plane_main.push((cut, Some(across)));
        rail_uses.push((cap_arc, Some(reverse(with(plane_turn)))));
        new_faces.push(crate::extrude::face(
            model,
            surface.clone(),
            plane_sense,
            vec![cap],
            None,
            origins.clone(),
        ));
        let cone_turn = turn(model, uses, axis, mouth.face, mouth.edge, dec)?;
        let cone_loop = chain(
            model,
            vec![
                (short, Some(with(cone_turn))),
                (long, Some(with(cone_turn))),
            ],
        )?;
        replaced.insert(
            mouth.face,
            rebuild(model, mouth.face, &[mouth.edge], cone_loop),
        );
    }
    let plane_main = chain(model, plane_main)?;
    let mut plane_drop = vec![edge];
    plane_drop.extend(mouths.iter().map(|mouth| mouth.edge));
    replaced.insert(
        ring.plane_face,
        rebuild(model, ring.plane_face, &plane_drop, plane_main),
    );

    // The cylinder keeps one ring, at the height of the ball centre.
    let input = model.edge(edge).clone();
    let cylinder_use = uses[&edge]
        .iter()
        .find(|(f, _)| *f == ring.cylinder_face)
        .map(|(_, s)| *s)
        .ok_or_else(|| bad("replaced boundary is not on its face"))?;
    let cylinder_ring = model.add_edge(Edge {
        curve: cylinder_rail,
        sense: input.sense,
        start: None,
        end: None,
    });
    replaced.insert(
        ring.cylinder_face,
        rebuild(
            model,
            ring.cylinder_face,
            &[edge],
            vec![(cylinder_ring, cylinder_use)],
        ),
    );
    let rail_loop = chain(model, rail_uses)?;
    new_faces.push(crate::extrude::face(
        model,
        definition.surface.clone(),
        Sense::Reversed,
        vec![rail_loop, vec![(cylinder_ring, reverse(cylinder_use))]],
        None,
        vec![definition.origin.clone()],
    ));

    let output = with_faces(model, body, ring.plane_face, &replaced, &new_faces);
    let output = model.add_body(output);
    checked(model, output).map(Some)
}
