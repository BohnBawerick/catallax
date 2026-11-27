//! Plane sections of a solid, face by face.
//!
//! OpenCASCADE's section builder (`BRepAlgoAPI_Section` over
//! `BOPAlgo_Section`) intersects every face with the section face, trims each
//! curve to both faces, and adds the edges of faces that coincide with it.
//! SMLib's `SmBrep::CreatePlanarSectionCurves` sections each face's surface
//! by the infinite plane (`SmSurface::CreatePlanarSectionCurves`, as
//! `SmCutter` uses it) and trims the curves to the face.
//! Against an unbounded plane the trim to the plane is empty, so each face
//! of the solid is cut on its own: the exact intersection of its surface
//! with the plane, cut where the face boundary crosses the plane, and each
//! piece kept or dropped by exact membership of one interior point, the way
//! [`super::contour_generators`] trims a generator to its face. A planar face
//! in the plane gives its boundary edges, and any boundary edge lying in the
//! plane is reported once, under its own id; a cut along such an edge is
//! that edge. Only the faces themselves are touched, never a Boolean of the
//! whole body.

use std::collections::BTreeSet;

use catallax_boolean::{BoolOp, boolean, face_contains_point};
use catallax_intersect::{Component, SurfacePair, intersect_surfaces};
use catallax_numbers::{Answer, Answered, Decider, Rational, Refusal};
use catallax_shapes::{Plane, PointId, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense, block_in_frame};

use super::{boundary_cuts, carrier_pieces, on_span, runs, same_point};
use crate::queries::{SectionArc, SectionCurves, SectionPlane, SectionSolid, perpendicular, sign};

pub(crate) fn cut_solid(
    model: &Model,
    body: BodyId,
    plane: &SectionPlane,
    dec: &mut Decider,
) -> Answered<SectionSolid> {
    let bounds = catallax_measure::bounding_box(model, body, &Default::default(), dec)?;
    let mut scratch = model.clone();
    let Some(bounds) = bounds.value else {
        return Ok(Answer::with_marks(
            SectionSolid {
                model: scratch,
                body,
            },
            model.body(body).marks.clone(),
        ));
    };
    let o = &plane.origin;
    // This L1 radius strictly contains the whole body. Each frame vector
    // has length at least one, so only the near face of the tool can cut it.
    let reach = [&o.x, &o.y, &o.z]
        .iter()
        .enumerate()
        .fold(Rational::one(), |sum, (i, &coord)| {
            sum + (&bounds.min[i] - coord)
                .abs()
                .max((&bounds.max[i] - coord).abs())
        });
    let rescale = |v: RVector3| {
        let component = v.x.abs().max(v.y.abs()).max(v.z.abs());
        v.scaled(&(Rational::one() / component))
    };
    let n = rescale(plane.normal.clone());
    let axis = if n.x.is_zero() && n.y.is_zero() {
        RVector3::new(1.into(), 0.into(), 0.into())
    } else {
        RVector3::new(0.into(), 0.into(), 1.into())
    };
    // The plane does not need an orthonormal frame. Keeping the box rational
    // avoids introducing square roots before any surfaces are intersected.
    let u = rescale(axis.cross(&n));
    let v = rescale(n.cross(&u));
    let a = u.scaled(&reach);
    let b = v.scaled(&reach);
    let corner = o.plus(&RVector3::zero().minus(&a).minus(&b));
    let tool = block_in_frame(
        &mut scratch,
        &corner,
        &[a.plus(&a), b.plus(&b), n.scaled(&reach)],
    )?;
    let cut = boolean(&mut scratch, body, tool, BoolOp::Subtract, dec).map_err(|e| e.refusal)?;
    Ok(Answer::with_marks(
        SectionSolid {
            model: scratch,
            body: cut.body,
        },
        cut.marks,
    ))
}

/// The section of `body` by the plane through `origin` with `normal`.
pub(crate) fn plane_section(
    model: &Model,
    body: BodyId,
    origin: &RPoint3,
    normal: &RVector3,
    dec: &mut Decider,
) -> Answered<SectionCurves> {
    let n = normal.to_real();
    let o = origin.to_real();
    let faces = model.faces_of_body(body);
    let mut scratch = model.clone();
    let mut marks = model.body(body).marks.clone();
    let ref_dir = perpendicular(&n, dec)?;
    let plane = scratch.geom_mut().add_surface(Surface::Plane(Plane {
        origin: o.clone(),
        normal: n.clone(),
        ref_dir,
    }));
    let mut in_plane: BTreeSet<EdgeId> = BTreeSet::new();
    let mut arcs: Vec<SectionArc> = Vec::new();
    let mut points: Vec<PointId> = Vec::new();
    for &face in &faces {
        dec.check_deadline()?;
        let support = scratch.face(face).surface;
        if let Surface::Plane(p) = scratch.geom().surface(support)
            && sign(&p.normal.cross(&n).norm_squared(), dec)? == 0
            && sign(&o.minus(&p.origin).dot(&n), dec)? == 0
        {
            for &lp in &scratch.face(face).loops {
                if let Loop::Coedges(coedges) = scratch.loop_(lp) {
                    in_plane.extend(coedges.iter().map(|&c| scratch.coedge(c).edge));
                }
            }
            continue;
        }
        let pair = intersect_surfaces(scratch.geom_mut(), support, plane, dec)?;
        marks.extend(&pair.marks);
        let components = match pair.value {
            SurfacePair::Disjoint => continue,
            SurfacePair::Coincident { .. } => {
                return Err(Refusal::new(
                    "app.section-coincident",
                    "a curved face coincides with the section plane",
                ));
            }
            SurfacePair::Meets(components) => components,
        };
        let boundary = boundary_cuts(&mut scratch, face, plane, &mut marks, dec)?;
        in_plane.extend(&boundary.in_plane);
        for component in components {
            match component {
                Component::Point { point } => {
                    if face_contains_point(&mut scratch, face, point, dec)? {
                        points.push(point);
                    }
                }
                Component::Curve {
                    curve, singular, ..
                } => {
                    let pieces = carrier_pieces(
                        &mut scratch,
                        &faces,
                        face,
                        curve,
                        &boundary,
                        &singular,
                        &mut marks,
                        dec,
                    )?;
                    // A piece along an edge in the plane is that edge.
                    let keys = pieces
                        .labels
                        .iter()
                        .map(|l| {
                            l.as_ref()
                                .filter(|along| along.edge.is_none())
                                .map(|along| along.faces.clone())
                        })
                        .collect();
                    for ((start, end), faces) in runs(&pieces.spans, keys, pieces.closed) {
                        arcs.push(SectionArc {
                            curve,
                            start,
                            end,
                            faces,
                            edge: None,
                        });
                    }
                }
            }
        }
    }
    let mut edges = Vec::new();
    for edge in model.edges_of_body(body) {
        if !in_plane.contains(&edge) {
            continue;
        }
        let stored = scratch.edge(edge);
        let (start, end) = match (stored.start, stored.end) {
            (Some(a), Some(b)) => {
                let (a, b) = (scratch.vertex(a).point, scratch.vertex(b).point);
                match stored.sense {
                    Sense::Same => (Some(a), Some(b)),
                    Sense::Reversed => (Some(b), Some(a)),
                }
            }
            _ => (None, None),
        };
        let mut on: Vec<FaceId> = Vec::new();
        for coedge in scratch.coedges_of_edge(edge) {
            if let Some((_, face)) = scratch.owner_of_coedge(coedge)
                && faces.contains(&face)
                && !on.contains(&face)
            {
                on.push(face);
            }
        }
        on.sort_by_key(|f| faces.iter().position(|g| g == f));
        edges.push(SectionArc {
            curve: stored.curve,
            start,
            end,
            faces: on,
            edge: Some(edge),
        });
    }
    edges.append(&mut arcs);
    let curves = edges;
    // Vertices that touch the plane off every curve are isolated contacts
    // too, as are tangent points; each is kept once.
    let mut touching: Vec<PointId> = Vec::new();
    for vertex in model.vertices_of_body(body) {
        let point = scratch.vertex(vertex).point;
        if sign(&scratch.geom().position(point).minus(&o).dot(&n), dec)? == 0 {
            touching.push(point);
        }
    }
    touching.extend(std::mem::take(&mut points));
    'candidates: for point in touching {
        dec.check_deadline()?;
        for &kept in &points {
            if same_point(&scratch, point, kept, dec)? {
                continue 'candidates;
            }
        }
        for arc in &curves {
            for end in [arc.start, arc.end].into_iter().flatten() {
                if same_point(&scratch, point, end, dec)? {
                    continue 'candidates;
                }
            }
            let ends = arc.start.zip(arc.end).map(|(a, b)| [a, b]);
            if on_span(&scratch, arc.curve, ends, point, dec)? {
                continue 'candidates;
            }
        }
        points.push(point);
    }
    Ok(Answer::with_marks(
        SectionCurves {
            model: scratch,
            curves,
            points,
        },
        marks,
    ))
}
