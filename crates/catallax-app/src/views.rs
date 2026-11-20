//! Drawing views from the kernel: contour generators of curved faces under
//! parallel projection, independent of the display mesh.
//!
//! A contour generator (a silhouette before hiding) is the set of face points
//! whose surface normal is perpendicular to the view direction. On a quadric
//! it is the quadric's section by the polar plane of the view direction's point
//! at infinity: for `x^T A x + 2 b^T x + c = 0`, the points where
//! `d^T (A x + b) = 0` (Salmon, "A Treatise on the Analytic Geometry of Three
//! Dimensions", polar planes; Hartley and Zisserman, "Multiple View Geometry",
//! section 8.3, the contour generator of a quadric). OpenCASCADE's
//! `Contap_ContAna` gives the same lines and circles in closed form for
//! cylinders, cones and spheres, and `Contap_Contour` trims them to the face;
//! SMLib's `SmSurfaceSilhouette` traces them instead.
//!
//! Catallax follows OpenCASCADE: the exact surface-plane intersection, or for
//! a cone the same closed form built directly, gives the carriers. The face's
//! boundary edges cut them where they cross the polar plane, and exact face
//! membership of one interior point keeps or drops each piece. Only the curved
//! face itself is touched, never the rest of the body.
//!
//! Tori and spline faces have no polar plane. Their generators need a fitted
//! curve, which is not built yet; those faces refuse with
//! `app.contour-surface` rather than returning a partial outline.

use catallax_boolean::face_contains;
use catallax_intersect::{
    Component, CurveSurface, Cyclic, SurfacePair, arc_point, curve_is_closed, cyclic_order,
    intersect_curve_surface, intersect_surfaces, on_curve, order_on_curve,
};
use catallax_numbers::{Answer, Decider, MarkLog, Outcome, Real, Refusal, StrictOrder};
use catallax_shapes::{
    Curve, CurveId, Line, Plane, Point3, PointDef, PointId, RVector3, Surface, SurfaceId, Vector3,
};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model};

use crate::ApplicationBody;
use crate::queries::{QueryRefusal, perpendicular, query, sign};

/// One trimmed contour generator.
#[derive(Clone, Debug)]
pub struct ContourArc {
    /// Curved face of the queried body whose generator this arc is.
    pub face: FaceId,
    /// Carrier in [`Contours::model`]: the exact face-plane intersection.
    pub curve: CurveId,
    /// Start of the arc; `None` with `end` for a whole closed curve.
    pub start: Option<PointId>,
    /// End of the arc; `None` with `start` for a whole closed curve.
    pub end: Option<PointId>,
    /// Every face of the queried body the arc lies on, in `faces_of_body`
    /// order, `face` included. More than one face means the generator runs
    /// along an existing edge, such as the edge between a cylinder and a
    /// tangent plane.
    pub on_faces: Vec<FaceId>,
}

/// Contour generators of every curved face of one body for one direction.
#[derive(Clone, Debug)]
pub struct Contours {
    /// Owns every curve and point handle below. The input model is untouched,
    /// and face handles in the arcs are the input body's.
    pub model: Model,
    /// Arcs in face order, then carrier order, then along each carrier.
    pub arcs: Vec<ContourArc>,
}

/// The polar plane of the point at infinity in direction `d` for the quadric
/// carrying `surface`, as a point and an unnormalised normal `A d`. `None`
/// when the face has no generators for this direction: a plane, or a
/// cylinder seen along its axis.
fn polar_plane(
    surface: &Surface,
    d: &Vector3,
    dec: &mut Decider,
) -> Result<Option<(Point3, Vector3)>, Refusal> {
    let nonzero = |normal: Vector3, dec: &mut Decider| -> Result<Option<Vector3>, Refusal> {
        Ok((sign(&normal.norm_squared(), dec)? != 0).then_some(normal))
    };
    Ok(match surface {
        Surface::Plane(_) => None,
        // A = I - a a^T / (a.a) about the axis point; A d is d's part
        // perpendicular to the axis, scaled by a.a.
        Surface::Cylinder(c) => nonzero(
            d.scaled(&c.axis.norm_squared())
                .minus(&c.axis.scaled(&c.axis.dot(d))),
            dec,
        )?
        .map(|normal| (c.origin.clone(), normal)),
        // About the apex, the nappe satisfies |x|^2 = (1 + k^2) (x.a)^2 / (a.a)
        // with k the radius gained per unit axial distance, so
        // A d = (a.a) d - (1 + k^2) (a.d) a, never zero for k > 0.
        Surface::Cone(c) => {
            let steep = Real::one() + c.slope.square();
            let normal = d
                .scaled(&c.axis.norm_squared())
                .minus(&c.axis.scaled(&(&steep * &c.axis.dot(d))));
            Some((c.apex.clone(), normal))
        }
        // A = I about the centre: the great circle facing the viewer.
        Surface::Sphere(s) => Some((s.center.clone(), d.clone())),
        // With m the major direction and w = a x m the minor direction,
        // A d = (d.m) / ((m.m) major^2) m + (d.w) / ((w.w) minor^2) w.
        Surface::EllipticCylinder(e) => {
            let minor_dir = e.axis.cross(&e.major_dir);
            let major = e.major_dir.dot(d) / (&e.major_dir.norm_squared() * &e.major_sq);
            let minor = minor_dir.dot(d) / (&minor_dir.norm_squared() * &e.minor_sq);
            nonzero(
                e.major_dir.scaled(&major).plus(&minor_dir.scaled(&minor)),
                dec,
            )?
            .map(|normal| (e.origin.clone(), normal))
        }
        _ => {
            return Err(Refusal::new(
                "app.contour-surface",
                "contour generators of tori and spline faces need a fitted curve; unfinished",
            ));
        }
    })
}

fn same_point(model: &Model, a: PointId, b: PointId, dec: &mut Decider) -> Result<bool, Refusal> {
    let (a, b) = (model.geom().position(a), model.geom().position(b));
    Ok(sign(&a.minus(&b).norm_squared(), dec)? == 0)
}

fn lies_on(
    model: &Model,
    curve: CurveId,
    point: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    match on_curve(model.geom(), curve, point, dec) {
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::Proven(_) => Ok(false),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    }
}

/// Sorts distinct points along a carrier: by parameter on an open curve, by
/// cyclic order from the first point on a closed one. The same insertion
/// order as the Boolean engine's cut sorting.
fn sort_along(
    model: &Model,
    curve: CurveId,
    closed: bool,
    points: &mut [PointId],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let anchor = points.first().copied();
    for i in 1..points.len() {
        let mut j = i;
        while j > 0 {
            let before = if !closed {
                match order_on_curve(model.geom(), curve, points[j], points[j - 1], dec) {
                    Outcome::Proven(order) => order == StrictOrder::Less,
                    Outcome::ProvenCoincident(_) => {
                        return Err(Refusal::new(
                            "app.contour-cuts",
                            "distinct generator cuts were not ordered",
                        ));
                    }
                    Outcome::NotExact(n) => return Err(Refusal::from(n)),
                    Outcome::Refused(r) => return Err(r),
                }
            } else if j == 1 {
                false
            } else {
                let anchor = anchor.expect("nonempty");
                match cyclic_order(model.geom(), curve, [anchor, points[j], points[j - 1]], dec) {
                    Outcome::Proven(order) => order == Cyclic::Forward,
                    Outcome::ProvenCoincident(_) => {
                        return Err(Refusal::new(
                            "app.contour-cuts",
                            "distinct generator cuts were not ordered",
                        ));
                    }
                    Outcome::NotExact(n) => return Err(Refusal::from(n)),
                    Outcome::Refused(r) => return Err(r),
                }
            };
            if !before {
                break;
            }
            points.swap(j, j - 1);
            j -= 1;
        }
    }
    Ok(())
}

fn normal_of(model: &Model, plane: SurfaceId) -> Vector3 {
    match model.geom().surface(plane) {
        Surface::Plane(p) => p.normal.clone(),
        _ => unreachable!("the polar plane is a plane"),
    }
}

/// The cone's generators in its polar plane through the apex: lines through
/// the apex along `p +- f t`, with `p` the axis projected into the plane,
/// `t = n x a` and `f^2 = -((a.a) - (1 + k^2)(p.p)) / ((a.a)(n.n))`. This is the
/// closed form OpenCASCADE's `Contap_ContAna` uses for cones and the one
/// Catallax's plane-cone intersection uses for algebraic inputs
/// (`catallax_intersect` `pair/conics.rs`, `plane_cone`). That intersection
/// refuses irrational generator directions for rational inputs, which an
/// oblique view of a rational cone gives, so the lines are built here.
fn cone_generators(
    scratch: &mut Model,
    cone: &catallax_shapes::Cone,
    normal: &Vector3,
    dec: &mut Decider,
) -> Result<Vec<CurveId>, Refusal> {
    let axis2 = cone.axis.norm_squared();
    let normal2 = normal.norm_squared();
    let projected = cone
        .axis
        .minus(&normal.scaled(&(cone.axis.dot(normal) / normal2.clone())));
    let transverse = normal.cross(&cone.axis);
    let steep = Real::one() + cone.slope.square();
    let quadratic = &axis2 - &(&steep * &projected.norm_squared());
    let directions = match sign(&quadratic, dec)? {
        // The polar plane meets the nappe only at the apex.
        1 => return Ok(Vec::new()),
        // The view runs along one generator, where the plane is tangent.
        0 => vec![projected],
        _ => {
            let factor = (-&quadratic / (&axis2 * &normal2)).sqrt();
            vec![
                projected.plus(&transverse.scaled(&factor)),
                projected.minus(&transverse.scaled(&factor)),
            ]
        }
    };
    Ok(directions
        .into_iter()
        .map(|dir| {
            scratch.geom_mut().add_curve(Curve::Line(Line {
                origin: cone.apex.clone(),
                dir,
            }))
        })
        .collect())
}

/// The generators of one quadric face, trimmed to the face.
#[allow(clippy::too_many_arguments)]
fn face_generators(
    scratch: &mut Model,
    body_faces: &[FaceId],
    face: FaceId,
    origin: Point3,
    normal: Vector3,
    arcs: &mut Vec<ContourArc>,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let ref_dir = perpendicular(&normal, dec)?;
    let plane: SurfaceId = scratch.geom_mut().add_surface(Surface::Plane(Plane {
        origin,
        normal,
        ref_dir,
    }));
    let support = scratch.face(face).surface;
    let carriers: Vec<CurveId> = if let Surface::Cone(cone) = scratch.geom().surface(support) {
        let cone = cone.clone();
        cone_generators(scratch, &cone, &normal_of(scratch, plane), dec)?
    } else {
        let pair = intersect_surfaces(scratch.geom_mut(), support, plane, dec)?;
        marks.extend(&pair.marks);
        match pair.value {
            SurfacePair::Disjoint => return Ok(()),
            SurfacePair::Coincident { .. } => {
                return Err(Refusal::new(
                    "app.contour-degenerate",
                    "a quadric face coincides with its polar plane",
                ));
            }
            SurfacePair::Meets(components) => components
                .into_iter()
                .filter_map(|component| match component {
                    Component::Curve { curve, .. } => Some(curve),
                    Component::Point { .. } => None,
                })
                .collect(),
        }
    };
    if carriers.is_empty() {
        return Ok(());
    }
    // Every place a boundary carrier crosses the polar plane. Points off the
    // edge's own span only split a carrier further; the membership test of
    // each piece decides, and kept neighbours are merged again.
    let mut cuts: Vec<PointId> = Vec::new();
    let mut in_plane: Vec<EdgeId> = Vec::new();
    for lp in scratch.face(face).loops.clone() {
        match scratch.loop_(lp).clone() {
            Loop::Vertex(vertex) => cuts.push(scratch.vertex(vertex).point),
            Loop::Coedges(coedges) => {
                for coedge in coedges {
                    let edge = scratch.coedge(coedge).edge;
                    let stored = scratch.edge(edge).clone();
                    let hits =
                        intersect_curve_surface(scratch.geom_mut(), stored.curve, plane, dec)?;
                    marks.extend(&hits.marks);
                    match hits.value {
                        CurveSurface::Disjoint => {}
                        CurveSurface::Contained => {
                            in_plane.push(edge);
                            for vertex in [stored.start, stored.end].into_iter().flatten() {
                                cuts.push(scratch.vertex(vertex).point);
                            }
                        }
                        CurveSurface::Meets(hits) => cuts.extend(hits.into_iter().map(|h| h.point)),
                    }
                }
            }
        }
    }
    for carrier in carriers {
        dec.check_deadline()?;
        let closed = curve_is_closed(scratch.geom(), carrier, dec)?;
        let mut on: Vec<PointId> = Vec::new();
        for &cut in &cuts {
            if !lies_on(scratch, carrier, cut, dec)? {
                continue;
            }
            let mut seen = false;
            for &kept in &on {
                if same_point(scratch, cut, kept, dec)? {
                    seen = true;
                    break;
                }
            }
            if !seen {
                on.push(cut);
            }
        }
        sort_along(scratch, carrier, closed, &mut on, dec)?;
        // Pieces between consecutive cuts. An open carrier is unbounded, so
        // the face, being finite, holds none of its two outer pieces.
        let pieces: Vec<(Option<PointId>, Option<PointId>)> = if closed && on.len() < 2 {
            vec![(None, None)]
        } else if closed {
            (0..on.len())
                .map(|i| (Some(on[i]), Some(on[(i + 1) % on.len()])))
                .collect()
        } else {
            on.windows(2).map(|w| (Some(w[0]), Some(w[1]))).collect()
        };
        let mut inside = Vec::with_capacity(pieces.len());
        let mut interiors = Vec::with_capacity(pieces.len());
        for &(start, end) in &pieces {
            let interior = arc_point(scratch.geom(), carrier, start, end, dec)?;
            marks.extend(&interior.marks);
            inside.push(face_contains(scratch, face, &interior.value, dec)?);
            interiors.push(interior.value);
        }
        // Merge runs of kept pieces. On a closed carrier a run may wrap
        // round; a carrier kept everywhere is one whole closed arc.
        let n = pieces.len();
        if inside.iter().all(|&kept| kept) && closed {
            push_arc(
                scratch,
                body_faces,
                face,
                carrier,
                (None, None),
                &interiors[0],
                &in_plane,
                arcs,
                dec,
            )?;
            continue;
        }
        let first = if closed {
            inside
                .iter()
                .position(|&kept| !kept)
                .map_or(0, |i| (i + 1) % n)
        } else {
            0
        };
        let mut run: Option<(usize, usize)> = None;
        for step in 0..n {
            let i = (first + step) % n;
            if inside[i] {
                run = Some(run.map_or((i, i), |(s, _)| (s, i)));
            }
            if (!inside[i] || step + 1 == n)
                && let Some((s, e)) = run.take()
            {
                let span = (pieces[s].0, pieces[e].1);
                push_arc(
                    scratch,
                    body_faces,
                    face,
                    carrier,
                    span,
                    &interiors[s],
                    &in_plane,
                    arcs,
                    dec,
                )?;
            }
        }
    }
    Ok(())
}

/// Records one kept arc and the faces it lies on: its own face, and the far
/// face of any in-plane boundary edge that carries its interior point.
#[allow(clippy::too_many_arguments)]
fn push_arc(
    scratch: &mut Model,
    body_faces: &[FaceId],
    face: FaceId,
    curve: CurveId,
    (start, end): (Option<PointId>, Option<PointId>),
    interior: &Point3,
    in_plane: &[EdgeId],
    arcs: &mut Vec<ContourArc>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut on = vec![face];
    if !in_plane.is_empty() {
        let probe = scratch
            .geom_mut()
            .add_point(PointDef::Explicit(interior.clone()));
        for &edge in in_plane {
            if !lies_on(scratch, scratch.edge(edge).curve, probe, dec)? {
                continue;
            }
            for coedge in scratch.coedges_of_edge(edge) {
                if let Some((_, other)) = scratch.owner_of_coedge(coedge)
                    && !on.contains(&other)
                    && body_faces.contains(&other)
                    && face_contains(scratch, other, interior, dec)?
                {
                    on.push(other);
                }
            }
        }
    }
    on.sort_by_key(|f| body_faces.iter().position(|g| g == f));
    arcs.push(ContourArc {
        face,
        curve,
        start,
        end,
        on_faces: on,
    });
    Ok(())
}

/// Contour generators of every curved face of `body` seen along `direction`,
/// which points toward the viewer. Only its line matters: the generators of a
/// direction and its reverse are the same. Each arc is exact: its carrier is
/// the exact intersection of the face's quadric with the polar plane, and its
/// ends are exact points where the face boundary crosses that plane.
/// Isolated generator points, such as a cone apex seen from inside the cone's
/// angle, draw nothing and are left out.
///
/// Arcs are not split where they cross in projection and carry no
/// visibility; hiding is a separate step. A torus or spline face refuses the
/// whole query with `app.contour-surface`.
pub fn contour_generators(
    model: &Model,
    body: BodyId,
    direction: &RVector3,
    dec: &mut Decider,
) -> Result<Answer<Contours>, QueryRefusal> {
    query(model, &[body], || {
        dec.check_deadline()?;
        if direction.is_zero() {
            return Err(Refusal::new(
                "app.contour-direction",
                "view direction is zero",
            ));
        }
        let d = direction.to_real();
        let faces = model.faces_of_body(body);
        let mut polar = Vec::new();
        for &face in &faces {
            dec.check_deadline()?;
            if let Some((origin, normal)) =
                polar_plane(model.geom().surface(model.face(face).surface), &d, dec)?
            {
                polar.push((face, origin, normal));
            }
        }
        let mut scratch = model.clone();
        let mut marks = model.body(body).marks.clone();
        let mut arcs = Vec::new();
        for (face, origin, normal) in polar {
            face_generators(
                &mut scratch,
                &faces,
                face,
                origin,
                normal,
                &mut arcs,
                &mut marks,
                dec,
            )?;
        }
        Ok(Answer::with_marks(
            Contours {
                model: scratch,
                arcs,
            },
            marks,
        ))
    })
}

impl ApplicationBody {
    /// Contour generators of this result's curved faces for a direction
    /// toward the viewer.
    pub fn contour_generators(
        &self,
        direction: &RVector3,
        dec: &mut Decider,
    ) -> Result<Answer<Contours>, QueryRefusal> {
        contour_generators(&self.built.model, self.built.result, direction, dec)
    }
}

#[cfg(test)]
mod tests;
