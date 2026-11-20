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
//! A ring torus has no polar plane. Seen along its axis its generators are
//! the two equator circles, and seen across it the two meridian circles in
//! the plane through the axis and the top and bottom circles; those are exact
//! sections and take the quadric path. Every other view gives two closed
//! branches of degree 8, which `torus.rs` fits and marks approximate
//! (decision Q14), with exact trims. Spline faces and horn or spindle torus
//! sheets refuse with `app.contour-surface` rather than returning a partial
//! outline.

use catallax_boolean::face_contains;
use catallax_intersect::{
    Component, CurveSurface, Cyclic, SurfacePair, arc_point, curve_is_closed, cyclic_order,
    intersect_curve_surface, intersect_surfaces, on_curve, order_on_curve,
};
use catallax_numbers::{
    Answer, BoundKind, Decider, Mark, MarkLog, Outcome, Rational, Real, Refusal, StrictOrder,
};
use catallax_shapes::{
    Circle, Curve, CurveId, Line, Plane, Point3, PointDef, PointId, RVector3, Surface, SurfaceId,
    Vector3,
};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

use crate::ApplicationBody;
use crate::queries::{QueryRefusal, perpendicular, query, sign};

mod fit;
mod torus;

/// Code of the mark a [`Contours`] answer carries when any arc is fitted. Its
/// bound is the largest [`ContourFit::tolerance_mm`]: a measured value, not a
/// proof.
pub const CONTOUR_FIT_MARK: &str = "app.contour-fitted";

/// How far a fitted arc may lie from the true contour generator. Decision
/// Q14 marks such pieces approximate; exact arcs have none.
#[derive(Clone, Debug, PartialEq)]
pub struct ContourFit {
    /// The largest distance found between the fitted carrier and the true
    /// generator at the same parameter, plus its measurement error, in
    /// millimetres: the reported gap of rule T4.
    pub deviation_mm: f64,
    /// The tolerance rule T4 stores for it: 1.1 times the deviation, rounded
    /// up.
    pub tolerance_mm: f64,
    /// Rule T2's fit target the fitter aimed at.
    pub target_mm: f64,
    /// Cubic Hermite spans in the fitted carrier.
    pub spans: usize,
}

/// One trimmed contour generator.
#[derive(Clone, Debug)]
pub struct ContourArc {
    /// Curved face of the queried body whose generator this arc is.
    pub face: FaceId,
    /// Carrier in [`Contours::model`]: the exact face-plane intersection, or
    /// for a fitted arc a cubic B-spline over exactly the arc.
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
    /// `None` for an exact arc. A fitted arc's ends are still the exact
    /// boundary crossings; its carrier passes within the deviation of them.
    pub fit: Option<ContourFit>,
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

/// A plane whose section of a face's surface holds the face's generators.
struct Section {
    origin: Point3,
    normal: Vector3,
    /// The generator itself where the plane touches the surface along it,
    /// built directly rather than intersected.
    carrier: Option<Curve>,
}

/// Where one face's generators come from.
enum Generators {
    /// Plane sections; none for a plane or a cylinder seen along its axis.
    Sections(Vec<Section>),
    /// A ring torus seen obliquely: two fitted closed branches.
    ObliqueTorus,
}

fn sections(surface: &Surface, d: &Vector3, dec: &mut Decider) -> Result<Generators, Refusal> {
    if let Surface::Torus(t) = surface {
        // Along the axis, N.d = cos v (a.d) vanishes on the equator circles
        // rho = R +- r, the section by the plane through the centre normal to
        // the axis.
        if sign(&t.axis.cross(d).norm_squared(), dec)? == 0 {
            return Ok(Generators::Sections(vec![Section {
                origin: t.center.clone(),
                normal: t.axis.clone(),
                carrier: None,
            }]));
        }
        if sign(&t.axis.dot(d), dec)? != 0 {
            return Ok(Generators::ObliqueTorus);
        }
        // Across the axis, N.d = cos v (e(u).d): the two meridian circles in
        // the plane through the axis normal to d, and the top and bottom
        // circles v = +-pi/2, along which the plane normal to the axis at
        // height +-r touches the torus.
        let lift = t.axis.scaled(&(&t.minor / &t.axis.norm_squared().sqrt()));
        let mut sections = vec![Section {
            origin: t.center.clone(),
            normal: d.clone(),
            carrier: None,
        }];
        for side in [1, -1] {
            let center = t.center.plus(&lift.scaled(&Real::from_i64(side)));
            sections.push(Section {
                origin: center.clone(),
                normal: t.axis.clone(),
                carrier: Some(Curve::Circle(Circle {
                    center,
                    normal: t.axis.clone(),
                    radius_sq: t.major.square(),
                    ref_dir: t.ref_dir.clone(),
                })),
            });
        }
        return Ok(Generators::Sections(sections));
    }
    Ok(Generators::Sections(
        polar_plane(surface, d, dec)?
            .map(|(origin, normal)| Section {
                origin,
                normal,
                carrier: None,
            })
            .into_iter()
            .collect(),
    ))
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
                "contour generators of spline faces and horn or spindle torus sheets need a \
                 traced curve; unfinished",
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

/// The generators of one face in one plane section, trimmed to the face.
fn face_generators(
    scratch: &mut Model,
    body_faces: &[FaceId],
    face: FaceId,
    section: Section,
    arcs: &mut Vec<ContourArc>,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let Section {
        origin,
        normal,
        carrier,
    } = section;
    let ref_dir = perpendicular(&normal, dec)?;
    let plane: SurfaceId = scratch.geom_mut().add_surface(Surface::Plane(Plane {
        origin,
        normal,
        ref_dir,
    }));
    let support = scratch.face(face).surface;
    let carriers: Vec<CurveId> = if let Some(curve) = carrier {
        vec![scratch.geom_mut().add_curve(curve)]
    } else if let Surface::Cone(cone) = scratch.geom().surface(support) {
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
                    "a curved face coincides with its generator plane",
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
        // Each piece kept by the face, with the faces it lies on. Every
        // in-plane edge's ends are cuts, so a piece lies along an edge's span
        // entirely or not at all, and one interior point decides both.
        let mut labels: Vec<Option<Vec<FaceId>>> = Vec::with_capacity(pieces.len());
        for &(start, end) in &pieces {
            let interior = arc_point(scratch.geom(), carrier, start, end, dec)?;
            marks.extend(&interior.marks);
            labels.push(if face_contains(scratch, face, &interior.value, dec)? {
                Some(lies_along(
                    scratch,
                    body_faces,
                    face,
                    &interior.value,
                    &in_plane,
                    dec,
                )?)
            } else {
                None
            });
        }
        // Merge runs of neighbouring kept pieces that lie on the same faces,
        // so that `on_faces` holds along the whole arc. On a closed carrier a
        // run may wrap round; a carrier kept everywhere on one set of faces
        // is one whole closed arc.
        let n = pieces.len();
        let arc = |span: (Option<PointId>, Option<PointId>), on: &Option<Vec<FaceId>>| ContourArc {
            face,
            curve: carrier,
            start: span.0,
            end: span.1,
            on_faces: on.clone().expect("a kept piece"),
            fit: None,
        };
        if closed && labels[0].is_some() && labels.iter().all(|l| *l == labels[0]) {
            arcs.push(arc((None, None), &labels[0]));
            continue;
        }
        // Start a closed carrier where one run ends and another begins.
        let first = if closed {
            (0..n)
                .find(|&i| labels[i] != labels[(i + n - 1) % n])
                .unwrap_or(0)
        } else {
            0
        };
        let mut run: Option<(usize, usize)> = None;
        for step in 0..n {
            let i = (first + step) % n;
            if let Some((s, e)) = run
                && labels[i] != labels[s]
            {
                arcs.push(arc((pieces[s].0, pieces[e].1), &labels[s]));
                run = None;
            }
            if labels[i].is_some() {
                run = Some(run.map_or((i, i), |(s, _)| (s, i)));
            }
        }
        if let Some((s, e)) = run {
            arcs.push(arc((pieces[s].0, pieces[e].1), &labels[s]));
        }
    }
    Ok(())
}

/// The faces a kept piece lies on: its own face, and both faces of every
/// in-plane boundary edge whose span holds the piece's interior point. The
/// span decides, not the edge's carrier alone: two edges on one line can
/// border different faces.
fn lies_along(
    scratch: &mut Model,
    body_faces: &[FaceId],
    face: FaceId,
    interior: &Point3,
    in_plane: &[EdgeId],
    dec: &mut Decider,
) -> Result<Vec<FaceId>, Refusal> {
    let mut on = vec![face];
    if !in_plane.is_empty() {
        let probe = scratch
            .geom_mut()
            .add_point(PointDef::Explicit(interior.clone()));
        for &edge in in_plane {
            if !on_edge_span(scratch, edge, probe, dec)? {
                continue;
            }
            for coedge in scratch.coedges_of_edge(edge) {
                if let Some((_, other)) = scratch.owner_of_coedge(coedge)
                    && !on.contains(&other)
                    && body_faces.contains(&other)
                {
                    on.push(other);
                }
            }
        }
    }
    on.sort_by_key(|f| body_faces.iter().position(|g| g == f));
    Ok(on)
}

/// Whether `point`, which is not a vertex, lies on the span of `edge`: on its
/// carrier, and on an open carrier between its two vertices, on a closed one
/// on the arc its sense runs from its start vertex to its end vertex.
fn on_edge_span(
    scratch: &Model,
    edge: EdgeId,
    point: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let stored = scratch.edge(edge);
    let curve = stored.curve;
    if !lies_on(scratch, curve, point, dec)? {
        return Ok(false);
    }
    let (Some(start), Some(end)) = (stored.start, stored.end) else {
        // An edge without vertices is its whole closed carrier.
        return Ok(true);
    };
    let (start, end) = (scratch.vertex(start).point, scratch.vertex(end).point);
    let geom = scratch.geom();
    if curve_is_closed(geom, curve, dec)? {
        if same_point(scratch, start, end, dec)? {
            return Ok(true);
        }
        let order = match stored.sense {
            Sense::Same => [start, point, end],
            Sense::Reversed => [end, point, start],
        };
        return match cyclic_order(geom, curve, order, dec) {
            Outcome::Proven(order) => Ok(order == Cyclic::Forward),
            Outcome::ProvenCoincident(_) => Ok(false),
            Outcome::NotExact(n) => Err(Refusal::from(n)),
            Outcome::Refused(r) => Err(r),
        };
    }
    let mut order = |a: PointId, b: PointId| match order_on_curve(geom, curve, a, b, dec) {
        Outcome::Proven(order) => Ok(Some(order)),
        Outcome::ProvenCoincident(_) => Ok(None),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    };
    let (before, after) = (order(start, point)?, order(point, end)?);
    Ok(before.is_some() && before == after)
}

/// Contour generators of every curved face of `body` seen along `direction`,
/// which points toward the viewer. Only its line matters: the generators of a
/// direction and its reverse are the same. On quadrics, and on a torus seen
/// along or across its axis, each arc is exact: its carrier is the exact
/// intersection of the face's surface with a plane, and its ends are exact
/// points where the face boundary crosses that plane. On a torus seen
/// obliquely each arc is fitted: its ends are exact, its carrier is a cubic
/// B-spline within [`ContourFit::deviation_mm`] of the true generator, and the
/// answer carries a [`CONTOUR_FIT_MARK`] mark. Isolated generator points, such
/// as a cone apex seen from inside the cone's angle, draw nothing and are left
/// out.
///
/// Arcs are not split where they cross in projection and carry no
/// visibility; hiding is a separate step. A spline face or a horn or spindle
/// torus sheet refuses the whole query with `app.contour-surface`; an oblique
/// torus face bounded by anything other than its parallels and meridians
/// refuses with `app.contour-trim`.
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
        let mut jobs = Vec::new();
        for &face in &faces {
            dec.check_deadline()?;
            match sections(model.geom().surface(model.face(face).surface), &d, dec)? {
                Generators::Sections(sections) => {
                    jobs.extend(sections.into_iter().map(|s| (face, Some(s))));
                }
                Generators::ObliqueTorus => jobs.push((face, None)),
            }
        }
        let mut scratch = model.clone();
        let mut marks = model.body(body).marks.clone();
        let mut arcs = Vec::new();
        for (face, section) in jobs {
            match section {
                Some(section) => face_generators(
                    &mut scratch,
                    &faces,
                    face,
                    section,
                    &mut arcs,
                    &mut marks,
                    dec,
                )?,
                None => {
                    let Surface::Torus(t) = scratch.geom().surface(scratch.face(face).surface)
                    else {
                        unreachable!("only a ring torus is fitted");
                    };
                    let t = t.clone();
                    torus::oblique_generators(
                        &mut scratch,
                        face,
                        &t,
                        &d,
                        &mut arcs,
                        &mut marks,
                        dec,
                    )?;
                }
            }
        }
        let worst = arcs
            .iter()
            .filter_map(|arc| arc.fit.as_ref().map(|f| f.tolerance_mm))
            .reduce(f64::max);
        if let Some(worst) = worst {
            let bound = Rational::from_f64_exact(worst).ok_or_else(|| {
                Refusal::new(
                    "app.contour-fit",
                    "a fitted contour generator has a nonfinite tolerance",
                )
            })?;
            marks.push(Mark::new(
                CONTOUR_FIT_MARK,
                BoundKind::Distance,
                bound,
                "fitted contour generators: rule T4 tolerance of the largest measured \
                 deviation, not a proof",
            ));
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
