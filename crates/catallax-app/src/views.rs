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
//! (decision Q14), with exact trims.
//!
//! On an extruded Bezier profile the generators are exact rulings; on a
//! complete rational patch, or a lofted face covered by such patches, they
//! are traced and fitted (`spline.rs`). Other spline bands refuse with
//! `app.contour-trim`, and horn or spindle torus sheets with
//! `app.contour-surface`, rather than return a partial outline.
//!
//! [`hidden_line_view`](crate::views::hidden_line_view) cuts the edges and generators where visibility can
//! change and marks each piece visible or hidden (`hidden.rs`).
//!
//! [`drawing_curve`](crate::views::drawing_curve) and its companions turn edges, face loops, hidden-line
//! strokes and section pieces into the binary64 lines, arcs, circles and
//! NURBS a drawing takes (`drawing.rs`). Plane sections, cut face by face
//! with the same trimming as the generators, are in `section.rs`.

use catallax_boolean::face_contains;
use catallax_intersect::{
    Component, CurveForm, CurveSurface, Cyclic, SurfacePair, arc_point, curve_form,
    curve_is_closed, cyclic_order, intersect_curve_surface, intersect_surfaces, on_curve,
    order_on_curve,
};
use catallax_numbers::{
    Answer, BoundKind, Decider, Mark, MarkLog, Outcome, Rational, Real, Refusal, StrictOrder,
};
use catallax_shapes::{
    Circle, Curve, CurveId, Line, Plane, Point3, PointDef, PointId, RVector3, Surface, SurfaceId,
    Vector3,
};
use catallax_topology::numerical::tolerance;
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

use crate::ApplicationBody;
use crate::queries::{QueryRefusal, drawing_query, perpendicular, sign};

mod drawing;
mod fit;
mod float;
mod hidden;
pub(crate) mod section;
mod spline;
mod torus;

pub use drawing::{
    DRAWING_FIT_MARK, DrawingCurve, DrawnCurve, DrawnEdge, drawing_curve, edge_drawing,
    face_drawing_loops, stroke_drawing,
};
pub use hidden::{
    HiddenLineView, StrokeSource, VIEW_RAY_MARK, VIEW_SHORT_PIECE_MARK, VIEW_SPLIT_MARK,
    ViewStroke, hidden_line_view,
};

/// Code of the mark a [`Contours`] answer carries when any arc is fitted. Its
/// bound is the largest [`ContourFit::tolerance_mm`]: a measured value, not a
/// proof.
pub const CONTOUR_FIT_MARK: &str = "app.contour-fitted";

/// How far a fitted arc may lie from the true contour generator, or a
/// fitted drawing curve from the true curve. Decision Q14 marks such pieces
/// approximate; exact arcs have none.
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
    /// `None` for an exact arc. On a torus a fitted arc's ends are the
    /// exact boundary crossings, and its carrier passes within the deviation
    /// of them. On a traced spline face its ends are where the fitted carrier
    /// starts and ends, within the deviation of the face boundary; where a
    /// generator crosses a chart seam of a lofted face, the arcs on the two
    /// charts share one end point there.
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

/// One face's work, in face order.
enum Job {
    Section(Section),
    ObliqueTorus,
    Spline,
}

/// Where one face's generators come from.
enum Generators {
    /// Plane sections; none for a plane or a cylinder seen along its axis.
    Sections(Vec<Section>),
    /// A ring torus seen obliquely: two fitted closed branches.
    ObliqueTorus,
    /// A spline face: exact rulings or traced generators.
    Spline,
}

fn sections(surface: &Surface, d: &Vector3, dec: &mut Decider) -> Result<Generators, Refusal> {
    if let Surface::BSpline(_) = surface {
        return Ok(Generators::Spline);
    }
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
                "contour generators of horn or spindle torus sheets need a traced curve; \
                 unfinished",
            ));
        }
    })
}

pub(crate) fn same_point(
    model: &Model,
    a: PointId,
    b: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
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
    let boundary = boundary_cuts(scratch, face, plane, marks, dec)?;
    for carrier in carriers {
        dec.check_deadline()?;
        let pieces = carrier_pieces(
            scratch,
            body_faces,
            face,
            carrier,
            &boundary,
            &[],
            marks,
            dec,
        )?;
        let keys = pieces
            .labels
            .iter()
            .map(|l| l.as_ref().map(|along| along.faces.clone()))
            .collect();
        for (span, on_faces) in runs(&pieces.spans, keys, pieces.closed) {
            arcs.push(ContourArc {
                face,
                curve: carrier,
                start: span.0,
                end: span.1,
                on_faces,
                fit: None,
            });
        }
    }
    Ok(())
}

/// Where a face's boundary meets a plane.
pub(crate) struct Boundary {
    /// Every place a boundary carrier crosses the plane, and the ends of
    /// boundary edges in it. Points off the edge's own span only split a
    /// carrier further; the membership test of each piece decides, and kept
    /// neighbours are merged again.
    cuts: Vec<PointId>,
    /// Boundary edges whose carriers lie in the plane.
    pub(crate) in_plane: Vec<EdgeId>,
}

/// The boundary of `face` against `plane`.
pub(crate) fn boundary_cuts(
    scratch: &mut Model,
    face: FaceId,
    plane: SurfaceId,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Boundary, Refusal> {
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
    Ok(Boundary { cuts, in_plane })
}

/// What a kept piece lies on.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Along {
    /// Its own face, and both faces of every in-plane boundary edge whose
    /// span holds it, in `faces_of_body` order.
    pub(crate) faces: Vec<FaceId>,
    /// The first in-plane boundary edge whose span holds it.
    pub(crate) edge: Option<EdgeId>,
}

/// The ends of a piece of a carrier: both `None` for the whole closed
/// carrier.
pub(crate) type Ends = (Option<PointId>, Option<PointId>);

/// One carrier cut where the face boundary crosses it.
pub(crate) struct Pieces {
    /// Pieces between consecutive cuts, in order along the carrier.
    pub(crate) spans: Vec<Ends>,
    /// For each piece, what it lies on when the face holds it.
    pub(crate) labels: Vec<Option<Along>>,
    /// The carrier is closed.
    pub(crate) closed: bool,
}

/// Cuts `carrier`, which lies in the plane of `boundary`, at the boundary's
/// cuts and at `extra` points of its own, and labels the pieces `face`
/// holds.
#[allow(clippy::too_many_arguments)]
pub(crate) fn carrier_pieces(
    scratch: &mut Model,
    body_faces: &[FaceId],
    face: FaceId,
    carrier: CurveId,
    boundary: &Boundary,
    extra: &[PointId],
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Pieces, Refusal> {
    let closed = curve_is_closed(scratch.geom(), carrier, dec)?;
    let mut on: Vec<PointId> = Vec::new();
    for &cut in boundary.cuts.iter().chain(extra) {
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
    let spans: Vec<Ends> = if closed && on.len() < 2 {
        vec![(None, None)]
    } else if closed {
        (0..on.len())
            .map(|i| (Some(on[i]), Some(on[(i + 1) % on.len()])))
            .collect()
    } else {
        on.windows(2).map(|w| (Some(w[0]), Some(w[1]))).collect()
    };
    // Each piece kept by the face, with what it lies on. Every in-plane
    // edge's ends are cuts, so a piece lies along an edge's span entirely
    // or not at all, and one interior point decides both.
    let mut labels: Vec<Option<Along>> = Vec::with_capacity(spans.len());
    for &span in &spans {
        let interior = piece_point(scratch, carrier, closed, span, marks, dec)?;
        labels.push(if face_contains(scratch, face, &interior, dec)? {
            Some(lies_along(
                scratch,
                body_faces,
                face,
                &interior,
                &boundary.in_plane,
                dec,
            )?)
        } else {
            None
        });
    }
    Ok(Pieces {
        spans,
        labels,
        closed,
    })
}

/// An exact point strictly inside the piece of `carrier` between `ends`,
/// whose face membership labels the piece. Cuts bound every piece, so any
/// point strictly inside gives the same label; only the work differs.
///
/// The intersection crate's arc point halves the piece in the carrier's
/// chart parameter, so between two algebraic cuts, such as the places an
/// intersection branch crosses the plane, its coordinates nest the cuts'
/// radicals, and exact face membership of such a point can cost tens of
/// millions of work units. On a line, circle or ellipse the carrier's point
/// at the middle binary64 parameter keeps short numbers instead: rational
/// on a rational line, and on a conic the half-angle point of
/// `float::conic_point`, whose axes carry one square root each. Exact
/// order confirms that it lies strictly inside the piece, as `hidden.rs`
/// confirms its sample points; failing that, and on other carriers, the arc
/// point decides.
fn piece_point(
    scratch: &mut Model,
    carrier: CurveId,
    closed: bool,
    (start, end): Ends,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    if let (Some(a), Some(b)) = (start, end) {
        for candidate in middles(scratch, carrier, a, b, dec)? {
            let c = scratch
                .geom_mut()
                .add_point(PointDef::Explicit(candidate.clone()));
            if strictly_inside(scratch, carrier, closed, [a, c, b], dec)? {
                return Ok(candidate);
            }
        }
    }
    let interior = arc_point(scratch.geom(), carrier, start, end, dec)?;
    marks.extend(&interior.marks);
    Ok(interior.value)
}

/// Exact points of a line, circle or ellipse carrier at the middle binary64
/// parameter between `a` and `b` in its closed form's direction; on a conic
/// also the middle of the rest of the conic, half a turn on, since a carrier
/// stored as an intersection component may run against its closed form.
/// None for other carriers.
fn middles(
    scratch: &Model,
    carrier: CurveId,
    a: PointId,
    b: PointId,
    dec: &mut Decider,
) -> Result<Vec<Point3>, Refusal> {
    // Rule T1's smallest floor: binary64 images only place the candidate,
    // and exact order then confirms it.
    let floor = tolerance::floor_minimum_mm();
    let image = |p: &Point3, dec: &mut Decider| -> Result<[f64; 3], Refusal> {
        Ok([
            float::double(&p.x, &floor, dec)?,
            float::double(&p.y, &floor, dec)?,
            float::double(&p.z, &floor, dec)?,
        ])
    };
    let (pa, pb) = (
        image(&scratch.geom().position(a), dec)?,
        image(&scratch.geom().position(b), dec)?,
    );
    let (center, u, w) = match curve_form(scratch.geom(), carrier, dec)? {
        CurveForm::Line { origin, dir } => {
            let o = image(&origin, dec)?;
            let v = float::doubles(&dir, &floor, dec)?;
            let along = |p: [f64; 3]| {
                float::dot(std::array::from_fn(|k| p[k] - o[k]), v) / float::dot(v, v)
            };
            let t = (along(pa) + along(pb)) / 2.0;
            return Ok(vec![
                origin.plus(&dir.scaled(&Real::from(float::rational(t)?))),
            ]);
        }
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } => {
            let side = normal.cross(&ref_dir);
            let u = ref_dir.scaled(&float::root(&(&radius_sq / &ref_dir.norm_squared())));
            let w = side.scaled(&float::root(&(&radius_sq / &side.norm_squared())));
            (center, u, w)
        }
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            let side = normal.cross(&major_dir);
            let u = major_dir.scaled(&float::root(&(&major_sq / &major_dir.norm_squared())));
            let w = side.scaled(&float::root(&(&minor_sq / &side.norm_squared())));
            (center, u, w)
        }
        _ => return Ok(Vec::new()),
    };
    // `p = c + cos t u + sin t w` with `u` and `w` orthogonal.
    let (c, u64, w64) = (
        image(&center, dec)?,
        float::doubles(&u, &floor, dec)?,
        float::doubles(&w, &floor, dec)?,
    );
    let angle = |p: [f64; 3]| {
        let r: [f64; 3] = std::array::from_fn(|k| p[k] - c[k]);
        (float::dot(r, w64) / float::dot(w64, w64)).atan2(float::dot(r, u64) / float::dot(u64, u64))
    };
    let start = angle(pa);
    let mut end = angle(pb);
    while end <= start {
        end += std::f64::consts::TAU;
    }
    let mid = (start + end) / 2.0;
    Ok(vec![
        float::conic_point(&center, [&u, &w], mid)?,
        float::conic_point(&center, [&u, &w], mid + std::f64::consts::PI)?,
    ])
}

/// Whether `c` lies strictly between `a` and `b` along `carrier`: on a
/// closed carrier going round from `a` in its natural direction. An order
/// the exact predicates cannot decide counts as not inside, and the arc
/// point decides instead.
fn strictly_inside(
    scratch: &Model,
    carrier: CurveId,
    closed: bool,
    [a, c, b]: [PointId; 3],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if closed {
        return match cyclic_order(scratch.geom(), carrier, [a, c, b], dec) {
            Outcome::Proven(order) => Ok(order == Cyclic::Forward),
            Outcome::ProvenCoincident(_) | Outcome::NotExact(_) => Ok(false),
            Outcome::Refused(r) => Err(r),
        };
    }
    for (p, q) in [(a, c), (c, b)] {
        match order_on_curve(scratch.geom(), carrier, p, q, dec) {
            Outcome::Proven(StrictOrder::Less) => {}
            Outcome::Proven(_) | Outcome::ProvenCoincident(_) | Outcome::NotExact(_) => {
                return Ok(false);
            }
            Outcome::Refused(r) => return Err(r),
        }
    }
    Ok(true)
}

/// Merges runs of neighbouring kept pieces with equal keys, so that a key
/// holds along the whole of each run; a piece without a key is dropped. On
/// a closed carrier a run may wrap round; a carrier kept everywhere with one
/// key is one whole closed run.
pub(crate) fn runs<K: Clone + PartialEq>(
    spans: &[Ends],
    keys: Vec<Option<K>>,
    closed: bool,
) -> Vec<(Ends, K)> {
    let n = spans.len();
    let mut out = Vec::new();
    if n == 0 {
        return out;
    }
    if closed && keys[0].is_some() && keys.iter().all(|k| *k == keys[0]) {
        out.push(((None, None), keys[0].clone().expect("kept")));
        return out;
    }
    // Start a closed carrier where one run ends and another begins.
    let first = if closed {
        (0..n)
            .find(|&i| keys[i] != keys[(i + n - 1) % n])
            .unwrap_or(0)
    } else {
        0
    };
    let mut run: Option<(usize, usize)> = None;
    for step in 0..n {
        let i = (first + step) % n;
        if let Some((s, e)) = run
            && keys[i] != keys[s]
        {
            out.push(((spans[s].0, spans[e].1), keys[s].clone().expect("kept")));
            run = None;
        }
        if keys[i].is_some() {
            run = Some(run.map_or((i, i), |(s, _)| (s, i)));
        }
    }
    if let Some((s, e)) = run {
        out.push(((spans[s].0, spans[e].1), keys[s].clone().expect("kept")));
    }
    out
}

/// What a kept piece lies on: its own face, and both faces of every
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
) -> Result<Along, Refusal> {
    let mut on = vec![face];
    let mut along = None;
    if !in_plane.is_empty() {
        let probe = scratch
            .geom_mut()
            .add_point(PointDef::Explicit(interior.clone()));
        for &edge in in_plane {
            if !on_edge_span(scratch, edge, probe, dec)? {
                continue;
            }
            along.get_or_insert(edge);
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
    Ok(Along {
        faces: on,
        edge: along,
    })
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
    // An edge without vertices is its whole closed carrier.
    let ends = match (stored.start, stored.end) {
        (Some(start), Some(end)) => {
            let (start, end) = (scratch.vertex(start).point, scratch.vertex(end).point);
            Some(match stored.sense {
                Sense::Same => [start, end],
                Sense::Reversed => [end, start],
            })
        }
        _ => None,
    };
    on_span(scratch, stored.curve, ends, point, dec)
}

/// Whether `point`, which is neither end, lies on the arc of `curve` from
/// `ends[0]` to `ends[1]` along the curve's natural direction: on the
/// carrier, and on an open carrier between the ends, on a closed one on the
/// arc between them. No ends, or two equal ends, mean the whole carrier.
pub(crate) fn on_span(
    scratch: &Model,
    curve: CurveId,
    ends: Option<[PointId; 2]>,
    point: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if !lies_on(scratch, curve, point, dec)? {
        return Ok(false);
    }
    let Some([start, end]) = ends else {
        return Ok(true);
    };
    let geom = scratch.geom();
    if curve_is_closed(geom, curve, dec)? {
        if same_point(scratch, start, end, dec)? {
            return Ok(true);
        }
        return match cyclic_order(geom, curve, [start, point, end], dec) {
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
/// answer carries a [`CONTOUR_FIT_MARK`] mark. On an extruded Bezier profile
/// each arc is an exact ruling. On a complete rational patch, or a lofted
/// face, each arc is fitted the same way, and its ends are where the fitted
/// carrier starts and ends, within its deviation of the face boundary; a
/// generator crossing a chart seam of a lofted face is one arc per chart,
/// and the two arcs at a seam share their end point there.
/// Isolated generator points, such as a cone apex seen from inside the
/// cone's angle, draw nothing and are left out.
///
/// Arcs are not split where they cross in projection and carry no
/// visibility; hiding is a separate step. A horn or spindle torus sheet
/// refuses the whole query with `app.contour-surface`; an oblique torus face
/// bounded by anything other than its parallels and meridians, a rational
/// patch face with further trims, and spline bands from fillets and drafts
/// refuse with `app.contour-trim`.
pub fn contour_generators(
    model: &Model,
    body: BodyId,
    direction: &RVector3,
    dec: &mut Decider,
) -> Result<Answer<Contours>, QueryRefusal> {
    drawing_query(model, &[body], || {
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
                    jobs.extend(sections.into_iter().map(|s| (face, Job::Section(s))));
                }
                Generators::ObliqueTorus => jobs.push((face, Job::ObliqueTorus)),
                Generators::Spline => jobs.push((face, Job::Spline)),
            }
        }
        let mut scratch = model.clone();
        let mut marks = model.body(body).marks.clone();
        let mut arcs = Vec::new();
        for (face, job) in jobs {
            match job {
                Job::Section(section) => face_generators(
                    &mut scratch,
                    &faces,
                    face,
                    section,
                    &mut arcs,
                    &mut marks,
                    dec,
                )?,
                Job::Spline => {
                    spline::generators(&mut scratch, &faces, face, &d, &mut arcs, dec)?;
                }
                Job::ObliqueTorus => {
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
