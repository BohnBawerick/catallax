//! Exact selected-edge definitions, independent of the display mesh.
//!
//! A finite circle or ellipse span is selected by its ordered endpoints, the edge's sense
//! relative to the carrier, and an exact interior point. No rounded angle or
//! polyline chooses the branch. Handles and meeting definitions belong to the
//! queried model; persist origin names and resolve them again after an edit.
//!
//! Query a result picking handle with [`ApplicationBody::edge_curve`] or an
//! incremental checkpoint with [`IncrementalApplicationBody::edge_curve`].
//! Inspect [`CurveSpan`] first: only lines, full circles and ellipses, and proper
//! circular and elliptical arcs and finite parabola and hyperbola arcs are
//! supported output. The other carriers retain their exact form and
//! trim with `Unsupported`. An underlying intersection-query refusal propagates.
//! Finite trim endpoints follow edge traversal; apply a face's coedge sense
//! separately when walking a wire. `crosses_seam` excludes seam endpoints and
//! does not replace the exact endpoints and interior point as the span definition.
//!
//! These queries borrow a single result and use its finite topology, the caller's
//! decider and any active arithmetic allowance. They never rebuild a recipe,
//! modify topology, read display samples or decide drawing visibility. Returned
//! geometry describes the retained model; any construction marks remain on that
//! model and are not erased by querying its definitions.

use catallax_intersect::{CurveForm, Cyclic, arc_point, curve_form, cyclic_order, order_on_curve};
use catallax_numbers::{Decider, Outcome, Real, Refusal, StrictOrder};
use catallax_recipe::{EdgeName, Origin};
use catallax_shapes::{Curve, CurveId, Point3, PointDef, PointId};
use catallax_topology::{BodyId, EdgeId, Loop, Model, Sense, VertexId};

use crate::{ApplicationBody, IncrementalApplicationBody};

/// An exact endpoint, retaining both its construction and its coordinates.
#[derive(Clone, Debug)]
pub struct CurveEndpoint {
    /// Model-local picking handle.
    pub vertex: VertexId,
    /// Model-local point handle.
    pub point: PointId,
    /// Original point definition, including meeting or isolated-root identity.
    pub definition: PointDef,
    /// Exact coordinates, never obtained from display positions.
    pub position: Point3,
}

/// Topological bounds in edge traversal order, before applying any coedge sense.
#[derive(Clone, Debug)]
pub enum CurveTrim {
    /// A ring edge, with no topological endpoints.
    Closed,
    /// Both endpoints, in edge traversal order even when the edge is reversed.
    Finite {
        /// First endpoint.
        start: Box<CurveEndpoint>,
        /// Last endpoint.
        end: Box<CurveEndpoint>,
    },
}

/// Supported finite span, or an explicit unsupported carrier.
#[derive(Clone, Debug)]
pub enum CurveSpan {
    /// The finite line segment between the trim endpoints.
    Line,
    /// A complete circle, including a closed edge with coincident endpoints.
    FullCircle,
    /// A proper circular arc, which may cross the carrier's parameter seam.
    CircularArc {
        /// Exact point strictly inside the selected branch; not an arc midpoint.
        interior: Point3,
        /// The seam lies strictly inside the arc. A seam endpoint is excluded.
        crosses_seam: bool,
    },
    /// A complete ellipse, including a closed edge with coincident endpoints.
    FullEllipse,
    /// A proper elliptical arc, which may cross the major-axis seam.
    EllipticArc {
        /// Exact point strictly inside the selected branch; not an arc midpoint.
        interior: Point3,
        /// The seam lies strictly inside the arc. A seam endpoint is excluded.
        crosses_seam: bool,
    },
    /// A finite parabola arc or single-branch hyperbola arc between distinct
    /// endpoints. The carrier kind and branch are on [`EdgeCurve::form`].
    OpenConicArc {
        /// Exact point strictly inside the selected arc, not a polyline vertex.
        interior: Point3,
        /// Transverse coordinate of the arc's endpoint that comes first along
        /// the carrier's natural direction, in the carrier frame. The frame
        /// origin is the parabola vertex or hyperbola centre and the axis is
        /// the unit of `normal x axis_dir`.
        first: Real,
        /// Transverse coordinate of the endpoint that comes last along the
        /// carrier's natural direction. The coordinate is strictly monotone
        /// on the arc, so the arc is exactly the carrier points whose
        /// coordinate lies between `first` and `last`. The edge's start is
        /// `first` when the sense is `Same` and `last` when it is `Reversed`.
        last: Real,
    },
    /// A spline or general intersection branch.
    /// The exact carrier and topological trim remain available on [`EdgeCurve`].
    Unsupported,
}

/// One selected edge's exact geometry and retained source identity.
#[derive(Clone, Debug)]
pub struct EdgeCurve {
    /// Model-local edge handle.
    pub edge: EdgeId,
    /// Model-local curve handle, including its original intersection component.
    pub curve: CurveId,
    /// Original stored curve definition, not a display fit or substituted carrier.
    pub definition: Curve,
    /// Exact carrier from the existing intersection query.
    pub form: CurveForm,
    /// Edge traversal relative to the carrier's natural direction.
    pub sense: Sense,
    /// Topological endpoints, independent of the carrier's parameter seam.
    pub trim: CurveTrim,
    /// Supported line/circle/ellipse span or explicit unsupported status.
    pub span: CurveSpan,
    /// Stable origin-pair name when the result admits one; not necessarily unique.
    pub name: Option<EdgeName>,
    /// All incident face origins, sorted and deduplicated, including merged sets.
    pub origins: Vec<Origin>,
    /// Original construction pairs retained through splitting and merging.
    pub constructions: Vec<EdgeName>,
}

impl ApplicationBody {
    /// Queries one result edge without reading or changing its display mesh.
    /// Uses the caller's decider and any active work allowance. Uncertain span
    /// decisions refuse without discarding this body or weakening exactness.
    pub fn edge_curve(&self, edge: EdgeId, decider: &mut Decider) -> Result<EdgeCurve, Refusal> {
        edge_curve(&self.built.model, self.built.result, edge, decider)
    }
}

impl IncrementalApplicationBody {
    /// Queries an edge of this checkpoint's result, with the same contract as
    /// [`ApplicationBody::edge_curve`]. Handles from another edit are not stable.
    pub fn edge_curve(&self, edge: EdgeId, decider: &mut Decider) -> Result<EdgeCurve, Refusal> {
        let built = self.checkpoint.built();
        edge_curve(&built.model, built.result, edge, decider)
    }
}

/// Queries one edge of `body`, using exact topology and existing curve queries.
///
/// Work is confined to this body's finite topology and one curve's definition
/// and span. It performs no recipe rebuild, tessellation or curve fitting.
/// As with other model queries, the source topology must already be valid.
/// An edge outside the result returns `app.curve-edge`; malformed endpoint
/// bounds return `app.curve-trim`. Foreign-model handles are programming errors,
/// as with other topology queries. Existing intersection refusals propagate.
pub fn edge_curve(
    model: &Model,
    body: BodyId,
    edge: EdgeId,
    decider: &mut Decider,
) -> Result<EdgeCurve, Refusal> {
    decider.check_deadline()?;
    if !model.edges_of_body(body).contains(&edge) {
        return Err(Refusal::new(
            "app.curve-edge",
            "edge is not in the queried body",
        ));
    }
    let source = model.edge(edge);
    let endpoint = |vertex| {
        let point = model.vertex(vertex).point;
        CurveEndpoint {
            vertex,
            point,
            definition: model.geom().point(point).clone(),
            position: model.geom().position(point),
        }
    };
    let trim = match (source.start, source.end) {
        (None, None) => CurveTrim::Closed,
        (Some(start), Some(end)) => CurveTrim::Finite {
            start: Box::new(endpoint(start)),
            end: Box::new(endpoint(end)),
        },
        _ => return Err(bad_trim()),
    };
    let form = curve_form(model.geom(), source.curve, decider)?;
    let span = match (&form, &trim) {
        (CurveForm::Line { .. }, CurveTrim::Finite { .. }) => CurveSpan::Line,
        (CurveForm::Line { .. }, CurveTrim::Closed) => return Err(bad_trim()),
        (CurveForm::Circle { .. }, CurveTrim::Closed) => CurveSpan::FullCircle,
        (CurveForm::Ellipse { .. }, CurveTrim::Closed) => CurveSpan::FullEllipse,
        (
            CurveForm::Circle {
                center,
                radius_sq,
                ref_dir,
                ..
            },
            CurveTrim::Finite { start, end },
        ) => {
            let seam = center.plus(&ref_dir.scaled(&(radius_sq / &ref_dir.norm_squared()).sqrt()));
            match finite_arc(model, source, start, end, seam, decider)? {
                None => CurveSpan::FullCircle,
                Some((interior, crosses_seam)) => CurveSpan::CircularArc {
                    interior,
                    crosses_seam,
                },
            }
        }
        (
            CurveForm::Ellipse {
                center,
                major_dir,
                major_sq,
                ..
            },
            CurveTrim::Finite { start, end },
        ) => {
            let seam =
                center.plus(&major_dir.scaled(&(major_sq / &major_dir.norm_squared()).sqrt()));
            match finite_arc(model, source, start, end, seam, decider)? {
                None => CurveSpan::FullEllipse,
                Some((interior, crosses_seam)) => CurveSpan::EllipticArc {
                    interior,
                    crosses_seam,
                },
            }
        }
        (CurveForm::Parabola { .. }, CurveTrim::Closed)
        | (CurveForm::Hyperbola { .. }, CurveTrim::Closed) => return Err(bad_trim()),
        (
            CurveForm::Parabola {
                vertex: origin,
                normal,
                axis_dir,
                ..
            }
            | CurveForm::Hyperbola {
                center: origin,
                normal,
                axis_dir,
                ..
            },
            CurveTrim::Finite { start, end },
        ) => {
            if exact(model.geom().same_point(start.point, end.point, decider))?.is_none() {
                return Err(bad_trim());
            }
            let (from, to) = match source.sense {
                Sense::Same => (start.point, end.point),
                Sense::Reversed => (end.point, start.point),
            };
            // The existing open-curve order decides the carrier direction; a
            // reversed pair is a malformed edge, not a span to flip silently.
            match order_on_curve(model.geom(), source.curve, from, to, decider) {
                Outcome::Proven(StrictOrder::Less) => {}
                Outcome::Proven(_) | Outcome::ProvenCoincident(_) => return Err(bad_trim()),
                Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
                Outcome::Refused(r) => return Err(r),
            }
            let interior =
                arc_point(model.geom(), source.curve, Some(from), Some(to), decider)?.value;
            let v = normal.cross(axis_dir);
            let v = v.scaled(&(Real::one() / v.norm_squared().sqrt()));
            let coordinate = |point| model.geom().position(point).minus(origin).dot(&v);
            CurveSpan::OpenConicArc {
                interior,
                first: coordinate(from),
                last: coordinate(to),
            }
        }
        _ => CurveSpan::Unsupported,
    };
    let mut origins = Vec::new();
    for face in model.faces_of_body(body) {
        decider.check_deadline()?;
        if model.face(face).loops.iter().any(|lp| {
            matches!(model.loop_(*lp), Loop::Coedges(cs) if cs.iter().any(|c| model.coedge(*c).edge == edge))
        }) {
            origins.extend_from_slice(model.face_origins(face));
        }
    }
    origins.sort();
    origins.dedup();
    let result = EdgeCurve {
        edge,
        curve: source.curve,
        definition: model.geom().curve(source.curve).clone(),
        form,
        sense: source.sense,
        trim,
        span,
        name: model.edge_name(body, edge),
        origins,
        constructions: model.edge_constructions(edge).to_vec(),
    };
    decider.check_deadline()?;
    Ok(result)
}

/// The proper arc between two distinct endpoints of a circle or ellipse, as its
/// exact interior point and strict seam crossing; `None` when they coincide.
fn finite_arc(
    model: &Model,
    source: &catallax_topology::Edge,
    start: &CurveEndpoint,
    end: &CurveEndpoint,
    seam: Point3,
    decider: &mut Decider,
) -> Result<Option<(Point3, bool)>, Refusal> {
    if exact(model.geom().same_point(start.point, end.point, decider))?.is_none() {
        return Ok(None);
    }
    let (from, to) = match source.sense {
        Sense::Same => (start.point, end.point),
        Sense::Reversed => (end.point, start.point),
    };
    let interior = arc_point(model.geom(), source.curve, Some(from), Some(to), decider)?.value;
    // A private persistent-store clone supplies the seam point to the existing
    // cyclic query without mutating the source model.
    let mut geom = model.geom().clone();
    let seam = geom.add_point(PointDef::Explicit(seam));
    let crosses_seam = matches!(
        exact(cyclic_order(&geom, source.curve, [from, seam, to], decider))?,
        Some(Cyclic::Forward)
    );
    Ok(Some((interior, crosses_seam)))
}

fn exact<T>(outcome: Outcome<T>) -> Result<Option<T>, Refusal> {
    match outcome {
        Outcome::Proven(value) => Ok(Some(value)),
        Outcome::ProvenCoincident(_) => Ok(None),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

fn bad_trim() -> Refusal {
    Refusal::new(
        "app.curve-trim",
        "line and open conic need distinct finite bounds in carrier order; endpoint bounds must be paired",
    )
}

#[cfg(test)]
mod tests;
