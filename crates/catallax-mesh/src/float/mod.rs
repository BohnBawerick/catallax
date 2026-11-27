//! Binary64 B-rep data for the unproved display path. No mesh or certificate.

use std::collections::BTreeMap;

use catallax_numbers::Refusal;
use catallax_shapes::{Surface, SurfaceId};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

mod edges;
mod evaluate;
mod numerical;
mod spline;

use evaluate::{Projector, align, anchor, describe, finite, point, positive, sub};
pub use numerical::describe_numerical_body;

type Xyz = [f64; 3];
type Uv = [f64; 2];

/// Version of the coordinate conventions and records below.
pub const DISPLAY_DESCRIPTION_VERSION: u32 = 1;

/// One body's shared boundaries and face charts, in model traversal order.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayDescription<F = FaceId, E = EdgeId> {
    /// [`DISPLAY_DESCRIPTION_VERSION`].
    pub version: u32,
    /// World-space origin in millimetres. All other positions are relative to it.
    pub origin: [f64; 3],
    /// Each topological edge appears once.
    pub edges: Vec<DisplayEdge<E>>,
    /// Revision-local handles, or stable definition keys on numerical bodies.
    pub faces: Vec<DisplayFace<F>>,
}

/// One polyline in topological edge direction. A closed edge repeats its endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayEdge<E = EdgeId> {
    /// Revision-local picking handle or stable definition key.
    pub edge: E,
    /// Body-relative positions in millimetres.
    pub points: Vec<[f64; 3]>,
}

/// One oriented face and its chart.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayFace<F = FaceId> {
    /// Revision-local picking handle or stable definition key.
    pub face: F,
    /// Reverse the normal from `du × dv` when true.
    pub reversed: bool,
    /// Parameterization of the retained carrier.
    pub surface: DisplaySurface,
    /// Inclusive UV bounds of this face's chart, as `[lower, upper]`.
    pub bounds: [[f64; 2]; 2],
    /// UV period per axis, if any. Trim UVs may be outside the principal period.
    pub periods: [Option<f64>; 2],
    /// Boundary loops in their topological orientation; no duplicated closing use.
    pub loops: Vec<DisplayLoop>,
}

/// Consecutive oriented uses around a boundary loop.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayLoop {
    /// Directed uses; empty only for a single-vertex boundary.
    pub uses: Vec<DisplayEdgeUse>,
    /// Singular point boundary (for example a cone apex), in body coordinates.
    pub vertex: Option<([f64; 3], [f64; 2])>,
}

/// One occurrence of a shared edge, in loop traversal direction.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayEdgeUse {
    /// Index into [`DisplayDescription::edges`].
    pub edge: usize,
    /// Samples in traversal order, including both ends.
    pub samples: Vec<DisplayTrimSample>,
}

/// A UV occurrence of an edge sample. At a seam or pole the same point may have
/// more than one UV occurrence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayTrimSample {
    /// Index into the referenced edge's [`DisplayEdge::points`].
    pub point: usize,
    /// Coordinates in the incident face's chart.
    pub uv: [f64; 2],
}

/// Right-handed unit frame, with `y = z × x`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayFrame {
    /// Body-relative parameter origin.
    pub origin: [f64; 3],
    /// Reference direction.
    pub x: [f64; 3],
    /// Transverse direction.
    pub y: [f64; 3],
    /// Axis or plane normal.
    pub z: [f64; 3],
}

/// Analytic angles use turns. Analytic linear parameters use millimetres.
/// Let `radial(u) = x cos(2πu) + y sin(2πu)`.
#[derive(Clone, Debug, PartialEq)]
pub enum DisplaySurface {
    /// `origin + u*x + v*y`.
    Plane(DisplayFrame),
    /// `origin + radius*radial(u) + v*z`.
    Cylinder {
        /// Cylinder axis frame.
        frame: DisplayFrame,
        /// Radius in millimetres.
        radius: f64,
    },
    /// `apex + v*slope*radial(u) + v*z`.
    Cone {
        /// Frame with origin at the apex.
        frame: DisplayFrame,
        /// Radius gained per unit axial distance.
        slope: f64,
    },
    /// `center + radius*(cos(2πv)*radial(u) + sin(2πv)*z)`.
    Sphere {
        /// Frame with origin at the centre.
        frame: DisplayFrame,
        /// Radius in millimetres.
        radius: f64,
    },
    /// `center + (major + minor*cos(2πv))*radial(u) + minor*sin(2πv)*z`.
    Torus {
        /// Frame with origin at the centre.
        frame: DisplayFrame,
        /// Radius of the core circle.
        major: f64,
        /// Radius of the tube.
        minor: f64,
    },
    /// `origin + radii[0]*cos(2πu)*x + radii[1]*sin(2πu)*y + v*z`.
    EllipticCylinder {
        /// Cylinder axis frame.
        frame: DisplayFrame,
        /// Semi-axis lengths in frame X and Y order.
        radii: [f64; 2],
    },
    /// Tensor-product net in the retained knot coordinates, including fitted faces.
    BSpline(DisplaySpline),
    /// Bilinearly interpolated samples for a carrier without a direct float form.
    Grid(DisplayGrid),
}

/// Nonhomogeneous control points and optional rational weights, U varying fastest.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplaySpline {
    /// Polynomial degrees in U, V order.
    pub degree: [usize; 2],
    /// Expanded knot vectors, including multiplicities.
    pub knots: [Vec<f64>; 2],
    /// Control-point counts in U, V order.
    pub size: [usize; 2],
    /// Body-relative control points.
    pub points: Vec<[f64; 3]>,
    /// One weight per control, or `None` for a polynomial net.
    pub weights: Option<Vec<f64>>,
}

/// Uniform samples over [`DisplayFace::bounds`], including both ends, U fastest.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayGrid {
    /// Sample counts in U, V order.
    pub size: [usize; 2],
    /// Body-relative sample positions.
    pub points: Vec<[f64; 3]>,
}

/// Sampling controls for the display description, not proof tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayOptions {
    /// Visual target for the sampled edge chord distance in millimetres.
    pub chord_tolerance: f64,
    /// Maximum total number of points stored across the shared edge polylines.
    pub max_edge_points: usize,
}

impl Default for DisplayOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: 0.01,
            max_edge_points: 1_000_000,
        }
    }
}

fn validate_options(options: &DisplayOptions) -> Result<(), Refusal> {
    if positive(options.chord_tolerance).is_none() || options.max_edge_points < 2 {
        return Err(Refusal::new(
            "mesh.display-description-options",
            "a positive finite chord target and at least two edge points are required",
        ));
    }
    Ok(())
}

/// Read one retained B-rep into a versioned, unproved display description.
///
/// Requires a valid retained body. Numeric definitions are converted once at
/// their read boundary; sampling and parameter recovery use binary64 only.
/// This call does not validate the model, construct a mesh, or request proof.
pub fn describe_body(
    model: &Model,
    body: BodyId,
    options: &DisplayOptions,
) -> Result<DisplayDescription, Refusal> {
    validate_options(options)?;
    let geom = model.geom();
    let faces = model.faces_of_body(body);
    let edges = model.edges_of_body(body);
    let vertex_ids = model.vertices_of_body(body);
    let vertex = |v| point(&geom.position(model.vertex(v).point));
    let origin = if let Some(&v) = vertex_ids.first() {
        vertex(v).ok_or_else(|| failed("body origin exceeds binary64"))?
    } else if let Some(&f) = faces.first() {
        anchor(geom.surface(model.face(f).surface))
            .ok_or_else(|| failed("face carrier has no finite display origin"))?
    } else {
        [0.0; 3]
    };
    let mut vertices = BTreeMap::new();
    for v in vertex_ids {
        let p = vertex(v).ok_or_else(|| failed(format!("vertex {v:?} exceeds binary64")))?;
        // Normalize zero so topological vertices have byte-identical copies.
        let [x, y, z] = sub(p, origin).map(finite);
        let overflow = || failed(format!("vertex {v:?} exceeds body-relative binary64"));
        vertices.insert(
            v,
            [
                x.ok_or_else(overflow)?,
                y.ok_or_else(overflow)?,
                z.ok_or_else(overflow)?,
            ],
        );
    }

    let mut surfaces = BTreeMap::new();
    let mut incident = BTreeMap::<EdgeId, Vec<SurfaceId>>::new();
    for &id in &faces {
        let face = model.face(id);
        read_surface(model, &mut surfaces, face.surface, origin)?;
        for &boundary in &face.loops {
            if let Loop::Coedges(uses) = model.loop_(boundary) {
                for &coedge in uses {
                    let entries = incident.entry(model.coedge(coedge).edge).or_default();
                    if !entries.contains(&face.surface) {
                        entries.push(face.surface);
                    }
                }
            }
        }
    }
    for &edge in &edges {
        if let catallax_shapes::Curve::Intersection(c) = geom.curve(model.edge(edge).curve) {
            for surface in c.surfaces {
                read_surface(model, &mut surfaces, surface, origin)?;
            }
        }
    }
    let reader = edges::Reader {
        model,
        origin,
        options,
        vertices: &vertices,
        surfaces: &surfaces,
    };
    let mut result = DisplayDescription {
        version: DISPLAY_DESCRIPTION_VERSION,
        origin,
        edges: Vec::with_capacity(edges.len()),
        faces: Vec::with_capacity(faces.len()),
    };
    let mut edge_index = BTreeMap::new();
    let mut remaining = options.max_edge_points;
    for edge in edges {
        let supports = incident.get(&edge).map_or(&[][..], Vec::as_slice);
        let points = reader.edge(edge, supports, remaining).ok_or_else(|| {
            failed(format!(
                "edge {edge:?} ({:?}) could not be sampled",
                geom.curve(model.edge(edge).curve).kind()
            ))
        })?;
        remaining -= points.len();
        edge_index.insert(edge, result.edges.len());
        result.edges.push(DisplayEdge { edge, points });
    }

    let mut remaining_uses = BTreeMap::<SurfaceId, usize>::new();
    for &f in &faces {
        *remaining_uses.entry(model.face(f).surface).or_default() += 1;
    }
    for face in faces {
        let f = model.face(face);
        let source = geom.surface(f.surface);
        let count = remaining_uses.get_mut(&f.surface).unwrap();
        *count -= 1;
        let surface = if *count == 0 {
            surfaces.remove(&f.surface).unwrap()
        } else {
            surfaces[&f.surface].clone()
        };
        let mut chart = Projector::new(&surface, Some(source))
            .ok_or_else(|| failed(format!("face {face:?} has no finite display chart")))?;
        let mut loops = Vec::with_capacity(f.loops.len());
        let mut trim_bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
        for &boundary in &f.loops {
            let mut previous = None;
            let mut project = |p| -> Result<Uv, Refusal> {
                let mut uv = chart.project(p, previous);
                for i in 0..2 {
                    if let (Some(near), Some(period)) = (previous, chart.periods[i]) {
                        uv[i] = align(uv[i], near[i], period);
                    }
                    uv[i] = finite(uv[i])
                        .ok_or_else(|| failed(format!("face {face:?} has nonfinite trim UV")))?;
                    trim_bounds[i][0] = trim_bounds[i][0].min(uv[i]);
                    trim_bounds[i][1] = trim_bounds[i][1].max(uv[i]);
                }
                previous = Some(uv);
                Ok(uv)
            };
            let loop_ = match model.loop_(boundary) {
                Loop::Vertex(v) => {
                    let p = vertices[v];
                    DisplayLoop {
                        uses: Vec::new(),
                        vertex: Some((p, project(p)?)),
                    }
                }
                Loop::Coedges(uses) => {
                    let mut samples = Vec::with_capacity(uses.len());
                    for &coedge in uses {
                        let coedge = model.coedge(coedge);
                        let edge = edge_index[&coedge.edge];
                        let points = &result.edges[edge].points;
                        let samples_for_use = (0..points.len())
                            .map(|i| {
                                let point = if coedge.sense == Sense::Same {
                                    i
                                } else {
                                    points.len() - 1 - i
                                };
                                Ok(DisplayTrimSample {
                                    point,
                                    uv: project(points[point])?,
                                })
                            })
                            .collect::<Result<Vec<_>, Refusal>>()?;
                        samples.push(DisplayEdgeUse {
                            edge,
                            samples: samples_for_use,
                        });
                    }
                    DisplayLoop {
                        uses: samples,
                        vertex: None,
                    }
                }
            };
            loops.push(loop_);
        }
        let bounds = std::array::from_fn(|i| chart.domain[i].unwrap_or(trim_bounds[i]));
        if bounds
            .iter()
            .any(|[lo, hi]| !lo.is_finite() || !hi.is_finite() || lo >= hi)
        {
            return Err(failed(format!(
                "face {face:?} has an empty or unbounded display chart"
            )));
        }
        let periods = chart.periods;
        let reversed = (f.sense == Sense::Reversed) ^ carrier_reversed(source);
        result.faces.push(DisplayFace {
            face,
            reversed,
            surface,
            bounds,
            periods,
            loops,
        });
    }
    Ok(result)
}

fn read_surface(
    model: &Model,
    surfaces: &mut BTreeMap<SurfaceId, DisplaySurface>,
    id: SurfaceId,
    origin: Xyz,
) -> Result<(), Refusal> {
    if let std::collections::btree_map::Entry::Vacant(e) = surfaces.entry(id) {
        let s = model.geom().surface(id);
        e.insert(describe(s, origin).ok_or_else(|| {
            failed(format!(
                "carrier {id:?} ({:?}) has no finite float form",
                s.kind()
            ))
        })?);
    }
    Ok(())
}

fn carrier_reversed(surface: &Surface) -> bool {
    matches!(surface, Surface::BSpline(catallax_shapes::SplineSupport::LinearFillet(s)) if s.orientation() < 0)
}

fn failed(detail: impl Into<String>) -> Refusal {
    Refusal::new("mesh.display-description", detail)
}
