//! Binary64 B-rep data for the unproved display path. No mesh or certificate.

use catallax_topology::{EdgeId, FaceId};

/// Version of the coordinate conventions and records below.
pub const DISPLAY_DESCRIPTION_VERSION: u32 = 1;

/// One body's shared boundaries and face charts, in model traversal order.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayDescription {
    /// [`DISPLAY_DESCRIPTION_VERSION`].
    pub version: u32,
    /// World-space origin in millimetres. All other positions are relative to it.
    pub origin: [f64; 3],
    /// Each topological edge appears once.
    pub edges: Vec<DisplayEdge>,
    /// Face handles resolve in the retained model; the app supplies picking names.
    pub faces: Vec<DisplayFace>,
}

/// One polyline in topological edge direction. A closed edge repeats its endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayEdge {
    /// Revision-local picking handle.
    pub edge: EdgeId,
    /// Body-relative positions in millimetres.
    pub points: Vec<[f64; 3]>,
}

/// One oriented face and its chart.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayFace {
    /// Revision-local picking handle.
    pub face: FaceId,
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
