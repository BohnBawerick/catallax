//! Piece 2 of catallax: shapes.
//!
//! Exact geometric definitions (points, poses, analytic surfaces and curves,
//! intersection curves and vertices by definition), the write-once geometry
//! store with its fact book, and decisions about single objects.
//!
//! The contract is `docs/pieces/02-shapes.md`.

mod bspline;
mod curve;
mod isometry;
mod math;
mod offset;
mod offset_chart;
mod point;
mod pose;
mod store;
mod surface;
mod vector;

#[cfg(test)]
mod tests;

pub use bspline::{BSplineJet, BSplineSurface};
pub use curve::{Circle, Curve, CurveKind, Ellipse, IntersectionCurve, Line};
pub use isometry::Isometry;
pub use offset::{EllipticOffsetDefinition, OffsetSide};
pub use offset_chart::EllipticOffsetChart;
pub use point::{Meet, MeetOf, PointDef, Position};
pub use pose::{Pose, Rotation};
pub use store::{Apart, Arena, CurveId, GeomKey, GeomStore, Handle, PointId, Side, SurfaceId};
pub use surface::{
    Cone, Cylinder, EllipticCylinder, ParamKind, Plane, Sphere, Surface, SurfaceKind, Torus,
};
pub use vector::{Point3, RPoint3, RVector3, Vector3};
