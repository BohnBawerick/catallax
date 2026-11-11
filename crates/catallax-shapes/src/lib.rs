//! Piece 2 of catallax: shapes.
//!
//! Exact geometric definitions (points, poses, analytic surfaces and curves,
//! intersection curves and vertices by definition), the write-once geometry
//! store with its fact book, and decisions about single objects.
//!
//! The contract is `docs/pieces/02-shapes.md`.

mod bezier_extrusion;
mod bspline;
mod curve;
/// Certified minimum-distance bounds on trimmed analytic carriers.
pub mod distance;
mod intent;
mod isometry;
mod linear_fillet;
mod math;
mod oblique_draft;
mod offset;
mod offset_chart;
mod offset_sweep;
mod parallel;
mod point;
mod pose;
mod rational_atlas;
mod rational_patch;
mod rational_sheet;
mod rolling_envelope;
mod rolling_mitre;
mod rolling_pipe;
mod rolling_saddle;
mod rolling_sphere;
mod rolling_torus;
mod spline_surface;
/// Exact exchange spline data and finite-domain knot operations.
pub mod spline_tools;
mod store;
mod surface;
mod vector;

#[cfg(test)]
mod tests;

pub use bezier_extrusion::BezierExtrusion;
pub use bspline::{
    BSplineCurve, BSplineCurveJet, BSplineJet, BSplineSurface, SplineBoundary, SplineCurveParameter,
};
pub use curve::{Circle, Curve, CurveKind, Ellipse, IntersectionCurve, Line};
pub use isometry::Isometry;
pub use linear_fillet::LinearFillet;
pub use oblique_draft::ObliqueDraft;
pub use offset::{EllipticOffsetDefinition, OffsetSide};
pub use offset_chart::EllipticOffsetChart;
pub use offset_sweep::{
    OffsetSweepProof, certify_convex_cubic_span, certify_elliptic_offset_sweep,
};
pub use point::{Meet, MeetOf, PointDef, Position};
pub use pose::{Pose, Rotation};
pub use rational_atlas::{AtlasChart, AtlasSection, RationalAtlas};
pub use rational_patch::{
    OffsetPatch, PairArc, PairCertificate, PairPoint, PatchPairBranch, PatchPairIntersection,
    PatchPairSearch, PatchSectionBranch, RationalGridEmbedding, RationalPatch, RationalPatchGrid,
    RationalPatchSection, SectionBranchApproximation, SectionCandidate, SectionHits, SectionSearch,
    SectionSurface,
};
pub use rolling_envelope::{
    BoundaryRow, ContactMoments, EnvelopeGrid, JunctionSection, ParameterTable,
    ROLLING_ENVELOPE_TARGET_NM, RollingEnvelope, RollingEnvelopeDefinition,
};
pub use rolling_mitre::RollingMitre;
pub use rolling_pipe::{PipeLaw, RollingPipe};
pub use rolling_saddle::{
    BallSide, ROLLING_SADDLE_TARGET_NM, RollingSaddle, RollingSaddleDefinition, SaddleLaw,
    SaddleRing,
};
pub use rolling_sphere::{ROLLING_SPHERE_TARGET_NM, RollingSphere, RollingSphereDefinition};
pub use rolling_torus::{InscribedTorus, RollingTorus};
pub use spline_surface::{SplinePointParameters, SplineSupport, SplineSurface};
pub use store::{Apart, Arena, CurveId, GeomKey, GeomStore, Handle, PointId, Side, SurfaceId};
pub use surface::{
    Cone, Cylinder, EllipticCylinder, ParamKind, Plane, Sphere, Surface, SurfaceKind, Torus,
    TorusSheet,
};
pub use vector::{Point3, RPoint3, RVector3, Vector3};
