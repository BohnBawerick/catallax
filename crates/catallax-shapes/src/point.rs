//! Points stored by definition.

use catallax_numbers::IsolatedRoot;

use crate::{CurveId, Point3, RPoint3, SurfaceId};

/// The definition of a point.
#[derive(Clone, Debug)]
pub enum PointDef {
    /// An exact rational point.
    Explicit(RPoint3),
    /// Where some objects meet.
    Meet(Meet),
}

/// Where these objects meet, and which point of the meeting.
#[derive(Clone, Debug)]
pub struct Meet {
    /// The objects that meet.
    pub of: MeetOf,
    /// Which point of the meeting.
    pub position: Position,
}

/// The objects a point is defined to lie on.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum MeetOf {
    /// Two or more surfaces.
    Surfaces(Vec<SurfaceId>),
    /// A curve and a surface.
    CurveSurface(CurveId, SurfaceId),
    /// Two curves.
    Curves(CurveId, CurveId),
}

/// Which point of a meeting.
#[derive(Clone, Debug)]
pub enum Position {
    /// Exact closed-form coordinates.
    Closed(Point3),
    /// A certified root of a system in `x, y, z`.
    Isolated(IsolatedRoot),
}
