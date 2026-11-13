//! Operations and their parameters. `docs/RECIPE.md` gives each one's definition.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::ids::{ElementId, StepId};
use crate::names::{EdgeRef, EdgeSelection, FaceRef, VertexName};
use crate::number::{Axis, Number, Plane, Point2, Point3, Pose, Vector3};

/// What a step's value is.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum ValueKind {
    /// A planar region, made by a sketch.
    Profile,
    /// A chain of curves in space, made by a path.
    Path,
    /// A solid.
    Body,
    /// An open scan: reference data with no inside or outside.
    Scan,
}

impl ValueKind {
    /// The kind's name, for messages.
    pub fn name(self) -> &'static str {
        match self {
            ValueKind::Profile => "profile",
            ValueKind::Path => "path",
            ValueKind::Body => "body",
            ValueKind::Scan => "scan",
        }
    }
}

/// The authored policy when a rolling-ball contact reaches a neighbour.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FilletOverflow {
    /// Continue across neighbours, consuming their trimmed regions when required.
    RollOver,
    /// Preserve unselected neighbours and close runs at explicit station planes.
    KeepNeighbours(Vec<FilletSetback>),
}

/// How an exposed side of a selected fillet junction closes.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FilletCorner {
    /// Close the active normal-cone cell against its exposed side plane.
    Cliff,
    /// Continue a radius-preserving contact sheet across the exposed side.
    Smooth,
}

/// A station measured along a selected input edge from one of its endpoints.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FilletSetback {
    /// Must resolve to exactly one selected input edge.
    pub edge: EdgeRef,
    /// The endpoint's canonical input vertex name.
    pub at: VertexName,
    /// Positive distance along the edge from that endpoint, in millimetres.
    pub setback: Number,
}

/// One operation with its parameters.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum Operation {
    /// The box from the local origin to `size`, placed by `pose`.
    Block {
        /// Positive side lengths in millimetres.
        size: [Number; 3],
        /// Rigid placement.
        pose: Pose,
    },
    /// A capped cylinder on local `z`, from zero to `height`.
    Cylinder {
        /// Positive radius.
        radius: Number,
        /// Positive height.
        height: Number,
        /// Rigid placement.
        pose: Pose,
    },
    /// A sphere centred at the local origin.
    Sphere {
        /// Positive radius.
        radius: Number,
        /// Rigid placement.
        pose: Pose,
    },
    /// A capped circular frustum on local `z`.
    Cone {
        /// Non-negative radius at zero.
        base_radius: Number,
        /// Non-negative radius at `height`; at least one radius is positive.
        top_radius: Number,
        /// Positive height.
        height: Number,
        /// Rigid placement.
        pose: Pose,
    },
    /// A ring torus about local `z`, centred at the local origin.
    Torus {
        /// Radius of the core circle.
        major: Number,
        /// Positive tube radius, strictly less than `major`.
        minor: Number,
        /// Rigid placement.
        pose: Pose,
    },
    /// A planar region bounded by closed loops.
    Sketch {
        /// Places the sketch: its plane is the local `z = 0` plane, its normal local `+z`.
        plane: Pose,
        /// The loops. The region is the set of points enclosed by an odd number of them.
        loops: Vec<Loop>,
    },
    /// A tangent-continuous chain of curves in space, for a sweep.
    Path {
        /// Where the first segment starts.
        start: Point3,
        /// The segments, each starting where the previous one ends.
        segments: Vec<PathSegment>,
    },
    /// A profile swept along a straight line.
    Extrude {
        /// The sketch step.
        profile: StepId,
        /// Signed length of the sweep along `direction`, not zero.
        distance: Number,
        /// Sweep direction, of any length; the sketch normal when absent.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        direction: Option<Vector3>,
    },
    /// Cuts a profile through a target in both directions, without end caps.
    ExtrudeThrough {
        /// The body to cut.
        target: StepId,
        /// The sketch step.
        profile: StepId,
        /// Tool direction; the sketch normal when absent.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        direction: Option<Vector3>,
        /// Which half of the sweep to cut; version 4 only.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        side: Option<ThroughSide>,
    },
    /// A body-valued prism capped by the first positive hit on selected trimmed faces.
    ExtrudeTo {
        /// Body supplying the limiting faces.
        target: StepId,
        /// Sketch region to sweep.
        profile: StepId,
        /// Sweep direction; the sketch normal when absent.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        direction: Option<Vector3>,
        /// The selected trimmed faces.
        limit: FaceRef,
    },
    /// A body-valued prism capped by the first positive material crossing.
    ExtrudeNext {
        /// Body supplying the limiting boundary.
        target: StepId,
        /// Sketch region to sweep.
        profile: StepId,
        /// Sweep direction; the sketch normal when absent.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        direction: Option<Vector3>,
    },
    /// A profile turned about an axis in its plane.
    Revolve {
        /// The sketch step.
        profile: StepId,
        /// The axis; positive angles turn counterclockwise about its direction.
        axis: Axis,
        /// Angle in degrees, in `(0, 360]`.
        angle: Number,
    },
    /// Holes drilled into a body.
    Hole {
        /// The body.
        target: StepId,
        /// Where each hole starts and which way it goes.
        at: Vec<HolePlacement>,
        /// Diameter of the main wall.
        diameter: Number,
        /// Depth of the main wall below the placement point.
        depth: HoleDepth,
        /// The bottom of a blind hole.
        tip: HoleTip,
        /// A counterbore or countersink at the entry.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        head: Option<HoleHead>,
    },
    /// A regularised set operation on two bodies.
    Boolean {
        /// Which operation.
        kind: BooleanKind,
        /// The first body.
        target: StepId,
        /// The second body.
        tool: StepId,
    },
    /// Edges rounded by a rolling ball.
    Fillet {
        /// The body.
        body: StepId,
        /// The edges.
        edges: EdgeSelection,
        /// Ball radius.
        radius: Number,
        /// Version-4 continuation policy; absent retains the recipe version's default.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        overflow: Option<FilletOverflow>,
        /// Version-4 junction closure; absence keeps the cliff rule.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        corner: Option<FilletCorner>,
    },
    /// Edges rounded by explicitly oriented radius laws.
    FilletVariable {
        /// The body.
        body: StepId,
        /// One law per selected open edge.
        edges: Vec<VariableFilletEdge>,
    },
    /// Edges bevelled by a ruled face.
    Chamfer {
        /// The body.
        body: StepId,
        /// The edges.
        edges: EdgeSelection,
        /// Setback on the face matching the first origin pattern of each reference.
        distance: Number,
        /// Setback on the other face; `distance` when absent.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        other_distance: Option<Number>,
    },
    /// A body hollowed to a wall of constant thickness.
    Shell {
        /// The body.
        body: StepId,
        /// Faces left open.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        open: Vec<FaceRef>,
        /// Wall thickness.
        thickness: Number,
    },
    /// A body grown or shrunk by a distance.
    Offset {
        /// The body.
        body: StepId,
        /// Positive grows, negative shrinks.
        distance: Number,
    },
    /// Faces of a body made into a new solid of constant thickness.
    Thicken {
        /// The body whose faces are thickened.
        body: StepId,
        /// The faces.
        faces: Vec<FaceRef>,
        /// Thickness.
        thickness: Number,
        /// Which side of the faces the solid lies on.
        side: ThickenSide,
    },
    /// Faces tilted away from a pull direction about a neutral plane.
    Draft {
        /// The body.
        body: StepId,
        /// The faces to tilt.
        faces: Vec<FaceRef>,
        /// The direction the part leaves its mould in.
        pull: Vector3,
        /// Where the drafted faces keep their place.
        neutral: Plane,
        /// Draft angle in degrees, in `(-90, 90)`.
        angle: Number,
        /// Saved interior points selecting retained components instead of the
        /// neutral-plane rule. An explicit empty list selects the empty result.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        keep_regions: Option<Vec<Point3>>,
    },
    /// A profile carried along a path.
    Sweep {
        /// The sketch step.
        profile: StepId,
        /// The path step.
        path: StepId,
        /// How the profile turns along the path.
        orientation: SweepOrientation,
    },
    /// A solid through a sequence of profiles.
    Loft {
        /// The profiles, at least two, in order.
        sections: Vec<LoftSection>,
        /// How the side faces pass between sections.
        shape: LoftShape,
    },
    /// Copies of a body, united.
    Pattern {
        /// The body.
        body: StepId,
        /// Where the copies go.
        layout: PatternLayout,
    },
    /// A body reflected in a plane.
    Mirror {
        /// The body.
        body: StepId,
        /// The mirror plane.
        plane: Plane,
        /// Unites the reflection with the original when true.
        keep_original: bool,
    },
    /// Beams along a repeated cell, inside a body.
    Lattice {
        /// The body the lattice fills.
        region: StepId,
        /// The cell.
        cell: LatticeCell,
        /// Edge length of the cubic cell.
        size: Number,
        /// Places the cell grid.
        frame: Pose,
        /// Beam radius.
        radius: Grade,
    },
    /// A gyroid inside a body.
    Gyroid {
        /// The body the gyroid fills.
        region: StepId,
        /// Length of one period.
        period: Number,
        /// The level of the gyroid function that bounds the solid.
        level: Grade,
        /// A sheet around the surface or one of the two networks it separates.
        kind: GyroidKind,
        /// Places the gyroid.
        frame: Pose,
    },
    /// A set operation with the crease rounded through the bodies' distance fields.
    SmoothBlend {
        /// Which operation.
        kind: BooleanKind,
        /// The first body.
        target: StepId,
        /// The second body.
        tool: StepId,
        /// Blend radius.
        radius: Number,
    },
    /// A STEP file entering the recipe as a fixed piece.
    ImportStep {
        /// The file.
        file: FileRef,
        /// Source solid entity key, bound to `file.sha256`. Omission imports all solids.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source_solid: Option<String>,
        /// Places the imported body.
        pose: Pose,
        /// Saved answers to import questions, such as whether a gap is intended.
        #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
        decisions: BTreeMap<String, String>,
    },
    /// An IGES file's manifold solids entering the recipe as a fixed piece,
    /// exactly as [`Operation::ImportStep`] admits a STEP file's.
    ImportIges {
        /// The file.
        file: FileRef,
        /// Source solid entity key, bound to `file.sha256`. Omission imports all solids.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source_solid: Option<String>,
        /// Places the imported body.
        pose: Pose,
        /// Saved answers to import questions, such as whether a gap is intended.
        #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
        decisions: BTreeMap<String, String>,
    },
    /// A closed triangle mesh entering the recipe as a fixed piece.
    ImportMesh {
        /// The file.
        file: FileRef,
        /// Places the imported body.
        pose: Pose,
        /// Millimetres per file unit.
        scale: Number,
    },
    /// An open scan entering the recipe as reference data, never a body.
    ScanImport {
        /// The file.
        file: FileRef,
        /// Places the scan.
        pose: Pose,
        /// Millimetres per file unit.
        scale: Number,
    },
}

impl Operation {
    /// The operation's name, as written in a recipe.
    pub fn name(&self) -> &'static str {
        match self {
            Operation::Block { .. } => "block",
            Operation::Cylinder { .. } => "cylinder",
            Operation::Sphere { .. } => "sphere",
            Operation::Cone { .. } => "cone",
            Operation::Torus { .. } => "torus",
            Operation::Sketch { .. } => "sketch",
            Operation::Path { .. } => "path",
            Operation::Extrude { .. } => "extrude",
            Operation::ExtrudeThrough { .. } => "extrude_through",
            Operation::ExtrudeTo { .. } => "extrude_to",
            Operation::ExtrudeNext { .. } => "extrude_next",
            Operation::Revolve { .. } => "revolve",
            Operation::Hole { .. } => "hole",
            Operation::Boolean { .. } => "boolean",
            Operation::Fillet { .. } => "fillet",
            Operation::FilletVariable { .. } => "fillet_variable",
            Operation::Chamfer { .. } => "chamfer",
            Operation::Shell { .. } => "shell",
            Operation::Offset { .. } => "offset",
            Operation::Thicken { .. } => "thicken",
            Operation::Draft { .. } => "draft",
            Operation::Sweep { .. } => "sweep",
            Operation::Loft { .. } => "loft",
            Operation::Pattern { .. } => "pattern",
            Operation::Mirror { .. } => "mirror",
            Operation::Lattice { .. } => "lattice",
            Operation::Gyroid { .. } => "gyroid",
            Operation::SmoothBlend { .. } => "smooth_blend",
            Operation::ImportStep { .. } => "import_step",
            Operation::ImportIges { .. } => "import_iges",
            Operation::ImportMesh { .. } => "import_mesh",
            Operation::ScanImport { .. } => "scan_import",
        }
    }

    /// What the operation makes.
    pub fn kind(&self) -> ValueKind {
        match self {
            Operation::Sketch { .. } => ValueKind::Profile,
            Operation::Path { .. } => ValueKind::Path,
            Operation::ScanImport { .. } => ValueKind::Scan,
            _ => ValueKind::Body,
        }
    }

    /// The steps the operation reads, each with the kind of value it needs, in
    /// parameter order.
    pub fn inputs(&self) -> Vec<(&StepId, ValueKind)> {
        use ValueKind::{Body, Path, Profile};
        match self {
            Operation::Sketch { .. }
            | Operation::Block { .. }
            | Operation::Cylinder { .. }
            | Operation::Sphere { .. }
            | Operation::Cone { .. }
            | Operation::Torus { .. }
            | Operation::Path { .. }
            | Operation::ImportStep { .. }
            | Operation::ImportIges { .. }
            | Operation::ImportMesh { .. }
            | Operation::ScanImport { .. } => vec![],
            Operation::Extrude { profile, .. } | Operation::Revolve { profile, .. } => {
                vec![(profile, Profile)]
            }
            Operation::ExtrudeThrough {
                target, profile, ..
            }
            | Operation::ExtrudeTo {
                target, profile, ..
            }
            | Operation::ExtrudeNext {
                target, profile, ..
            } => vec![(target, Body), (profile, Profile)],
            Operation::Hole { target, .. } => vec![(target, Body)],
            Operation::Boolean { target, tool, .. }
            | Operation::SmoothBlend { target, tool, .. } => {
                vec![(target, Body), (tool, Body)]
            }
            Operation::Fillet { body, .. }
            | Operation::FilletVariable { body, .. }
            | Operation::Chamfer { body, .. }
            | Operation::Shell { body, .. }
            | Operation::Offset { body, .. }
            | Operation::Thicken { body, .. }
            | Operation::Draft { body, .. }
            | Operation::Pattern { body, .. }
            | Operation::Mirror { body, .. } => vec![(body, Body)],
            Operation::Sweep { profile, path, .. } => vec![(profile, Profile), (path, Path)],
            Operation::Loft { sections, .. } => {
                sections.iter().map(|s| (&s.profile, Profile)).collect()
            }
            Operation::Lattice { region, .. } | Operation::Gyroid { region, .. } => {
                vec![(region, Body)]
            }
        }
    }

    /// The ids of the elements inside the operation, in parameter order: sketch
    /// segments, path segments or hole placements.
    pub fn elements(&self) -> Vec<&ElementId> {
        match self {
            Operation::Sketch { loops, .. } => loops
                .iter()
                .flat_map(|l| l.0.iter().map(Segment::id))
                .collect(),
            Operation::Path { segments, .. } => segments.iter().map(PathSegment::id).collect(),
            Operation::Hole { at, .. } => at.iter().map(|p| &p.id).collect(),
            _ => vec![],
        }
    }
}

/// A radius law on the normalized arc length of an oriented input edge.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum RadiusLaw {
    /// Linear interpolation between endpoint radii in millimetres.
    /// Both radii are non-negative and at least one is positive.
    Linear {
        /// Radius at the named starting vertex.
        start: Number,
        /// Radius at the opposite endpoint.
        end: Number,
    },
}

/// A selected open edge, its orientation and its radius law.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VariableFilletEdge {
    /// Must resolve to exactly one open input edge.
    pub edge: EdgeRef,
    /// Must resolve to one endpoint of that edge.
    pub from: VertexName,
    /// Radius along the edge, starting at `from`.
    pub radius: RadiusLaw,
}

/// The side of an unbounded profile tool, relative to its authored sketch plane.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ThroughSide {
    /// Both directions, with no physical end cap.
    Both,
    /// Nonnegative line parameters along the sweep direction.
    Forward,
    /// Nonpositive line parameters along the sweep direction.
    Backward,
}

/// A closed chain of sketch segments. Each segment starts where the previous one ends,
/// and the first starts where the last ends, so a loop is closed by construction.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Loop(pub Vec<Segment>);

/// One piece of a loop, in sketch coordinates.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Segment {
    /// A polynomial Bezier span from the previous endpoint to `to`.
    Bezier {
        /// The segment's id.
        id: ElementId,
        /// One interior control for a quadratic, two for a cubic.
        controls: Vec<Point2>,
        /// End point.
        to: Point2,
    },
    /// A straight segment to `to`.
    Line {
        /// The segment's id.
        id: ElementId,
        /// End point.
        to: Point2,
    },
    /// A circular arc about `centre` to `to`.
    Arc {
        /// The segment's id.
        id: ElementId,
        /// Centre; both ends lie exactly at the same distance from it.
        centre: Point2,
        /// End point.
        to: Point2,
        /// The way round, seen from the side the sketch normal points to.
        turn: Turn,
    },
    /// A full circle, alone in its loop.
    Circle {
        /// The segment's id.
        id: ElementId,
        /// Centre.
        centre: Point2,
        /// Radius.
        radius: Number,
    },
    /// A full ellipse, alone in its loop. Its second axis is `axis.perp() * ratio`.
    Ellipse {
        /// The segment's id.
        id: ElementId,
        /// Centre.
        centre: Point2,
        /// Nonzero first semi-axis vector, in sketch coordinates.
        axis: Point2,
        /// Positive ratio of the second semi-axis to the first.
        ratio: Number,
    },
    /// An elliptical arc. Both endpoints lie exactly on its ellipse.
    EllipticArc {
        /// The segment's id.
        id: ElementId,
        /// Centre.
        centre: Point2,
        /// Nonzero first semi-axis vector.
        axis: Point2,
        /// Positive ratio of the second semi-axis to the first.
        ratio: Number,
        /// End point.
        to: Point2,
        /// Direction seen from the sketch normal.
        turn: Turn,
    },
}

impl Segment {
    /// The segment's id.
    pub fn id(&self) -> &ElementId {
        match self {
            Segment::Line { id, .. }
            | Segment::Bezier { id, .. }
            | Segment::Arc { id, .. }
            | Segment::Circle { id, .. }
            | Segment::Ellipse { id, .. }
            | Segment::EllipticArc { id, .. } => id,
        }
    }
}

/// The way an arc runs round its centre.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Turn {
    /// Counterclockwise.
    Ccw,
    /// Clockwise.
    Cw,
}

/// One piece of a path, in model space.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PathSegment {
    /// A straight segment to `to`.
    Line {
        /// The segment's id.
        id: ElementId,
        /// End point.
        to: Point3,
        /// Positive half-turn axis at a tangent reversal, in version 4.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        corner_axis: Option<Vector3>,
    },
    /// A circular arc about `centre` to `to`, counterclockwise about `axis`.
    Arc {
        /// The segment's id.
        id: ElementId,
        /// Centre.
        centre: Point3,
        /// End point.
        to: Point3,
        /// Normal of the arc's plane; the arc turns counterclockwise about it.
        axis: Vector3,
        /// Positive half-turn axis at a tangent reversal, in version 4.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        corner_axis: Option<Vector3>,
    },
    /// A clamped rational B-spline in model space, in versions 4 and 5.
    Spline {
        /// The authored segment id, retained across internal knot spans.
        id: ElementId,
        /// Polynomial degree, at least one.
        degree: usize,
        /// Exact controls, including the preceding endpoint and this segment's end.
        control_points: Vec<Point3>,
        /// One strictly positive rational weight per control.
        weights: Vec<Number>,
        /// Expanded rational knots on [0,1], clamped at both endpoints.
        knots: Vec<Number>,
        /// Positive half-turn axis at a tangent reversal before this segment.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        corner_axis: Option<Vector3>,
    },
}

impl PathSegment {
    /// The segment's id.
    pub fn id(&self) -> &ElementId {
        match self {
            Self::Line { id, .. } | Self::Arc { id, .. } | Self::Spline { id, .. } => id,
        }
    }
}

impl PathSegment {
    /// Optional saved axis for a stationary turn at the preceding join.
    pub fn corner_axis(&self) -> Option<&Vector3> {
        match self {
            Self::Line { corner_axis, .. }
            | Self::Arc { corner_axis, .. }
            | Self::Spline { corner_axis, .. } => corner_axis.as_ref(),
        }
    }
}

/// Where one hole starts and which way it goes.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HolePlacement {
    /// The placement's id.
    pub id: ElementId,
    /// The point on the hole's axis where depth zero is.
    pub point: Point3,
    /// The drilling direction, into the material, of any length.
    pub direction: Vector3,
}

/// How deep a hole's main wall goes.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HoleDepth {
    /// To this depth below the placement point.
    Blind(Number),
    /// Without end.
    Through,
}

/// The bottom of a blind hole.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HoleTip {
    /// A flat floor.
    Flat,
    /// A drill-point cone with this full angle in degrees, in `(0, 180)`.
    Cone(Number),
}

/// A counterbore or countersink at a hole's entry.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum HoleHead {
    /// A wider cylinder at the entry.
    Counterbore {
        /// Its diameter, more than the hole's.
        diameter: Number,
        /// Its depth.
        depth: Number,
    },
    /// A cone at the entry.
    Countersink {
        /// Its diameter at depth zero, more than the hole's.
        diameter: Number,
        /// Its full angle in degrees, in `(0, 180)`.
        angle: Number,
    },
}

/// A regularised set operation.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BooleanKind {
    /// Material in either body.
    Union,
    /// Material in the target and not in the tool.
    Subtract,
    /// Material in both bodies.
    Intersect,
}

/// Which side of thickened faces the new solid lies on.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ThickenSide {
    /// Along the faces' outward normals, away from the body.
    Front,
    /// Against them, into the body.
    Back,
}

/// How a sweep's profile turns along the path.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SweepOrientation {
    /// It turns only as much as the path's tangent turns.
    RotationMinimising,
    /// It does not turn; it is only translated.
    Fixed,
}

/// One profile of a loft.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LoftSection {
    /// The sketch step.
    pub profile: StepId,
    /// One segment in each loop, where the loop's correspondence starts. Loops of
    /// different sections correspond by their place in this list.
    pub start: Vec<ElementId>,
}

/// How a loft's side faces pass between sections.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LoftShape {
    /// Straight lines between consecutive sections.
    Ruled,
    /// A natural cubic spline through all sections.
    Smooth,
}

/// Where a pattern's copies go. Copy `k` runs from 0, the body itself.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum PatternLayout {
    /// Copy `k` is translated by `k` times `vector`.
    Linear {
        /// The step between copies.
        vector: Vector3,
        /// Number of copies, including the body itself.
        count: u32,
    },
    /// Copy `k` is turned by `k` times `angle` about `axis`.
    Circular {
        /// The axis.
        axis: Axis,
        /// Number of copies, including the body itself.
        count: u32,
        /// The angle between copies, in degrees.
        angle: Number,
    },
}

impl PatternLayout {
    /// Number of copies, including the body itself.
    pub fn count(&self) -> u32 {
        match self {
            PatternLayout::Linear { count, .. } | PatternLayout::Circular { count, .. } => *count,
        }
    }
}

/// The cell of a lattice, in coordinates where the cell is the unit cube.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum LatticeCell {
    /// The 12 edges of the cube.
    Cubic,
    /// The centre joined to the 8 corners.
    BodyCentred,
    /// Each face centre joined to the 4 corners of its face.
    FaceCentred,
    /// The face-centred beams and each face centre joined to its 4 neighbours.
    Octet,
    /// Explicit nodes and beams.
    Custom {
        /// Node positions.
        nodes: Vec<Point3>,
        /// Beams, each a pair of indices into `nodes`.
        beams: Vec<[u32; 2]>,
    },
}

/// A value that may vary over space.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Grade {
    /// The same everywhere.
    Constant(Number),
    /// `start` at `from`, `end` at `to`, linear in between along the line through them,
    /// constant beyond.
    Linear {
        /// Where the value is `start`.
        from: Point3,
        /// Where the value is `end`.
        to: Point3,
        /// The value at `from`.
        start: Number,
        /// The value at `to`.
        end: Number,
    },
    /// `near` on the faces, `far` at `span` from them and beyond, linear in between.
    Distance {
        /// Faces of the region body.
        faces: Vec<FaceRef>,
        /// The value on the faces.
        near: Number,
        /// The value at `span` and beyond.
        far: Number,
        /// The distance over which the value changes.
        span: Number,
    },
}

/// Whether a gyroid is a sheet or a network.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GyroidKind {
    /// Points where the absolute value of the gyroid function is at most the level.
    Sheet,
    /// Points where the gyroid function is at most the level.
    Network,
}

/// A file the recipe reads, fixed by its content.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FileRef {
    /// The file's name, for people and for finding it.
    pub name: String,
    /// SHA-256 of the file's bytes, 64 lowercase hexadecimal digits.
    pub sha256: String,
}
