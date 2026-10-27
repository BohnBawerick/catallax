use serde::{Deserialize, Serialize};

/// Format tag of [`ProposalSet`].
pub const PROPOSALS_FORMAT: &str = "catallax-import-proposals/1";
/// Format tag of [`RepairRecord`].
pub const REPAIR_FORMAT: &str = "catallax-import-repair/1";
/// Largest gap repaired without asking, in millimetres.
pub const AUTO_CAP_MM: &str = "0.01";

/// What a measured gap is.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GapKind {
    /// Two boundary features of one shell that the file declares distinct and
    /// that closing would join into one seam.
    SeamClose,
    /// Features of different shells or bodies. Never merged automatically.
    BodyGap,
    /// A boundary loop with no neighbouring face. Never closed automatically.
    Opening,
}

/// Whether the policy table settles a proposal without the user.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Disposition {
    /// Within the automatic cap. Applied and recorded without a question.
    Auto,
    /// Needs an explicit decision before anything moves.
    Ask,
}

/// One measured gap and what the policy says about it.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Proposal {
    /// Stable within the source bytes, such as `"g0"`.
    pub id: String,
    /// Recipe decision key: `repair:<first 12 hex of source sha256>:<id>`.
    pub question: String,
    /// `ADVANCED_FACE` or edge entity names, such as `"#12"`.
    pub entities: Vec<String>,
    /// Kind of gap.
    pub kind: GapKind,
    /// Certified lower bound of the gap distance.
    pub distance_min_mm: String,
    /// Certified upper bound of the gap distance. Equal to the lower bound when exact.
    pub distance_max_mm: String,
    /// Tolerance the move would be validated against.
    pub tolerance_mm: String,
    /// `auto` only when the policy table allows closure without asking.
    pub disposition: Disposition,
    /// Why this disposition, in plain words.
    pub reason: String,
}

/// What an import returns: the proposals for one source file.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct ProposalSet {
    /// [`PROPOSALS_FORMAT`].
    pub format: String,
    /// Lowercase hex SHA-256 of the source bytes.
    pub source_sha256: String,
    /// Every measured gap, ordered by id. Empty when the file is already closed.
    pub proposals: Vec<Proposal>,
}

/// An application's answer to one `ask` proposal.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(tag = "decision", rename_all = "snake_case")]
pub enum Decision {
    /// Allow the move, validated against `tolerance_mm`.
    Accept {
        /// Proposal id.
        id: String,
        /// Approved tolerance, at least the measured upper bound.
        tolerance_mm: String,
        /// Identity the application supplies for the record.
        approver: String,
    },
    /// Leave the geometry as in the file.
    Decline {
        /// Proposal id.
        id: String,
    },
}

/// How a recorded change was approved.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Approval {
    /// Within the cap, no question asked.
    Automatic,
    /// Accepted by a user above the cap or where intent was uncertain.
    User,
    /// Declined. The geometry is unchanged.
    Declined,
}

/// Axis-aligned bounds, `[[min x, min y, min z], [max x, max y, max z]]`.
pub type Bounds = [[String; 3]; 2];

/// One recorded change.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Change {
    /// Proposal id.
    pub id: String,
    /// Affected entity names.
    pub entities: Vec<String>,
    /// Kind of gap.
    pub kind: GapKind,
    /// Body bounds before the change.
    pub bounds_before: Bounds,
    /// Body bounds after the change.
    pub bounds_after: Bounds,
    /// Maximum movement measured on the repaired result, not copied from the proposal.
    pub max_movement_mm: String,
    /// Tolerance the change was validated against.
    pub tolerance_mm: String,
    /// Why the change was made or declined.
    pub reason: String,
    /// How it was approved.
    pub approval: Approval,
    /// Approver identity for [`Approval::User`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub approver: Option<String>,
}

/// Source reference of a record.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct SourceRef {
    /// File name.
    pub name: String,
    /// Lowercase hex SHA-256 of the immutable source bytes.
    pub sha256: String,
}

/// Derivative reference of a record.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct DerivativeRef {
    /// Lowercase hex SHA-256 of the derivative's canonical bytes.
    pub sha256: String,
}

/// The persisted change record. Stored beside the recipe and rebuilt into the
/// derivative from the source and this record alone.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct RepairRecord {
    /// [`REPAIR_FORMAT`].
    pub format: String,
    /// The immutable source.
    pub source: SourceRef,
    /// The repaired derivative. Absent when every change was declined.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub derivative: Option<DerivativeRef>,
    /// Changes and declines, ordered by id.
    pub changes: Vec<Change>,
}

/// Typed refusals of an import. None of these is an impossible-geometry claim
/// unless named so.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(tag = "code", rename_all = "snake_case")]
pub enum ImportRefusal {
    /// Not a readable Part 21 file.
    Malformed {
        /// What failed to parse.
        detail: String,
    },
    /// An entity outside the supported first path.
    Unsupported {
        /// Entity name, such as `"#40"`.
        entity: String,
        /// STEP type name, such as `"B_SPLINE_SURFACE_WITH_KNOTS"`.
        type_name: String,
    },
    /// Length unit missing or not one of millimetre, inch, metre.
    UnknownUnit {
        /// The unit text found, if any.
        found: String,
    },
    /// Supplied bytes do not match the recipe hash. A recipe error.
    HashMismatch {
        /// Hash in the recipe.
        expected: String,
        /// Hash of the bytes.
        actual: String,
    },
    /// The faces do not bound one closed solid after the decisions. Geometry.
    NotSolid {
        /// Why.
        detail: String,
    },
    /// A decision names an unknown or auto proposal, or a tolerance below the measured bound.
    InvalidDecision {
        /// Proposal id.
        id: String,
        /// Why.
        detail: String,
    },
    /// The moved result failed closed, non-self-intersecting, oriented or movement checks.
    ValidationFailed {
        /// Proposal id.
        id: String,
        /// Which check.
        detail: String,
    },
    /// Two faces of one solid pass through each other. Geometry.
    SelfIntersecting {
        /// Solid entity name.
        solid: String,
        /// Entity names of the first pair of faces whose interiors intersect.
        faces: [String; 2],
    },
    /// The work allowance ended first. Incomplete, not impossible.
    ResourceLimit,
}
