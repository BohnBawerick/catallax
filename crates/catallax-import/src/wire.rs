use serde::{Deserialize, Serialize};

/// Format tag of [`ProposalSet`].
pub const PROPOSALS_FORMAT: &str = "catallax-import-proposals/1";
/// Format tag of [`RepairRecord`].
pub const REPAIR_FORMAT: &str = "catallax-import-repair/2";
/// Repair records with an explicit interpretation of inconsistent planar source data.
pub const INCIDENCE_REPAIR_FORMAT: &str = "catallax-import-repair/3";
/// Repair records with an explicit complete-cylinder source interpretation.
pub const CYLINDER_REPAIR_FORMAT: &str = "catallax-import-repair/4";
/// Proposals that include bounded planar source-incidence repair.
pub const INCIDENCE_PROPOSALS_FORMAT: &str = "catallax-import-proposals/2";
/// Largest gap repaired without asking, in millimetres.
pub const AUTO_CAP_MM: &str = "0.01";

/// What a measured gap is.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GapKind {
    /// A planar carrier or straight boundary disagrees with its source vertices.
    PlanarIncidence,
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

/// Exact decimal or rational axis-aligned bounds,
/// `[[min x, min y, min z], [max x, max y, max z]]`.
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

/// A whole trimmed-face distance enclosure against the retained source definition.
/// Convex polygon covers and complete cylinder bands prove both directions.
/// `/3` records identify their projected source in [`SourceIncidence`]; `/4`
/// records identify their cylindrical source in [`CylinderRadiusRepair`].
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct FaceTolerance {
    /// Source solid entity.
    pub solid: String,
    /// Source face entity, also the imported origin key.
    pub face: String,
    /// Certified lower bound on the whole-face Hausdorff deviation, in mm.
    pub deviation_lower_mm: String,
    /// Proven Hausdorff upper bound in millimetres, retained as local tolerance.
    pub deviation_upper_mm: String,
    /// File approval limit, separate from the attained deviation.
    pub approved_tolerance_mm: String,
    /// Conservative nearest-feature size squared, an exact rational in mm².
    pub feature_size_squared_mm2: String,
}

/// How an inconsistent planar source was made exactly incident.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PlanarRepairMethod {
    /// Keep shared vertices and refit planar carriers through their boundary wires.
    RefitPlanes,
    /// Keep carriers and project each shared vertex onto all its incident planes.
    SnapSharedVertices,
}

/// Source interpretation and complete face bounds for one planar incidence repair.
///
/// The source face is its vertex boundary projected onto its declared plane.
/// The measured discrepancy covers the supplied vertices and finite line-carrier
/// pieces separately from geometric deviation against that projected patch.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct SourceIncidence {
    /// Source solid entity.
    pub solid: String,
    /// Source face entity.
    pub face: String,
    /// Candidate construction, selected deterministically for the solid.
    pub method: PlanarRepairMethod,
    /// Versioned source interpretation, always `projected_vertex_boundary/1`.
    pub source_definition: String,
    /// Lower bound on maximum source-to-projection correspondence displacement.
    pub source_discrepancy_lower_mm: String,
    /// Upper bound over the complete source fan and finite boundary pieces.
    pub source_discrepancy_upper_mm: String,
    /// Whole-face Hausdorff lower bound from the projected source to the result.
    pub geometric_deviation_lower_mm: String,
    /// Whole-face Hausdorff upper bound from the projected source to the result.
    pub geometric_deviation_upper_mm: String,
    /// Direct bound against source constraints, not a sum of the preceding bounds.
    pub source_to_result_upper_mm: String,
    /// Maximum movement of a supplied vertex, separate from face deviation.
    pub vertex_movement_upper_mm: String,
}

/// Measured, typed explanation of an unfinished planar incidence repair.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct PlanarRepairDiagnostic {
    /// Source solid entity.
    pub solid: String,
    /// Source face that first prevented the repair.
    pub face: String,
    /// Stable stage, such as `source_cover`, `import_cap` or `nearest_feature`.
    pub stage: String,
    /// Quantity enclosed; a point residual is never called a whole-face proof.
    pub measurement: String,
    /// Certified lower endpoint when available.
    pub deviation_lower_mm: Option<String>,
    /// Certified upper endpoint when available.
    pub deviation_upper_mm: Option<String>,
    /// The unchanged automatic import cap, in millimetres.
    pub limit_mm: String,
    /// Proven neighbouring-feature size when that check determined the refusal.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_size_squared_mm2: Option<String>,
    /// Why a complete repair was not returned.
    pub reason: String,
}

/// Certified distance enclosure for one failed source incidence check.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct IncidenceGap {
    /// Source edge entity.
    pub edge: String,
    /// Source face entity for an edge/surface check.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub face: Option<String>,
    /// The failed certification step.
    pub step: String,
    /// The geometric quantity enclosed; support distance does not certify a trim.
    pub metric: String,
    /// Proven lower bound in millimetres.
    pub distance_lower_mm: String,
    /// Proven upper bound in millimetres.
    pub distance_upper_mm: String,
    /// Whether both bounds describe the finite edge or an enclosing carrier.
    pub coverage: String,
    /// Whether the requested measurement width was achieved.
    pub converged: bool,
}

/// Source interpretation for an automatic radius repair. The circular rims and
/// cap planes stay unchanged; the source is the complete cylindrical band at
/// its declared radius between those planes, as named by `source_definition`.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct CylinderRadiusRepair {
    /// Source solid entity.
    pub solid: String,
    /// Source cylindrical face entity.
    pub face: String,
    /// Versioned interpretation, `cylinder_between_cap_planes/1`.
    pub source_definition: String,
    /// Declared radius in millimetres, as an exact rational.
    pub source_radius_mm: String,
    /// Common radius of the two unchanged rims, in millimetres.
    pub repaired_radius_mm: String,
    /// Source edges defining the two complete rims, in face order.
    pub unchanged_rim_edges: Vec<String>,
}

/// The persisted change record. Stored beside the recipe and rebuilt into the
/// derivative from the source and this record alone.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct RepairRecord {
    /// [`REPAIR_FORMAT`], [`INCIDENCE_REPAIR_FORMAT`] or [`CYLINDER_REPAIR_FORMAT`].
    pub format: String,
    /// The immutable source.
    pub source: SourceRef,
    /// The repaired derivative. Absent when every change was declined.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub derivative: Option<DerivativeRef>,
    /// Changes and declines, in deterministic application order.
    pub changes: Vec<Change>,
    /// Nonzero whole-face exceptions against source geometry, not output error.
    /// Empty for clean sources and declined repairs.
    pub face_tolerances: Vec<FaceTolerance>,
    /// Source interpretation and measured discrepancy for planar repairs.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_incidence: Vec<SourceIncidence>,
    /// Explicit cylindrical source interpretation, present only in `/4` records.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cylinder_radius: Vec<CylinderRadiusRepair>,
}

/// Typed refusals of an import. None of these is an impossible-geometry claim
/// unless named so.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(tag = "code", rename_all = "snake_case")]
pub enum ImportRefusal {
    /// A measured source incidence has no certified repair within the import rules.
    IncidenceGap {
        /// Attained distance enclosure, never a sampled acceptance.
        gap: Box<IncidenceGap>,
        /// Which repair or admission certificate is missing.
        reason: String,
        /// The strict solid check reached by a bounded exact derivative, if any.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        next_refusal: Option<Box<ImportRefusal>>,
    },
    /// Bounded planar repair could not meet its source or complete-proof requirements.
    PlanarRepairNotVerified {
        /// Face, stage, measurement enclosure and unchanged cap.
        #[serde(flatten)]
        detail: Box<PlanarRepairDiagnostic>,
    },
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
    /// A pair of faces could not be decided, so the solid is not verified. Never a pass.
    SelfIntersectionNotVerified {
        /// Solid entity name.
        solid: String,
        /// Entity names of the first pair of faces that could not be decided.
        faces: [String; 2],
        /// Why the engine could not decide, as `code: detail`.
        reason: String,
        /// The engine's refusal: a `numbers.*` code that keeps the undecided
        /// decision and its bound when a decision ran out of work or
        /// precision, `import.pair-undecided` when the pair is outside what
        /// the check covers. Kept in process; on the wire only as `reason`.
        #[serde(skip)]
        refusal: Option<Box<catallax_numbers::Refusal>>,
    },
    /// The work allowance ended first. Incomplete, not impossible.
    ResourceLimit,
}
