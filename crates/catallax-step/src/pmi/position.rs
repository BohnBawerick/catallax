//! Position controls that a plain geometric tolerance cannot state: a
//! composite frame of two ordered segments, a projected tolerance zone, or a
//! datum reference with its own material boundary.

use catallax_numbers::Rational;

use super::{Feature, MaterialCondition};

/// One datum reference of a position segment.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct DatumRef {
    /// The datum letter, such as `A`.
    pub letter: String,
    /// The datum feature's material boundary: maximum (MMB) or least (LMB).
    /// `None` is regardless of material boundary. It belongs to this
    /// reference alone: a later segment naming the same letter does not
    /// inherit it.
    pub material: Option<MaterialCondition>,
}

impl DatumRef {
    /// A reference regardless of material boundary.
    pub fn new(letter: &str) -> Self {
        Self {
            letter: letter.to_owned(),
            material: None,
        }
    }
}

/// One segment of a position control's frame.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct PositionSegment {
    /// The positive total size of the zone in millimetres: its diameter, not
    /// its radius, when `diameter_zone` is set.
    pub value_mm: Rational,
    /// The zone is a cylinder. The axis controls of this subset require it.
    pub diameter_zone: bool,
    /// The controlled features' material condition for this segment only.
    /// `None` is regardless of feature size.
    pub material: Option<MaterialCondition>,
    /// The datum references, primary first. At most three.
    pub datums: Vec<DatumRef>,
}

/// The frame of a position control.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum PositionFrame<F> {
    /// One segment, with or without a projected zone.
    Single {
        /// The segment.
        segment: PositionSegment,
        /// The zone projected beyond the part, if any.
        projected_zone: Option<ProjectedZone<F>>,
    },
    /// A composite frame of exactly two segments, in drawing order. The upper
    /// one locates the pattern from its datums; the lower one refines the
    /// features' relation to each other. It is not a second location control.
    Composite {
        /// The pattern-locating segment.
        upper: PositionSegment,
        /// The feature-relating segment.
        lower: PositionSegment,
    },
}

impl<F> PositionFrame<F> {
    /// The segments in drawing order, upper first.
    pub fn segments(&self) -> Vec<&PositionSegment> {
        match self {
            PositionFrame::Single { segment, .. } => vec![segment],
            PositionFrame::Composite { upper, lower } => vec![upper, lower],
        }
    }
}

/// A tolerance zone projected out of the part from where each controlled
/// feature enters it.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ProjectedZone<F> {
    /// The positive height of the zone in millimetres, along `direction`.
    pub height_mm: Rational,
    /// The zone's direction in part axes, pointing out of the part from the
    /// entry faces into free space. Not a length; only its sense and line
    /// matter.
    pub direction: [Rational; 3],
    /// One planar entry face per controlled feature, in the order of the
    /// control's features. Where the feature's axis meets it, the zone starts.
    pub entry_faces: Vec<F>,
}

/// A position control on one or more cylindrical features. Every segment
/// governs the same features.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct PositionControl<F> {
    /// The controlled features. Several are a pattern.
    pub features: Vec<Feature<F>>,
    /// The frame.
    pub frame: PositionFrame<F>,
}

impl<F> PositionControl<F> {
    /// Whether any segment refers to a datum with a material boundary.
    pub fn has_datum_boundary(&self) -> bool {
        self.frame
            .segments()
            .iter()
            .any(|segment| segment.datums.iter().any(|datum| datum.material.is_some()))
    }
}
