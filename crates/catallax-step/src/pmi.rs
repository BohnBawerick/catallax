//! Manufacturing information on a part, written to STEP as AP242 semantic PMI.
//!
//! A [`Pmi`] is a list of [`Item`]s: datums, dimensions with their limits,
//! geometric tolerances and surface textures, each attached to the faces it
//! applies to, plus the notes that belong to the whole part. The type is
//! generic over the face key: a writer takes model faces, a reader returns
//! the entity names of the file's faces.
//!
//! A dimension carries no nominal value on the way out. The writer measures
//! it from the exact geometry ([`nominals`]), so a file never states a size
//! the solid does not have. A reader returns the value the file states.
//!
//! The entity forms follow the CAx-IF recommended practices for the
//! representation of PMI in AP242, version 4.1.
//! Width, chamfer and edge break use dimensional sizes with explicit
//! side roles. Chamfer angle and edge break bounds use named properties.

mod callout;
mod check;
mod measure;
pub mod patch;
mod position;
mod target;
mod view;
mod write;

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use catallax_numbers::{Rational, Refusal};

pub use check::check_on_model;
pub use measure::nominals;
pub use position::{DatumRef, PositionControl, PositionFrame, PositionSegment, ProjectedZone};
pub use target::{DatumTarget, TargetArea, TargetKind, off_face};
pub use view::{Camera, SavedView};
pub(crate) use write::emit;
pub use write::{
    BEVEL_ANGLE, CHAMFER, COMPOSITE, COMPOSITE_STANDARD, COUNTERBORE, COUNTERSINK, DATUM_TARGET,
    EDGE_BREAK, FIRST_SIDE, HEAD_DEPTH, HEAD_DIAMETER, HOLE_DEPTH, HOLE_PILOT, INCLUDED_ANGLE,
    KNURL, MAX_SETBACK, MIN_SETBACK, ORIENTATION, SECOND_SIDE, SEMANTIC_TEXT, TARGET_DIAMETER,
    TARGET_LENGTH, TARGET_WIDTH, THREAD, THREAD_DEPTH, THROUGH, WIDTH,
};

/// The faces of one feature: a bore's cylinder, or the two flats of a slot.
pub type Feature<F> = Vec<F>;

/// The manufacturing information of one part.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Pmi<F> {
    /// The items attached to faces, in the order the part lists them.
    pub items: Vec<Item<F>>,
    /// What applies to the whole part.
    pub part: PartNotes,
    /// The saved views, in the order the part lists them.
    pub views: Vec<SavedView>,
}

impl<F> Default for Pmi<F> {
    fn default() -> Self {
        Self {
            items: Vec::new(),
            part: PartNotes::default(),
            views: Vec::new(),
        }
    }
}

impl<F> Pmi<F> {
    /// True when there is nothing to write.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty() && self.part.is_empty() && self.views.is_empty()
    }

    /// The same information with every face key replaced. The first face
    /// `map` has no key for refuses with that face's item.
    pub fn try_map<G, E>(
        &self,
        mut map: impl FnMut(&Item<F>, &F) -> Result<G, E>,
    ) -> Result<Pmi<G>, E> {
        let mut items = Vec::with_capacity(self.items.len());
        for item in &self.items {
            let mut feature = |faces: &Feature<F>| -> Result<Feature<G>, E> {
                faces.iter().map(|face| map(item, face)).collect()
            };
            let kind = match &item.kind {
                ItemKind::Datum(datum) => ItemKind::Datum(Datum {
                    letter: datum.letter.clone(),
                    feature: feature(&datum.feature)?,
                }),
                ItemKind::Size(size) => ItemKind::Size(Size {
                    kind: size.kind.clone(),
                    features: size
                        .features
                        .iter()
                        .map(&mut feature)
                        .collect::<Result<_, _>>()?,
                    limits: size.limits.clone(),
                    through: size.through,
                    thread: size.thread.clone(),
                    nominal: size.nominal.clone(),
                }),
                ItemKind::Location(location) => {
                    let mut anchor = |anchor: &Anchor<F>| -> Result<Anchor<G>, E> {
                        Ok(Anchor {
                            feature: feature(&anchor.feature)?,
                            at: anchor.at,
                        })
                    };
                    ItemKind::Location(Location {
                        from: anchor(&location.from)?,
                        to: anchor(&location.to)?,
                        limits: location.limits.clone(),
                        nominal: location.nominal.clone(),
                    })
                }
                ItemKind::Width {
                    flats: [a, b],
                    limits,
                    nominal,
                } => ItemKind::Width {
                    flats: [feature(a)?, feature(b)?],
                    limits: limits.clone(),
                    nominal: nominal.clone(),
                },
                ItemKind::Chamfer {
                    bevel,
                    sides: [a, b],
                    nominal,
                } => ItemKind::Chamfer {
                    bevel: feature(bevel)?,
                    sides: [feature(a)?, feature(b)?],
                    nominal: nominal.clone(),
                },
                ItemKind::EdgeBreak {
                    bevel,
                    sides: [a, b],
                    min_mm,
                    max_mm,
                    nominal,
                } => ItemKind::EdgeBreak {
                    bevel: feature(bevel)?,
                    sides: [feature(a)?, feature(b)?],
                    min_mm: min_mm.clone(),
                    max_mm: max_mm.clone(),
                    nominal: nominal.clone(),
                },
                ItemKind::BasicAngle {
                    planes: [a, b],
                    nominal,
                } => ItemKind::BasicAngle {
                    planes: [feature(a)?, feature(b)?],
                    nominal: nominal.clone(),
                },
                ItemKind::Tolerance(tolerance) => ItemKind::Tolerance(GeometricTolerance {
                    characteristic: tolerance.characteristic,
                    value: tolerance.value.clone(),
                    diameter_zone: tolerance.diameter_zone,
                    material: tolerance.material,
                    datums: tolerance.datums.clone(),
                    features: tolerance
                        .features
                        .iter()
                        .map(&mut feature)
                        .collect::<Result<_, _>>()?,
                }),
                ItemKind::Texture(texture) => ItemKind::Texture(SurfaceTexture {
                    parameter: texture.parameter.clone(),
                    value: texture.value.clone(),
                    process: texture.process,
                    feature: feature(&texture.feature)?,
                }),
                ItemKind::FaceNote {
                    text,
                    feature: faces,
                } => ItemKind::FaceNote {
                    text: text.clone(),
                    feature: feature(faces)?,
                },
                ItemKind::Knurl {
                    designation,
                    feature: faces,
                } => ItemKind::Knurl {
                    designation: designation.clone(),
                    feature: feature(faces)?,
                },
                ItemKind::DatumTarget(target) => ItemKind::DatumTarget(DatumTarget {
                    datum: target.datum.clone(),
                    number: target.number,
                    face: map(item, &target.face)?,
                    point_mm: target.point_mm.clone(),
                    kind: target.kind.clone(),
                }),
                ItemKind::Position(control) => ItemKind::Position(PositionControl {
                    features: control
                        .features
                        .iter()
                        .map(&mut feature)
                        .collect::<Result<_, _>>()?,
                    frame: match &control.frame {
                        PositionFrame::Single {
                            segment,
                            projected_zone,
                        } => PositionFrame::Single {
                            segment: segment.clone(),
                            projected_zone: match projected_zone {
                                Some(zone) => Some(ProjectedZone {
                                    height_mm: zone.height_mm.clone(),
                                    direction: zone.direction.clone(),
                                    entry_faces: zone
                                        .entry_faces
                                        .iter()
                                        .map(|face| map(item, face))
                                        .collect::<Result<_, _>>()?,
                                }),
                                None => None,
                            },
                        },
                        PositionFrame::Composite { upper, lower } => PositionFrame::Composite {
                            upper: upper.clone(),
                            lower: lower.clone(),
                        },
                    },
                }),
            };
            items.push(Item {
                id: item.id.clone(),
                kind,
            });
        }
        Ok(Pmi {
            items,
            part: self.part.clone(),
            views: self.views.clone(),
        })
    }
}

/// One piece of manufacturing information on some faces.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Item<F> {
    /// The item's name within the part, such as `D1`. Unique, and kept by
    /// the file, so a reader can tell which item is which.
    pub id: String,
    /// What the item says.
    pub kind: ItemKind<F>,
}

/// What an item says.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum ItemKind<F> {
    /// A datum feature and its letter.
    Datum(Datum<F>),
    /// The size of a feature: a diameter or a radius.
    Size(Size<F>),
    /// The distance between two features.
    Location(Location<F>),
    /// The width between opposed planar face groups.
    Width {
        /// The two opposed flats.
        flats: [Feature<F>; 2],
        /// The width tolerance, if stated.
        limits: Option<Limits>,
        /// The width in millimetres stated by an imported file.
        nominal: Option<Rational>,
    },
    /// The setback on the reference side and its angle to the bevel.
    Chamfer {
        /// The planar bevel.
        bevel: Feature<F>,
        /// The two adjacent planar sides, reference side first.
        sides: [Feature<F>; 2],
        /// The setback and angle stated by an imported file.
        nominal: Option<ChamferNominal>,
    },
    /// Explicit bounds on the larger setback of a planar broken edge.
    EdgeBreak {
        /// The planar break face.
        bevel: Feature<F>,
        /// The two adjacent planar sides.
        sides: [Feature<F>; 2],
        /// Minimum setback in millimetres.
        min_mm: Rational,
        /// Maximum setback in millimetres.
        max_mm: Rational,
        /// The larger setback in millimetres stated by an imported file.
        nominal: Option<Rational>,
    },
    /// The basic included angle between two planes.
    BasicAngle {
        /// The two planar face groups.
        planes: [Feature<F>; 2],
        /// The angle in degrees stated by an imported file.
        nominal: Option<Rational>,
    },
    /// A geometric tolerance.
    Tolerance(GeometricTolerance<F>),
    /// A surface texture requirement.
    Texture(SurfaceTexture<F>),
    /// A datum target on one face of a declared datum.
    DatumTarget(DatumTarget<F>),
    /// Free text a person reads, attached to its faces. Written as CAx-IF
    /// semantic text: an editable string, not a computable requirement.
    FaceNote {
        /// The note, one or more lines.
        text: String,
        /// The faces it is attached to.
        feature: Feature<F>,
    },
    /// A knurl designation on its faces, such as `RGE 0.8 DIN 82`. Written
    /// as CAx-IF semantic text under the variable name [`KNURL`]: AP242 has
    /// no knurl entity, so no reader learns a pitch or pattern from it.
    Knurl {
        /// The designation as a drawing writes it.
        designation: String,
        /// The knurled faces.
        feature: Feature<F>,
    },
    /// A position control with a composite frame, a projected zone or a
    /// datum reference at a material boundary. A plain single position is a
    /// [`GeometricTolerance`], and this variant refuses one.
    Position(PositionControl<F>),
}

impl<F> ItemKind<F> {
    /// Every feature the item names, in order.
    pub fn features(&self) -> Vec<&Feature<F>> {
        match self {
            ItemKind::Datum(datum) => vec![&datum.feature],
            ItemKind::Size(size) => size.features.iter().collect(),
            ItemKind::Location(location) => vec![&location.from.feature, &location.to.feature],
            ItemKind::Width { flats: [a, b], .. } | ItemKind::BasicAngle { planes: [a, b], .. } => {
                vec![a, b]
            }
            ItemKind::Chamfer {
                bevel,
                sides: [a, b],
                ..
            }
            | ItemKind::EdgeBreak {
                bevel,
                sides: [a, b],
                ..
            } => vec![bevel, a, b],
            ItemKind::Tolerance(tolerance) => tolerance.features.iter().collect(),
            ItemKind::Texture(texture) => vec![&texture.feature],
            ItemKind::Position(control) => control.features.iter().collect(),
            ItemKind::FaceNote { feature, .. } | ItemKind::Knurl { feature, .. } => vec![feature],
            ItemKind::DatumTarget(_) => Vec::new(),
        }
    }

    /// Whether the item is text a person reads rather than semantic PMI. Such
    /// items are counted apart (recommended practices, section 7.4.3).
    pub fn is_text(&self) -> bool {
        matches!(self, ItemKind::FaceNote { .. } | ItemKind::Knurl { .. })
    }

    /// Every face the item names: its features' faces, a datum target's face
    /// and a projected zone's entry faces.
    pub fn faces(&self) -> Vec<&F> {
        let mut faces: Vec<&F> = self.features().into_iter().flatten().collect();
        match self {
            ItemKind::DatumTarget(target) => faces.push(&target.face),
            ItemKind::Position(PositionControl {
                frame:
                    PositionFrame::Single {
                        projected_zone: Some(zone),
                        ..
                    },
                ..
            }) => faces.extend(&zone.entry_faces),
            _ => {}
        }
        faces
    }
}

/// The two values of an imported chamfer dimension.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ChamferNominal {
    /// Setback in millimetres along the first side.
    pub setback_mm: Rational,
    /// Acute angle in degrees between the first side and the bevel.
    pub angle_deg: Rational,
}

/// A datum feature.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Datum<F> {
    /// The datum letter, such as `A`.
    pub letter: String,
    /// The feature that establishes it.
    pub feature: Feature<F>,
}

/// The size of one feature, or of every feature of a pattern.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Size<F> {
    /// Which size.
    pub kind: SizeKind,
    /// For [`SizeKind::Hole`], the pilot followed by the head if present.
    /// Otherwise each feature is an equal member of the dimension's pattern.
    pub features: Vec<Feature<F>>,
    /// The limits, or `None` for a size under the part's general tolerance.
    pub limits: Option<Limits>,
    /// The size is a hole's diameter and the hole goes through the part.
    pub through: bool,
    /// The feature is threaded: its designation as a drawing writes it, such as
    /// `M6x1 - 6H`. The diameter is then the thread's nominal diameter or the
    /// diameter of the hole it is cut in, whichever the faces have.
    pub thread: Option<String>,
    /// The nominal in millimetres a file states (the pilot diameter for a
    /// hole). `None` on the way out: the writer measures it.
    pub nominal: Option<Rational>,
}

/// Which size of a feature a dimension gives.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum SizeKind {
    /// The diameter of a cylindrical feature.
    Diameter,
    /// The radius of a cylindrical feature.
    Radius,
    /// A single hole with a compound head, a known drill depth or a known
    /// usable thread depth. Simple diameters retain their existing encoding.
    Hole(Box<HoleSize>),
}

/// Design intent not captured by the pilot diameter. Lengths are in millimetres.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct HoleSize {
    /// The compound head, if the pilot has one.
    pub head: Option<HoleHead>,
    /// Axial depth from the entry to the specified drill bottom.
    pub drill_depth: Option<Rational>,
    /// Usable axial thread length, independent of the drill depth.
    pub thread_depth: Option<Rational>,
}

/// The design size of a head attached to its own face set.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum HoleHead {
    /// A cylindrical enlargement at the entry.
    Counterbore {
        /// The head diameter.
        diameter: Rational,
        /// The axial depth from the entry to the counterbore floor.
        depth: Rational,
    },
    /// A conical enlargement at the entry.
    Countersink {
        /// The major diameter at the entry.
        diameter: Rational,
        /// The full included angle, in degrees.
        included_angle: Rational,
    },
}

/// The distance between two features.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Location<F> {
    /// The feature it is measured from.
    pub from: Anchor<F>,
    /// The feature it is measured to.
    pub to: Anchor<F>,
    /// The limits, or `None` for a distance under the general tolerance.
    pub limits: Option<Limits>,
    /// The nominal value in millimetres a file states. `None` on the way
    /// out: the writer measures it.
    pub nominal: Option<Rational>,
}

/// One end of a distance.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Anchor<F> {
    /// The feature.
    pub feature: Feature<F>,
    /// Where on the feature the distance ends.
    pub at: At,
}

/// Where on a feature a distance ends.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum At {
    /// The feature itself: a plane's surface, a cylinder's axis.
    Feature,
    /// Midway between the feature's two faces: the centre plane of a slot's
    /// flats, or the axis midway between two equal parallel cylinders.
    Middle,
    /// The side of a cylindrical feature nearest the other end.
    Inner,
    /// The side of a cylindrical feature farthest from the other end.
    Outer,
}

/// The limits of a dimension.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Limits {
    /// Deviations from the nominal value in millimetres, such as +0.1 and 0.
    Deviations {
        /// The upper deviation.
        upper: Rational,
        /// The lower deviation, below the upper one.
        lower: Rational,
    },
    /// Theoretically exact (basic): the value places a tolerance zone and
    /// takes no limits of its own. A geometric tolerance governs the
    /// feature it locates.
    Basic,
    /// An ISO 286 tolerance class, such as H7.
    Fit {
        /// The fundamental deviation: capital letters for a hole, small
        /// letters for a shaft.
        deviation: String,
        /// The tolerance grade number, such as `7`.
        grade: String,
    },
}

/// A geometric tolerance.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct GeometricTolerance<F> {
    /// The characteristic.
    pub characteristic: Characteristic,
    /// The width of the tolerance zone in millimetres.
    pub value: Rational,
    /// The zone is a cylinder and the value its diameter.
    pub diameter_zone: bool,
    /// The material condition the tolerance applies at.
    pub material: Option<MaterialCondition>,
    /// The datum letters it refers to: primary, secondary, tertiary.
    pub datums: Vec<String>,
    /// The features it applies to. Several features are a pattern.
    pub features: Vec<Feature<F>>,
}

/// A geometric characteristic, ISO 1101.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Characteristic {
    /// Angularity.
    Angularity,
    /// Circular run-out.
    CircularRunout,
    /// Coaxiality.
    Coaxiality,
    /// Concentricity.
    Concentricity,
    /// Cylindricity.
    Cylindricity,
    /// Flatness.
    Flatness,
    /// Profile of a line.
    LineProfile,
    /// Parallelism.
    Parallelism,
    /// Perpendicularity.
    Perpendicularity,
    /// Position.
    Position,
    /// Roundness.
    Roundness,
    /// Straightness.
    Straightness,
    /// Profile of a surface.
    SurfaceProfile,
    /// Symmetry.
    Symmetry,
    /// Total run-out.
    TotalRunout,
}

/// Whether a characteristic refers to datums.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DatumUse {
    /// Never: a form tolerance.
    None,
    /// With or without.
    Optional,
    /// Always, one to three.
    Required,
}

impl Characteristic {
    /// Every characteristic.
    pub const ALL: [Characteristic; 15] = [
        Characteristic::Angularity,
        Characteristic::CircularRunout,
        Characteristic::Coaxiality,
        Characteristic::Concentricity,
        Characteristic::Cylindricity,
        Characteristic::Flatness,
        Characteristic::LineProfile,
        Characteristic::Parallelism,
        Characteristic::Perpendicularity,
        Characteristic::Position,
        Characteristic::Roundness,
        Characteristic::Straightness,
        Characteristic::SurfaceProfile,
        Characteristic::Symmetry,
        Characteristic::TotalRunout,
    ];

    /// The AP242 entity that carries it.
    pub fn entity(self) -> &'static str {
        match self {
            Characteristic::Angularity => "ANGULARITY_TOLERANCE",
            Characteristic::CircularRunout => "CIRCULAR_RUNOUT_TOLERANCE",
            Characteristic::Coaxiality => "COAXIALITY_TOLERANCE",
            Characteristic::Concentricity => "CONCENTRICITY_TOLERANCE",
            Characteristic::Cylindricity => "CYLINDRICITY_TOLERANCE",
            Characteristic::Flatness => "FLATNESS_TOLERANCE",
            Characteristic::LineProfile => "LINE_PROFILE_TOLERANCE",
            Characteristic::Parallelism => "PARALLELISM_TOLERANCE",
            Characteristic::Perpendicularity => "PERPENDICULARITY_TOLERANCE",
            Characteristic::Position => "POSITION_TOLERANCE",
            Characteristic::Roundness => "ROUNDNESS_TOLERANCE",
            Characteristic::Straightness => "STRAIGHTNESS_TOLERANCE",
            Characteristic::SurfaceProfile => "SURFACE_PROFILE_TOLERANCE",
            Characteristic::Symmetry => "SYMMETRY_TOLERANCE",
            Characteristic::TotalRunout => "TOTAL_RUNOUT_TOLERANCE",
        }
    }

    /// The characteristic an AP242 entity name carries.
    pub fn of_entity(name: &str) -> Option<Characteristic> {
        Characteristic::ALL
            .into_iter()
            .find(|characteristic| characteristic.entity() == name)
    }

    /// Whether this characteristic is supported under ASME Y14.5-2018.
    pub fn supports_asme_2018(self) -> bool {
        !matches!(
            self,
            Self::Concentricity | Self::Coaxiality | Self::Symmetry
        )
    }

    /// Whether the characteristic refers to datums.
    pub fn datum_use(self) -> DatumUse {
        match self {
            Characteristic::Cylindricity
            | Characteristic::Flatness
            | Characteristic::Roundness
            | Characteristic::Straightness => DatumUse::None,
            Characteristic::LineProfile
            | Characteristic::Position
            | Characteristic::SurfaceProfile => DatumUse::Optional,
            Characteristic::Angularity
            | Characteristic::CircularRunout
            | Characteristic::Coaxiality
            | Characteristic::Concentricity
            | Characteristic::Parallelism
            | Characteristic::Perpendicularity
            | Characteristic::Symmetry
            | Characteristic::TotalRunout => DatumUse::Required,
        }
    }
}

/// The material condition a geometric tolerance applies at.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum MaterialCondition {
    /// Maximum material requirement.
    Maximum,
    /// Least material requirement.
    Least,
}

/// A surface texture requirement, ISO 1302.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SurfaceTexture<F> {
    /// The roughness parameter, such as `Ra`.
    pub parameter: String,
    /// Its limit in micrometres.
    pub value: Rational,
    /// What the process may do to the surface.
    pub process: Process,
    /// The feature it applies to.
    pub feature: Feature<F>,
}

/// What a surface texture requirement allows the process to do.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub enum Process {
    /// Any process.
    #[default]
    Any,
    /// Material must be removed.
    Removal,
    /// Material must not be removed.
    NoRemoval,
}

impl Process {
    /// The phrase ISO 10303-1110 gives the condition.
    pub fn phrase(self) -> &'static str {
        match self {
            Process::Any => "any process allowed",
            Process::Removal => "material removal required",
            Process::NoRemoval => "no material removal",
        }
    }

    /// The condition a phrase names.
    pub fn of_phrase(phrase: &str) -> Option<Process> {
        [Process::Any, Process::Removal, Process::NoRemoval]
            .into_iter()
            .find(|process| process.phrase() == phrase)
    }
}

/// The explicitly selected dimensioning standard for a part.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DimensioningStandard {
    /// ASME Y14.5-2018.
    AsmeY14_5_2018,
}

/// What applies to the whole part.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct PartNotes {
    /// The dimensioning standard. Omission uses ISO 1101:2017.
    pub standard: Option<DimensioningStandard>,
    /// The general tolerance, such as `ISO 2768-m`.
    pub general_tolerance: Option<String>,
    /// The material, such as `6061-T6`.
    pub material: Option<String>,
    /// Other notes, such as `Break sharp edges`.
    pub notes: Vec<String>,
}

impl PartNotes {
    /// True when the part has no notes or explicit standard.
    pub fn is_empty(&self) -> bool {
        self.standard.is_none()
            && self.general_tolerance.is_none()
            && self.material.is_none()
            && self.notes.is_empty()
    }
}

fn refuse(item: &str, what: impl std::fmt::Display) -> Refusal {
    Refusal::new("step.pmi-invalid", format!("{item}: {what}"))
}

fn positive(value: &Rational) -> bool {
    *value > Rational::zero()
}

/// Why one item did not pass validation. The caller supplies its name or
/// source entity when reporting the error.
#[derive(Debug)]
pub enum ItemError<'a> {
    /// The name already belongs to an item in this part.
    DuplicateName,
    /// The letter already belongs to an accepted datum.
    DuplicateDatum(&'a str),
    /// An accepted datum already names this set of faces.
    DuplicateDatumFeature,
    /// The tolerance references a letter that no accepted datum defines.
    UndefinedDatum(&'a str),
    /// An accepted datum target already has this datum letter and number.
    DuplicateTarget(&'a str, u32),
    /// Another item rule did not pass.
    Rule(&'static str),
}

impl fmt::Display for ItemError<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateName => f.write_str("the name is used twice"),
            Self::DuplicateDatum(letter) => write!(f, "datum {letter} is defined twice"),
            Self::DuplicateDatumFeature => f.write_str("the feature already carries a datum"),
            Self::UndefinedDatum(letter) => write!(f, "datum {letter} is not defined"),
            Self::DuplicateTarget(letter, number) => {
                write!(f, "datum target {letter}{number} is defined twice")
            }
            Self::Rule(reason) => f.write_str(reason),
        }
    }
}

impl<F> Datum<F> {
    /// An authored datum letter is one or two capital letters.
    pub fn validate_letter(&self) -> Result<(), ItemError<'static>> {
        if (1..=2).contains(&self.letter.len())
            && self.letter.bytes().all(|b| b.is_ascii_uppercase())
        {
            Ok(())
        } else {
            Err(ItemError::Rule(
                "a datum letter is one or two capital letters",
            ))
        }
    }
}

/// Item validation against the names and datums accepted for one part.
pub struct PartValidator<'a, F> {
    names: BTreeSet<Cow<'a, str>>,
    datum_features: BTreeSet<Vec<&'a F>>,
    datum_faces: BTreeMap<&'a str, Vec<&'a F>>,
    targets: BTreeSet<(String, u32)>,
}

impl<F> Default for PartValidator<'_, F> {
    fn default() -> Self {
        Self {
            names: BTreeSet::new(),
            datum_features: BTreeSet::new(),
            datum_faces: BTreeMap::new(),
            targets: BTreeSet::new(),
        }
    }
}

impl<'a, F: Ord> PartValidator<'a, F> {
    /// Gives an accepted import item a name, adding a suffix on collision.
    pub fn assign_name(&mut self, base: String) -> String {
        let name = if self.names.contains(base.as_str()) {
            (2..)
                .map(|copy| format!("{base} ({copy})"))
                .find(|name| !self.names.contains(name.as_str()))
                .expect("a finite part has an unused name")
        } else {
            base
        };
        self.names.insert(Cow::Owned(name.clone()));
        name
    }

    /// Checks for conflicts with accepted datums. The letter and face set are
    /// registered together or not at all. The caller checks the datum's local
    /// rules first, since a reader can retain a letter the writer cannot emit.
    pub fn accept_datum(&mut self, datum: &'a Datum<F>) -> Result<(), ItemError<'a>> {
        let letter = datum.letter.as_str();
        if self.contains_datum(letter) {
            return Err(ItemError::DuplicateDatum(letter));
        }
        if !self.datum_features.insert(feature_key(&datum.feature)) {
            return Err(ItemError::DuplicateDatumFeature);
        }
        self.datum_faces
            .insert(letter, datum.feature.iter().collect());
        Ok(())
    }

    /// Checks a datum target against the accepted datums and targets, and
    /// accepts it. Its datum must be accepted first, its face must be one of
    /// that datum's faces, and its letter and number must be new.
    pub fn accept_target<'t>(&mut self, target: &'t DatumTarget<F>) -> Result<(), ItemError<'t>> {
        let letter = target.datum.as_str();
        if target.number == 0 {
            return Err(ItemError::Rule("a datum target's number is positive"));
        }
        let Some(faces) = self.datum_faces.get(letter) else {
            return Err(ItemError::UndefinedDatum(letter));
        };
        if !faces.contains(&&target.face) {
            return Err(ItemError::Rule(
                "the target's face is not on its datum's feature",
            ));
        }
        if !self.targets.insert((letter.to_owned(), target.number)) {
            return Err(ItemError::DuplicateTarget(letter, target.number));
        }
        Ok(())
    }

    /// Whether the part has accepted a datum with this letter.
    pub fn contains_datum(&self, letter: &str) -> bool {
        self.datum_faces.contains_key(letter)
    }

    /// Checks an item's features and values against the accepted datums.
    pub fn check_kind<'item>(&self, kind: &'item ItemKind<F>) -> Result<(), ItemError<'item>> {
        check_features(kind)?;
        self.check_values(kind)
    }

    fn check_values<'item>(&self, kind: &'item ItemKind<F>) -> Result<(), ItemError<'item>> {
        use ItemError::Rule;
        match kind {
            ItemKind::Datum(_) => {}
            ItemKind::Size(size) => {
                if size.features.is_empty() {
                    return Err(Rule("the size names no feature"));
                }
                if size.through && size.kind == SizeKind::Radius {
                    return Err(Rule("only a diameter goes through"));
                }
                if let Some(thread) = &size.thread {
                    if size.kind == SizeKind::Radius {
                        return Err(Rule("only a diameter is threaded"));
                    }
                    // A designation such as `M6x1 - 6H`: plain text a STEP
                    // string holds as it stands.
                    let plain = !thread.trim().is_empty()
                        && thread.len() <= 40
                        && thread.bytes().all(|b| b.is_ascii_graphic() || b == b' ')
                        && !thread.contains(['\'', '\\']);
                    if !plain {
                        return Err(Rule("the thread's designation is not plain text"));
                    }
                }
                if let SizeKind::Hole(hole) = &size.kind {
                    check_hole(size, hole)?;
                }
                check_limits(size.limits.as_ref())?;
            }
            ItemKind::Location(location) => {
                let from = &location.from.feature;
                let to = &location.to.feature;
                if location.from.at == location.to.at
                    && from.len() == to.len()
                    && (from == to || feature_key(from) == feature_key(to))
                {
                    return Err(Rule("the distance starts and ends on one feature"));
                }
                for anchor in [&location.from, &location.to] {
                    if anchor.at == At::Middle && anchor.feature.len() != 2 {
                        return Err(Rule("the middle is taken between two faces"));
                    }
                }
                check_limits(location.limits.as_ref())?;
            }
            ItemKind::Width {
                limits, nominal, ..
            } => {
                if nominal.as_ref().is_some_and(|value| !positive(value)) {
                    return Err(Rule("the width must be more than zero"));
                }
                check_limits(limits.as_ref())?;
            }
            ItemKind::Chamfer { nominal, .. } => {
                if nominal.as_ref().is_some_and(|nominal| {
                    !positive(&nominal.setback_mm)
                        || !positive(&nominal.angle_deg)
                        || nominal.angle_deg > Rational::from_i64(90)
                }) {
                    return Err(Rule(
                        "the chamfer needs a positive setback and an angle up to 90 degrees",
                    ));
                }
            }
            ItemKind::EdgeBreak {
                min_mm,
                max_mm,
                nominal,
                ..
            } => {
                if *min_mm < Rational::zero() || max_mm <= min_mm {
                    return Err(Rule("edge break bounds must be nonnegative and ordered"));
                }
                if nominal.as_ref().is_some_and(|value| !positive(value)) {
                    return Err(Rule("the edge break setback must be more than zero"));
                }
            }
            ItemKind::BasicAngle { nominal, .. } => {
                if nominal
                    .as_ref()
                    .is_some_and(|value| !positive(value) || *value >= Rational::from_i64(180))
                {
                    return Err(Rule("the basic angle must lie between 0 and 180 degrees"));
                }
            }
            ItemKind::Tolerance(tolerance) => {
                if tolerance.features.is_empty() {
                    return Err(Rule("the tolerance names no feature"));
                }
                if !positive(&tolerance.value) {
                    return Err(Rule("the tolerance value must be more than zero"));
                }
                if tolerance.datums.len() > 3 {
                    return Err(Rule("a tolerance refers to at most three datums"));
                }
                for (index, letter) in tolerance.datums.iter().enumerate() {
                    if tolerance.datums[..index].contains(letter) {
                        return Err(Rule("a datum is referred to twice"));
                    }
                }
                for letter in &tolerance.datums {
                    if !self.contains_datum(letter) {
                        return Err(ItemError::UndefinedDatum(letter));
                    }
                }
                match (tolerance.characteristic.datum_use(), tolerance.datums.len()) {
                    (DatumUse::None, 1..) => {
                        return Err(Rule("a form tolerance refers to no datum"));
                    }
                    (DatumUse::Required, 0) => {
                        return Err(Rule("the tolerance needs a datum"));
                    }
                    _ => {}
                }
            }
            ItemKind::Texture(texture) => {
                if texture.parameter.is_empty()
                    || !texture.parameter.bytes().all(|b| b.is_ascii_alphanumeric())
                {
                    return Err(Rule("the roughness parameter is letters and digits"));
                }
                if !positive(&texture.value) {
                    return Err(Rule("the roughness value must be more than zero"));
                }
            }
            ItemKind::DatumTarget(target) => {
                if !target.kind.sizes().into_iter().all(positive) {
                    return Err(Rule("a datum target's sizes must be more than zero"));
                }
                if target.kind.direction().is_some_and(is_zero_vector) {
                    return Err(Rule("a datum target's direction must not be zero"));
                }
            }
            ItemKind::Position(control) => self.check_position(control)?,
            ItemKind::FaceNote { text, .. } => {
                if text.trim().is_empty() {
                    return Err(Rule("a face note needs text"));
                }
            }
            ItemKind::Knurl { designation, .. } => {
                if designation.trim().is_empty() {
                    return Err(Rule("a knurl designation needs text"));
                }
            }
        }
        Ok(())
    }

    fn check_position<'item>(
        &self,
        control: &'item PositionControl<F>,
    ) -> Result<(), ItemError<'item>> {
        use ItemError::Rule;
        if control.features.is_empty() {
            return Err(Rule("the position control names no feature"));
        }
        for segment in control.frame.segments() {
            if !positive(&segment.value_mm) {
                return Err(Rule("a position zone must be more than zero"));
            }
            if !segment.diameter_zone {
                return Err(Rule("a position control's zone is cylindrical"));
            }
            if segment.datums.len() > 3 {
                return Err(Rule("a segment refers to at most three datums"));
            }
            for (index, datum) in segment.datums.iter().enumerate() {
                if segment.datums[..index]
                    .iter()
                    .any(|other| other.letter == datum.letter)
                {
                    return Err(Rule("a datum is referred to twice"));
                }
                if !self.contains_datum(&datum.letter) {
                    return Err(ItemError::UndefinedDatum(&datum.letter));
                }
            }
        }
        match &control.frame {
            PositionFrame::Composite { upper, lower } => {
                if upper.datums.is_empty() {
                    return Err(Rule("the upper segment of a composite needs datums"));
                }
                if lower.datums.len() > upper.datums.len()
                    || lower
                        .datums
                        .iter()
                        .zip(&upper.datums)
                        .any(|(low, up)| low.letter != up.letter)
                {
                    return Err(Rule(
                        "the lower segment's datums must repeat the upper segment's first datums in order",
                    ));
                }
                if lower.value_mm > upper.value_mm {
                    return Err(Rule(
                        "the lower segment's zone cannot be larger than the upper one",
                    ));
                }
            }
            PositionFrame::Single {
                projected_zone: Some(zone),
                ..
            } => {
                if !positive(&zone.height_mm) {
                    return Err(Rule("a projected zone needs a positive height"));
                }
                if is_zero_vector(&zone.direction) {
                    return Err(Rule("a projected zone needs a direction"));
                }
                if zone.entry_faces.len() != control.features.len() {
                    return Err(Rule(
                        "a projected zone needs one entry face per controlled feature",
                    ));
                }
                if zone
                    .entry_faces
                    .iter()
                    .any(|face| control.features.iter().flatten().any(|f| f == face))
                {
                    return Err(Rule("an entry face is one of the controlled faces"));
                }
            }
            PositionFrame::Single {
                projected_zone: None,
                ..
            } => {
                if !control.has_datum_boundary() {
                    return Err(Rule(
                        "a single position with no projected zone and no datum boundary is a geometric tolerance",
                    ));
                }
            }
        }
        Ok(())
    }
}

fn is_zero_vector(vector: &[Rational; 3]) -> bool {
    vector.iter().all(Rational::is_zero)
}

fn feature_key<F: Ord>(faces: &[F]) -> Vec<&F> {
    let mut key: Vec<_> = faces.iter().collect();
    key.sort_unstable();
    key
}

fn check_features<F: Ord>(kind: &ItemKind<F>) -> Result<(), ItemError<'static>> {
    fn check<F: Ord>(feature: &Feature<F>) -> Result<(), ItemError<'static>> {
        match feature.len() {
            0 => Err(ItemError::Rule("a feature has no face")),
            1 => Ok(()),
            _ => {
                let mut unique = BTreeSet::new();
                if feature.iter().all(|face| unique.insert(face)) {
                    Ok(())
                } else {
                    Err(ItemError::Rule("a feature names a face twice"))
                }
            }
        }
    }
    match kind {
        ItemKind::Datum(datum) => check(&datum.feature),
        ItemKind::Size(size) => size.features.iter().try_for_each(check),
        ItemKind::Location(location) => {
            check(&location.from.feature)?;
            check(&location.to.feature)
        }
        ItemKind::Width { .. }
        | ItemKind::Chamfer { .. }
        | ItemKind::EdgeBreak { .. }
        | ItemKind::BasicAngle { .. } => {
            let mut used = BTreeSet::new();
            for feature in kind.features() {
                check(feature)?;
                if !feature.iter().all(|face| used.insert(face)) {
                    return Err(ItemError::Rule("planar dimension roles share a face"));
                }
            }
            Ok(())
        }
        ItemKind::Tolerance(tolerance) => tolerance.features.iter().try_for_each(check),
        ItemKind::Texture(texture) => check(&texture.feature),
        ItemKind::Position(control) => control.features.iter().try_for_each(check),
        ItemKind::FaceNote { feature, .. } | ItemKind::Knurl { feature, .. } => check(feature),
        ItemKind::DatumTarget(_) => Ok(()),
    }
}

impl<F: Ord> Pmi<F> {
    /// Checks that the information is consistent in itself, before any face
    /// is looked up: unique names and datum letters, no empty feature, every
    /// referenced datum defined, values in range. Refuses `step.pmi-invalid`
    /// with the item's name.
    pub fn validate(&self) -> Result<(), Refusal> {
        let mut validator = PartValidator::default();
        // Check names, features and datum declarations before item values.
        // A tolerance can name a datum that appears later in the list.
        for item in &self.items {
            let id = item.id.as_str();
            if id.is_empty() {
                return Err(refuse("an item", "has no name"));
            }
            if !validator.names.insert(Cow::Borrowed(id)) {
                return Err(refuse(id, ItemError::DuplicateName));
            }
            check_features(&item.kind).map_err(|reason| refuse(id, reason))?;
            if let ItemKind::Datum(datum) = &item.kind {
                datum
                    .validate_letter()
                    .map_err(|reason| refuse(id, reason))?;
                validator
                    .accept_datum(datum)
                    .map_err(|reason| refuse(id, reason))?;
            }
        }
        // A target needs its datum, which can come later in the list.
        for item in &self.items {
            if let ItemKind::DatumTarget(target) = &item.kind {
                validator
                    .accept_target(target)
                    .map_err(|reason| refuse(&item.id, reason))?;
            }
        }
        for item in &self.items {
            validator
                .check_values(&item.kind)
                .map_err(|reason| refuse(&item.id, reason))?;
        }
        self.validate_views()?;
        for (what, text) in self
            .part
            .general_tolerance
            .iter()
            .map(|text| ("general tolerance", text))
            .chain(self.part.material.iter().map(|text| ("material", text)))
            .chain(self.part.notes.iter().map(|text| ("note", text)))
        {
            if text.trim().is_empty() {
                return Err(refuse(what, "is empty"));
            }
        }
        Ok(())
    }
}

impl<F> Pmi<F> {
    fn validate_views(&self) -> Result<(), Refusal> {
        let mut names = BTreeSet::new();
        for view in &self.views {
            let name = view.name.as_str();
            if name.trim().is_empty() {
                return Err(refuse("a saved view", "has no name"));
            }
            if !names.insert(name) {
                return Err(refuse(name, "the view name is used twice"));
            }
            let camera = &view.camera;
            let dot = (0..3).fold(Rational::zero(), |sum, i| {
                &sum + &(&camera.up[i] * &camera.direction[i])
            });
            if is_zero_vector(&camera.direction) || is_zero_vector(&camera.up) || !dot.is_zero() {
                return Err(refuse(
                    name,
                    "the camera needs a direction and an up direction square to it",
                ));
            }
            if !positive(&camera.distance) || !camera.window.iter().all(positive) {
                return Err(refuse(
                    name,
                    "the camera needs a positive distance and window",
                ));
            }
            let mut shown = BTreeSet::new();
            for id in &view.items {
                if !shown.insert(id.as_str()) {
                    return Err(refuse(name, format!("item {id} is listed twice")));
                }
                let item = self
                    .items
                    .iter()
                    .find(|item| &item.id == id)
                    .ok_or_else(|| refuse(name, format!("item {id} is not on the part")))?;
                if !item.kind.has_callout() {
                    return Err(refuse(
                        name,
                        format!("item {id} has no callout this writer draws"),
                    ));
                }
            }
        }
        Ok(())
    }
}

impl<F> ItemKind<F> {
    /// Whether this item has a supported saved-view drawing.
    pub fn has_callout(&self) -> bool {
        matches!(
            self,
            ItemKind::Datum(_)
                | ItemKind::Size(Size {
                    kind: SizeKind::Diameter,
                    thread: None,
                    ..
                })
                | ItemKind::Tolerance(GeometricTolerance {
                    characteristic: Characteristic::Position,
                    ..
                })
                | ItemKind::Texture(_)
        )
    }
}

fn check_hole<F: Ord>(size: &Size<F>, hole: &HoleSize) -> Result<(), ItemError<'static>> {
    use ItemError::Rule;

    if size.features.len() != 1 + usize::from(hole.head.is_some()) {
        return Err(Rule("a hole needs one pilot and, if compound, one head"));
    }
    if let [pilot, head] = size.features.as_slice()
        && pilot.iter().any(|face| head.contains(face))
    {
        return Err(Rule("the pilot and head share a face"));
    }
    if hole.head.is_none() && hole.drill_depth.is_none() && hole.thread_depth.is_none() {
        return Err(Rule("a hole needs a head, a drill depth or a thread depth"));
    }
    if let Some(depth) = &hole.drill_depth
        && (size.through || !positive(depth))
    {
        return Err(Rule("a drill depth must be positive and cannot be through"));
    }
    if let Some(thread) = &hole.thread_depth
        && (size.thread.is_none()
            || !positive(thread)
            || hole
                .drill_depth
                .as_ref()
                .is_some_and(|drill| thread > drill))
    {
        return Err(Rule(
            "a thread depth must be positive, designated and no longer than the drill depth",
        ));
    }
    match &hole.head {
        Some(HoleHead::Counterbore { diameter, depth })
            if !positive(diameter)
                || !positive(depth)
                || hole
                    .drill_depth
                    .as_ref()
                    .is_some_and(|drill| depth >= drill) =>
        {
            Err(Rule(
                "the counterbore needs a positive diameter and a depth short of the drill bottom",
            ))
        }
        Some(HoleHead::Countersink {
            diameter,
            included_angle,
        }) if !positive(diameter)
            || !positive(included_angle)
            || *included_angle >= Rational::from_i64(180) =>
        {
            Err(Rule(
                "the countersink needs a positive diameter and an included angle between 0 and 180 degrees",
            ))
        }
        _ => Ok(()),
    }
}

fn check_limits(limits: Option<&Limits>) -> Result<(), ItemError<'static>> {
    match limits {
        None | Some(Limits::Basic) => Ok(()),
        Some(Limits::Deviations { upper, lower }) => {
            if upper > lower {
                Ok(())
            } else {
                Err(ItemError::Rule(
                    "the upper deviation must be above the lower one",
                ))
            }
        }
        Some(Limits::Fit { deviation, grade }) => {
            let letters = (1..=2).contains(&deviation.len())
                && (deviation.bytes().all(|b| b.is_ascii_uppercase())
                    || deviation.bytes().all(|b| b.is_ascii_lowercase()));
            let digits =
                (1..=2).contains(&grade.len()) && grade.bytes().all(|b| b.is_ascii_digit());
            if letters && digits {
                Ok(())
            } else {
                Err(ItemError::Rule(
                    "a tolerance class is one or two letters and a grade number",
                ))
            }
        }
    }
}
