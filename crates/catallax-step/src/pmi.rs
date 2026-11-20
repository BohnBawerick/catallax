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

mod measure;
mod write;

use std::collections::BTreeSet;

use catallax_numbers::{Rational, Refusal};

pub use measure::nominals;
pub(crate) use write::emit;
pub use write::{HOLE_DEPTH, THREAD, THROUGH};

/// The faces of one feature: a bore's cylinder, or the two flats of a slot.
pub type Feature<F> = Vec<F>;

/// The manufacturing information of one part.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Pmi<F> {
    /// The items attached to faces, in the order the part lists them.
    pub items: Vec<Item<F>>,
    /// What applies to the whole part.
    pub part: PartNotes,
}

impl<F> Default for Pmi<F> {
    fn default() -> Self {
        Self {
            items: Vec::new(),
            part: PartNotes::default(),
        }
    }
}

impl<F> Pmi<F> {
    /// True when there is nothing to write.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty() && self.part.is_empty()
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
                    kind: size.kind,
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
            };
            items.push(Item {
                id: item.id.clone(),
                kind,
            });
        }
        Ok(Pmi {
            items,
            part: self.part.clone(),
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
    /// A geometric tolerance.
    Tolerance(GeometricTolerance<F>),
    /// A surface texture requirement.
    Texture(SurfaceTexture<F>),
}

impl<F> ItemKind<F> {
    /// Every feature the item names, in order.
    pub fn features(&self) -> Vec<&Feature<F>> {
        match self {
            ItemKind::Datum(datum) => vec![&datum.feature],
            ItemKind::Size(size) => size.features.iter().collect(),
            ItemKind::Location(location) => vec![&location.from.feature, &location.to.feature],
            ItemKind::Tolerance(tolerance) => tolerance.features.iter().collect(),
            ItemKind::Texture(texture) => vec![&texture.feature],
        }
    }
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
    /// The features it applies to. Several features are a pattern, such as
    /// four equal holes, and the size applies to each of them.
    pub features: Vec<Feature<F>>,
    /// The limits, or `None` for a size under the part's general tolerance.
    pub limits: Option<Limits>,
    /// The size is a hole's diameter and the hole goes through the part.
    pub through: bool,
    /// The feature is threaded: its designation as a drawing writes it, such as
    /// `M6x1 - 6H`. The diameter is then the thread's nominal diameter or the
    /// diameter of the hole it is cut in, whichever the faces have.
    pub thread: Option<String>,
    /// The nominal value in millimetres a file states. `None` on the way
    /// out: the writer measures it.
    pub nominal: Option<Rational>,
}

/// Which size of a feature a dimension gives.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum SizeKind {
    /// The diameter of a cylindrical feature.
    Diameter,
    /// The radius of a cylindrical feature.
    Radius,
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

/// What applies to the whole part.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct PartNotes {
    /// The general tolerance, such as `ISO 2768-m`.
    pub general_tolerance: Option<String>,
    /// The material, such as `6061-T6`.
    pub material: Option<String>,
    /// Other notes, such as `Break sharp edges`.
    pub notes: Vec<String>,
}

impl PartNotes {
    /// True when the part has no note at all.
    pub fn is_empty(&self) -> bool {
        self.general_tolerance.is_none() && self.material.is_none() && self.notes.is_empty()
    }
}

fn refuse(item: &str, what: impl std::fmt::Display) -> Refusal {
    Refusal::new("step.pmi-invalid", format!("{item}: {what}"))
}

fn positive(value: &Rational) -> bool {
    *value > Rational::zero()
}

/// A datum letter is one or two capital letters.
pub(crate) fn datum_letter(letter: &str) -> bool {
    (1..=2).contains(&letter.len()) && letter.bytes().all(|b| b.is_ascii_uppercase())
}

impl<F: Ord + Clone> Pmi<F> {
    /// Checks that the information is consistent in itself, before any face
    /// is looked up: unique names and datum letters, no empty feature, every
    /// referenced datum defined, values in range. Refuses `step.pmi-invalid`
    /// with the item's name.
    pub fn validate(&self) -> Result<(), Refusal> {
        let mut ids = BTreeSet::new();
        let mut letters = BTreeSet::new();
        let mut datum_features = BTreeSet::new();
        for item in &self.items {
            let id = item.id.as_str();
            if id.is_empty() {
                return Err(refuse("an item", "has no name"));
            }
            if !ids.insert(id) {
                return Err(refuse(id, "the name is used twice"));
            }
            for feature in item.kind.features() {
                if feature.is_empty() {
                    return Err(refuse(id, "a feature has no face"));
                }
                let unique: BTreeSet<&F> = feature.iter().collect();
                if unique.len() != feature.len() {
                    return Err(refuse(id, "a feature names a face twice"));
                }
            }
            if let ItemKind::Datum(datum) = &item.kind {
                if !datum_letter(&datum.letter) {
                    return Err(refuse(id, "a datum letter is one or two capital letters"));
                }
                if !letters.insert(datum.letter.as_str()) {
                    return Err(refuse(
                        id,
                        format!("datum {} is defined twice", datum.letter),
                    ));
                }
                let mut key = datum.feature.clone();
                key.sort();
                if !datum_features.insert(key) {
                    return Err(refuse(id, "the feature already carries a datum"));
                }
            }
        }
        for item in &self.items {
            let id = item.id.as_str();
            match &item.kind {
                ItemKind::Datum(_) => {}
                ItemKind::Size(size) => {
                    if size.features.is_empty() {
                        return Err(refuse(id, "the size names no feature"));
                    }
                    if size.through && size.kind != SizeKind::Diameter {
                        return Err(refuse(id, "only a diameter goes through"));
                    }
                    if let Some(thread) = &size.thread {
                        if size.kind != SizeKind::Diameter {
                            return Err(refuse(id, "only a diameter is threaded"));
                        }
                        // A designation such as `M6x1 - 6H`: plain text a STEP string
                        // holds as it stands.
                        let plain = !thread.trim().is_empty()
                            && thread.len() <= 40
                            && thread.bytes().all(|b| b.is_ascii_graphic() || b == b' ')
                            && !thread.contains(['\'', '\\']);
                        if !plain {
                            return Err(refuse(id, "the thread's designation is not plain text"));
                        }
                    }
                    check_limits(id, size.limits.as_ref())?;
                }
                ItemKind::Location(location) => {
                    let mut from = location.from.feature.clone();
                    let mut to = location.to.feature.clone();
                    from.sort();
                    to.sort();
                    if from == to && location.from.at == location.to.at {
                        return Err(refuse(id, "the distance starts and ends on one feature"));
                    }
                    for anchor in [&location.from, &location.to] {
                        if anchor.at == At::Middle && anchor.feature.len() != 2 {
                            return Err(refuse(id, "the middle is taken between two faces"));
                        }
                    }
                    check_limits(id, location.limits.as_ref())?;
                }
                ItemKind::Tolerance(tolerance) => {
                    if tolerance.features.is_empty() {
                        return Err(refuse(id, "the tolerance names no feature"));
                    }
                    if !positive(&tolerance.value) {
                        return Err(refuse(id, "the tolerance value must be more than zero"));
                    }
                    if tolerance.datums.len() > 3 {
                        return Err(refuse(id, "a tolerance refers to at most three datums"));
                    }
                    let distinct: BTreeSet<&str> =
                        tolerance.datums.iter().map(String::as_str).collect();
                    if distinct.len() != tolerance.datums.len() {
                        return Err(refuse(id, "a datum is referred to twice"));
                    }
                    for letter in &tolerance.datums {
                        if !letters.contains(letter.as_str()) {
                            return Err(refuse(id, format!("datum {letter} is not defined")));
                        }
                    }
                    match (tolerance.characteristic.datum_use(), tolerance.datums.len()) {
                        (DatumUse::None, 1..) => {
                            return Err(refuse(id, "a form tolerance refers to no datum"));
                        }
                        (DatumUse::Required, 0) => {
                            return Err(refuse(id, "the tolerance needs a datum"));
                        }
                        _ => {}
                    }
                }
                ItemKind::Texture(texture) => {
                    if texture.parameter.is_empty()
                        || !texture.parameter.bytes().all(|b| b.is_ascii_alphanumeric())
                    {
                        return Err(refuse(id, "the roughness parameter is letters and digits"));
                    }
                    if !positive(&texture.value) {
                        return Err(refuse(id, "the roughness value must be more than zero"));
                    }
                }
            }
        }
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

fn check_limits(id: &str, limits: Option<&Limits>) -> Result<(), Refusal> {
    match limits {
        None | Some(Limits::Basic) => Ok(()),
        Some(Limits::Deviations { upper, lower }) => {
            if upper > lower {
                Ok(())
            } else {
                Err(refuse(
                    id,
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
                Err(refuse(
                    id,
                    "a tolerance class is one or two letters and a grade number",
                ))
            }
        }
    }
}
