//! Manufacturing information read from a STEP file's AP242 semantic PMI.
//!
//! [`read_step_pmi`] returns the datums, dimensions, geometric tolerances and
//! surface textures a file attaches to faces, keyed by the entity names of its
//! `ADVANCED_FACE`s (`#17`), which are the face keys an import of the same
//! bytes reports. It reads the forms `catallax_step::write_step_with_pmi`
//! writes and the common forms of the CAx-IF recommended practices.
//!
//! Nothing is guessed. An item this reader cannot represent in full, such as
//! an angular dimension, a datum target or a tolerance with a modifier it does
//! not know, is left out and listed in [`StepPmi::skipped`] with its entity
//! and the reason.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::Rational;
use catallax_step::pmi::{
    Anchor, At, Characteristic, Datum, Feature, GeometricTolerance, Item, ItemKind, Limits,
    Location, MaterialCondition, PartNotes, Pmi, Process, Size, SizeKind, SurfaceTexture,
};

use crate::ImportRefusal;
use crate::part21::{Entity, Graph, Param, parse};
use crate::recognize::unit_millimetres;

/// The manufacturing information of a STEP file.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct StepPmi {
    /// What was read, with faces named by their entity, such as `#17`.
    pub pmi: Pmi<String>,
    /// Items in the file that were not read.
    pub skipped: Vec<SkippedPmi>,
}

/// An item of manufacturing information that was not read.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SkippedPmi {
    /// The entity, such as `#412`.
    pub entity: String,
    /// Its STEP type.
    pub kind: String,
    /// Why it was left out.
    pub reason: String,
}

const ASPECT_PARTS: [&str; 9] = [
    "SHAPE_ASPECT",
    "DATUM_FEATURE",
    "COMPOSITE_GROUP_SHAPE_ASPECT",
    "COMPOSITE_SHAPE_ASPECT",
    "COMPOSITE_UNIT_SHAPE_ASPECT",
    "ALL_AROUND_SHAPE_ASPECT",
    "CENTRE_OF_SYMMETRY",
    "DERIVED_SHAPE_ASPECT",
    "CONTINUOUS_SHAPE_ASPECT",
];

/// Dimension and datum entities this reader knows of and does not read.
const UNREAD: [(&str, &str); 9] = [
    ("ANGULAR_LOCATION", "an angle between features is not read"),
    ("ANGULAR_SIZE", "an angular size is not read"),
    (
        "DIMENSIONAL_LOCATION_WITH_PATH",
        "a distance along a path is not read",
    ),
    (
        "DIMENSIONAL_SIZE_WITH_PATH",
        "a size along a path is not read",
    ),
    (
        "DIMENSIONAL_LOCATION_WITH_DATUM_FEATURE",
        "a distance that is also a datum feature is not read",
    ),
    (
        "DIMENSIONAL_SIZE_WITH_DATUM_FEATURE",
        "a size that is also a datum feature is not read",
    ),
    ("DATUM_TARGET", "a datum target is not read"),
    ("PLACED_DATUM_TARGET_FEATURE", "a datum target is not read"),
    ("COMMON_DATUM", "a common datum is not read"),
];

fn reference(param: Option<&Param>) -> Option<u64> {
    match param {
        Some(Param::Ref(id)) => Some(*id),
        _ => None,
    }
}

fn references(param: Option<&Param>) -> Vec<u64> {
    match param {
        Some(Param::Ref(id)) => vec![*id],
        Some(Param::List(items)) => items.iter().flat_map(|p| references(Some(p))).collect(),
        Some(Param::Typed(_, items)) => items.iter().flat_map(|p| references(Some(p))).collect(),
        _ => Vec::new(),
    }
}

fn collect_references(params: &[Param], out: &mut Vec<u64>) {
    for param in params {
        match param {
            Param::Ref(id) => out.push(*id),
            Param::List(items) | Param::Typed(_, items) => collect_references(items, out),
            _ => {}
        }
    }
}

/// A STEP string with its control directives decoded.
fn decode(raw: &str) -> String {
    let mut out = String::new();
    let mut rest = raw;
    while let Some(at) = rest.find('\\') {
        out.push_str(&rest[..at]);
        let tail = &rest[at..];
        let hex = |width: usize, body: &str, out: &mut String| {
            let digits: Vec<char> = body.chars().collect();
            for chunk in digits.chunks(width) {
                let text: String = chunk.iter().collect();
                if let Some(c) = u32::from_str_radix(&text, 16).ok().and_then(char::from_u32) {
                    out.push(c);
                }
            }
        };
        if let Some(body) = tail.strip_prefix("\\\\") {
            out.push('\\');
            rest = body;
        } else if let Some(body) = tail.strip_prefix("\\X2\\")
            && let Some(end) = body.find("\\X0\\")
        {
            // UTF-16 code units, so a pair of surrogates is one character.
            let units: Vec<u16> = body.as_bytes()[..end]
                .chunks(4)
                .filter_map(|chunk| u16::from_str_radix(std::str::from_utf8(chunk).ok()?, 16).ok())
                .collect();
            out.push_str(&String::from_utf16_lossy(&units));
            rest = &body[end + 4..];
        } else if let Some(body) = tail.strip_prefix("\\X4\\")
            && let Some(end) = body.find("\\X0\\")
        {
            hex(8, &body[..end], &mut out);
            rest = &body[end + 4..];
        } else if let Some(body) = tail.strip_prefix("\\X\\")
            && body.len() >= 2
            && body.is_char_boundary(2)
        {
            hex(2, &body[..2], &mut out);
            rest = &body[2..];
        } else if let Some(body) = tail.strip_prefix("\\S\\")
            && let Some(c) = body.chars().next()
        {
            // ISO 8859-1 upper half: the character's code plus 128.
            if let Some(c) = char::from_u32(c as u32 + 128) {
                out.push(c);
            }
            rest = &body[c.len_utf8()..];
        } else {
            out.push('\\');
            rest = &tail[1..];
        }
    }
    out.push_str(rest);
    out
}

fn string(param: Option<&Param>) -> Option<String> {
    match param {
        Some(Param::Str(text)) => Some(decode(text)),
        _ => None,
    }
}

fn number(param: Option<&Param>) -> Option<Rational> {
    let text = match param {
        Some(Param::Number(text)) => text,
        Some(Param::Typed(_, items)) => return number(items.first()),
        _ => return None,
    };
    let fixed = text.replace(".E", ".0E").replace(".e", ".0e");
    let fixed = if fixed.ends_with('.') {
        format!("{fixed}0")
    } else {
        fixed
    };
    Rational::parse(&fixed).ok()
}

struct Reader<'a> {
    g: &'a Graph,
    /// Every entity that refers to an entity, by the entity referred to.
    sharing: BTreeMap<u64, Vec<u64>>,
}

/// One item as read, before names are made unique.
struct Found {
    order: u64,
    name: Option<String>,
    kind: ItemKind<String>,
}

impl<'a> Reader<'a> {
    fn new(g: &'a Graph) -> Self {
        let mut sharing: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        let mut found = Vec::new();
        for (id, entity) in g {
            found.clear();
            for (_, params) in &entity.parts {
                collect_references(params, &mut found);
            }
            found.sort_unstable();
            found.dedup();
            for target in &found {
                sharing.entry(*target).or_default().push(*id);
            }
        }
        Self { g, sharing }
    }

    fn entity(&self, id: u64) -> Option<&'a Entity> {
        self.g.get(&id)
    }

    fn has(&self, id: u64, part: &str) -> bool {
        self.entity(id).is_some_and(|e| e.part(part).is_some())
    }

    /// The entities of type `kind` that refer to `id` in parameter `index`.
    fn sharing_at(&self, id: u64, kind: &str, index: usize) -> Vec<(u64, &'a [Param])> {
        let mut out = Vec::new();
        for user in self.sharing.get(&id).map_or(&[][..], Vec::as_slice) {
            if let Some(entity) = self.entity(*user)
                && let Some(params) = entity.part(kind)
                && references(params.get(index)).contains(&id)
            {
                out.push((*user, params));
            }
        }
        out
    }

    /// Name, description, relating and related aspect of a relationship,
    /// whichever subtype or complex form it is written in.
    fn relationship(&self, id: u64) -> Option<(String, u64, u64)> {
        let entity = self.entity(id)?;
        let params = match entity.part("SHAPE_ASPECT_RELATIONSHIP") {
            Some(params) if params.len() >= 4 => params,
            _ if entity.parts.len() == 1 => entity.parts[0].1.as_slice(),
            _ => return None,
        };
        Some((
            string(params.first()).unwrap_or_default(),
            reference(params.get(2))?,
            reference(params.get(3))?,
        ))
    }

    /// The parameters every shape aspect has: name, description, shape, flag.
    fn aspect_params(&self, id: u64) -> Option<&'a [Param]> {
        let entity = self.entity(id)?;
        if let Some(params) = entity.part("SHAPE_ASPECT").filter(|p| p.len() >= 4) {
            return Some(params);
        }
        (entity.parts.len() == 1 && ASPECT_PARTS.contains(&entity.type_name()))
            .then(|| entity.parts[0].1.as_slice())
    }

    fn id_attribute(&self, id: u64) -> Option<String> {
        self.sharing_at(id, "ID_ATTRIBUTE", 1)
            .first()
            .and_then(|(_, params)| string(params.first()))
            .filter(|name| !name.is_empty())
    }

    /// The faces a usage ties straight to an aspect, if it has one.
    fn used_faces(&self, aspect: u64) -> Result<Option<Feature<String>>, String> {
        let mut faces = Vec::new();
        let mut any = false;
        for kind in [
            "GEOMETRIC_ITEM_SPECIFIC_USAGE",
            "ITEM_IDENTIFIED_REPRESENTATION_USAGE",
        ] {
            for (_, params) in self.sharing_at(aspect, kind, 2) {
                any = true;
                for item in references(params.get(4)) {
                    let kind = self.entity(item).map_or("", Entity::type_name);
                    if kind != "ADVANCED_FACE" && kind != "FACE_SURFACE" {
                        return Err(format!("it is attached to {kind}, not to a face"));
                    }
                    let name = format!("#{item}");
                    if !faces.contains(&name) {
                        faces.push(name);
                    }
                }
            }
        }
        Ok(any.then_some(faces).filter(|faces| !faces.is_empty()))
    }

    /// The features an aspect stands for: one, or the members of a pattern.
    fn features(&self, aspect: u64, depth: u32) -> Result<Vec<Feature<String>>, String> {
        if depth > 4 {
            return Err("its features are nested too deeply".into());
        }
        if let Some(faces) = self.used_faces(aspect)? {
            return Ok(vec![faces]);
        }
        let mut members = Vec::new();
        for user in self.sharing.get(&aspect).map_or(&[][..], Vec::as_slice) {
            let Some(entity) = self.entity(*user) else {
                continue;
            };
            if entity.type_name() != "SHAPE_ASPECT_RELATIONSHIP" || entity.parts.len() != 1 {
                continue;
            }
            if let Some((_, relating, related)) = self.relationship(*user)
                && relating == aspect
                && !self.has(related, "DATUM")
            {
                members.push(related);
            }
        }
        if members.is_empty() {
            return Err("it is attached to no face".into());
        }
        let mut features = Vec::new();
        for member in members {
            features.extend(self.features(member, depth + 1)?);
        }
        let params = self.aspect_params(aspect).unwrap_or(&[]);
        let merged = [string(params.first()), string(params.get(1))]
            .into_iter()
            .flatten()
            .any(|text| text == "multiple elements");
        if merged {
            let mut faces: Feature<String> = Vec::new();
            for face in features.into_iter().flatten() {
                if !faces.contains(&face) {
                    faces.push(face);
                }
            }
            return Ok(vec![faces]);
        }
        Ok(features)
    }

    fn one_feature(&self, aspect: u64) -> Result<Feature<String>, String> {
        let mut features = self.features(aspect, 0)?;
        if features.len() != 1 {
            return Err("it applies to a pattern where one feature is expected".into());
        }
        Ok(features.remove(0))
    }

    fn round(&self, face: &str) -> bool {
        face.strip_prefix('#')
            .and_then(|id| id.parse().ok())
            .and_then(|id| self.entity(id))
            .and_then(|face| reference(face.parts[0].1.get(2)))
            .and_then(|surface| self.entity(surface))
            .is_some_and(|surface| surface.type_name() == "CYLINDRICAL_SURFACE")
    }

    /// A length in millimetres from a measure in any of its written forms,
    /// with the name it carries as a representation item.
    fn length(&self, id: u64) -> Option<(String, Rational)> {
        let entity = self.entity(id)?;
        let (name, value, unit) = if let Some(measure) = entity.part("MEASURE_WITH_UNIT") {
            let name = entity
                .part("REPRESENTATION_ITEM")
                .and_then(|params| string(params.first()));
            (name, measure.first(), measure.get(1))
        } else {
            let params = entity.parts[0].1.as_slice();
            match (entity.type_name(), params.len()) {
                ("MEASURE_REPRESENTATION_ITEM", 3) => {
                    (string(params.first()), params.get(1), params.get(2))
                }
                ("LENGTH_MEASURE_WITH_UNIT" | "MEASURE_WITH_UNIT", 2) => {
                    (None, params.first(), params.get(1))
                }
                _ => return None,
            }
        };
        let scale = self.millimetres(reference(unit)?, 0)?;
        Some((name.unwrap_or_default(), &number(value)? * &scale))
    }

    /// Millimetres in one length unit. The importer's own conversion reads a
    /// factor written as `LENGTH_MEASURE_WITH_UNIT`; files also write it as a
    /// plain `MEASURE_WITH_UNIT` holding a length measure.
    fn millimetres(&self, unit: u64, depth: u32) -> Option<Rational> {
        if let Ok(scale) = unit_millimetres(self.g, unit, 0) {
            return Some(scale);
        }
        if depth >= 8 {
            return None;
        }
        let entity = self.entity(unit)?;
        entity.part("LENGTH_UNIT")?;
        let factor = reference(entity.part("CONVERSION_BASED_UNIT")?.get(1))?;
        let factor = self.entity(factor)?;
        if factor.parts.len() != 1 || factor.type_name() != "MEASURE_WITH_UNIT" {
            return None;
        }
        let params = factor.parts[0].1.as_slice();
        if !matches!(params.first(), Some(Param::Typed(kind, _)) if kind.contains("LENGTH_MEASURE"))
        {
            return None;
        }
        let value = number(params.first()).filter(|value| *value > Rational::zero())?;
        Some(&value * &self.millimetres(reference(params.get(1))?, depth + 1)?)
    }

    fn datum_letter(&self, datum: u64) -> Option<String> {
        let entity = self.entity(datum)?;
        if entity.type_name() != "DATUM" || entity.parts.len() != 1 {
            return None;
        }
        string(entity.parts[0].1.get(4))
    }

    fn datum(&self, id: u64) -> Result<Found, String> {
        let letter = self
            .datum_letter(id)
            .ok_or("it has no identification letter")?;
        let mut features = Vec::new();
        let mut name = None;
        for user in self.sharing.get(&id).map_or(&[][..], Vec::as_slice) {
            let Some((_, relating, related)) = self.relationship(*user) else {
                continue;
            };
            if related != id {
                continue;
            }
            if self.has(relating, "DATUM_TARGET")
                || self.has(relating, "PLACED_DATUM_TARGET_FEATURE")
            {
                return Err("it is established by datum targets".into());
            }
            if !self.has(relating, "DATUM_FEATURE") {
                continue;
            }
            name = name.or_else(|| self.id_attribute(relating));
            // A datum feature is one feature however many faces it has. As a
            // group it stands for each member on its own, which only one can.
            let members = self.features(relating, 0)?;
            if self.has(relating, "COMPOSITE_GROUP_SHAPE_ASPECT") && members.len() != 1 {
                return Err("it is on each feature of a pattern".into());
            }
            let mut faces: Feature<String> = Vec::new();
            for face in members.into_iter().flatten() {
                if !faces.contains(&face) {
                    faces.push(face);
                }
            }
            features.push(faces);
        }
        let [feature] = <[Feature<String>; 1]>::try_from(features)
            .map_err(|_| "it does not have exactly one datum feature")?;
        Ok(Found {
            order: id,
            name,
            kind: ItemKind::Datum(Datum { letter, feature }),
        })
    }

    /// The nominal value and limits a dimension states.
    fn dimension_value(&self, id: u64) -> Result<(Option<Rational>, Option<Limits>), String> {
        let (mut nominal, mut upper, mut lower) = (None, None, None);
        for (_, params) in self.sharing_at(id, "DIMENSIONAL_CHARACTERISTIC_REPRESENTATION", 0) {
            let Some(representation) = reference(params.get(1)).and_then(|r| self.entity(r)) else {
                continue;
            };
            for item in references(representation.parts[0].1.get(1)) {
                let Some(entity) = self.entity(item) else {
                    continue;
                };
                if let Some(qualified) = entity.part("QUALIFIED_REPRESENTATION_ITEM") {
                    for qualifier in references(qualified.first()) {
                        if self.has(qualifier, "TYPE_QUALIFIER") {
                            return Err(
                                "its value has a maximum, minimum or average qualifier".into()
                            );
                        }
                    }
                }
                if let Some((name, value)) = self.length(item) {
                    if name.contains("upper") {
                        upper = Some(value);
                    } else if name.contains("lower") {
                        lower = Some(value);
                    } else {
                        nominal = Some(value);
                    }
                    continue;
                }
                match entity.type_name() {
                    "DESCRIPTIVE_REPRESENTATION_ITEM" => {
                        let params = entity.parts[0].1.as_slice();
                        // A basic or reference dimension, or any other note, changes
                        // what the dimension means, and there is nowhere to keep it.
                        let note = string(params.get(1)).unwrap_or_default();
                        return Err(match note.as_str() {
                            "theoretical" => "it is a basic (theoretical) dimension".into(),
                            "auxiliary" => "it is a reference (auxiliary) dimension".into(),
                            _ => format!("it carries the note `{note}`"),
                        });
                    }
                    "COMPOUND_REPRESENTATION_ITEM" => {
                        return Err("it has dimension modifiers".into());
                    }
                    _ => {}
                }
            }
        }
        let mut limits = None;
        for (_, params) in self.sharing_at(id, "PLUS_MINUS_TOLERANCE", 1) {
            let Some(range) = reference(params.first()).and_then(|r| self.entity(r)) else {
                continue;
            };
            let fields = range.parts[0].1.as_slice();
            limits = Some(match range.type_name() {
                "TOLERANCE_VALUE" => {
                    let bound = |index: usize| {
                        reference(fields.get(index))
                            .and_then(|measure| self.length(measure))
                            .map(|(_, value)| value)
                            .ok_or("a deviation is not a length")
                    };
                    Limits::Deviations {
                        lower: bound(0)?,
                        upper: bound(1)?,
                    }
                }
                "LIMITS_AND_FITS" => {
                    // Written as `H` and `7`, or whole as `H7` in either field.
                    let class = format!(
                        "{}{}",
                        string(fields.first()).unwrap_or_default().trim(),
                        string(fields.get(2)).unwrap_or_default().trim()
                    );
                    let split = class
                        .find(|c: char| c.is_ascii_digit())
                        .unwrap_or(class.len());
                    Limits::Fit {
                        deviation: class[..split].to_owned(),
                        grade: class[split..].to_owned(),
                    }
                }
                other => return Err(format!("its limits are a {other}")),
            });
        }
        if limits.is_none()
            && let (Some(nominal), Some(upper), Some(lower)) = (&nominal, &upper, &lower)
        {
            limits = Some(Limits::Deviations {
                upper: upper - nominal,
                lower: lower - nominal,
            });
        }
        Ok((nominal, limits))
    }

    fn size(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let kind = match string(params.get(1)).unwrap_or_default().as_str() {
            "diameter" => SizeKind::Diameter,
            "radius" => SizeKind::Radius,
            other => return Err(format!("a `{other}` size is not read")),
        };
        let aspect = reference(params.first()).ok_or("it names no feature")?;
        let (nominal, limits) = self.dimension_value(id)?;
        Ok(Found {
            order: id,
            name: self.id_attribute(id),
            kind: ItemKind::Size(Size {
                kind,
                features: self.features(aspect, 0)?,
                limits,
                nominal,
            }),
        })
    }

    /// One end of a distance: its faces, and whether it is their middle.
    fn end(&self, aspect: u64) -> Result<(Feature<String>, bool), String> {
        if self.has(aspect, "CENTRE_OF_SYMMETRY") {
            let mut faces = Vec::new();
            for user in self.sharing.get(&aspect).map_or(&[][..], Vec::as_slice) {
                if self.entity(*user).map(Entity::type_name)
                    != Some("SHAPE_ASPECT_DERIVING_RELATIONSHIP")
                {
                    continue;
                }
                if let Some((_, relating, related)) = self.relationship(*user)
                    && relating == aspect
                {
                    faces.extend(self.one_feature(related)?);
                }
            }
            if faces.len() != 2 {
                return Err("its centre is not taken between two faces".into());
            }
            return Ok((faces, true));
        }
        for derived in [
            "DERIVED_SHAPE_ASPECT",
            "APEX",
            "GEOMETRIC_ALIGNMENT",
            "PERPENDICULAR_TO",
            "EXTENSION",
            "TANGENT",
            "PARALLEL_OFFSET",
        ] {
            if self.has(aspect, derived) {
                return Err(format!("it ends on a derived shape ({derived})"));
            }
        }
        Ok((self.one_feature(aspect)?, false))
    }

    fn location(&self, id: u64) -> Result<Found, String> {
        let (name, relating, related) = self.relationship(id).ok_or("it names no two features")?;
        let sides = name
            .strip_prefix("linear distance")
            .ok_or_else(|| format!("a `{name}` is not read"))?;
        let sides: Vec<&str> = sides.split_whitespace().collect();
        let (from_faces, from_middle) = self.end(relating)?;
        let (to_faces, to_middle) = self.end(related)?;
        let from_round = from_faces.first().is_some_and(|face| self.round(face));
        let to_round = to_faces.first().is_some_and(|face| self.round(face));
        let (from_side, to_side) = match (sides.as_slice(), from_round, to_round) {
            ([], _, _) => ("centre", "centre"),
            ([a, b], _, _) => (*a, *b),
            ([a], true, false) => (*a, "centre"),
            ([a], false, true) => ("centre", *a),
            _ => {
                return Err(format!(
                    "the distance name `{name}` does not fit its features"
                ));
            }
        };
        let at = |side: &str, middle: bool| match (side, middle) {
            ("centre", true) => Ok(At::Middle),
            ("centre", false) => Ok(At::Feature),
            ("inner", false) => Ok(At::Inner),
            ("outer", false) => Ok(At::Outer),
            _ => Err(format!(
                "the distance name `{name}` does not fit its features"
            )),
        };
        let (nominal, limits) = self.dimension_value(id)?;
        Ok(Found {
            order: id,
            name: self.id_attribute(id),
            kind: ItemKind::Location(Location {
                from: Anchor {
                    feature: from_faces,
                    at: at(from_side, from_middle)?,
                },
                to: Anchor {
                    feature: to_faces,
                    at: at(to_side, to_middle)?,
                },
                limits,
                nominal,
            }),
        })
    }

    /// The datum letters a tolerance refers to, primary first.
    fn datum_letters(&self, systems: &[u64]) -> Result<Vec<String>, String> {
        let mut letters = Vec::new();
        let mut ranked = Vec::new();
        for system in systems {
            let entity = self.entity(*system).ok_or("a datum reference is missing")?;
            let params = entity.parts[0].1.as_slice();
            match entity.type_name() {
                "DATUM_SYSTEM" => {
                    for compartment in references(params.get(4)) {
                        let fields = self
                            .entity(compartment)
                            .filter(|e| e.type_name() == "DATUM_REFERENCE_COMPARTMENT")
                            .map(|e| e.parts[0].1.as_slice())
                            .ok_or("a datum reference is not a single compartment")?;
                        if !matches!(fields.get(5), Some(Param::Unset) | None) {
                            return Err("a datum reference has a modifier".into());
                        }
                        letters.push(
                            reference(fields.get(4))
                                .and_then(|datum| self.datum_letter(datum))
                                .ok_or("it refers to a common datum")?,
                        );
                    }
                }
                "DATUM_REFERENCE" => {
                    let rank = number(params.first()).ok_or("a datum has no precedence")?;
                    let letter = reference(params.get(1))
                        .and_then(|datum| self.datum_letter(datum))
                        .ok_or("a datum reference names no datum")?;
                    ranked.push((rank, letter));
                }
                other => return Err(format!("it refers to datums through a {other}")),
            }
        }
        ranked.sort();
        letters.extend(ranked.into_iter().map(|(_, letter)| letter));
        Ok(letters)
    }

    fn tolerance(&self, id: u64, characteristic: Characteristic) -> Result<Found, String> {
        let entity = self.entity(id).ok_or("it is missing")?;
        let own = characteristic.entity();
        for (part, _) in &entity.parts {
            let known = [
                own,
                "GEOMETRIC_TOLERANCE",
                "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE",
                "GEOMETRIC_TOLERANCE_WITH_MODIFIERS",
                "MODIFIED_GEOMETRIC_TOLERANCE",
            ];
            if !known.contains(&part.as_str()) {
                return Err(format!("it is also a {part}"));
            }
        }
        let simple = entity.parts.len() == 1;
        let base = match entity.part("GEOMETRIC_TOLERANCE") {
            Some(params) => params,
            None if simple => entity.parts[0].1.as_slice(),
            None => return Err("it has no tolerance values".into()),
        };
        let systems = match entity.part("GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE") {
            Some(params) => references(params.first()),
            None if simple => references(base.get(4)),
            None => Vec::new(),
        };
        let mut material = None;
        for part in [
            "GEOMETRIC_TOLERANCE_WITH_MODIFIERS",
            "MODIFIED_GEOMETRIC_TOLERANCE",
        ] {
            let Some(params) = entity.part(part) else {
                continue;
            };
            let modifiers: Vec<&Param> = match params.first() {
                Some(Param::List(items)) => items.iter().collect(),
                Some(other) => vec![other],
                None => Vec::new(),
            };
            for modifier in modifiers {
                let Param::Enum(name) = modifier else {
                    return Err("it has an unreadable modifier".into());
                };
                material = Some(match name.as_str() {
                    "MAXIMUM_MATERIAL_REQUIREMENT" | "MAXIMUM_MATERIAL_CONDITION" => {
                        MaterialCondition::Maximum
                    }
                    "LEAST_MATERIAL_REQUIREMENT" | "LEAST_MATERIAL_CONDITION" => {
                        MaterialCondition::Least
                    }
                    "REGARDLESS_OF_FEATURE_SIZE" => continue,
                    other => return Err(format!("it has the modifier {other}")),
                });
            }
        }
        let value = reference(base.get(2))
            .and_then(|measure| self.length(measure))
            .map(|(_, value)| value)
            .ok_or("it has no length value")?;
        let target = reference(base.get(3)).ok_or("it names no feature")?;
        let target = match self.entity(target) {
            Some(entity) if entity.type_name() == "DIMENSIONAL_SIZE" => {
                reference(entity.parts[0].1.first()).ok_or("it names no feature")?
            }
            Some(entity) if entity.type_name() == "PRODUCT_DEFINITION_SHAPE" => {
                return Err("it applies to the whole part".into());
            }
            Some(_) if self.aspect_params(target).is_some() => target,
            _ => return Err("it is not attached to a feature".into()),
        };
        let mut diameter_zone = false;
        for (_, params) in self.sharing_at(id, "TOLERANCE_ZONE", 4) {
            let form = reference(params.get(5))
                .and_then(|form| self.entity(form))
                .and_then(|form| string(form.parts[0].1.first()))
                .unwrap_or_default();
            match form.as_str() {
                "cylindrical or circular" | "within a cylinder" | "within a circle" => {
                    diameter_zone = true;
                }
                "spherical" | "within a sphere" => {
                    return Err("its zone is a sphere".into());
                }
                _ => {}
            }
        }
        let name = string(base.first())
            .filter(|name| !name.is_empty())
            .or_else(|| self.id_attribute(id));
        Ok(Found {
            order: id,
            name,
            kind: ItemKind::Tolerance(GeometricTolerance {
                characteristic,
                value,
                diameter_zone,
                material,
                datums: self.datum_letters(&systems)?,
                features: self.features(target, 0)?,
            }),
        })
    }

    fn texture(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let mut measures = Vec::new();
        for item in references(params.get(1)) {
            match self.entity(item).map(Entity::type_name) {
                Some("VALUE_RANGE") => return Err("its value is a range".into()),
                Some("DESCRIPTIVE_REPRESENTATION_ITEM") => {}
                _ => measures.push(self.length(item).ok_or("its value is not a length")?),
            }
        }
        let [(parameter, millimetres)] = <[(String, Rational); 1]>::try_from(measures)
            .map_err(|_| "it does not have exactly one value")?;
        let uses = self.sharing_at(id, "PROPERTY_DEFINITION_REPRESENTATION", 1);
        let property = uses
            .first()
            .and_then(|(_, params)| reference(params.first()))
            .ok_or("it belongs to no property")?;
        let fields = self
            .entity(property)
            .filter(|e| e.type_name() == "PROPERTY_DEFINITION")
            .map(|e| e.parts[0].1.as_slice())
            .ok_or("it belongs to no property")?;
        let aspect = reference(fields.get(2))
            .filter(|aspect| self.aspect_params(*aspect).is_some())
            .ok_or("it is not attached to a feature")?;
        let process = string(fields.get(1))
            .and_then(|phrase| Process::of_phrase(&phrase))
            .unwrap_or_default();
        Ok(Found {
            order: property,
            name: self.id_attribute(property),
            kind: ItemKind::Texture(SurfaceTexture {
                parameter,
                value: &millimetres * &Rational::from_i64(1000),
                process,
                feature: self.one_feature(aspect)?,
            }),
        })
    }

    /// The notes written as user defined attributes of the product.
    fn part_notes(&self) -> PartNotes {
        let mut notes = PartNotes::default();
        let mut numbered = Vec::new();
        for (id, entity) in self.g {
            if entity.type_name() != "PROPERTY_DEFINITION" || entity.parts.len() != 1 {
                continue;
            }
            let params = entity.parts[0].1.as_slice();
            if string(params.get(1)).as_deref() != Some("user defined attribute") {
                continue;
            }
            let on_product = reference(params.get(2)).is_some_and(|target| {
                self.has(target, "PRODUCT_DEFINITION")
                    || self.has(target, "PRODUCT_DEFINITION_SHAPE")
            });
            let Some(name) = string(params.first()).filter(|_| on_product) else {
                continue;
            };
            let mut text = None;
            for (_, uses) in self.sharing_at(*id, "PROPERTY_DEFINITION_REPRESENTATION", 0) {
                let Some(representation) = reference(uses.get(1)).and_then(|r| self.entity(r))
                else {
                    continue;
                };
                for item in references(representation.parts[0].1.get(1)) {
                    if let Some(item) = self
                        .entity(item)
                        .filter(|e| e.type_name() == "DESCRIPTIVE_REPRESENTATION_ITEM")
                    {
                        text = text.or_else(|| string(item.parts[0].1.get(1)));
                    }
                }
            }
            let Some(text) = text.filter(|text| !text.trim().is_empty()) else {
                continue;
            };
            if name == "general tolerance" {
                notes.general_tolerance.get_or_insert(text);
            } else if name == "material" {
                notes.material.get_or_insert(text);
            } else if let Some(rest) = name.strip_prefix("note") {
                let rank = rest.trim().parse::<u64>().unwrap_or(u64::MAX);
                numbered.push((rank, *id, text));
            }
        }
        numbered.sort();
        notes.notes = numbered.into_iter().map(|(_, _, text)| text).collect();
        notes
    }
}

/// Reads the manufacturing information of a STEP file.
///
/// Faces are named by entity, such as `#17`: the keys an import of the same
/// bytes gives its faces. A dimension's [`Size::nominal`] or
/// [`Location::nominal`] is the value the file states, in millimetres; the
/// caller decides what to do when the geometry disagrees.
///
/// Refuses only a file that is not readable Part 21 text. What cannot be
/// represented is listed in [`StepPmi::skipped`].
pub fn read_step_pmi(bytes: &[u8]) -> Result<StepPmi, ImportRefusal> {
    let graph = parse(bytes)?;
    let reader = Reader::new(&graph);
    let mut found = Vec::new();
    let mut skipped = Vec::new();
    let mut skip = |id: u64, kind: &str, reason: String| {
        skipped.push(SkippedPmi {
            entity: format!("#{id}"),
            kind: kind.to_owned(),
            reason,
        });
    };
    for (id, entity) in &graph {
        let id = *id;
        let kind = entity.type_name();
        if let Some((part, reason)) = UNREAD.iter().find(|(part, _)| entity.part(part).is_some()) {
            skip(id, part, (*reason).to_owned());
            continue;
        }
        let read = if kind == "DATUM" && entity.parts.len() == 1 {
            reader.datum(id)
        } else if kind == "DIMENSIONAL_SIZE" && entity.parts.len() == 1 {
            reader.size(id, entity.parts[0].1.as_slice())
        } else if entity.part("DIMENSIONAL_LOCATION").is_some() {
            reader.location(id)
        } else if kind == "SURFACE_TEXTURE_REPRESENTATION" && entity.parts.len() == 1 {
            reader.texture(id, entity.parts[0].1.as_slice())
        } else if let Some(characteristic) = entity
            .parts
            .iter()
            .find_map(|(part, _)| Characteristic::of_entity(part))
        {
            reader.tolerance(id, characteristic)
        } else {
            continue;
        };
        match read {
            Ok(item) => found.push((kind.to_owned(), item)),
            Err(reason) => skip(id, kind, reason),
        }
    }
    found.sort_by_key(|(_, item)| item.order);

    // A tolerance needs its datums, so datums are settled first.
    let mut pmi = Pmi::<String> {
        items: Vec::new(),
        part: reader.part_notes(),
    };
    let mut names = BTreeSet::new();
    let mut letters = BTreeSet::new();
    let mut datum_features = BTreeSet::new();
    let mut later = Vec::new();
    let name_of = |found: &Found, names: &mut BTreeSet<String>| {
        let base = found.name.clone().unwrap_or_else(|| match &found.kind {
            ItemKind::Datum(datum) => format!("Datum {}", datum.letter),
            _ => format!("#{}", found.order),
        });
        let mut name = base.clone();
        let mut copy = 2;
        while !names.insert(name.clone()) {
            name = format!("{base} ({copy})");
            copy += 1;
        }
        name
    };
    for (kind, item) in found {
        let ItemKind::Datum(datum) = &item.kind else {
            later.push((kind, item));
            continue;
        };
        let mut feature = datum.feature.clone();
        feature.sort();
        if !letters.insert(datum.letter.clone()) {
            skip(
                item.order,
                &kind,
                format!("datum {} is defined twice", datum.letter),
            );
            continue;
        }
        if !datum_features.insert(feature) {
            skip(
                item.order,
                &kind,
                "its feature already carries a datum".into(),
            );
            continue;
        }
        let id = name_of(&item, &mut names);
        pmi.items.push(Item {
            id,
            kind: item.kind,
        });
    }
    for (kind, item) in later {
        if let ItemKind::Tolerance(tolerance) = &item.kind
            && let Some(letter) = tolerance.datums.iter().find(|l| !letters.contains(*l))
        {
            skip(
                item.order,
                &kind,
                format!("it refers to datum {letter}, which was not read"),
            );
            continue;
        }
        let order = item.order;
        let candidate = Item {
            id: "item".to_owned(),
            kind: item.kind.clone(),
        };
        // Each item is checked on its own, with the datums already settled.
        let mut alone = Pmi::<String> {
            items: pmi
                .items
                .iter()
                .filter(|existing| matches!(existing.kind, ItemKind::Datum(_)))
                .cloned()
                .collect(),
            part: PartNotes::default(),
        };
        alone.items.push(candidate);
        if let Err(refusal) = alone.validate() {
            let reason = refusal
                .detail
                .strip_prefix("item: ")
                .unwrap_or(&refusal.detail)
                .to_owned();
            skip(order, &kind, reason);
            continue;
        }
        let id = name_of(&item, &mut names);
        pmi.items.push(Item {
            id,
            kind: item.kind,
        });
    }
    Ok(StepPmi { pmi, skipped })
}

#[cfg(test)]
mod tests {
    use super::decode;

    #[test]
    fn strings_decode_their_control_directives() {
        assert_eq!(decode("plain"), "plain");
        assert_eq!(decode("a\\\\b"), "a\\b");
        assert_eq!(decode("\\X2\\00D8\\X0\\8 H7"), "\u{d8}8 H7");
        assert_eq!(decode("\\X\\D8"), "\u{d8}");
        assert_eq!(decode("\\S\\X"), "\u{d8}");
        assert_eq!(decode("\\X4\\0001F600\\X0\\"), "\u{1f600}");
    }
}
