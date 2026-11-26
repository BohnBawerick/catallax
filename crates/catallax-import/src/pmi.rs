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
use catallax_shapes::{RPoint3, RVector3};
use catallax_step::pmi::{
    Anchor, At, COUNTERBORE, COUNTERSINK, Characteristic, Datum, Feature, GeometricTolerance,
    HEAD_DEPTH, HEAD_DIAMETER, HOLE_DEPTH, HOLE_PILOT, HoleHead, HoleSize, INCLUDED_ANGLE, Item,
    ItemError, ItemKind, Limits, Location, MaterialCondition, PartNotes, PartValidator, Pmi,
    Process, Size, SizeKind, SurfaceTexture, THREAD, THREAD_DEPTH, THROUGH,
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

pub(crate) fn string(param: Option<&Param>) -> Option<String> {
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

fn coordinates(param: Option<&Param>) -> Option<[Rational; 3]> {
    let Some(Param::List(values)) = param else {
        return None;
    };
    let [x, y, z] = values.as_slice() else {
        return None;
    };
    Some([number(Some(x))?, number(Some(y))?, number(Some(z))?])
}

struct RoundSurface {
    placement: u64,
    radius: Rational,
    angle: Option<Rational>,
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

    fn owner_shapes(&self, id: u64) -> BTreeSet<u64> {
        let mut pending = vec![id];
        // A surface texture representation is owned through its property.
        pending.extend(
            self.sharing_at(id, "PROPERTY_DEFINITION_REPRESENTATION", 1)
                .into_iter()
                .filter_map(|(_, params)| reference(params.first())),
        );
        let mut seen = BTreeSet::new();
        let mut shapes = BTreeSet::new();
        while let Some(id) = pending.pop() {
            let Some(entity) = self.entity(id).filter(|_| seen.insert(id)) else {
                continue;
            };
            if entity.type_name() == "PRODUCT_DEFINITION_SHAPE" {
                shapes.insert(id);
                continue;
            }
            // Ownership is through the aspect or property, never geometry.
            if matches!(entity.type_name(), "PRODUCT_DEFINITION" | "ADVANCED_FACE")
                || entity.type_name().ends_with("SHAPE_REPRESENTATION")
            {
                continue;
            }
            for (_, params) in &entity.parts {
                collect_references(params, &mut pending);
            }
        }
        shapes
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

    fn face_surface(&self, face: &str) -> Option<u64> {
        face.strip_prefix('#')
            .and_then(|id| id.parse().ok())
            .and_then(|id| self.entity(id))
            .and_then(|face| reference(face.parts[0].1.get(2)))
    }

    fn round(&self, face: &str) -> bool {
        self.face_surface(face)
            .and_then(|surface| self.entity(surface))
            .is_some_and(|surface| surface.type_name() == "CYLINDRICAL_SURFACE")
    }

    fn coaxial(&self, a: u64, b: u64) -> Option<bool> {
        if a == b {
            return Some(true);
        }
        let axis = |id| {
            let placement = self.entity(id)?.part("AXIS2_PLACEMENT_3D")?;
            let point = self
                .entity(reference(placement.get(1))?)?
                .part("CARTESIAN_POINT")?;
            let [x, y, z] = coordinates(point.get(1))?;
            let point = RPoint3::new(x, y, z);
            let direction = match placement.get(2) {
                Some(Param::Unset) => {
                    RVector3::new(Rational::zero(), Rational::zero(), Rational::one())
                }
                param => {
                    let dir = self.entity(reference(param)?)?.part("DIRECTION")?;
                    let [x, y, z] = coordinates(dir.get(1))?;
                    RVector3::new(x, y, z)
                }
            };
            (!direction.is_zero()).then_some((point, direction))
        };
        let (point, direction) = axis(a)?;
        let (other_point, other_direction) = axis(b)?;
        Some(
            direction.cross(&other_direction).is_zero()
                && other_point.minus(&point).cross(&direction).is_zero(),
        )
    }

    fn round_feature(&self, faces: &[String], cone: bool) -> Result<RoundSurface, String> {
        let mut surfaces = faces.iter().map(|face| {
            self.face_surface(face)
                .ok_or("a size face has no readable surface")
        });
        let first = surfaces.next().ok_or("the size has no face")??;
        let read = |id| {
            let entity = self.entity(id)?;
            let params = entity.part(if cone {
                "CONICAL_SURFACE"
            } else {
                "CYLINDRICAL_SURFACE"
            })?;
            Some(RoundSurface {
                placement: reference(params.get(1))?,
                radius: number(params.get(2))?,
                angle: if cone {
                    Some(number(params.get(3))?)
                } else {
                    None
                },
            })
        };
        let bad_surface = if cone {
            "the head is not conical"
        } else {
            "a size face is not cylindrical"
        };
        let shape = read(first).ok_or(bad_surface)?;
        for surface in surfaces {
            let surface = surface?;
            if surface == first {
                continue;
            }
            let other = read(surface).ok_or(bad_surface)?;
            if shape.radius != other.radius
                || shape.angle != other.angle
                || (if cone {
                    shape.placement != other.placement
                } else {
                    self.coaxial(shape.placement, other.placement) != Some(true)
                })
            {
                return Err("the size groups non-equivalent faces".into());
            }
        }
        Ok(shape)
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
        let mut basic = false;
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
                        // A reference dimension, or any other note, changes what the
                        // dimension means, and there is nowhere to keep it.
                        let note = string(params.get(1)).unwrap_or_default();
                        if note == "theoretical" {
                            basic = true;
                            continue;
                        }
                        return Err(match note.as_str() {
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
        if basic {
            if limits.is_some() {
                return Err("it is basic (theoretical) and has limits".into());
            }
            limits = Some(Limits::Basic);
        }
        Ok((nominal, limits))
    }

    fn hole_attribute(&self, aspect: u64, name: &str) -> Result<Option<u64>, String> {
        let properties = self.sharing_at(aspect, "PROPERTY_DEFINITION", 2);
        let mut matching = properties.into_iter().filter(|(_, params)| {
            matches!(params.first(), Some(Param::Str(value)) if value == name)
                && matches!(params.get(1), Some(Param::Str(value)) if value == "user defined attribute")
        });
        let Some((id, _)) = matching.next() else {
            return Ok(None);
        };
        let links = self.sharing_at(id, "PROPERTY_DEFINITION_REPRESENTATION", 0);
        if matching.next().is_none()
            && let [(_, link)] = links.as_slice()
            && let Some(rep) = reference(link.get(1)).and_then(|id| self.entity(id))
            && let [item] = references(rep.parts[0].1.get(1)).as_slice()
        {
            return Ok(Some(*item));
        }
        Err(format!("the `{name}` attribute has no unique value"))
    }

    fn hole_text(&self, item: u64) -> Option<String> {
        let params = self.entity(item)?.part("DESCRIPTIVE_REPRESENTATION_ITEM")?;
        string(params.get(1))
    }

    fn hole_length(&self, item: u64, name: &str) -> Result<Rational, String> {
        self.length(item)
            .map(|(_, value)| value)
            .ok_or_else(|| format!("the `{name}` attribute is not a length with a known unit"))
    }

    fn hole_angle(&self, item: u64) -> Option<Rational> {
        let measure = self.entity(item)?;
        measure.part("PLANE_ANGLE_MEASURE_WITH_UNIT")?;
        let params = measure.part("MEASURE_WITH_UNIT")?;
        let unit = self.entity(reference(params.get(1))?)?;
        unit.part("PLANE_ANGLE_UNIT")?;
        let conversion = unit.part("CONVERSION_BASED_UNIT")?;
        if !matches!(conversion.first(), Some(Param::Str(name)) if name.eq_ignore_ascii_case("degree"))
        {
            return None;
        }
        let factor = self
            .entity(reference(conversion.get(1))?)?
            .part("PLANE_ANGLE_MEASURE_WITH_UNIT")?;
        let radian = self.entity(reference(factor.get(1))?)?;
        radian.part("PLANE_ANGLE_UNIT")?;
        let si = radian.part("SI_UNIT")?;
        if !matches!(si, [Param::Unset, Param::Enum(name)] if name == "RADIAN")
            || number(factor.first())? != Rational::parse("0.017453292519943295").ok()?
        {
            return None;
        }
        match params.first() {
            Some(Param::Typed(kind, values)) if kind == "PLANE_ANGLE_MEASURE" => {
                number(values.first())
            }
            _ => None,
        }
    }

    fn hole_head(
        &self,
        pilot: u64,
        pilot_radius: &Rational,
        head: Option<(&str, u64)>,
    ) -> Result<Option<(HoleHead, Feature<String>)>, String> {
        let diameter = self.hole_attribute(pilot, HEAD_DIAMETER)?;
        let depth = self.hole_attribute(pilot, HEAD_DEPTH)?;
        let angle = self.hole_attribute(pilot, INCLUDED_ANGLE)?;
        let head = match (head, diameter, depth, angle) {
            (Some((COUNTERBORE, aspect)), Some(diameter), Some(depth), None) => (
                HoleHead::Counterbore {
                    diameter: self.hole_length(diameter, HEAD_DIAMETER)?,
                    depth: self.hole_length(depth, HEAD_DEPTH)?,
                },
                aspect,
            ),
            (Some((COUNTERSINK, aspect)), Some(diameter), None, Some(angle)) => (
                HoleHead::Countersink {
                    diameter: self.hole_length(diameter, HEAD_DIAMETER)?,
                    included_angle: self
                        .hole_angle(angle)
                        .ok_or("the included angle is not in the supported degree unit")?,
                },
                aspect,
            ),
            (None, None, None, None) => return Ok(None),
            _ => return Err("the hole has an unknown or incomplete head".into()),
        };
        if self
            .aspect_params(head.1)
            .and_then(|params| reference(params.get(2)))
            != self
                .aspect_params(pilot)
                .and_then(|params| reference(params.get(2)))
        {
            return Err("the hole's head belongs to another part".into());
        }
        let faces = self
            .used_faces(head.1)?
            .ok_or("the hole head names no face")?;
        let cone = matches!(&head.0, HoleHead::Countersink { .. });
        let shape = self.round_feature(&faces, cone)?;
        // The explicit relationship identifies the head. In a STEP export
        // the two placements round independently, so exact coaxiality is
        // checked against the model on export, not these written coordinates.
        if !cone && shape.radius <= *pilot_radius {
            return Err("the counterbore head is not wider than the pilot".into());
        }
        Ok(Some((head.0, faces)))
    }

    /// The text of a user defined attribute.
    fn attribute_text(&self, property: u64) -> Option<String> {
        let mut text = None;
        for (_, uses) in self.sharing_at(property, "PROPERTY_DEFINITION_REPRESENTATION", 0) {
            let Some(representation) = reference(uses.get(1)).and_then(|r| self.entity(r)) else {
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
        text
    }

    fn size(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let mut kind = match string(params.get(1)).unwrap_or_default().as_str() {
            "diameter" => SizeKind::Diameter,
            "radius" => SizeKind::Radius,
            other => return Err(format!("a `{other}` size is not read")),
        };
        let aspect = reference(params.first()).ok_or("it names no feature")?;
        let (nominal, limits) = self.dimension_value(id)?;
        let is_hole = self
            .aspect_params(aspect)
            .and_then(|params| string(params.first()))
            .as_deref()
            == Some(HOLE_PILOT);
        let relationships = self.sharing_at(aspect, "SHAPE_ASPECT_RELATIONSHIP", 2);
        let mut heads = relationships
            .iter()
            .filter_map(|(_, params)| match params.first() {
                Some(Param::Str(kind)) if matches!(kind.as_str(), COUNTERBORE | COUNTERSINK) => {
                    Some((kind.as_str(), reference(params.get(3))))
                }
                _ => None,
            });
        let depth = self.hole_attribute(aspect, HOLE_DEPTH)?;
        let through = depth.is_some_and(|value| self.hole_text(value).as_deref() == Some(THROUGH));
        let drill = depth.filter(|_| !through);
        let thread = self
            .hole_attribute(aspect, THREAD)?
            .map(|value| {
                self.hole_text(value)
                    .ok_or("the thread designation is not text")
            })
            .transpose()?;
        let thread_depth = self.hole_attribute(aspect, THREAD_DEPTH)?;
        let mut features;
        if is_hole {
            if kind != SizeKind::Diameter || nominal.is_none() {
                return Err("the hole has no readable pilot diameter".into());
            }
            let head = match (heads.next(), relationships.len()) {
                (None, 0) => None,
                (Some((kind, Some(id))), 1) => Some((kind, id)),
                _ => return Err("the hole has an unknown or ambiguous head relationship".into()),
            };
            let pilot = self
                .used_faces(aspect)?
                .ok_or("the hole pilot names no face")?;
            let shape = self.round_feature(&pilot, false)?;
            let head = self.hole_head(aspect, &shape.radius, head)?;
            if let Some((head, _)) = &head {
                let diameter = match head {
                    HoleHead::Counterbore { diameter, .. }
                    | HoleHead::Countersink { diameter, .. } => diameter,
                };
                if nominal.as_ref().is_some_and(|pilot| diameter <= pilot) {
                    return Err("the head diameter must be larger than the pilot".into());
                }
            }
            features = vec![pilot];
            let head = head.map(|(head, faces)| {
                features.push(faces);
                head
            });
            kind = SizeKind::Hole(HoleSize {
                head,
                drill_depth: drill
                    .map(|id| self.hole_length(id, HOLE_DEPTH))
                    .transpose()?,
                thread_depth: thread_depth
                    .map(|id| self.hole_length(id, THREAD_DEPTH))
                    .transpose()?,
            });
        } else {
            if heads.next().is_some() || drill.is_some() || thread_depth.is_some() {
                return Err(
                    "the size carries hole information without a supported pilot aspect".into(),
                );
            }
            features = self.features(aspect, 0)?;
            let mut radius = None;
            for feature in &features {
                let next = self.round_feature(feature, false)?.radius;
                if radius.as_ref().is_some_and(|radius| *radius != next) {
                    return Err("the size groups features of different diameters".into());
                }
                radius = Some(next);
            }
        }
        Ok(Found {
            order: id,
            name: self.id_attribute(id),
            kind: ItemKind::Size(Size {
                kind,
                features,
                limits,
                through,
                thread,
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
    fn part_notes(&self, scope: Option<(u64, u64)>) -> PartNotes {
        let mut notes = PartNotes::default();
        let mut numbered = Vec::new();
        let candidates: Box<dyn Iterator<Item = (&u64, &Entity)>> = match scope {
            Some((definition, shape)) => Box::new(
                [definition, shape]
                    .into_iter()
                    .flat_map(|target| self.sharing.get(&target).into_iter().flatten())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .filter_map(|id| self.entity(*id).map(|entity| (id, entity))),
            ),
            None => Box::new(self.g.iter()),
        };
        for (id, entity) in candidates {
            if entity.type_name() != "PROPERTY_DEFINITION" || entity.parts.len() != 1 {
                continue;
            }
            let params = entity.parts[0].1.as_slice();
            if string(params.get(1)).as_deref() != Some("user defined attribute") {
                continue;
            }
            let on_product = reference(params.get(2)).is_some_and(|target| match scope {
                Some((definition, shape)) => target == definition || target == shape,
                None => {
                    self.has(target, "PRODUCT_DEFINITION")
                        || self.has(target, "PRODUCT_DEFINITION_SHAPE")
                }
            });
            let Some(name) = string(params.first()).filter(|_| on_product) else {
                continue;
            };
            let text = self.attribute_text(*id);
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
    Ok(read_items(
        &reader,
        graph
            .iter()
            .filter_map(|(&id, entity)| candidate(id, entity)),
        None,
    ))
}

pub(crate) fn read_assembly_pmi(
    graph: &Graph,
    shapes: &BTreeMap<u64, u64>,
) -> Result<BTreeMap<u64, StepPmi>, ImportRefusal> {
    let reader = Reader::new(graph);
    let mut items = BTreeMap::<u64, Vec<Candidate<'_>>>::new();
    for item in graph
        .iter()
        .filter_map(|(&id, entity)| candidate(id, entity))
    {
        let owners = reader.owner_shapes(item.id);
        let owner = owners
            .first()
            .filter(|id| owners.len() == 1 && shapes.contains_key(id));
        let Some(&owner) = owner else {
            return Err(ImportRefusal::Unsupported {
                entity: format!("#{}", item.id),
                type_name: "PMI without a unique part attachment".into(),
            });
        };
        items.entry(owner).or_default().push(item);
    }
    Ok(shapes
        .iter()
        .map(|(&shape, &definition)| {
            let pmi = read_items(
                &reader,
                items.remove(&shape).unwrap_or_default(),
                Some((definition, shape)),
            );
            (shape, pmi)
        })
        .collect())
}

enum PmiKind {
    Unread(&'static str, &'static str),
    Datum,
    Size,
    Location,
    Texture,
    Tolerance(Characteristic),
}

struct Candidate<'a> {
    id: u64,
    entity: &'a Entity,
    kind: PmiKind,
}

fn candidate(id: u64, entity: &Entity) -> Option<Candidate<'_>> {
    let name = entity.type_name();
    let kind = if let Some(&(part, reason)) =
        UNREAD.iter().find(|(part, _)| entity.part(part).is_some())
    {
        PmiKind::Unread(part, reason)
    } else if name == "DATUM" && entity.parts.len() == 1 {
        PmiKind::Datum
    } else if name == "DIMENSIONAL_SIZE" && entity.parts.len() == 1 {
        PmiKind::Size
    } else if entity.part("DIMENSIONAL_LOCATION").is_some() {
        PmiKind::Location
    } else if name == "SURFACE_TEXTURE_REPRESENTATION" && entity.parts.len() == 1 {
        PmiKind::Texture
    } else {
        PmiKind::Tolerance(
            entity
                .parts
                .iter()
                .find_map(|(part, _)| Characteristic::of_entity(part))?,
        )
    };
    Some(Candidate { id, entity, kind })
}

fn read_items<'a>(
    reader: &Reader<'a>,
    candidates: impl IntoIterator<Item = Candidate<'a>>,
    scope: Option<(u64, u64)>,
) -> StepPmi {
    let mut datums = Vec::new();
    let mut later = Vec::new();
    let mut skipped = Vec::new();
    let mut skip = |id: u64, kind: &str, reason: String| {
        skipped.push(SkippedPmi {
            entity: format!("#{id}"),
            kind: kind.to_owned(),
            reason,
        });
    };
    for Candidate { id, entity, kind } in candidates {
        let read = match kind {
            PmiKind::Unread(part, reason) => {
                skip(id, part, reason.to_owned());
                continue;
            }
            PmiKind::Datum => reader.datum(id),
            PmiKind::Size => reader.size(id, &entity.parts[0].1),
            PmiKind::Location => reader.location(id),
            PmiKind::Texture => reader.texture(id, &entity.parts[0].1),
            PmiKind::Tolerance(characteristic) => reader.tolerance(id, characteristic),
        };
        let kind = entity.type_name();
        match read {
            Ok(item) if matches!(item.kind, ItemKind::Datum(_)) => datums.push((kind, item)),
            Ok(item) => later.push((kind, item)),
            Err(reason) => skip(id, kind, reason),
        }
    }
    // A tolerance needs its datums, so datums are settled first.
    datums.sort_by_key(|(_, item)| item.order);
    later.sort_by_key(|(_, item)| item.order);
    let mut validator = PartValidator::default();
    let mut conflicting_letters = BTreeSet::new();
    let mut recovered_names = Vec::new();
    let mut invalid_datum = None;
    let mut datum_named_item = false;
    for (kind, item) in &mut datums {
        let ItemKind::Datum(datum) = &item.kind else {
            unreachable!("datums were collected separately");
        };
        let name = item.name.take();
        let letter = datum.letter.as_str();
        let follows_conflict = conflicting_letters.contains(letter);
        let letter_error = datum.validate_letter().err();
        if follows_conflict && letter_error.is_some() {
            skip(
                item.order,
                kind,
                ItemError::DuplicateDatum(letter).to_string(),
            );
            continue;
        }
        if let Err(error) = validator.accept_datum(datum) {
            let reason = match error {
                // Remember these letters for diagnostic precedence only.
                // A later datum can still claim the letter on a free feature.
                ItemError::DuplicateDatumFeature if conflicting_letters.insert(letter) => {
                    "its feature already carries a datum".into()
                }
                ItemError::DuplicateDatumFeature => ItemError::DuplicateDatum(letter).to_string(),
                error => error.to_string(),
            };
            skip(item.order, kind, reason);
            continue;
        }
        let base = name.unwrap_or_else(|| format!("Datum {letter}"));
        if follows_conflict {
            // The recovered datum can satisfy references now. Give it a
            // name after the other items so it does not take an existing name.
            recovered_names.push((&mut item.name, base));
            continue;
        }
        let id = validator.assign_name(base);
        // Older imports kept such a datum, but its authoring error was the
        // first refusal for each later item. Keep that diagnostic here.
        if invalid_datum.is_none()
            && let Some(error) = letter_error
        {
            invalid_datum = Some(if id == "item" {
                error.to_string()
            } else {
                let name = id.strip_prefix("item: ").unwrap_or(&id);
                format!("{name}: {error}")
            });
        }
        datum_named_item |= id == "item";
        item.name = Some(id);
    }
    // A datum named `item` collided with the candidate's temporary name.
    let datum_error =
        invalid_datum.or_else(|| datum_named_item.then(|| ItemError::DuplicateName.to_string()));
    for (kind, item) in &mut later {
        let name = item.name.take();
        if let ItemKind::Tolerance(tolerance) = &item.kind
            && let Some(letter) = tolerance.datums.iter().find(|letter| {
                !validator.contains_datum(letter) && !conflicting_letters.contains(letter.as_str())
            })
        {
            skip(
                item.order,
                kind,
                format!("it refers to datum {letter}, which was not read"),
            );
            continue;
        }
        if let Some(reason) = &datum_error {
            skip(item.order, kind, reason.clone());
            continue;
        }
        if let Err(error) = validator.check_kind(&item.kind) {
            skip(item.order, kind, error.to_string());
            continue;
        }
        let base = name.unwrap_or_else(|| format!("#{}", item.order));
        item.name = Some(validator.assign_name(base));
    }
    for (name, base) in recovered_names {
        *name = Some(validator.assign_name(base));
    }
    drop(conflicting_letters);
    drop(validator);
    let items = datums
        .into_iter()
        .chain(later)
        .filter_map(|(_, item)| {
            Some(Item {
                id: item.name?,
                kind: item.kind,
            })
        })
        .collect();
    StepPmi {
        pmi: Pmi {
            items,
            part: reader.part_notes(scope),
        },
        skipped,
    }
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
