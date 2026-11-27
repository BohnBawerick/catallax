//! Manufacturing information read from a STEP file's AP242 semantic PMI.
//!
//! [`read_step_pmi`] returns the datums, dimensions, geometric tolerances and
//! surface textures a file attaches to faces, keyed by the entity names of its
//! `ADVANCED_FACE`s (`#17`), which are the face keys an import of the same
//! bytes reports. It reads the forms `catallax_step::write_step_with_pmi`
//! writes and the common forms of the CAx-IF recommended practices.
//!
//! Nothing is guessed. An item this reader cannot represent in full, such as
//! an angular size, an explicit area datum target or a tolerance with a
//! modifier it does not know, is left out and listed in [`StepPmi::skipped`]
//! with its entity and the reason. A datum whose target is left out is left
//! out too, with every item that refers to it.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Rational, Real};
use catallax_shapes::distance::{PlanarBoundary, PlanarEdge};
use catallax_shapes::{Circle, Plane, Point3, RPoint3, RVector3, Vector3};
use catallax_step::pmi::{
    Anchor, At, BEVEL_ANGLE, CHAMFER, COMPOSITE, COUNTERBORE, COUNTERSINK, Camera, ChamferNominal,
    Characteristic, Datum, DatumRef, DatumTarget, EDGE_BREAK, FIRST_SIDE, Feature,
    GeometricTolerance, HEAD_DEPTH, HEAD_DIAMETER, HOLE_DEPTH, HOLE_PILOT, HoleHead, HoleSize,
    INCLUDED_ANGLE, Item, ItemError, ItemKind, KNURL, Limits, Location, MAX_SETBACK, MIN_SETBACK,
    MaterialCondition, ORIENTATION, PartNotes, PartValidator, Pmi, PositionControl, PositionFrame,
    PositionSegment, Process, ProjectedZone, SECOND_SIDE, SEMANTIC_TEXT, SavedView, Size, SizeKind,
    SurfaceTexture, TARGET_DIAMETER, TARGET_LENGTH, TARGET_WIDTH, THREAD, THREAD_DEPTH, THROUGH,
    TargetArea, TargetKind, WIDTH, off_face,
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
    /// Graphic callouts that present no item that was read: graphical-only
    /// PMI, kept apart so nothing is counted twice.
    pub graphics: Vec<GraphicCallout>,
}

/// A graphic callout with no semantic item read from the file.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct GraphicCallout {
    /// The callout entity, such as `#412`.
    pub entity: String,
    /// The callout's name.
    pub name: String,
    /// The presented PMI type its geometry is named with, such as
    /// `position` (recommended practices, table 18), or empty.
    pub kind: String,
    /// The text it shows, from its equivalent Unicode string, if stated.
    pub text: Option<String>,
    /// The faces of the shape aspect it is associated with, if any.
    pub faces: Option<Feature<String>>,
    /// The semantic entity it presents when that entity was not read.
    pub presents: Option<String>,
    /// The saved views that show it.
    pub views: Vec<String>,
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
const UNREAD: [(&str, &str); 7] = [
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
    (
        "DATUM_TARGET",
        "a datum target with explicit geometry is not read",
    ),
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

/// A datum reference compartment's modifier: its material boundary, or none.
fn datum_boundary(param: Option<&Param>) -> Result<Option<MaterialCondition>, String> {
    let items = match param {
        None | Some(Param::Unset) => return Ok(None),
        Some(Param::List(items)) => items,
        Some(_) => return Err("a datum reference has an unreadable modifier".into()),
    };
    let mut found = None;
    for item in items {
        let name = match item {
            Param::Typed(kind, values) if kind == "SIMPLE_DATUM_REFERENCE_MODIFIER" => {
                match values.as_slice() {
                    [Param::Enum(name)] => name.as_str(),
                    _ => return Err("a datum reference has an unreadable modifier".into()),
                }
            }
            Param::Enum(name) => name.as_str(),
            Param::Typed(kind, _) => return Err(format!("a datum reference has a {kind}")),
            _ => return Err("a datum reference has an unreadable modifier".into()),
        };
        let material = match name {
            "MAXIMUM_MATERIAL_REQUIREMENT" => MaterialCondition::Maximum,
            "LEAST_MATERIAL_REQUIREMENT" => MaterialCondition::Least,
            other => return Err(format!("a datum reference has the modifier {other}")),
        };
        if found.replace(material).is_some() {
            return Err("a datum reference has more than one modifier".into());
        }
    }
    Ok(found)
}

fn sub3(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

fn neg3(a: &[Rational; 3]) -> [Rational; 3] {
    sub3(&[Rational::zero(), Rational::zero(), Rational::zero()], a)
}

fn dot3(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}

fn cross3(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

fn is_zero3(a: &[Rational; 3]) -> bool {
    a.iter().all(Rational::is_zero)
}

/// A nonzero vector square to `axis`.
fn perpendicular3(axis: &[Rational; 3]) -> [Rational; 3] {
    let [x, y, _] = axis;
    if x.is_zero() && y.is_zero() {
        [Rational::one(), Rational::zero(), Rational::zero()]
    } else {
        [&Rational::zero() - y, x.clone(), Rational::zero()]
    }
}

/// An `AXIS2_PLACEMENT_3D` with its defaults filled in, in file units.
struct Placement {
    name: String,
    origin: [Rational; 3],
    axis: [Rational; 3],
    reference: Option<[Rational; 3]>,
}

/// One geometric tolerance entity: a plain tolerance or one segment of a
/// position control.
struct Segment {
    value: Rational,
    diameter_zone: bool,
    material: Option<MaterialCondition>,
    datums: Vec<(DatumRef, u64)>,
    target: u64,
    projected: Option<(Rational, u64)>,
    name: Option<String>,
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
        let mut targeted = false;
        for user in self.sharing.get(&id).map_or(&[][..], Vec::as_slice) {
            let Some((_, relating, related)) = self.relationship(*user) else {
                continue;
            };
            if related != id {
                continue;
            }
            if self.has(relating, "DATUM_TARGET") {
                return Err("it is established by datum targets with explicit geometry".into());
            }
            if self.has(relating, "PLACED_DATUM_TARGET_FEATURE") {
                targeted = true;
                continue;
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
        if features.is_empty() && targeted {
            return Err(
                "it is set up only by datum targets, with no datum feature to name its faces"
                    .into(),
            );
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
        self.dimension_value_as(id, Self::length)
    }

    fn dimension_value_as(
        &self,
        id: u64,
        measure: fn(&Self, u64) -> Option<(String, Rational)>,
    ) -> Result<(Option<Rational>, Option<Limits>), String> {
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
                if let Some((name, value)) = measure(self, item) {
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

    fn degrees(&self, item: u64) -> Option<Rational> {
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
                        .degrees(angle)
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

    fn planar_feature(&self, aspect: u64, owner: u64) -> Result<Feature<String>, String> {
        if self.aspect_params(aspect).and_then(|p| reference(p.get(2))) != Some(owner) {
            return Err("a planar dimension's face group belongs to another part".into());
        }
        let (feature, derived) = self.end(aspect)?;
        if derived
            || !feature.iter().all(|face| {
                self.face_surface(face)
                    .is_some_and(|id| self.has(id, "PLANE"))
            })
        {
            return Err("the dimension needs planar faces, not derived shapes".into());
        }
        Ok(feature)
    }

    fn planar_size(&self, id: u64, aspect: u64, name: &str) -> Result<Found, String> {
        let owner = self
            .aspect_params(aspect)
            .and_then(|p| reference(p.get(2)))
            .ok_or("the planar size has no owning part")?;
        let links = self.sharing_at(aspect, "SHAPE_ASPECT_RELATIONSHIP", 2);
        let side = |role| {
            let mut found = links
                .iter()
                .filter(|(_, p)| matches!(p.first(), Some(Param::Str(name)) if name == role));
            let (_, params) = found
                .next()
                .ok_or_else(|| format!("the planar size has no {role}"))?;
            let target = reference(params.get(3))
                .filter(|target| *target != aspect)
                .ok_or("a planar size side names no separate face group")?;
            if found.next().is_some() {
                return Err(format!("the planar size has more than one {role}"));
            }
            self.planar_feature(target, owner)
        };
        let sides = [side(FIRST_SIDE)?, side(SECOND_SIDE)?];
        let (nominal, limits) = self.dimension_value(id)?;
        let nominal = nominal.ok_or("the planar size has no nominal length with a known unit")?;
        let attribute = |name| {
            self.hole_attribute(aspect, name)?
                .ok_or_else(|| format!("the planar size has no `{name}` attribute"))
        };
        let kind = match name {
            WIDTH => ItemKind::Width {
                flats: sides,
                limits,
                nominal: Some(nominal),
            },
            CHAMFER | EDGE_BREAK => {
                if limits.is_some() {
                    return Err("the bevel size carries unsupported dimensional limits".into());
                }
                let bevel = self.planar_feature(aspect, owner)?;
                if name == CHAMFER {
                    ItemKind::Chamfer {
                        bevel,
                        sides,
                        nominal: Some(ChamferNominal {
                            setback_mm: nominal,
                            angle_deg: self
                                .degrees(attribute(BEVEL_ANGLE)?)
                                .ok_or("the bevel angle is not in the supported degree unit")?,
                        }),
                    }
                } else {
                    ItemKind::EdgeBreak {
                        bevel,
                        sides,
                        min_mm: self.hole_length(attribute(MIN_SETBACK)?, MIN_SETBACK)?,
                        max_mm: self.hole_length(attribute(MAX_SETBACK)?, MAX_SETBACK)?,
                        nominal: Some(nominal),
                    }
                }
            }
            _ => unreachable!("planar size was selected by its name"),
        };
        Ok(Found {
            order: id,
            name: self.id_attribute(id),
            kind,
        })
    }

    fn size(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let name = string(params.get(1)).unwrap_or_default();
        let mut kind = match name.as_str() {
            WIDTH | CHAMFER | EDGE_BREAK => {
                let aspect = reference(params.first()).ok_or("it names no feature")?;
                return self.planar_size(id, aspect, &name);
            }
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
            kind = SizeKind::Hole(Box::new(HoleSize {
                head,
                drill_depth: drill
                    .map(|id| self.hole_length(id, HOLE_DEPTH))
                    .transpose()?,
                thread_depth: thread_depth
                    .map(|id| self.hole_length(id, THREAD_DEPTH))
                    .transpose()?,
            }));
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

    fn basic_angle(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let (name, a, b) = self.relationship(id).ok_or("it names no two planes")?;
        if name != "angular distance" {
            return Err(format!("a `{name}` angular location is not read"));
        }
        let (nominal, limits) = self.dimension_value_as(id, |reader, measure| {
            let params = reader.entity(measure)?.part("REPRESENTATION_ITEM")?;
            Some((string(params.first())?, reader.degrees(measure)?))
        })?;
        if limits != Some(Limits::Basic) {
            return Err("the angle is not a basic (theoretical) dimension".into());
        }
        let nominal = nominal.ok_or("the basic angle is not in the supported degree unit")?;
        let selection = if params.len() == 1 {
            params.first()
        } else {
            params.get(4)
        };
        if !matches!(selection, Some(Param::Enum(relator)) if match relator.as_str() {
            "SMALL" => nominal < Rational::from_i64(90),
            "EQUAL" => nominal == Rational::from_i64(90),
            "LARGE" => nominal > Rational::from_i64(90),
            _ => false,
        }) {
            return Err("the angle selection does not agree with the stated angle".into());
        }
        let owner = self
            .aspect_params(a)
            .and_then(|p| reference(p.get(2)))
            .ok_or("the angular location has no owning part")?;
        Ok(Found {
            order: id,
            name: self.id_attribute(id),
            kind: ItemKind::BasicAngle {
                planes: [
                    self.planar_feature(a, owner)?,
                    self.planar_feature(b, owner)?,
                ],
                nominal: Some(nominal),
            },
        })
    }

    /// The datum references a tolerance makes, primary first, each with the
    /// `DATUM` entity it names.
    fn datum_refs(&self, systems: &[u64]) -> Result<Vec<(DatumRef, u64)>, String> {
        let mut found = Vec::new();
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
                        let material = datum_boundary(fields.get(5))?;
                        let datum = reference(fields.get(4))
                            .filter(|datum| self.datum_letter(*datum).is_some())
                            .ok_or("it refers to a common datum")?;
                        let letter = self.datum_letter(datum).unwrap_or_default();
                        found.push((DatumRef { letter, material }, datum));
                    }
                }
                "DATUM_REFERENCE" => {
                    let rank = number(params.first()).ok_or("a datum has no precedence")?;
                    let datum = reference(params.get(1))
                        .filter(|datum| self.datum_letter(*datum).is_some())
                        .ok_or("a datum reference names no datum")?;
                    ranked.push((rank, datum));
                }
                other => return Err(format!("it refers to datums through a {other}")),
            }
        }
        ranked.sort();
        found.extend(ranked.into_iter().map(|(_, datum)| {
            let letter = self.datum_letter(datum).unwrap_or_default();
            (DatumRef::new(&letter), datum)
        }));
        Ok(found)
    }

    /// One geometric tolerance entity as the file states it.
    fn segment(&self, id: u64, characteristic: Characteristic) -> Result<Segment, String> {
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
        let mut projected = None;
        for (zone, params) in self.sharing_at(id, "TOLERANCE_ZONE", 4) {
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
            for (_, definition) in self.sharing_at(zone, "PROJECTED_ZONE_DEFINITION", 0) {
                if !references(definition.get(1)).is_empty() {
                    return Err("its projected zone names boundaries, which are not read".into());
                }
                let start =
                    reference(definition.get(2)).ok_or("its projected zone has no start")?;
                let height = reference(definition.get(3))
                    .and_then(|measure| self.length(measure))
                    .map(|(_, value)| value)
                    .ok_or("its projected zone's height is not a length with a known unit")?;
                if projected.replace((height, start)).is_some() {
                    return Err("it has more than one projected zone".into());
                }
            }
        }
        let name = string(base.first())
            .filter(|name| !name.is_empty())
            .or_else(|| self.id_attribute(id));
        Ok(Segment {
            value,
            diameter_zone,
            material,
            datums: self.datum_refs(&systems)?,
            target,
            projected,
            name,
        })
    }

    /// The upper segment of the composite frame `id` is the lower segment of.
    fn composite_upper(&self, id: u64) -> Option<u64> {
        self.sharing_at(id, "GEOMETRIC_TOLERANCE_RELATIONSHIP", 3)
            .first()
            .and_then(|(_, params)| reference(params.get(2)))
    }

    /// The lower segment of the composite frame whose upper segment is `id`.
    fn composite_lower(&self, id: u64) -> Result<Option<u64>, String> {
        let mut lowers = Vec::new();
        for (_, params) in self.sharing_at(id, "GEOMETRIC_TOLERANCE_RELATIONSHIP", 2) {
            let name = string(params.first()).unwrap_or_default();
            // CAx-IF writes `composite`; AP242's own example names it in full.
            if name != COMPOSITE && name != "composite tolerance" {
                return Err(format!(
                    "it is in a `{name}` tolerance relationship, which is not read"
                ));
            }
            lowers.extend(reference(params.get(3)));
        }
        match lowers.as_slice() {
            [] => Ok(None),
            [lower] => Ok(Some(*lower)),
            _ => Err("it has more than one lower composite segment".into()),
        }
    }

    /// Whether a datum's feature is a feature of size that can take a
    /// material boundary: one cylinder, or two opposed parallel planes. A
    /// datum that targets set up is not.
    fn datum_feature_of_size(&self, datum: u64) -> bool {
        if !self.placed_targets(datum).is_empty() {
            return false;
        }
        let Ok(Found {
            kind: ItemKind::Datum(Datum { feature, .. }),
            ..
        }) = self.datum(datum)
        else {
            return false;
        };
        if feature.iter().all(|face| self.round(face)) {
            return self.round_feature(&feature, false).is_ok();
        }
        let [a, b] = feature.as_slice() else {
            return false;
        };
        let (Some(a), Some(b)) = (self.outward_plane(a), self.outward_plane(b)) else {
            return false;
        };
        is_zero3(&cross3(&a.1, &b.1))
            && dot3(&a.1, &b.1) < Rational::zero()
            && !dot3(&sub3(&b.0, &a.0), &a.1).is_zero()
    }

    fn position_segment(&self, segment: Segment) -> Result<PositionSegment, String> {
        let mut datums = Vec::with_capacity(segment.datums.len());
        for (reference, datum) in segment.datums {
            if reference.material.is_some() && !self.datum_feature_of_size(datum) {
                return Err(format!(
                    "datum {} is not a feature of size and has no material boundary",
                    reference.letter
                ));
            }
            datums.push(reference);
        }
        Ok(PositionSegment {
            value_mm: segment.value,
            diameter_zone: segment.diameter_zone,
            material: segment.material,
            datums,
        })
    }

    fn position_features(&self, target: u64) -> Result<Vec<Feature<String>>, String> {
        let features = self.features(target, 0)?;
        for feature in &features {
            if self.round_feature(feature, false).is_err() {
                return Err("a position control's features are each one cylinder".into());
            }
        }
        Ok(features)
    }

    /// A projected zone from where each controlled feature meets the start
    /// faces (recommended practices, section 6.9.2.2). Each feature shares an
    /// edge with exactly one planar start face, whose outward normal gives the
    /// zone's direction.
    fn projection(
        &self,
        start: u64,
        height_mm: Rational,
        features: &[Feature<String>],
    ) -> Result<ProjectedZone<String>, String> {
        let starts = self
            .used_faces(start)?
            .ok_or("its projected zone's start names no face")?;
        let mut planes = Vec::with_capacity(starts.len());
        for face in &starts {
            let (_, normal) = self
                .outward_plane(face)
                .ok_or("its projected zone starts on a face that is not planar")?;
            planes.push((face, normal, self.face_edge_curves(face)));
        }
        let mut direction: Option<[Rational; 3]> = None;
        let mut entry_faces = Vec::with_capacity(features.len());
        for feature in features {
            let edges: BTreeSet<u64> = feature
                .iter()
                .flat_map(|face| self.face_edge_curves(face))
                .collect();
            let meeting: Vec<_> = planes
                .iter()
                .filter(|(_, _, own)| !own.is_disjoint(&edges))
                .collect();
            let [(face, normal, _)] = meeting.as_slice() else {
                return Err(
                    "each controlled feature needs exactly one projected zone start face it meets"
                        .into(),
                );
            };
            match &direction {
                None => direction = Some(normal.clone()),
                Some(first) => {
                    if !is_zero3(&cross3(first, normal)) || dot3(first, normal) <= Rational::zero()
                    {
                        return Err("the projected zone's start faces do not face one way".into());
                    }
                }
            }
            let axis = feature
                .first()
                .and_then(|face| self.face_surface(face))
                .and_then(|surface| self.entity(surface))
                .and_then(|surface| surface.part("CYLINDRICAL_SURFACE"))
                .and_then(|params| reference(params.get(1)))
                .and_then(|placement| self.placement(placement))
                .ok_or("a controlled feature has no readable axis")?
                .axis;
            if !is_zero3(&cross3(&axis, normal)) {
                return Err("a controlled axis is not square to its start face".into());
            }
            entry_faces.push((*face).clone());
        }
        Ok(ProjectedZone {
            height_mm,
            direction: direction.ok_or("the position control names no feature")?,
            entry_faces,
        })
    }

    fn tolerance(&self, id: u64, characteristic: Characteristic) -> Result<Found, String> {
        let lower = self.composite_lower(id)?;
        let upper = self.segment(id, characteristic)?;
        if let Some(lower) = lower {
            return self.composite(id, characteristic, upper, lower);
        }
        let boundary = upper
            .datums
            .iter()
            .any(|(reference, _)| reference.material.is_some());
        if characteristic == Characteristic::Position && (upper.projected.is_some() || boundary) {
            let features = self.position_features(upper.target)?;
            let projected_zone = match &upper.projected {
                Some((height, start)) => {
                    Some(self.projection(*start, height.clone(), &features)?)
                }
                None => None,
            };
            let name = upper.name.clone();
            return Ok(Found {
                order: id,
                name,
                kind: ItemKind::Position(PositionControl {
                    features,
                    frame: PositionFrame::Single {
                        segment: self.position_segment(upper)?,
                        projected_zone,
                    },
                }),
            });
        }
        if boundary {
            return Err("a datum reference has a modifier".into());
        }
        if upper.projected.is_some() {
            return Err("only a position tolerance's projected zone is read".into());
        }
        Ok(Found {
            order: id,
            kind: ItemKind::Tolerance(GeometricTolerance {
                characteristic,
                value: upper.value,
                diameter_zone: upper.diameter_zone,
                material: upper.material,
                datums: upper
                    .datums
                    .into_iter()
                    .map(|(reference, _)| reference.letter)
                    .collect(),
                features: self.features(upper.target, 0)?,
            }),
            name: upper.name,
        })
    }

    /// A composite position frame: the upper segment `id` and the lower one
    /// its `composite` relationship names (section 6.9.9).
    fn composite(
        &self,
        id: u64,
        characteristic: Characteristic,
        upper: Segment,
        lower: u64,
    ) -> Result<Found, String> {
        if characteristic != Characteristic::Position
            || self
                .entity(lower)
                .and_then(|entity| entity.part("POSITION_TOLERANCE"))
                .is_none()
        {
            return Err("only a composite position frame is read".into());
        }
        if self.composite_lower(lower)?.is_some() {
            return Err("its composite frame has more than two segments".into());
        }
        let lower = self
            .segment(lower, Characteristic::Position)
            .map_err(|reason| format!("its lower segment #{lower}: {reason}"))?;
        if upper.projected.is_some() || lower.projected.is_some() {
            return Err("a composite frame with a projected zone is not read".into());
        }
        let features = self.position_features(upper.target)?;
        if lower.target != upper.target && self.features(lower.target, 0)? != features {
            return Err("the composite frame's segments control different features".into());
        }
        let name = upper.name.clone();
        Ok(Found {
            order: id,
            name,
            kind: ItemKind::Position(PositionControl {
                features,
                frame: PositionFrame::Composite {
                    upper: self.position_segment(upper)?,
                    lower: self.position_segment(lower)?,
                },
            }),
        })
    }

    /// The placed datum targets that set up a datum.
    fn placed_targets(&self, datum: u64) -> Vec<u64> {
        let mut targets = Vec::new();
        for user in self.sharing.get(&datum).map_or(&[][..], Vec::as_slice) {
            if let Some((_, relating, related)) = self.relationship(*user)
                && related == datum
                && self.has(relating, "PLACED_DATUM_TARGET_FEATURE")
            {
                targets.push(relating);
            }
        }
        targets
    }

    /// An `AXIS2_PLACEMENT_3D` as the file writes it.
    fn placement(&self, id: u64) -> Option<Placement> {
        let params = self.entity(id)?.part("AXIS2_PLACEMENT_3D")?;
        let origin = self
            .entity(reference(params.get(1))?)?
            .part("CARTESIAN_POINT")?;
        let direction = |param: Option<&Param>| -> Option<Option<[Rational; 3]>> {
            match param {
                Some(Param::Unset) | None => Some(None),
                param => {
                    let found =
                        coordinates(self.entity(reference(param)?)?.part("DIRECTION")?.get(1))?;
                    (!is_zero3(&found)).then_some(Some(found))
                }
            }
        };
        Some(Placement {
            name: string(params.first()).unwrap_or_default(),
            origin: coordinates(origin.get(1))?,
            axis: direction(params.get(2))?.unwrap_or([
                Rational::zero(),
                Rational::zero(),
                Rational::one(),
            ]),
            reference: direction(params.get(3))?,
        })
    }

    /// The plane of a planar face and its normal pointing out of the material.
    fn outward_plane(&self, face: &str) -> Option<([Rational; 3], [Rational; 3])> {
        let id: u64 = face.strip_prefix('#')?.parse().ok()?;
        let params = self.entity(id)?.part("ADVANCED_FACE")?;
        let plane = self.entity(reference(params.get(2))?)?.part("PLANE")?;
        let placement = self.placement(reference(plane.get(1))?)?;
        let normal = match params.get(3) {
            Some(Param::Enum(sense)) if sense == "T" => placement.axis,
            Some(Param::Enum(sense)) if sense == "F" => neg3(&placement.axis),
            _ => return None,
        };
        Some((placement.origin, normal))
    }

    /// The `EDGE_CURVE`s that bound a face.
    fn face_edge_curves(&self, face: &str) -> BTreeSet<u64> {
        let mut edges = BTreeSet::new();
        let Some(params) = face
            .strip_prefix('#')
            .and_then(|id| id.parse().ok())
            .and_then(|id| self.entity(id))
            .and_then(|face| face.part("ADVANCED_FACE"))
        else {
            return edges;
        };
        for bound in references(params.get(1)) {
            let Some(cycle) = self
                .entity(bound)
                .and_then(|bound| reference(bound.parts[0].1.get(1)))
                .and_then(|cycle| self.entity(cycle))
                .and_then(|cycle| cycle.part("EDGE_LOOP"))
            else {
                continue;
            };
            for oriented in references(cycle.get(1)) {
                if let Some(edge) = self
                    .entity(oriented)
                    .and_then(|oriented| oriented.part("ORIENTED_EDGE"))
                    .and_then(|params| reference(params.get(3)))
                {
                    edges.insert(edge);
                }
            }
        }
        edges
    }

    /// Millimetres in the length unit of a representation context.
    fn context_millimetres(&self, context: u64) -> Option<Rational> {
        let units = self
            .entity(context)?
            .part("GLOBAL_UNIT_ASSIGNED_CONTEXT")?
            .first();
        references(units)
            .into_iter()
            .find_map(|unit| self.millimetres(unit, 0))
    }

    /// The complete trim of a planar face from the file's own edge records,
    /// in millimetres, or `None` when an edge is not a line or a circle.
    fn face_boundary(&self, face: &str, scale: &Rational) -> Option<PlanarBoundary> {
        let id: u64 = face.strip_prefix('#')?.parse().ok()?;
        let params = self.entity(id)?.part("ADVANCED_FACE")?;
        let plane = self.entity(reference(params.get(2))?)?.part("PLANE")?;
        let frame = self.placement(reference(plane.get(1))?)?;
        let at = |p: &[Rational; 3]| {
            Point3::new(
                Real::from_rational(&p[0] * scale),
                Real::from_rational(&p[1] * scale),
                Real::from_rational(&p[2] * scale),
            )
        };
        let along = |v: &[Rational; 3]| {
            Vector3::new(
                Real::from_rational(v[0].clone()),
                Real::from_rational(v[1].clone()),
                Real::from_rational(v[2].clone()),
            )
        };
        let mut vertices = Vec::new();
        let mut indices = BTreeMap::new();
        let mut edges = Vec::new();
        for bound in references(params.get(1)) {
            let cycle = self.entity(reference(self.entity(bound)?.parts[0].1.get(1))?)?;
            if cycle.part("VERTEX_LOOP").is_some() {
                continue;
            }
            for oriented in references(cycle.part("EDGE_LOOP")?.get(1)) {
                let edge = reference(self.entity(oriented)?.part("ORIENTED_EDGE")?.get(3))?;
                let edge = self.entity(edge)?.part("EDGE_CURVE")?;
                let mut ends = [0usize; 2];
                for (slot, vertex) in ends.iter_mut().zip([edge.get(1), edge.get(2)]) {
                    let vertex = reference(vertex)?;
                    *slot = match indices.get(&vertex) {
                        Some(index) => *index,
                        None => {
                            let point =
                                reference(self.entity(vertex)?.part("VERTEX_POINT")?.get(1))?;
                            let point =
                                coordinates(self.entity(point)?.part("CARTESIAN_POINT")?.get(1))?;
                            vertices.push(at(&point));
                            indices.insert(vertex, vertices.len() - 1);
                            vertices.len() - 1
                        }
                    };
                }
                let same_sense = match edge.get(4) {
                    Some(Param::Enum(sense)) if sense == "T" => true,
                    Some(Param::Enum(sense)) if sense == "F" => false,
                    _ => return None,
                };
                let mut curve = self.entity(reference(edge.get(3))?)?;
                for wrapper in ["SURFACE_CURVE", "SEAM_CURVE"] {
                    if let Some(params) = curve.part(wrapper) {
                        curve = self.entity(reference(params.get(1))?)?;
                    }
                }
                if curve.part("LINE").is_some() {
                    edges.push(PlanarEdge::Segment(ends));
                } else {
                    // A face bounded by another curve is not checked.
                    let circle = curve.part("CIRCLE")?;
                    let placement = self.placement(reference(circle.get(1))?)?;
                    let radius = &number(circle.get(2))? * scale;
                    let ref_dir = placement
                        .reference
                        .clone()
                        .unwrap_or_else(|| perpendicular3(&placement.axis));
                    edges.push(PlanarEdge::Arc {
                        circle: Circle {
                            center: at(&placement.origin),
                            normal: along(&placement.axis),
                            radius_sq: Real::from_rational(&radius * &radius),
                            ref_dir: along(&ref_dir),
                        },
                        ends: (ends[0] != ends[1])
                            .then(|| if same_sense { ends } else { [ends[1], ends[0]] }),
                    });
                }
            }
        }
        let ref_dir = frame
            .reference
            .clone()
            .unwrap_or_else(|| perpendicular3(&frame.axis));
        (!edges.is_empty()).then(|| PlanarBoundary {
            plane: Plane {
                origin: at(&frame.origin),
                normal: along(&frame.axis),
                ref_dir: along(&ref_dir),
            },
            vertices,
            edges,
        })
    }

    /// The equivalent Unicode string a callout's validation property states.
    fn callout_text(&self, callout: u64) -> Option<String> {
        for (within, _) in self.sharing_at(callout, "CHARACTERIZED_ITEM_WITHIN_REPRESENTATION", 2) {
            for (property, _) in self.sharing_at(within, "PROPERTY_DEFINITION", 2) {
                for (_, link) in self.sharing_at(property, "PROPERTY_DEFINITION_REPRESENTATION", 0)
                {
                    let Some(representation) =
                        reference(link.get(1)).and_then(|id| self.entity(id))
                    else {
                        continue;
                    };
                    for item in references(representation.parts[0].1.get(1)) {
                        if let Some(fields) = self
                            .entity(item)
                            .and_then(|entity| entity.part("DESCRIPTIVE_REPRESENTATION_ITEM"))
                            && string(fields.first()).as_deref()
                                == Some("equivalent unicode string")
                        {
                            return string(fields.get(1));
                        }
                    }
                }
            }
        }
        None
    }

    /// The presented PMI type a callout's geometry is named with.
    fn callout_kind(&self, callout: u64) -> String {
        let Some(entity) = self.entity(callout) else {
            return String::new();
        };
        let params = entity
            .part("DRAUGHTING_CALLOUT")
            .unwrap_or(entity.parts[0].1.as_slice());
        for occurrence in references(params.get(1)) {
            let Some(occurrence) = self.entity(occurrence) else {
                continue;
            };
            let fields = occurrence.parts[0].1.as_slice();
            if let Some(name) = reference(fields.get(2))
                .and_then(|set| self.entity(set))
                .and_then(|set| string(set.parts[0].1.first()))
                .filter(|name| !name.is_empty())
            {
                return name;
            }
        }
        String::new()
    }

    /// A saved view's camera (recommended practices, section 9.4.2.6), with
    /// lengths in millimetres.
    fn camera(&self, id: u64, scale: &Rational) -> Result<Camera, String> {
        let entity = self.entity(id).ok_or("the camera is missing")?;
        if entity.parts.len() != 1 || entity.type_name() != "CAMERA_MODEL_D3" {
            return Err(format!("a {} is not read", entity.type_name()));
        }
        let params = entity.parts[0].1.as_slice();
        let frame = reference(params.get(1))
            .and_then(|frame| self.placement(frame))
            .ok_or("the camera has no readable placement")?;
        let across = frame
            .reference
            .ok_or("the camera's placement has no x axis")?;
        let norm = dot3(&frame.axis, &frame.axis);
        let up = cross3(&frame.axis, &across).map(|value| &value / &norm);
        let volume = reference(params.get(2))
            .and_then(|volume| self.entity(volume))
            .and_then(|volume| volume.part("VIEW_VOLUME"))
            .ok_or("the camera has no view volume")?;
        let perspective = match volume.first() {
            Some(Param::Enum(kind)) if kind == "CENTRAL" => true,
            Some(Param::Enum(kind)) if kind == "PARALLEL" => false,
            _ => return Err("the camera's projection is neither central nor parallel".into()),
        };
        let distance = number(volume.get(2)).ok_or("the camera has no view distance")?;
        let window = reference(volume.get(8))
            .and_then(|window| self.entity(window))
            .and_then(|window| window.part("PLANAR_BOX"))
            .ok_or("the camera has no view window")?;
        let size = |index: usize| {
            number(window.get(index))
                .map(|value| &value * scale)
                .ok_or("the camera's window has no size")
        };
        Ok(Camera {
            eye: frame.origin.map(|value| &value * scale),
            direction: frame.axis,
            up,
            distance: &distance * scale,
            window: [size(1)?, size(2)?],
            perspective,
        })
    }

    /// The saved views and graphic callouts of the file. A callout whose
    /// association names a read item, or the private aspect of a read surface
    /// texture, presents that item and is not listed again as graphics.
    fn presentation(
        &self,
        items: &[(u64, Item<String>)],
        lowers: &BTreeMap<u64, u64>,
        skipped: &mut Vec<SkippedPmi>,
    ) -> (Vec<SavedView>, Vec<GraphicCallout>) {
        let mut by_entity: BTreeMap<u64, &str> = items
            .iter()
            .map(|(entity, item)| (*entity, item.id.as_str()))
            .collect();
        for (lower, upper) in lowers {
            if let Some(name) = by_entity.get(upper).copied() {
                by_entity.insert(*lower, name);
            }
        }
        let is_callout = |entity: &Entity| {
            entity
                .parts
                .iter()
                .any(|(part, _)| part == "DRAUGHTING_CALLOUT" || part.ends_with("_CALLOUT"))
        };
        let mut presented = BTreeMap::new();
        let mut graphics = BTreeMap::new();
        for (&callout, entity) in self.g {
            if !is_callout(entity) {
                continue;
            }
            let mut item = None;
            let mut presents = None;
            let mut faces = None;
            for (_, association) in self.sharing_at(callout, "DRAUGHTING_MODEL_ITEM_ASSOCIATION", 4)
            {
                let Some(definition) = reference(association.get(2)) else {
                    continue;
                };
                if let Some(name) = by_entity.get(&definition) {
                    item = Some((*name).to_owned());
                    break;
                }
                if let Some(params) = self.aspect_params(definition) {
                    let feature: Option<Feature<String>> = self
                        .features(definition, 0)
                        .ok()
                        .map(|features| features.into_iter().flatten().collect());
                    let name = string(params.first()).unwrap_or_default();
                    if let Some((_, texture)) = items.iter().find(|(_, texture)| {
                        texture.id == name
                            && matches!(&texture.kind, ItemKind::Texture(t)
                            if feature.as_ref().is_some_and(|faces| {
                                let mut a = faces.clone();
                                let mut b = t.feature.clone();
                                a.sort();
                                b.sort();
                                a == b
                            }))
                    }) {
                        item = Some(texture.id.clone());
                        break;
                    }
                    faces = feature;
                } else {
                    presents = Some(format!("#{definition}"));
                }
            }
            match item {
                Some(name) => {
                    presented.insert(callout, name);
                }
                None => {
                    let params = entity.parts[0].1.as_slice();
                    graphics.insert(
                        callout,
                        GraphicCallout {
                            entity: format!("#{callout}"),
                            name: string(params.first()).unwrap_or_default(),
                            kind: self.callout_kind(callout),
                            text: self.callout_text(callout),
                            faces,
                            presents,
                            views: Vec::new(),
                        },
                    );
                }
            }
        }
        let mut views: Vec<SavedView> = Vec::new();
        for (&model, entity) in self.g {
            if entity.type_name() != "DRAUGHTING_MODEL" || entity.parts.len() != 1 {
                continue;
            }
            let params = entity.parts[0].1.as_slice();
            let members = references(params.get(1));
            let cameras: Vec<u64> = members
                .iter()
                .copied()
                .filter(|member| {
                    self.entity(*member)
                        .is_some_and(|e| e.type_name().starts_with("CAMERA_MODEL"))
                })
                .collect();
            if cameras.is_empty() {
                continue;
            }
            let mut shown = Vec::new();
            for member in &members {
                let Some(entity) = self.entity(*member) else {
                    continue;
                };
                if is_callout(entity) {
                    shown.push(*member);
                } else if let Some(plane) = entity.part("ANNOTATION_PLANE") {
                    shown.extend(
                        references(plane.get(3))
                            .into_iter()
                            .filter(|element| self.entity(*element).is_some_and(is_callout)),
                    );
                }
            }
            let Some(scale) =
                reference(params.get(2)).and_then(|context| self.context_millimetres(context))
            else {
                skipped.push(SkippedPmi {
                    entity: format!("#{model}"),
                    kind: "DRAUGHTING_MODEL".into(),
                    reason: "the saved view has no known length unit".into(),
                });
                continue;
            };
            for camera in cameras {
                let camera_name = self
                    .entity(camera)
                    .and_then(|camera| string(camera.parts[0].1.first()))
                    .filter(|name| !name.is_empty());
                let mut name = camera_name
                    .or_else(|| string(params.first()))
                    .unwrap_or_default();
                if views.iter().any(|view| view.name == name) {
                    name = format!("{name} ({})", string(params.first()).unwrap_or_default());
                }
                let camera_record = match self.camera(camera, &scale) {
                    Ok(camera) => camera,
                    Err(reason) => {
                        skipped.push(SkippedPmi {
                            entity: format!("#{camera}"),
                            kind: "CAMERA_MODEL_D3".into(),
                            reason,
                        });
                        continue;
                    }
                };
                let mut view_items: Vec<String> = Vec::new();
                for callout in &shown {
                    if let Some(item) = presented.get(callout) {
                        if !view_items.contains(item) {
                            view_items.push(item.clone());
                        }
                    } else if let Some(graphic) = graphics.get_mut(callout)
                        && !graphic.views.contains(&name)
                    {
                        graphic.views.push(name.clone());
                    }
                }
                views.push(SavedView {
                    name,
                    camera: camera_record,
                    items: view_items,
                });
            }
        }
        (views, graphics.into_values().collect())
    }

    /// Editable text on faces (recommended practices, section 7.4): a face
    /// note whose lines have no variable name, or a knurl designation held
    /// under [`KNURL`].
    fn text(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let target = reference(params.get(2)).ok_or("it is attached to nothing")?;
        match self.entity(target).map(Entity::type_name) {
            Some("PRODUCT_DEFINITION" | "PRODUCT_DEFINITION_SHAPE") => {
                return Err("a text note on the whole part is not read".into());
            }
            _ if self.aspect_params(target).is_none() => {
                return Err("a text note on a PMI element is not read".into());
            }
            _ => {}
        }
        let mut lines = Vec::new();
        for (_, link) in self.sharing_at(id, "PROPERTY_DEFINITION_REPRESENTATION", 0) {
            let representation = reference(link.get(1))
                .and_then(|id| self.entity(id))
                .ok_or("its text is missing")?;
            for item in references(representation.parts[0].1.get(1)) {
                let fields = self
                    .entity(item)
                    .and_then(|entity| entity.part("DESCRIPTIVE_REPRESENTATION_ITEM"))
                    .ok_or("its text holds an item that is not a string")?;
                lines.push((
                    string(fields.first()).unwrap_or_default(),
                    string(fields.get(1)).unwrap_or_default(),
                ));
            }
        }
        let feature = self.one_feature(target)?;
        let kind = match lines.as_slice() {
            [] => return Err("it holds no text".into()),
            [(name, designation)] if name == KNURL => ItemKind::Knurl {
                designation: designation.clone(),
                feature,
            },
            _ => {
                if let Some((name, _)) = lines.iter().find(|(name, _)| !name.is_empty()) {
                    return Err(format!("a text variable `{name}` is not read"));
                }
                let text: Vec<&str> = lines.iter().map(|(_, line)| line.as_str()).collect();
                ItemKind::FaceNote {
                    text: text.join("\n"),
                    feature,
                }
            }
        };
        Ok(Found {
            order: id,
            name: self.id_attribute(id),
            kind,
        })
    }

    /// A placed datum target (recommended practices, section 6.6): its datum,
    /// its one face, its `orientation` and sizes. The placement's axis must
    /// point out of the face, the point must lie on the face's plane, and the
    /// whole target must lie on the face within its trim.
    fn target(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let description = string(params.get(1)).unwrap_or_default();
        let number = string(params.get(4))
            .and_then(|text| text.trim().parse::<u32>().ok())
            .filter(|number| *number > 0)
            .ok_or("its target number is not a positive whole number")?;
        let mut datums = Vec::new();
        for user in self.sharing.get(&id).map_or(&[][..], Vec::as_slice) {
            if let Some((_, relating, related)) = self.relationship(*user)
                && relating == id
                && self.datum_letter(related).is_some()
            {
                datums.push(related);
            }
        }
        let [datum] = datums.as_slice() else {
            return Err("it does not belong to exactly one datum".into());
        };
        let letter = self.datum_letter(*datum).unwrap_or_default();
        let links = self.sharing_at(id, "FEATURE_FOR_DATUM_TARGET_RELATIONSHIP", 3);
        let [(_, link)] = links.as_slice() else {
            return Err(if links.is_empty() {
                "it is tied to no face".into()
            } else {
                "it is tied to more than one face feature".into()
            });
        };
        let aspect = reference(link.get(2)).ok_or("its face feature is missing")?;
        if self.aspect_params(aspect).and_then(|p| reference(p.get(2))) != reference(params.get(2))
        {
            return Err("its face feature belongs to another part".into());
        }
        let faces = self
            .used_faces(aspect)?
            .ok_or("its face feature names no face")?;
        let [face] = faces.as_slice() else {
            return Err("its face feature names more than one face".into());
        };
        let mut representations = Vec::new();
        for (property, _) in self.sharing_at(id, "PROPERTY_DEFINITION", 2) {
            for (_, link) in self.sharing_at(property, "SHAPE_DEFINITION_REPRESENTATION", 0) {
                representations.extend(reference(link.get(1)));
            }
        }
        let [representation] = representations.as_slice() else {
            return Err("its parameters are not one shape representation".into());
        };
        let representation = self
            .entity(*representation)
            .and_then(|entity| entity.part("SHAPE_REPRESENTATION_WITH_PARAMETERS"))
            .ok_or("its parameters are not a shape representation with parameters")?;
        let scale = reference(representation.get(2))
            .and_then(|context| self.context_millimetres(context))
            .ok_or("its parameters have no known length unit")?;
        let (mut orientation, mut length, mut width, mut diameter) = (None, None, None, None);
        for item in references(representation.get(1)) {
            if let Some(placement) = self.placement(item) {
                if placement.name != ORIENTATION {
                    return Err(format!(
                        "its parameters include a placement named `{}`",
                        placement.name
                    ));
                }
                if orientation.replace(placement).is_some() {
                    return Err("it has more than one orientation".into());
                }
                continue;
            }
            let Some((name, value)) = self.length(item) else {
                let kind = self
                    .entity(item)
                    .map_or("missing entity", Entity::type_name);
                return Err(format!("its parameters include a {kind}"));
            };
            let slot = match name.as_str() {
                TARGET_LENGTH => &mut length,
                TARGET_WIDTH => &mut width,
                TARGET_DIAMETER => &mut diameter,
                other => return Err(format!("its parameters include `{other}`")),
            };
            if slot.replace(value).is_some() {
                return Err(format!("it states `{name}` twice"));
            }
        }
        let orientation = orientation.ok_or("it has no orientation")?;
        let direction = || {
            orientation
                .reference
                .clone()
                .ok_or_else(|| format!("the {description} target has no direction"))
        };
        let kind = match (description.as_str(), length, width, diameter) {
            ("point", None, None, None) => TargetKind::Point,
            ("line", Some(length_mm), None, None) => TargetKind::Line {
                length_mm,
                direction: direction()?,
            },
            ("circle", None, None, Some(diameter_mm)) => {
                TargetKind::Area(TargetArea::Circular { diameter_mm })
            }
            ("rectangle", Some(length_mm), Some(width_mm), None) => {
                TargetKind::Area(TargetArea::Rectangular {
                    length_mm,
                    width_mm,
                    direction: direction()?,
                })
            }
            ("point" | "line" | "circle" | "rectangle", ..) => {
                return Err(format!("its sizes do not fit a {description} target"));
            }
            (other, ..) => return Err(format!("a `{other}` datum target is not read")),
        };
        let (origin, outward) = self
            .outward_plane(face)
            .ok_or("a datum target is read only on a planar face")?;
        if !is_zero3(&cross3(&orientation.axis, &outward))
            || dot3(&orientation.axis, &outward) <= Rational::zero()
        {
            return Err("its orientation does not point out of its face".into());
        }
        if !dot3(&sub3(&orientation.origin, &origin), &outward).is_zero() {
            return Err("its point is not on its face's plane".into());
        }
        let target = DatumTarget {
            datum: letter.clone(),
            number,
            face: face.clone(),
            point_mm: orientation.origin.clone().map(|value| &value * &scale),
            kind,
        };
        let boundary = self.face_boundary(face, &scale).ok_or(
            "its face's boundary is not lines and circles, so the target cannot be checked on it",
        )?;
        match off_face(&boundary, &target, &mut Decider::default()) {
            Ok(None) => {}
            Ok(Some(reason)) => return Err(reason.into()),
            Err(refusal) => {
                return Err(format!(
                    "the target could not be checked on its face: {}",
                    refusal.detail
                ));
            }
        }
        Ok(Found {
            order: id,
            name: Some(
                self.id_attribute(id)
                    .unwrap_or_else(|| format!("Datum target {letter}{number}")),
            ),
            kind: ItemKind::DatumTarget(target),
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
/// bytes gives its faces. A dimension's nominal is the value the file states,
/// in millimetres or degrees; the caller decides what to do when the geometry
/// disagrees.
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
    Target,
    Text,
    Size,
    Location,
    BasicAngle,
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
    } else if name == "PLACED_DATUM_TARGET_FEATURE" && entity.parts.len() == 1 {
        PmiKind::Target
    } else if name == "PROPERTY_DEFINITION"
        && entity.parts.len() == 1
        && string(entity.parts[0].1.first()).as_deref() == Some(SEMANTIC_TEXT)
    {
        PmiKind::Text
    } else if name == "DIMENSIONAL_SIZE" && entity.parts.len() == 1 {
        PmiKind::Size
    } else if entity.part("ANGULAR_LOCATION").is_some() {
        PmiKind::BasicAngle
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
    // A lower composite segment is read with its upper one.
    let mut lowers = Vec::new();
    let mut failed = BTreeSet::new();
    for Candidate { id, entity, kind } in candidates {
        let read = match kind {
            PmiKind::Unread(part, reason) => {
                skip(id, part, reason.to_owned());
                continue;
            }
            PmiKind::Datum => reader.datum(id),
            PmiKind::Target => reader.target(id, &entity.parts[0].1),
            PmiKind::Text => reader.text(id, &entity.parts[0].1),
            PmiKind::Size => reader.size(id, &entity.parts[0].1),
            PmiKind::Location => reader.location(id),
            PmiKind::BasicAngle => reader.basic_angle(id, entity.part("ANGULAR_LOCATION").unwrap()),
            PmiKind::Texture => reader.texture(id, &entity.parts[0].1),
            PmiKind::Tolerance(characteristic) => {
                if let Some(upper) = reader.composite_upper(id) {
                    lowers.push((id, entity.type_name(), upper));
                    continue;
                }
                reader.tolerance(id, characteristic)
            }
        };
        let kind = entity.type_name();
        match read {
            Ok(item) if matches!(item.kind, ItemKind::Datum(_)) => datums.push((kind, item)),
            Ok(item) => later.push((kind, item)),
            Err(reason) => {
                failed.insert(id);
                skip(id, kind, reason);
            }
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
        // A datum its targets set up is not realized from its whole feature
        // when one of them is missing.
        if let Some(target) = reader
            .placed_targets(item.order)
            .into_iter()
            .find(|target| failed.contains(target))
        {
            item.name = None;
            skip(
                item.order,
                kind,
                format!("its datum target #{target} was not read"),
            );
            continue;
        }
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
        let referenced: Vec<&str> = match &item.kind {
            ItemKind::Tolerance(tolerance) => tolerance.datums.iter().map(String::as_str).collect(),
            ItemKind::Position(control) => control
                .frame
                .segments()
                .into_iter()
                .flat_map(|segment| segment.datums.iter().map(|datum| datum.letter.as_str()))
                .collect(),
            ItemKind::DatumTarget(target) => vec![target.datum.as_str()],
            _ => Vec::new(),
        };
        if let Some(letter) = referenced.into_iter().find(|letter| {
            !validator.contains_datum(letter) && !conflicting_letters.contains(*letter)
        }) {
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
        if let ItemKind::DatumTarget(target) = &item.kind
            && let Err(error) = validator.accept_target(target)
        {
            skip(item.order, kind, error.to_string());
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
    let mut paired_lowers = BTreeMap::new();
    for (id, kind, upper) in lowers {
        paired_lowers.insert(id, upper);
        let read = later.iter().any(|(_, item)| {
            item.order == upper
                && item.name.is_some()
                && matches!(
                    &item.kind,
                    ItemKind::Position(PositionControl {
                        frame: PositionFrame::Composite { .. },
                        ..
                    })
                )
        });
        if !read {
            skip(
                id,
                kind,
                format!(
                    "it is the lower segment of the composite frame #{upper}, which was not read"
                ),
            );
        }
    }
    drop(conflicting_letters);
    drop(validator);
    // Each read item with the entity it was read from.
    let items: Vec<(u64, Item<String>)> = datums
        .into_iter()
        .chain(later)
        .filter_map(|(_, item)| {
            Some((
                item.order,
                Item {
                    id: item.name?,
                    kind: item.kind,
                },
            ))
        })
        .collect();
    let (views, graphics) = match scope {
        None => reader.presentation(&items, &paired_lowers, &mut skipped),
        // Saved views and graphics of an assembly's parts are not read.
        Some(_) => (Vec::new(), Vec::new()),
    };
    StepPmi {
        pmi: Pmi {
            items: items.into_iter().map(|(_, item)| item).collect(),
            part: reader.part_notes(scope),
            views,
        },
        skipped,
        graphics,
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
