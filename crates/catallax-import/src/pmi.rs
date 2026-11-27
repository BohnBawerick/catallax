//! AP242 manufacturing information in a STEP file.
//!
//! Use [`read_step_pmi`] to get a [`StepPmi`] report and [`patch_step_pmi`]
//! to apply source-record edits.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::distance::{PlanarBoundary, PlanarEdge};
use catallax_shapes::{Circle, Plane, Point3, Vector3};
use catallax_step::pmi::{
    Anchor, At, BEVEL_ANGLE, CHAMFER, COMPOSITE, COMPOSITE_STANDARD, COUNTERBORE, COUNTERSINK,
    Camera, ChamferNominal, Characteristic, Datum, DatumRef, DatumTarget, EDGE_BREAK, FIRST_SIDE,
    Feature, GeometricTolerance, HEAD_DEPTH, HEAD_DIAMETER, HOLE_DEPTH, HOLE_PILOT, HoleHead,
    HoleSize, INCLUDED_ANGLE, Item, ItemError, ItemKind, KNURL, Limits, Location, MAX_SETBACK,
    MIN_SETBACK, MaterialCondition, ORIENTATION, PartNotes, PartValidator, Pmi, PositionControl,
    PositionFrame, PositionSegment, Process, ProjectedZone, SECOND_SIDE, SEMANTIC_TEXT, SavedView,
    Size, SizeKind, SurfaceTexture, TARGET_DIAMETER, TARGET_LENGTH, TARGET_WIDTH, THREAD,
    THREAD_DEPTH, THROUGH, TargetArea, TargetKind, WIDTH, off_face, projection_circle_behind,
    projection_misses_cylinder, projection_plane_hit,
};

use crate::ImportRefusal;
use crate::part21::{Entity, Graph, Param, parse};
use crate::recognize::{Vec3, length_measure, placement_reference, unit_millimetres};

/// The manufacturing information of a STEP file.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct StepPmi {
    /// What was read, with faces named by their entity, such as `#17`.
    pub pmi: Pmi<String>,
    /// Diagnostics for unsupported or malformed items and presentation.
    pub skipped: Vec<SkippedPmi>,
    /// Graphics without a semantic item and graphics whose associations include
    /// an unread semantic item.
    pub graphics: Vec<GraphicCallout>,
    /// The entity each item was read from, such as `#412`, by item name: a
    /// tolerance, dimension, datum, target or property. A source patch starts
    /// from it.
    pub sources: BTreeMap<String, String>,
}

/// A graphic callout with no semantic item or with an unread association.
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
    /// The semantic entities it presents that were not read.
    pub presents: Vec<String>,
    /// Polyline strokes in part coordinates, in millimetres.
    pub polylines_mm: Vec<Vec<[Rational; 3]>>,
    /// Annotation planes that contain the callout.
    pub planes: Vec<GraphicPlane>,
    /// The saved views that show it.
    pub views: Vec<String>,
}

/// An annotation plane in part coordinates.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct GraphicPlane {
    /// The plane's origin in millimetres.
    pub origin_mm: [Rational; 3],
    /// Its normal.
    pub normal: [Rational; 3],
    /// Its reference direction.
    pub reference: [Rational; 3],
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

fn reference_list(param: Option<&Param>) -> Option<Vec<u64>> {
    let Some(Param::List(items)) = param else {
        return None;
    };
    if items.is_empty() {
        return None;
    }
    let refs: Vec<_> = items
        .iter()
        .map(|item| reference(Some(item)))
        .collect::<Option<_>>()?;
    (refs.iter().collect::<BTreeSet<_>>().len() == refs.len()).then_some(refs)
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
fn decode(raw: &str) -> Option<String> {
    let mut out = String::new();
    let mut rest = raw;
    while let Some(at) = rest.find('\\') {
        out.push_str(&rest[..at]);
        let tail = &rest[at..];
        if let Some(body) = tail.strip_prefix("\\\\") {
            out.push('\\');
            rest = body;
        } else if let Some(body) = tail.strip_prefix("\\PA\\") {
            rest = body;
        } else if let Some(body) = tail.strip_prefix("\\S\\") {
            let character = *body.as_bytes().first()?;
            if !(b' '..=b'~').contains(&character) {
                return None;
            }
            out.push(char::from_u32(u32::from(character) + 128)?);
            rest = &body[1..];
        } else if let Some(body) = tail.strip_prefix("\\X\\") {
            let code = u32::from_str_radix(body.get(..2)?, 16).ok()?;
            out.push(char::from_u32(code)?);
            rest = &body[2..];
        } else {
            let (width, body) = if let Some(body) = tail.strip_prefix("\\X2\\") {
                (4, body)
            } else {
                (8, tail.strip_prefix("\\X4\\")?)
            };
            let (encoded, after) = body.split_once("\\X0\\")?;
            if encoded.is_empty() || encoded.len() % width != 0 {
                return None;
            }
            let codes = encoded
                .as_bytes()
                .chunks_exact(width)
                .map(|chunk| u32::from_str_radix(std::str::from_utf8(chunk).ok()?, 16).ok());
            if width == 4 {
                let units = codes
                    .map(|code| code.map(|code| code as u16))
                    .collect::<Option<Vec<_>>>()?;
                for character in char::decode_utf16(units) {
                    out.push(character.ok()?);
                }
            } else {
                for code in codes {
                    out.push(char::from_u32(code?)?);
                }
            }
            rest = after;
        }
    }
    out.push_str(rest);
    Some(out)
}

pub(crate) fn assembly_name(param: Option<&Param>) -> Option<String> {
    if let Some(Param::Str(raw)) = param {
        Some(decode(raw).unwrap_or_else(|| raw.clone()))
    } else {
        None
    }
}

fn string(param: Option<&Param>) -> Option<String> {
    match param {
        Some(Param::Str(text)) => decode(text),
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
fn datum_boundary(param: &Param) -> Result<Option<MaterialCondition>, String> {
    let items = match param {
        Param::Unset => return Ok(None),
        Param::List(items) if !items.is_empty() => items,
        _ => return Err("a datum reference has an unreadable modifier".into()),
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

/// An `AXIS2_PLACEMENT_3D` with its defaults filled in, in file units.
struct Placement {
    name: String,
    origin: [Rational; 3],
    axis: [Rational; 3],
    reference: [Rational; 3],
}

enum ZoneForm {
    Cylinder,
    Circle,
    Other,
}

/// One geometric tolerance entity: a plain tolerance or one segment of a
/// position control.
struct Segment {
    value: Rational,
    zone_form: ZoneForm,
    material: Option<MaterialCondition>,
    datums: Vec<(DatumRef, u64)>,
    target: u64,
    projected: Option<(Rational, u64)>,
    name: Option<String>,
}

struct RoundSurface {
    placement: u64,
    axis: Placement,
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
    datum_sources: Vec<u64>,
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

    fn read_round(&self, id: u64, cone: bool) -> Option<RoundSurface> {
        let kind = if cone {
            "CONICAL_SURFACE"
        } else {
            "CYLINDRICAL_SURFACE"
        };
        let params = self.entity(id)?.part(kind)?;
        let placement = reference(params.get(1))?;
        let axis = self.placement(placement)?;
        let radius = number(params.get(2))?;
        let angle = if cone {
            Some(number(params.get(3))?)
        } else {
            None
        };
        (cone || radius > Rational::zero()).then_some(RoundSurface {
            placement,
            axis,
            radius,
            angle,
        })
    }

    fn round_feature(&self, faces: &[String], cone: bool) -> Result<RoundSurface, String> {
        let mut surfaces = faces.iter().map(|face| {
            self.face_surface(face)
                .ok_or("a size face has no readable surface")
        });
        let first = surfaces.next().ok_or("the size has no face")??;
        let bad_surface = if cone {
            "the head is not conical"
        } else {
            "a size face is not cylindrical"
        };
        let shape = self.read_round(first, cone).ok_or(bad_surface)?;
        for surface in surfaces {
            let surface = surface?;
            if surface == first {
                continue;
            }
            let other = self.read_round(surface, cone).ok_or(bad_surface)?;
            if shape.radius != other.radius
                || shape.angle != other.angle
                || (if cone {
                    shape.placement != other.placement
                } else {
                    !is_zero3(&cross3(&shape.axis.axis, &other.axis.axis))
                        || !is_zero3(&cross3(
                            &sub3(&other.axis.origin, &shape.axis.origin),
                            &shape.axis.axis,
                        ))
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
            let [value, unit] = measure else { return None };
            let name = entity
                .part("REPRESENTATION_ITEM")
                .and_then(|params| string(params.first()));
            (name, Some(value), Some(unit))
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
        let scale = unit_millimetres(self.g, reference(unit)?, 0).ok()?;
        let (_, value) = length_measure(value?)?;
        Some((name.unwrap_or_default(), &value * &scale))
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
        let targeted = !self.placed_targets(id)?.is_empty();
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
            datum_sources: Vec::new(),
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
            datum_sources: Vec::new(),
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
            datum_sources: Vec::new(),
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
            datum_sources: Vec::new(),
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
            datum_sources: Vec::new(),
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

    fn datum_refs(&self, systems: &[u64], owner: u64) -> Result<Vec<(DatumRef, u64)>, String> {
        let resolve = |param| {
            let datum = reference(param).ok_or("a datum reference names no datum")?;
            let letter = self
                .datum_letter(datum)
                .ok_or("it refers to a common or unreadable datum")?;
            if self
                .entity(datum)
                .and_then(|e| reference(e.parts[0].1.get(2)))
                != Some(owner)
            {
                return Err("a datum reference belongs to another part");
            }
            Ok((letter, datum))
        };
        if let [system] = systems
            && let Some(fields) = self
                .entity(*system)
                .filter(|entity| entity.parts.len() == 1)
                .and_then(|entity| entity.part("DATUM_SYSTEM"))
        {
            if fields.len() != 5
                || reference(fields.get(2)) != Some(owner)
                || !matches!(fields.get(3), Some(Param::Enum(flag)) if flag == "F")
            {
                return Err("its datum system has an invalid shape or attributes".into());
            }
            let compartments = reference_list(fields.get(4))
                .filter(|refs| refs.len() <= 3)
                .ok_or("its datum system needs one to three distinct ordered compartments")?;
            return compartments
                .into_iter()
                .map(|compartment| {
                    let fields = self
                        .entity(compartment)
                        .filter(|entity| entity.parts.len() == 1)
                        .and_then(|entity| entity.part("DATUM_REFERENCE_COMPARTMENT"))
                        .ok_or("a datum reference is not a single compartment")?;
                    let [_, _, shape, Param::Enum(flag), datum, modifiers] = fields else {
                        return Err("a datum reference compartment has invalid attributes".into());
                    };
                    if reference(Some(shape)) != Some(owner)
                        || !matches!(flag.as_str(), "T" | "F" | "U")
                        || self.sharing_at(compartment, "DATUM_SYSTEM", 4).len() != 1
                    {
                        return Err(
                            "a datum reference compartment has invalid attributes or ownership"
                                .into(),
                        );
                    }
                    let material = datum_boundary(modifiers)?;
                    let (letter, datum) = resolve(Some(datum))?;
                    Ok((DatumRef { letter, material }, datum))
                })
                .collect();
        }
        let mut ranked = Vec::new();
        for system in systems {
            let fields = self
                .entity(*system)
                .filter(|entity| entity.parts.len() == 1)
                .and_then(|entity| entity.part("DATUM_REFERENCE"))
                .ok_or("a tolerance needs one datum system or ranked datum references")?;
            let [Param::Number(rank), datum] = fields else {
                return Err("a legacy datum reference has invalid attributes".into());
            };
            let rank = rank
                .parse::<usize>()
                .map_err(|_| "a datum's precedence is not an integer")?;
            let (letter, datum) = resolve(Some(datum))?;
            ranked.push((rank, DatumRef::new(&letter), datum));
        }
        ranked.sort_by_key(|(rank, _, _)| *rank);
        if ranked
            .iter()
            .enumerate()
            .any(|(index, (rank, _, _))| *rank != index + 1)
        {
            return Err("datum precedence must be distinct and start at one".into());
        }
        Ok(ranked
            .into_iter()
            .map(|(_, reference, datum)| (reference, datum))
            .collect())
    }

    /// One geometric tolerance entity as the file states it.
    fn segment(&self, id: u64, characteristic: Characteristic) -> Result<Segment, String> {
        let entity = self.entity(id).ok_or("it is missing")?;
        let own = characteristic.entity();
        let mut seen = BTreeSet::new();
        for (part, _) in &entity.parts {
            let known = [
                own,
                "GEOMETRIC_TOLERANCE",
                "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE",
                "GEOMETRIC_TOLERANCE_WITH_MODIFIERS",
                "MODIFIED_GEOMETRIC_TOLERANCE",
            ];
            if !seen.insert(part.as_str()) {
                return Err(format!("it states {part} more than once"));
            }
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
        let datum_part = entity.part("GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE");
        let (base_count, datum_field) = match (simple, datum_part, characteristic.datum_use()) {
            (true, None, catallax_step::pmi::DatumUse::Required) => (5, base.get(4)),
            (true, None, _) => (4, None),
            (false, Some([field]), _) => (4, Some(field)),
            (false, None, catallax_step::pmi::DatumUse::Required) => {
                return Err("it lacks its datum reference attributes".into());
            }
            (false, None, _) => (4, None),
            _ => return Err("its datum reference has invalid attributes".into()),
        };
        if base.len() != base_count || (!simple && entity.part(own) != Some(&[])) {
            return Err("its geometric tolerance has invalid attributes".into());
        }
        let systems = match datum_field {
            Some(field) => {
                reference_list(Some(field)).ok_or("it has an invalid datum reference aggregate")?
            }
            None => Vec::new(),
        };
        let mut material = None;
        let mut has_modifier = false;
        for part in [
            "GEOMETRIC_TOLERANCE_WITH_MODIFIERS",
            "MODIFIED_GEOMETRIC_TOLERANCE",
        ] {
            let Some(params) = entity.part(part) else {
                continue;
            };
            let [modifiers] = params else {
                return Err("it has an unreadable modifier".into());
            };
            let modifiers = match modifiers {
                Param::List(items)
                    if part == "GEOMETRIC_TOLERANCE_WITH_MODIFIERS" && !items.is_empty() =>
                {
                    items.as_slice()
                }
                Param::Enum(_) if part == "MODIFIED_GEOMETRIC_TOLERANCE" => {
                    std::slice::from_ref(modifiers)
                }
                _ => return Err("it has an unreadable modifier".into()),
            };
            for modifier in modifiers {
                let Param::Enum(name) = modifier else {
                    return Err("it has an unreadable modifier".into());
                };
                if std::mem::replace(&mut has_modifier, true) {
                    return Err("its material modifier is stated more than once".into());
                }
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
        let datums = if systems.is_empty() {
            Vec::new()
        } else {
            let owner = self
                .aspect_params(target)
                .and_then(|fields| reference(fields.get(2)))
                .ok_or("the controlled feature has no owning part")?;
            self.datum_refs(&systems, owner)?
        };
        let mut zones = Vec::new();
        for user in self.sharing.get(&id).into_iter().flatten() {
            let Some(zone) = self.entity(*user) else {
                continue;
            };
            let simple = zone.parts.len() == 1;
            let specialized = [
                "DIRECTED_TOLERANCE_ZONE",
                "ORIENTED_TOLERANCE_ZONE",
                "DIRECTION_FEATURE_TOLERANCE_ZONE",
                "TOLERANCE_ZONE_WITH_DATUM",
            ];
            let fields = match zone.part("TOLERANCE_ZONE") {
                Some(fields) => fields,
                None if simple && specialized.contains(&zone.type_name()) => {
                    zone.parts[0].1.as_slice()
                }
                None if specialized.iter().any(|kind| zone.part(kind).is_some()) => {
                    return Err("its tolerance zone has unsupported shape attributes".into());
                }
                _ => continue,
            };
            let defining = if simple { 4 } else { 0 };
            let shape = if simple && zone.type_name() == "TOLERANCE_ZONE" {
                fields.get(..4)
            } else if !simple && zone.parts.len() == 2 {
                zone.part("SHAPE_ASPECT")
            } else {
                None
            };
            let Some(
                [
                    Param::Str(_),
                    Param::Str(_) | Param::Unset,
                    owner,
                    Param::Enum(flag),
                ],
            ) = shape
            else {
                return Err("its tolerance zone has unsupported shape attributes".into());
            };
            if !matches!(flag.as_str(), "T" | "F" | "U")
                || reference(Some(owner)).filter(|id| self.has(*id, "PRODUCT_DEFINITION_SHAPE"))
                    != self.aspect_params(target).and_then(|p| reference(p.get(2)))
                || reference(Some(owner)).is_none()
            {
                return Err("its tolerance zone has invalid ownership or flag".into());
            }
            let members = reference_list(fields.get(defining))
                .ok_or("its zone has an invalid defining tolerance list")?;
            if fields.len() != defining + 2
                || !members.contains(&id)
                || members.into_iter().any(|member| {
                    self.entity(member).is_none_or(|entity| {
                        !entity.parts.iter().any(|(kind, _)| {
                            kind == "GEOMETRIC_TOLERANCE"
                                || Characteristic::of_entity(kind).is_some()
                        })
                    })
                })
            {
                return Err("its zone has an invalid defining tolerance list".into());
            }
            zones.push((*user, reference(fields.get(defining + 1))));
        }
        if zones.len() > 1 {
            return Err("it has more than one tolerance zone".into());
        }
        let mut zone_form = ZoneForm::Other;
        let mut projected = None;
        if let [(zone, form)] = zones.as_slice() {
            let form = form
                .and_then(|form| self.entity(form))
                .and_then(|form| form.part("TOLERANCE_ZONE_FORM"))
                .and_then(|form| match form {
                    [Param::Str(label)] => decode(label),
                    _ => None,
                })
                .ok_or("its zone has no readable form")?;
            zone_form = match form.as_str() {
                "cylindrical or circular" | "within a cylinder" => ZoneForm::Cylinder,
                "within a circle" => ZoneForm::Circle,
                "spherical" | "within a sphere" => return Err("its zone is a sphere".into()),
                _ => ZoneForm::Other,
            };
            projected = self.projected_definition(*zone)?;
            if projected.is_some() && !matches!(zone_form, ZoneForm::Cylinder) {
                return Err("its projected zone needs a cylindrical zone".into());
            }
        }
        let name = string(base.first())
            .filter(|name| !name.is_empty())
            .or_else(|| self.id_attribute(id));
        Ok(Segment {
            value,
            zone_form,
            material,
            datums,
            target,
            projected,
            name,
        })
    }

    fn projected_definition(&self, zone: u64) -> Result<Option<(Rational, u64)>, String> {
        let mut projected = None;
        for user in self.sharing.get(&zone).into_iter().flatten() {
            let Some(entity) = self.entity(*user) else {
                continue;
            };
            let kind = entity.type_name();
            let simple = entity.parts.len() == 1;
            let base = if simple
                && matches!(
                    kind,
                    "TOLERANCE_ZONE_DEFINITION"
                        | "PROJECTED_ZONE_DEFINITION"
                        | "PROJECTED_ZONE_DEFINITION_WITH_OFFSET"
                        | "NON_UNIFORM_ZONE_DEFINITION"
                        | "RUNOUT_ZONE_DEFINITION"
                ) {
                Some(entity.parts[0].1.as_slice())
            } else {
                entity.part("TOLERANCE_ZONE_DEFINITION")
            };
            let Some(base) = base else {
                if entity.part("PROJECTED_ZONE_DEFINITION").is_some() {
                    return Err("its projected zone lacks a tolerance zone definition".into());
                }
                continue;
            };
            if reference(base.first()) != Some(zone) {
                return Err("its projected zone must name one tolerance zone".into());
            }
            let fields = match (simple, kind, entity.part("PROJECTED_ZONE_DEFINITION")) {
                (true, "PROJECTED_ZONE_DEFINITION", _) => base.get(2..),
                (false, _, Some(fields)) if entity.parts.len() == 2 && base.len() == 2 => {
                    Some(fields)
                }
                _ => {
                    return Err(format!(
                        "its zone definition uses {kind}, which is not read"
                    ));
                }
            };
            let Some([start, height]) = fields else {
                return Err("its projected zone has unreadable parameters".into());
            };
            if !matches!(base.get(1), Some(Param::List(bounds)) if bounds.is_empty()) {
                return Err("its projected zone names boundaries, which are not read".into());
            }
            let start = reference(Some(start)).ok_or("its projected zone has no start")?;
            let height = reference(Some(height))
                .and_then(|measure| self.length(measure))
                .map(|(_, value)| value)
                .ok_or("its projected zone's height is not a length with a known unit")?;
            if projected.replace((height, start)).is_some() {
                return Err("it has more than one projected zone".into());
            }
        }
        Ok(projected)
    }

    fn tolerance_relationships(&self, id: u64) -> impl Iterator<Item = &'a [Param]> + '_ {
        self.sharing
            .get(&id)
            .into_iter()
            .flatten()
            .filter_map(|user| self.entity(*user)?.part("GEOMETRIC_TOLERANCE_RELATIONSHIP"))
    }

    fn composite_link(&self, id: u64) -> Result<Option<(u64, u64)>, String> {
        let mut links = self.tolerance_relationships(id);
        let Some(params) = links.next() else {
            return Ok(None);
        };
        let [Param::Str(raw), Param::Str(_) | Param::Unset, upper, lower] = params else {
            return Err("its tolerance relationship has invalid attributes".into());
        };
        let name = decode(raw).ok_or("its tolerance relationship has no readable name")?;
        if name != COMPOSITE {
            return Err(format!(
                "it is in a `{name}` tolerance relationship, which is not read"
            ));
        }
        let upper = reference(Some(upper))
            .filter(|id| self.has(*id, "POSITION_TOLERANCE"))
            .ok_or("its composite frame has no upper position segment")?;
        let lower = reference(Some(lower))
            .filter(|id| self.has(*id, "POSITION_TOLERANCE"))
            .ok_or("its composite frame has no lower position segment")?;
        let other = if id == upper { lower } else { upper };
        if upper == lower
            || (id != upper && id != lower)
            || links.next().is_some()
            || self.tolerance_relationships(other).count() != 1
        {
            return Err("its composite frame does not have two uniquely paired segments".into());
        }
        Ok(Some((upper, lower)))
    }

    /// Whether a datum's feature is a feature of size that can take a
    /// material boundary: one cylinder, or two opposed parallel planes. A
    /// datum that targets set up is not.
    fn datum_feature_of_size(&self, datum: u64) -> bool {
        if !matches!(self.placed_targets(datum), Ok(targets) if targets.is_empty()) {
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
        if !matches!(segment.zone_form, ZoneForm::Cylinder) {
            return Err("a position control's zone is cylindrical".into());
        }
        Ok(PositionSegment {
            value_mm: segment.value,
            diameter_zone: true,
            material: segment.material,
            datums,
        })
    }

    fn position_features(
        &self,
        target: u64,
    ) -> Result<(Vec<Feature<String>>, Vec<Placement>), String> {
        let features = self.features(target, 0)?;
        let axes = features
            .iter()
            .map(|feature| {
                self.round_feature(feature, false)
                    .map(|round| round.axis)
                    .map_err(|_| "a position control's features are each one cylinder")
            })
            .collect::<Result<_, _>>()?;
        Ok((features, axes))
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
        axes: &[Placement],
    ) -> Result<ProjectedZone<String>, String> {
        let starts = self
            .used_faces(start)?
            .ok_or("its projected zone's start names no face")?;
        let mut planes = Vec::with_capacity(starts.len());
        for face in &starts {
            let (point, normal) = self
                .outward_plane(face)
                .ok_or("its projected zone starts on a face that is not planar")?;
            planes.push((face, point, normal, self.face_edge_curves(face)));
        }
        let mut direction: Option<[Rational; 3]> = None;
        let mut entry_faces = Vec::with_capacity(features.len());
        for (feature, axis) in features.iter().zip(axes) {
            let edges: BTreeSet<u64> = feature
                .iter()
                .flat_map(|face| self.face_edge_curves(face))
                .collect();
            let meeting: Vec<_> = planes
                .iter()
                .filter(|(_, _, _, own)| !own.is_disjoint(&edges))
                .collect();
            let [(face, point, normal, _)] = meeting.as_slice() else {
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
            if !is_zero3(&cross3(&axis.axis, normal)) {
                return Err("a controlled axis is not square to its start face".into());
            }
            let offset = &dot3(&sub3(point, &axis.origin), normal) / &dot3(&axis.axis, normal);
            let intersection: [Rational; 3] =
                std::array::from_fn(|i| &axis.origin[i] + &(&offset * &axis.axis[i]));
            self.check_inward_trim(&edges, point, normal)?;
            self.check_open_entry(face, &intersection, normal)?;
            entry_faces.push((*face).clone());
        }
        if starts.iter().any(|face| !entry_faces.contains(face)) {
            return Err("its projected zone names an unused start face".into());
        }
        Ok(ProjectedZone {
            height_mm,
            direction: direction.ok_or("the position control names no feature")?,
            entry_faces,
        })
    }

    fn check_inward_trim(
        &self,
        edges: &BTreeSet<u64>,
        origin: &[Rational; 3],
        normal: &[Rational; 3],
    ) -> Result<(), String> {
        let unknown = "the cylinder's axial trim cannot be checked";
        let mut inward = false;
        let mut check = |point: &[Rational; 3]| -> Result<(), String> {
            let distance = dot3(&sub3(point, origin), normal);
            if distance > Rational::zero() {
                return Err("the projected zone points back into its controlled cylinder".into());
            }
            inward |= distance < Rational::zero();
            Ok(())
        };
        for id in edges {
            let edge = self
                .entity(*id)
                .and_then(|entity| entity.part("EDGE_CURVE"))
                .ok_or(unknown)?;
            let curve = self
                .trim_curve(reference(edge.get(3)).ok_or(unknown)?)
                .ok_or(unknown)?;
            if curve.part("LINE").is_some() {
                for vertex in [edge.get(1), edge.get(2)] {
                    let point = reference(vertex)
                        .and_then(|id| self.entity(id)?.part("VERTEX_POINT"))
                        .and_then(|fields| reference(fields.get(1)))
                        .and_then(|id| self.point_coordinates(id))
                        .ok_or(unknown)?;
                    check(&point)?;
                }
            } else {
                let circle = curve
                    .part("CIRCLE")
                    .and_then(|p| reference(p.get(1)))
                    .and_then(|id| self.placement(id))
                    .ok_or(unknown)?;
                if !is_zero3(&cross3(&circle.axis, normal)) {
                    return Err(unknown.into());
                }
                check(&circle.origin)?;
            }
        }
        if !inward {
            return Err(unknown.into());
        }
        Ok(())
    }

    fn check_open_entry(
        &self,
        entry: &str,
        intersection: &[Rational; 3],
        normal: &[Rational; 3],
    ) -> Result<(), String> {
        let entry_id: u64 = entry[1..].parse().expect("a used face names its entity");
        let shells = self.sharing_at(entry_id, "CLOSED_SHELL", 1);
        let [(_, shell)] = shells.as_slice() else {
            return Err("its projected zone entry has no unique closed shell".into());
        };
        let unknown = "the projected zone's exterior start cannot be proved";
        let point = |p: [Rational; 3]| {
            let [x, y, z] = p.map(Real::from_rational);
            Point3::new(x, y, z)
        };
        let vector = |p: [Rational; 3]| {
            let [x, y, z] = p.map(Real::from_rational);
            Vector3::new(x, y, z)
        };
        let at = point(intersection.clone());
        let direction = vector(normal.clone());
        let mut dec = Decider::default();
        for id in reference_list(shell.get(1)).ok_or(unknown)? {
            let face = format!("#{id}");
            let surface = self.face_surface(&face).ok_or(unknown)?;
            if self.has(surface, "PLANE") {
                let (origin, axis) = self.outward_plane(&face).ok_or(unknown)?;
                let plane = Plane {
                    origin: point(origin),
                    normal: vector(axis),
                    ref_dir: direction.clone(),
                };
                if let Ok(hit) = projection_plane_hit(&at, &direction, &plane, &mut dec) {
                    if let Some(hit) = hit
                        && self
                            .face_boundary(&face, &Rational::one())
                            .ok_or(unknown)?
                            .contains(&hit, &mut dec)
                            .map_err(|error| error.detail)?
                    {
                        return Err("its projected zone's axis meets a closed end".into());
                    }
                    continue;
                }
            } else if self.has(surface, "CYLINDRICAL_SURFACE") {
                let round = self.read_round(surface, false).ok_or(unknown)?;
                let cylinder = catallax_shapes::Cylinder {
                    origin: point(round.axis.origin),
                    axis: vector(round.axis.axis),
                    ref_dir: vector(round.axis.reference),
                    radius_sq: Real::from_rational(&round.radius * &round.radius),
                };
                if projection_misses_cylinder(&at, &direction, &cylinder, &mut dec)
                    .map_err(|error| error.detail)?
                {
                    continue;
                }
            } else if !self.has(surface, "CONICAL_SURFACE") {
                return Err(unknown.into());
            }
            let edges = self.face_edge_curves(&face);
            if edges.is_empty() {
                return Err(unknown.into());
            }
            for id in edges {
                let edge = self
                    .entity(id)
                    .and_then(|edge| edge.part("EDGE_CURVE"))
                    .ok_or(unknown)?;
                let curve = self
                    .trim_curve(reference(edge.get(3)).ok_or(unknown)?)
                    .ok_or(unknown)?;
                if curve.part("LINE").is_some() {
                    for vertex in [edge.get(1), edge.get(2)] {
                        let end = reference(vertex)
                            .and_then(|id| self.entity(id)?.part("VERTEX_POINT"))
                            .and_then(|fields| reference(fields.get(1)))
                            .and_then(|id| self.point_coordinates(id))
                            .ok_or(unknown)?;
                        if dot3(&sub3(&end, intersection), normal) > Rational::zero() {
                            return Err(unknown.into());
                        }
                    }
                } else {
                    let circle = curve.part("CIRCLE").ok_or(unknown)?;
                    let frame = reference(circle.get(1))
                        .and_then(|id| self.placement(id))
                        .ok_or(unknown)?;
                    let radius = number(circle.get(2))
                        .filter(|radius| *radius > Rational::zero())
                        .ok_or(unknown)?;
                    if !projection_circle_behind(
                        &at,
                        &direction,
                        &point(frame.origin),
                        &vector(frame.axis),
                        &Real::from_rational(&radius * &radius),
                        &mut dec,
                    )
                    .map_err(|error| error.detail)?
                    {
                        return Err(unknown.into());
                    }
                }
            }
        }
        Ok(())
    }

    fn tolerance(
        &self,
        id: u64,
        characteristic: Characteristic,
        lower: Option<u64>,
    ) -> Result<Found, String> {
        let upper = self.segment(id, characteristic)?;
        if let Some(lower) = lower {
            return self.composite(id, characteristic, upper, lower);
        }
        if !characteristic.supports_asme_2018() && self.aspect_has_asme_standard(upper.target) {
            return Err(format!(
                "this tolerance is not supported under {COMPOSITE_STANDARD}"
            ));
        }
        let datum_sources = upper.datums.iter().map(|(_, source)| *source).collect();
        let boundary = upper
            .datums
            .iter()
            .any(|(reference, _)| reference.material.is_some());
        if characteristic == Characteristic::Position && (upper.projected.is_some() || boundary) {
            let (features, axes) = self.position_features(upper.target)?;
            let projected_zone = match &upper.projected {
                Some((height, start)) => {
                    Some(self.projection(*start, height.clone(), &features, &axes)?)
                }
                None => None,
            };
            let name = upper.name.clone();
            return Ok(Found {
                order: id,
                datum_sources,
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
            datum_sources,
            kind: ItemKind::Tolerance(GeometricTolerance {
                characteristic,
                value: upper.value,
                diameter_zone: matches!(upper.zone_form, ZoneForm::Cylinder | ZoneForm::Circle),
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

    fn aspect_has_asme_standard(&self, target: u64) -> bool {
        let Some(product) = self
            .aspect_params(target)
            .and_then(|p| reference(p.get(2)))
            .and_then(|shape| self.entity(shape)?.part("PRODUCT_DEFINITION_SHAPE"))
            .and_then(|shape| reference(shape.get(2)))
        else {
            return false;
        };
        self.has_asme_standard(product)
    }

    fn has_asme_standard(&self, product: u64) -> bool {
        let mut standards = BTreeSet::new();
        for (_, fields) in self.sharing_at(product, "PRODUCT_DEFINITION_CONTEXT_ASSOCIATION", 0) {
            let Some(context) = reference(fields.get(1)) else {
                continue;
            };
            let application = self
                .entity(context)
                .and_then(|e| e.part("PRODUCT_DEFINITION_CONTEXT"))
                .and_then(|p| reference(p.get(1)))
                .and_then(|id| self.entity(id)?.part("APPLICATION_CONTEXT"))
                .and_then(|p| string(p.first()));
            if application.as_deref()
                != Some("geometrical dimensioning and tolerancing representation")
            {
                continue;
            }
            for (_, fields) in self.sharing_at(context, "APPLIED_DOCUMENT_REFERENCE", 2) {
                let Some(name) = reference(fields.first())
                    .and_then(|id| self.entity(id)?.part("DOCUMENT"))
                    .and_then(|p| string(p.first()))
                else {
                    return false;
                };
                standards.insert(name);
            }
        }
        standards.len() == 1 && standards.contains(COMPOSITE_STANDARD)
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
        let lower = self
            .segment(lower, Characteristic::Position)
            .map_err(|reason| format!("its lower segment #{lower}: {reason}"))?;
        if upper.projected.is_some() || lower.projected.is_some() {
            return Err("a composite frame with a projected zone is not read".into());
        }
        if !self.aspect_has_asme_standard(upper.target)
            || !self.aspect_has_asme_standard(lower.target)
        {
            return Err(format!(
                "a composite position frame requires a part declaration of {COMPOSITE_STANDARD}"
            ));
        }
        let (features, _) = self.position_features(upper.target)?;
        if lower.target != upper.target && self.features(lower.target, 0)? != features {
            return Err("the composite frame's segments control different features".into());
        }
        let datum_sources = upper
            .datums
            .iter()
            .chain(&lower.datums)
            .map(|(_, source)| *source)
            .collect();
        let name = upper.name.clone();
        Ok(Found {
            order: id,
            datum_sources,
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
    fn placed_targets(&self, datum: u64) -> Result<Vec<u64>, String> {
        self.sharing
            .get(&datum)
            .into_iter()
            .flatten()
            .filter_map(|id| self.entity(*id)?.part("SHAPE_ASPECT_RELATIONSHIP"))
            .filter_map(|fields| {
                let from = fields.get(2);
                let declared =
                    string(fields.get(1)).as_deref() == Some(catallax_step::pmi::DATUM_TARGET);
                let placed = fields.iter().any(|field| {
                    references(Some(field))
                        .into_iter()
                        .any(|id| self.has(id, "PLACED_DATUM_TARGET_FEATURE"))
                });
                (declared || placed).then(|| {
                    if !declared || fields.len() != 4 || reference(fields.get(3)) != Some(datum) {
                        return Err("its datum target relationship is unreadable".into());
                    }
                    reference(from)
                        .ok_or_else(|| "its datum target relationship has no target".into())
                })
            })
            .collect()
    }

    /// An `AXIS2_PLACEMENT_3D` as the file writes it.
    fn placement(&self, id: u64) -> Option<Placement> {
        let [Param::Str(name), origin, axis, reference_dir] =
            self.entity(id)?.part("AXIS2_PLACEMENT_3D")?
        else {
            return None;
        };
        let direction = |param: &Param| -> Option<Option<[Rational; 3]>> {
            match param {
                Param::Unset => Some(None),
                param => {
                    let [Param::Str(_), ratios] =
                        self.entity(reference(Some(param))?)?.part("DIRECTION")?
                    else {
                        return None;
                    };
                    let found = coordinates(Some(ratios))?;
                    (!is_zero3(&found)).then_some(Some(found))
                }
            }
        };
        let axis =
            direction(axis)?.unwrap_or([Rational::zero(), Rational::zero(), Rational::one()]);
        let reference_dir = placement_reference(&axis, direction(reference_dir)?)?;
        Some(Placement {
            name: decode(name).unwrap_or_default(),
            origin: self.point_coordinates(reference(Some(origin))?)?,
            axis,
            reference: reference_dir,
        })
    }

    /// The plane of a planar face and its normal pointing out of the material.
    fn outward_plane(&self, face: &str) -> Option<([Rational; 3], [Rational; 3])> {
        let id: u64 = face.strip_prefix('#')?.parse().ok()?;
        let params = self.entity(id)?.part("ADVANCED_FACE")?;
        let plane = self.entity(reference(params.get(2))?)?.part("PLANE")?;
        let placement = self.placement(reference(plane.get(1))?)?;
        let Some(Param::Enum(sense)) = params.get(3) else {
            return None;
        };
        let mut same = match sense.as_str() {
            "T" => true,
            "F" => false,
            _ => return None,
        };
        let shells = self.sharing_at(id, "CLOSED_SHELL", 1);
        if shells.len() > 1 {
            return None;
        }
        if let [(shell, _)] = shells.as_slice() {
            let oriented = self.sharing_at(*shell, "ORIENTED_CLOSED_SHELL", 2);
            match oriented.as_slice() {
                [] => {}
                [(_, [_, Param::Derived, _, Param::Enum(flag)])] if flag == "T" => {}
                [(_, [_, Param::Derived, _, Param::Enum(flag)])] if flag == "F" => same = !same,
                _ => return None,
            }
        }
        let normal = if same {
            placement.axis
        } else {
            neg3(&placement.axis)
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
    fn context_millimetres(&self, context: u64) -> Result<Rational, String> {
        let units = self
            .entity(context)
            .and_then(|entity| entity.part("GLOBAL_UNIT_ASSIGNED_CONTEXT"))
            .and_then(|fields| reference_list(fields.first()))
            .ok_or("the context has no known length unit")?;
        let mut scale = None;
        for unit in units {
            let entity = self
                .entity(unit)
                .ok_or("the context refers to a missing unit")?;
            if entity.part("LENGTH_UNIT").is_none() {
                continue;
            }
            let length = unit_millimetres(self.g, unit, 0)
                .map_err(|_| "the context has an unreadable length unit")?;
            if scale.as_ref().is_some_and(|known| known != &length) {
                return Err("the context has conflicting length units".into());
            }
            scale = Some(length);
        }
        scale.ok_or_else(|| "the context has no known length unit".into())
    }

    fn representation_contains_face(&self, representation: u64, face: u64) -> bool {
        let Some(rep) = self.entity(representation) else {
            return false;
        };
        let mut pending = references(rep.parts[0].1.get(1));
        let mut seen = BTreeSet::new();
        while let Some(id) = pending.pop() {
            if id == face {
                return true;
            }
            if !seen.insert(id) {
                continue;
            }
            let Some(entity) = self.entity(id) else {
                continue;
            };
            if matches!(
                entity.type_name(),
                "MANIFOLD_SOLID_BREP"
                    | "BREP_WITH_VOIDS"
                    | "BREP_WITH_VOIDS_AND_BOUNDS"
                    | "ORIENTED_CLOSED_SHELL"
                    | "CLOSED_SHELL"
                    | "OPEN_SHELL"
                    | "SHELL_BASED_SURFACE_MODEL"
            ) {
                collect_references(&entity.parts[0].1, &mut pending);
            }
        }
        false
    }

    fn representation_scale(&self, representation: u64) -> Result<Rational, String> {
        let context = self
            .entity(representation)
            .and_then(|entity| {
                entity
                    .part("REPRESENTATION")
                    .or_else(|| (entity.parts.len() == 1).then(|| entity.parts[0].1.as_slice()))
            })
            .and_then(|fields| reference(fields.get(2)))
            .ok_or("the representation has no unit context")?;
        self.context_millimetres(context)
    }

    fn presentation_scale(&self, model: u64) -> Result<Rational, String> {
        let scale = self.representation_scale(model)?;
        let members = self
            .entity(model)
            .and_then(Self::draughting_model)
            .and_then(|fields| reference_list(fields.get(1)))
            .ok_or("the draughting model has unreadable members")?;
        let maps: BTreeSet<_> = members
            .into_iter()
            .filter_map(|id| {
                let item = self.entity(id)?;
                let mapped = item
                    .part("STYLED_ITEM")
                    .and_then(|fields| reference(fields.get(2)))
                    .unwrap_or(id);
                self.has(mapped, "MAPPED_ITEM").then_some(mapped)
            })
            .collect();
        let mapped = match maps.len() {
            0 => return Ok(scale),
            1 => *maps.first().expect("one mapped item"),
            _ => return Err("the draughting model has ambiguous geometry mappings".into()),
        };
        let frames = (|| {
            let [Param::Str(_), source, target] = self.entity(mapped)?.part("MAPPED_ITEM")? else {
                return None;
            };
            let [origin, geometry] = self
                .entity(reference(Some(source))?)?
                .part("REPRESENTATION_MAP")?
            else {
                return None;
            };
            let geometry = reference(Some(geometry))?;
            if !matches!(
                self.entity(geometry)?.type_name(),
                "SHAPE_REPRESENTATION"
                    | "ADVANCED_BREP_SHAPE_REPRESENTATION"
                    | "MANIFOLD_SURFACE_SHAPE_REPRESENTATION"
                    | "EDGE_BASED_WIREFRAME_SHAPE_REPRESENTATION"
            ) || self.representation_scale(geometry).ok()? != scale
            {
                return None;
            }
            Some((
                self.placement(reference(Some(origin))?)?,
                self.placement(reference(Some(target))?)?,
            ))
        })();
        let (source, target) =
            frames.ok_or("the draughting model has an unreadable geometry mapping")?;
        let same_direction = |a, b| is_zero3(&cross3(a, b)) && dot3(a, b) > Rational::zero();
        if source.origin != target.origin
            || !same_direction(&source.axis, &target.axis)
            || !same_direction(&source.reference, &target.reference)
        {
            return Err("the draughting model has a nonidentity geometry mapping".into());
        }
        Ok(scale)
    }

    fn trim_curve(&self, id: u64) -> Option<&'a Entity> {
        let mut curve = self.entity(id)?;
        for wrapper in ["SURFACE_CURVE", "SEAM_CURVE"] {
            if let Some(params) = curve.part(wrapper) {
                curve = self.entity(reference(params.get(1))?)?;
            }
        }
        Some(curve)
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
                let curve = self.trim_curve(reference(edge.get(3))?)?;
                if curve.part("LINE").is_some() {
                    edges.push(PlanarEdge::Segment(ends));
                } else {
                    // A face bounded by another curve is not checked.
                    let circle = curve.part("CIRCLE")?;
                    let placement = self.placement(reference(circle.get(1))?)?;
                    let radius = &number(circle.get(2))? * scale;
                    edges.push(PlanarEdge::Arc {
                        circle: Circle {
                            center: at(&placement.origin),
                            normal: along(&placement.axis),
                            radius_sq: Real::from_rational(&radius * &radius),
                            ref_dir: along(&placement.reference),
                        },
                        ends: (ends[0] != ends[1])
                            .then(|| if same_sense { ends } else { [ends[1], ends[0]] }),
                    });
                }
            }
        }
        (!edges.is_empty()).then(|| PlanarBoundary {
            plane: Plane {
                origin: at(&frame.origin),
                normal: along(&frame.axis),
                ref_dir: along(&frame.reference),
            },
            vertices,
            edges,
        })
    }

    /// The equivalent Unicode string a callout's validation property states.
    fn callout_text(&self, callout: u64) -> Result<Option<String>, String> {
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
                            return string(fields.get(1)).map(Some).ok_or_else(|| "the callout equivalent text has an unsupported STEP encoding or is missing".into());
                        }
                    }
                }
            }
        }
        Ok(None)
    }

    fn point_coordinates(&self, id: u64) -> Option<[Rational; 3]> {
        let [Param::Str(_), point] = self.entity(id)?.part("CARTESIAN_POINT")? else {
            return None;
        };
        coordinates(Some(point))
    }

    fn graphic_geometry(
        &self,
        occurrences: &[u64],
        scale: &Rational,
    ) -> Result<(String, Vec<Vec<Vec3>>), String> {
        let mut kind = String::new();
        let mut strokes = Vec::new();
        for id in occurrences {
            let occurrence = self
                .entity(*id)
                .ok_or("an annotation occurrence is missing")?;
            if !self.has(*id, "ANNOTATION_CURVE_OCCURRENCE") {
                return Err(format!(
                    "a {} has no supported strokes",
                    occurrence.type_name()
                ));
            }
            let fields = occurrence
                .part("STYLED_ITEM")
                .unwrap_or(occurrence.parts[0].1.as_slice());
            let item = reference(fields.get(2))
                .and_then(|id| self.entity(id))
                .ok_or("an annotation occurrence has no geometry")?;
            let mut curves = vec![item];
            if let Some(set) = item.part("GEOMETRIC_CURVE_SET") {
                if kind.is_empty() {
                    kind = string(set.first()).unwrap_or_default();
                }
                curves = references(set.get(1))
                    .into_iter()
                    .map(|id| self.entity(id).ok_or("an annotation curve is missing"))
                    .collect::<Result<_, _>>()?;
            }
            for curve in curves {
                let polyline = curve.part("POLYLINE").ok_or_else(|| {
                    format!(
                        "an annotation curve is a {}, which is not read",
                        curve.type_name()
                    )
                })?;
                let points = references(polyline.get(1))
                    .into_iter()
                    .map(|id| {
                        self.point_coordinates(id)
                            .map(|point| point.map(|value| &value * scale))
                            .ok_or("an annotation polyline has an unreadable point")
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                if points.len() < 2 {
                    return Err("an annotation polyline needs at least two points".into());
                }
                strokes.push(points);
            }
        }
        if strokes.is_empty() {
            return Err("the callout has no supported strokes".into());
        }
        Ok((kind, strokes))
    }

    fn annotation_plane(&self, id: u64, scale: &Rational) -> Option<GraphicPlane> {
        let plane = self.entity(id)?.part("ANNOTATION_PLANE")?;
        let carrier = self.entity(reference(plane.get(2))?)?;
        let placement = if let Some(surface) = carrier.part("PLANE") {
            reference(surface.get(1))?
        } else {
            reference(carrier.part("PLANAR_BOX")?.get(3))?
        };
        let frame = self.placement(placement)?;
        Some(GraphicPlane {
            origin_mm: frame.origin.map(|value| &value * scale),
            normal: frame.axis,
            reference: frame.reference,
        })
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
            .ok_or("the camera has no readable placement or its placement axes are parallel")?;
        let norm = dot3(&frame.axis, &frame.axis);
        let up = cross3(&frame.axis, &frame.reference).map(|value| &value / &norm);
        let volume = reference(params.get(2))
            .and_then(|volume| self.entity(volume))
            .and_then(|volume| volume.part("VIEW_VOLUME"))
            .ok_or("the camera has no view volume")?;
        let perspective = match volume.first() {
            Some(Param::Enum(kind)) if kind == "CENTRAL" => true,
            Some(Param::Enum(kind)) if kind == "PARALLEL" => false,
            _ => return Err("the camera's projection is neither central nor parallel".into()),
        };
        let projection = reference(volume.get(1))
            .and_then(|id| self.point_coordinates(id))
            .ok_or("the camera has no projection point")?;
        if !is_zero3(&projection) {
            return Err("the camera's projection point is outside the placement origin".into());
        }
        if [4, 6, 7]
            .iter()
            .any(|&index| !matches!(volume.get(index), Some(Param::Enum(flag)) if flag == "F"))
        {
            return Err("the camera has active or unreadable clipping".into());
        }
        let distance = number(volume.get(2))
            .filter(|value| *value > Rational::zero())
            .ok_or("the camera needs a positive view distance")?;
        let window = reference(volume.get(8))
            .and_then(|window| self.entity(window))
            .and_then(|window| window.part("PLANAR_BOX"))
            .ok_or("the camera has no view window")?;
        let size = |index: usize| {
            number(window.get(index))
                .filter(|value| *value > Rational::zero())
                .ok_or("the camera's window needs a positive size")
        };
        let dimensions = [size(1)?, size(2)?];
        let placement = reference(window.get(3))
            .and_then(|id| self.placement(id))
            .ok_or("the camera's window has no 3D placement")?;
        let negative_half = Rational::from_ratio(-1, 2).expect("a nonzero denominator");
        if placement.origin[0] != &dimensions[0] * &negative_half
            || placement.origin[1] != &dimensions[1] * &negative_half
        {
            return Err("the camera's window is off centre".into());
        }
        if placement.origin[2] != distance {
            return Err("the camera's window depth differs from its view distance".into());
        }
        let across = placement.reference;
        if !placement.axis[0].is_zero()
            || !placement.axis[1].is_zero()
            || placement.axis[2] <= Rational::zero()
            || across[0] <= Rational::zero()
            || !across[1].is_zero()
            || !across[2].is_zero()
        {
            return Err("the camera's window is rotated".into());
        }
        Ok(Camera {
            eye: frame.origin.map(|value| &value * scale),
            direction: frame.axis,
            up,
            distance: &distance * scale,
            window: dimensions.map(|value| &value * scale),
            perspective,
        })
    }

    fn graphic_callout(entity: &Entity) -> Option<&str> {
        const CALLOUTS: &[&str] = &[
            "ANGULAR_DIMENSION",
            "CURVE_DIMENSION",
            "DATUM_FEATURE_CALLOUT",
            "DATUM_TARGET_CALLOUT",
            "DIAMETER_DIMENSION",
            "DIMENSION_CALLOUT",
            "DIMENSION_CURVE_DIRECTED_CALLOUT",
            "DRAUGHTING_CALLOUT",
            "DRAUGHTING_ELEMENTS",
            "GEOMETRICAL_TOLERANCE_CALLOUT",
            "LEADER_DIRECTED_CALLOUT",
            "LEADER_DIRECTED_DIMENSION",
            "LINEAR_DIMENSION",
            "ORDINATE_DIMENSION",
            "PROJECTION_DIRECTED_CALLOUT",
            "RADIUS_DIMENSION",
            "STRUCTURED_DIMENSION_CALLOUT",
            "SURFACE_CONDITION_CALLOUT",
        ];
        entity
            .parts
            .iter()
            .rev()
            .find_map(|(kind, _)| CALLOUTS.contains(&kind.as_str()).then_some(kind.as_str()))
    }

    fn draughting_model(entity: &Entity) -> Option<&[Param]> {
        let fields = entity.part("DRAUGHTING_MODEL")?;
        let fields = if entity.parts.len() == 1 {
            fields
        } else {
            entity.part("REPRESENTATION")?
        };
        matches!(fields, [Param::Str(_), Param::List(_), Param::Ref(_)]).then_some(fields)
    }

    fn draughting_models_containing(&self, member: u64) -> impl Iterator<Item = u64> + '_ {
        self.sharing
            .get(&member)
            .into_iter()
            .flatten()
            .copied()
            .filter(move |id| {
                self.entity(*id)
                    .and_then(Self::draughting_model)
                    .is_some_and(|fields| references(fields.get(1)).contains(&member))
            })
    }

    fn graphic_faces(&self, id: u64) -> Result<Vec<String>, String> {
        let entity = self.entity(id).ok_or("its linked item is missing")?;
        let mut aspects = Vec::new();
        if let Some(fields) = entity.part("DIMENSIONAL_SIZE") {
            aspects.extend(reference(fields.first()));
        } else if let Some((_, a, b)) = self.relationship(id)
            && (entity.part("DIMENSIONAL_LOCATION").is_some()
                || entity.part("ANGULAR_LOCATION").is_some())
        {
            aspects.extend([a, b]);
        } else if let Some(fields) = entity.part("PROPERTY_DEFINITION") {
            aspects.extend(reference(fields.get(2)));
        } else if entity.type_name() == "DATUM" {
            aspects.extend(
                self.sharing
                    .get(&id)
                    .into_iter()
                    .flatten()
                    .filter_map(|link| self.relationship(*link))
                    .filter_map(|(_, feature, datum)| {
                        (datum == id && self.has(feature, "DATUM_FEATURE")).then_some(feature)
                    }),
            );
        } else if let Some(fields) = entity.part("GEOMETRIC_TOLERANCE").or_else(|| {
            Characteristic::of_entity(entity.type_name()).map(|_| entity.parts[0].1.as_slice())
        }) {
            let target = reference(fields.get(3)).ok_or("its linked tolerance names no feature")?;
            aspects.push(
                self.entity(target)
                    .and_then(|item| item.part("DIMENSIONAL_SIZE"))
                    .and_then(|fields| reference(fields.first()))
                    .unwrap_or(target),
            );
        } else if self.aspect_params(id).is_some() {
            aspects.push(id);
        }
        if aspects.is_empty() {
            return Err("its linked item names no face feature".into());
        }
        let mut faces = Vec::new();
        for aspect in aspects {
            faces.extend(self.features(aspect, 0)?.into_iter().flatten());
        }
        Ok(faces)
    }

    fn legacy_texture<'b>(
        &self,
        name: &str,
        text: Option<&str>,
        aspect: u64,
        items: &'b [(u64, Item<String>)],
    ) -> Result<Option<&'b str>, String> {
        let Some([Param::Str(aspect_name), _, owner, Param::Enum(flag)]) = self
            .entity(aspect)
            .filter(|entity| entity.parts.len() == 1)
            .and_then(|entity| entity.part("SHAPE_ASPECT"))
        else {
            return Ok(None);
        };
        let Some(owner) =
            reference(Some(owner)).filter(|id| self.has(*id, "PRODUCT_DEFINITION_SHAPE"))
        else {
            return Ok(None);
        };
        if flag != "T" || decode(aspect_name).as_deref() != Some(name) {
            return Ok(None);
        }
        let Some((prefix, requirement)) = text.and_then(|text| text.split_once(r"\w\w\w\w")) else {
            return Ok(None);
        };
        let process = match prefix {
            "SRFAPAN" => Process::Any,
            "SRFMRRN" => Process::Removal,
            "SRFNMRN" => Process::NoRemoval,
            _ => return Ok(None),
        };
        let Some((parameter, value)) = requirement.split_once(' ') else {
            return Ok(None);
        };
        let Ok(value) = Rational::parse(value) else {
            return Ok(None);
        };
        let Some(mut faces) = self.used_faces(aspect)? else {
            return Ok(None);
        };
        faces.sort();
        let mut matched = items.iter().filter(|(property, item)| {
            let ItemKind::Texture(texture) = &item.kind else {
                return false;
            };
            let semantic_aspect = self
                .entity(*property)
                .and_then(|e| e.part("PROPERTY_DEFINITION"))
                .and_then(|fields| reference(fields.get(2)));
            self.id_attribute(*property).as_deref() == Some(name)
                && semantic_aspect.is_some_and(|id| {
                    id != aspect
                        && self
                            .aspect_params(id)
                            .and_then(|fields| reference(fields.get(2)))
                            == Some(owner)
                })
                && texture.parameter == parameter
                && texture.value == value
                && texture.process == process
                && texture
                    .feature
                    .iter()
                    .cloned()
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .eq(faces.iter().cloned())
        });
        let first = matched.next();
        if matched.next().is_some() {
            return Err("its legacy surface texture association matches several finishes".into());
        }
        Ok(first.map(|(_, item)| item.id.as_str()))
    }

    /// The saved views and graphic callouts of the file. A callout whose
    /// association names a read item presents that item and is not listed
    /// again as graphics.
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
        for (datum, item) in items {
            if matches!(item.kind, ItemKind::Datum(_)) {
                for link in self.sharing.get(datum).into_iter().flatten() {
                    if let Some((_, feature, related)) = self.relationship(*link)
                        && related == *datum
                        && self.has(feature, "DATUM_FEATURE")
                    {
                        by_entity.insert(feature, item.id.as_str());
                    }
                }
            }
        }
        for (lower, upper) in lowers {
            if let Some(name) = by_entity.get(upper).copied() {
                by_entity.insert(*lower, name);
            }
        }
        let mut presented = BTreeMap::new();
        let mut graphics = BTreeMap::new();
        for (&callout, entity) in self.g {
            let Some(kind) = Self::graphic_callout(entity) else {
                continue;
            };
            let params = &entity.parts[0].1;
            if kind != "DRAUGHTING_CALLOUT"
                || entity.parts.len() != 1
                || !matches!(params.as_slice(), [Param::Str(_), Param::List(_)])
            {
                skipped.push(SkippedPmi {
                    entity: format!("#{callout}"),
                    kind: kind.into(),
                    reason: "this draughting callout form is not supported".into(),
                });
                continue;
            }
            let name = string(params.first()).unwrap_or_default();
            let occurrences = references(params.get(1));
            let subjects: Vec<_> = std::iter::once(callout)
                .chain(occurrences.iter().copied())
                .collect();
            let mut linked = Vec::new();
            let mut presents = BTreeSet::new();
            let mut legacy_aspects = BTreeSet::new();
            let mut faces = Vec::new();
            let mut models = BTreeSet::new();
            let mut planes = BTreeSet::new();
            for &subject in &subjects {
                for (_, association) in
                    self.sharing_at(subject, "DRAUGHTING_MODEL_ITEM_ASSOCIATION", 4)
                {
                    models.extend(reference(association.get(3)));
                    let Some(definition) = reference(association.get(2)) else {
                        continue;
                    };
                    if let Some(name) = by_entity.get(&definition) {
                        linked.push((*name).to_owned());
                        continue;
                    }
                    match self.graphic_faces(definition) {
                        Ok(attached) => faces.extend(attached),
                        Err(reason) => skipped.push(SkippedPmi {
                            entity: format!("#{callout}"),
                            kind: entity.type_name().into(),
                            reason: format!(
                                "cannot resolve graphic attachment #{definition}: {reason}"
                            ),
                        }),
                    }
                    if self.aspect_params(definition).is_none()
                        || self.has(definition, "DATUM_FEATURE")
                        || self.has(definition, "PLACED_DATUM_TARGET_FEATURE")
                    {
                        presents.insert(definition);
                    } else if !name.is_empty()
                        && string(association.first()).as_deref() == Some(&name)
                    {
                        legacy_aspects.insert(definition);
                    }
                }
                models.extend(self.draughting_models_containing(subject));
                planes.extend(
                    self.sharing_at(subject, "ANNOTATION_PLANE", 3)
                        .into_iter()
                        .map(|(id, _)| id),
                );
            }
            let mut equivalent = None;
            for &id in &subjects {
                match self.callout_text(id) {
                    Ok(Some(text)) => {
                        equivalent.get_or_insert(text);
                    }
                    Ok(None) => {}
                    Err(reason) => skipped.push(SkippedPmi {
                        entity: format!("#{callout}"),
                        kind: entity.type_name().into(),
                        reason,
                    }),
                }
            }
            for &plane in &planes {
                models.extend(self.draughting_models_containing(plane));
            }
            let geometry = (|| -> Result<_, String> {
                let mut scale = None;
                for model in models {
                    let unit = self.presentation_scale(model)?;
                    if scale.as_ref().is_some_and(|known| known != &unit) {
                        return Err("the callout has conflicting length units".into());
                    }
                    scale = Some(unit);
                }
                let scale = scale.ok_or("the callout has no known length unit")?;
                let (kind, polylines) = self.graphic_geometry(&occurrences, &scale)?;
                let planes = planes
                    .into_iter()
                    .map(|id| {
                        self.annotation_plane(id, &scale)
                            .ok_or("the callout's annotation plane is unreadable")
                    })
                    .collect::<Result<_, _>>()?;
                Ok((kind, polylines, planes))
            })();
            let (kind, polylines_mm, planes) = match geometry {
                Ok(geometry) => geometry,
                Err(reason) => {
                    skipped.push(SkippedPmi {
                        entity: format!("#{callout}"),
                        kind: entity.type_name().into(),
                        reason,
                    });
                    continue;
                }
            };
            if !linked.is_empty() {
                presented.insert(callout, linked);
                if presents.is_empty() {
                    continue;
                }
            }
            if kind == "surface roughness" && presents.is_empty() {
                for aspect in legacy_aspects {
                    match self.legacy_texture(&name, equivalent.as_deref(), aspect, items) {
                        Ok(Some(item)) => {
                            presented.entry(callout).or_default().push(item.to_owned())
                        }
                        Ok(None) => {}
                        Err(reason) => skipped.push(SkippedPmi {
                            entity: format!("#{callout}"),
                            kind: entity.type_name().into(),
                            reason,
                        }),
                    }
                }
                if presented.contains_key(&callout) {
                    continue;
                }
            }
            faces.sort();
            faces.dedup();
            graphics.insert(
                callout,
                GraphicCallout {
                    entity: format!("#{callout}"),
                    name,
                    kind,
                    text: equivalent,
                    faces: (!faces.is_empty()).then_some(faces),
                    presents: presents.into_iter().map(|id| format!("#{id}")).collect(),
                    polylines_mm,
                    planes,
                    views: Vec::new(),
                },
            );
        }
        let mut views: Vec<SavedView> = Vec::new();
        let mut view_names = BTreeSet::new();
        for (&model, entity) in self.g {
            if entity.part("DRAUGHTING_MODEL").is_none() {
                continue;
            }
            let Some(params) = Self::draughting_model(entity) else {
                skipped.push(SkippedPmi {
                    entity: format!("#{model}"),
                    kind: "DRAUGHTING_MODEL".into(),
                    reason: "the draughting model has no readable representation".into(),
                });
                continue;
            };
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
                if Self::graphic_callout(entity).is_some() {
                    shown.push(*member);
                } else if let Some(plane) = entity.part("ANNOTATION_PLANE") {
                    for element in references(plane.get(3)) {
                        if self
                            .entity(element)
                            .and_then(Self::graphic_callout)
                            .is_some()
                        {
                            shown.push(element);
                        } else {
                            shown.extend(
                                self.sharing
                                    .get(&element)
                                    .into_iter()
                                    .flatten()
                                    .filter(|id| {
                                        self.entity(**id).and_then(Self::graphic_callout).is_some()
                                    })
                                    .copied(),
                            );
                        }
                    }
                }
            }
            let scale = match self.presentation_scale(model) {
                Ok(scale) => scale,
                Err(reason) => {
                    skipped.push(SkippedPmi {
                        entity: format!("#{model}"),
                        kind: "DRAUGHTING_MODEL".into(),
                        reason,
                    });
                    continue;
                }
            };
            for camera in cameras {
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
                let base = self
                    .entity(camera)
                    .and_then(|entity| string(entity.parts[0].1.first()))
                    .filter(|name| !name.trim().is_empty())
                    .or_else(|| string(params.first()).filter(|name| !name.trim().is_empty()))
                    .unwrap_or_else(|| format!("View #{camera}"));
                let mut name = base.clone();
                for index in 2.. {
                    if view_names.insert(name.clone()) {
                        break;
                    }
                    name = format!("{base} ({index})");
                }
                let mut view_items: Vec<String> = Vec::new();
                for callout in &shown {
                    if let Some(items) = presented.get(callout) {
                        for item in items {
                            if !view_items.contains(item) {
                                view_items.push(item.clone());
                            }
                        }
                    }
                    if let Some(graphic) = graphics.get_mut(callout)
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
        let mut items = Vec::new();
        for (_, link) in self.sharing_at(id, "PROPERTY_DEFINITION_REPRESENTATION", 0) {
            let representation = reference(link.get(1))
                .and_then(|id| self.entity(id))
                .ok_or("its text is missing")?;
            items.extend(references(representation.parts[0].1.get(1)));
        }
        let item = match items.as_slice() {
            [] => return Err("it holds no text".into()),
            [item] => *item,
            _ => return Err("its text has no defined order".into()),
        };
        let fields = self
            .entity(item)
            .and_then(|entity| entity.part("DESCRIPTIVE_REPRESENTATION_ITEM"))
            .ok_or("its text holds an item that is not a string")?;
        let name = string(fields.first()).ok_or("its text has no variable name")?;
        let text = string(fields.get(1))
            .ok_or("its text has an unsupported STEP encoding or is missing")?;
        let feature = self.one_feature(target)?;
        let kind = match name.as_str() {
            "" => ItemKind::FaceNote { text, feature },
            KNURL => ItemKind::Knurl {
                designation: text,
                feature,
            },
            _ => return Err(format!("a text variable `{name}` is not read")),
        };
        Ok(Found {
            order: id,
            datum_sources: Vec::new(),
            name: self.id_attribute(id),
            kind,
        })
    }

    /// A placed datum target (recommended practices, section 6.6): its datum,
    /// its one face, its `orientation` and sizes. The placement's axis must
    /// point out of the face, the point must lie on the face's plane, and the
    /// whole target must lie on the face within its trim.
    fn target(&self, id: u64, params: &[Param]) -> Result<Found, String> {
        let [
            Param::Str(_),
            Param::Str(raw),
            owner,
            Param::Enum(flag),
            number,
        ] = params
        else {
            return Err("its placed target has invalid attributes".into());
        };
        if !matches!(flag.as_str(), "T" | "F" | "U") {
            return Err("its placed target has an invalid definitional flag".into());
        }
        let description = decode(raw).ok_or("its target description is unreadable")?;
        let owner = reference(Some(owner))
            .filter(|id| self.has(*id, "PRODUCT_DEFINITION_SHAPE"))
            .ok_or("its target has no owning part")?;
        let number = string(Some(number))
            .and_then(|text| text.trim().parse::<u32>().ok())
            .filter(|number| *number > 0)
            .ok_or("its target number is not a positive whole number")?;
        let mut datums = Vec::new();
        for user in self.sharing.get(&id).map_or(&[][..], Vec::as_slice) {
            if let Some(fields) = self
                .entity(*user)
                .and_then(|e| e.part("SHAPE_ASPECT_RELATIONSHIP"))
                && (string(fields.get(1)).as_deref() == Some(catallax_step::pmi::DATUM_TARGET)
                    || fields.iter().any(|field| {
                        references(Some(field))
                            .into_iter()
                            .any(|id| self.has(id, "DATUM"))
                    }))
            {
                let [Param::Str(_), description, from, to] = fields else {
                    return Err("its datum relationship has invalid attributes".into());
                };
                if reference(Some(from)) != Some(id)
                    || string(Some(description)).as_deref()
                        != Some(catallax_step::pmi::DATUM_TARGET)
                {
                    return Err("its datum relationship is unreadable".into());
                }
                let related = reference(Some(to))
                    .filter(|id| self.has(*id, "DATUM"))
                    .ok_or("its datum relationship has no datum")?;
                datums.push(related);
            }
        }
        let [datum] = datums.as_slice() else {
            return Err("it does not belong to exactly one datum".into());
        };
        if self
            .entity(*datum)
            .and_then(|datum| datum.part("DATUM"))
            .and_then(|fields| reference(fields.get(2)))
            != Some(owner)
        {
            return Err("its datum belongs to another part".into());
        }
        let letter = self.datum_letter(*datum).ok_or("its datum has no letter")?;
        let links = self.sharing_at(id, "FEATURE_FOR_DATUM_TARGET_RELATIONSHIP", 3);
        let [(_, link)] = links.as_slice() else {
            return Err(if links.is_empty() {
                "it is tied to no face".into()
            } else {
                "it is tied to more than one face feature".into()
            });
        };
        let [Param::Str(_), Param::Str(_) | Param::Unset, aspect, placed] = *link else {
            return Err("its face relationship has invalid attributes".into());
        };
        if reference(Some(placed)) != Some(id) {
            return Err("its face relationship does not name one target".into());
        }
        let aspect = reference(Some(aspect)).ok_or("its face feature is missing")?;
        if self.aspect_params(aspect).and_then(|p| reference(p.get(2))) != Some(owner) {
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
        let context =
            reference(representation.get(2)).ok_or("its parameters have no unit context")?;
        let scale = self.context_millimetres(context)?;
        let (mut orientation, mut length, mut width, mut diameter) = (None, None, None, None);
        let [Param::Str(_), parameters, _] = representation else {
            return Err("its parameter representation has invalid attributes".into());
        };
        for item in reference_list(Some(parameters))
            .ok_or("its parameter list is not a flat set of distinct items")?
        {
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
        let direction = orientation.reference;
        let kind = match (description.as_str(), length, width, diameter) {
            ("point", None, None, None) => TargetKind::Point,
            ("line", Some(length_mm), None, None) => TargetKind::Line {
                length_mm,
                direction: direction.clone(),
            },
            ("circle", None, None, Some(diameter_mm)) => {
                TargetKind::Area(TargetArea::Circular { diameter_mm })
            }
            ("rectangle", Some(length_mm), Some(width_mm), None) => {
                TargetKind::Area(TargetArea::Rectangular {
                    length_mm,
                    width_mm,
                    direction,
                })
            }
            ("point" | "line" | "circle" | "rectangle", ..) => {
                return Err(format!("its sizes do not fit a {description} target"));
            }
            (other, ..) => return Err(format!("a `{other}` datum target is not read")),
        };
        let (_, outward) = self
            .outward_plane(face)
            .ok_or("a datum target is read only on a planar face")?;
        if !is_zero3(&cross3(&orientation.axis, &outward))
            || dot3(&orientation.axis, &outward) <= Rational::zero()
        {
            return Err("its orientation does not point out of its face".into());
        }
        let target = DatumTarget {
            datum: letter.clone(),
            number,
            face: face.clone(),
            point_mm: orientation.origin.clone().map(|value| &value * &scale),
            kind,
        };
        let face_id = face[1..].parse().expect("a used face names its entity");
        let mut face_scale = None;
        for kind in [
            "GEOMETRIC_ITEM_SPECIFIC_USAGE",
            "ITEM_IDENTIFIED_REPRESENTATION_USAGE",
        ] {
            for (_, usage) in self.sharing_at(aspect, kind, 2) {
                if !references(usage.get(4)).contains(&face_id) {
                    continue;
                }
                let [
                    Param::Str(_),
                    Param::Str(_) | Param::Unset,
                    definition,
                    representation,
                    identified,
                ] = usage
                else {
                    return Err("its face usage has invalid attributes".into());
                };
                let one_face = match identified {
                    Param::Ref(id) => *id == face_id,
                    Param::Typed(kind, fields)
                        if kind == "SET_REPRESENTATION_ITEM" && fields.len() == 1 =>
                    {
                        reference_list(fields.first()).as_deref() == Some(&[face_id])
                    }
                    _ => false,
                };
                if !one_face || reference(Some(definition)) != Some(aspect) {
                    return Err("its face usage must name one face and one feature".into());
                }
                let representation = reference(Some(representation))
                    .filter(|rep| {
                        self.sharing_at(*rep, "SHAPE_DEFINITION_REPRESENTATION", 1)
                            .iter()
                            .any(|(_, link)| reference(link.first()) == Some(owner))
                            && self.representation_contains_face(*rep, face_id)
                    })
                    .ok_or(
                        "its face usage does not identify a face in the owning representation",
                    )?;
                let scale = self.representation_scale(representation)?;
                if face_scale.as_ref().is_some_and(|known| known != &scale) {
                    return Err("its face has conflicting length units".into());
                }
                face_scale = Some(scale);
            }
        }
        let face_scale = face_scale.ok_or("its face has no known length unit")?;
        let boundary = self.face_boundary(face, &face_scale).ok_or(
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
            datum_sources: vec![*datum],
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
            datum_sources: Vec::new(),
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
    fn assembly_presentation(&self, shapes: &BTreeMap<u64, u64>) -> BTreeMap<u64, Vec<SkippedPmi>> {
        let mut out = BTreeMap::<_, Vec<_>>::new();
        let mut representation_owners = BTreeMap::<_, BTreeSet<_>>::new();
        let mut models = BTreeMap::<_, BTreeSet<_>>::new();
        for (&id, entity) in self.g {
            if let Some(fields) = entity.part("SHAPE_DEFINITION_REPRESENTATION")
                && let (Some(shape), Some(rep)) =
                    (reference(fields.first()), reference(fields.get(1)))
                && shapes.contains_key(&shape)
            {
                representation_owners.entry(rep).or_default().insert(shape);
            }
            if let Some(fields) = entity.part("DRAUGHTING_MODEL_ITEM_ASSOCIATION")
                && let (Some(definition), Some(model)) =
                    (reference(fields.get(2)), reference(fields.get(3)))
            {
                models
                    .entry(model)
                    .or_default()
                    .extend(self.owner_shapes(definition));
            }
            if entity.part("DRAUGHTING_MODEL").is_some() {
                models.entry(id).or_default();
            }
        }
        for (model, mut owners) in models {
            let Some(entity) = self
                .entity(model)
                .filter(|entity| entity.part("DRAUGHTING_MODEL").is_some())
            else {
                continue;
            };
            let fields = Self::draughting_model(entity);
            let members = references(fields.and_then(|fields| fields.get(1)));
            for &member in &members {
                let mapped = self
                    .entity(member)
                    .and_then(|e| e.part("MAPPED_ITEM"))
                    .and_then(|fields| reference(fields.get(1)))
                    .and_then(|id| self.entity(id)?.part("REPRESENTATION_MAP"))
                    .and_then(|fields| reference(fields.get(1)));
                if let Some(shapes) = mapped.and_then(|rep| representation_owners.get(&rep)) {
                    owners.extend(shapes);
                }
            }
            let has_view = members.iter().any(|id| {
                self.entity(*id)
                    .is_some_and(|e| e.type_name().starts_with("CAMERA_MODEL"))
            });
            let has_graphics = members.iter().any(|id| {
                self.entity(*id).is_some_and(|e| {
                    e.part("ANNOTATION_PLANE").is_some() || Self::graphic_callout(e).is_some()
                })
            });
            if fields.is_some() && !(has_view || has_graphics) {
                continue;
            }
            let diagnostic = SkippedPmi {
                entity: format!("#{model}"),
                kind: "DRAUGHTING_MODEL".into(),
                reason: if fields.is_none() {
                    "the assembly draughting model has no readable representation"
                } else if has_view {
                    "assembly saved-view presentation is not read"
                } else {
                    "assembly graphical presentation is not read"
                }
                .into(),
            };
            owners.retain(|owner| shapes.contains_key(owner));
            if owners.is_empty() {
                owners.extend(shapes.keys().copied());
            }
            for owner in owners {
                out.entry(owner).or_default().push(diagnostic.clone());
            }
        }
        out
    }

    fn part_notes(&self, scope: Option<(u64, u64)>) -> PartNotes {
        let mut notes = PartNotes::default();
        let product = scope.map(|(product, _)| product).or_else(|| {
            let mut products = self
                .g
                .keys()
                .copied()
                .filter(|&id| self.has(id, "PRODUCT_DEFINITION"));
            let product = products.next()?;
            products.next().is_none().then_some(product)
        });
        if product.is_some_and(|id| self.has_asme_standard(id)) {
            notes.standard = Some(catallax_step::pmi::DimensioningStandard::AsmeY14_5_2018);
        }
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
        let mut owners = reader.owner_shapes(item.id);
        if matches!(item.kind, PmiKind::Text)
            && let Some(definition) = item
                .entity
                .part("PROPERTY_DEFINITION")
                .and_then(|fields| reference(fields.get(2)))
                .filter(|definition| reader.has(*definition, "PRODUCT_DEFINITION"))
        {
            owners.extend(
                shapes
                    .iter()
                    .filter_map(|(&shape, &leaf)| (leaf == definition).then_some(shape)),
            );
        }
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
    let mut presentation = reader.assembly_presentation(shapes);
    Ok(shapes
        .iter()
        .map(|(&shape, &definition)| {
            let mut pmi = read_items(
                &reader,
                items.remove(&shape).unwrap_or_default(),
                Some((definition, shape)),
            );
            pmi.skipped
                .extend(presentation.remove(&shape).unwrap_or_default());
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
            PmiKind::Tolerance(characteristic) => match reader.composite_link(id) {
                Ok(Some((upper, lower))) if id == lower => {
                    lowers.push((id, entity.type_name(), upper));
                    continue;
                }
                Ok(link) => reader.tolerance(id, characteristic, link.map(|(_, lower)| lower)),
                Err(reason) => Err(reason),
            },
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
    let mut accepted_datums = BTreeSet::new();
    let mut rejected_targets = BTreeSet::new();
    let targets: BTreeMap<_, _> = later
        .iter()
        .filter_map(|(_, item)| {
            if let ItemKind::DatumTarget(target) = &item.kind {
                Some((item.order, (&item.kind, target)))
            } else {
                None
            }
        })
        .collect();
    for (kind, item) in &mut datums {
        let ItemKind::Datum(datum) = &item.kind else {
            unreachable!("datums were collected separately");
        };
        let name = item.name.take();
        let letter = datum.letter.as_str();
        if let Err(error) = datum.validate_letter() {
            skip(item.order, kind, error.to_string());
            continue;
        }
        let mut local = PartValidator::default();
        local
            .accept_datum(datum)
            .expect("the first datum cannot conflict");
        let mut bad_target = None;
        let Ok(linked_targets) = reader.placed_targets(item.order) else {
            unreachable!("the datum's target relationships were read");
        };
        for target in linked_targets {
            match targets.get(&target) {
                Some((kind, target_item)) => {
                    if let Err(error) = local
                        .check_kind(kind)
                        .and_then(|()| local.accept_target(target_item))
                    {
                        bad_target.get_or_insert(target);
                        rejected_targets.insert(target);
                        skip(target, "PLACED_DATUM_TARGET_FEATURE", error.to_string());
                    }
                }
                None => {
                    bad_target.get_or_insert(target);
                }
            }
        }
        if let Some(target) = bad_target {
            skip(
                item.order,
                kind,
                format!("its datum target #{target} was not read"),
            );
            continue;
        }
        if let Err(error) = validator.accept_datum(datum) {
            let reason = match error {
                ItemError::DuplicateDatumFeature if conflicting_letters.insert(letter) => {
                    "its feature already carries a datum".into()
                }
                ItemError::DuplicateDatumFeature => ItemError::DuplicateDatum(letter).to_string(),
                error => error.to_string(),
            };
            skip(item.order, kind, reason);
            continue;
        }
        accepted_datums.insert(item.order);
        item.name = Some(validator.assign_name(name.unwrap_or_else(|| format!("Datum {letter}"))));
    }
    drop(targets);
    for (kind, item) in &mut later {
        let name = item.name.take();
        if rejected_targets.contains(&item.order) {
            continue;
        }
        if let Some(source) = item
            .datum_sources
            .iter()
            .find(|source| !accepted_datums.contains(source))
        {
            let letter = reader.datum_letter(*source).unwrap_or_default();
            let reason = if conflicting_letters.contains(letter.as_str())
                && !validator.contains_datum(&letter)
            {
                ItemError::UndefinedDatum(&letter).to_string()
            } else {
                format!("it refers to datum {letter}, which was not read")
            };
            skip(item.order, kind, reason);
            continue;
        }
        if let Err(error) = validator.check_kind(&item.kind) {
            skip(item.order, kind, error.to_string());
            continue;
        }
        let base = name.unwrap_or_else(|| format!("#{}", item.order));
        item.name = Some(validator.assign_name(base));
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
    let sources = items
        .iter()
        .map(|(entity, item)| (item.id.clone(), format!("#{entity}")))
        .collect();
    StepPmi {
        sources,
        pmi: Pmi {
            items: items.into_iter().map(|(_, item)| item).collect(),
            part: reader.part_notes(scope),
            views,
        },
        skipped,
        graphics,
    }
}

/// One edit of a source file's PMI.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum PmiEdit {
    /// Gives a geometric tolerance a new zone value in millimetres.
    ToleranceValue {
        /// The item's source entity (such as `#412`) from [`StepPmi::sources`].
        source: String,
        /// The new value.
        value: Rational,
    },
}

fn patch_refusal(what: impl std::fmt::Display) -> Refusal {
    Refusal::new("step.pmi-patch", what.to_string())
}

/// A positive millimetre value in a file's length unit, as a STEP real.
fn file_real(value: &Rational, scale: &Rational) -> Option<String> {
    let text = (value / scale).to_decimal_exact()?;
    Some(if text.contains('.') {
        text
    } else {
        format!("{text}.")
    })
}

/// Edits the PMI of a STEP file in its own records, keeping every other
/// record as it was written: the assembly, colours, other annotations and
/// records this reader does not know. No edits give back the same bytes.
///
/// A tolerance value is a new length measure under a fresh number. The patch
/// refuses when a dependent record could retain the old meaning or when the
/// edited file fails read-back.
pub fn patch_step_pmi(bytes: &[u8], edits: &[PmiEdit]) -> Result<Vec<u8>, Refusal> {
    use catallax_step::pmi::patch::{Patch, Records};
    if edits.is_empty() {
        return Ok(bytes.to_vec());
    }
    let records = Records::parse(bytes)?;
    let graph = parse(bytes).map_err(|refusal| patch_refusal(format!("{refusal:?}")))?;
    let reader = Reader::new(&graph);
    let mut selected = BTreeMap::new();
    let mut edits_by_id = Vec::with_capacity(edits.len());
    let mut edited = BTreeSet::new();
    for PmiEdit::ToleranceValue { source, value } in edits {
        let root: u64 = source
            .strip_prefix('#')
            .and_then(|id| id.parse().ok())
            .ok_or_else(|| patch_refusal(format!("{source}: expected a source entity number")))?;
        if !edited.insert(root) {
            return Err(patch_refusal(format!("{source}: it is edited twice")));
        }
        let owners = reader.owner_shapes(root);
        let owner = owners
            .first()
            .copied()
            .filter(|_| owners.len() == 1)
            .ok_or_else(|| patch_refusal(format!("{source}: no item on one part was read")))?;
        let product = reader
            .entity(owner)
            .and_then(|shape| shape.part("PRODUCT_DEFINITION_SHAPE"))
            .and_then(|shape| reference(shape.get(2)))
            .filter(|id| reader.has(*id, "PRODUCT_DEFINITION"))
            .ok_or_else(|| patch_refusal(format!("{source}: no item on one part was read")))?;
        selected.insert(owner, product);
        edits_by_id.push((root, owner, source, value));
    }
    let scoped = |reader: &Reader<'_>| {
        let mut candidates = BTreeMap::<_, Vec<_>>::new();
        for item in reader
            .g
            .iter()
            .filter_map(|(&id, entity)| candidate(id, entity))
        {
            let owners = reader.owner_shapes(item.id);
            if let Some(owner) = owners
                .first()
                .copied()
                .filter(|owner| owners.len() == 1 && selected.contains_key(owner))
            {
                candidates.entry(owner).or_default().push(item);
            }
        }
        selected
            .iter()
            .map(|(&owner, &product)| {
                (
                    owner,
                    read_items(
                        reader,
                        candidates.remove(&owner).unwrap_or_default(),
                        Some((product, owner)),
                    ),
                )
            })
            .collect::<BTreeMap<_, _>>()
    };
    let mut patch = Patch::new(&records);
    let mut expected = scoped(&reader);
    for (root, owner, source, value) in edits_by_id {
        let part = expected.get_mut(&owner).expect("selected part was read");
        let (item, _) = part
            .sources
            .iter()
            .find(|(_, id)| *id == source)
            .ok_or_else(|| patch_refusal(format!("{source}: no such item was read")))?;
        let changed = part
            .pmi
            .items
            .iter_mut()
            .find(|found| &found.id == item)
            .expect("read parts have an item for each source");
        let ItemKind::Tolerance(tolerance) = &mut changed.kind else {
            return Err(patch_refusal(format!(
                "{item}: only a geometric tolerance's value is patched"
            )));
        };
        if *value <= Rational::zero() {
            return Err(patch_refusal(format!(
                "{item}: the value must be more than zero"
            )));
        }
        let refuse_dependency = |id: u64| {
            let kind = graph.get(&id).map_or("unknown", Entity::type_name);
            patch_refusal(format!(
                "{item}: #{id} ({kind}) may still depend on the old value"
            ))
        };
        for &user in records.users(root) {
            let entity = graph.get(&user).ok_or_else(|| refuse_dependency(user))?;
            let fields = entity.parts[0].1.as_slice();
            match (entity.parts.len(), entity.type_name()) {
                (1, "ID_ATTRIBUTE") if reference(fields.get(1)) == Some(root) => {}
                (1, "TOLERANCE_ZONE") if references(fields.get(4)) == [root] => {
                    if let Some(&dependent) = records.users(user).first() {
                        return Err(refuse_dependency(dependent));
                    }
                }
                _ => return Err(refuse_dependency(user)),
            }
        }
        let entity = &graph[&root];
        let base = entity
            .part("GEOMETRIC_TOLERANCE")
            .unwrap_or(entity.parts[0].1.as_slice());
        let magnitude = reference(base.get(2)).expect("a read tolerance has a value");
        if let Some(&user) = records.users(magnitude).iter().find(|&&user| user != root) {
            return Err(refuse_dependency(user));
        }
        let measure = &graph[&magnitude];
        let (true, [select, Param::Ref(unit)]) = (
            measure.parts.len() == 1 && measure.type_name() == "LENGTH_MEASURE_WITH_UNIT",
            measure.parts[0].1.as_slice(),
        ) else {
            return Err(patch_refusal(format!(
                "{item}: its value #{magnitude} is not a plain length measure"
            )));
        };
        let (kind, _) = length_measure(select).ok_or_else(|| {
            patch_refusal(format!("{item}: its value #{magnitude} is not a length"))
        })?;
        let scale = unit_millimetres(&graph, *unit, 0)
            .map_err(|_| patch_refusal(format!("{item}: its value's unit is not known")))?;
        let written = file_real(value, &scale).ok_or_else(|| {
            patch_refusal(format!(
                "{item}: the value is not a decimal in the file's unit"
            ))
        })?;
        let fresh = patch.add(format!(
            "LENGTH_MEASURE_WITH_UNIT({kind}({written}),#{unit})"
        ))?;
        patch.rewrite(root, records.replace_reference(root, magnitude, fresh)?);
        patch.remove(magnitude);
        tolerance.value = value.clone();
    }
    let out = patch.apply(&records);
    let graph = parse(&out).map_err(|refusal| patch_refusal(format!("{refusal:?}")))?;
    if scoped(&Reader::new(&graph)) != expected {
        return Err(patch_refusal(
            "the patched file does not read back as the edited PMI",
        ));
    }
    let patched = Records::parse(&out)?;
    let touched = patch.touched();
    for record in records.all() {
        if !touched.contains(&record.id) && patched.text(record.id) != records.text(record.id) {
            return Err(patch_refusal(format!(
                "#{}: an untouched record changed",
                record.id
            )));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::decode;

    #[test]
    fn strings_decode_their_control_directives() {
        for (encoded, decoded) in [
            ("plain", "plain"),
            (r"a\\b", r"a\b"),
            (r"\X2\00D8\X0\8 H7", "Ø8 H7"),
            (r"\X\D8", "Ø"),
            (r"\S\X", "Ø"),
            (r"\PA\\S\X", "Ø"),
            (r"\X4\0001F600\X0\", "😀"),
            (r"\X2\0041D83DDE000042\X0\", "A😀B"),
        ] {
            assert_eq!(decode(encoded).as_deref(), Some(decoded));
        }
        assert_eq!(decode(r"\PE\\S\*\S\U\S\b"), None);
        assert_eq!(decode(r"\X2\D83D0041\X0\"), None);
        assert_eq!(decode(r"\X2\DE00\X0\"), None);
        assert_eq!(decode(r"\X2\D83D\X0\"), None);
    }
}
