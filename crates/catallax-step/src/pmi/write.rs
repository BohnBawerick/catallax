//! AP242 entities for the manufacturing information of one part.
//!
//! Each feature becomes one `SHAPE_ASPECT` tied to its written faces, shared
//! by every item on that feature. A datum's feature is a `DATUM_FEATURE`. A
//! pattern is a `COMPOSITE_GROUP_SHAPE_ASPECT` of its features. Sizes and
//! distances are `DIMENSIONAL_SIZE` and `DIMENSIONAL_LOCATION` with a
//! nominal value and, where given, limits. Geometric tolerances refer to one
//! shared `DATUM_SYSTEM` per datum sequence. A surface texture is a
//! `SURFACE_TEXTURE_REPRESENTATION` on the feature. Part notes are user defined attributes of the product
//! definition.

use std::collections::BTreeMap;

use catallax_numbers::{Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Model};

use super::{
    Anchor, At, Feature, GeometricTolerance, ItemKind, Limits, MaterialCondition, Pmi, SizeKind,
    SurfaceTexture,
};
use crate::graph::{Field, Graph, Id, boolean, refs, string, text};
use crate::number::decimal;

/// The fixed entities of the file header the information hangs on.
pub(crate) const PRODUCT_DEFINITION: &str = "#7";
pub(crate) const PRODUCT_DEFINITION_SHAPE: &str = "#8";
pub(crate) const MILLIMETRE: &str = "#10";
pub(crate) const CONTEXT: &str = "#14";
pub(crate) const REPRESENTATION: &str = "#16";

fn quoted(value: &str) -> Field {
    text(string(value))
}

fn header(reference: &'static str) -> Field {
    text(reference)
}

/// A terminating decimal as people write it, such as `0.05`.
pub(crate) fn plain(value: &Rational, item: &str) -> Result<String, Refusal> {
    value.to_decimal_exact().ok_or_else(|| {
        Refusal::new(
            "step.pmi-invalid",
            format!("{item}: a value is not a decimal number"),
        )
    })
}

fn real(value: &Rational, item: &str) -> Result<String, Refusal> {
    Ok(decimal(&plain(value, item)?))
}

struct Emitter<'a> {
    graph: &'a mut Graph,
    model: &'a Model,
    faces: &'a BTreeMap<FaceId, Vec<Id>>,
    roots: Vec<Id>,
    aspects: BTreeMap<Vec<FaceId>, Id>,
    groups: BTreeMap<Vec<Vec<FaceId>>, Id>,
    centres: BTreeMap<Vec<FaceId>, Id>,
    datums: BTreeMap<String, Id>,
    systems: BTreeMap<Vec<String>, Id>,
    micrometre: Option<Id>,
}

fn key(feature: &[FaceId]) -> Vec<FaceId> {
    let mut key = feature.to_vec();
    key.sort();
    key
}

impl Emitter<'_> {
    fn root(&mut self, name: &'static str, fields: Vec<Field>) -> Id {
        let id = self.graph.add(name, fields);
        self.roots.push(id);
        id
    }

    /// Ties a shape aspect to the written faces of a feature.
    fn usage(&mut self, aspect: Id, feature: &[FaceId], item: &str) -> Result<(), Refusal> {
        let mut written = Vec::new();
        for face in key(feature) {
            let ids = self.faces.get(&face).ok_or_else(|| {
                Refusal::new(
                    "step.pmi-face-not-written",
                    format!("{item}: a face it is attached to is not in the file"),
                )
            })?;
            written.extend(ids.iter().copied());
        }
        let common = vec![
            text("''"),
            text("''"),
            aspect.into(),
            header(REPRESENTATION),
        ];
        if let [face] = written.as_slice() {
            let mut fields = common;
            fields.push((*face).into());
            self.root("GEOMETRIC_ITEM_SPECIFIC_USAGE", fields);
        } else {
            let mut fields = common;
            fields.push(Field::Typed("SET_REPRESENTATION_ITEM", vec![refs(written)]));
            self.root("ITEM_IDENTIFIED_REPRESENTATION_USAGE", fields);
        }
        Ok(())
    }

    /// The one shape aspect of a feature.
    fn aspect(&mut self, feature: &[FaceId], item: &str) -> Result<Id, Refusal> {
        let key = key(feature);
        if let Some(id) = self.aspects.get(&key) {
            return Ok(*id);
        }
        let id = self.graph.add(
            "SHAPE_ASPECT",
            vec![
                text("''"),
                text("''"),
                header(PRODUCT_DEFINITION_SHAPE),
                boolean(true),
            ],
        );
        self.usage(id, feature, item)?;
        self.aspects.insert(key, id);
        Ok(id)
    }

    fn datum(&mut self, letter: &str, feature: &[FaceId], item: &str) -> Result<(), Refusal> {
        let aspect = self.graph.add(
            "DATUM_FEATURE",
            vec![
                text("''"),
                text("''"),
                header(PRODUCT_DEFINITION_SHAPE),
                boolean(true),
            ],
        );
        self.usage(aspect, feature, item)?;
        self.aspects.insert(key(feature), aspect);
        let datum = self.graph.add(
            "DATUM",
            vec![
                text("''"),
                text("''"),
                header(PRODUCT_DEFINITION_SHAPE),
                boolean(false),
                quoted(letter),
            ],
        );
        self.root(
            "SHAPE_ASPECT_RELATIONSHIP",
            vec![text("''"), text("''"), aspect.into(), datum.into()],
        );
        self.datums.insert(letter.to_owned(), datum);
        self.root("ID_ATTRIBUTE", vec![quoted(item), aspect.into()]);
        Ok(())
    }

    /// The aspect a size or tolerance applies to: the feature's own, or the
    /// group of a pattern's features.
    fn target(&mut self, features: &[Feature<FaceId>], item: &str) -> Result<Id, Refusal> {
        if let [feature] = features {
            return self.aspect(feature, item);
        }
        let mut members: Vec<Vec<FaceId>> = features.iter().map(|feature| key(feature)).collect();
        members.sort();
        if let Some(id) = self.groups.get(&members) {
            return Ok(*id);
        }
        let group = self.graph.add(
            "COMPOSITE_GROUP_SHAPE_ASPECT",
            vec![
                text("''"),
                text("'pattern of features'"),
                header(PRODUCT_DEFINITION_SHAPE),
                boolean(true),
            ],
        );
        for feature in features {
            let member = self.aspect(feature, item)?;
            self.root(
                "SHAPE_ASPECT_RELATIONSHIP",
                vec![text("''"), text("''"), group.into(), member.into()],
            );
        }
        self.groups.insert(members, group);
        Ok(group)
    }

    /// The aspect a distance ends on.
    fn anchor(&mut self, anchor: &Anchor<FaceId>, item: &str) -> Result<Id, Refusal> {
        if anchor.at != At::Middle {
            return self.aspect(&anchor.feature, item);
        }
        let key = key(&anchor.feature);
        if let Some(id) = self.centres.get(&key) {
            return Ok(*id);
        }
        let centre = self.graph.add(
            "CENTRE_OF_SYMMETRY",
            vec![
                text("''"),
                text("''"),
                header(PRODUCT_DEFINITION_SHAPE),
                boolean(false),
            ],
        );
        for face in &key {
            let side = self.aspect(&[*face], item)?;
            self.root(
                "SHAPE_ASPECT_DERIVING_RELATIONSHIP",
                vec![text("''"), text("''"), centre.into(), side.into()],
            );
        }
        self.centres.insert(key, centre);
        Ok(centre)
    }

    fn round(&self, feature: &[FaceId]) -> bool {
        feature.first().is_some_and(|face| {
            matches!(
                self.model.geom().surface(self.model.face(*face).surface),
                Surface::Cylinder(_)
            )
        })
    }

    /// The name AP242 gives a distance by where it ends on round features.
    fn location_name(&self, from: &Anchor<FaceId>, to: &Anchor<FaceId>) -> String {
        let side = |anchor: &Anchor<FaceId>| {
            self.round(&anchor.feature).then_some(match anchor.at {
                At::Feature | At::Middle => "centre",
                At::Inner => "inner",
                At::Outer => "outer",
            })
        };
        match (side(from), side(to)) {
            (Some("centre"), Some("centre")) | (None, None) => "linear distance".to_owned(),
            (Some(a), Some(b)) => format!("linear distance {a} {b}"),
            (Some("centre"), None) | (None, Some("centre")) => "linear distance".to_owned(),
            (Some(a), None) | (None, Some(a)) => format!("linear distance {a}"),
        }
    }

    fn length(&mut self, value: &Rational, item: &str) -> Result<Id, Refusal> {
        Ok(self.graph.add(
            "LENGTH_MEASURE_WITH_UNIT",
            vec![
                text(format!("LENGTH_MEASURE({})", real(value, item)?)),
                header(MILLIMETRE),
            ],
        ))
    }

    /// The value, limits and name of a dimension.
    fn dimension(
        &mut self,
        dimension: Id,
        nominal: &str,
        limits: Option<&Limits>,
        item: &str,
    ) -> Result<(), Refusal> {
        let value = self.graph.add_complex(vec![
            ("LENGTH_MEASURE_WITH_UNIT", vec![]),
            ("MEASURE_REPRESENTATION_ITEM", vec![]),
            (
                "MEASURE_WITH_UNIT",
                vec![
                    text(format!("POSITIVE_LENGTH_MEASURE({nominal})")),
                    header(MILLIMETRE),
                ],
            ),
            ("REPRESENTATION_ITEM", vec![text("'nominal value'")]),
        ]);
        let items = vec![value];
        let representation = self.graph.add(
            "SHAPE_DIMENSION_REPRESENTATION",
            vec![text("''"), refs(items), header(CONTEXT)],
        );
        self.root(
            "DIMENSIONAL_CHARACTERISTIC_REPRESENTATION",
            vec![dimension.into(), representation.into()],
        );
        match limits {
            None => {}
            Some(Limits::Deviations { upper, lower }) => {
                let lower = self.length(lower, item)?;
                let upper = self.length(upper, item)?;
                let range = self
                    .graph
                    .add("TOLERANCE_VALUE", vec![lower.into(), upper.into()]);
                self.root("PLUS_MINUS_TOLERANCE", vec![range.into(), dimension.into()]);
            }
            Some(Limits::Fit { deviation, grade }) => {
                let zone = if deviation.bytes().all(|b| b.is_ascii_uppercase()) {
                    "'hole'"
                } else {
                    "'shaft'"
                };
                let range = self.graph.add(
                    "LIMITS_AND_FITS",
                    vec![quoted(deviation), text(zone), quoted(grade), text("''")],
                );
                self.root("PLUS_MINUS_TOLERANCE", vec![range.into(), dimension.into()]);
            }
        }
        self.root("ID_ATTRIBUTE", vec![quoted(item), dimension.into()]);
        Ok(())
    }

    fn system(&mut self, letters: &[String], item: &str) -> Result<Id, Refusal> {
        if let Some(id) = self.systems.get(letters) {
            return Ok(*id);
        }
        let mut compartments = Vec::with_capacity(letters.len());
        for letter in letters {
            let datum = *self.datums.get(letter).ok_or_else(|| {
                Refusal::new(
                    "step.pmi-invalid",
                    format!("{item}: datum {letter} is not defined"),
                )
            })?;
            compartments.push(self.graph.add(
                "DATUM_REFERENCE_COMPARTMENT",
                vec![
                    text("''"),
                    text("$"),
                    header(PRODUCT_DEFINITION_SHAPE),
                    boolean(false),
                    datum.into(),
                    text("$"),
                ],
            ));
        }
        let id = self.graph.add(
            "DATUM_SYSTEM",
            vec![
                quoted(&letters.join("|")),
                text("''"),
                header(PRODUCT_DEFINITION_SHAPE),
                boolean(false),
                refs(compartments),
            ],
        );
        self.systems.insert(letters.to_vec(), id);
        Ok(id)
    }

    fn tolerance(
        &mut self,
        tolerance: &GeometricTolerance<FaceId>,
        item: &str,
    ) -> Result<(), Refusal> {
        let target = self.target(&tolerance.features, item)?;
        let magnitude = self.length(&tolerance.value, item)?;
        let base = vec![quoted(item), text("''"), magnitude.into(), target.into()];
        let entity = tolerance.characteristic.entity();
        // An entity with one leaf type is a plain instance: every form
        // tolerance, and a characteristic that always has datums, whose
        // entity is a subtype of the one that holds them. Several leaf
        // types make a complex instance.
        let datums_inherited = tolerance.characteristic.datum_use() == super::DatumUse::Required;
        let id = if tolerance.datums.is_empty() && tolerance.material.is_none() {
            self.root(entity, base)
        } else if datums_inherited && tolerance.material.is_none() {
            let system = self.system(&tolerance.datums, item)?;
            let mut fields = base;
            fields.push(refs([system]));
            self.root(entity, fields)
        } else {
            let mut parts = vec![("GEOMETRIC_TOLERANCE", base), (entity, vec![])];
            if !tolerance.datums.is_empty() {
                let system = self.system(&tolerance.datums, item)?;
                parts.push((
                    "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE",
                    vec![refs([system])],
                ));
            }
            if let Some(material) = tolerance.material {
                let modifier = match material {
                    MaterialCondition::Maximum => ".MAXIMUM_MATERIAL_REQUIREMENT.",
                    MaterialCondition::Least => ".LEAST_MATERIAL_REQUIREMENT.",
                };
                parts.push((
                    "GEOMETRIC_TOLERANCE_WITH_MODIFIERS",
                    vec![Field::List(vec![text(modifier)])],
                ));
            }
            parts.sort_by_key(|part| part.0);
            let id = self.graph.add_complex(parts);
            self.roots.push(id);
            id
        };
        if tolerance.diameter_zone {
            let form = self.graph.add(
                "TOLERANCE_ZONE_FORM",
                vec![text("'cylindrical or circular'")],
            );
            self.root(
                "TOLERANCE_ZONE",
                vec![
                    text("''"),
                    text("''"),
                    header(PRODUCT_DEFINITION_SHAPE),
                    boolean(false),
                    refs([id]),
                    form.into(),
                ],
            );
        }
        Ok(())
    }

    fn texture(&mut self, texture: &SurfaceTexture<FaceId>, item: &str) -> Result<(), Refusal> {
        let aspect = self.aspect(&texture.feature, item)?;
        let micrometre = match self.micrometre {
            Some(unit) => unit,
            None => {
                let unit = self.graph.add_complex(vec![
                    ("LENGTH_UNIT", vec![]),
                    ("NAMED_UNIT", vec![text("*")]),
                    ("SI_UNIT", vec![text(".MICRO."), text(".METRE.")]),
                ]);
                self.micrometre = Some(unit);
                unit
            }
        };
        let property = self.graph.add(
            "PROPERTY_DEFINITION",
            vec![
                text("'surface_condition'"),
                quoted(texture.process.phrase()),
                aspect.into(),
            ],
        );
        let measure = self.graph.add(
            "MEASURE_REPRESENTATION_ITEM",
            vec![
                quoted(&texture.parameter),
                text(format!("LENGTH_MEASURE({})", real(&texture.value, item)?)),
                micrometre.into(),
            ],
        );
        let method = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![text("'measuring method'"), text("'ISO 4287'")],
        );
        let representation = self.graph.add(
            "SURFACE_TEXTURE_REPRESENTATION",
            vec![text("''"), refs([measure, method]), header(CONTEXT)],
        );
        self.root(
            "PROPERTY_DEFINITION_REPRESENTATION",
            vec![property.into(), representation.into()],
        );
        let general = self.graph.add(
            "GENERAL_PROPERTY",
            vec![quoted(item), text("'surface_condition'"), text("''")],
        );
        self.root(
            "GENERAL_PROPERTY_ASSOCIATION",
            vec![text("''"), text("''"), general.into(), property.into()],
        );
        self.root("ID_ATTRIBUTE", vec![quoted(item), property.into()]);
        Ok(())
    }

    /// A note on the whole part, as a user defined attribute.
    fn note(&mut self, name: &str, value: &str) {
        let property = self.graph.add(
            "PROPERTY_DEFINITION",
            vec![
                quoted(name),
                text("'user defined attribute'"),
                header(PRODUCT_DEFINITION),
            ],
        );
        let item = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![quoted(name), quoted(value)],
        );
        let representation = self.graph.add(
            "REPRESENTATION",
            vec![text("''"), refs([item]), header(CONTEXT)],
        );
        self.root(
            "PROPERTY_DEFINITION_REPRESENTATION",
            vec![property.into(), representation.into()],
        );
        let general = self.graph.add(
            "GENERAL_PROPERTY",
            vec![quoted(name), quoted(name), text("'user defined attribute'")],
        );
        self.root(
            "GENERAL_PROPERTY_ASSOCIATION",
            vec![text("''"), text("''"), general.into(), property.into()],
        );
    }

    /// States that the tolerances read by ISO 1101.
    fn standard(&mut self) {
        let application = self.graph.add(
            "APPLICATION_CONTEXT",
            vec![text(
                "'geometrical dimensioning and tolerancing representation'",
            )],
        );
        let context = self.graph.add(
            "PRODUCT_DEFINITION_CONTEXT",
            vec![text("''"), application.into(), text("'design'")],
        );
        let role = self.graph.add(
            "PRODUCT_DEFINITION_CONTEXT_ROLE",
            vec![text("'additional context'"), text("''")],
        );
        self.root(
            "PRODUCT_DEFINITION_CONTEXT_ASSOCIATION",
            vec![header(PRODUCT_DEFINITION), context.into(), role.into()],
        );
        let kind = self.graph.add(
            "DOCUMENT_TYPE",
            vec![text("'configuration controlled document version'")],
        );
        let document = self.graph.add(
            "DOCUMENT",
            vec![
                text("'ISO 1101:2017'"),
                text("'Geometrical tolerancing'"),
                text("''"),
                kind.into(),
            ],
        );
        self.root(
            "APPLIED_DOCUMENT_REFERENCE",
            vec![document.into(), text("''"), refs([context])],
        );
    }
}

/// The names part notes are written under, in order.
pub(crate) const GENERAL_TOLERANCE: &str = "general tolerance";
pub(crate) const MATERIAL: &str = "material";
pub(crate) const NOTE: &str = "note";

/// Adds the entities of `pmi` to `graph` and returns the ones nothing else
/// refers to, for the writer to number after the solids. `faces` maps each
/// model face to its written `ADVANCED_FACE` entities, and `nominals` holds
/// the measured value of each dimension as STEP text, by item.
pub(crate) fn emit(
    graph: &mut Graph,
    model: &Model,
    faces: &BTreeMap<FaceId, Vec<Id>>,
    pmi: &Pmi<FaceId>,
    nominals: &[Option<String>],
) -> Result<Vec<Id>, Refusal> {
    pmi.validate()?;
    let mut emitter = Emitter {
        graph,
        model,
        faces,
        roots: Vec::new(),
        aspects: BTreeMap::new(),
        groups: BTreeMap::new(),
        centres: BTreeMap::new(),
        datums: BTreeMap::new(),
        systems: BTreeMap::new(),
        micrometre: None,
    };
    if !pmi.items.is_empty() {
        emitter.standard();
    }
    // A datum's feature is a DATUM_FEATURE whatever else is attached to it,
    // so datums come first.
    for item in &pmi.items {
        if let ItemKind::Datum(datum) = &item.kind {
            emitter.datum(&datum.letter, &datum.feature, &item.id)?;
        }
    }
    for (item, nominal) in pmi.items.iter().zip(nominals) {
        let id = item.id.as_str();
        let measured = || {
            nominal.as_deref().ok_or_else(|| {
                Refusal::new(
                    "step.pmi-not-measurable",
                    format!("{id}: the dimension has no measured value"),
                )
            })
        };
        match &item.kind {
            ItemKind::Datum(_) => {}
            ItemKind::Size(size) => {
                let target = emitter.target(&size.features, id)?;
                let name = match size.kind {
                    SizeKind::Diameter => "'diameter'",
                    SizeKind::Radius => "'radius'",
                };
                let dimension = emitter
                    .graph
                    .add("DIMENSIONAL_SIZE", vec![target.into(), text(name)]);
                emitter.dimension(dimension, measured()?, size.limits.as_ref(), id)?;
            }
            ItemKind::Location(location) => {
                let from = emitter.anchor(&location.from, id)?;
                let to = emitter.anchor(&location.to, id)?;
                let name = emitter.location_name(&location.from, &location.to);
                let dimension = emitter.graph.add(
                    "DIMENSIONAL_LOCATION",
                    vec![quoted(&name), text("''"), from.into(), to.into()],
                );
                emitter.dimension(dimension, measured()?, location.limits.as_ref(), id)?;
            }
            ItemKind::Tolerance(tolerance) => emitter.tolerance(tolerance, id)?,
            ItemKind::Texture(texture) => emitter.texture(texture, id)?,
        }
    }
    if let Some(value) = &pmi.part.general_tolerance {
        emitter.note(GENERAL_TOLERANCE, value);
    }
    if let Some(value) = &pmi.part.material {
        emitter.note(MATERIAL, value);
    }
    for (index, value) in pmi.part.notes.iter().enumerate() {
        emitter.note(&format!("{NOTE} {}", index + 1), value);
    }
    Ok(emitter.roots)
}
