//! AP242 entities for the manufacturing information of one part.
//!
//! Each feature becomes one `SHAPE_ASPECT` tied to its written faces, shared
//! by every item on that feature. A datum's feature is a `DATUM_FEATURE`. A
//! pattern is a `COMPOSITE_GROUP_SHAPE_ASPECT` of its features. Sizes and
//! distances are `DIMENSIONAL_SIZE` and `DIMENSIONAL_LOCATION` with a
//! nominal value and, where given, limits. Geometric tolerances refer to one
//! shared `DATUM_SYSTEM` per datum sequence. A surface texture is a
//! `SURFACE_TEXTURE_REPRESENTATION` on the feature, and is also drawn as a
//! graphic callout (`super::callout`) on an aspect of its own, the one form of
//! it OpenCASCADE reads. Part notes are user defined attributes of the product
//! definition.
//!
//! The file carries the PMI validation properties of the CAx-IF recommended
//! practices (version 4.1, section 10) that AP242 edition 1 allows: the counts
//! of the part's items and annotations, the equivalent Unicode string of each
//! datum and distance, and the string, length and centre of each callout. A
//! size or a geometric tolerance takes a property only from edition 2 on.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Loop, Model};

use super::callout::{self, presented};
use super::{
    Anchor, At, Camera, DatumRef, DatumTarget, Feature, GeometricTolerance, HoleHead, HoleSize,
    ItemKind, Limits, MaterialCondition, Pmi, PositionControl, PositionFrame, PositionSegment,
    SavedView, Size, SizeKind, SurfaceTexture, TargetArea, TargetKind, check, measure,
};
use crate::document::PartReferences;
use crate::graph::{Field, Graph, Id, boolean, refs, string, text};
use crate::number::{self, decimal};

/// The attribute of a hole feature that says how deep it goes.
pub const HOLE_DEPTH: &str = "hole depth";
/// The [`HOLE_DEPTH`] of a hole that goes through the part.
pub const THROUGH: &str = "through";
/// The name of the user defined attribute that carries a thread's designation.
pub const THREAD: &str = "thread";
/// AP242 edition 1 fallback for a typed Catallax hole: the pilot's aspect
/// carries its diameter as semantic PMI. The head relationship and measured
/// properties retain hole intent for readers that know this convention.
pub const HOLE_PILOT: &str = "catallax hole pilot";
/// The relationship from a pilot aspect to a cylindrical head.
pub const COUNTERBORE: &str = "counterbore";
/// The relationship from a pilot aspect to a conical head.
pub const COUNTERSINK: &str = "countersink";
/// The measured property for the head's diameter.
pub const HEAD_DIAMETER: &str = "head diameter";
/// The measured property for axial counterbore depth.
pub const HEAD_DEPTH: &str = "head depth";
/// The measured property for the full countersink angle.
pub const INCLUDED_ANGLE: &str = "included angle";
/// The measured property for usable axial thread depth.
pub const THREAD_DEPTH: &str = "thread depth";
/// The dimensional size name for an opposed width.
pub const WIDTH: &str = "width";
/// The dimensional size name for the Catallax chamfer convention.
pub const CHAMFER: &str = "chamfer";
/// The dimensional size name for the Catallax bounded edge break convention.
pub const EDGE_BREAK: &str = "edge break";
/// The relationship to the first planar side of a planar size.
pub const FIRST_SIDE: &str = "first side";
/// The relationship to the second planar side of a planar size.
pub const SECOND_SIDE: &str = "second side";
/// The chamfer's measured angle to the reference side.
pub const BEVEL_ANGLE: &str = "bevel angle";
/// The lower setback bound of an edge break.
pub const MIN_SETBACK: &str = "minimum setback";
/// The upper setback bound of an edge break.
pub const MAX_SETBACK: &str = "maximum setback";
/// The name of the relationship from an upper composite segment to its lower
/// one (recommended practices, section 6.9.9).
pub const COMPOSITE: &str = "composite";

/// The supported ASME standard for composite position frames.
pub const COMPOSITE_STANDARD: &str = "ASME Y14.5-2018";
/// The description of the relationship from a datum target to its datum.
pub const DATUM_TARGET: &str = "datum target";
/// The name of a datum target's placement (section 6.6.1, table 10).
pub const ORIENTATION: &str = "orientation";
/// The length of a line or rectangular datum target.
pub const TARGET_LENGTH: &str = "target length";
/// The width of a rectangular datum target.
pub const TARGET_WIDTH: &str = "target width";
/// The diameter of a circular datum target.
pub const TARGET_DIAMETER: &str = "target diameter";
/// The name of the property that carries editable text (recommended
/// practices, section 7.4.2).
pub const SEMANTIC_TEXT: &str = "semantic text";
/// The variable name of a knurl designation held as semantic text.
pub const KNURL: &str = "knurl designation";

/// The AP242 modifier of a material condition or boundary.
fn modifier(material: MaterialCondition) -> &'static str {
    match material {
        MaterialCondition::Maximum => ".MAXIMUM_MATERIAL_REQUIREMENT.",
        MaterialCondition::Least => ".LEAST_MATERIAL_REQUIREMENT.",
    }
}

/// The letter a drawing circles for a material condition or boundary.
fn letter(material: MaterialCondition) -> &'static str {
    match material {
        MaterialCondition::Maximum => "M",
        MaterialCondition::Least => "L",
    }
}

/// A nonzero exact vector square to `normal`, for the reference direction of
/// a point or circle target, which has no direction of its own.
fn perpendicular(normal: &[Rational; 3]) -> [Rational; 3] {
    let [x, y, _] = normal;
    if x.is_zero() && y.is_zero() {
        [Rational::one(), Rational::zero(), Rational::zero()]
    } else {
        [&Rational::zero() - y, x.clone(), Rational::zero()]
    }
}

enum Nominal<'a> {
    Length(&'a str),
    Angle(&'a str),
}

fn quoted(value: &str) -> Field {
    text(string(value))
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

/// A deviation as a drawing writes it: `+0.1`, `0` or `-0.1`.
fn signed(value: &Rational, item: &str) -> Result<String, Refusal> {
    let text = plain(value, item)?;
    Ok(if *value > Rational::zero() {
        format!("+{text}")
    } else {
        text
    })
}

/// The equivalent Unicode string of a distance (CAx-IF PMI Unicode string
/// practice, version 1.0, sections 4.2, 5.3 and 6.1), such as `DIM\w15±0.1`,
/// `DIM\w[25]` for a basic distance or `DIM\w6\u+0.1\n0\u` for stacked limits.
fn distance_string(nominal: &str, limits: Option<&Limits>, item: &str) -> Result<String, Refusal> {
    let value = nominal.strip_suffix('.').unwrap_or(nominal);
    Ok(match limits {
        None => format!("DIM\\w{value}"),
        Some(Limits::Basic) => format!("DIM\\w[{value}]"),
        Some(Limits::Fit { deviation, grade }) => format!("DIM\\w{value} {deviation}{grade}"),
        Some(Limits::Deviations { upper, lower }) => {
            if *upper == &Rational::zero() - lower {
                format!("DIM\\w{value}±{}", plain(upper, item)?)
            } else {
                format!(
                    "DIM\\w{value}\\u{}\\n{}\\u",
                    signed(upper, item)?,
                    signed(lower, item)?
                )
            }
        }
    })
}

/// Whether a refusal met while finding a point for a leader stops the export:
/// the decider's deadline or budget, as opposed to a curve or surface that
/// has no point to give there.
fn stops(refusal: &Refusal) -> bool {
    refusal.code.starts_with("numbers.")
}

/// The global draughting model of the part's callouts, its annotation plane
/// and the curve style they share.
#[derive(Clone, Copy)]
struct Presentation {
    model: Id,
    fill: Id,
    style: Id,
}

struct CalloutLinks {
    semantic: Id,
    geometry: Option<Id>,
    semantic_on_occurrence: bool,
}

struct Emitter<'a> {
    graph: &'a mut Graph,
    model: &'a Model,
    faces: &'a BTreeMap<FaceId, Vec<Id>>,
    references: PartReferences,
    roots: Vec<Id>,
    aspects: BTreeMap<Vec<FaceId>, Id>,
    groups: BTreeMap<Vec<Vec<FaceId>>, Id>,
    centres: BTreeMap<Vec<FaceId>, Id>,
    datums: BTreeMap<String, Id>,
    systems: BTreeMap<Vec<(String, u8)>, Id>,
    micrometre: Option<Id>,
    degree: Option<Id>,
    presentation: Option<Presentation>,
    front_plane: Option<Id>,
    annotations: usize,
    /// The semantic entity each item is written as, by item name, for the
    /// association of its callout.
    semantic: BTreeMap<String, Id>,
    reversed: BTreeSet<Id>,
    dec: &'a mut Decider,
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
        self.written_usage(aspect, &written);
        Ok(())
    }

    fn written_usage(&mut self, aspect: Id, written: &[Id]) {
        let common = vec![
            text("''"),
            text("''"),
            aspect.into(),
            self.references.product.representation.into(),
        ];
        if let [face] = written {
            let mut fields = common;
            fields.push((*face).into());
            self.root("GEOMETRIC_ITEM_SPECIFIC_USAGE", fields);
        } else {
            let mut fields = common;
            fields.push(Field::Typed(
                "SET_REPRESENTATION_ITEM",
                vec![refs(written.iter().copied())],
            ));
            self.root("ITEM_IDENTIFIED_REPRESENTATION_USAGE", fields);
        }
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
                self.references.product.shape.into(),
                boolean(true),
            ],
        );
        self.usage(id, feature, item)?;
        self.aspects.insert(key, id);
        Ok(id)
    }

    fn hole(&mut self, size: &Size<FaceId>, hole: &HoleSize, item: &str) -> Result<Id, Refusal> {
        // A private aspect keeps per-item intent off any datum or independent
        // size dimension that happens to name the same face.
        let pilot = self.graph.add(
            "SHAPE_ASPECT",
            vec![
                quoted(HOLE_PILOT),
                text("''"),
                self.references.product.shape.into(),
                boolean(true),
            ],
        );
        self.usage(pilot, &size.features[0], item)?;
        if let Some(head) = &hole.head {
            let (kind, diameter) = match head {
                HoleHead::Counterbore { diameter, .. } => (COUNTERBORE, diameter),
                HoleHead::Countersink { diameter, .. } => (COUNTERSINK, diameter),
            };
            let head_aspect = self.aspect(&size.features[1], item)?;
            self.root(
                "SHAPE_ASPECT_RELATIONSHIP",
                vec![quoted(kind), text("''"), pilot.into(), head_aspect.into()],
            );
            self.attribute_length(HEAD_DIAMETER, diameter, pilot, item)?;
            match head {
                HoleHead::Counterbore { depth, .. } => {
                    self.attribute_length(HEAD_DEPTH, depth, pilot, item)?;
                }
                HoleHead::Countersink { included_angle, .. } => {
                    let value = real(included_angle, item)?;
                    let measure = self.angle_measure(INCLUDED_ANGLE, &value);
                    self.attribute_item(INCLUDED_ANGLE, measure, pilot.into());
                }
            }
        }
        if let Some(depth) = &hole.drill_depth {
            self.attribute_length(HOLE_DEPTH, depth, pilot, item)?;
        }
        if let Some(depth) = &hole.thread_depth {
            self.attribute_length(THREAD_DEPTH, depth, pilot, item)?;
        }
        Ok(pilot)
    }

    fn attribute_length(
        &mut self,
        name: &str,
        value: &Rational,
        aspect: Id,
        item: &str,
    ) -> Result<(), Refusal> {
        let measure = self.named_length(name, value, item)?;
        self.attribute_item(name, measure, aspect.into());
        Ok(())
    }

    fn angle_measure(&mut self, name: &str, value: &str) -> Id {
        let unit = self.degree();
        self.graph.add_complex(vec![
            ("MEASURE_REPRESENTATION_ITEM", vec![]),
            (
                "MEASURE_WITH_UNIT",
                vec![text(format!("PLANE_ANGLE_MEASURE({value})")), unit.into()],
            ),
            ("PLANE_ANGLE_MEASURE_WITH_UNIT", vec![]),
            ("REPRESENTATION_ITEM", vec![quoted(name)]),
        ])
    }

    fn planar_size(
        &mut self,
        name: &str,
        bevel: Option<&Feature<FaceId>>,
        sides: &[Feature<FaceId>; 2],
        item: &str,
    ) -> Result<(Id, Id), Refusal> {
        // The side roles and properties belong to this dimension, not to
        // another control (or datum) that names one of the same faces.
        let aspect = self.graph.add(
            if bevel.is_some() {
                "SHAPE_ASPECT"
            } else {
                "COMPOSITE_SHAPE_ASPECT"
            },
            vec![
                text("''"),
                text("''"),
                self.references.product.shape.into(),
                boolean(false),
            ],
        );
        if let Some(bevel) = bevel {
            self.usage(aspect, bevel, item)?;
        }
        for (role, side) in [FIRST_SIDE, SECOND_SIDE].into_iter().zip(sides) {
            let side = self.aspect(side, item)?;
            self.root(
                "SHAPE_ASPECT_RELATIONSHIP",
                vec![quoted(role), text("''"), aspect.into(), side.into()],
            );
        }
        let dimension = self
            .graph
            .add("DIMENSIONAL_SIZE", vec![aspect.into(), quoted(name)]);
        Ok((aspect, dimension))
    }

    fn degree(&mut self) -> Id {
        *self.degree.get_or_insert_with(|| {
            let radians = self.graph.add_complex(vec![
                ("NAMED_UNIT", vec![text("*")]),
                ("PLANE_ANGLE_UNIT", vec![]),
                ("SI_UNIT", vec![text("$"), text(".RADIAN.")]),
            ]);
            let conversion = self.graph.add(
                "PLANE_ANGLE_MEASURE_WITH_UNIT",
                vec![
                    text("PLANE_ANGLE_MEASURE(0.017453292519943295)"),
                    radians.into(),
                ],
            );
            let dimension = self.graph.add("DIMENSIONAL_EXPONENTS", vec![text("0."); 7]);
            self.graph.add_complex(vec![
                (
                    "CONVERSION_BASED_UNIT",
                    vec![text("'DEGREE'"), conversion.into()],
                ),
                ("NAMED_UNIT", vec![dimension.into()]),
                ("PLANE_ANGLE_UNIT", vec![]),
            ])
        })
    }

    fn datum(&mut self, letter: &str, feature: &[FaceId], item: &str) -> Result<(), Refusal> {
        let aspect = self.graph.add(
            "DATUM_FEATURE",
            vec![
                text("''"),
                text("''"),
                self.references.product.shape.into(),
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
                self.references.product.shape.into(),
                boolean(false),
                quoted(letter),
            ],
        );
        self.root(
            "SHAPE_ASPECT_RELATIONSHIP",
            vec![text("''"), text("''"), aspect.into(), datum.into()],
        );
        self.datums.insert(letter.to_owned(), datum);
        self.semantic.insert(item.to_owned(), aspect);
        self.root("ID_ATTRIBUTE", vec![quoted(item), aspect.into()]);
        self.equivalent(datum, &format!("DTM\\w{letter}"));
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
                self.references.product.shape.into(),
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
                self.references.product.shape.into(),
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
                self.references.length_unit.into(),
            ],
        ))
    }

    /// The value, limits and name of a dimension.
    fn dimension(
        &mut self,
        dimension: Id,
        nominal: Nominal<'_>,
        limits: Option<&Limits>,
        item: &str,
    ) -> Result<(), Refusal> {
        let value = match nominal {
            Nominal::Angle(value) => self.angle_measure("nominal value", value),
            Nominal::Length(value) => self.graph.add_complex(vec![
                ("LENGTH_MEASURE_WITH_UNIT", vec![]),
                ("MEASURE_REPRESENTATION_ITEM", vec![]),
                (
                    "MEASURE_WITH_UNIT",
                    vec![
                        text(format!("POSITIVE_LENGTH_MEASURE({value})")),
                        self.references.length_unit.into(),
                    ],
                ),
                ("REPRESENTATION_ITEM", vec![text("'nominal value'")]),
            ]),
        };
        let mut items = vec![value];
        if limits == Some(&Limits::Basic) {
            items.push(self.graph.add(
                "DESCRIPTIVE_REPRESENTATION_ITEM",
                vec![text("'dimensional note'"), text("'theoretical'")],
            ));
        }
        let representation = self.graph.add(
            "SHAPE_DIMENSION_REPRESENTATION",
            vec![text("''"), refs(items), self.references.context.into()],
        );
        self.root(
            "DIMENSIONAL_CHARACTERISTIC_REPRESENTATION",
            vec![dimension.into(), representation.into()],
        );
        match limits {
            None | Some(Limits::Basic) => {}
            Some(Limits::Deviations { upper, lower }) => {
                let lower = self.length(lower, item)?;
                let upper = self.length(upper, item)?;
                let range = self
                    .graph
                    .add("TOLERANCE_VALUE", vec![lower.into(), upper.into()]);
                self.root("PLUS_MINUS_TOLERANCE", vec![range.into(), dimension.into()]);
            }
            Some(Limits::Fit { deviation, grade }) => {
                let zone = if super::is_bore_fit(deviation) {
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

    /// The one datum system of a sequence of references. A reference with a
    /// material boundary is a compartment of its own, so `A|B` and `A|B(M)`
    /// are different systems (recommended practices, section 6.9.7).
    fn system(&mut self, datums: &[DatumRef], item: &str) -> Result<Id, Refusal> {
        let key: Vec<(String, u8)> = datums
            .iter()
            .map(|datum| {
                let boundary = match datum.material {
                    None => 0,
                    Some(MaterialCondition::Maximum) => 1,
                    Some(MaterialCondition::Least) => 2,
                };
                (datum.letter.clone(), boundary)
            })
            .collect();
        if let Some(id) = self.systems.get(&key) {
            return Ok(*id);
        }
        let mut compartments = Vec::with_capacity(datums.len());
        for reference in datums {
            let letter = &reference.letter;
            let datum = *self.datums.get(letter).ok_or_else(|| {
                Refusal::new(
                    "step.pmi-invalid",
                    format!("{item}: datum {letter} is not defined"),
                )
            })?;
            let modifiers = match reference.material {
                None => text("$"),
                Some(material) => Field::List(vec![Field::Typed(
                    "SIMPLE_DATUM_REFERENCE_MODIFIER",
                    vec![text(modifier(material))],
                )]),
            };
            compartments.push(self.graph.add(
                "DATUM_REFERENCE_COMPARTMENT",
                vec![
                    text("''"),
                    text("$"),
                    self.references.product.shape.into(),
                    boolean(false),
                    datum.into(),
                    modifiers,
                ],
            ));
        }
        let name: Vec<String> = datums
            .iter()
            .map(|datum| match datum.material {
                None => datum.letter.clone(),
                Some(material) => format!("{}({})", datum.letter, letter(material)),
            })
            .collect();
        let id = self.graph.add(
            "DATUM_SYSTEM",
            vec![
                quoted(&name.join("|")),
                text("''"),
                self.references.product.shape.into(),
                boolean(false),
                refs(compartments),
            ],
        );
        self.systems.insert(key, id);
        Ok(id)
    }

    fn tolerance(
        &mut self,
        tolerance: &GeometricTolerance<FaceId>,
        item: &str,
    ) -> Result<(), Refusal> {
        let target = self.target(&tolerance.features, item)?;
        let datums: Vec<DatumRef> = tolerance
            .datums
            .iter()
            .map(|letter| DatumRef::new(letter))
            .collect();
        let (id, _) = self.geometric(
            tolerance.characteristic,
            &tolerance.value,
            tolerance.diameter_zone,
            tolerance.material,
            &datums,
            target,
            item,
        )?;
        self.semantic.insert(item.to_owned(), id);
        Ok(())
    }

    /// One geometric tolerance entity on `target`, and its cylindrical zone
    /// if it has one.
    #[allow(clippy::too_many_arguments)]
    fn geometric(
        &mut self,
        characteristic: super::Characteristic,
        value: &Rational,
        diameter_zone: bool,
        material: Option<MaterialCondition>,
        datums: &[DatumRef],
        target: Id,
        item: &str,
    ) -> Result<(Id, Option<Id>), Refusal> {
        let magnitude = self.length(value, item)?;
        let base = vec![quoted(item), text("''"), magnitude.into(), target.into()];
        let entity = characteristic.entity();
        // An entity with one leaf type is a plain instance: every form
        // tolerance, and a characteristic that always has datums, whose
        // entity is a subtype of the one that holds them. Several leaf
        // types make a complex instance.
        let datums_inherited = characteristic.datum_use() == super::DatumUse::Required;
        let id = if datums.is_empty() && material.is_none() {
            self.root(entity, base)
        } else if datums_inherited && material.is_none() {
            let system = self.system(datums, item)?;
            let mut fields = base;
            fields.push(refs([system]));
            self.root(entity, fields)
        } else {
            let mut parts = vec![("GEOMETRIC_TOLERANCE", base), (entity, vec![])];
            if !datums.is_empty() {
                let system = self.system(datums, item)?;
                parts.push((
                    "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE",
                    vec![refs([system])],
                ));
            }
            if let Some(material) = material {
                parts.push((
                    "GEOMETRIC_TOLERANCE_WITH_MODIFIERS",
                    vec![Field::List(vec![text(modifier(material))])],
                ));
            }
            parts.sort_by_key(|part| part.0);
            let id = self.graph.add_complex(parts);
            self.roots.push(id);
            id
        };
        let zone = if diameter_zone {
            let form = self.graph.add(
                "TOLERANCE_ZONE_FORM",
                vec![text("'cylindrical or circular'")],
            );
            Some(self.root(
                "TOLERANCE_ZONE",
                vec![
                    text("''"),
                    text("''"),
                    self.references.product.shape.into(),
                    boolean(false),
                    refs([id]),
                    form.into(),
                ],
            ))
        } else {
            None
        };
        Ok((id, zone))
    }

    fn segment(
        &mut self,
        segment: &PositionSegment,
        target: Id,
        item: &str,
    ) -> Result<(Id, Option<Id>), Refusal> {
        self.geometric(
            super::Characteristic::Position,
            &segment.value_mm,
            segment.diameter_zone,
            segment.material,
            &segment.datums,
            target,
            item,
        )
    }

    /// A position control (recommended practices, sections 6.9.2.2 and
    /// 6.9.9). Both segments of a composite are position tolerances on one
    /// aspect, related upper to lower by a `composite` relationship. A
    /// projected zone starts on an aspect of the entry faces.
    fn position(
        &mut self,
        control: &PositionControl<FaceId>,
        datum_features: &BTreeMap<&str, &Vec<FaceId>>,
        item: &str,
    ) -> Result<(), Refusal> {
        let mut axis = |feature: &[FaceId]| {
            check::written_axis(
                self.graph,
                feature
                    .iter()
                    .flat_map(|face| self.faces[face].iter().copied()),
                item,
                self.dec,
            )
        };
        let axes = control
            .features
            .iter()
            .map(|feature| axis(feature))
            .collect::<Result<Vec<_>, _>>()?;
        let segments: &[_] = match &control.frame {
            PositionFrame::Single { segment, .. } => &[segment],
            PositionFrame::Composite { upper, lower } => &[upper, lower],
        };
        for feature in segments
            .iter()
            .flat_map(|segment| &segment.datums)
            .filter(|datum| datum.material.is_some())
            .filter_map(|datum| datum_features.get(datum.letter.as_str()))
        {
            check::written_datum_boundary(
                self.graph,
                feature
                    .iter()
                    .flat_map(|face| self.faces[face].iter().copied())
                    .collect(),
                &self.reversed,
                item,
                self.dec,
            )?;
        }
        let target = self.target(&control.features, item)?;
        match &control.frame {
            PositionFrame::Single {
                segment,
                projected_zone,
            } => {
                let (_, zone) = self.segment(segment, target, item)?;
                if let Some(projection) = projected_zone {
                    check::check_written_projection(
                        self.graph,
                        self.faces,
                        &self.reversed,
                        control.features.iter().zip(&axes),
                        projection,
                        item,
                        self.dec,
                    )?;
                    let zone = zone.ok_or_else(|| {
                        Refusal::new(
                            "step.pmi-invalid",
                            format!("{item}: a projected zone needs a cylindrical zone"),
                        )
                    })?;
                    let mut entries = projection.entry_faces.clone();
                    entries.sort();
                    entries.dedup();
                    let end = self.aspect(&entries, item)?;
                    let height = self.length(&projection.height_mm, item)?;
                    self.root(
                        "PROJECTED_ZONE_DEFINITION",
                        vec![
                            zone.into(),
                            Field::List(Vec::new()),
                            end.into(),
                            height.into(),
                        ],
                    );
                }
            }
            PositionFrame::Composite { upper, lower } => {
                let (upper, _) = self.segment(upper, target, item)?;
                let (lower, _) = self.segment(lower, target, item)?;
                self.root(
                    "GEOMETRIC_TOLERANCE_RELATIONSHIP",
                    vec![quoted(COMPOSITE), text("''"), upper.into(), lower.into()],
                );
            }
        }
        Ok(())
    }

    fn named_length(&mut self, name: &str, value: &Rational, item: &str) -> Result<Id, Refusal> {
        let value = real(value, item)?;
        Ok(self.graph.add_complex(vec![
            ("LENGTH_MEASURE_WITH_UNIT", vec![]),
            ("MEASURE_REPRESENTATION_ITEM", vec![]),
            (
                "MEASURE_WITH_UNIT",
                vec![
                    text(format!("LENGTH_MEASURE({value})")),
                    self.references.length_unit.into(),
                ],
            ),
            ("REPRESENTATION_ITEM", vec![quoted(name)]),
        ]))
    }

    fn exact_point(&mut self, at: &[Rational; 3], item: &str) -> Result<Id, Refusal> {
        let xyz = at
            .iter()
            .map(|value| real(value, item).map(text))
            .collect::<Result<_, _>>()?;
        Ok(self
            .graph
            .add("CARTESIAN_POINT", vec![text("''"), Field::List(xyz)]))
    }

    fn direction(&mut self, at: &[Rational; 3], item: &str) -> Result<Id, Refusal> {
        let xyz = super::decimal_direction(at)
            .iter()
            .map(|value| real(value, item).map(text))
            .collect::<Result<_, _>>()?;
        Ok(self
            .graph
            .add("DIRECTION", vec![text("''"), Field::List(xyz)]))
    }

    fn face_text(
        &mut self,
        variable: &str,
        value: &str,
        feature: &[FaceId],
        item: &str,
    ) -> Result<(), Refusal> {
        let aspect = self.graph.add(
            "SHAPE_ASPECT",
            vec![
                quoted(item),
                text("''"),
                self.references.product.shape.into(),
                boolean(true),
            ],
        );
        self.usage(aspect, feature, item)?;
        let value = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![quoted(variable), quoted(value)],
        );
        let property = self.graph.add(
            "PROPERTY_DEFINITION",
            vec![quoted(SEMANTIC_TEXT), text("''"), aspect.into()],
        );
        let representation = self.graph.add(
            "REPRESENTATION",
            vec![text("''"), refs([value]), self.references.context.into()],
        );
        self.root(
            "PROPERTY_DEFINITION_REPRESENTATION",
            vec![property.into(), representation.into()],
        );
        self.root("ID_ATTRIBUTE", vec![quoted(item), property.into()]);
        Ok(())
    }

    /// A placed datum target (recommended practices, section 6.6): its
    /// `orientation` placement and sizes as parameters, its relationship to
    /// its datum, and its relationship to an aspect of its face (6.6.3).
    fn datum_target(&mut self, target: &DatumTarget<FaceId>, item: &str) -> Result<(), Refusal> {
        let datum = *self.datums.get(&target.datum).ok_or_else(|| {
            Refusal::new(
                "step.pmi-invalid",
                format!("{item}: datum {} is not defined", target.datum),
            )
        })?;
        let not_exact = || {
            Refusal::new(
                "step.pmi-invalid",
                format!("{item}: the target face's normal is not an exact vector"),
            )
        };
        let (_, normal) = check::outward_plane(self.model, target.face).ok_or_else(not_exact)?;
        check::check_written_target(
            self.graph,
            &self.faces[&target.face],
            &self.reversed,
            target,
            &normal,
            item,
            self.dec,
        )?;
        let normal = [&normal.x, &normal.y, &normal.z]
            .map(|value| value.as_rational().ok_or_else(not_exact));
        let [x, y, z] = normal;
        let normal = [x?, y?, z?];
        let reference = match target.kind.direction() {
            Some(direction) => direction.clone(),
            None => perpendicular(&normal),
        };
        let placed = self.graph.add(
            "PLACED_DATUM_TARGET_FEATURE",
            vec![
                quoted(item),
                quoted(target.kind.description()),
                self.references.product.shape.into(),
                boolean(true),
                quoted(&target.number.to_string()),
            ],
        );
        self.root(
            "SHAPE_ASPECT_RELATIONSHIP",
            vec![
                text("''"),
                quoted(DATUM_TARGET),
                placed.into(),
                datum.into(),
            ],
        );
        let face = self.graph.add(
            "SHAPE_ASPECT",
            vec![
                text("''"),
                text("''"),
                self.references.product.shape.into(),
                boolean(true),
            ],
        );
        self.usage(face, &[target.face], item)?;
        self.root(
            "FEATURE_FOR_DATUM_TARGET_RELATIONSHIP",
            vec![text("''"), text("$"), face.into(), placed.into()],
        );
        let origin = self.exact_point(&target.point_mm, item)?;
        let axis = self.direction(&normal, item)?;
        let along = self.direction(&reference, item)?;
        let placement = self.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![
                quoted(ORIENTATION),
                origin.into(),
                axis.into(),
                along.into(),
            ],
        );
        let mut items = vec![placement];
        match &target.kind {
            TargetKind::Point => {}
            TargetKind::Line { length_mm, .. } => {
                items.push(self.named_length(TARGET_LENGTH, length_mm, item)?);
            }
            TargetKind::Area(TargetArea::Circular { diameter_mm }) => {
                items.push(self.named_length(TARGET_DIAMETER, diameter_mm, item)?);
            }
            TargetKind::Area(TargetArea::Rectangular {
                length_mm,
                width_mm,
                ..
            }) => {
                items.push(self.named_length(TARGET_LENGTH, length_mm, item)?);
                items.push(self.named_length(TARGET_WIDTH, width_mm, item)?);
            }
        }
        let property = self.graph.add(
            "PROPERTY_DEFINITION",
            vec![text("''"), text("$"), placed.into()],
        );
        let representation = self.graph.add(
            "SHAPE_REPRESENTATION_WITH_PARAMETERS",
            vec![text("''"), refs(items), self.references.context.into()],
        );
        self.root(
            "SHAPE_DEFINITION_REPRESENTATION",
            vec![property.into(), representation.into()],
        );
        self.root("ID_ATTRIBUTE", vec![quoted(item), placed.into()]);
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
            vec![
                text("''"),
                refs([measure, method]),
                self.references.context.into(),
            ],
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
        self.semantic.insert(item.to_owned(), property);
        Ok(())
    }

    fn curve_style(&mut self) -> Id {
        let font = self.graph.add(
            "DRAUGHTING_PRE_DEFINED_CURVE_FONT",
            vec![text("'continuous'")],
        );
        let colour = self
            .graph
            .add("DRAUGHTING_PRE_DEFINED_COLOUR", vec![text("'black'")]);
        let curve = self.graph.add(
            "CURVE_STYLE",
            vec![
                text("''"),
                font.into(),
                Field::Typed("POSITIVE_LENGTH_MEASURE", vec![text("0.25")]),
                colour.into(),
            ],
        );
        self.graph
            .add("PRESENTATION_STYLE_ASSIGNMENT", vec![refs([curve])])
    }

    fn presentation(&mut self) -> Presentation {
        if let Some(presentation) = self.presentation {
            return presentation;
        }
        let colour = self
            .graph
            .add("DRAUGHTING_PRE_DEFINED_COLOUR", vec![text("'black'")]);
        let solid = self
            .graph
            .add("FILL_AREA_STYLE_COLOUR", vec![text("''"), colour.into()]);
        let fill = self
            .graph
            .add("FILL_AREA_STYLE", vec![text("''"), refs([solid])]);
        let fill = self
            .graph
            .add("PRESENTATION_STYLE_ASSIGNMENT", vec![refs([fill])]);
        let model = self.graph.add(
            "DRAUGHTING_MODEL",
            vec![text("''"), refs([]), self.references.context.into()],
        );
        let presentation = Presentation {
            model,
            fill,
            style: self.curve_style(),
        };
        self.presentation = Some(presentation);
        presentation
    }

    /// A coordinate of the model as the drawing places it.
    fn place(&self, value: &Real) -> Result<Rational, Refusal> {
        Ok(presented(&match value.as_rational() {
            Some(value) => value,
            None => self.dec.enclose_at(value, 64)?.mid(),
        }))
    }

    /// A point on one of the feature's faces for a callout's leader to end
    /// on: a vertex of the faces, else a point of an edge's curve, else of a
    /// face's surface. `None` when none of them can be evaluated.
    fn anchor_point(&self, feature: &[FaceId]) -> Result<Option<callout::Point>, Refusal> {
        let faces = key(feature);
        for face in &faces {
            for boundary in &self.model.face(*face).loops {
                let vertex = match self.model.loop_(*boundary) {
                    Loop::Vertex(vertex) => Some(*vertex),
                    Loop::Coedges(coedges) => coedges.iter().find_map(|coedge| {
                        let edge = self.model.edge(self.model.coedge(*coedge).edge);
                        edge.start.or(edge.end)
                    }),
                };
                if let Some(vertex) = vertex {
                    let at = self.model.geom().position(self.model.vertex(vertex).point);
                    return Ok(Some([
                        self.place(&at.x)?,
                        self.place(&at.y)?,
                        self.place(&at.z)?,
                    ]));
                }
            }
        }
        // Faces bounded only by ring edges, or by none.
        let start = Ball::from_rational(&Rational::zero(), 64);
        let mid = |point: [Ball; 3]| point.map(|ball| presented(&ball.mid()));
        for face in &faces {
            for boundary in &self.model.face(*face).loops {
                if let Loop::Coedges(coedges) = self.model.loop_(*boundary) {
                    for coedge in coedges {
                        let edge = self.model.edge(self.model.coedge(*coedge).edge);
                        let curve = self.model.geom().curve(edge.curve);
                        match curve.evaluate(&start, self.dec) {
                            Ok(Some(point)) => return Ok(Some(mid(point))),
                            Err(refusal) if stops(&refusal) => return Err(refusal),
                            Ok(None) | Err(_) => {}
                        }
                    }
                }
            }
            let surface = self.model.geom().surface(self.model.face(*face).surface);
            match surface.evaluate(&start, &start, self.dec) {
                Ok(point) => return Ok(Some(mid(point))),
                Err(refusal) if stops(&refusal) => return Err(refusal),
                Err(_) => {}
            }
        }
        Ok(None)
    }

    fn point(&mut self, name: &str, at: &callout::Point) -> Id {
        let xyz = at
            .iter()
            .map(|c| {
                text(decimal(
                    &c.to_decimal_exact()
                        .expect("a drawing coordinate is a decimal"),
                ))
            })
            .collect();
        self.graph
            .add("CARTESIAN_POINT", vec![quoted(name), Field::List(xyz)])
    }

    fn texture_callout(
        &mut self,
        texture: &SurfaceTexture<FaceId>,
        item: &str,
        slot: usize,
        view: Option<(&Camera, usize)>,
    ) -> Result<Id, Refusal> {
        let anchor = self.anchor_point(&texture.feature)?;
        let value = plain(&texture.value, item)?;
        let mut drawing = callout::draw(
            texture,
            &value,
            anchor.as_ref().filter(|_| view.is_none()),
            if view.is_some() { 0 } else { slot },
        );
        let equivalent = callout::equivalent_string(texture.process, &drawing.text);
        let plane = self.callout_plane(
            view.map(|(camera, count)| (camera, slot, count)),
            anchor.as_ref(),
            &mut drawing,
            item,
        )?;
        let aspect = self.graph.add(
            "SHAPE_ASPECT",
            vec![
                quoted(item),
                text("''"),
                self.references.product.shape.into(),
                boolean(true),
            ],
        );
        self.usage(aspect, &texture.feature, item)?;
        self.graphic_callout(
            item,
            "surface roughness",
            &drawing,
            &equivalent,
            plane,
            CalloutLinks {
                semantic: self.semantic[item],
                geometry: Some(aspect),
                semantic_on_occurrence: false,
            },
        )
    }

    fn graphic_attachment(&mut self, item: &super::Item<FaceId>) -> Id {
        let written: BTreeSet<_> = item
            .kind
            .faces()
            .into_iter()
            .flat_map(|face| self.faces[face].iter().copied())
            .collect();
        let owner = self.references.product.shape;
        let aspect_fields = || vec![text("''"), text("''"), owner.into(), boolean(true)];
        // Each written face gets its own usage because Edition 1 permits one per aspect.
        let members: Vec<_> = written
            .into_iter()
            .map(|face| {
                let aspect = self.graph.add("SHAPE_ASPECT", aspect_fields());
                self.written_usage(aspect, &[face]);
                aspect
            })
            .collect();
        if let [member] = members.as_slice() {
            return *member;
        }
        let group = self
            .graph
            .add("COMPOSITE_GROUP_SHAPE_ASPECT", aspect_fields());
        for member in members {
            self.root(
                "SHAPE_ASPECT_RELATIONSHIP",
                vec![text("''"), text("''"), group.into(), member.into()],
            );
        }
        group
    }

    fn label_callout(
        &mut self,
        item: &super::Item<FaceId>,
        label: Label<'_>,
        view: (&Camera, usize, usize),
        repeated: bool,
    ) -> Result<Id, Refusal> {
        let id = item.id.as_str();
        if !callout::drawable(&label.line) {
            return Err(Refusal::new(
                "step.pmi-invalid",
                format!("{id}: the callout font cannot draw {}", label.line),
            ));
        }
        let anchor = self.anchor_point(label.feature)?;
        let mut drawing = callout::label(&label.line);
        let plane = self.callout_plane(Some(view), anchor.as_ref(), &mut drawing, id)?;
        let links = CalloutLinks {
            semantic: self.semantic[id],
            geometry: repeated.then(|| self.graphic_attachment(item)),
            semantic_on_occurrence: repeated,
        };
        self.graphic_callout(id, label.kind, &drawing, &label.equivalent, plane, links)
    }

    fn callout_plane(
        &mut self,
        view: Option<(&Camera, usize, usize)>,
        anchor: Option<&callout::Point>,
        drawing: &mut callout::Drawing,
        item: &str,
    ) -> Result<Id, Refusal> {
        let presentation = self.presentation();
        let (origin, normal, across) = if let Some((camera, slot, count)) = view {
            let anchor = anchor.ok_or_else(|| {
                Refusal::new(
                    "step.pmi-invalid",
                    format!("{item}: the callout has no face anchor"),
                )
            })?;
            let [across, normal] = drawing
                .in_view(anchor, camera, slot, count)
                .map_err(|reason| Refusal::new("step.pmi-invalid", format!("{item}: {reason}")))?;
            (
                anchor.clone(),
                normal.map(|value| &Rational::zero() - &value),
                across,
            )
        } else {
            if let Some(plane) = self.front_plane {
                return Ok(plane);
            }
            let [zero, one, minus] = [0, 1, -1].map(Rational::from_i64);
            (
                [zero.clone(), zero.clone(), zero.clone()],
                [zero.clone(), minus, zero.clone()],
                [one, zero.clone(), zero],
            )
        };
        let origin = self.point("", &origin);
        let normal = self.direction(&normal, item)?;
        let across = self.direction(&across, item)?;
        let frame = self.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), origin.into(), normal.into(), across.into()],
        );
        let surface = self.graph.add("PLANE", vec![text("''"), frame.into()]);
        let plane = self.graph.add(
            "ANNOTATION_PLANE",
            vec![
                text("''"),
                refs([presentation.fill]),
                surface.into(),
                refs([]),
            ],
        );
        self.graph.extend_list(presentation.model, 1, [plane]);
        if view.is_none() {
            self.front_plane = Some(plane);
        }
        Ok(plane)
    }

    fn graphic_callout(
        &mut self,
        item: &str,
        set_name: &str,
        drawing: &callout::Drawing,
        equivalent: &str,
        plane: Id,
        links: CalloutLinks,
    ) -> Result<Id, Refusal> {
        let (length, centre) = drawing.length_and_centre();
        let length = presented(&length);
        if length.is_zero() {
            return Err(Refusal::new(
                "step.pmi-invalid",
                format!("{item}: the callout length rounds to zero"),
            ));
        }
        let presentation = self.presentation();
        let mut curves = Vec::with_capacity(drawing.polylines.len());
        for line in &drawing.polylines {
            let points: Vec<Id> = line.iter().map(|at| self.point("", at)).collect();
            curves.push(self.graph.add("POLYLINE", vec![text("''"), refs(points)]));
        }
        let set = self
            .graph
            .add("GEOMETRIC_CURVE_SET", vec![quoted(set_name), refs(curves)]);
        let occurrence = self.graph.add(
            "ANNOTATION_CURVE_OCCURRENCE",
            vec![quoted(item), refs([presentation.style]), set.into()],
        );
        let callout = self
            .graph
            .add("DRAUGHTING_CALLOUT", vec![quoted(item), refs([occurrence])]);
        self.graph.extend_list(plane, 3, [callout]);
        // Open CASCADE uses the first association to find a graphic callout's faces.
        // Its semantic label can store only the first view's presentation.
        let subject = if links.semantic_on_occurrence {
            occurrence
        } else {
            callout
        };
        let geometry = links
            .geometry
            .map(|id| ("PMI presentation to geometry link", id, callout));
        let semantic = (
            "PMI representation to presentation link",
            links.semantic,
            subject,
        );
        for (name, definition, subject) in geometry.into_iter().chain([semantic]) {
            self.root(
                "DRAUGHTING_MODEL_ITEM_ASSOCIATION",
                vec![
                    quoted(name),
                    text("''"),
                    definition.into(),
                    presentation.model.into(),
                    subject.into(),
                ],
            );
        }
        // What the callout shows, how long its lines are and where they
        // stand (section 10.2.3).
        let within = self.graph.add(
            "CHARACTERIZED_ITEM_WITHIN_REPRESENTATION",
            vec![
                text("''"),
                text("''"),
                callout.into(),
                presentation.model.into(),
            ],
        );
        let string = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![text("'equivalent unicode string'"), quoted(equivalent)],
        );
        let length = self.graph.add(
            "MEASURE_REPRESENTATION_ITEM",
            vec![
                text("'polyline curve length'"),
                Field::Typed(
                    "POSITIVE_LENGTH_MEASURE",
                    vec![text(decimal(
                        &length
                            .to_decimal_exact()
                            .expect("a rounded length is a decimal"),
                    ))],
                ),
                self.references.length_unit.into(),
            ],
        );
        let centre = self.point("polyline centre point", &centre.map(|c| presented(&c)));
        self.validation(within.into(), vec![string, length, centre]);
        self.annotations += 1;
        Ok(callout)
    }

    /// A saved view (recommended practices, section 9.4): a draughting model
    /// named for the view holding its callouts, part geometry and camera.
    fn view(&mut self, view: &SavedView, mut items: Vec<Id>) -> Result<(), Refusal> {
        let name = view.name.as_str();
        let camera = &view.camera;
        let across = callout::cross(&camera.up, &camera.direction);
        let eye = self.exact_point(&camera.eye, name)?;
        let axis = self.direction(&camera.direction, name)?;
        let across = self.direction(&across, name)?;
        let frame = self.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), eye.into(), axis.into(), across.into()],
        );
        let projection_point = self
            .graph
            .add("CARTESIAN_POINT", vec![text("''"), text("(0.,0.,0.)")]);
        let half = Rational::from_ratio(-1, 2).expect("a nonzero denominator");
        let corner = self.exact_point(
            &[
                &camera.window[0] * &half,
                &camera.window[1] * &half,
                camera.distance.clone(),
            ],
            name,
        )?;
        let window_frame = self.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), corner.into(), text("$"), text("$")],
        );
        let window = self.graph.add(
            "PLANAR_BOX",
            vec![
                text("''"),
                text(real(&camera.window[0], name)?),
                text(real(&camera.window[1], name)?),
                window_frame.into(),
            ],
        );
        let distance = real(&camera.distance, name)?;
        let volume = self.graph.add(
            "VIEW_VOLUME",
            vec![
                text(if camera.perspective {
                    ".CENTRAL."
                } else {
                    ".PARALLEL."
                }),
                projection_point.into(),
                text(distance.clone()),
                text(distance.clone()),
                boolean(false),
                text(distance),
                boolean(false),
                boolean(false),
                window.into(),
            ],
        );
        let camera = self.graph.add(
            "CAMERA_MODEL_D3",
            vec![quoted(name), frame.into(), volume.into()],
        );
        self.root(
            "ITEM_DEFINED_TRANSFORMATION",
            vec![text("''"), text("$"), camera.into(), frame.into()],
        );
        let map = self.graph.add(
            "REPRESENTATION_MAP",
            vec![
                self.references.product.frame.into(),
                self.references.product.representation.into(),
            ],
        );
        let target_point = self
            .graph
            .add("CARTESIAN_POINT", vec![text("''"), text("(0.,0.,0.)")]);
        let target = self.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), target_point.into(), text("$"), text("$")],
        );
        let mapped = self
            .graph
            .add("MAPPED_ITEM", vec![text("''"), map.into(), target.into()]);
        let colour = self.graph.add("COLOUR", vec![]);
        let curve = self.graph.add(
            "CURVE_STYLE",
            vec![text("''"), text("$"), text("$"), colour.into()],
        );
        let assignment = self
            .graph
            .add("PRESENTATION_STYLE_ASSIGNMENT", vec![refs([curve])]);
        let styled = self.graph.add(
            "STYLED_ITEM",
            vec![text("''"), refs([assignment]), mapped.into()],
        );
        items.extend([mapped, styled, camera]);
        let model = self.graph.add(
            "DRAUGHTING_MODEL",
            vec![quoted(name), refs(items), self.references.context.into()],
        );
        if let Some(presentation) = self.presentation {
            self.root(
                "MECHANICAL_DESIGN_AND_DRAUGHTING_RELATIONSHIP",
                vec![
                    text("''"),
                    text("''"),
                    model.into(),
                    presentation.model.into(),
                ],
            );
        } else {
            self.roots.push(model);
        }
        Ok(())
    }

    /// A `pmi validation property` of `of` that holds `items`, all in one
    /// representation (recommended practices, sections 10.1 and 10.3).
    fn validation(&mut self, of: Field, items: Vec<Id>) {
        let property = self.graph.add(
            "PROPERTY_DEFINITION",
            vec![text("'pmi validation property'"), text("''"), of],
        );
        let representation = self.graph.add(
            "REPRESENTATION",
            vec![text("''"), refs(items), self.references.context.into()],
        );
        self.root(
            "PROPERTY_DEFINITION_REPRESENTATION",
            vec![property.into(), representation.into()],
        );
    }

    /// The equivalent Unicode string of a datum or a distance (section
    /// 10.1.3.3).
    fn equivalent(&mut self, of: Id, string: &str) {
        let item = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![text("'equivalent unicode string'"), quoted(string)],
        );
        self.validation(of.into(), vec![item]);
    }

    /// The counts of the part's semantic items by type and of its
    /// annotations (sections 10.1.1 and 10.2.1). Every type of table 22 is
    /// given, so a zero says the part has none.
    fn counts(&mut self, pmi: &Pmi<FaceId>) {
        let of = |kind: fn(&ItemKind<FaceId>) -> bool| {
            pmi.items.iter().filter(|item| kind(&item.kind)).count()
        };
        let counts = [
            (
                "number of semantic pmi elements",
                of(|kind| !kind.is_text()),
            ),
            (
                "number of dimensional locations",
                of(|kind| matches!(kind, ItemKind::Location(_) | ItemKind::BasicAngle { .. })),
            ),
            (
                "number of dimensional sizes",
                of(|kind| {
                    matches!(
                        kind,
                        ItemKind::Size(_)
                            | ItemKind::Width { .. }
                            | ItemKind::Chamfer { .. }
                            | ItemKind::EdgeBreak { .. }
                    )
                }),
            ),
            (
                "number of geometric tolerances",
                pmi.items
                    .iter()
                    .map(|item| match &item.kind {
                        ItemKind::Tolerance(_) => 1,
                        ItemKind::Position(control) => control.frame.segments().len(),
                        _ => 0,
                    })
                    .sum(),
            ),
            (
                "number of composite tolerances",
                of(|kind| {
                    matches!(
                        kind,
                        ItemKind::Position(PositionControl {
                            frame: PositionFrame::Composite { .. },
                            ..
                        })
                    )
                }),
            ),
            (
                "number of datums",
                of(|kind| matches!(kind, ItemKind::Datum(_))),
            ),
            (
                "number of datum targets",
                of(|kind| matches!(kind, ItemKind::DatumTarget(_))),
            ),
            ("number of annotations", self.annotations),
        ];
        // Text notes are counted apart, only by a part that has any
        // (section 7.4.3); all of them are on geometry.
        let notes = of(ItemKind::is_text);
        let text_counts = [
            ("number of semantic text notes", notes),
            ("number of semantic text notes on part", 0),
            ("number of semantic text notes on geometry", notes),
            ("number of semantic text notes on PMI", 0),
        ];
        let counts = counts
            .into_iter()
            .chain(text_counts.into_iter().filter(|_| notes > 0));
        let items = counts
            .into_iter()
            .map(|(name, count)| {
                // AP242 declares INTEGER_REPRESENTATION_ITEM.the_value as REAL.
                self.graph.add(
                    "INTEGER_REPRESENTATION_ITEM",
                    vec![quoted(name), text(format!("{count}."))],
                )
            })
            .collect();
        self.validation(self.references.product.shape.into(), items);
    }

    /// A note on the whole part, as a user defined attribute.
    fn note(&mut self, name: &str, value: &str) {
        self.attribute(name, value, self.references.product.definition.into());
    }

    /// A user defined attribute of the part or of one of its features.
    fn attribute(&mut self, name: &str, value: &str, of: Field) {
        let item = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![quoted(name), quoted(value)],
        );
        self.attribute_item(name, item, of);
    }

    fn attribute_item(&mut self, name: &str, item: Id, of: Field) {
        let property = self.graph.add(
            "PROPERTY_DEFINITION",
            vec![quoted(name), text("'user defined attribute'"), of],
        );
        let representation = self.graph.add(
            "REPRESENTATION",
            vec![text("''"), refs([item]), self.references.context.into()],
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

    fn standard(&mut self, name: &str) {
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
            vec![
                self.references.product.definition.into(),
                context.into(),
                role.into(),
            ],
        );
        let kind = self.graph.add(
            "DOCUMENT_TYPE",
            vec![text("'configuration controlled document version'")],
        );
        let document = self.graph.add(
            "DOCUMENT",
            vec![
                quoted(name),
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

/// A material condition as a label writes it: `(M)`, `(L)` or nothing.
fn circled(material: Option<MaterialCondition>) -> String {
    material.map_or_else(String::new, |material| format!("({})", letter(material)))
}

fn frame(tolerance: &GeometricTolerance<FaceId>, item: &str) -> Result<(String, String), Refusal> {
    let value = plain(&tolerance.value, item)?;
    let diameter = if tolerance.diameter_zone { "⌀" } else { "" };
    let visible_modifier = circled(tolerance.material);
    let modifier = match tolerance.material {
        Some(MaterialCondition::Maximum) => "Ⓜ",
        Some(MaterialCondition::Least) => "Ⓛ",
        None => "",
    };
    let mut visible = vec![
        "POS".to_owned(),
        format!("{diameter}{value}{visible_modifier}"),
    ];
    visible.extend(tolerance.datums.iter().cloned());
    let mut equivalent = vec![
        "FCF".to_owned(),
        "⌖".to_owned(),
        format!("{diameter}{value}{modifier}"),
    ];
    equivalent.extend(tolerance.datums.iter().cloned());
    Ok((visible.join("|"), equivalent.join(r"\w")))
}

struct Label<'a> {
    kind: &'static str,
    line: String,
    equivalent: String,
    feature: &'a [FaceId],
}

fn label<'a>(item: &'a super::Item<FaceId>, nominal: Option<&str>) -> Result<Label<'a>, Refusal> {
    let id = item.id.as_str();
    let (kind, line, equivalent, feature) = match &item.kind {
        ItemKind::Datum(datum) => (
            "datum",
            datum.letter.clone(),
            format!(r"DTM\w{}", datum.letter),
            datum.feature.as_slice(),
        ),
        ItemKind::Size(Size {
            kind: SizeKind::Diameter,
            limits: Some(Limits::Fit { deviation, grade }),
            thread: None,
            through,
            features,
            ..
        }) if super::is_bore_fit(deviation) => {
            let value = nominal.expect("a diameter has a measured value");
            let value = value.strip_suffix('.').unwrap_or(value);
            let mut main = format!("⌀{value} {deviation}{grade}");
            if *through {
                main.push_str(" THRU");
            }
            let equivalent = format!(r"DIM\w{main}");
            (
                "diameter dimension",
                main,
                equivalent,
                features[0].as_slice(),
            )
        }
        ItemKind::Tolerance(tolerance)
            if tolerance.characteristic == super::Characteristic::Position =>
        {
            let (line, equivalent) = frame(tolerance, id)?;
            (
                "position",
                line,
                equivalent,
                tolerance.features[0].as_slice(),
            )
        }
        _ => {
            return Err(Refusal::new(
                "step.pmi-invalid",
                format!("{id}: the writer draws no callout for this item"),
            ));
        }
    };
    Ok(Label {
        kind,
        line,
        equivalent,
        feature,
    })
}

/// The names part notes are written under, in order.
pub(crate) const GENERAL_TOLERANCE: &str = "general tolerance";
pub(crate) const MATERIAL: &str = "material";
pub(crate) const NOTE: &str = "note";

/// Adds the entities of `pmi` to `graph` and returns the ones nothing else
/// refers to, for the writer to number after the solids. `faces` maps each
/// model face to its written `ADVANCED_FACE` entities. The caller validates
/// part before calling this. `dec` measures dimension nominals and places a
/// callout's leader on a face with no rational vertex.
pub(crate) fn emit(
    graph: &mut Graph,
    part: (&Model, BodyId),
    faces: &BTreeMap<FaceId, Vec<Id>>,
    references: PartReferences,
    pmi: &Pmi<FaceId>,
    digits: u32,
    dec: &mut Decider,
) -> Result<Vec<Id>, Refusal> {
    for item in &pmi.items {
        if item
            .kind
            .faces()
            .into_iter()
            .any(|face| !faces.contains_key(face))
        {
            return Err(Refusal::new(
                "step.pmi-face-not-written",
                format!("{}: a face it is attached to is not in the file", item.id),
            ));
        }
    }
    let (model, body) = part;
    let reversed = check::reversed_faces(graph);
    let mut emitter = Emitter {
        graph,
        model,
        faces,
        references,
        roots: Vec::new(),
        aspects: BTreeMap::new(),
        groups: BTreeMap::new(),
        centres: BTreeMap::new(),
        datums: BTreeMap::new(),
        systems: BTreeMap::new(),
        micrometre: None,
        degree: None,
        presentation: None,
        front_plane: None,
        annotations: 0,
        semantic: BTreeMap::new(),
        reversed,
        dec,
    };
    check::check_on_model(model, body, pmi, emitter.dec)?;
    if !pmi.items.is_empty() || pmi.part.standard.is_some() {
        let composite = pmi.items.iter().any(|item| {
            matches!(
                &item.kind,
                ItemKind::Position(PositionControl {
                    frame: PositionFrame::Composite { .. },
                    ..
                })
            )
        });
        let asme = pmi.part.standard == Some(super::DimensioningStandard::AsmeY14_5_2018);
        if composite && !asme {
            return Err(Refusal::new(
                "step.pmi-invalid",
                format!("composite position requires explicit selection of {COMPOSITE_STANDARD}"),
            ));
        }
        if asme
            && let Some(item) = pmi.items.iter().find(|item| {
                matches!(&item.kind, ItemKind::Tolerance(tolerance)
                    if !tolerance.characteristic.supports_asme_2018())
            })
        {
            return Err(Refusal::new(
                "step.pmi-invalid",
                format!(
                    "{}: this tolerance is not supported under {COMPOSITE_STANDARD}",
                    item.id
                ),
            ));
        }
        emitter.standard(if asme {
            COMPOSITE_STANDARD
        } else {
            "ISO 1101:2017"
        });
    }
    // A datum's feature is a DATUM_FEATURE whatever else is attached to it,
    // so datums come first.
    let mut datum_features = BTreeMap::new();
    for item in &pmi.items {
        if let ItemKind::Datum(datum) = &item.kind {
            emitter.datum(&datum.letter, &datum.feature, &item.id)?;
            datum_features.insert(datum.letter.as_str(), &datum.feature);
        }
    }
    let shown: BTreeSet<_> = pmi.views.iter().flat_map(|view| &view.items).collect();
    let mut texture_slot = 0;
    let mut sizes = BTreeMap::new();
    for item in &pmi.items {
        let id = item.id.as_str();
        let measured = measure::dimension(model, item, emitter.dec)?;
        let nominal = measured
            .as_ref()
            .map(|value| {
                if matches!(&item.kind, ItemKind::BasicAngle { .. }) {
                    real(
                        value
                            .value
                            .rational_ref()
                            .expect("measured angle is decimal"),
                        id,
                    )
                } else {
                    Ok(decimal(
                        &number::real(&value.value, digits, emitter.dec)?.text,
                    ))
                }
            })
            .transpose()?;
        let nominal = || nominal.as_deref().expect("dimension has a measured value");
        match &item.kind {
            ItemKind::Datum(_) => {}
            ItemKind::Size(size) => {
                let (target, name) = match &size.kind {
                    SizeKind::Diameter => (emitter.target(&size.features, id)?, "'diameter'"),
                    SizeKind::Radius => (emitter.target(&size.features, id)?, "'radius'"),
                    SizeKind::Hole(hole) => (emitter.hole(size, hole, id)?, "'diameter'"),
                };
                let dimension = emitter
                    .graph
                    .add("DIMENSIONAL_SIZE", vec![target.into(), text(name)]);
                emitter.semantic.insert(id.to_owned(), dimension);
                sizes.insert(id, nominal().to_owned());
                emitter.dimension(
                    dimension,
                    Nominal::Length(nominal()),
                    size.limits.as_ref(),
                    id,
                )?;
                if size.through {
                    // AP242 edition 1 has no hole entity, so the hole's depth is an
                    // attribute of the feature the diameter is on.
                    emitter.attribute(HOLE_DEPTH, THROUGH, target.into());
                }
                if let Some(thread) = &size.thread {
                    // Nor has it a thread entity: the designation is an attribute too.
                    emitter.attribute(THREAD, thread, target.into());
                }
            }
            ItemKind::Location(location) => {
                let from = emitter.anchor(&location.from, id)?;
                let to = emitter.anchor(&location.to, id)?;
                let name = emitter.location_name(&location.from, &location.to);
                let dimension = emitter.graph.add(
                    "DIMENSIONAL_LOCATION",
                    vec![quoted(&name), text("''"), from.into(), to.into()],
                );
                emitter.dimension(
                    dimension,
                    Nominal::Length(nominal()),
                    location.limits.as_ref(),
                    id,
                )?;
                let string = distance_string(nominal(), location.limits.as_ref(), id)?;
                emitter.equivalent(dimension, &string);
            }
            ItemKind::Width { flats, limits, .. } => {
                let (_, dimension) = emitter.planar_size(WIDTH, None, flats, id)?;
                emitter.dimension(dimension, Nominal::Length(nominal()), limits.as_ref(), id)?;
            }
            ItemKind::Chamfer { bevel, sides, .. } => {
                let (aspect, dimension) = emitter.planar_size(CHAMFER, Some(bevel), sides, id)?;
                emitter.dimension(dimension, Nominal::Length(nominal()), None, id)?;
                let angle = measured
                    .as_ref()
                    .and_then(|m| m.bevel_angle.as_ref())
                    .expect("chamfer has a measured decimal angle");
                let angle = emitter.angle_measure(BEVEL_ANGLE, &real(angle, id)?);
                emitter.attribute_item(BEVEL_ANGLE, angle, aspect.into());
            }
            ItemKind::EdgeBreak {
                bevel,
                sides,
                min_mm,
                max_mm,
                ..
            } => {
                let (aspect, dimension) =
                    emitter.planar_size(EDGE_BREAK, Some(bevel), sides, id)?;
                emitter.dimension(dimension, Nominal::Length(nominal()), None, id)?;
                emitter.attribute_length(MIN_SETBACK, min_mm, aspect, id)?;
                emitter.attribute_length(MAX_SETBACK, max_mm, aspect, id)?;
            }
            ItemKind::BasicAngle { planes: [a, b], .. } => {
                let a = emitter.aspect(a, id)?;
                let b = emitter.aspect(b, id)?;
                let angle = measured
                    .as_ref()
                    .and_then(|m| m.value.rational_ref())
                    .expect("basic angle is measured in decimal degrees");
                let relator = match angle.cmp(&Rational::from_i64(90)) {
                    std::cmp::Ordering::Less => ".SMALL.",
                    std::cmp::Ordering::Equal => ".EQUAL.",
                    std::cmp::Ordering::Greater => ".LARGE.",
                };
                let dimension = emitter.graph.add(
                    "ANGULAR_LOCATION",
                    vec![
                        text("'angular distance'"),
                        text("''"),
                        a.into(),
                        b.into(),
                        text(relator),
                    ],
                );
                emitter.dimension(
                    dimension,
                    Nominal::Angle(nominal()),
                    Some(&Limits::Basic),
                    id,
                )?;
                let value = nominal().strip_suffix('.').unwrap_or_else(nominal);
                emitter.equivalent(dimension, &format!("DIM\\w[{value}°]"));
            }
            ItemKind::Tolerance(tolerance) => emitter.tolerance(tolerance, id)?,
            ItemKind::DatumTarget(target) => emitter.datum_target(target, id)?,
            ItemKind::FaceNote { text, feature } => emitter.face_text("", text, feature, id)?,
            ItemKind::Knurl {
                designation,
                feature,
            } => emitter.face_text(KNURL, designation, feature, id)?,
            ItemKind::Position(control) => emitter.position(control, &datum_features, id)?,
            ItemKind::Texture(texture) => {
                emitter.texture(texture, id)?;
                if !shown.contains(&item.id) {
                    emitter.texture_callout(texture, id, texture_slot, None)?;
                    texture_slot += 1;
                }
            }
        }
    }
    let by_name: BTreeMap<_, _> = pmi.items.iter().map(|item| (&item.id, item)).collect();
    let mut shown_labels = BTreeSet::new();
    for view in &pmi.views {
        let mut callouts = Vec::with_capacity(view.items.len());
        for (slot, id) in view.items.iter().enumerate() {
            let item = by_name[id];
            let callout = if let ItemKind::Texture(texture) = &item.kind {
                emitter.texture_callout(
                    texture,
                    id,
                    slot,
                    Some((&view.camera, view.items.len())),
                )?
            } else {
                let label = label(item, sizes.get(id.as_str()).map(String::as_str))?;
                emitter.label_callout(
                    item,
                    label,
                    (&view.camera, slot, view.items.len()),
                    !shown_labels.insert(id),
                )?
            };
            callouts.push(callout);
        }
        emitter.view(view, callouts)?;
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
    if !pmi.items.is_empty() {
        emitter.counts(pmi);
    }
    Ok(emitter.roots)
}
