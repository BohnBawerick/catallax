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

use std::collections::BTreeMap;

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Loop, Model};

use super::callout::{self, presented};
use super::{
    Anchor, At, Feature, GeometricTolerance, ItemKind, Limits, MaterialCondition, Pmi, SizeKind,
    SurfaceTexture,
};
use crate::document::PartReferences;
use crate::graph::{Field, Graph, Id, boolean, refs, string, text};
use crate::number::decimal;

/// The attribute of a hole feature that says how deep it goes.
pub const HOLE_DEPTH: &str = "hole depth";
/// The [`HOLE_DEPTH`] of a hole that goes through the part.
pub const THROUGH: &str = "through";
/// The name of the user defined attribute that carries a thread's designation.
pub const THREAD: &str = "thread";

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
    plane: Id,
    style: Id,
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
    systems: BTreeMap<Vec<String>, Id>,
    micrometre: Option<Id>,
    presentation: Option<Presentation>,
    annotations: usize,
    dec: &'a Decider,
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
            self.references.product.representation.into(),
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
                self.references.product.shape.into(),
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
                    self.references.length_unit.into(),
                ],
            ),
            ("REPRESENTATION_ITEM", vec![text("'nominal value'")]),
        ]);
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
                    self.references.product.shape.into(),
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
                self.references.product.shape.into(),
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
                    self.references.product.shape.into(),
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
        self.callout(texture, item)
    }

    /// The draughting model, made with the first callout: one annotation
    /// plane facing the front, whose x axis reads along +x and whose y axis
    /// is +z (recommended practices, section 9.1), and one black continuous
    /// curve style.
    fn presentation(&mut self) -> Presentation {
        if let Some(presentation) = self.presentation {
            return presentation;
        }
        let origin = self
            .graph
            .add("CARTESIAN_POINT", vec![text("''"), text("(0.,0.,0.)")]);
        let towards = self
            .graph
            .add("DIRECTION", vec![text("''"), text("(0.,-1.,0.)")]);
        let along = self
            .graph
            .add("DIRECTION", vec![text("''"), text("(1.,0.,0.)")]);
        let frame = self.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), origin.into(), towards.into(), along.into()],
        );
        let surface = self.graph.add("PLANE", vec![text("''"), frame.into()]);
        let unstyled = self.graph.add(
            "PRESENTATION_STYLE_ASSIGNMENT",
            vec![Field::List(vec![Field::Typed(
                "NULL_STYLE",
                vec![text(".NULL.")],
            )])],
        );
        // The callouts join the plane's elements as they are drawn.
        let plane = self.graph.add(
            "ANNOTATION_PLANE",
            vec![
                text("''"),
                refs([unstyled]),
                surface.into(),
                Field::List(Vec::new()),
            ],
        );
        let model = self.graph.add(
            "DRAUGHTING_MODEL",
            vec![text("''"), refs([plane]), self.references.context.into()],
        );
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
        let style = self
            .graph
            .add("PRESENTATION_STYLE_ASSIGNMENT", vec![refs([curve])]);
        let presentation = Presentation {
            model,
            plane,
            style,
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

    /// The graphic callout of a surface texture (recommended practices,
    /// sections 8.1 and 9.3.1), linked to an aspect of its own on the
    /// texture's faces. On the feature's shared aspect OpenCASCADE merges the
    /// callout into a dimension there; on its own it reads a presentation on
    /// the faces.
    fn callout(&mut self, texture: &SurfaceTexture<FaceId>, item: &str) -> Result<(), Refusal> {
        let anchor = self.anchor_point(&texture.feature)?;
        let value = plain(&texture.value, item)?;
        let drawing = callout::draw(texture, &value, anchor.as_ref());
        let presentation = self.presentation();
        let mut curves = Vec::with_capacity(drawing.polylines.len());
        for line in &drawing.polylines {
            let points: Vec<Id> = line.iter().map(|at| self.point("", at)).collect();
            curves.push(self.graph.add("POLYLINE", vec![text("''"), refs(points)]));
        }
        let set = self.graph.add(
            "GEOMETRIC_CURVE_SET",
            vec![text("'surface roughness'"), refs(curves)],
        );
        let occurrence = self.graph.add(
            "ANNOTATION_CURVE_OCCURRENCE",
            vec![quoted(item), refs([presentation.style]), set.into()],
        );
        let callout = self
            .graph
            .add("DRAUGHTING_CALLOUT", vec![quoted(item), refs([occurrence])]);
        self.graph.extend_list(presentation.plane, 3, [callout]);
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
        self.root(
            "DRAUGHTING_MODEL_ITEM_ASSOCIATION",
            vec![
                quoted(item),
                text("''"),
                aspect.into(),
                presentation.model.into(),
                callout.into(),
            ],
        );
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
        let (length, centre) = drawing.length_and_centre();
        let string = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![
                text("'equivalent unicode string'"),
                quoted(&callout::equivalent_string(texture.process, &drawing.text)),
            ],
        );
        let length = self.graph.add(
            "MEASURE_REPRESENTATION_ITEM",
            vec![
                text("'polyline curve length'"),
                Field::Typed(
                    "POSITIVE_LENGTH_MEASURE",
                    vec![text(decimal(
                        &presented(&length)
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
            ("number of semantic pmi elements", pmi.items.len()),
            (
                "number of dimensional locations",
                of(|kind| matches!(kind, ItemKind::Location(_))),
            ),
            (
                "number of dimensional sizes",
                of(|kind| matches!(kind, ItemKind::Size(_))),
            ),
            (
                "number of geometric tolerances",
                of(|kind| matches!(kind, ItemKind::Tolerance(_))),
            ),
            ("number of composite tolerances", 0),
            (
                "number of datums",
                of(|kind| matches!(kind, ItemKind::Datum(_))),
            ),
            ("number of datum targets", 0),
            ("number of annotations", self.annotations),
        ];
        let items = counts
            .into_iter()
            .map(|(name, count)| {
                // The practices write an integer item's value as a real.
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
        let property = self.graph.add(
            "PROPERTY_DEFINITION",
            vec![quoted(name), text("'user defined attribute'"), of],
        );
        let item = self.graph.add(
            "DESCRIPTIVE_REPRESENTATION_ITEM",
            vec![quoted(name), quoted(value)],
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
/// model face to its written `ADVANCED_FACE` entities. The caller validates
/// the part and measures its dimension nominals before calling this. `dec`
/// places a callout's leader on a face with no rational vertex.
pub(crate) fn emit(
    graph: &mut Graph,
    model: &Model,
    faces: &BTreeMap<FaceId, Vec<Id>>,
    references: PartReferences,
    pmi: &Pmi<FaceId>,
    nominals: &[Option<String>],
    dec: &Decider,
) -> Result<Vec<Id>, Refusal> {
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
        presentation: None,
        annotations: 0,
        dec,
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
                emitter.dimension(dimension, measured()?, location.limits.as_ref(), id)?;
                let string = distance_string(measured()?, location.limits.as_ref(), id)?;
                emitter.equivalent(dimension, &string);
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
    if !pmi.items.is_empty() {
        emitter.counts(pmi);
    }
    Ok(emitter.roots)
}
