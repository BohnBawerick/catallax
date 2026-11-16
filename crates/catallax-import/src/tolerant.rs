//! Tolerance-bearing source topology for the tolerant import fallback.
//!
//! Reading preserves geometry, topology and source identity; it grants no solid
//! validity or exact certification. Healing consumes this graph separately.
//! Like OCCT BRep_TVertex/TEdge::UpdateTolerance, growth is monotone. Unlike an
//! unconstrained setter, each growth here is capped, metered and recorded.

use std::collections::BTreeMap;

use catallax_numbers::{
    Decider, Rational, Refusal, Sign, operation_budget_active, with_operation_budget,
};
use serde::{Deserialize, Serialize};

pub use crate::numerical::ImportError as Error;

/// Measured, transactional repair of source topology.
pub mod healing;
/// Per-body admission with retained surface representations on solid refusal.
pub mod opened;
use crate::{Curve, Face, FaceLoop, ImportRefusal, OrientedEdge, StepFile, Surface, Vertex};

/// Initial numerical uncertainty and the maximum healing tolerance, in mm.
/// These are topology controls, not a source-to-result deviation certificate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ToleranceOptions {
    /// Starting tolerance on every read vertex and edge.
    pub initial_mm: Rational,
    /// Healing cannot grow an entity beyond this cap.
    pub maximum_mm: Rational,
}

impl Default for ToleranceOptions {
    fn default() -> Self {
        Self {
            initial_mm: Rational::parse("0.0000001").expect("constant"),
            maximum_mm: Rational::parse("0.01").expect("constant"),
        }
    }
}

/// A retained source vertex with its own current tolerance.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantVertex {
    /// Original source entity and position, in millimetres.
    pub source: Vertex,
    tolerance_mm: Rational,
}
impl TolerantVertex {
    /// Current local topology tolerance, not a movement or global bound.
    pub fn tolerance_mm(&self) -> &Rational {
        &self.tolerance_mm
    }
}

/// One shared source edge, independent of its oriented uses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantEdge {
    /// Original EDGE_CURVE name.
    pub name: String,
    /// Unmodified carrier.
    pub curve: Curve,
    /// Direction against the carrier's natural direction.
    pub curve_same_sense: bool,
    /// Start and end indices in this body's shared vertex table.
    pub vertices: [usize; 2],
    tolerance_mm: Rational,
}
impl TolerantEdge {
    /// Current curve/representation tolerance, in millimetres.
    pub fn tolerance_mm(&self) -> &Rational {
        &self.tolerance_mm
    }
}

/// One source use of a shared edge.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantUse {
    /// Index in the body's edge table.
    pub edge: usize,
    /// Traversal direction against EDGE_CURVE.
    pub forward: bool,
}

/// Original face-bound declaration, without ordering or closure assumptions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantLoop {
    /// Original outer-bound declaration.
    pub outer: bool,
    /// Original bound orientation.
    pub orientation: bool,
    /// Uses in source order, including repeated seam uses.
    pub edges: Vec<TolerantUse>,
    /// A VERTEX_LOOP, when present.
    pub vertex: Option<usize>,
}

/// A face with unchanged carrier and source bounds.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantFace {
    /// Source face identity.
    pub name: String,
    /// Original supported surface data.
    pub surface: Surface,
    /// Original normal sense.
    pub same_sense: bool,
    /// Source loops, without an assertion of finite embedding.
    pub loops: Vec<TolerantLoop>,
}

/// One source shell. Declared closure does not establish material.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantShell {
    /// Source shell key, or the containing solid's key for a solid shell.
    pub name: String,
    /// Source closure declaration.
    pub declared_closed: bool,
    /// Unmodified faces in source order.
    pub faces: Vec<TolerantFace>,
}

/// The kind of an entity whose local tolerance grew.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "name", rename_all = "snake_case")]
pub enum ToleranceEntity {
    /// A source VERTEX_POINT.
    Vertex(String),
    /// A source EDGE_CURVE.
    Edge(String),
}

/// One actual tolerance increase. Values are rational millimetre strings.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ToleranceChange {
    /// Source identity, scoped to its body and immutable source hash.
    pub entity: ToleranceEntity,
    /// Value immediately before this change.
    pub before_mm: String,
    /// Value immediately after this change.
    pub after_mm: String,
    /// Named measurement or repair that required the increase.
    pub cause: String,
    /// Edge whose increase required its endpoint's tolerance to grow.
    pub propagated_from_edge: Option<String>,
}

/// An imported source body, before validity checks and healing.
/// This type cannot be extracted as an exact model or counted as a valid solid.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantBody {
    name: String,
    declared_solid: bool,
    shells: Vec<TolerantShell>,
    vertices: Vec<TolerantVertex>,
    edges: Vec<TolerantEdge>,
    vertex_indices: BTreeMap<String, usize>,
    edge_indices: BTreeMap<String, usize>,
    options: ToleranceOptions,
    maximum_attained_mm: Rational,
    changes: Vec<ToleranceChange>,
}

impl TolerantBody {
    /// Source body name.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Whether the source declares a solid; this is not a validity check.
    pub fn declared_solid(&self) -> bool {
        self.declared_solid
    }
    /// Source shells, in source order.
    pub fn shells(&self) -> &[TolerantShell] {
        &self.shells
    }
    /// Shared vertices, in first-use source order.
    pub fn vertices(&self) -> &[TolerantVertex] {
        &self.vertices
    }
    /// Shared edges, in first-use source order.
    pub fn edges(&self) -> &[TolerantEdge] {
        &self.edges
    }
    /// Controls fixed before reading or healing this body.
    pub fn options(&self) -> &ToleranceOptions {
        &self.options
    }
    /// Largest current entity tolerance; not a geometric displacement bound.
    pub fn maximum_tolerance_mm(&self) -> &Rational {
        &self.maximum_attained_mm
    }
    /// Ordered changes, including propagated endpoint increases.
    pub fn tolerance_changes(&self) -> &[ToleranceChange] {
        &self.changes
    }
    /// Whether tolerance growth has made this body explicitly repaired.
    /// Absence of repairs does not imply validity or exact certification.
    pub fn has_tolerance_repairs(&self) -> bool {
        !self.changes.is_empty()
    }

    /// Grow a source vertex's tolerance. Failure leaves the body unchanged.
    pub fn grow_vertex_tolerance(
        &mut self,
        name: &str,
        required_mm: Rational,
        cause: &str,
        dec: &mut Decider,
    ) -> Result<bool, Error> {
        self.validate_growth(name, &required_mm, cause, dec)?;
        let index = *self
            .vertex_indices
            .get(name)
            .ok_or_else(|| invalid("unknown vertex"))?;
        if required_mm <= self.vertices[index].tolerance_mm {
            return Ok(false);
        }
        self.change_vertex(index, required_mm, cause, None);
        Ok(true)
    }

    /// Grow an edge and, when necessary, both endpoint tolerances atomically.
    /// A vertex's tolerance remains at least that of every incident edge.
    pub fn grow_edge_tolerance(
        &mut self,
        name: &str,
        required_mm: Rational,
        cause: &str,
        dec: &mut Decider,
    ) -> Result<bool, Error> {
        self.validate_growth(name, &required_mm, cause, dec)?;
        let index = *self
            .edge_indices
            .get(name)
            .ok_or_else(|| invalid("unknown edge"))?;
        if required_mm <= self.edges[index].tolerance_mm {
            return Ok(false);
        }
        // Prepay all possible changes before the first write. No budget or cap
        // failure can leave an edge changed without its endpoint propagation.
        dec.charge_work(3)?;
        let edge = &mut self.edges[index];
        self.changes.push(ToleranceChange {
            entity: ToleranceEntity::Edge(name.into()),
            before_mm: edge.tolerance_mm.to_string(),
            after_mm: required_mm.to_string(),
            cause: cause.into(),
            propagated_from_edge: None,
        });
        edge.tolerance_mm = required_mm.clone();
        let endpoints = edge.vertices;
        for vertex in endpoints {
            if self.vertices[vertex].tolerance_mm < required_mm {
                self.change_vertex(vertex, required_mm.clone(), cause, Some(name));
            }
        }
        self.maximum_attained_mm = self.maximum_attained_mm.clone().max(required_mm);
        Ok(true)
    }

    fn validate_growth(
        &self,
        name: &str,
        required: &Rational,
        cause: &str,
        dec: &mut Decider,
    ) -> Result<(), Error> {
        dec.charge_work(1 + name.len() as u64 + cause.len() as u64)?;
        if name.is_empty() || cause.trim().is_empty() || required.sign() == Sign::Negative {
            return Err(invalid(
                "growth needs an entity, named cause and nonnegative tolerance",
            ));
        }
        if required > &self.options.maximum_mm {
            return Err(Refusal::new(
                "import.tolerance-cap",
                format!(
                    "{name}: {required} mm exceeds {} mm",
                    self.options.maximum_mm
                ),
            )
            .into());
        }
        Ok(())
    }

    fn change_vertex(&mut self, index: usize, value: Rational, cause: &str, edge: Option<&str>) {
        let vertex = &mut self.vertices[index];
        self.changes.push(ToleranceChange {
            entity: ToleranceEntity::Vertex(vertex.source.name.clone()),
            before_mm: vertex.tolerance_mm.to_string(),
            after_mm: value.to_string(),
            cause: cause.into(),
            propagated_from_edge: edge.map(str::to_owned),
        });
        vertex.tolerance_mm = value.clone();
        self.maximum_attained_mm = self.maximum_attained_mm.clone().max(value);
    }

    fn vertex(&mut self, source: &Vertex, dec: &mut Decider) -> Result<usize, Error> {
        dec.charge_work(1 + source.name.len() as u64)?;
        if let Some(&index) = self.vertex_indices.get(&source.name) {
            if self.vertices[index].source != *source {
                return Err(invalid("one vertex name describes different geometry"));
            }
            return Ok(index);
        }
        let index = self.vertices.len();
        self.vertex_indices.insert(source.name.clone(), index);
        self.vertices.push(TolerantVertex {
            source: source.clone(),
            tolerance_mm: self.options.initial_mm.clone(),
        });
        self.maximum_attained_mm = self
            .options
            .initial_mm
            .clone()
            .max(self.maximum_attained_mm.clone());
        Ok(index)
    }

    fn edge(&mut self, source: &OrientedEdge, dec: &mut Decider) -> Result<usize, Error> {
        dec.charge_work(1 + source.edge.len() as u64)?;
        if let Curve::BSpline(spline) = &source.curve {
            dec.charge_work(
                (spline.geometry.knots().len() + 4 * spline.geometry.controls().len()) as u64,
            )?;
        }
        let vertices = [
            self.vertex(&source.start, dec)?,
            self.vertex(&source.end, dec)?,
        ];
        if let Some(&index) = self.edge_indices.get(&source.edge) {
            let edge = &self.edges[index];
            if edge.vertices != vertices
                || edge.curve != source.curve
                || edge.curve_same_sense != source.curve_same_sense
            {
                return Err(invalid(
                    "one edge name describes different geometry or endpoints",
                ));
            }
            return Ok(index);
        }
        let index = self.edges.len();
        self.edge_indices.insert(source.edge.clone(), index);
        self.edges.push(TolerantEdge {
            name: source.edge.clone(),
            curve: source.curve.clone(),
            curve_same_sense: source.curve_same_sense,
            vertices,
            tolerance_mm: self.options.initial_mm.clone(),
        });
        Ok(index)
    }

    fn boundary(&mut self, source: &FaceLoop, dec: &mut Decider) -> Result<TolerantLoop, Error> {
        let mut edges = Vec::new();
        for edge in &source.edges {
            edges.push(TolerantUse {
                edge: self.edge(edge, dec)?,
                forward: edge.forward,
            });
        }
        Ok(TolerantLoop {
            outer: source.outer,
            orientation: source.orientation,
            edges,
            vertex: source
                .vertex
                .as_ref()
                .map(|v| self.vertex(v, dec))
                .transpose()?,
        })
    }

    fn shell(
        &mut self,
        name: &str,
        closed: bool,
        faces: &[Face],
        dec: &mut Decider,
    ) -> Result<(), Error> {
        let mut result = Vec::new();
        for face in faces {
            dec.charge_work(1)?;
            if let Surface::BSpline(spline) = &face.surface {
                dec.charge_work(
                    (spline
                        .geometry
                        .knots()
                        .iter()
                        .map(|k| k.len())
                        .sum::<usize>()
                        + 4 * spline.geometry.controls().len()) as u64,
                )?;
            }
            let mut loops = Vec::new();
            for boundary in &face.loops {
                loops.push(self.boundary(boundary, dec)?);
            }
            result.push(TolerantFace {
                name: face.name.clone(),
                surface: face.surface.clone(),
                same_sense: face.same_sense,
                loops,
            });
        }
        self.shells.push(TolerantShell {
            name: name.into(),
            declared_closed: closed,
            faces: result,
        });
        Ok(())
    }

    fn new(name: &str, declared_solid: bool, options: &ToleranceOptions) -> Self {
        Self {
            name: name.into(),
            declared_solid,
            options: options.clone(),
            shells: Vec::new(),
            vertices: Vec::new(),
            edges: Vec::new(),
            vertex_indices: BTreeMap::new(),
            edge_indices: BTreeMap::new(),
            maximum_attained_mm: Rational::zero(),
            changes: Vec::new(),
        }
    }
}

/// Recognized geometry with tolerance-bearing topology. It is not yet admitted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TolerantFile {
    source_sha256: String,
    mm_per_unit: Rational,
    options: ToleranceOptions,
    bodies: Vec<TolerantBody>,
}
impl TolerantFile {
    /// Hash of the immutable source bytes; entity numbers are scoped to it.
    pub fn source_sha256(&self) -> &str {
        &self.source_sha256
    }
    /// Source solids followed by source surface models, as in StepFile.
    pub fn bodies(&self) -> &[TolerantBody] {
        &self.bodies
    }
    fn one_body(&self, index: usize, dec: &mut Decider) -> Result<Self, Error> {
        let body = &self.bodies[index];
        charge_text(&self.source_sha256, dec)?;
        healing::charge_value(&self.mm_per_unit, dec)?;
        healing::charge_body(body, dec)?;
        Ok(Self {
            source_sha256: self.source_sha256.clone(),
            mm_per_unit: self.mm_per_unit.clone(),
            options: body.options.clone(),
            bodies: vec![body.clone()],
        })
    }
    /// Reconstruct the retained source definition for geometry consumers.
    /// Local tolerance changes do not move or replace any source geometry.
    pub fn definition(&self, dec: &mut Decider) -> Result<StepFile, Error> {
        let mut file = StepFile {
            mm_per_unit: self.mm_per_unit.clone(),
            solids: Vec::new(),
            surface_models: Vec::new(),
        };
        for body in &self.bodies {
            let mut shells = Vec::new();
            for shell in &body.shells {
                let mut faces = Vec::new();
                for face in &shell.faces {
                    dec.charge_work(1)?;
                    let mut loops = Vec::new();
                    for boundary in &face.loops {
                        let mut edges = Vec::new();
                        for coedge in &boundary.edges {
                            dec.charge_work(1)?;
                            let edge = &body.edges[coedge.edge];
                            edges.push(OrientedEdge {
                                edge: edge.name.clone(),
                                curve: edge.curve.clone(),
                                curve_same_sense: edge.curve_same_sense,
                                forward: coedge.forward,
                                start: body.vertices[edge.vertices[0]].source.clone(),
                                end: body.vertices[edge.vertices[1]].source.clone(),
                            });
                        }
                        loops.push(FaceLoop {
                            outer: boundary.outer,
                            orientation: boundary.orientation,
                            edges,
                            vertex: boundary.vertex.map(|v| body.vertices[v].source.clone()),
                        });
                    }
                    faces.push(Face {
                        name: face.name.clone(),
                        surface: face.surface.clone(),
                        same_sense: face.same_sense,
                        loops,
                    });
                }
                shells.push(crate::SurfaceShell {
                    name: shell.name.clone(),
                    closed: shell.declared_closed,
                    faces,
                });
            }
            if body.declared_solid {
                file.solids.push(crate::Solid {
                    name: body.name.clone(),
                    faces: shells.remove(0).faces,
                });
            } else {
                file.surface_models.push(crate::SurfaceModel {
                    name: body.name.clone(),
                    shells,
                });
            }
        }
        Ok(file)
    }
    /// Grow one source vertex's tolerance without changing geometry or bindings.
    pub fn grow_vertex_tolerance(
        &mut self,
        body: &str,
        vertex: &str,
        required_mm: Rational,
        cause: &str,
        dec: &mut Decider,
    ) -> Result<bool, Error> {
        self.body_mut(body, dec)?
            .grow_vertex_tolerance(vertex, required_mm, cause, dec)
    }
    /// Grow one source edge and its endpoints without changing source bindings.
    pub fn grow_edge_tolerance(
        &mut self,
        body: &str,
        edge: &str,
        required_mm: Rational,
        cause: &str,
        dec: &mut Decider,
    ) -> Result<bool, Error> {
        self.body_mut(body, dec)?
            .grow_edge_tolerance(edge, required_mm, cause, dec)
    }
    fn body_mut(&mut self, name: &str, dec: &mut Decider) -> Result<&mut TolerantBody, Error> {
        dec.charge_work(self.bodies.len() as u64 + name.len() as u64)?;
        self.bodies
            .iter_mut()
            .find(|b| b.name == name)
            .ok_or_else(|| invalid("unknown body"))
    }
    /// Serializable state, verified by replay when reopened against the source.
    pub fn record(&self) -> ToleranceRecord {
        ToleranceRecord {
            format: "catallax-import-tolerances/1".into(),
            source_sha256: self.source_sha256.clone(),
            initial_mm: self.options.initial_mm.to_string(),
            maximum_mm: self.options.maximum_mm.to_string(),
            bodies: self
                .bodies
                .iter()
                .map(|body| BodyToleranceRecord {
                    name: body.name.clone(),
                    vertices: body
                        .vertices
                        .iter()
                        .map(|v| (v.source.name.clone(), v.tolerance_mm.to_string()))
                        .collect(),
                    edges: body
                        .edges
                        .iter()
                        .map(|e| (e.name.clone(), e.tolerance_mm.to_string()))
                        .collect(),
                    changes: body.changes.clone(),
                })
                .collect(),
        }
    }
}

/// Versioned persistence record. Geometry is reread from the matching source.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ToleranceRecord {
    /// Always catallax-import-tolerances/1.
    pub format: String,
    /// Exact source hash to which entity identities belong.
    pub source_sha256: String,
    /// Starting tolerance, in rational millimetres.
    pub initial_mm: String,
    /// Healing cap, in rational millimetres.
    pub maximum_mm: String,
    /// Source bodies and all their entity values and growth histories.
    pub bodies: Vec<BodyToleranceRecord>,
}

/// All current entity tolerances and the history needed to reproduce them.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BodyToleranceRecord {
    /// Source body identity.
    pub name: String,
    /// Source vertex names and current tolerances, in first-use order.
    pub vertices: Vec<(String, String)>,
    /// Source edge names and current tolerances, in first-use order.
    pub edges: Vec<(String, String)>,
    /// Complete ordered growth history, including endpoint propagation.
    pub changes: Vec<ToleranceChange>,
}

fn invalid(detail: &str) -> Error {
    Refusal::new("import.invalid-tolerances", detail).into()
}

fn budget<T>(f: impl FnOnce() -> Result<T, Error>) -> Result<T, Error> {
    if operation_budget_active() {
        return f();
    }
    let (result, work) = with_operation_budget(catallax_topology::numerical::DEFAULT_WORK, f);
    if work.exhausted {
        return Err(ImportRefusal::ResourceLimit.into());
    }
    result
}

/// Read supported source geometry without requiring exact incidence or closure.
/// No repair, validity tier or new geometry-reading support is implied.
pub fn read_step(
    bytes: &[u8],
    expected_sha256: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<TolerantFile, Error> {
    budget(|| read_inner(bytes, expected_sha256, options, dec))
}

fn read_inner(
    bytes: &[u8],
    expected_sha256: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<TolerantFile, Error> {
    dec.charge_work(bytes.len() as u64)?;
    let actual = crate::sha256_hex(bytes);
    if actual != expected_sha256 {
        return Err(ImportRefusal::HashMismatch {
            expected: expected_sha256.into(),
            actual,
        }
        .into());
    }
    if options.initial_mm.sign() == Sign::Negative
        || options.maximum_mm.sign() != Sign::Positive
        || options.initial_mm > options.maximum_mm
    {
        return Err(invalid(
            "initial tolerance must be nonnegative and no greater than a positive cap",
        ));
    }
    let file = crate::read_step_with_decider(bytes, dec)?;
    from_file(file, expected_sha256, options, dec)
}

fn from_file(
    file: StepFile,
    hash: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<TolerantFile, Error> {
    let mut bodies = Vec::new();
    for solid in &file.solids {
        let mut body = TolerantBody::new(&solid.name, true, &options);
        body.shell(&solid.name, true, &solid.faces, dec)?;
        bodies.push(body);
    }
    for model in &file.surface_models {
        let mut body = TolerantBody::new(&model.name, false, &options);
        for shell in &model.shells {
            body.shell(&shell.name, shell.closed, &shell.faces, dec)?;
        }
        bodies.push(body);
    }
    Ok(TolerantFile {
        source_sha256: hash.into(),
        mm_per_unit: file.mm_per_unit,
        options,
        bodies,
    })
}

/// Reopen a tolerance record only against its exact immutable source. Replaying
/// all requested changes verifies values, caps, causes and endpoint propagation.
/// Renumbered or changed source bytes never inherit the old entity bindings.
pub fn reopen(
    bytes: &[u8],
    record: &ToleranceRecord,
    dec: &mut Decider,
) -> Result<TolerantFile, Error> {
    budget(|| {
        charge_record(record, dec)?;
        if record.format != "catallax-import-tolerances/1" {
            return Err(invalid("unsupported tolerance record format"));
        }
        // charge_record prepays parsing and formatting before either occurs.
        // Records written by this API use canonical integers or fractions;
        // exponent syntax is excluded before it can expand a short input.
        let parse = |s: &str| {
            let value =
                Rational::parse(s).map_err(|_| invalid("invalid rational millimetre value"))?;
            if value.to_string() != s {
                return Err(invalid("noncanonical rational millimetre value"));
            }
            Ok(value)
        };
        let mut file = read_inner(
            bytes,
            &record.source_sha256,
            ToleranceOptions {
                initial_mm: parse(&record.initial_mm)?,
                maximum_mm: parse(&record.maximum_mm)?,
            },
            dec,
        )?;
        if file.bodies.len() != record.bodies.len() {
            return Err(invalid("record body inventory differs from source"));
        }
        for (body, saved) in file.bodies.iter_mut().zip(&record.bodies) {
            if body.name != saved.name {
                return Err(invalid("record body identity differs from source"));
            }
            for change in &saved.changes {
                dec.charge_work(1)?;
                if change.propagated_from_edge.is_some() {
                    continue;
                }
                let value = parse(&change.after_mm)?;
                match &change.entity {
                    ToleranceEntity::Vertex(name) => {
                        body.grow_vertex_tolerance(name, value, &change.cause, dec)?;
                    }
                    ToleranceEntity::Edge(name) => {
                        body.grow_edge_tolerance(name, value, &change.cause, dec)?;
                    }
                }
            }
        }
        if !matches_record(&file, record, dec)? {
            return Err(invalid(
                "record values or history do not match the replayed changes",
            ));
        }
        dec.check_deadline()?;
        Ok(file)
    })
}

fn charge_text(value: &str, dec: &mut Decider) -> Result<(), Error> {
    dec.charge_work(1_u64.saturating_add(value.len() as u64))?;
    Ok(())
}

fn charge_number(value: &str, dec: &mut Decider) -> Result<(), Error> {
    let size = value.len() as u64;
    // Conservative quadratic prepayment covers big-integer parsing, fraction
    // reduction, comparisons, decimal formatting and the replay's copies.
    dec.charge_work(
        1_u64.saturating_add(
            size.saturating_mul(size)
                .saturating_add(size)
                .saturating_mul(8),
        ),
    )?;
    if value.is_empty()
        || !value
            .bytes()
            .all(|c| c.is_ascii_digit() || matches!(c, b'-' | b'/'))
    {
        return Err(invalid(
            "saved millimetre values must be canonical integers or fractions",
        ));
    }
    Ok(())
}

fn charge_change(change: &ToleranceChange, dec: &mut Decider) -> Result<(), Error> {
    let name = match &change.entity {
        ToleranceEntity::Vertex(name) | ToleranceEntity::Edge(name) => name,
    };
    charge_text(name, dec)?;
    charge_text(&change.cause, dec)?;
    charge_number(&change.before_mm, dec)?;
    charge_number(&change.after_mm, dec)?;
    if let Some(edge) = &change.propagated_from_edge {
        charge_text(edge, dec)?;
    }
    Ok(())
}

fn charge_record(record: &ToleranceRecord, dec: &mut Decider) -> Result<(), Error> {
    charge_text(&record.format, dec)?;
    charge_text(&record.source_sha256, dec)?;
    charge_number(&record.initial_mm, dec)?;
    charge_number(&record.maximum_mm, dec)?;
    for body in &record.bodies {
        charge_text(&body.name, dec)?;
        for (name, value) in body.vertices.iter().chain(&body.edges) {
            charge_text(name, dec)?;
            charge_number(value, dec)?;
        }
        for change in &body.changes {
            charge_change(change, dec)?;
        }
    }
    Ok(())
}

fn matches_record(
    file: &TolerantFile,
    record: &ToleranceRecord,
    dec: &mut Decider,
) -> Result<bool, Error> {
    // Compare the replay directly. Do not allocate a second complete record or
    // hide inventory traversal, number formatting or log equality in derive Eq.
    for (body, saved) in file.bodies.iter().zip(&record.bodies) {
        dec.charge_work(1)?;
        if body.vertices.len() != saved.vertices.len()
            || body.edges.len() != saved.edges.len()
            || body.changes.len() != saved.changes.len()
        {
            return Ok(false);
        }
        let values = body
            .vertices
            .iter()
            .map(|v| (&v.source.name, &v.tolerance_mm))
            .chain(body.edges.iter().map(|e| (&e.name, &e.tolerance_mm)));
        for ((name, value), (saved_name, saved_value)) in
            values.zip(saved.vertices.iter().chain(&saved.edges))
        {
            charge_text(name, dec)?;
            // All values came from the prepaid initial value or growth history.
            // Reserve formatting explicitly too, using a digit-count upper
            // bound derived without first formatting the integer operands.
            let bits = value
                .numerator()
                .bit_length()
                .saturating_add(value.denominator().bit_length());
            dec.charge_work(bits.saturating_mul(bits).saturating_add(1))?;
            if name != saved_name || value.to_string() != *saved_value {
                return Ok(false);
            }
        }
        for (change, saved) in body.changes.iter().zip(&saved.changes) {
            charge_change(change, dec)?;
            if change != saved {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
