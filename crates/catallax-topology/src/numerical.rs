//! Numerical admission of finite edges and explicit parameter boundaries.
//! Exact constructors, incidence facts and checkpoints are not used for admission.

use std::collections::BTreeSet;

use catallax_numbers::{
    Ball, Decider, MarkLog, NumericalCheck, NumericalEvidence, NumericalRepair, NumericalSample,
    Rational, Refusal, Sign, operation_budget_active, with_operation_budget,
};
use catallax_shapes::{Curve, Point3, PointDef, SplineSupport, Surface};

use crate::{
    Body, BodyId, BodyKind, Coedge, Edge, EntityCounts, Face, Loop, Lump, Model, Sense, Shell,
    Vertex,
};

mod analytic;
mod atlas;
mod curved;
mod handover;
mod native;
pub use native::{NumericalNativeTrim, build_with_native_trims};
mod source;
pub use source::definition_from_model;
/// Editable finite spline sheets with trim, extension and normal thickening.
pub mod surface;

/// Shared default work limit for a numerical build or continuation.
pub const DEFAULT_WORK: u64 = 20_000_000;
const SUBDIVISIONS: i64 = 8;
const BITS: u32 = 128;

/// A retained vertex definition in millimetres.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct NumericalVertex {
    /// Stable key, unique among this definition's vertices.
    pub key: String,
    /// Position before any approximation by an evaluator.
    pub position: [Rational; 3],
}

/// A finite edge, with parameters in vertex order.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct NumericalEdge {
    /// Stable key, unique among edges.
    pub key: String,
    /// Retained line, finite circular arc or native B-spline carrier.
    pub curve: Curve,
    /// Start and end indices in the vertex table.
    pub vertices: [usize; 2],
    /// Parameters corresponding to those vertices. Descending ranges are allowed.
    pub parameters: [Rational; 2],
}

/// One edge use and its explicit finite trim on a face.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct NumericalTrim {
    /// Index in the edge table.
    pub edge: usize,
    /// Direction relative to the edge's vertex order.
    pub sense: Sense,
    /// U/v endpoints in traversal order. Patch trims interpolate linearly unless
    /// this use has a retained [`NumericalNativeTrim`].
    /// A native spline on a plane retains its projected curve between them.
    pub uv: [[Rational; 2]; 2],
}

/// One oriented surface and its finite parameter boundaries.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct NumericalFace {
    /// Stable key, unique among faces. It survives a surface replacement.
    pub key: String,
    /// A plane, finite cylinder or cone chart, or embedded rational patch.
    /// The trim need not be a patch row.
    pub surface: Surface,
    /// Material normal relative to the support normal.
    pub sense: Sense,
    /// Ordered finite loops: one simple outer loop followed by interior holes.
    /// Holes run opposite to the outer loop, remain disjoint and may not nest.
    /// Planar native-spline walks receive sampled parameter-polygon checks.
    pub loops: Vec<Vec<NumericalTrim>>,
}

/// Immutable generating data retained by a checked numerical body.
/// Indices reference these tables; stable keys identify entities across rebuilds.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct NumericalDefinition {
    /// Definition name for diagnostics.
    pub name: String,
    /// Shared vertices.
    pub vertices: Vec<NumericalVertex>,
    /// Shared finite edges.
    pub edges: Vec<NumericalEdge>,
    /// Faces of one connected shell. `build` requires closure; `build_surface`
    /// permits boundary edges without asserting enclosed material.
    pub faces: Vec<NumericalFace>,
}

/// Stated controls, fixed before a numerical build starts.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub struct NumericalOptions {
    /// Maximum accepted distance at each checked location, in millimetres.
    /// This is not a promised global deviation bound.
    pub tolerance_mm: Rational,
}

impl Default for NumericalOptions {
    fn default() -> Self {
        Self {
            tolerance_mm: Rational::parse("0.01").expect("constant"),
        }
    }
}

/// An editable numerical B-rep with its generating definition and actual checks.
///
/// Later numerical consumers may rely on the retained definition, stable keys,
/// stated tolerance, completed check records and ordered repair history. Each
/// continuation rebuilds and rechecks a new body; failure leaves this body intact.
/// No exact Incident facts are inferred from samples. The numerical build and
/// continuation paths return this type with [`NumericalEvidence`], never an exact
/// or bounded answer. There is no conversion to an exact evaluation checkpoint.
///
/// This guarantee applies to these producers, not to arbitrary caller code.
/// The shared generic `Answer` and `BoundedEvidence` constructors still trust
/// their callers. An explicit caller can mislabel a numerical body as exact or
/// a sampled value as a complete bound. Preventing that misuse requires a separate
/// change to producer-owned result authority; this API does not claim that fence.
///
/// Acceptance examples for this first path:
/// - A structurally closed, consistently oriented shell with simple finite u/v
///   loops and all endpoint/trim samples within tolerance returns this type.
/// - Zero sampled residual still returns a numerical body, never an exact proof.
/// - A sampled maximum of 1/4500 mm passes a 0.01 mm tolerance; 1/45 mm refuses.
/// - A missing coedge, crossing checked u/v polygon, reversed loop, missing trim,
///   unsupported carrier, or exhausted work refuses, even if other samples pass.
/// - A vertex repair records its before/after positions and point displacement;
///   no repair is an explicit empty history. Historical distances do not become
///   a bound on the current body.
///
/// Nine equally spaced locations per coedge, including its endpoints, are checked
/// against its retained face parameters. Planar circular/spline edges use projected curve
/// samples; other ordinary trims use their linear u/v segment. Retained native
/// UV curves use endpoints and eight subdivisions of each clipped 3D/UV knot
/// span at the same finite native parameter. Their endpoint correspondence is
/// checked in millimetres against the declared loop endpoints. Every edge endpoint is
/// checked against its vertex. These are pointwise distance enclosures, not bounds
/// between samples. Projected or retained native UV walks keep their sampled
/// simplicity checks. Complete planar circular boundaries additionally receive
/// continuous projected-disk clearance checks, not chord-polygon pair tests.
/// Other curved boundary containment remains sampled and explicitly unresolved
/// between nodes.
/// The checks do not establish whole-body separation, global self-intersection
/// freedom or a whole-body approximation bound. They grant no STEP export authority.
/// Approximate export requires a numerical writer and its own output checks.
/// Existing exact builds and exporters retain their interfaces and behavior.
///
/// ```compile_fail
/// use catallax_numbers::Answer;
/// use catallax_topology::numerical::NumericalBody;
/// fn promote(body: NumericalBody) -> Answer<NumericalBody> { body.into() }
/// ```
#[derive(Clone, Debug)]
pub struct NumericalBody {
    source: NumericalDefinition,
    source_native_trims: Vec<NumericalNativeTrim>,
    native_trims: Vec<NumericalNativeTrim>,
    history: Vec<NumericalHistory>,
    definition: NumericalDefinition,
    options: NumericalOptions,
    model: Model,
    body: BodyId,
    evidence: NumericalEvidence,
}

/// A previous checked state. Its checks apply to this definition and tolerance only.
#[derive(Clone, Debug)]
pub struct NumericalHistory {
    /// The generating definition at this stage.
    pub definition: NumericalDefinition,
    /// Native parameter curves checked at this stage.
    pub native_trims: Vec<NumericalNativeTrim>,
    /// The controls used at this stage.
    pub options: NumericalOptions,
    /// The checks and cumulative repair history at this stage.
    pub evidence: NumericalEvidence,
}

impl NumericalBody {
    /// Original generating input, including endpoints removed by explicit sewing.
    pub fn source_definition(&self) -> &NumericalDefinition {
        &self.source
    }
    /// Earlier checked states in order. These checks do not describe the current body.
    pub fn history(&self) -> &[NumericalHistory] {
        &self.history
    }
    /// Rechecks an edited definition under explicitly stated current controls.
    /// Prior checks and repairs remain historical; no old sample certifies this body.
    pub fn rebuild(
        &self,
        definition: NumericalDefinition,
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        self.rebuild_with_native_trims(definition, self.native_trims.clone(), options, dec)
    }
    fn inherit(&mut self, previous: &Self) {
        self.source = previous.source.clone();
        self.source_native_trims = previous.source_native_trims.clone();
        self.history = previous.history.clone();
        self.history.push(NumericalHistory {
            definition: previous.definition.clone(),
            native_trims: previous.native_trims.clone(),
            options: previous.options.clone(),
            evidence: previous.evidence.clone(),
        });
    }

    /// Retained generating data, including the finite parameter trims.
    pub fn definition(&self) -> &NumericalDefinition {
        &self.definition
    }
    /// The tolerance used for this result and its continuations.
    pub fn options(&self) -> &NumericalOptions {
        &self.options
    }
    /// Actual checks and repair history, with no implicit complete bound.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
    /// Stored topological kind, without exposing the model to exact consumers.
    pub fn kind(&self) -> BodyKind {
        self.model.body(self.body).kind
    }
    /// Counts of the stored B-rep. No mesh stands in for this body.
    pub fn counts(&self) -> EntityCounts {
        self.model.counts(self.body)
    }

    /// Diagnostics from the unchanged exact topology checker. Missing exact
    /// incidence remains missing even when every numerical sample is zero.
    /// The raw model stays private so this result cannot enter exact consumers.
    pub fn exact_topology_faults(&self) -> Vec<crate::TopologyFault> {
        self.model.check(self.body)
    }

    /// Replaces one support and rechecks its retained trims and the complete body.
    /// This is a parameter edit, not a repair. The old result remains unchanged.
    pub fn replace_surface(
        &self,
        face: &str,
        surface: Surface,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        guarded(|| {
            let mut definition = self.definition.clone();
            let selected = definition
                .faces
                .iter_mut()
                .find(|f| f.key == face)
                .ok_or_else(|| bad("unknown-face", format!("unknown face {face}")))?;
            selected.surface = surface;
            self.rebuild(definition, self.options.clone(), dec)
        })
    }

    /// Moves one vertex under the unchanged tolerance, records its actual point
    /// displacement, and rechecks every required condition before returning it.
    pub fn repair_vertex(
        &self,
        vertex: &str,
        position: [Rational; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        guarded(|| {
            let mut definition = self.definition.clone();
            let selected = definition
                .vertices
                .iter_mut()
                .find(|v| v.key == vertex)
                .ok_or_else(|| bad("unknown-vertex", format!("unknown vertex {vertex}")))?;
            precision_allowed(dec)?;
            dec.charge_work(1)?;
            let before = selected.position.clone();
            let movement = distance(&point_box(&before), &point_box(&position));
            if movement > self.options.tolerance_mm {
                return Err(bad(
                    "repair-tolerance-exceeded",
                    "vertex movement exceeds the unchanged tolerance",
                ));
            }
            selected.position = position.clone();
            let mut repairs = self.evidence.repairs().to_vec();
            if before != position {
                repairs.push(NumericalRepair::VertexMoved {
                    vertex: vertex.into(),
                    before,
                    after: position,
                    distance_upper_mm: movement,
                });
            }
            let mut next = build_inner_kind(
                self.kind(),
                definition,
                self.native_trims.clone(),
                self.options.clone(),
                repairs,
                dec,
            )?;
            next.inherit(self);
            Ok(next)
        })
    }
}

/// Builds a numerical result under the caller's active work limit, or 20M work
/// units when no limit is active. Point enclosures use 128 bits; a lower caller
/// precision cap refuses. It never retries by increasing a tolerance or cap.
/// See [`NumericalBody`] for the completed checks and their precise limits.
pub fn build(
    definition: NumericalDefinition,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    guarded(|| build_inner(definition, Vec::new(), options, Vec::new(), dec))
}

/// Builds an open or closed surface shell under stated numerical tolerances.
/// Boundary edges have one use; shared edges require two opposite uses.
/// The existing surface topology checker validates boundary cycles, connectivity,
/// vertex fans and Euler characteristic. No material or solid volume is inferred.
pub fn build_surface(
    definition: NumericalDefinition,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    build_surface_with_native_trims(definition, Vec::new(), options, dec)
}

/// Builds a surface shell with retained native parameter curves.
/// Its sample checks and limitations match `build_with_native_trims`.
pub fn build_surface_with_native_trims(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    guarded(|| {
        build_inner_kind(
            BodyKind::Surface,
            definition,
            native_trims,
            options,
            Vec::new(),
            dec,
        )
    })
}

fn guarded<T>(run: impl FnOnce() -> Result<T, Refusal>) -> Result<T, Refusal> {
    if operation_budget_active() {
        return run();
    }
    let (result, work) = with_operation_budget(DEFAULT_WORK, run);
    if work.exhausted {
        Err(Refusal::new(
            "numbers.operation-budget",
            "numerical build exhausted its work allowance",
        ))
    } else {
        result
    }
}
fn bad(code: &'static str, detail: impl Into<String>) -> Refusal {
    let code = match code {
        "definition" => "numerical.invalid-definition",
        "unsupported" => "numerical.unsupported-carrier",
        "trim" => "numerical.invalid-trim",
        "unknown-face" => "numerical.unknown-face",
        "unknown-vertex" => "numerical.unknown-vertex",
        "repair-tolerance-exceeded" => "numerical.repair-tolerance-exceeded",
        _ => unreachable!(),
    };
    Refusal::new(code, detail)
}
fn precision_allowed(dec: &Decider) -> Result<(), Refusal> {
    if dec.policy().max_precision_bits < BITS {
        return Err(Refusal::new(
            "numbers.precision-limit",
            "numerical admission requires 128-bit point enclosures",
        ));
    }
    Ok(())
}
// Rational sample addresses stay exact even when their binary enclosure
// extends beyond a finite spline endpoint.
fn curve_at(
    curve: &Curve,
    parameter: &Rational,
    dec: &Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    match curve {
        Curve::BSpline(curve) => curve
            .evaluate_box([parameter, parameter], BITS, dec)
            .map(Some),
        _ => curve.evaluate(&Ball::from_rational(parameter, BITS), dec),
    }
}

fn point_box(p: &[Rational; 3]) -> [Ball; 3] {
    p.each_ref().map(|x| Ball::from_rational(x, BITS))
}
fn distance(a: &[Ball; 3], b: &[Ball; 3]) -> Rational {
    let d = std::array::from_fn::<_, 3, _>(|i| &a[i] - &b[i]);
    (d[0].square() + d[1].square() + d[2].square())
        .sqrt()
        .expect("sum of squares")
        .hi()
}
fn keys<'a>(values: impl Iterator<Item = &'a str>, dec: &mut Decider) -> Result<(), Refusal> {
    let mut seen = BTreeSet::new();
    for key in values {
        dec.charge_work(1)?;
        if key.trim().is_empty() || !seen.insert(key) {
            return Err(bad(
                "definition",
                "entity keys must be nonempty and unique within each table",
            ));
        }
    }
    Ok(())
}

fn build_inner(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    repairs: Vec<NumericalRepair>,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    build_inner_kind(
        BodyKind::Solid,
        definition,
        native_trims,
        options,
        repairs,
        dec,
    )
}

fn build_inner_kind(
    kind: BodyKind,
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    repairs: Vec<NumericalRepair>,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    precision_allowed(dec)?;
    dec.charge_work(1)?;
    if definition.name.trim().is_empty() || options.tolerance_mm.sign() != Sign::Positive {
        return Err(bad(
            "definition",
            "a named definition and positive tolerance are required",
        ));
    }
    keys(definition.vertices.iter().map(|v| v.key.as_str()), dec)?;
    keys(definition.edges.iter().map(|e| e.key.as_str()), dec)?;
    keys(definition.faces.iter().map(|f| f.key.as_str()), dec)?;
    let native = native::Index::new(&definition, &native_trims, dec)?;
    let mut model = Model::new();
    let vertices: Vec<_> = definition
        .vertices
        .iter()
        .map(|v| {
            let [x, y, z] = v.position.clone();
            let point = model.geom_mut().add_point(PointDef::Explicit(Point3::new(
                x.into(),
                y.into(),
                z.into(),
            )));
            model.add_vertex(Vertex { point })
        })
        .collect();
    let mut edges = Vec::new();
    for edge in &definition.edges {
        dec.charge_work(1)?;
        match &edge.curve {
            Curve::Line(line) => {
                catallax_shapes::Line::new(line.origin.clone(), line.dir.clone(), dec)?;
            }
            Curve::BSpline(curve) => {
                let [lo, hi] = curve.domain();
                if edge.parameters.iter().any(|t| t < lo || t > hi) {
                    return Err(bad(
                        "definition",
                        "edge parameters leave the retained spline domain",
                    ));
                }
            }
            Curve::Circle(circle) => analytic::validate_circle(
                circle,
                &edge.parameters,
                edge.vertices[0] == edge.vertices[1],
                dec,
            )?,
            _ => {
                return Err(bad(
                    "unsupported",
                    "numerical admission requires a line or finite native B-spline edge",
                ));
            }
        }
        if edge.vertices.iter().any(|i| *i >= vertices.len())
            || edge.parameters[0] == edge.parameters[1]
        {
            return Err(bad(
                "definition",
                "edge has an invalid vertex index or empty parameter interval",
            ));
        }
        let curve = model.geom_mut().add_curve(edge.curve.clone());
        edges.push(model.add_edge(Edge {
            curve,
            sense: if edge.parameters[0] < edge.parameters[1] {
                Sense::Same
            } else {
                Sense::Reversed
            },
            start: Some(vertices[edge.vertices[0]]),
            end: Some(vertices[edge.vertices[1]]),
        }));
    }
    let mut faces = Vec::new();
    for (face_index, face) in definition.faces.iter().enumerate() {
        dec.charge_work(1)?;
        if !analytic::supports(&face.surface) {
            return Err(bad(
                "unsupported",
                "numerical admission supports planes, finite cylinder or cone charts and embedded rational patches",
            ));
        }
        analytic::validate_surface(&face.surface, dec)?;
        if let Surface::Plane(plane) = &face.surface {
            catallax_shapes::Plane::new(
                plane.origin.clone(),
                plane.normal.clone(),
                plane.ref_dir.clone(),
                dec,
            )?;
        }
        if face.loops.is_empty() {
            return Err(bad("trim", "this admission path requires an outer loop"));
        }
        // Chart-union containment is checked at sampled native UV locations,
        // so an atlas face accepts no straight parameter trim, and its charts
        // must fill their box so that the filled loop stays on them.
        if let Surface::BSpline(SplineSupport::RationalAtlas(source)) = &face.surface {
            atlas::check_filled(source, dec)?;
            for (row, boundary) in face.loops.iter().enumerate() {
                for use_index in 0..boundary.len() {
                    dec.charge_work(1)?;
                    if native.get(face_index, row, use_index).is_none() {
                        return Err(bad(
                            "trim",
                            "an atlas face needs a native UV curve on every edge use",
                        ));
                    }
                }
            }
        }
        // Validate every interior boundary before building its coedges.
        // Duplicated, nested or crossing boundaries remain invalid trims.
        if face.loops.len() > 1 {
            native::check_loops(
                face_index,
                face,
                &definition.edges,
                &native,
                &options.tolerance_mm,
                dec,
            )?;
            check_distinct_boundary_vertices(face, &definition.edges, dec)?;
        }
        let mut loops = Vec::new();
        for boundary in &face.loops {
            let mut uses = Vec::new();
            for trim in boundary {
                dec.charge_work(1)?;
                let edge = *edges
                    .get(trim.edge)
                    .ok_or_else(|| bad("definition", "trim names an invalid edge index"))?;
                uses.push(model.add_coedge(Coedge {
                    edge,
                    sense: trim.sense,
                }));
            }
            loops.push(model.add_loop(Loop::Coedges(uses)));
        }
        let surface = model.geom_mut().add_surface(face.surface.clone());
        faces.push(model.add_face(Face {
            surface,
            sense: face.sense,
            loops,
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind,
        lumps: vec![lump],
        wires: Vec::new(),
        marks: MarkLog::new(),
    });
    if let Some(fault) = crate::check::structure(&model, body).into_iter().next() {
        return Err(Refusal::new(fault.code, fault.detail));
    }
    let counts = model.counts(body);
    if counts.vertices != definition.vertices.len() || counts.edges != definition.edges.len() {
        return Err(bad(
            "definition",
            "unused vertices or edges are not part of this result",
        ));
    }
    let mut checks = vec![NumericalCheck::Passed {
        name: "topology.structure".into(),
        scope: definition.name.clone(),
    }];
    for edge in &definition.edges {
        let mut samples = Vec::new();
        for i in 0..2 {
            dec.charge_work(1)?;
            let t = &edge.parameters[i];
            let value = crate::numerical::curve_at(&edge.curve, t, dec)?
                .ok_or_else(|| bad("unsupported", "edge has no parameter evaluator"))?;
            samples.push(NumericalSample {
                parameters: vec![t.clone()],
                distance_upper_mm: distance(
                    &value,
                    &point_box(&definition.vertices[edge.vertices[i]].position),
                ),
            });
        }
        checks.push(NumericalCheck::Samples {
            name: "edge.vertex-distance".into(),
            scope: edge.key.clone(),
            samples,
        });
    }
    for (face_index, face) in definition.faces.iter().enumerate() {
        if face.loops.len() == 1 {
            native::check_loops(
                face_index,
                face,
                &definition.edges,
                &native,
                &options.tolerance_mm,
                dec,
            )?;
        }
        checks.extend(native::join_samples(
            face_index,
            face,
            &definition.edges,
            &native,
            dec,
        )?);
        let sampled = native.has_face(face_index) || analytic::projected(face, &definition.edges);
        checks.push(NumericalCheck::Passed {
            name: if sampled {
                "trim.sampled-simple-oriented-parameter-loop"
            } else {
                "trim.simple-oriented-parameter-loop"
            }
            .into(),
            scope: face.key.clone(),
        });
        if matches!(
            face.surface,
            Surface::BSpline(SplineSupport::RationalAtlas(_))
        ) {
            checks.push(NumericalCheck::Passed {
                name: "surface.atlas-charts-fill-box".into(),
                scope: face.key.clone(),
            });
            checks.push(NumericalCheck::Passed {
                name: "trim.sampled-atlas-chart-union".into(),
                scope: face.key.clone(),
            });
        }
        if face.loops.len() > 1 {
            // This family necessarily ran analytic::circle_clearance during
            // loop validation. Native curves on rational patches did not.
            let planar_circle_clearance = matches!(face.surface, Surface::Plane(_))
                && face.loops.iter().any(|boundary| {
                    analytic::complete_circle(boundary, &definition.edges).is_some()
                });
            checks.push(NumericalCheck::Passed {
                name: if planar_circle_clearance {
                    "trim.sampled-boundaries-with-continuous-circle-clearance"
                } else if sampled {
                    "trim.sampled-hole-contained-disjoint"
                } else {
                    "trim.hole-contained-disjoint"
                }
                .into(),
                scope: face.key.clone(),
            });
            if planar_circle_clearance
                && face.loops.iter().skip(1).any(|boundary| {
                    analytic::complete_circle(boundary, &definition.edges).is_some()
                })
            {
                checks.push(NumericalCheck::Passed {
                    name: "trim.continuous-circle-hole-clearance".into(),
                    scope: face.key.clone(),
                });
            }
            if planar_circle_clearance
                && analytic::complete_circle(&face.loops[0], &definition.edges).is_some()
            {
                checks.push(NumericalCheck::Passed {
                    name: "trim.continuous-circle-outer-clearance".into(),
                    scope: face.key.clone(),
                });
            }
            checks.push(NumericalCheck::Passed {
                name: "trim.boundary-endpoints-disjoint".into(),
                scope: face.key.clone(),
            });
        }
        for (index, (boundary_index, use_index, trim)) in face
            .loops
            .iter()
            .enumerate()
            .flat_map(|(b, row)| row.iter().enumerate().map(move |(u, t)| (b, u, t)))
            .enumerate()
        {
            let edge = &definition.edges[trim.edge];
            if let Some(curve) = native.get(face_index, boundary_index, use_index) {
                let (samples, endpoints) = native::samples(face, trim, edge, curve, dec)?;
                checks.push(NumericalCheck::Samples {
                    name: "trim.native-endpoint-position".into(),
                    scope: format!("face {:?}, edge {:?}, use {index}", face.key, edge.key),
                    samples: endpoints,
                });
                checks.push(NumericalCheck::Samples {
                    name: "trim.edge-surface-distance".into(),
                    scope: format!("face {:?}, edge {:?}, use {index}", face.key, edge.key),
                    samples,
                });
                continue;
            }
            let [a, b] = if trim.sense == Sense::Same {
                edge.parameters.each_ref()
            } else {
                [&edge.parameters[1], &edge.parameters[0]]
            };
            let mut samples = Vec::new();
            for i in 0..=SUBDIVISIONS {
                dec.charge_work(1)?;
                let s = Rational::from_ratio(i, SUBDIVISIONS).expect("positive denominator");
                let t = a + &((&s) * &(b - a));
                let mut uv: [Rational; 2] = std::array::from_fn(|axis| {
                    &trim.uv[0][axis] + &(&s * &(&trim.uv[1][axis] - &trim.uv[0][axis]))
                });
                let p = crate::numerical::curve_at(&edge.curve, &t, dec)?
                    .ok_or_else(|| bad("unsupported", "edge has no parameter evaluator"))?;
                if i != 0 && i != SUBDIVISIONS && analytic::projected_edge(face, edge) {
                    uv = curved::parameters(&face.surface, &p, dec)?;
                }
                let q = if matches!(
                    face.surface,
                    Surface::BSpline(SplineSupport::RationalPatch(_))
                ) {
                    native::point(
                        &face.surface,
                        &uv.clone().map(catallax_numbers::Real::from),
                        dec,
                    )?
                } else {
                    face.surface.evaluate(
                        &Ball::from_rational(&uv[0], BITS),
                        &Ball::from_rational(&uv[1], BITS),
                        dec,
                    )?
                };
                samples.push(NumericalSample {
                    parameters: vec![t, uv[0].clone(), uv[1].clone()],
                    distance_upper_mm: distance(&p, &q),
                });
            }
            checks.push(NumericalCheck::Samples {
                name: "trim.edge-surface-distance".into(),
                scope: format!("face {:?}, edge {:?}, use {index}", face.key, edge.key),
                samples,
            });
        }
    }
    let evidence = NumericalEvidence::new(options.tolerance_mm.clone(), checks, repairs)?;
    dec.check_deadline()?;
    Ok(NumericalBody {
        source: definition.clone(),
        source_native_trims: native_trims.clone(),
        native_trims,
        history: Vec::new(),
        definition,
        options,
        model,
        body,
        evidence,
    })
}

fn orient(a: &[Rational; 2], b: &[Rational; 2], c: &[Rational; 2]) -> Rational {
    (&b[0] - &a[0]) * (&c[1] - &a[1]) - (&b[1] - &a[1]) * (&c[0] - &a[0])
}
fn on_segment(a: &[Rational; 2], b: &[Rational; 2], p: &[Rational; 2]) -> bool {
    orient(a, b, p).is_zero()
        && (0..2).all(|i| {
            p[i] >= a[i].clone().min(b[i].clone()) && p[i] <= a[i].clone().max(b[i].clone())
        })
}
fn intersects(a: &[Rational; 2], b: &[Rational; 2], c: &[Rational; 2], d: &[Rational; 2]) -> bool {
    if on_segment(a, b, c) || on_segment(a, b, d) || on_segment(c, d, a) || on_segment(c, d, b) {
        return true;
    }
    let signs = [
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    ]
    .map(|v| v.sign());
    signs.iter().all(|s| *s != Sign::Zero) && signs[0] != signs[1] && signs[2] != signs[3]
}
fn check_loop(
    face: &NumericalFace,
    boundary: &[NumericalTrim],
    hole: bool,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let n = boundary.len();
    if n < 3 {
        return Err(bad(
            "trim",
            "a parameter polygon needs at least three edges",
        ));
    }
    let mut area = Rational::zero();
    for (i, trim) in boundary.iter().enumerate() {
        dec.charge_work(1)?;
        let [a, b] = &trim.uv;
        if a == b || *b != boundary[(i + 1) % n].uv[0] {
            return Err(bad(
                "trim",
                format!("{} has an empty or open parameter edge", face.key),
            ));
        }
        if let Surface::BSpline(SplineSupport::RationalPatch(patch)) = &face.surface {
            for p in [a, b] {
                if (0..2).any(|axis| {
                    p[axis] < patch.domain()[axis][0] || p[axis] > patch.domain()[axis][1]
                }) {
                    return Err(bad(
                        "trim",
                        "parameter trim leaves the retained patch domain",
                    ));
                }
            }
        }
        // An atlas face evaluates its whole retained source, which can hold
        // other components. Its trims must stay on its own charts.
        if let Surface::BSpline(SplineSupport::RationalAtlas(atlas)) = &face.surface {
            for p in [a, b] {
                dec.charge_work(atlas.charts().len() as u64)?;
                if !atlas.charts().iter().any(|chart| {
                    let domain = chart.patch().domain();
                    (0..2).all(|axis| p[axis] >= domain[axis][0] && p[axis] <= domain[axis][1])
                }) {
                    return Err(bad(
                        "trim",
                        "parameter trim leaves the charts of its atlas face",
                    ));
                }
            }
        }
        area = area + (&a[0] * &b[1] - &a[1] * &b[0]);
        let next = &boundary[(i + 1) % n].uv[1];
        if on_segment(a, b, next) || on_segment(b, next, a) {
            return Err(bad("trim", "adjacent parameter edges overlap or reverse"));
        }
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            dec.charge_work(1)?;
            let [c, d] = &boundary[j].uv;
            if intersects(a, b, c, d) {
                return Err(bad("trim", "parameter loop crosses or touches itself"));
            }
        }
    }
    let wanted = if (face.sense == Sense::Same) != hole {
        Sign::Positive
    } else {
        Sign::Negative
    };
    if area.sign() != wanted {
        return Err(bad(
            "trim",
            "parameter loop orientation disagrees with material normal",
        ));
    }
    Ok(())
}

// UV polygons alone do not establish separation of their topological walks:
// explicit sewing can unify endpoint IDs while leaving the UV trims unchanged.
fn check_distinct_boundary_vertices(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut previous = BTreeSet::new();
    for boundary in &face.loops {
        let mut walked = BTreeSet::new();
        for trim in boundary {
            dec.charge_work(1)?;
            let edge = edges
                .get(trim.edge)
                .ok_or_else(|| bad("definition", "trim names an invalid edge index"))?;
            // Sense reverses the traversal, but both endpoints are walked.
            walked.extend(edge.vertices);
        }
        if !previous.is_disjoint(&walked) {
            return Err(bad("trim", "outer and hole boundaries share an endpoint"));
        }
        previous.extend(walked);
    }
    Ok(())
}

fn inside_loop(
    point: &[Rational; 2],
    boundary: &[NumericalTrim],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut inside = false;
    for edge in boundary {
        dec.charge_work(1)?;
        let [a, b] = &edge.uv;
        if (a[1] > point[1]) != (b[1] > point[1])
            && (orient(a, b, point).sign() == Sign::Positive) == (b[1] > a[1])
        {
            inside = !inside;
        }
    }
    Ok(inside)
}

fn check_loops(
    face: &NumericalFace,
    cleared_circles: Option<(&NumericalFace, &[NumericalEdge])>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for (i, boundary) in face.loops.iter().enumerate() {
        check_loop(face, boundary, i != 0, dec)?;
    }
    for (i, hole) in face.loops.iter().enumerate().skip(1) {
        let hole_circle = cleared_circles.is_some_and(|(source, edges)| {
            analytic::complete_circle(&source.loops[i], edges).is_some()
        });
        for (j, other) in face.loops[..i].iter().enumerate() {
            let other_circle = cleared_circles.is_some_and(|(source, edges)| {
                analytic::complete_circle(&source.loops[j], edges).is_some()
            });
            if !hole_circle && !other_circle {
                for a in other {
                    for b in hole {
                        dec.charge_work(1)?;
                        if intersects(&a.uv[0], &a.uv[1], &b.uv[0], &b.uv[1]) {
                            return Err(bad("trim", "hole crosses or touches another boundary"));
                        }
                    }
                }
            }
            // Inscribed chord polygons can cross even when two real circular
            // rims are disjoint. Those pairs already passed disk clearance.
            if other_circle && (hole_circle || j == 0) {
                continue;
            }
            // Simple disjoint boundaries need one strict witness per component.
            let contained = inside_loop(&hole[0].uv[0], other, dec)?;
            if j == 0 {
                if !contained {
                    return Err(bad("trim", "hole is not strictly inside the outer loop"));
                }
            } else if contained || inside_loop(&other[0].uv[0], hole, dec)? {
                return Err(bad("trim", "interior holes overlap or nest"));
            }
        }
    }
    Ok(())
}

/// Sews explicitly paired endpoint vertices, then admits a separate numerical body.
/// Each pair is `[keeper, removed]` in the input vertex table. Pairs must be
/// disjoint; there is no transitive, nearest-neighbour or automatic gap search.
/// Removed endpoints are relocated onto the keeper only within the stated local
/// tolerance, recorded with N01's actual point-movement history, and their edge
/// references are unified. The removed key remains in history; the keeper key
/// identifies the retained vertex. Carriers and trims are never moved.
///
/// All checks and repairs share one active request allowance. Rebuilding the
/// sewn shell checks topology, every finite loop and every required sample.
/// Excessive gaps, crossing loops, collapsed edges or invalid pairs return a
/// refusal with no partial body. The borrowed input is never mutated.
/// Even a zero-gap result is numerical, with no new proof or export authority.
pub fn sew(
    input: &NumericalDefinition,
    pairs: &[[usize; 2]],
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    guarded(|| {
        precision_allowed(dec)?;
        if options.tolerance_mm.sign() != Sign::Positive {
            return Err(bad("definition", "sewing needs a positive local tolerance"));
        }
        keys(input.vertices.iter().map(|v| v.key.as_str()), dec)?;
        let mut definition = input.clone();
        let mut seen = BTreeSet::new();
        let mut remove = BTreeSet::new();
        let mut repairs = Vec::new();
        let mut mapping: Vec<_> = (0..definition.vertices.len()).collect();
        for &[keeper, removed] in pairs {
            dec.charge_work(1)?;
            if keeper >= mapping.len()
                || removed >= mapping.len()
                || keeper == removed
                || !seen.insert(keeper)
                || !seen.insert(removed)
            {
                return Err(bad(
                    "definition",
                    "sewing pairs must name disjoint distinct vertices",
                ));
            }
            let before = definition.vertices[removed].position.clone();
            let after = definition.vertices[keeper].position.clone();
            let movement = distance(&point_box(&before), &point_box(&after));
            if movement > options.tolerance_mm {
                return Err(bad(
                    "repair-tolerance-exceeded",
                    "sewing gap exceeds the local tolerance",
                ));
            }
            if before != after {
                repairs.push(NumericalRepair::VertexMoved {
                    vertex: definition.vertices[removed].key.clone(),
                    before,
                    after,
                    distance_upper_mm: movement,
                });
            }
            mapping[removed] = keeper;
            remove.insert(removed);
        }
        let mut compact = vec![0; mapping.len()];
        let mut retained = Vec::new();
        for (index, vertex) in definition.vertices.iter().enumerate() {
            dec.charge_work(1)?;
            if !remove.contains(&index) {
                compact[index] = retained.len();
                retained.push(vertex.clone());
            }
        }
        for edge in &mut definition.edges {
            dec.charge_work(1)?;
            for vertex in &mut edge.vertices {
                let target = *mapping
                    .get(*vertex)
                    .ok_or_else(|| bad("definition", "edge names an invalid vertex"))?;
                *vertex = compact[target];
            }
        }
        definition.vertices = retained;
        let mut body = build_inner(definition, Vec::new(), options, repairs, dec)?;
        body.source = input.clone();
        let mut checks = body.evidence.checks().to_vec();
        for &[keeper, removed] in pairs {
            dec.charge_work(1)?;
            checks.push(NumericalCheck::Passed {
                name: "sewing.explicit-vertex-pair".into(),
                scope: format!(
                    "{} -> {}",
                    input.vertices[removed].key, input.vertices[keeper].key
                ),
            });
        }
        body.evidence = NumericalEvidence::new(
            body.options.tolerance_mm.clone(),
            checks,
            body.evidence.repairs().to_vec(),
        )?;
        dec.check_deadline()?;
        Ok(body)
    })
}
