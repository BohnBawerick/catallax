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
mod edge_gap;
pub use edge_gap::{DerivedCopyAudit, audit_derived_copies};
/// The edge-gap checker: the largest distance between an edge's two copies.
pub mod gap;
mod handover;
mod import_gap;
pub use import_gap::{ImportGapOutcome, ImportGapRow};
mod native;
pub use native::{NumericalNativeTrim, build_with_native_trims, build_with_sampled_carriers};
mod source;
pub use source::definition_from_model;
/// Editable finite spline sheets with trim, extension and normal thickening.
pub mod surface;
/// Tolerance rules T0 to T8 (decision Q7, section 4.3).
pub mod tolerance;
pub use tolerance::{ToleranceReport, ToleranceWarning};

#[cfg(test)]
mod known_answers;
#[cfg(test)]
mod segment_tests;

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

/// What the face copy of an edge use was made from (decision Q8, section 4.4
/// of the Checkpoint 2 decisions).
///
/// Every edge use has two copies of its curve: the 3D curve and the face
/// copy, which is the u/v endpoints of a [`NumericalTrim`] or the u/v curve of
/// a [`NumericalNativeTrim`]. The tag says which one is authoritative.
/// Admission measures the gap of an independent copy and records the bound of
/// a derived one.
///
/// The tag is not written to a stored body. A body that is read back has
/// every copy `Independent`, as an imported one has, and a rebuild from the
/// recipe sets the tags again.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum CopyTag {
    /// The two copies were made separately: fitted, traced, projected, read
    /// from a file, or edited on one side.
    #[default]
    Independent,
    /// One copy was made from the other by a rule whose gap has a known
    /// bound.
    Derived {
        /// The copy the other one was made from.
        from: CopySource,
        /// The rule that made it.
        rule: CopyRule,
        /// The largest distance between the two copies that the rule allows,
        /// in millimetres. Zero for a rule that is exact.
        bound: Rational,
    },
}

/// The authoritative copy of a derived edge use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CopySource {
    /// The 3D curve. The face copy was computed from it.
    Edge3d,
    /// The face's surface. The 3D curve is a curve of that surface, and the
    /// face copy names it.
    Carrier,
}

/// The rule that made a derived copy.
///
/// A line or a circle of an analytic surface, embedded by the surface's
/// chart, has no rule here: its u/v ends are rounded, so it is
/// [`CopyTag::Independent`] and admission measures it by closed form (owner
/// ruling of 20 November 2026 on section 4.4).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CopyRule {
    /// The 3D curve is a row of the blend surface's control net, and the
    /// face copy is that row's parameter line.
    BlendRow,
    /// A u/v curve composed exactly from another copy.
    ExactComposition,
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
    /// What these endpoints were made from. A retained
    /// [`NumericalNativeTrim`] on the same use carries its own tag, and that
    /// one is read.
    #[transfer(skip)]
    pub tag: CopyTag,
}

/// One oriented surface and its finite parameter boundaries.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct NumericalFace {
    /// Stable key, unique among faces. It survives a surface replacement.
    pub key: String,
    /// A finite analytic chart or retained regular B-spline carrier.
    /// The trim need not follow a carrier parameter row.
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

impl NumericalDefinition {
    /// Sets every edge use back to [`CopyTag::Independent`].
    ///
    /// A tag is a statement by the code that built both copies of a use. Code
    /// that takes a definition over from a body and then changes a curve, a
    /// surface, a parameter range or a u/v end calls this first, unless it
    /// applies one exact operation to both copies and that operation keeps the
    /// stated bound. Admission then measures the uses again.
    pub fn reset_copy_tags(&mut self) {
        for trim in self
            .faces
            .iter_mut()
            .flat_map(|face| face.loops.iter_mut().flatten())
        {
            trim.tag = CopyTag::Independent;
        }
    }
}

/// Stated controls, fixed before a numerical build starts.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub struct NumericalOptions {
    /// The largest tolerance an edge use may store, and the largest distance
    /// accepted at every other checked location, in millimetres. An edge use
    /// stores 1.1 times its measured gap (rule T4 of the tolerance rules).
    /// This is not a promised global deviation bound.
    pub tolerance_mm: Rational,
    /// The largest distance from an edge end to its vertex, where the class
    /// gives vertices their own cap: rule T3's 0.02 mm for imported vertices
    /// ([`tolerance::import_options`]). `None`: vertices share
    /// `tolerance_mm`. A stored body keeps it in its tolerance report.
    #[transfer(skip)]
    pub vertex_tolerance_mm: Option<Rational>,
}

impl Default for NumericalOptions {
    fn default() -> Self {
        Self {
            tolerance_mm: Rational::parse("0.01").expect("constant"),
            vertex_tolerance_mm: None,
        }
    }
}

impl NumericalOptions {
    /// The cap of an edge end's distance to its vertex.
    pub fn vertex_cap_mm(&self) -> &Rational {
        self.vertex_tolerance_mm
            .as_ref()
            .unwrap_or(&self.tolerance_mm)
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
/// - A measured gap of 1/4500 mm passes a 0.01 mm tolerance; 1/45 mm refuses.
/// - A missing coedge, crossing checked u/v polygon, reversed loop, missing trim,
///   unsupported carrier, or exhausted work refuses, even if other samples pass.
/// - A vertex repair records its before/after positions and point displacement;
///   no repair is an explicit empty history. Historical distances do not become
///   a bound on the current body.
///
/// Every edge use is measured by [`gap::measure`]: the largest distance
/// between the edge's 3D curve and its face copy at the same parameter. A
/// retained native UV curve is the face copy where one exists. A circle or
/// spline on a plane is measured against the plane, and its two declared end
/// values are measured as well. Other ordinary trims use their linear u/v
/// segment. The record of a use holds one sample: the parameter of the
/// largest gap and the tolerance stored for it, 1.1 times the gap plus its
/// measurement error. A stored tolerance above the body's limit refuses. A
/// value within its binary64 rounding error of the limit is decided on a
/// 128-bit enclosure at that parameter; nothing else raises the precision. This is a search for the maximum, not a proof
/// of one: see `docs/milestones/3.0A-edge-gap-checker.md` for what it can miss.
/// Native UV curves also have their endpoint correspondence checked in
/// millimetres against the declared loop endpoints. Every edge endpoint is
/// checked against its vertex. Projected or retained native UV walks keep
/// their sampled simplicity checks. Complete planar circular boundaries
/// additionally receive continuous projected-disk clearance checks, not
/// chord-polygon pair tests. Other curved boundary containment remains
/// sampled and explicitly unresolved between nodes.
/// Explicitly declared numerical atlas carriers use a 9 by 9 chart-box grid
/// and nine UV-walk locations per coedge for face-local checks. Their evidence
/// says numerical samples, never continuous filled-box or curved-loop proof.
/// Edge uses on such a carrier are measured like every other use and keep
/// their own labelled check names.
/// Other face-local checks keep their original path; incidence stays unchanged.
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
    source_sampled_carriers: BTreeSet<String>,
    sampled_carriers: BTreeSet<String>,
    native_trims: Vec<NumericalNativeTrim>,
    history: Vec<NumericalHistory>,
    definition: NumericalDefinition,
    options: NumericalOptions,
    model: Model,
    body: BodyId,
    evidence: NumericalEvidence,
    tolerances: ToleranceReport,
    import_gaps: Vec<ImportGapRow>,
}

/// A previous checked state. Its checks apply to this definition and tolerance only.
#[derive(Clone, Debug)]
pub struct NumericalHistory {
    /// The generating definition at this stage.
    pub definition: NumericalDefinition,
    /// Native parameter curves checked at this stage.
    pub native_trims: Vec<NumericalNativeTrim>,
    /// Numerical carriers admitted with sampled face-local checks.
    pub sampled_carriers: BTreeSet<String>,
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
    ///
    /// The retained native curves come from this body, not from the caller.
    /// Copy tags are kept only where both copies are unchanged; see
    /// [`NumericalBody::rebuild_with_native_trims`].
    pub fn rebuild(
        &self,
        definition: NumericalDefinition,
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        self.rebuild_with_native_trims(definition, self.native_trims.clone(), options, dec)
    }

    /// Sets back to [`CopyTag::Independent`] every derived tag of an edit that
    /// this body cannot show still holds (decision Q8, section 4.4 of the
    /// Checkpoint 2 decisions: any change to one copy resets it, unless the
    /// same exact operation applies to both).
    ///
    /// The tags of an edited definition are the caller's, and a caller can
    /// clone a tagged definition and change one copy without touching the
    /// tag. So a use keeps a derived tag only when this body has a use at the
    /// same face, loop and place with the same tag, the same face surface, the
    /// same sense and u/v ends, an edge with the same curve and parameter
    /// range, and the same face copy: the linear u/v ends on both sides, or
    /// the same retained u/v curve on both sides, so dropping or replacing a
    /// retained curve resets the use. A retained u/v curve keeps its tag only
    /// when, in addition, this
    /// body retains the same curve with the same tag for that use, and the
    /// use still states that tag after the reset above, so resetting a use
    /// resets its curve as well. Identity is [`Curve::same_definition`] and
    /// [`Surface::same_definition`]: a false answer only costs a measurement.
    fn keep_unchanged_copy_tags(
        &self,
        definition: &mut NumericalDefinition,
        native_trims: &mut [NumericalNativeTrim],
    ) {
        let before = &self.definition;
        let changed: Vec<(usize, usize, usize)> = definition
            .faces
            .iter()
            .enumerate()
            .flat_map(|(f, face)| {
                face.loops.iter().enumerate().flat_map(move |(b, row)| {
                    row.iter().enumerate().map(move |(u, trim)| (f, b, u, trim))
                })
            })
            .filter(|(f, b, u, trim)| {
                let key = definition.faces[*f].key.as_str();
                trim.tag != CopyTag::Independent
                    && (same_use(before, definition, *f, *b, *u)
                        .is_none_or(|old| old.tag != trim.tag)
                        || !same_face_copy(
                            native_at(&self.native_trims, key, *b, *u),
                            native_at(&*native_trims, key, *b, *u),
                        ))
            })
            .map(|(f, b, u, _)| (f, b, u))
            .collect();
        for (f, b, u) in changed {
            definition.faces[f].loops[b][u].tag = CopyTag::Independent;
        }
        for binding in native_trims {
            if binding.tag == CopyTag::Independent {
                continue;
            }
            let kept = definition
                .faces
                .iter()
                .position(|face| face.key == binding.face)
                .filter(|&f| {
                    same_use(before, definition, f, binding.boundary, binding.use_index).is_some()
                        && definition.faces[f].loops[binding.boundary][binding.use_index].tag
                            == binding.tag
                })
                .is_some()
                && self.native_trims.iter().any(|old| {
                    old.face == binding.face
                        && old.boundary == binding.boundary
                        && old.use_index == binding.use_index
                        && old.edge == binding.edge
                        && old.tag == binding.tag
                        && old.curve.same_definition(&binding.curve)
                });
            if !kept {
                binding.tag = CopyTag::Independent;
            }
        }
    }
    /// Refits the face copy of every measured edge use whose stored
    /// tolerance is over `above_mm`, aiming under it (rule T8: "refit edges
    /// above the receiver tolerance"). A refit is kept only where the checker
    /// reads a stored tolerance under `above_mm`; the body is then admitted
    /// again as a continuation of this one. Returns the body and one row per
    /// kept refit. A body with nothing to lower is returned as it is.
    pub fn refit_face_copies(
        &self,
        above_mm: &Rational,
        dec: &mut Decider,
    ) -> Result<(Self, Vec<RefitRow>), Refusal> {
        guarded(|| {
            let records = tolerance::use_records(&self.evidence);
            let mut chosen = Vec::new();
            let mut index = 0;
            for (face_index, face) in self.definition.faces.iter().enumerate() {
                for (boundary, row) in face.loops.iter().enumerate() {
                    for use_index in 0..row.len() {
                        dec.charge_work(1)?;
                        let (stored, derived) = &records[index];
                        index += 1;
                        if !derived && stored > above_mm {
                            chosen.push(((face_index, boundary, use_index), stored.clone()));
                        }
                    }
                }
            }
            let target = above_mm.to_f64_nearest();
            match self.refit_uses(&chosen, target, dec)? {
                None => Ok((self.clone(), Vec::new())),
                Some((native_trims, rows)) => {
                    let mut next = build_inner_kind(
                        self.kind(),
                        self.definition.clone(),
                        native_trims,
                        self.sampled_carriers.clone(),
                        self.options.clone(),
                        self.evidence.repairs().to_vec(),
                        dec,
                    )?;
                    next.inherit(self, dec)?;
                    Ok((next, rows))
                }
            }
        })
    }

    /// Refits the chosen uses and returns the native curves with every refit
    /// that brings its use's stored tolerance under `target_mm`, or `None`
    /// when none does.
    #[allow(clippy::type_complexity)]
    fn refit_uses(
        &self,
        chosen: &[((usize, usize, usize), Rational)],
        target_mm: f64,
        dec: &mut Decider,
    ) -> Result<Option<(Vec<NumericalNativeTrim>, Vec<RefitRow>)>, Refusal> {
        if chosen.is_empty() {
            return Ok(None);
        }
        let index = native::Index::new(
            &self.definition,
            &self.native_trims,
            &self.sampled_carriers,
            dec,
        )?;
        let mut native_trims = self.native_trims.clone();
        let mut rows = Vec::new();
        for ((face_index, boundary, use_index), before) in chosen {
            let face = &self.definition.faces[*face_index];
            let trim = &face.loops[*boundary][*use_index];
            let edge = &self.definition.edges[trim.edge];
            let native = index.get(*face_index, *boundary, *use_index);
            let Some(found) = edge_gap::refit_use(face, trim, edge, native, target_mm, dec)? else {
                continue;
            };
            let after = Rational::from_f64_exact(found.stored_mm)
                .ok_or_else(|| bad("trim", "a refit gap is not a finite number"))?;
            // A refit is kept where it reaches its target: a copy that
            // only moves the reading in its last bits is no lower gap.
            if after >= *before || after.to_f64_nearest() > target_mm {
                continue;
            }
            native_trims.retain(|binding| {
                !(binding.face == face.key
                    && binding.boundary == *boundary
                    && binding.use_index == *use_index)
            });
            native_trims.push(NumericalNativeTrim {
                face: face.key.clone(),
                boundary: *boundary,
                use_index: *use_index,
                edge: edge.key.clone(),
                curve: found.curve,
                tag: CopyTag::Independent,
            });
            rows.push(RefitRow {
                face: face.key.clone(),
                edge: edge.key.clone(),
                boundary: *boundary,
                use_index: *use_index,
                before_mm: before.clone(),
                after_mm: after,
                spans: found.spans,
            });
        }
        Ok((!rows.is_empty()).then_some((native_trims, rows)))
    }

    /// Rule T7 after a continuation: an edge the input already had whose
    /// tolerance rose more than ten times above its own earlier value has the
    /// face copies of its measured uses refitted, aiming under ten times that
    /// value, and the body is admitted again where that lowers a stored
    /// tolerance. What still rose stays a warning. A new edge was just fitted
    /// by the operation that made it, so it is not fitted twice; its growth
    /// over the input's largest tolerance is reported. Every use was measured
    /// afresh by the continuation's own admission (rule T4, never stacked),
    /// which is where a tolerance comes down when its gap does.
    fn lower_growth(self, previous: &Self, dec: &mut Decider) -> Result<Self, Refusal> {
        let grown: BTreeSet<&str> = self
            .tolerances
            .warnings
            .iter()
            .filter_map(|warning| match warning {
                ToleranceWarning::Growth {
                    kind: tolerance::EntityKind::Edge,
                    key,
                    same_entity: true,
                    ..
                } => Some(key.as_str()),
                _ => None,
            })
            .collect();
        if grown.is_empty() {
            return Ok(self);
        }
        let records = tolerance::use_records(&self.evidence);
        let earlier = tolerance::by_key(&previous.tolerances, dec)?;
        let mut chosen = Vec::new();
        let mut index = 0;
        let mut target = f64::INFINITY;
        for (face_index, face) in self.definition.faces.iter().enumerate() {
            for (boundary, row) in face.loops.iter().enumerate() {
                for (use_index, trim) in row.iter().enumerate() {
                    dec.charge_work(1)?;
                    let (stored, derived) = &records[index];
                    index += 1;
                    let key = self.definition.edges[trim.edge].key.as_str();
                    if !derived && grown.contains(key) {
                        chosen.push(((face_index, boundary, use_index), stored.clone()));
                        dec.charge_work(1)?;
                        let input = earlier
                            .get(&(tolerance::EntityKind::Edge, key))
                            .map(|value| (*value).clone())
                            .unwrap_or_default();
                        target = target.min(
                            (input * Rational::from_i64(tolerance::GROWTH_FACTOR)).to_f64_nearest(),
                        );
                    }
                }
            }
        }
        let Some((native_trims, _)) = self.refit_uses(&chosen, target, dec)? else {
            return Ok(self);
        };
        let mut next = build_inner_kind(
            self.kind(),
            self.definition.clone(),
            native_trims,
            self.sampled_carriers.clone(),
            self.options.clone(),
            self.evidence.repairs().to_vec(),
            dec,
        )?;
        next.inherit(previous, dec)?;
        Ok(next)
    }

    fn inherit(&mut self, previous: &Self, dec: &mut Decider) -> Result<(), Refusal> {
        // Rule T7: this body was measured again by its own admission, with
        // rule T4's margin on fresh measurements only, so no tolerance is
        // carried over or stacked. What rose more than ten times above its
        // input is reported with its absolute values.
        let growth = tolerance::growth_warnings(&previous.tolerances, &self.tolerances, dec)?;
        self.tolerances.warnings.extend(growth);
        self.source = previous.source.clone();
        self.source_native_trims = previous.source_native_trims.clone();
        self.source_sampled_carriers = previous.source_sampled_carriers.clone();
        self.history = previous.history.clone();
        self.history.push(NumericalHistory {
            definition: previous.definition.clone(),
            native_trims: previous.native_trims.clone(),
            sampled_carriers: previous.sampled_carriers.clone(),
            options: previous.options.clone(),
            evidence: previous.evidence.clone(),
        });
        Ok(())
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
    /// What the tolerance rules found at admission: every vertex and edge
    /// tolerance, rule T1's floor and frame size, and the warnings of rules
    /// T1, T5 and T7. See [`tolerance`].
    pub fn tolerances(&self) -> &ToleranceReport {
        &self.tolerances
    }
    /// What the gap rule for imports read (decision Q9): one row per edge
    /// use that [`admit_import`] read geometrically, with the two gaps and
    /// what was done. Empty for an own build, for a continuation and for a
    /// body that was read back.
    pub fn import_gaps(&self) -> &[ImportGapRow] {
        &self.import_gaps
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
            // One copy of every edge use on this face changed and the other
            // did not, so none of them is derived any more.
            for trim in selected.loops.iter_mut().flatten() {
                trim.tag = CopyTag::Independent;
            }
            let mut native_trims = self.native_trims.clone();
            for binding in native_trims.iter_mut().filter(|b| b.face == face) {
                binding.tag = CopyTag::Independent;
            }
            self.rebuild_with_native_trims(definition, native_trims, self.options.clone(), dec)
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
                self.sampled_carriers.clone(),
                self.options.clone(),
                repairs,
                dec,
            )?;
            next.inherit(self, dec)?;
            next.lower_growth(self, dec)
        })
    }
}

/// One edge use whose face copy a refit replaced (rules T7 and T8).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RefitRow {
    /// Key of the face of the use.
    pub face: String,
    /// Key of the edge.
    pub edge: String,
    /// Loop index on the face.
    pub boundary: usize,
    /// Use index in the loop.
    pub use_index: usize,
    /// Stored tolerance before, in millimetres.
    pub before_mm: Rational,
    /// Stored tolerance of the refit copy as the checker read it.
    pub after_mm: Rational,
    /// Spans of the refit copy.
    pub spans: usize,
}

/// The refusal of a stored tolerance over the body's cap.
const TOLERANCE_EXCEEDED: &str = "numbers.numerical-tolerance-exceeded";

/// Admits an imported definition under decision Q9 and rule T3.
///
/// Every entity is measured; the caller's `options` hold the cap: the
/// automatic 0.01 mm, or a value the user opted in to, at most 0.1 mm. A
/// file's declared uncertainty plays no part here: rule T3 makes it a hint
/// only.
///
/// Decision Q9, the gap rule for imports: an edge use whose same-parameter
/// gap is much larger than its geometric gap has copies that run at
/// different speeds. Its face copy is given the parameter of its 3D curve,
/// and the tolerance stored is that of the new copy's same-parameter gap; no
/// tolerance is raised to cover a difference of parameter. See
/// [`NumericalBody::import_gaps`] for what was read and done, and
/// `docs/milestones/3.0A-import-gap-rule.md` for the method.
///
/// Rule T3: where a stored tolerance is still over the cap, the face copies
/// over it are refitted from their 3D curves first ("above that, refit
/// first"), on a probe admitted at the opt-in limit, and the definition is
/// admitted again under the caller's cap with the refitted copies. Only what
/// is still over the cap refuses, and the refusal says the user may opt in.
///
/// A body that is admitted as imported and is not admitted with a
/// reparametrised copy is returned as imported.
pub fn admit_import(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    if options.tolerance_mm > tolerance::import_opt_in_limit_mm() {
        return Err(Refusal::new(
            "numerical.import-opt-in",
            format!(
                "an import cap of {} mm is over the {} mm a user may opt in to (rule T3)",
                options.tolerance_mm,
                tolerance::import_opt_in_limit_mm()
            ),
        ));
    }
    let options = tolerance::import_options(options.tolerance_mm);
    let first = build_with_native_trims(
        definition.clone(),
        native_trims.clone(),
        options.clone(),
        dec,
    );
    let known = match &first {
        Ok(body) => Some(tolerance::use_records(&body.evidence)),
        Err(refusal) if refusal.code == TOLERANCE_EXCEEDED => None,
        Err(_) => return first,
    };
    let (reparametrised, mut rows) =
        guarded(|| import_gap::settle(&definition, &native_trims, known.as_deref(), dec))?;
    let with_new_copies = reparametrised.map(|native_trims| {
        admit_under_cap(definition.clone(), native_trims, options.clone(), dec)
    });
    let admitted = import_gap::settled_body(first, with_new_copies, &mut rows, |refusal| {
        refit_under_cap(definition, native_trims, options, refusal, dec)
    });
    admitted.map(|mut body| {
        body.import_gaps = rows;
        body
    })
}

/// Admits a definition under the caller's cap, with rule T3's refit where a
/// stored tolerance is over it.
fn admit_under_cap(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    match build_with_native_trims(
        definition.clone(),
        native_trims.clone(),
        options.clone(),
        dec,
    ) {
        Err(refusal) if refusal.code == TOLERANCE_EXCEEDED => {
            refit_under_cap(definition, native_trims, options, refusal, dec)
        }
        other => other,
    }
}

/// Rule T3 after `refusal` at the cap: refits the face copies over the cap
/// on a probe admitted at the opt-in limit, and admits the definition again
/// under the caller's cap with the refitted copies.
fn refit_under_cap(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    refusal: Refusal,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    let probe = build_with_native_trims(
        definition.clone(),
        native_trims,
        tolerance::import_options(tolerance::import_opt_in_limit_mm()),
        dec,
    )
    .map_err(|_| refusal.clone())?;
    let (refitted, rows) = probe.refit_face_copies(&options.tolerance_mm, dec)?;
    if rows.is_empty() {
        return Err(opt_in_refusal(refusal));
    }
    build_with_native_trims(definition, refitted.native_trims, options, dec).map_err(|again| {
        if again.code == TOLERANCE_EXCEEDED {
            opt_in_refusal(again)
        } else {
            again
        }
    })
}

fn opt_in_refusal(mut refusal: Refusal) -> Refusal {
    refusal.detail = format!(
        "{}; refitted first; the user may opt in to a cap of up to {} mm (rule T3)",
        refusal.detail,
        tolerance::import_opt_in_limit_mm()
    );
    refusal
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
    guarded(|| {
        build_inner(
            definition,
            Vec::new(),
            BTreeSet::new(),
            options,
            Vec::new(),
            dec,
        )
    })
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
    build_surface_with_sampled_carriers(definition, native_trims, BTreeSet::new(), options, dec)
}

/// Builds a numerical surface shell while retaining declared sampled carriers.
pub fn build_surface_with_sampled_carriers(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    sampled_carriers: BTreeSet<String>,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    guarded(|| {
        build_inner_kind(
            BodyKind::Surface,
            definition,
            native_trims,
            sampled_carriers,
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
/// The use of `before` at the same face key, loop and place as use `use_index`
/// of loop `boundary` of face `face` of `after`, when both of its copies are
/// the same definition there: the face surface, the sense and u/v ends, and
/// the edge's curve and parameter range. Tags are not compared here.
fn same_use<'a>(
    before: &'a NumericalDefinition,
    after: &NumericalDefinition,
    face: usize,
    boundary: usize,
    use_index: usize,
) -> Option<&'a NumericalTrim> {
    let face = after.faces.get(face)?;
    let trim = face.loops.get(boundary)?.get(use_index)?;
    let edge = after.edges.get(trim.edge)?;
    let face_before = before.faces.iter().find(|f| f.key == face.key)?;
    let old = face_before.loops.get(boundary)?.get(use_index)?;
    let edge_before = before.edges.get(old.edge)?;
    (face_before.surface.same_definition(&face.surface)
        && old.sense == trim.sense
        && old.uv == trim.uv
        && edge_before.curve.same_definition(&edge.curve)
        && edge_before.parameters == edge.parameters)
        .then_some(old)
}
/// The retained u/v curve bound to use `use_index` of loop `boundary` of face
/// `face`, if any.
fn native_at<'a>(
    bindings: &'a [NumericalNativeTrim],
    face: &str,
    boundary: usize,
    use_index: usize,
) -> Option<&'a NumericalNativeTrim> {
    bindings
        .iter()
        .find(|b| b.face == face && b.boundary == boundary && b.use_index == use_index)
}
/// Whether a use's face copy is represented the same way before and after an
/// edit: by the linear u/v ends on both sides, or by the same retained u/v
/// curve on both sides. Adding, removing or replacing a retained curve changes
/// the face copy even when the u/v ends stay.
fn same_face_copy(
    before: Option<&NumericalNativeTrim>,
    after: Option<&NumericalNativeTrim>,
) -> bool {
    match (before, after) {
        (None, None) => true,
        (Some(a), Some(b)) => a.edge == b.edge && a.curve.same_definition(&b.curve),
        _ => false,
    }
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
    sampled_carriers: BTreeSet<String>,
    options: NumericalOptions,
    repairs: Vec<NumericalRepair>,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    build_inner_kind(
        BodyKind::Solid,
        definition,
        native_trims,
        sampled_carriers,
        options,
        repairs,
        dec,
    )
}

fn build_inner_kind(
    kind: BodyKind,
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    sampled_carriers: BTreeSet<String>,
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
    for key in &sampled_carriers {
        let face = definition
            .faces
            .iter()
            .find(|face| &face.key == key)
            .ok_or_else(|| bad("definition", "sampled carrier names an unknown face"))?;
        if !matches!(
            face.surface,
            Surface::BSpline(SplineSupport::RationalAtlas(_))
        ) {
            return Err(bad(
                "unsupported",
                "sampled carrier coverage requires a retained rational atlas",
            ));
        }
    }
    let native = native::Index::new(&definition, &native_trims, &sampled_carriers, dec)?;
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
                "numerical admission needs a supported finite analytic chart or retained regular spline carrier",
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
        // Every atlas edge needs native UV data. Exact faces also establish
        // complete chart-box coverage; explicitly numerical carriers record
        // sampled coverage only.
        if let Surface::BSpline(SplineSupport::RationalAtlas(source)) = &face.surface {
            if sampled_carriers.contains(&face.key) {
                atlas::check_sampled(source, dec)?;
            } else {
                atlas::check_filled(source, dec)?;
            }
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
    for (edge_index, edge) in definition.edges.iter().enumerate() {
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
            name: if native.numerical_pair(edge_index) {
                "edge.numerical-pair-vertex-distance"
            } else {
                "edge.vertex-distance"
            }
            .into(),
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
            name: if sampled_carriers.contains(&face.key) {
                "trim.numerical-carrier-sampled-simple-oriented-parameter-loop-9-per-use"
            } else if sampled {
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
                name: if sampled_carriers.contains(&face.key) {
                    "surface.sampled-atlas-chart-box-coverage-9x9"
                } else {
                    "surface.atlas-charts-fill-box"
                }
                .into(),
                scope: face.key.clone(),
            });
            checks.push(NumericalCheck::Passed {
                name: if sampled_carriers.contains(&face.key) {
                    "trim.numerical-carrier-sampled-atlas-chart-union-9-per-use"
                } else {
                    "trim.sampled-atlas-chart-union"
                }
                .into(),
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
                name: if sampled_carriers.contains(&face.key) {
                    "trim.numerical-carrier-sampled-hole-contained-disjoint"
                } else if planar_circle_clearance {
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
            let numerical_pair = native.numerical_pair(trim.edge);
            let curve = native.get(face_index, boundary_index, use_index);
            let scope = format!("face {:?}, edge {:?}, use {index}", face.key, edge.key);
            if let Some(curve) = curve {
                checks.push(NumericalCheck::Samples {
                    name: if numerical_pair {
                        "trim.numerical-pair-native-endpoint-position"
                    } else {
                        "trim.native-endpoint-position"
                    }
                    .into(),
                    scope: scope.clone(),
                    samples: native::endpoints(face, trim, edge, curve, dec)?,
                });
            }
            // One record per use. A measured use says where its two copies
            // are furthest apart and the tolerance stored for that gap. A
            // derived use says the bound its rule states, under its own name,
            // so that a bound is never read as a measurement.
            let (derived, sample) = edge_gap::stored_sample(
                face,
                trim,
                edge,
                curve,
                // A retained native curve is the face copy of its use, so its
                // tag is the one that counts.
                native
                    .tag(face_index, boundary_index, use_index)
                    .unwrap_or(&trim.tag),
                &options.tolerance_mm,
                dec,
            )?;
            checks.push(NumericalCheck::Samples {
                name: match (numerical_pair, derived) {
                    (false, false) => "trim.edge-surface-distance",
                    (true, false) => "trim.numerical-pair-edge-surface-distance",
                    (false, true) => "trim.edge-surface-derived-bound",
                    (true, true) => "trim.numerical-pair-edge-surface-derived-bound",
                }
                .into(),
                scope,
                samples: vec![sample],
            });
        }
    }
    // Rule T3 caps a vertex apart from its edges and faces. Every record of an
    // edge end's distance to its vertex is held to the vertex cap with T4's
    // margin, and every other record to the class cap.
    let evidence_tolerance = if let Some(vertex) = &options.vertex_tolerance_mm {
        tolerance::check_by_class(&checks, &options.tolerance_mm, vertex)?;
        vertex.clone().max(options.tolerance_mm.clone())
    } else {
        options.tolerance_mm.clone()
    };
    let evidence = NumericalEvidence::new(evidence_tolerance, checks, repairs)?;
    let tolerances = tolerance_report(&definition, &native, &evidence, &options, dec)?;
    dec.check_deadline()?;
    Ok(NumericalBody {
        source: definition.clone(),
        source_native_trims: native_trims.clone(),
        source_sampled_carriers: sampled_carriers.clone(),
        sampled_carriers,
        native_trims,
        history: Vec::new(),
        definition,
        options,
        model,
        body,
        evidence,
        tolerances,
        import_gaps: Vec::new(),
    })
}

/// The tolerance rules on an admitted definition: every entity's tolerance
/// over rule T1's floor, rule T1's size warning, and rule T5, which refuses a
/// pair whose tolerances add up to more than half their distance. That is the
/// local term the class cap (`options.tolerance_mm`) is compared with: the
/// cap that applies is the smaller of the two.
fn tolerance_report(
    definition: &NumericalDefinition,
    native: &native::Index<'_>,
    evidence: &NumericalEvidence,
    options: &NumericalOptions,
    dec: &mut Decider,
) -> Result<ToleranceReport, Refusal> {
    let frame = tolerance::BodyFrame::of(definition, dec)?;
    let floor_mm = frame.floor_mm();
    dec.charge_work(
        (definition.vertices.len() + definition.edges.len() + evidence.checks().len()) as u64,
    )?;
    let entities = tolerance::entity_tolerances(
        definition,
        evidence,
        &floor_mm,
        options.vertex_tolerance_mm.as_ref(),
    )?;
    let mut warnings = Vec::new();
    if frame.extent_mm > tolerance::body_size_warning_mm() {
        warnings.push(ToleranceWarning::LargeBody {
            extent_mm: frame.extent_mm.clone(),
        });
    }
    warnings.extend(
        tolerance::pair_rule_check(definition, native, &entities, &frame.centre, dec)?
            .into_iter()
            .map(ToleranceWarning::Pair),
    );
    Ok(ToleranceReport {
        coordinate_size_mm: frame.coordinate_size_mm,
        floor_mm,
        cap_mm: options.tolerance_mm.clone(),
        vertex_cap_mm: options.vertex_tolerance_mm.clone(),
        entities,
        warnings,
    })
}

fn orient(a: &[Rational; 2], b: &[Rational; 2], c: &[Rational; 2]) -> Rational {
    (&b[0] - &a[0]) * (&c[1] - &a[1]) - (&b[1] - &a[1]) * (&c[0] - &a[0])
}
fn in_segment_box(a: &[Rational; 2], b: &[Rational; 2], p: &[Rational; 2]) -> bool {
    (0..2).all(|i| &p[i] >= (&a[i]).min(&b[i]) && &p[i] <= (&a[i]).max(&b[i]))
}
fn on_segment(a: &[Rational; 2], b: &[Rational; 2], p: &[Rational; 2]) -> bool {
    in_segment_box(a, b, p) && orient(a, b, p).is_zero()
}
fn intersects(a: &[Rational; 2], b: &[Rational; 2], c: &[Rational; 2], d: &[Rational; 2]) -> bool {
    // Disjoint exact boxes certify separation without reducing any determinant
    // fractions. Native pcurve samples can have large unrelated denominators.
    if (0..2)
        .any(|i| (&a[i]).min(&b[i]) > (&c[i]).max(&d[i]) || (&c[i]).min(&d[i]) > (&a[i]).max(&b[i]))
    {
        return false;
    }
    let signs = [
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    ]
    .map(|v| v.sign());
    if (signs[0] == Sign::Zero && in_segment_box(a, b, c))
        || (signs[1] == Sign::Zero && in_segment_box(a, b, d))
        || (signs[2] == Sign::Zero && in_segment_box(c, d, a))
        || (signs[3] == Sign::Zero && in_segment_box(c, d, b))
    {
        return true;
    }
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
        if let Surface::BSpline(source) = &face.surface
            && !matches!(
                source,
                SplineSupport::RationalPatch(_) | SplineSupport::RationalAtlas(_)
            )
        {
            let domain = source.carrier().domain();
            if [a, b].iter().any(|point| {
                (0..2).any(|axis| point[axis] < *domain[axis][0] || point[axis] > *domain[axis][1])
            }) {
                return Err(bad(
                    "trim",
                    "parameter trim leaves its retained spline carrier domain",
                ));
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
    check_boundaries(face, &face.loops, cleared_circles, dec)
}

fn check_boundaries(
    face: &NumericalFace,
    boundaries: &[Vec<NumericalTrim>],
    cleared_circles: Option<(&NumericalFace, &[NumericalEdge])>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for (i, boundary) in boundaries.iter().enumerate() {
        check_loop(face, boundary, i != 0, dec)?;
    }
    for (i, hole) in boundaries.iter().enumerate().skip(1) {
        let hole_circle = cleared_circles.is_some_and(|(source, edges)| {
            analytic::complete_circle(&source.loops[i], edges).is_some()
        });
        for (j, other) in boundaries[..i].iter().enumerate() {
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
        let mut body = build_inner(
            definition,
            Vec::new(),
            BTreeSet::new(),
            options,
            repairs,
            dec,
        )?;
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
