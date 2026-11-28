//! Piece 5 of catallax: the boolean engine.
//!
//! Unite, subtract and intersect through five reusable stages: intersect,
//! imprint, classify, keep, sew. Plus exact containment of rational points
//! and rational points inside faces.
//!
//! [`check_body`] explicitly checks structure and opposite-sheet overlaps on
//! matching carriers, with conservative spatial pruning before exact proofs.
//!
//! The contract is `docs/pieces/05-boolean-engine.md`.
//!
//! Face imprints are oriented halfedge arrangements. Equal tangents are
//! ordered by exact curvature from curve equations. Face membership uses
//! closed-form boundary predicates, exact surface probes or certified trace
//! homotopies in rational charts; containment counts transverse rational rays.
//! Rational face witnesses follow exact chart paths strictly left of a
//! boundary arc; Bernstein coefficient signs certify separation from every edge.
//!
//! `cargo bench -p catallax-boolean --bench speed` prints non-gating stage timings,
//! decision shares, containment timings and reachable stored-number growth.
//! Add `-- --case 10` to run one numbered fixture, or `-- --case 44` for a chain.

mod arrangement;
mod bounds;
mod broad_phase;
mod carrier_bounds;
mod containment;
mod curve_pair;
mod face_queries;
mod fast_region;
mod geometry;
mod germ;
mod graph;
mod imprinting;
pub mod numerical;
mod point_identity;
mod points;
mod probe_path;
mod probed_region;
mod project;
mod region;
mod sewing;
mod spline_fillet;
mod stages;
mod subgap;
mod torus_region;
mod trim_bounds;
mod validation;

#[cfg(test)]
mod tests;

pub use arrangement::{PatchSelection, arrange, trim_patch};
pub use stages::{
    Classification, FaceClass, Imprinted, IntersectionGraph, SectionEdge, Selection, classify,
    imprint, intersect, keep, make_maximal, sew,
};

use std::collections::BTreeMap;

use catallax_numbers::{Answer, Answered, Decider, DecisionStats, MarkLog, Outcome, Refusal};
use catallax_shapes::{Point3, RPoint3};
use catallax_topology::{Body, BodyId, BodyKind, EdgeId, FaceId, Model};

pub use face_queries::PreparedFaceQueries;
pub use subgap::SUB_TOLERANCE_MERGE;
pub use validation::check_body;

/// A boolean operation.
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum BoolOp {
    /// `A ∪ B`.
    Unite,
    /// `A − B`.
    Subtract,
    /// `A ∩ B`.
    Intersect,
}

/// The five stages of the engine.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Stage {
    /// Find where the boundaries meet.
    Intersect,
    /// Split faces and edges along the meetings.
    Imprint,
    /// Label face pieces against the other body.
    Classify,
    /// Choose the face pieces to keep.
    Keep,
    /// Join kept pieces into the result body.
    Sew,
}

/// What one stage did.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct StageRecord {
    /// The stage.
    pub stage: Stage,
    /// Stage counter increases; `max_bits` is the decider's running maximum at stage end.
    pub stats: DecisionStats,
    /// Marks the stage added.
    pub marks: usize,
}

/// Per-stage record of one boolean.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct StageReport {
    /// Records in stage order.
    pub stages: Vec<StageRecord>,
}

/// A finished boolean.
#[derive(Clone, Debug)]
pub struct BooleanAnswer {
    /// The new body.
    pub body: BodyId,
    /// Marks of this operation followed by those of both operands.
    pub marks: MarkLog,
    /// Per-stage record.
    pub report: StageReport,
}

/// Immediate operand IDs contributing to the result of a subtraction.
/// Lists are sorted and deduplicated. Newly created intersection edges have
/// an empty source list unless they coincide with an operand boundary edge.
#[derive(Clone, Debug, Default)]
pub struct SubtractionHistory {
    /// Result face to the operand faces from which it was built.
    pub faces: BTreeMap<FaceId, Vec<FaceId>>,
    /// Result edge to the operand edges from which it was built.
    pub edges: BTreeMap<EdgeId, Vec<EdgeId>>,
}

/// A refused boolean, with the stage that broke.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BooleanRefusal {
    /// The first stage that broke.
    pub stage: Stage,
    /// Why.
    pub refusal: Refusal,
}

/// Where a point is relative to a body; the boundary is `ProvenCoincident`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Containment {
    /// Strictly inside the material.
    Inside,
    /// Strictly outside the material.
    Outside,
}

/// The refusal of a boolean on a body that holds a tolerant entity.
pub const TOLERANT_PARTICIPANT: &str = "boolean.tolerant-participant";

/// Adds the regularised result body, using set rules for repeated operands.
/// Both operands are exact bodies: one that stores a tolerance on a vertex,
/// an edge or a face is refused with [`TOLERANT_PARTICIPANT`] before any
/// decision is made.
pub fn boolean(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    dec: &mut Decider,
) -> Result<BooleanAnswer, BooleanRefusal> {
    boolean_impl(model, a, b, op, dec, false, None)
}

/// [`boolean`] for constructing a part, which merges sub-tolerance gaps
/// between the operands once. A plane of `b` parallel to a nearby
/// straight-edged face plane of `a`, with its face corners closer than the
/// 1 nm modelling tolerance to it, is merged with it first, on a copy of `b`; the result then carries a [`SUB_TOLERANCE_MERGE`]
/// mark bounding how far that copy's boundary moved. Without such a plane the
/// answer is the one [`boolean`] gives.
pub fn boolean_merging_gaps(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    dec: &mut Decider,
) -> Result<BooleanAnswer, BooleanRefusal> {
    boolean_impl(model, a, b, op, dec, true, None)
}

/// Subtract `tool` from `target`, retaining the source IDs of sewn topology.
pub fn subtract_with_history(
    model: &mut Model,
    target: BodyId,
    tool: BodyId,
    dec: &mut Decider,
) -> Result<(BooleanAnswer, SubtractionHistory), BooleanRefusal> {
    let mut history = SubtractionHistory::default();
    let answer = boolean_impl(
        model,
        target,
        tool,
        BoolOp::Subtract,
        dec,
        false,
        Some(&mut history),
    )?;
    Ok((answer, history))
}

fn boolean_impl(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    dec: &mut Decider,
    merge: bool,
    history: Option<&mut SubtractionHistory>,
) -> Result<BooleanAnswer, BooleanRefusal> {
    if [a, b]
        .iter()
        .any(|body| !matches!(model.body(*body).kind, BodyKind::Solid | BodyKind::Empty))
    {
        return Err(BooleanRefusal {
            stage: Stage::Intersect,
            refusal: geometry::refusal(
                "boolean.not-solid",
                "boolean operands must be solid or empty",
            ),
        });
    }
    // Checkpoint 2 decisions, section 4.5: "Tolerant entities never reach
    // the exact intersector." Everything below decides by exact predicates
    // on exact entities, and this is the one entry to it. A body with a
    // tolerance on any entity is cut by a construction that decides distance
    // against the tolerance instead, or it is refused here by name.
    if let Some((body, entity)) = [a, b]
        .into_iter()
        .find_map(|body| Some((body, model.tolerant_entity(body)?)))
    {
        return Err(BooleanRefusal {
            stage: Stage::Intersect,
            refusal: geometry::refusal(
                TOLERANT_PARTICIPANT,
                format!(
                    "{entity} of body {body:?} stores a tolerance; the exact boolean takes exact entities only"
                ),
            ),
        });
    }
    let before = dec.stats().clone();
    let merged = if merge
        && a != b
        && [a, b]
            .iter()
            .all(|body| model.body(*body).kind == BodyKind::Solid)
    {
        subgap::merge(model, a, b, dec).map_err(|refusal| BooleanRefusal {
            stage: Stage::Intersect,
            refusal,
        })?
    } else {
        None
    };
    // The search belongs to finding where the boundaries meet.
    let searched = statistics_delta(&before, dec.stats());
    let with_search = |mut answer: BooleanAnswer| {
        if let Some(first) = answer.report.stages.first_mut() {
            first.stats = statistics_sum(&searched, &first.stats);
        }
        answer
    };
    // A history would name the merged copy's faces, so no caller asks for both.
    debug_assert!(!(merge && history.is_some()));
    let b = merged.unwrap_or(b);
    combine(model, a, b, op, dec, history).map(with_search)
}

fn combine(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    dec: &mut Decider,
    mut history: Option<&mut SubtractionHistory>,
) -> Result<BooleanAnswer, BooleanRefusal> {
    // Appending is all a boolean does, so this snapshot shares the history.
    let snapshot = model.clone();
    let mut apart_proof = None;
    let refused = match staged(
        model,
        a,
        b,
        op,
        dec,
        &mut apart_proof,
        history.as_deref_mut(),
    ) {
        Ok(answer) => return Ok(answer),
        Err(refused) => refused,
    };
    // A carried certificate was paid by the main imprint call. Otherwise, try
    // the existing fitted-body fallback. Strictly apart lump boxes exclude
    // every possible surface or boundary meeting.
    if a == b
        || (!apart_proof.is_some_and(|proof| proof.matches(a, b))
            && !trim_bounds::TrimBounds::for_bodies(&snapshot, a, b)
                .bodies_apart(&snapshot, a, b, dec))
    {
        return Err(refused);
    }
    *model = snapshot;
    if let Some(history) = history {
        // A disjoint subtraction keeps precisely the target's topology.
        *history = SubtractionHistory {
            faces: model
                .faces_of_body(a)
                .into_iter()
                .map(|id| (id, vec![id]))
                .collect(),
            edges: model
                .edges_of_body(a)
                .into_iter()
                .map(|id| (id, vec![id]))
                .collect(),
        };
    }
    Ok(apart(model, a, b, op))
}

/// The answer for two solids proved apart: the target, both, or nothing.
fn apart(model: &mut Model, a: BodyId, b: BodyId, op: BoolOp) -> BooleanAnswer {
    let mut marks = model.body(a).marks.clone();
    marks.extend(&model.body(b).marks);
    let (kind, lumps) = match op {
        BoolOp::Subtract => (BodyKind::Solid, model.body(a).lumps.clone()),
        BoolOp::Unite => {
            let mut lumps = model.body(a).lumps.clone();
            lumps.extend_from_slice(&model.body(b).lumps);
            (BodyKind::Solid, lumps)
        }
        BoolOp::Intersect => (BodyKind::Empty, Vec::new()),
    };
    let body = model.add_body(Body {
        kind,
        lumps,
        wires: Vec::new(),
        marks: marks.clone(),
    });
    let (body, marks) = fit_marked(model, body, marks);
    // The difference is the target's own lumps, so a fitted blend keeps the
    // unfilleted base its certificate is written against.
    if op == BoolOp::Subtract
        && let Some(base) = model.blend_base(a)
    {
        model.record_blend_base(body, base);
    }
    BooleanAnswer {
        body,
        marks,
        report: StageReport::default(),
    }
}

fn staged(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    dec: &mut Decider,
    apart_proof: &mut Option<trim_bounds::BodyApartProof>,
    history: Option<&mut SubtractionHistory>,
) -> Result<BooleanAnswer, BooleanRefusal> {
    let mut operand_marks = model.body(a).marks.clone();
    operand_marks.extend(&model.body(b).marks);
    let a_empty = model.body(a).kind == BodyKind::Empty;
    let b_empty = model.body(b).kind == BodyKind::Empty;
    let empty = match op {
        BoolOp::Unite => a_empty && b_empty,
        BoolOp::Intersect => a_empty || b_empty,
        BoolOp::Subtract => a_empty || a == b,
    };
    if empty {
        if [a, b]
            .iter()
            .any(|body| !matches!(model.body(*body).kind, BodyKind::Solid | BodyKind::Empty))
        {
            return Err(BooleanRefusal {
                stage: Stage::Intersect,
                refusal: geometry::refusal(
                    "boolean.not-solid",
                    "boolean operands must be solid or empty",
                ),
            });
        }
        let body = model.add_body(Body {
            kind: BodyKind::Empty,
            lumps: Vec::new(),
            wires: Vec::new(),
            marks: operand_marks.clone(),
        });
        let (body, marks) = fit_marked(model, body, operand_marks);
        return Ok(BooleanAnswer {
            body,
            marks,
            report: StageReport::default(),
        });
    }
    let mut report = StageReport::default();
    let mut marks = MarkLog::new();
    macro_rules! run {
        ($stage:expr, $operation:expr) => {{
            #[cfg(feature = "work-span")]
            let _stage = catallax_numbers::work_span::region(format!("{:?}", $stage));
            let before = dec.stats().clone();
            let answer = $operation.map_err(|refusal| BooleanRefusal {
                stage: $stage,
                refusal,
            })?;
            report.stages.push(StageRecord {
                stage: $stage,
                stats: statistics_delta(&before, dec.stats()),
                marks: answer.marks.len(),
            });
            marks.extend(&answer.marks);
            answer.value
        }};
    }
    let body = if a == b {
        if model.body(a).kind != BodyKind::Solid {
            return Err(BooleanRefusal {
                stage: Stage::Intersect,
                refusal: geometry::refusal(
                    "boolean.not-solid",
                    "boolean operands must be solid or empty",
                ),
            });
        }
        run!(Stage::Sew, make_maximal(model, a, dec))
    } else {
        let graph = run!(
            Stage::Intersect,
            if op == BoolOp::Unite && bounds::bodies_separated(model, a, b, dec) {
                Ok(Answer::exact(None))
            } else {
                intersect(model, a, b, dec).map(|answer| answer.map(Some))
            }
        );
        if let Some(graph) = graph {
            let mut sources = history.is_some().then(imprinting::ImprintSources::default);
            let imprinted = run!(
                Stage::Imprint,
                stages::imprint_for_boolean(model, &graph, dec, apart_proof, sources.as_mut())
            );
            let classes = run!(Stage::Classify, classify(model, &imprinted, dec));
            let before = dec.stats().clone();
            let selection = keep(&classes, &imprinted, op);
            report.stages.push(StageRecord {
                stage: Stage::Keep,
                stats: statistics_delta(&before, dec.stats()),
                marks: 0,
            });
            run!(
                Stage::Sew,
                sewing::sew(
                    model,
                    &imprinted,
                    &selection,
                    dec,
                    sources.as_ref().zip(history)
                )
            )
        } else {
            let mut lumps = model.body(a).lumps.clone();
            lumps.extend_from_slice(&model.body(b).lumps);
            let body = model.add_body(Body {
                kind: BodyKind::Solid,
                lumps,
                wires: Vec::new(),
                marks: MarkLog::new(),
            });
            run!(Stage::Sew, make_maximal(model, body, dec))
        }
    };
    marks.extend(&operand_marks);
    let (body, marks) = fit_marked(model, body, marks);
    Ok(BooleanAnswer {
        body,
        marks,
        report,
    })
}

// A result keeps the operands' marks, except that its fit mark is the one its
// own fitted faces call for: the worst measured bound of the fitted faces it
// keeps, or none. An operand's fit bound for faces the result does not have
// is not carried.
fn fit_marked(model: &mut Model, body: BodyId, marks: MarkLog) -> (BodyId, MarkLog) {
    let marks = catallax_topology::with_fit_mark(model, body, &marks);
    if model.body(body).marks == marks {
        return (body, marks);
    }
    let mut result = model.body(body).clone();
    result.marks = marks.clone();
    (model.add_body(result), marks)
}

// Counters of `earlier` then `stage`; the running precision maximum is the stage's.
fn statistics_sum(earlier: &DecisionStats, stage: &DecisionStats) -> DecisionStats {
    DecisionStats {
        decisions: earlier.decisions + stage.decisions,
        filter: earlier.filter + stage.filter,
        refined: earlier.refined + stage.refined,
        exact: earlier.exact + stage.exact,
        coincident_exact: earlier.coincident_exact + stage.coincident_exact,
        coincident_separation: earlier.coincident_separation + stage.coincident_separation,
        coincident_fact: earlier.coincident_fact + stage.coincident_fact,
        not_exact: earlier.not_exact + stage.not_exact,
        refused: earlier.refused + stage.refused,
        max_bits: stage.max_bits,
        work: earlier.work + stage.work,
    }
}

fn statistics_delta(before: &DecisionStats, after: &DecisionStats) -> DecisionStats {
    DecisionStats {
        decisions: after.decisions - before.decisions,
        filter: after.filter - before.filter,
        refined: after.refined - before.refined,
        exact: after.exact - before.exact,
        coincident_exact: after.coincident_exact - before.coincident_exact,
        coincident_separation: after.coincident_separation - before.coincident_separation,
        coincident_fact: after.coincident_fact - before.coincident_fact,
        not_exact: after.not_exact - before.not_exact,
        refused: after.refused - before.refused,
        max_bits: after.max_bits,
        work: after.work - before.work,
    }
}

/// Whether an exact point is inside, outside or on the boundary of a body.
/// May add auxiliary geometry to the model's store; never changes topology.
pub fn classify_point<P: Clone + Into<Point3>>(
    model: &mut Model,
    body: BodyId,
    point: &P,
    dec: &mut Decider,
) -> Outcome<Containment> {
    containment::classify_point(model, body, &point.clone().into(), dec)
}

/// An exact point strictly inside a face, retaining its chart definition.
pub fn face_point(model: &mut Model, face: FaceId, dec: &mut Decider) -> Answered<Point3> {
    points::face_point(model, face, dec)
}

/// Membership of an exact point in a trimmed face. Rational callers are unchanged.
/// This query may append auxiliary geometry, but never changes topology.
/// Unsettled support or trimming predicates return a refusal, not membership.
/// A point proved off a finite spline carrier is no member, with no side query.
pub fn face_contains<P: Clone + Into<Point3>>(
    model: &mut Model,
    face: FaceId,
    point: &P,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    face_contains_with(model, face, point, None, dec, region::locate)
}

/// Membership of a stored exact point, preserving its support incidence facts.
/// The point must belong to `model`; the face trim is still checked.
pub fn face_contains_point(
    model: &mut Model,
    face: FaceId,
    point: catallax_shapes::PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let position = model.geom().position(point);
    face_contains_with(model, face, &position, Some(point), dec, region::locate)
}

fn face_contains_with<P: Clone + Into<Point3>>(
    model: &mut Model,
    face: FaceId,
    point: &P,
    stored: Option<catallax_shapes::PointId>,
    dec: &mut Decider,
    locate: impl FnOnce(
        &mut Model,
        FaceId,
        catallax_shapes::PointId,
        &mut Decider,
    ) -> Result<region::Location, Refusal>,
) -> Result<bool, Refusal> {
    let surface = model.face(face).surface;
    let point: Point3 = point.clone().into();
    if matches!(
        model.geom().surface(surface),
        catallax_shapes::Surface::BSpline(catallax_shapes::SplineSupport::RationalPatch(_))
    ) && !spline_fillet::borders_fillet(model, face)
    {
        let checked =
            catallax_topology::RationalPatchFace::of_face(model, face, dec)?.ok_or_else(|| {
                Refusal::new(
                    "boolean.rational-patch-domain",
                    "finite rational face membership needs a checked complete trim domain",
                )
            })?;
        return Ok(checked.locate(&point, dec)? != catallax_numbers::Sign::Positive);
    }
    if matches!(
        model.geom().surface(surface),
        catallax_shapes::Surface::BSpline(catallax_shapes::SplineSupport::BezierExtrusion(_))
    ) {
        let band = catallax_topology::BezierExtrusionBand::of_face(model, face, dec)?.ok_or_else(
            || {
                Refusal::new(
                    "boolean.bezier-band-required",
                    "finite Bezier face membership needs a checked trim domain",
                )
            },
        )?;
        return Ok(band.locate(&point, dec)? != catallax_numbers::Sign::Positive);
    }
    if let catallax_shapes::Surface::Cone(cone) = model.geom().surface(surface)
        && face_sign(&point.minus(&cone.apex).dot(&cone.axis), dec)?
            == catallax_numbers::Sign::Negative
    {
        return Ok(false);
    }
    if let catallax_shapes::Surface::BSpline(support) = model.geom().surface(surface)
        && support.misses_point(&point, dec)?
    {
        return Ok(false);
    }
    let p = stored.unwrap_or_else(|| {
        model
            .geom_mut()
            .add_point(catallax_shapes::PointDef::Explicit(point.clone()))
    });
    if matches!(
        model.geom().surface(surface),
        catallax_shapes::Surface::BSpline(catallax_shapes::SplineSupport::RationalAtlas(_))
    ) {
        if !geometry::on_surface(model, surface, p, dec, "boolean.face-unresolved")? {
            return Ok(false);
        }
    } else {
        match model.geom().side(surface, p, dec) {
            Outcome::ProvenCoincident(proof)
                if proof.proof != catallax_numbers::ZeroProof::SeparationBound => {}
            Outcome::Proven(_) => return Ok(false),
            Outcome::Refused(error) => return Err(error),
            other => {
                return Err(Refusal::new(
                    "boolean.face-unresolved",
                    "face support incidence is not exact",
                )
                .keeping_undecided(other));
            }
        }
    }
    if let Some(inside) = comparison_cone_band(model, face, &point, dec)? {
        return Ok(inside);
    }
    if let Some(inside) = region::prism_face(model, face, p, dec)? {
        return Ok(inside);
    }
    // The incidence proves the point is on this support. The ray that
    // found it is not a useful seed for constructing a trim probe.
    let p = if stored.is_some() {
        model
            .geom_mut()
            .add_point(catallax_shapes::PointDef::Explicit(point))
    } else {
        p
    };
    Ok(!matches!(
        locate(model, face, p, dec)?,
        region::Location::Outside
    ))
}

/// Constructs an exact rational supporting-surface witness near a target.
/// This is a rational chart projection, not necessarily the closest point.
/// Callers must measure its distance and check trimmed-face membership.
/// `None` means the surface frame has no supported rational parametrisation.
/// `width` must be positive. It controls coordinate/parameter refinement, not
/// the final physical distance to the target.
pub fn surface_witness(
    surface: &catallax_shapes::Surface,
    target: &catallax_shapes::Point3,
    width: &catallax_numbers::Rational,
    dec: &mut Decider,
) -> Result<Option<RPoint3>, Refusal> {
    points::comparison_witness(surface, target, width, dec)
}
fn face_sign(
    value: &catallax_numbers::Real,
    dec: &mut Decider,
) -> Result<catallax_numbers::Sign, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.sign());
    }
    geometry::sign(value, dec, "boolean.face-unresolved")
}

// A complete circular band has exact height limits. This query-only path avoids
// asking the intersection engine for a conical chart it does not implement.
fn comparison_cone_band(
    model: &Model,
    face: FaceId,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    use catallax_numbers::{Real, Sign};
    use catallax_shapes::{Curve, PointDef, Surface};
    use catallax_topology::{Loop, Sense};
    let f = model.face(face);
    let Surface::Cone(s) = model.geom().surface(f.surface) else {
        return Ok(None);
    };
    if f.loops.len() != 2 {
        return Ok(None);
    }
    let mut bounds: [(Real, Option<Real>); 2] = std::array::from_fn(|_| (Real::zero(), None));
    for (index, &lp) in f.loops.iter().enumerate() {
        match model.loop_(lp) {
            Loop::Vertex(v) => {
                let PointDef::Explicit(p) = model.geom().point(model.vertex(*v).point) else {
                    return Ok(None);
                };
                if face_sign(&p.minus(&s.apex).norm_squared(), dec)? != Sign::Zero {
                    return Ok(None);
                }
            }
            Loop::Coedges(cs) if cs.len() == 1 => {
                let coedge = model.coedge(cs[0]);
                let edge = model.edge(coedge.edge);
                if edge.start.is_some() || edge.end.is_some() {
                    return Ok(None);
                }
                let Curve::Circle(c) = model.geom().curve(edge.curve) else {
                    return Ok(None);
                };
                let delta = c.center.minus(&s.apex);
                if face_sign(&delta.cross(&s.axis).norm_squared(), dec)? != Sign::Zero
                    || face_sign(&c.normal.cross(&s.axis).norm_squared(), dec)? != Sign::Zero
                {
                    return Ok(None);
                }
                let h = delta.dot(&s.axis);
                if face_sign(&h, dec)? == Sign::Negative {
                    return Ok(None);
                }
                let radius = &c.radius_sq * &s.axis.norm_squared();
                let expected = (&s.slope * &s.slope) * (&h * &h);
                if face_sign(&(radius - expected), dec)? != Sign::Zero {
                    return Ok(None);
                }
                let mut orientation = c.normal.dot(&s.axis);
                for sense in [coedge.sense, edge.sense, f.sense] {
                    if sense == Sense::Reversed {
                        orientation = -orientation;
                    }
                }
                bounds[index] = (h, Some(orientation));
            }
            _ => return Ok(None),
        }
    }
    if face_sign(&(&bounds[0].0 - &bounds[1].0), dec)? == Sign::Positive {
        bounds.swap(0, 1);
    }
    if let Some(lower) = &bounds[0].1
        && face_sign(lower, dec)? != Sign::Positive
    {
        return Ok(None);
    }
    let Some(upper) = &bounds[1].1 else {
        return Ok(None);
    };
    if face_sign(upper, dec)? != Sign::Negative {
        return Ok(None);
    }
    let h = point.minus(&s.apex).dot(&s.axis);
    Ok(Some(
        face_sign(&(&h - &bounds[0].0), dec)? != Sign::Negative
            && face_sign(&(&h - &bounds[1].0), dec)? != Sign::Positive,
    ))
}
