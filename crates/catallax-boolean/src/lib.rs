//! Piece 5 of catallax: the boolean engine.
//!
//! Unite, subtract and intersect through five reusable stages: intersect,
//! imprint, classify, keep, sew. Plus exact containment of rational points
//! and rational points inside faces.
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

mod containment;
mod curve_pair;
mod fast_region;
mod geometry;
mod germ;
mod graph;
mod imprinting;
mod points;
mod probe_path;
mod probed_region;
mod project;
mod region;
mod sewing;
mod stages;
mod torus_region;

#[cfg(test)]
mod tests;

pub use stages::{
    Classification, FaceClass, Imprinted, IntersectionGraph, SectionEdge, Selection, classify,
    imprint, intersect, keep, sew,
};

use catallax_numbers::{Answered, Decider, DecisionStats, MarkLog, Outcome, Refusal};
use catallax_shapes::RPoint3;
use catallax_topology::{Body, BodyId, BodyKind, FaceId, Model};

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

/// Runs all five stages and adds the result body to the model.
pub fn boolean(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    dec: &mut Decider,
) -> Result<BooleanAnswer, BooleanRefusal> {
    let mut operand_marks = model.body(a).marks.clone();
    operand_marks.extend(&model.body(b).marks);
    let a_empty = model.body(a).kind == BodyKind::Empty;
    let b_empty = model.body(b).kind == BodyKind::Empty;
    let empty = match op {
        BoolOp::Unite => a_empty && b_empty,
        BoolOp::Intersect => a_empty || b_empty,
        BoolOp::Subtract => a_empty,
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
            marks: operand_marks.clone(),
        });
        return Ok(BooleanAnswer {
            body,
            marks: operand_marks,
            report: StageReport::default(),
        });
    }
    let mut report = StageReport::default();
    let mut marks = MarkLog::new();
    macro_rules! run {
        ($stage:expr, $operation:expr) => {{
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
    let graph = run!(Stage::Intersect, intersect(model, a, b, dec));
    let imprinted = run!(Stage::Imprint, imprint(model, &graph, dec));
    let classes = run!(Stage::Classify, classify(model, &imprinted, dec));
    let before = dec.stats().clone();
    let selection = keep(&classes, &imprinted, op);
    report.stages.push(StageRecord {
        stage: Stage::Keep,
        stats: statistics_delta(&before, dec.stats()),
        marks: 0,
    });
    let body = run!(Stage::Sew, sew(model, &imprinted, &selection, dec));
    marks.extend(&operand_marks);
    let body = if model.body(body).marks == marks {
        body
    } else {
        let mut result = model.body(body).clone();
        result.marks = marks.clone();
        model.add_body(result)
    };
    Ok(BooleanAnswer {
        body,
        marks,
        report,
    })
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

/// Whether a rational point is inside, outside or on the boundary of a body.
/// May add auxiliary geometry to the model's store; never changes topology.
pub fn classify_point(
    model: &mut Model,
    body: BodyId,
    point: &RPoint3,
    dec: &mut Decider,
) -> Outcome<Containment> {
    containment::classify_point(model, body, point, dec)
}

/// A rational point strictly inside a face.
pub fn face_point(model: &mut Model, face: FaceId, dec: &mut Decider) -> Answered<RPoint3> {
    points::face_point(model, face, dec)
}
