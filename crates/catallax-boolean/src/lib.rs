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
    imprint, intersect, keep, make_maximal, sew,
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

/// Membership of an exact rational point in a trimmed face.
/// This query may append auxiliary geometry, but never changes topology.
/// Unsettled support or trimming predicates return a refusal, not membership.
pub fn face_contains(
    model: &mut Model,
    face: FaceId,
    point: &RPoint3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let surface = model.face(face).surface;
    if let catallax_shapes::Surface::Cone(cone) = model.geom().surface(surface)
        && point.minus(&cone.apex).dot(&cone.axis) < catallax_numbers::Rational::zero()
    {
        return Ok(false);
    }
    let p = model
        .geom_mut()
        .add_point(catallax_shapes::PointDef::Explicit(point.clone()));
    match model.geom().side(surface, p, dec) {
        Outcome::ProvenCoincident(proof)
            if proof.proof != catallax_numbers::ZeroProof::SeparationBound => {}
        Outcome::Proven(_) => return Ok(false),
        Outcome::Refused(error) => return Err(error),
        _ => {
            return Err(Refusal::new(
                "boolean.face-unresolved",
                "face support incidence is not exact",
            ));
        }
    }
    if let Some(inside) = comparison_cone_band(model, face, point) {
        return Ok(inside);
    }
    Ok(!matches!(
        region::locate(model, face, p, dec)?,
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

// A complete circular band has exact height limits. This query-only path avoids
// asking the intersection engine for a conical chart it does not implement.
fn comparison_cone_band(model: &Model, face: FaceId, point: &RPoint3) -> Option<bool> {
    use catallax_shapes::{Curve, PointDef, Surface};
    use catallax_topology::{Loop, Sense};
    let f = model.face(face);
    let Surface::Cone(s) = model.geom().surface(f.surface) else {
        return None;
    };
    if f.loops.len() != 2 {
        return None;
    }
    let mut bounds = std::array::from_fn::<_, 2, _>(|_| (catallax_numbers::Rational::zero(), None));
    for (index, &lp) in f.loops.iter().enumerate() {
        match model.loop_(lp) {
            Loop::Vertex(v) => {
                let PointDef::Explicit(p) = model.geom().point(model.vertex(*v).point) else {
                    return None;
                };
                if p != &s.apex {
                    return None;
                }
                bounds[index] = (catallax_numbers::Rational::zero(), None);
            }
            Loop::Coedges(cs) if cs.len() == 1 => {
                let coedge = model.coedge(cs[0]);
                let edge = model.edge(coedge.edge);
                if edge.start.is_some() || edge.end.is_some() {
                    return None;
                }
                let Curve::Circle(c) = model.geom().curve(edge.curve) else {
                    return None;
                };
                let delta = c.center.minus(&s.apex);
                if !delta.cross(&s.axis).is_zero() || !c.normal.cross(&s.axis).is_zero() {
                    return None;
                }
                let h = delta.dot(&s.axis);
                if h < catallax_numbers::Rational::zero() {
                    return None;
                }
                if &c.radius_sq * &s.axis.norm_squared() != &(&s.slope * &s.slope) * &(&h * &h) {
                    return None;
                }
                let mut orientation = c.normal.dot(&s.axis);
                for sense in [coedge.sense, edge.sense, f.sense] {
                    if sense == Sense::Reversed {
                        orientation = -orientation;
                    }
                }
                bounds[index] = (h, Some(orientation));
            }
            _ => return None,
        }
    }
    bounds.sort_by(|a, b| a.0.cmp(&b.0));
    if bounds[0]
        .1
        .as_ref()
        .is_some_and(|x| x <= &catallax_numbers::Rational::zero())
        || bounds[1]
            .1
            .as_ref()
            .is_none_or(|x| x >= &catallax_numbers::Rational::zero())
    {
        return None;
    }
    let h = point.minus(&s.apex).dot(&s.axis);
    Some(bounds[0].0 <= h && h <= bounds[1].0)
}
