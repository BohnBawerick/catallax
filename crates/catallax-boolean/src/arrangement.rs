//! One shared surface arrangement, classified before any boundary is sewn.

use std::collections::BTreeMap;

use catallax_numbers::{Answer, Decider, MarkLog, Refusal, current_operation_work};
use catallax_shapes::Point3;
use catallax_topology::{BodyId, BodyKind, FaceId, Model, Sense};

use crate::{BooleanAnswer, BooleanRefusal, Selection, Stage, StageRecord, StageReport};

/// A proved final material transition on one arrangement cell.
#[derive(Clone, Copy, Debug)]
pub enum PatchSelection {
    /// No boundary of the final material remains on this cell.
    Discard,
    /// Keep the cell with this sense relative to its input patch.
    Keep(Sense),
}

/// Arranges the source boundary and finite candidate patches together.
///
/// `classify` receives the original patch and a strict interior point of each
/// imprinted cell. It must prove the final material transition there. `None`
/// asks for another certified interior witness; `Discard` removes the cell.
/// Candidate/source bounds only exclude impossible meetings.
/// The ordinary sewing checks validate the final selected boundary once.
pub fn arrange(
    model: &mut Model,
    source: BodyId,
    candidates: &[FaceId],
    mut classify: impl FnMut(
        &mut Model,
        FaceId,
        &Point3,
        &mut Decider,
    ) -> Result<Option<PatchSelection>, Refusal>,
    dec: &mut Decider,
) -> Result<BooleanAnswer, BooleanRefusal> {
    if !matches!(model.body(source).kind, BodyKind::Solid | BodyKind::Empty) {
        return Err(BooleanRefusal {
            stage: Stage::Intersect,
            refusal: crate::geometry::refusal("boolean.not-solid", "arrangement source is not solid or empty"),
        });
    }
    if model.body(source).kind == BodyKind::Empty {
        return Ok(BooleanAnswer {
            body: source,
            marks: model.body(source).marks.clone(),
            report: StageReport::default(),
        });
    }
    let mut faces = model.faces_of_body(source);
    faces.extend_from_slice(candidates);
    let mut report = StageReport::default();
    let mut marks = MarkLog::new();
    macro_rules! run {
        ($stage:expr, $operation:expr) => {{
            let before = dec.stats().clone();
            let work_before = current_operation_work().used;
            let result = $operation;
            eprintln!("MULTIWAY phase={:?} work={} cumulative={} outcome={:?}", $stage, current_operation_work().used - work_before, current_operation_work().used, result.as_ref().map(|_| "ok").map_err(|e| e.code));
            let answer = result.map_err(|refusal| BooleanRefusal { stage: $stage, refusal })?;
            report.stages.push(StageRecord {
                stage: $stage,
                stats: crate::statistics_delta(&before, dec.stats()),
                marks: answer.marks.len(),
            });
            marks.extend(&answer.marks);
            answer.value
        }};
    }
    let graph = run!(Stage::Intersect, crate::graph::arrangement(model, source, &faces, dec));
    let (imprinted, owners) = run!(Stage::Imprint, crate::imprinting::arrangement(model, &graph, &faces, dec));
    let selected = run!(Stage::Classify, {
        let mut selected = Vec::new();
        let result: Result<_, Refusal> = (|| {
            for &piece in &imprinted.first {
                let owner = owners[&piece];
                let choice = crate::points::find(model, piece, dec, |model, point, _, dec| {
                    let point = model.geom().position(point);
                    classify(model, owner, &point, dec)
                })?;
                if let PatchSelection::Keep(sense) = choice {
                    selected.push((piece, sense));
                }
            }
            Ok(Answer::exact(selected))
        })();
        result
    });
    let selection = Selection {
        faces: selected,
        coincident: BTreeMap::new(),
        intersection: false,
        marks: marks.clone(),
    };
    let body = run!(Stage::Sew, crate::sewing::sew(model, &imprinted, &selection, dec));
    marks.extend(&model.body(source).marks);
    Ok(BooleanAnswer { body, marks, report })
}

/// Splits a curved patch against a bounded clipping solid without sewing an
/// intermediate body. The caller certifies each retained cell with `keep`;
/// `None` requests another strict interior witness. All marks are returned.
pub fn trim_patch(
    model: &mut Model,
    source: BodyId,
    patch: FaceId,
    mut keep: impl FnMut(&Point3, &mut Decider) -> Result<Option<bool>, Refusal>,
    dec: &mut Decider,
) -> Result<Answer<Vec<FaceId>>, BooleanRefusal> {
    if !matches!(model.body(source).kind, BodyKind::Solid | BodyKind::Empty) {
        return Err(BooleanRefusal {
            stage: Stage::Intersect,
            refusal: crate::geometry::refusal("boolean.not-solid", "patch clipping source is not solid or empty"),
        });
    }
    let mut marks = model.body(source).marks.clone();
    if model.body(source).kind == BodyKind::Empty {
        return Ok(Answer::with_marks(Vec::new(), marks));
    }
    let mut faces = model.faces_of_body(source);
    faces.push(patch);
    macro_rules! run {
        ($stage:expr, $operation:expr) => {{
            let before = current_operation_work().used;
            let result = $operation;
            eprintln!("NEAREST_TRIM phase={:?} work={} cumulative={} outcome={:?}", $stage, current_operation_work().used-before, current_operation_work().used, result.as_ref().map(|_|"ok").map_err(|e|e.code));
            let answer = result.map_err(|refusal| BooleanRefusal { stage: $stage, refusal })?;
            marks.extend(&answer.marks);
            answer.value
        }};
    }
    let graph = run!(Stage::Intersect, crate::graph::arrangement(model, source, &faces, dec));
    let (imprinted, owners) = run!(Stage::Imprint, crate::imprinting::arrangement(model, &graph, &faces, dec));
    let mut retained = Vec::new();
    let before = current_operation_work().used;
    let classified: Result<(), Refusal> = (|| {
        for &piece in &imprinted.first {
            if owners[&piece] != patch {
                continue;
            }
            let selected = crate::points::find(model, piece, dec, |model, point, _, dec| {
                keep(&model.geom().position(point), dec)
            })?;
            if selected {
                retained.push(piece);
            }
        }
        Ok(())
    })();
    eprintln!("NEAREST_TRIM phase=Classify work={} cumulative={} outcome={:?}", current_operation_work().used-before, current_operation_work().used, classified.as_ref().map(|_|"ok").map_err(|e|e.code));
    classified.map_err(|refusal| BooleanRefusal { stage: Stage::Classify, refusal })?;
    Ok(Answer::with_marks(retained, marks))
}
