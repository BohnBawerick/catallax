//! Rendering answers with a preset bounded spline rule.

use super::*;
use catallax_numbers::{
    Answer, BoundedAnswer, BoundedEvidence, BoundedRule, CertificateIdentity, CertificateScope,
    EvidenceBound, EvidenceUnit, OperationResult, OperationTask, operation_budget_active,
    with_operation_budget,
};
use catallax_shapes::{SplineSupport, Surface};

/// Fixed policy for rendering spline boundaries. It grants no solid or export proof.
pub const SPLINE_DISPLAY_RULE: BoundedRule = BoundedRule {
    case_kind: "spline display mesh",
    name: "spline-display-evaluator-v1",
    level: "1/100",
    unit: EvidenceUnit::Millimetres,
};

/// How the complete rendering boundary was compared with the retained surfaces.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DisplayMeasurement {
    /// Complete-cell method, including shared-edge replacement and binary64 rounding.
    pub method: &'static str,
    /// Largest complete face or physical-edge deviation, in millimetres.
    pub maximum_deviation: Rational,
    /// Number of original finite faces covered by the measurement.
    pub covered_faces: usize,
}

/// Rendering data with picking handles and a distance measurement.
/// There is no conversion to a certified [`Mesh`] or an export writer.
///
/// ```compile_fail
/// use catallax_mesh::{DisplayMesh, Mesh};
/// fn as_export_mesh(display: DisplayMesh) -> Mesh { display.into() }
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayMesh {
    /// Positions, triangles and shared physical edges for rendering and selection.
    pub geometry: PreviewMesh,
    /// Whole-boundary measurement against the kept model.
    pub measurement: DisplayMeasurement,
    /// Marks retained from the body and sampled exact boundary traces.
    pub source_marks: catallax_numbers::MarkLog,
}

fn rendering(
    mesh: Mesh,
    source_marks: catallax_numbers::MarkLog,
    method: &'static str,
) -> DisplayMesh {
    let maximum_deviation = mesh
        .faces
        .iter()
        .map(|f| &f.max_chord_error)
        .chain(mesh.edges.iter().map(|e| &e.max_chord_error))
        .max()
        .cloned()
        .unwrap_or_else(Rational::zero);
    let covered_faces = mesh.faces.len();
    DisplayMesh {
        geometry: PreviewMesh {
            positions: mesh.positions,
            triangles: mesh.triangles,
            faces: mesh
                .faces
                .into_iter()
                .map(|f| PreviewFace {
                    face: f.face,
                    triangles: f.triangles,
                })
                .collect(),
            edges: mesh
                .edges
                .into_iter()
                .map(|e| PreviewEdge {
                    edge: e.edge,
                    vertices: e.vertices,
                })
                .collect(),
        },
        source_marks,
        measurement: DisplayMeasurement {
            method,
            maximum_deviation,
            covered_faces,
        },
    }
}

/// Tries certification first, then the fixed 0.01 mm spline rendering rule.
/// The certified attempt receives one quarter of the remaining caller allowance;
/// its actual work is charged before the positional attempt spends what remains.
/// No retry creates a fresh allowance. Other tolerances and optional controls use
/// the existing certified entry. `revision` must identify the immutable model.
pub fn display_mesh_body(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    revision: &str,
    dec: &mut Decider,
) -> OperationResult<DisplayMesh> {
    if !operation_budget_active() {
        let (result, work) = with_operation_budget(20_000_000, || {
            display_mesh_body(model, body, options, revision, dec)
        });
        return if work.exhausted && !matches!(result, OperationResult::Refused(_)) {
            OperationResult::Refused(Refusal::new(
                "mesh.spline-display-budget",
                "spline display exhausted the shared 20000000-unit allowance",
            ))
        } else {
            result
        };
    }
    let faces = model.faces_of_body(body);
    let eligible = options.chord_tolerance == Rational::from_ratio(1, 100).unwrap()
        && options.angular_deviation.is_none()
        && options.max_edge_length.is_none()
        && options.max_aspect_ratio.is_none()
        && faces.iter().any(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::BSpline(
                    SplineSupport::SplineFillet(_)
                        | SplineSupport::RationalPatch(_)
                        | SplineSupport::RationalAtlas(_)
                )
            )
        })
        && faces.iter().all(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Plane(_) | Surface::BSpline(_)
            )
        });
    let strict = if eligible {
        let task = OperationTask::capture()
            .split(std::num::NonZeroUsize::new(4).unwrap())
            .remove(0);
        match task
            .run(|| mesh_body(model, body, options, dec))
            .commit_share()
        {
            Ok(Some(result)) => result,
            Ok(None) => Err(Refusal::new(
                "numbers.operation-budget",
                "certified spline display exhausted its reserved attempt",
            )),
            Err(error) => return OperationResult::Refused(error),
        }
    } else {
        mesh_body(model, body, options, dec)
    };
    if let Ok(answer) = strict {
        return OperationResult::Proven(Answer {
            value: rendering(
                answer.value,
                answer.marks.clone(),
                "certified whole-cell chord and embedding bounds",
            ),
            marks: answer.marks,
        });
    }
    if !eligible {
        return OperationResult::Refused(strict.unwrap_err());
    }
    let result = (|| {
        let (mesh, marks) = crate::triangulate::build_position(model, body, options, dec)?;
        let display = rendering(
            mesh,
            marks,
            "3x3 direct surface evaluation, or endpoint-curve stencils for checked parallel rulings; complete source-span derivative residual enclosures between samples; cut-domain displacement, shared-edge replacement and binary64 rounding",
        );
        if display.measurement.covered_faces != faces.len() {
            return Err(failure("positional display omits an original face"));
        }
        dec.check_deadline()?;
        let evidence = BoundedEvidence::new(
            &SPLINE_DISPLAY_RULE,
            CertificateIdentity {
                revision: revision.to_owned(),
                representation: format!("body:{}", body.index()),
                target: "retained exact finite B-rep boundary".to_owned(),
                scope: CertificateScope::Entities(
                    faces
                        .iter()
                        .map(|f| format!("face:{}", f.index()))
                        .collect(),
                ),
            },
            EvidenceBound {
                upper: display.measurement.maximum_deviation.clone(),
                unit: EvidenceUnit::Millimetres,
            },
        )?;
        Ok(BoundedAnswer::new(display, evidence))
    })();
    match result {
        Ok(answer) => OperationResult::Bounded(answer),
        Err(error)
            if catallax_numbers::current_operation_work().exhausted
                || error.code == "numbers.operation-budget"
                || error.code == "numbers.budget-exhausted" =>
        {
            OperationResult::Refused(Refusal::new(
                "mesh.spline-display-budget",
                format!(
                    "spline display mesh at level 1/100 mm: direct evaluation exhausted its share after {} total request units; {}: {}",
                    catallax_numbers::current_operation_work().used,
                    error.code,
                    error.detail
                ),
            ))
        }
        Err(mut error) => {
            error.detail = format!(
                "spline display mesh level 1/100 mm, vertex cap {}, {} total request units: {}",
                options.max_vertices,
                catallax_numbers::current_operation_work().used,
                error.detail
            );
            OperationResult::Refused(error)
        }
    }
}
