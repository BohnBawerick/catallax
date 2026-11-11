//! Piece 7 of catallax: STEP export.
//!
//! Writes bodies as STEP AP242 files and proves the largest distance between
//! the exported and the exact geometry.
//!
//! The contract is `docs/pieces/07-step-export.md`.

use std::io::Write;

use catallax_numbers::{Answered, Decider, Rational};
use catallax_topology::{BodyId, Model};

mod band;
mod certified_spline;
mod chart;
mod cut;
mod equator;
mod export;
mod fit;
mod frames;
mod graph;
mod joint;
mod junction;
mod number;
mod partial;
mod pcurve;
mod poles;
mod rational_atlas;
mod regularize;
mod seams;
mod trim;
mod uv;
mod uv_loop;
mod wire;

pub use certified_spline::{CertifiedSplineStepOptions, write_spline_prism_step};

/// Export options.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct StepOptions {
    /// Product name written to the file.
    pub product_name: String,
    /// Requested maximum deviation, in millimetres, including rounding and spline fitting.
    pub curve_tolerance: Rational,
    /// Significant digits of every written number.
    pub significant_digits: u32,
    /// Which exchange contract governs seam curves and pcurves.
    pub seams: Seams,
}

/// The seam-curve rule a file must meet.
///
/// Both settings write the same bytes. `Required` additionally checks the
/// finished file and refuses `step.seam-unsupported` instead of writing one
/// whose curved faces a reader would have to complete itself.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub enum Seams {
    /// The original contract: pcurves and seam curves where an encoder writes
    /// them, and otherwise edges with their 3D curve only.
    #[default]
    Legacy,
    /// Every edge of a face that is not planar carries a `PCURVE` on that face
    /// through a `SURFACE_CURVE` or `SEAM_CURVE`, and an edge the face uses
    /// twice is a `SEAM_CURVE` with a pcurve for each use. A `VERTEX_LOOP`
    /// appears only on a sphere.
    Required,
}

impl Default for StepOptions {
    fn default() -> Self {
        Self {
            product_name: "catallax".into(),
            curve_tolerance: Rational::from_ratio(1, 10_000_000).unwrap(),
            significant_digits: 17,
            seams: Seams::Legacy,
        }
    }
}

/// A repair made only in the exported topology.
///
/// The exact body's face count equals `StepReport::faces` plus every record's
/// `removed_faces`, and its explicit edge count, before the writer's seams and
/// rings, equals `StepReport::edges` plus every record's `removed_edges`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct StepRegularization {
    /// Model face holding the ribbon: the face whose two decimal-rounded
    /// boundary wires overlap, the narrow step face removed from the file, or
    /// the face whose wire closes over a narrow wedge or loses a redundant
    /// straight boundary subdivision. A repair that follows
    /// another can name a face or edge of the writer's copy instead.
    pub face: catallax_topology::FaceId,
    /// The two straight edges bounding the collapsed ribbon: outer and
    /// contained inner, or the moved wall's and the kept wall's. For a boundary
    /// subdivision, these are the short edge and its collinear continuation.
    pub edges: [catallax_topology::EdgeId; 2],
    /// Neighbouring faces joined across that ribbon, or the two faces whose
    /// boundary subdivision was removed without joining faces.
    pub merged_faces: [catallax_topology::FaceId; 2],
    /// Proven displacement bound in millimetres, zero for a redundant subdivision.
    pub max_deviation: Rational,
    /// Faces this repair took out of the file.
    pub removed_faces: usize,
    /// Net explicit edges this repair took out of the file.
    pub removed_edges: usize,
}

/// One written face and the model face it covers.
///
/// A file may split one model face into several charts. Each record names a
/// reachable `ADVANCED_FACE`; together they cover the model's faces once.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct StepFaceRegion {
    /// Model face this file face belongs to.
    pub face: catallax_topology::FaceId,
    /// The `ADVANCED_FACE` name in the file.
    pub file_region: String,
    /// Exact `[[u_lo, u_hi], [v_lo, v_hi]]` chart of a spline patch.
    pub file_domain: Option<[[Rational; 2]; 2]>,
}

/// What an export wrote and how far it is from the exact model.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct StepReport {
    /// Proven upper bound, in millimetres, on the distance between exported and exact geometry.
    pub max_deviation: Rational,
    /// Intersection branches written as fitted splines.
    pub approximated_curves: usize,
    /// Vertices added for ring origins, natural poles, or periodic chart cuts.
    pub added_vertices: usize,
    /// Net explicit edges added by seams and shared-edge splits for periodic charts.
    pub added_edges: usize,
    /// Solids written.
    pub solids: usize,
    /// Faces written.
    pub faces: usize,
    /// Explicit EDGE_CURVE entities written, excluding implicit natural boundaries.
    pub edges: usize,
    /// Explicit VERTEX_POINT entities written.
    pub vertices: usize,
    /// Export-only topology changes, each with its own certified displacement.
    pub regularizations: Vec<StepRegularization>,
    /// Every written face, for writers that map file faces to model regions.
    /// Other writers leave it empty.
    pub face_map: Vec<StepFaceRegion>,
}

/// Checks a finished STEP file against the seam contract of [`Seams::Required`],
/// refusing `step.seam-unsupported` with the count of faces that break it.
/// Writers that take no [`StepOptions`] can be held to the contract this way.
pub fn check_seams(file: &[u8]) -> Result<(), catallax_numbers::Refusal> {
    seams::check(&String::from_utf8_lossy(file))
}

/// Writes `bodies` as one STEP AP242 file to `out`.
///
/// With [`Seams::Required`] the file is checked before any byte leaves; a file
/// that breaks the seam contract is refused with `step.seam-unsupported`.
pub fn write_step<W: Write>(
    model: &Model,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    if options.seams == Seams::Legacy {
        return write_any(model, bodies, options, dec, out);
    }
    let mut file = Vec::new();
    let answer = write_any(model, bodies, options, dec, &mut file)?;
    check_seams(&file)?;
    out.write_all(&file)
        .map_err(|error| catallax_numbers::Refusal::new("step.io-error", error.to_string()))?;
    Ok(answer)
}

fn write_any<W: Write>(
    model: &Model,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    if bodies.iter().any(|body| {
        model
            .chart_expansion(*body)
            .is_some_and(catallax_topology::ChartExpansion::needs_cells)
    }) {
        return rational_atlas::write(model, bodies, options, dec, out);
    }
    if let [body] = bodies
        && catallax_topology::DraftedCylinder::of_body(model, *body, dec)?.is_some()
    {
        let spline_options = CertifiedSplineStepOptions {
            product_name: options.product_name.clone(),
            file_tolerance: options.curve_tolerance.clone(),
            decimal_places: options.significant_digits,
        };
        return certified_spline::oblique::write(model, *body, &spline_options, dec, out);
    }
    write_prepared(model, bodies, options, dec, out, false)
}

fn write_prepared<W: Write>(
    model: &Model,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
    face_regions: bool,
) -> Answered<StepReport> {
    let prepared = regularize::prepare(model, bodies, options, dec)?;
    let answer = if let Some(prepared) = prepared {
        export::Exporter::new(&prepared.model, options, dec)
            .face_regions(face_regions)
            .regularizations(prepared.records, prepared.face_errors, prepared.edge_errors)
            .build(&prepared.bodies)?
    } else {
        export::Exporter::new(model, options, dec)
            .face_regions(face_regions)
            .build(bodies)?
    };
    let (file, report) = answer.value;
    out.write_all(file.as_bytes())
        .map_err(|error| catallax_numbers::Refusal::new("step.io-error", error.to_string()))?;
    Ok(catallax_numbers::Answer {
        value: report,
        marks: answer.marks,
    })
}

/// Writes a certified body with rolling-envelope faces. `max_deviation` bounds
/// the file against the represented model; compose it with the model's
/// definition certificate. A bound over the requested tolerance writes nothing.
pub fn write_rolling_junction_step<W: Write>(
    model: &Model,
    body: BodyId,
    options: &CertifiedSplineStepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    let mut bytes = Vec::new();
    let answer = write_step(
        model,
        &[body],
        &StepOptions {
            product_name: options.product_name.clone(),
            curve_tolerance: options.file_tolerance.clone(),
            significant_digits: options.decimal_places.saturating_add(2),
            seams: Seams::Legacy,
        },
        dec,
        &mut bytes,
    )?;
    if answer.value.max_deviation > options.file_tolerance {
        return Err(catallax_numbers::Refusal::new(
            "step.tolerance-not-reached",
            format!(
                "rounded geometry reaches {} mm, exceeding the requested {} mm file tolerance",
                answer.value.max_deviation, options.file_tolerance
            ),
        ));
    }
    out.write_all(&bytes)
        .map_err(|error| catallax_numbers::Refusal::new("step.io-error", error.to_string()))?;
    Ok(answer)
}

/// Writes a massless surface envelope and separately selected unit-density
/// surrogate solids. This opt-in representation contains no material metadata.
/// Ordinary [`write_step`] remains the solid-only exchange.
pub fn write_stand_in_step<W: Write>(
    model: &Model,
    envelope: BodyId,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    let answer = export::Exporter::new(model, options, dec).build_stand_in(envelope, bodies)?;
    let (file, report) = answer.value;
    if options.seams == Seams::Required {
        seams::check(&file)?;
    }
    out.write_all(file.as_bytes())
        .map_err(|error| catallax_numbers::Refusal::new("step.io-error", error.to_string()))?;
    Ok(catallax_numbers::Answer {
        value: report,
        marks: answer.marks,
    })
}

#[cfg(test)]
mod tests;
