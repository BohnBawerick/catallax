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
mod export;
mod fit;
mod frames;
mod graph;
mod number;
mod partial;
mod pcurve;
mod poles;
mod regularize;
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
}

impl Default for StepOptions {
    fn default() -> Self {
        Self {
            product_name: "catallax".into(),
            curve_tolerance: Rational::from_ratio(1, 10_000_000).unwrap(),
            significant_digits: 17,
        }
    }
}

/// A straight wire gap eliminated only in the exported topology.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct StepRegularization {
    /// Model face whose two decimal-rounded boundary wires overlap.
    pub face: catallax_topology::FaceId,
    /// Outer and contained inner straight edges bounding the collapsed ribbon.
    pub edges: [catallax_topology::EdgeId; 2],
    /// Parallel neighbouring faces joined across that ribbon.
    pub merged_faces: [catallax_topology::FaceId; 2],
    /// Proven whole-ribbon displacement bound in millimetres.
    pub max_deviation: Rational,
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
}

/// Writes `bodies` as one STEP AP242 file to `out`.
pub fn write_step<W: Write>(
    model: &Model,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    let prepared = regularize::prepare(model, bodies, options, dec)?;
    let answer = if let Some(prepared) = prepared {
        export::Exporter::new(&prepared.model, options, dec)
            .regularizations(prepared.records, prepared.face_errors, prepared.edge_errors)
            .build(&prepared.bodies)?
    } else {
        export::Exporter::new(model, options, dec).build(bodies)?
    };
    let (file, report) = answer.value;
    out.write_all(file.as_bytes())
        .map_err(|error| catallax_numbers::Refusal::new("step.io-error", error.to_string()))?;
    Ok(catallax_numbers::Answer {
        value: report,
        marks: answer.marks,
    })
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
    out.write_all(file.as_bytes())
        .map_err(|error| catallax_numbers::Refusal::new("step.io-error", error.to_string()))?;
    Ok(catallax_numbers::Answer {
        value: report,
        marks: answer.marks,
    })
}

#[cfg(test)]
mod tests;
