//! Piece 7 of catallax: STEP export.
//!
//! Writes bodies as STEP AP242 files and proves the largest distance between
//! the exported and the exact geometry.
//!
//! The contract is `docs/pieces/07-step-export.md`.

use std::io::Write;

use catallax_numbers::{Answered, Decider, Rational};
use catallax_topology::{BodyId, Model};

mod export;
mod fit;
mod graph;
mod number;

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

/// What an export wrote and how far it is from the exact model.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct StepReport {
    /// Proven upper bound, in millimetres, on the distance between exported and exact geometry.
    pub max_deviation: Rational,
    /// Intersection branches written as fitted splines.
    pub approximated_curves: usize,
    /// Vertices added because STEP needs them (ring edges, poles).
    pub added_vertices: usize,
    /// Edges added because STEP needs them (seams of loopless faces).
    pub added_edges: usize,
    /// Solids written.
    pub solids: usize,
    /// Faces written.
    pub faces: usize,
    /// Edges written.
    pub edges: usize,
    /// Vertices written.
    pub vertices: usize,
}

/// Writes `bodies` as one STEP AP242 file to `out`.
pub fn write_step<W: Write>(
    model: &Model,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    let answer = export::Exporter::new(model, options, dec).build(bodies)?;
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
