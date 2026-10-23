//! Independent whole-domain verification of polynomial elliptic-offset sweeps.

use crate::measurement::MeasurementOptions;
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{BSplineSurface, EllipticOffsetDefinition, certify_elliptic_offset_sweep};

/// Two-way boundary-correspondence enclosure for a complete finite offset sweep.
/// This is not a body certificate: cap domains, common edges, material orientation
/// and the association with a particular body's faces remain publication checks.
#[derive(Clone, Debug)]
pub struct OffsetSweepMeasurement {
    /// Non-negative lower bound on the worst two-way boundary distance, in mm.
    pub lower: Rational,
    /// Uniform upper bound, in mm, including every knot span and the closure.
    pub upper: Rational,
    /// Whether the requested width was attained. No bound is silently enlarged.
    pub converged: bool,
    /// Number of nonzero spans checked, in complete parameter order.
    pub spans: usize,
}

/// Measures the complete signed elliptic-offset sweep using the shared checked
/// geometric theorem, independently of the producer's fit/error records.
/// The result does not prove a body's cap domains, incidence or material sense.
pub fn measure_elliptic_offset_sweep(
    surface: &BSplineSurface,
    definition: &EllipticOffsetDefinition,
    axial: [&Real; 2],
    options: &MeasurementOptions,
    dec: &mut Decider,
) -> Result<OffsetSweepMeasurement, Refusal> {
    if options.width <= Rational::zero() || options.max_cells == 0 {
        dec.check_deadline()?;
        return Err(Refusal::new(
            "field.offset-spline",
            "positive measurement width and one streaming cell required",
        ));
    }
    let proof =
        certify_elliptic_offset_sweep(surface, definition, axial, dec).map_err(|error| {
            if error.code == "shapes.offset-sweep" {
                Refusal::new("field.offset-spline", error.detail)
            } else {
                error
            }
        })?;
    let spans = proof.span_count();
    let upper = proof.into_bound().upper;
    Ok(OffsetSweepMeasurement {
        lower: Rational::zero(),
        converged: upper <= options.width,
        upper,
        spans,
    })
}
