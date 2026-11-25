//! Separate approximate STEP artifacts for checked numerical bodies.
//! No exact answer, complete bound, model extraction or healing is provided.
//!
//! Check once (roadmap item 3.0A, step 5): [`build`] writes the numbers and
//! records the body's largest entity tolerance plus the rounding error of
//! the written numbers. It does not admit the written values again.
//! [`audit_written`] is that second admission, for tests and CI.
use crate::graph::{Field, Graph, Id, boolean, refs, string, text};
use catallax_numbers::{
    Decider, NumericalCheck, NumericalEvidence, Rational, Real, Refusal, Sign,
    operation_budget_active, with_operation_budget,
};
use catallax_shapes::{Curve, Point3, Surface, Vector3};
use catallax_topology::{
    BodyKind, Sense,
    numerical::{
        DEFAULT_WORK, NumericalBody, NumericalDefinition, NumericalNativeTrim, NumericalOptions,
        RefitRow,
        tolerance::{self, EntityTolerance, FeatureBand, FeaturePair},
    },
};

pub use catallax_topology::numerical::tolerance::ReceiverProfile;

mod analytic;
mod curved;
mod rounding;
pub use rounding::{RoundingError, RoundingErrors, RoundingLevers};

/// Controls fixed before a numerical export. They do not change source tolerance.
#[derive(Clone, Debug)]
pub struct NumericalStepOptions {
    /// STEP product name, independent of the retained source name.
    pub product_name: String,
    /// Significant digits, from 1 through 17 in this first writer.
    pub significant_digits: u32,
    /// Maximum decimal error of each written length scalar, in millimetres.
    /// This is not a whole-body approximation bound.
    pub max_scalar_rounding_error_mm: Rational,
    /// The receiving system, for rule T8's export with a receiver profile.
    pub receiver: Option<ReceiverProfile>,
    /// A scaled export (rule T8, for micro parts): the same numbers written
    /// under a larger length unit, so a receiver reads a part this many
    /// times larger.
    pub scale: Option<LengthScale>,
}
impl Default for NumericalStepOptions {
    fn default() -> Self {
        Self {
            product_name: "catallax numerical approximate".into(),
            significant_digits: 17,
            max_scalar_rounding_error_mm: Rational::from_ratio(1, 10_000_000).unwrap(),
            receiver: None,
            scale: None,
        }
    }
}

/// The factor of a scaled export (rule T8). The written numbers stay the
/// millimetre values; the file's length unit becomes the centimetre, the
/// decimetre or the metre, so a receiver reads every length, and the file's
/// uncertainty, that many times larger. Nothing is rounded again.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum LengthScale {
    /// Ten times: centimetres.
    Ten,
    /// A hundred times: decimetres.
    Hundred,
    /// A thousand times: metres.
    Thousand,
}

impl LengthScale {
    /// The factor.
    pub fn factor(self) -> i64 {
        match self {
            Self::Ten => 10,
            Self::Hundred => 100,
            Self::Thousand => 1000,
        }
    }
    /// The STEP `SI_UNIT` prefix of the written length unit.
    fn prefix(self) -> &'static str {
        match self {
            Self::Ten => ".CENTI.",
            Self::Hundred => ".DECI.",
            Self::Thousand => "$",
        }
    }
    /// Every scale, smallest first.
    pub const ALL: [Self; 3] = [Self::Ten, Self::Hundred, Self::Thousand];
}

/// A small feature as a receiver will see it (rule T8).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReceiverFeature {
    /// Two entities that share no topology, with their distance in the
    /// model's millimetres.
    pub pair: FeaturePair,
    /// The distance the receiver reads: the model's distance times the
    /// export's scale, in the receiver's millimetres.
    pub receiver_distance_mm: Rational,
    /// Under 10 times, or under 4 times, the receiver's tolerance.
    pub band: FeatureBand,
}

/// An entity whose tolerance stays over the receiver's after the refit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReceiverResidual {
    /// The entity, with its tolerance in the model's millimetres.
    pub entity: EntityTolerance,
    /// The tolerance the receiver reads: the model's times the export's
    /// scale.
    pub receiver_tolerance_mm: Rational,
}

/// What rule T8 did for a receiver profile.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReceiverReport {
    /// The profile.
    pub profile: ReceiverProfile,
    /// Edge uses over the receiver's tolerance whose face copies were
    /// refitted, with their stored tolerances before and after in the
    /// receiver's millimetres.
    pub refits: Vec<RefitRow>,
    /// The scale the receiver reads the file at: 1, or the export's scale.
    pub scale: i64,
    /// Entities still over the receiver's tolerance after the refit, read at
    /// the export's scale.
    pub residuals: Vec<ReceiverResidual>,
    /// Features under 10 times and 4 times the receiver's tolerance.
    pub features: Vec<ReceiverFeature>,
    /// Offered for a micro part, one with a feature under 4 times the
    /// receiver's tolerance: the smallest scale that lifts every feature to
    /// at least 10 times it, or the largest scale where none does.
    pub scaled_export_offer: Option<LengthScale>,
}
/// Rounding of native spline parameters and weights, separate from lengths.
#[derive(Clone, Debug, Default)]
pub struct NumericalSplineRounding {
    /// Maximum error in a knot, finite curve parameter or patch UV coordinate,
    /// measured in that native parameter's units. No spatial bound is inferred.
    pub max_parameter_error: Rational,
    /// Maximum error of a dimensionless positive control weight.
    pub max_weight_error: Rational,
}
/// Source and file checks with separate scopes; always a numerical result.
#[derive(Clone, Debug)]
pub struct NumericalStepReport {
    /// The retained body kind. An open surface is never reported as a solid.
    pub kind: BodyKind,
    /// Original generating input, before any explicit sewing or later edits.
    pub original_definition: catallax_topology::numerical::NumericalDefinition,
    /// Native parameter curves in the original generating input.
    pub original_native_trims: Vec<catallax_topology::numerical::NumericalNativeTrim>,
    /// Native parameter curves checked on the current source.
    pub source_native_trims: Vec<catallax_topology::numerical::NumericalNativeTrim>,
    /// Native parameter curves reconstructed from the emitted decimal values.
    pub emitted_native_trims: Vec<catallax_topology::numerical::NumericalNativeTrim>,
    /// Earlier checked states. Their checks apply only to their own definitions.
    pub source_history: Vec<catallax_topology::numerical::NumericalHistory>,
    /// Retained generating definition name.
    pub source_name: String,
    /// Immutable copy of the generating tables, trims and stable source keys.
    pub source_definition: catallax_topology::numerical::NumericalDefinition,
    /// N01/N02 evidence, including tolerance, actual samples and repair history.
    /// Sample maxima are never promoted into a complete bound.
    pub source_evidence: NumericalEvidence,
    /// How far the rounding of the written numbers can move a point of the
    /// written model, to first order. The written values are not measured.
    pub rounding: RoundingError,
    /// Requested error limit on written length scalars only.
    pub scalar_rounding_limit_mm: Rational,
    /// Actual maximum rounding error among written length scalars only.
    pub max_scalar_rounding_error_mm: Rational,
    /// Maximum rounding error of normalized DIRECTION ratios, dimensionless.
    /// This is not an angular or spatial deviation bound.
    pub max_direction_ratio_rounding_error: Rational,
    /// Stored faces, explicit edges and vertices, with no inferred topology.
    pub counts: [usize; 3],
    /// Native spline rounding, when the curved writer is used.
    pub spline_rounding: Option<NumericalSplineRounding>,
    /// Checks performed by this writer, distinct from external reference checks.
    pub checks: Vec<String>,
    /// The tolerance the file records, in millimetres before any scale: the
    /// largest entity tolerance of the body plus the rounding error, rounded
    /// up to the written digits. The stamp carries it as
    /// `recorded_tolerance_mm`, and reopening trusts it. Rule T8 also writes
    /// it as the file's uncertainty when the export has a receiver profile.
    pub uncertainty_mm: Rational,
    /// Rule T8's report, when a receiver profile was given.
    pub receiver: Option<ReceiverReport>,
    /// The scale written, if any.
    pub scale: Option<LengthScale>,
}
/// An approximate file, with no implicit conversion to a proven/bounded result.
/// The private bytes were checked before this value was returned.
/// Raw reference import and edit acceptance require an independent check.
///
/// ```compile_fail
/// use catallax_step::{StepReport, numerical::NumericalStep};
/// use catallax_numbers::Answer;
/// fn promote(file: NumericalStep) -> Answer<StepReport> { file.into() }
/// ```
#[derive(Clone, Debug)]
pub struct NumericalStep {
    bytes: Vec<u8>,
    report: NumericalStepReport,
    /// The full check of the written values, kept only by [`audit_written`].
    written: Option<NumericalEvidence>,
}

/// The full check of the values a file holds, against the tolerance the
/// file records. For tests and CI: an export does not run it.
#[derive(Clone, Debug)]
pub struct WrittenAudit {
    /// The file, byte for byte what [`build`] writes for the same body.
    pub file: NumericalStep,
    /// Admission's checks of the written decimal values, with the checks of
    /// the written u/v curves of an analytic face. Each stored record is
    /// 1.1 times a gap measured on the written values (rule T4). None is
    /// over 1.1 times the recorded tolerance, or the audit would have
    /// refused: no measured gap is over what the file records.
    pub evidence: NumericalEvidence,
}
impl NumericalStep {
    /// Approximate STEP bytes; the file itself states the numerical category.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
    /// Source provenance and scoped file checks, never an exact STEP report.
    pub fn report(&self) -> &NumericalStepReport {
        &self.report
    }
    /// Checks bytes reread from the actual output, including complete byte identity
    /// to this checked artifact and the seam syntax contract.
    /// This is not an OCCT import, healing, geometric certification or proof.
    pub fn check_written(&self, reread: &[u8], dec: &mut Decider) -> Result<(), Refusal> {
        guarded(|| {
            if reread.len() != self.bytes.len() {
                return Err(bad("file-changed", "reread file length differs"));
            }
            dec.charge_work(reread.len() as u64)?;
            if reread != self.bytes {
                return Err(bad("file-changed", "reread file bytes differ"));
            }
            crate::check_seams(reread)?;
            dec.check_deadline()?;
            Ok(())
        })
    }
}
fn bad(kind: &str, detail: impl Into<String>) -> Refusal {
    let code = match kind {
        "file-changed" => "step.numerical-file-changed",
        "file-number" => "step.numerical-file-number",
        "rounding-topology" => "step.numerical-rounding-topology",
        "rounding-limit" => "step.numerical-rounding-limit",
        "unsupported" => "step.numerical-unsupported",
        "options" => "step.numerical-options",
        "gap" => "step.numerical-gap",
        _ => unreachable!(),
    };
    Refusal::new(code, detail)
}
fn guarded<T>(run: impl FnOnce() -> Result<T, Refusal>) -> Result<T, Refusal> {
    if operation_budget_active() {
        return run();
    }
    let (result, work) = with_operation_budget(DEFAULT_WORK, run);
    if work.exhausted {
        Err(Refusal::new(
            "numbers.operation-budget",
            "numerical STEP exhausted its work allowance",
        ))
    } else {
        result
    }
}
fn rational_point(p: &Point3) -> Result<[Rational; 3], Refusal> {
    [&p.x, &p.y, &p.z]
        .map(|v| {
            v.as_rational().ok_or_else(|| {
                bad(
                    "unsupported",
                    "first numerical writer needs rational coordinates",
                )
            })
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| unreachable!())
}
fn rational_vector(p: &Vector3) -> Result<[Rational; 3], Refusal> {
    rational_point(&Point3::new(p.x.clone(), p.y.clone(), p.z.clone()))
}
fn axis(v: &[Rational; 3], unit: bool) -> Option<usize> {
    let nonzero: Vec<_> = v.iter().enumerate().filter(|(_, q)| !q.is_zero()).collect();
    if nonzero.len() != 1 || (unit && nonzero[0].1.abs() != Rational::one()) {
        return None;
    }
    Some(nonzero[0].0)
}
fn cross(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
// Parse the selected STEP decimal token independently. Rational::parse needs a
// digit after '.', while STEP permits a trailing decimal point before E/end.
fn written(q: &Rational, options: &NumericalStepOptions) -> Result<Rational, Refusal> {
    let mut token = crate::number::rational(q, options.significant_digits)
        .text
        .replace(".E", ".0E");
    if token.ends_with('.') {
        token.push('0');
    }
    Rational::parse(&token).map_err(|_| bad("file-number", "written decimal cannot be reread"))
}
fn point_from(p: &[Rational; 3]) -> Point3 {
    Point3::new(
        p[0].clone().into(),
        p[1].clone().into(),
        p[2].clone().into(),
    )
}

fn arithmetic(
    a: &Rational,
    b: &Rational,
    op: char,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, matches!(op, '*' | '/'))?;
    Ok(match op {
        '+' => a + b,
        '-' => a - b,
        '*' => a * b,
        '/' => a / b,
        _ => unreachable!(),
    })
}
fn dot(a: &[Rational; 3], b: &[Rational; 3], dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut result = Rational::zero();
    for i in 0..3 {
        let product = arithmetic(&a[i], &b[i], '*', dec)?;
        result = arithmetic(&result, &product, '+', dec)?;
    }
    Ok(result)
}
fn difference(
    a: &[Rational; 3],
    b: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    Ok([
        arithmetic(&a[0], &b[0], '-', dec)?,
        arithmetic(&a[1], &b[1], '-', dec)?,
        arithmetic(&a[2], &b[2], '-', dec)?,
    ])
}
fn charged_cross(
    a: &[Rational; 3],
    b: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    let mut result = std::array::from_fn(|_| Rational::zero());
    for (i, q) in result.iter_mut().enumerate() {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        let x = arithmetic(&a[j], &b[k], '*', dec)?;
        let y = arithmetic(&a[k], &b[j], '*', dec)?;
        *q = arithmetic(&x, &y, '-', dec)?;
    }
    Ok(result)
}
struct DirectionValues {
    values: [Rational; 3],
    tokens: [String; 3],
    error: Rational,
}
// STEP DIRECTION ratios are dimensionless and normalized by the reader. Scale
// them before decimal selection; never put their error in the length report.
fn direction_values(
    v: &[Rational; 3],
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<DirectionValues, Refusal> {
    if axis(v, false).is_some() {
        let values = v.each_ref().map(|q| {
            if q.is_zero() {
                Rational::zero()
            } else if q.sign() == Sign::Positive {
                Rational::one()
            } else {
                -Rational::one()
            }
        });
        let tokens = values.each_ref().map(|q| {
            if q.is_zero() {
                "0."
            } else if q.sign() == Sign::Positive {
                "1."
            } else {
                "-1."
            }
            .to_owned()
        });
        return Ok(DirectionValues {
            values,
            tokens,
            error: Rational::zero(),
        });
    }
    let mut scale = Rational::zero();
    for q in v {
        dec.reserve_uncounted_arithmetic(q, &scale, true)?;
        scale = scale.max(q.abs());
    }
    if scale.is_zero() {
        return Err(bad("rounding-topology", "direction vanished"));
    }
    let mut result = DirectionValues {
        values: std::array::from_fn(|_| Rational::zero()),
        tokens: std::array::from_fn(|_| String::new()),
        error: Rational::zero(),
    };
    for (i, q) in v.iter().enumerate() {
        let ratio = arithmetic(q, &scale, '/', dec)?;
        dec.charge_work(1)?;
        let number = crate::number::rational(&ratio, options.significant_digits);
        dec.reserve_uncounted_arithmetic(&result.error, &number.error, true)?;
        result.error = result.error.max(number.error);
        result.values[i] = written(&ratio, options)?;
        result.tokens[i] = number.text;
    }
    Ok(result)
}
fn inverse_uv(
    plane: &catallax_shapes::Plane,
    p: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let origin = rational_point(&plane.origin)?;
    let x = rational_vector(&plane.ref_dir)?;
    let y = charged_cross(&rational_vector(&plane.normal)?, &x, dec)?;
    let delta = difference(p, &origin, dec)?;
    let mut uv = [Rational::zero(), Rational::zero()];
    for (i, axis) in [&x, &y].into_iter().enumerate() {
        let squared = dot(axis, axis, dec)?;
        if squared.sign() != Sign::Positive {
            return Err(bad("rounding-topology", "written plane frame vanished"));
        }
        let value = &Real::from(dot(&delta, axis, dec)?) / &Real::from(squared).sqrt();
        uv[i] = dec.enclose_budgeted_at(&value, 128)?.mid();
    }
    Ok(uv)
}
// Same IDs/loops, and only values actually selected for serialization. The UVs
// below are inverse coordinates on emitted planar frames, not inferred topology.
// An audit's work: an export does not admit the written values.
fn planar_payload(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalDefinition, Refusal> {
    let mut d = body.definition().clone();
    // STEP carries no copy tag and both copies are rounded below, so the
    // check of the written values measures every use.
    d.reset_copy_tags();
    d.name.push_str(" STEP decimal payload");
    for vertex in &mut d.vertices {
        dec.charge_work(1)?;
        for q in &mut vertex.position {
            *q = written(q, options)?;
        }
    }
    for (edge, source) in d.edges.iter_mut().zip(&body.definition().edges) {
        dec.charge_work(1)?;
        let Curve::Line(line) = &mut edge.curve else {
            return Err(bad("unsupported", "finite line required"));
        };
        let mut origin = rational_point(&line.origin)?;
        let mut dir = rational_vector(&line.dir)?;
        for q in &mut origin {
            *q = written(q, options)?;
        }
        if let Some(a) = axis(&dir, false) {
            let magnitude = written(&dir[a].abs(), options)?;
            dir[a] = if dir[a].sign() == Sign::Negative {
                -magnitude
            } else {
                magnitude
            };
            if dir[a].is_zero() {
                return Err(bad("rounding-topology", "line direction vanished"));
            }
            edge.parameters = std::array::from_fn(|i| {
                (&d.vertices[edge.vertices[i]].position[a] - &origin[a]) / dir[a].clone()
            });
        } else {
            // Unit magnitude in STEP and these retained rational ratios define
            // the same oriented line carrier. Its finite parameters are local
            // inverse coordinates; EDGE_CURVE keeps the explicit vertex IDs.
            dir = direction_values(&dir, options, dec)?.values;
            let squared = dot(&dir, &dir, dec)?;
            for i in 0..2 {
                let delta = difference(&d.vertices[edge.vertices[i]].position, &origin, dec)?;
                edge.parameters[i] = arithmetic(&dot(&delta, &dir, dec)?, &squared, '/', dec)?;
            }
        }
        let written_order = edge.parameters[0].cmp(&edge.parameters[1]);
        if written_order != std::cmp::Ordering::Equal
            && written_order != source.parameters[0].cmp(&source.parameters[1])
        {
            return Err(bad(
                "rounding-topology",
                "written vertices reverse the retained finite line",
            ));
        }
        line.origin = point_from(&origin);
        line.dir = Vector3::new(
            dir[0].clone().into(),
            dir[1].clone().into(),
            dir[2].clone().into(),
        );
    }
    for face in &mut d.faces {
        dec.charge_work(1)?;
        let Surface::Plane(plane) = &mut face.surface else {
            return Err(bad("unsupported", "plane required"));
        };
        let mut origin = rational_point(&plane.origin)?;
        for q in &mut origin {
            *q = written(q, options)?;
        }
        plane.origin = point_from(&origin);
        let mut x = rational_vector(&plane.ref_dir)?;
        let normal = rational_vector(&plane.normal)?;
        let coordinate_frame = axis(&normal, true).is_some() && axis(&x, true).is_some();
        if !coordinate_frame {
            let normal = direction_values(&normal, options, dec)?.values;
            let reference = direction_values(&x, options, dec)?.values;
            // gp_Ax2 projects the written reference onto the written normal's
            // perpendicular plane. Reconstruct that actual STEP frame.
            let y = charged_cross(&normal, &reference, dec)?;
            x = charged_cross(&y, &normal, dec)?;
            if dot(&x, &x, dec)?.is_zero() {
                return Err(bad("rounding-topology", "written plane axes are parallel"));
            }
            plane.normal = Vector3::new(
                normal[0].clone().into(),
                normal[1].clone().into(),
                normal[2].clone().into(),
            );
            plane.ref_dir = Vector3::new(
                x[0].clone().into(),
                x[1].clone().into(),
                x[2].clone().into(),
            );
        }
        let y = cross(&rational_vector(&plane.normal)?, &x);
        for trim in face.loops.iter_mut().flatten() {
            let edge = &d.edges[trim.edge];
            let endpoints = if trim.sense == Sense::Same {
                edge.vertices
            } else {
                [edge.vertices[1], edge.vertices[0]]
            };
            for (i, endpoint) in endpoints.into_iter().enumerate() {
                let p = &d.vertices[endpoint].position;
                trim.uv[i] = if coordinate_frame {
                    [&x, &y].map(|axis| {
                        (0..3).fold(Rational::zero(), |sum, k| {
                            sum + (&p[k] - &origin[k]) * axis[k].clone()
                        })
                    })
                } else {
                    inverse_uv(plane, p, dec)?
                };
            }
        }
    }
    Ok(d)
}

/// What an audit admits: the decimal values a writer wrote.
enum Payload {
    /// Lines and planes: the values are rebuilt from the body when asked.
    Planar,
    /// The definition and u/v curves a curved writer wrote from, with the
    /// checks of what it wrote beside them.
    Curved {
        definition: NumericalDefinition,
        native_trims: Vec<NumericalNativeTrim>,
        checks: Vec<NumericalCheck>,
    },
}

/// The tolerance a file records and the rounding error inside it.
struct Recorded {
    rounding: RoundingError,
    /// The largest entity tolerance plus the rounding error, rounded up to
    /// the written digits.
    tolerance_mm: Rational,
    /// Its decimal text, as the stamp holds it.
    text: String,
    /// The rounding error's decimal text, rounded up.
    rounding_text: String,
}

/// Check once: the body's largest entity tolerance plus the rounding error
/// of the written numbers. Rule T4's margin is in the entity tolerance
/// already and is not applied again.
fn recorded(
    body: &NumericalBody,
    written: (&NumericalDefinition, &[NumericalNativeTrim]),
    errors: RoundingErrors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<Recorded, Refusal> {
    let levers = RoundingLevers::of(written.0, written.1, dec)?;
    let rounding = RoundingError::bound(errors, levers)?;
    let exact = largest_entity_tolerance(body, dec)? + rounding.bound_mm.clone();
    let text = exact.to_decimal_bounds(options.significant_digits).1;
    let tolerance_mm = Rational::parse(&text)
        .map_err(|_| bad("file-number", "recorded tolerance cannot be reread"))?;
    let rounding_text = rounding
        .bound_mm
        .to_decimal_bounds(options.significant_digits)
        .1;
    Ok(Recorded {
        rounding,
        tolerance_mm,
        text,
        rounding_text,
    })
}

/// The body's largest entity tolerance. A body read back from a handover
/// written before the tolerance report has no report, so its largest is read
/// from its checks and rule T1's floor, as the report would have read it.
fn largest_entity_tolerance(body: &NumericalBody, dec: &mut Decider) -> Result<Rational, Refusal> {
    let report = body.tolerances();
    if !report.entities.is_empty() {
        return Ok(report.largest_mm());
    }
    let floor = tolerance::BodyFrame::of(body.definition(), dec)?.floor_mm();
    Ok(tolerance::largest_recorded(body.evidence()).max(floor))
}

/// The record of the tolerance: a description string of its own, after the
/// stamp.
///
/// The stamp names the body and the receiver, which are a caller's text, so
/// a reader cannot tell a field there from a name that looks like one. This
/// string holds two numbers this writer formats and nothing else. A reader
/// takes the record from it alone (`catallax_import::recorded`). The values
/// are in the file's length unit: a scaled export writes the model's
/// millimetre values under a larger unit, and a reader takes the record that
/// many times larger with every other length.
fn recorded_stamp(recorded: &Recorded) -> String {
    format!(
        "catallax recorded tolerance;rounding_error_mm={};recorded_tolerance_mm={}",
        recorded.rounding_text, recorded.text
    )
}

/// An audit's limit. A stored record is 1.1 times a measured gap (rule T4),
/// so a record over 1.1 times the recorded tolerance is a gap over it.
fn audit_limit(recorded_mm: &Rational) -> NumericalOptions {
    // Source: rule T4 of the Checkpoint 2 decisions, section 4.3, as
    // `catallax_topology::numerical::gap::MARGIN` holds it.
    let margin = Rational::from_ratio(11, 10).expect("nonzero");
    NumericalOptions {
        tolerance_mm: margin * recorded_mm.clone(),
        vertex_tolerance_mm: None,
    }
}

/// An audit: admits the written decimal values under the audit's limit.
fn admit_written(
    body: &NumericalBody,
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    extra: Vec<NumericalCheck>,
    limit: &NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalEvidence, Refusal> {
    let build = if body.kind() == BodyKind::Surface {
        catallax_topology::numerical::build_surface_with_sampled_carriers
    } else {
        catallax_topology::numerical::build_with_sampled_carriers
    };
    let admitted = build(
        definition,
        native_trims,
        body.sampled_carriers().clone(),
        limit.clone(),
        dec,
    )?;
    let mut checks = admitted.evidence().checks().to_vec();
    checks.extend(extra);
    NumericalEvidence::new(limit.tolerance_mm.clone(), checks, Vec::new())
}

struct Writer<'a> {
    graph: Graph,
    options: &'a NumericalStepOptions,
    dec: &'a mut Decider,
    error: Rational,
    direction_error: Rational,
    /// The largest rounding error of a line's magnitude, which is the rate
    /// of the line's parameter as well as a length.
    magnitude_error: Rational,
    /// Whether this is [`audit_written`]: the written values are admitted.
    audit: bool,
}
impl Writer<'_> {
    fn scalar(&mut self, q: &Rational) -> Result<Field, Refusal> {
        self.dec.charge_work(1)?;
        let n = crate::number::rational(q, self.options.significant_digits);
        if n.error > self.options.max_scalar_rounding_error_mm {
            return Err(bad(
                "rounding-limit",
                "written scalar exceeds its separate rounding limit",
            ));
        }
        self.error = self.error.clone().max(n.error);
        Ok(text(n.text))
    }
    fn point(&mut self, q: &[Rational; 3]) -> Result<Id, Refusal> {
        let xyz = q
            .iter()
            .map(|q| self.scalar(q))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self
            .graph
            .add("CARTESIAN_POINT", vec![text("''"), Field::List(xyz)]))
    }
    fn direction(&mut self, q: &[Rational; 3]) -> Result<Id, Refusal> {
        let values = direction_values(q, self.options, self.dec)?;
        self.direction_error = self.direction_error.clone().max(values.error);
        Ok(self.graph.add(
            "DIRECTION",
            vec![
                text("''"),
                Field::List(values.tokens.into_iter().map(text).collect()),
            ],
        ))
    }
}
/// Builds a separate approximate STEP artifact for a checked numerical body.
/// Planar lines, native spline edges and rational patches retain their carriers
/// and explicit physical loops. Unsupported carriers refuse. No topology is
/// inferred and no vertex is healed.
///
/// Check once: the numbers are written and the file records one tolerance,
/// the body's largest entity tolerance plus the rounding error of the
/// written numbers ([`RoundingError`]). The stamp holds it as
/// `recorded_tolerance_mm`. The written values are not admitted again, so a
/// file is written whatever its rounding does to a gap; the recorded
/// tolerance says how much that can be. [`audit_written`] measures it.
/// With a receiver profile the file's global uncertainty is that same
/// tolerance (rule T8); without one the file has none.
/// The caller's allowance covers all work; no global policy is changed.
///
/// With a receiver profile (rule T8), the face copies of edge uses over the
/// receiver's tolerance are refitted first, the entities still over it are
/// listed, features under 10 times and 4 times its tolerance are listed, and
/// a micro part is offered a scaled export.
pub fn build(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    write(body, options, false, dec)
}

/// Writes the file as [`build`] does and then admits the values it wrote: the
/// full check of every vertex, edge and edge use of the decimal model, and of
/// the written u/v curves of an analytic face. It refuses with
/// `numbers.numerical-tolerance-exceeded` where a gap measured on the written
/// values is over the tolerance the file records. It is a full admission, so
/// any other rule can refuse it too: rule T5 does on a body whose coarse
/// rounding leaves its entities too loose for its features.
///
/// This is the admission an export ran before check once. It costs more than
/// the export itself and is for tests and CI, as the audit of derived copies
/// is ([`catallax_topology::numerical::audit_derived_copies`]).
pub fn audit_written(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<WrittenAudit, Refusal> {
    let mut file = write(body, options, true, dec)?;
    let evidence = file
        .written
        .take()
        .expect("an audit admits the written values");
    Ok(WrittenAudit { file, evidence })
}

fn write(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    audit: bool,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    guarded(|| {
        let (prepared, receiver) = match &options.receiver {
            Some(profile) => {
                let scale = options.scale.map_or(1, LengthScale::factor);
                let (prepared, report) = for_receiver(body, profile, scale, dec)?;
                (std::borrow::Cow::Owned(prepared), Some(report))
            }
            None => (std::borrow::Cow::Borrowed(body), None),
        };
        let mut artifact = if analytic::applies(&prepared) {
            analytic::build(&prepared, options, audit, dec)?
        } else {
            build_inner(&prepared, options, audit, dec)?
        };
        artifact.report.receiver = receiver;
        Ok(artifact)
    })
}

/// Rule T8 before writing: refit edges above the receiver's tolerance, list
/// what stays above it and the features under 10 and 4 times it, and offer
/// a scaled export for a micro part.
///
/// The receiver reads the file at `scale`: every length the model has, times
/// the scale. So its limits are taken in the model's millimetres divided by
/// the scale for the refit and the feature search, and what is reported is
/// read at the scale, in the receiver's millimetres. The offer is the
/// smallest scale that lifts every feature to 10 times the receiver's
/// tolerance, read from the model's own distances.
fn for_receiver(
    body: &NumericalBody,
    profile: &ReceiverProfile,
    scale: i64,
    dec: &mut Decider,
) -> Result<(NumericalBody, ReceiverReport), Refusal> {
    let factor = Rational::from_i64(scale);
    let model_limit = &profile.tolerance_mm / &factor;
    let (refitted, mut refits) = body.refit_face_copies(&model_limit, dec)?;
    for row in &mut refits {
        dec.charge_work(1)?;
        row.before_mm = &row.before_mm * &factor;
        row.after_mm = &row.after_mm * &factor;
    }
    let residuals = refitted
        .tolerances()
        .entities
        .iter()
        .filter(|entity| entity.tolerance_mm > model_limit)
        .map(|entity| ReceiverResidual {
            entity: entity.clone(),
            receiver_tolerance_mm: &entity.tolerance_mm * &factor,
        })
        .collect();
    let reach = &model_limit * &Rational::from_i64(tolerance::RECEIVER_WIDE_BAND);
    let mut features = Vec::new();
    for pair in tolerance::features_within(&refitted, &reach, dec)? {
        let receiver_distance_mm = &pair.distance_mm * &factor;
        if let Some(band) = tolerance::feature_band(&receiver_distance_mm, profile) {
            features.push(ReceiverFeature {
                pair,
                receiver_distance_mm,
                band,
            });
        }
    }
    let smallest = features
        .iter()
        .filter(|feature| feature.band == FeatureBand::UnderFour)
        .map(|feature| feature.pair.distance_mm.clone())
        .min();
    let scaled_export_offer = smallest.map(|size| {
        LengthScale::ALL
            .into_iter()
            .find(|scale| {
                tolerance::feature_band(&(&size * &Rational::from_i64(scale.factor())), profile)
                    .is_none()
            })
            .unwrap_or(LengthScale::Thousand)
    });
    Ok((
        refitted,
        ReceiverReport {
            profile: profile.clone(),
            refits,
            scale,
            residuals,
            features,
            scaled_export_offer,
        },
    ))
}

/// Serialize a numerical file with its stamp, recorded tolerance and units.
fn render_file(
    graph: &mut Graph,
    body: &NumericalBody,
    options: &NumericalStepOptions,
    solid: Id,
    description: [&str; 2],
    uncertainty: &str,
) -> String {
    use crate::document::{self, BREP, MILLI, Preamble};
    let prefix = options.scale.map_or(MILLI, LengthScale::prefix);
    // Rule T8 writes the uncertainty as part of an export with a receiver
    // profile. Without one the header stays as it was.
    let uncertainty = options.receiver.as_ref().map(|_| uncertainty);
    let representation = if body.kind() == BodyKind::Surface {
        "MANIFOLD_SURFACE_SHAPE_REPRESENTATION"
    } else {
        BREP
    };
    let preamble = Preamble::add(
        graph,
        &options.product_name,
        representation,
        &[solid],
        uncertainty,
        prefix,
        false,
    );
    let layout = preamble.layout(graph, &[solid], &[]);
    document::render(graph, &layout, &options.product_name, &description)
}

/// The stamp's note of a scale and a receiver.
fn stamp_extras(options: &NumericalStepOptions) -> String {
    let mut out = String::new();
    if let Some(scale) = options.scale {
        out.push_str(&format!(
            ";scaled_export=x{};written_values=millimetres",
            scale.factor()
        ));
    }
    if let Some(profile) = &options.receiver {
        out.push_str(&format!(
            ";receiver={};receiver_tolerance_mm={}",
            profile.name, profile.tolerance_mm
        ));
    }
    out
}
fn build_inner(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    audit: bool,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    if options.product_name.trim().is_empty()
        || !(1..=17).contains(&options.significant_digits)
        || options.max_scalar_rounding_error_mm.sign() != Sign::Positive
    {
        return Err(bad(
            "options",
            "named product, 1..17 digits and positive scalar rounding limit required",
        ));
    }
    let d = body.definition();
    let count = body.counts();
    if count.shells != 1 || count.lumps != 1 {
        return Err(bad("unsupported", "one retained closed shell required"));
    }
    if !body.native_trims().is_empty()
        || d.edges
            .iter()
            .any(|edge| !matches!(edge.curve, Curve::Line(_)))
        || d.faces
            .iter()
            .any(|face| !matches!(face.surface, Surface::Plane(_)))
    {
        return curved::build(body, options, audit, dec);
    }
    let mut w = Writer {
        graph: Graph::default(),
        options,
        dec,
        error: Rational::zero(),
        direction_error: Rational::zero(),
        magnitude_error: Rational::zero(),
        audit,
    };
    let mut vertices = Vec::new();
    for vertex in &d.vertices {
        let point = w.point(&vertex.position)?;
        vertices.push(w.graph.add(
            "VERTEX_POINT",
            vec![text(string(&vertex.key)), point.into()],
        ));
    }
    let mut edges = Vec::new();
    for edge in &d.edges {
        w.dec.charge_work(1)?;
        let Curve::Line(line) = &edge.curve else {
            return Err(bad(
                "unsupported",
                "only retained finite line edges are supported",
            ));
        };
        let origin = rational_point(&line.origin)?;
        let dir = rational_vector(&line.dir)?;
        let point = w.point(&origin)?;
        let direction = w.direction(&dir)?;
        let magnitude = if let Some(a) = axis(&dir, false) {
            w.scalar(&dir[a].abs())?
        } else {
            w.scalar(&Rational::one())?
        };
        let vector = w
            .graph
            .add("VECTOR", vec![text("''"), direction.into(), magnitude]);
        let line = w.graph.add(
            "LINE",
            vec![text(string(&edge.key)), point.into(), vector.into()],
        );
        edges.push(w.graph.add(
            "EDGE_CURVE",
            vec![
                text(string(&edge.key)),
                vertices[edge.vertices[0]].into(),
                vertices[edge.vertices[1]].into(),
                line.into(),
                boolean(edge.parameters[0] < edge.parameters[1]),
            ],
        ));
    }
    let mut faces = Vec::new();
    for face in &d.faces {
        w.dec.charge_work(1)?;
        let Surface::Plane(plane) = &face.surface else {
            return Err(bad(
                "unsupported",
                "numerical curved-face export is outside this bounded family",
            ));
        };
        let origin = rational_point(&plane.origin)?;
        let normal = rational_vector(&plane.normal)?;
        let x = rational_vector(&plane.ref_dir)?;
        let coordinate_frame = axis(&normal, true).is_some() && axis(&x, true).is_some();
        if coordinate_frame && axis(&normal, true) == axis(&x, true) {
            return Err(bad("unsupported", "plane axes must be orthogonal"));
        }
        let point = w.point(&origin)?;
        let n = w.direction(&normal)?;
        let reference = w.direction(&x)?;
        let placement = w.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), point.into(), n.into(), reference.into()],
        );
        let surface = w
            .graph
            .add("PLANE", vec![text(string(&face.key)), placement.into()]);
        let mut bounds = Vec::new();
        for (index, boundary) in face.loops.iter().enumerate() {
            let mut uses = Vec::new();
            for trim in boundary {
                w.dec.charge_work(1)?;
                uses.push(w.graph.add(
                    "ORIENTED_EDGE",
                    vec![
                        text("''"),
                        text("*"),
                        text("*"),
                        edges[trim.edge].into(),
                        boolean(trim.sense == Sense::Same),
                    ],
                ));
            }
            let loop_ = w.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]);
            bounds.push(w.graph.add(
                if index == 0 {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                },
                vec![text("''"), loop_.into(), boolean(true)],
            ));
        }
        faces.push(w.graph.add(
            "ADVANCED_FACE",
            vec![
                text(string(&face.key)),
                refs(bounds),
                surface.into(),
                boolean(face.sense == Sense::Same),
            ],
        ));
    }
    let solid = w.representation(faces, body)?;
    w.finish(
        solid,
        body,
        Payload::Planar,
        None,
        vec![
            "source.retained-line-and-plane-carriers".into(),
            "file.rounding-error-added-to-tolerance".into(),
            "file.canonical-byte-identity".into(),
            "file.planar-seam-contract".into(),
            "file.approximate-category-no-global-bound".into(),
        ],
    )
}

impl Writer<'_> {
    fn representation(&mut self, faces: Vec<Id>, body: &NumericalBody) -> Result<Id, Refusal> {
        if body.kind() == BodyKind::Surface {
            self.dec.charge_work(body.definition().edges.len() as u64)?;
            let mut uses = vec![0usize; body.definition().edges.len()];
            for trim in body
                .definition()
                .faces
                .iter()
                .flat_map(|f| f.loops.iter().flatten())
            {
                self.dec.charge_work(1)?;
                uses[trim.edge] += 1;
            }
            let shell = self.graph.add(
                if uses.iter().all(|n| *n == 2) {
                    "CLOSED_SHELL"
                } else {
                    "OPEN_SHELL"
                },
                vec![text("''"), refs(faces)],
            );
            Ok(self.graph.add(
                "SHELL_BASED_SURFACE_MODEL",
                vec![text("'numerical approximate surface'"), refs([shell])],
            ))
        } else {
            let shell = self
                .graph
                .add("CLOSED_SHELL", vec![text("''"), refs(faces)]);
            Ok(self.graph.add(
                "MANIFOLD_SOLID_BREP",
                vec![text("'numerical approximate'"), shell.into()],
            ))
        }
    }

    fn finish(
        mut self,
        solid: Id,
        body: &NumericalBody,
        payload: Payload,
        spline_rounding: Option<NumericalSplineRounding>,
        checks: Vec<String>,
    ) -> Result<NumericalStep, Refusal> {
        let options = self.options;
        let d = body.definition();
        let count = body.counts();
        let mut description = format!(
            "catallax numerical approximate;source={};source_tolerance_mm={};source_checks={};source_repairs={};complete_body_bound=unestablished;max_scalar_rounding_error_mm={}",
            d.name,
            body.options().tolerance_mm,
            body.evidence().checks().len(),
            body.evidence().repairs().len(),
            self.error
        );
        if !self.direction_error.is_zero() {
            description.push_str(&format!(
                ";max_direction_ratio_rounding_error={}",
                self.direction_error
            ));
        }
        if let Some(rounding) = &spline_rounding {
            description.push_str(&format!(
                ";max_native_parameter_rounding_error={};max_weight_rounding_error={}",
                rounding.max_parameter_error, rounding.max_weight_error
            ));
        }
        let errors = RoundingErrors {
            length_mm: self.error.clone(),
            line_magnitude: self.magnitude_error.clone(),
            direction: self.direction_error.clone(),
            parameter: spline_rounding
                .as_ref()
                .map_or_else(Rational::zero, |rounding| {
                    rounding.max_parameter_error.clone()
                }),
            weight: spline_rounding
                .as_ref()
                .map_or_else(Rational::zero, |rounding| rounding.max_weight_error.clone()),
        };
        // A planar file's numbers are its vertices and carriers' own, so the
        // body's tables hold the same sizes as the file.
        let written = match &payload {
            Payload::Planar => (body.definition(), &[][..]),
            Payload::Curved {
                definition,
                native_trims,
                ..
            } => (definition, native_trims.as_slice()),
        };
        let recorded = recorded(body, written, errors, options, self.dec)?;
        description.push_str(&stamp_extras(options));
        let record = recorded_stamp(&recorded);
        let uncertainty_text = crate::number::decimal(&recorded.text);
        let file = render_file(
            &mut self.graph,
            body,
            options,
            solid,
            [&description, &record],
            &uncertainty_text,
        );
        let written = if self.audit {
            let limit = audit_limit(&recorded.tolerance_mm);
            Some(match payload {
                Payload::Planar => admit_written(
                    body,
                    planar_payload(body, options, self.dec)?,
                    Vec::new(),
                    Vec::new(),
                    &limit,
                    self.dec,
                )?,
                Payload::Curved {
                    definition,
                    native_trims,
                    checks,
                } => admit_written(body, definition, native_trims, checks, &limit, self.dec)?,
            })
        } else {
            None
        };
        let artifact = NumericalStep {
            written,
            bytes: file.into_bytes(),
            report: NumericalStepReport {
                kind: body.kind(),
                original_definition: body.source_definition().clone(),
                original_native_trims: body.source_native_trims().to_vec(),
                source_native_trims: body.native_trims().to_vec(),
                emitted_native_trims: Vec::new(),
                source_history: body.history().to_vec(),
                source_name: d.name.clone(),
                source_definition: d.clone(),
                source_evidence: body.evidence().clone(),
                rounding: recorded.rounding,
                scalar_rounding_limit_mm: options.max_scalar_rounding_error_mm.clone(),
                max_scalar_rounding_error_mm: self.error,
                max_direction_ratio_rounding_error: self.direction_error,
                counts: [count.faces, count.edges, count.vertices],
                checks,
                spline_rounding,
                uncertainty_mm: recorded.tolerance_mm,
                receiver: None,
                scale: options.scale,
            },
        };
        artifact.check_written(artifact.bytes(), self.dec)?;
        Ok(artifact)
    }
}
