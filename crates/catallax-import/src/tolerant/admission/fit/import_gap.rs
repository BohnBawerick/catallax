//! The gap rule for imports on one edge use of a native tolerant body
//! (decision Q9, and section 4.2 step 2 of the Checkpoint 2 decisions).
//!
//! "Imports whose copies run at different speeds: measure the geometric gap
//! both ways, reparametrise, then continue with the same-parameter gap."
//!
//! A file's u/v curve is a second copy of its edge, and ISO 10303-42 ties its
//! parameter to the 3D curve's by nothing. This reader evaluates the two at
//! one parameter, so a u/v curve that traces the edge at another speed is far
//! from it at equal parameters with nothing wrong in its geometry.
//!
//! The rule is `catallax_topology::numerical::read_import_use`, the code the
//! numerical import runs on every use. The A+ checker reads the
//! same-parameter gap of the file's copy. Where that could be lowered, the
//! geometric gap is read both ways. Where the same-parameter gap is much
//! larger, the copy is given the parameter of its 3D curve, and the checker's
//! reading of the new copy is what the use continues with.
//!
//! The 3D curve and the surface are the file's own carriers, which the
//! kernel has no form for in general. The checker reads them through binary64
//! evaluators of their native charts
//! ([`Binary64Carrier`](super::super::super::native::Binary64Carrier)), in
//! millimetres about the file's origin. The u/v curve is read as the file
//! wrote it: a B-spline with its own controls and knots, or the straight
//! segment between the two ends of a u/v line.
//!
//! A file's u/v curve is kept as the use's face copy only where the rule read
//! it. These are not read, and the caller never keeps them: it rebuilds the
//! trace from the 3D curve, which is a copy of this reader's own in the
//! edge's parameter, and its warning says why.
//!
//! - a u/v curve that is neither a B-spline with written knots nor a line;
//! - a u/v curve whose interval does not hold the edge's, which needs the
//!   affine map of the interval first (OpenCASCADE's `SameRange`);
//! - a u/v curve whose controls leave the surface's charts;
//! - a use binary64 cannot place a sample of.

use catallax_numbers::Real;
use catallax_shapes::{BSplineCurve, Point3};
use catallax_topology::numerical::{self as n, gap};

use super::super::super::native::Binary64Carrier;
use super::super::ReparametrisedTrace;
use super::*;
use crate::geometry_source::{KnotData, Surface as SourceSurface};

/// What the gap rule read of one edge use.
pub(super) struct UseReading {
    /// The rule's own reading.
    pub rule: n::ImportUseReading,
    /// The copy the rule gave the use in place of the file's, in binary64.
    pub copy: Option<ReparametrisedTrace>,
}

/// The file's u/v curve of an edge use, in a form the checker reads.
enum FileCopy {
    /// A B-spline, with the controls, weights and knots the file wrote.
    Spline(BSplineCurve),
    /// A u/v line: its two ends at the edge's two parameters.
    Straight(Box<[[Rational; 2]; 2]>),
}

fn exact(value: f64) -> Result<Rational, Error> {
    Rational::from_f64_exact(value).ok_or_else(|| {
        invalid(
            "import.tolerant-broken-geometry",
            "non-finite parameter of an edge use",
        )
    })
}

/// Whether a refusal says that work or time ran out. It is returned as it
/// is: the hard work cap of an edge "refuses by name" (section 4.2 step 8),
/// and no other reading stands in for it.
fn resource_limit(refusal: &catallax_numbers::Refusal) -> bool {
    refusal.code == "numerical.gap-work-cap"
        || refusal.code.contains("budget")
        || refusal.code.contains("deadline")
}

/// The polynomial degree of a source curve's pieces, for the checker's point
/// count (section 4.2 step 4). Zero for a curve that has none: the checker
/// then takes its least count, and its coefficient tail adds points.
fn curve_degree(parsed: &ParsedStep, mut entity: u64, dec: &mut Decider) -> Result<usize, Error> {
    for _ in 0..32 {
        match parsed.geometry(entity, dec)? {
            Geometry::Curve(SourceCurve::BSpline(spline)) => return Ok(spline.degree),
            Geometry::Curve(SourceCurve::Line(..) | SourceCurve::Polyline(_)) => return Ok(1),
            Geometry::Curve(
                SourceCurve::SurfaceCurve(parent, ..)
                | SourceCurve::Replica(parent, _)
                | SourceCurve::Trimmed(parent, ..),
            ) => entity = parent,
            _ => return Ok(0),
        }
    }
    Ok(0)
}

/// The sum of a source surface's two polynomial degrees: one for a plane,
/// as the checker counts a kernel plane, and zero for a surface that has
/// none.
fn surface_degree(parsed: &ParsedStep, mut entity: u64, dec: &mut Decider) -> Result<usize, Error> {
    for _ in 0..32 {
        match parsed.geometry(entity, dec)? {
            Geometry::Surface(SourceSurface::BSpline(spline)) => {
                return Ok(spline.degree[0] + spline.degree[1]);
            }
            Geometry::Surface(SourceSurface::Plane(_)) => return Ok(1),
            Geometry::Surface(
                SourceSurface::RectangularTrim(parent, ..)
                | SourceSurface::Replica(parent, _)
                | SourceSurface::CurveBounded(parent, ..),
            ) => entity = parent,
            _ => return Ok(0),
        }
    }
    Ok(0)
}

impl Context<'_> {
    /// The binary64 reading of a source curve, compiled once.
    fn curve64(
        &mut self,
        curve: &NativeCurve,
        dec: &mut Decider,
    ) -> Result<Arc<Binary64Carrier>, Error> {
        if let Some(found) = self.curves64.get(&curve.source()) {
            return Ok(Arc::clone(found));
        }
        let compiled = Arc::new(curve.binary64(dec)?);
        self.curves64.insert(curve.source(), Arc::clone(&compiled));
        Ok(compiled)
    }

    /// The binary64 reading of a source surface, compiled once.
    fn surface64(
        &mut self,
        surface: &NativeSurface,
        dec: &mut Decider,
    ) -> Result<Arc<Binary64Carrier>, Error> {
        if let Some(found) = self.surfaces64.get(&surface.source()) {
            return Ok(Arc::clone(found));
        }
        let compiled = Arc::new(surface.binary64(dec)?);
        self.surfaces64
            .insert(surface.source(), Arc::clone(&compiled));
        Ok(compiled)
    }
}

/// The file's u/v curve as the checker reads it, or why it does not.
fn file_copy(
    parsed: &ParsedStep,
    entity: u64,
    pcurve: &NativeUvCurve,
    range: [f64; 2],
    dec: &mut Decider,
) -> Result<Result<FileCopy, &'static str>, Error> {
    let Geometry::Curve(SourceCurve::Pcurve(_, representation)) = parsed.geometry(entity, dec)?
    else {
        return Ok(Err("its face curve is not a PCURVE"));
    };
    let items = parsed.representation(representation)?.items;
    let [item] = items.as_slice() else {
        return Ok(Err("its PCURVE holds more than one curve"));
    };
    match parsed.geometry(*item, dec)? {
        Geometry::Curve(SourceCurve::BSpline(spline)) => {
            let KnotData::Explicit(multiplicities, values, _) = &spline.knots else {
                return Ok(Err("its u/v spline has no written knots"));
            };
            let knots: Vec<Rational> = multiplicities
                .iter()
                .zip(values)
                .flat_map(|(count, value)| std::iter::repeat_n(value.clone(), *count))
                .collect();
            let mut controls = Vec::with_capacity(spline.controls.len());
            for control in &spline.controls {
                dec.charge_work(1)?;
                let Geometry::Point(uv) = parsed.geometry(*control, dec)? else {
                    return Ok(Err("a control of its u/v spline is not a point"));
                };
                let [u, v] = uv.as_slice() else {
                    return Ok(Err("a control of its u/v spline is not a u/v point"));
                };
                controls.push(Point3::new(
                    Real::from(u.clone()),
                    Real::from(v.clone()),
                    Real::zero(),
                ));
            }
            let weights = match &spline.weights {
                Some(weights) => weights.iter().cloned().map(Real::from).collect(),
                None => vec![Real::one(); controls.len()],
            };
            match BSplineCurve::new(spline.degree, knots, controls, weights, dec) {
                Ok(curve) => Ok(Ok(FileCopy::Spline(curve))),
                Err(refusal) if resource_limit(&refusal) => Err(refusal.into()),
                Err(_) => Ok(Err("its u/v spline is not a clamped B-spline")),
            }
        }
        Geometry::Curve(SourceCurve::Line(..)) => {
            let mut ends = [
                [Rational::zero(), Rational::zero()],
                [Rational::zero(), Rational::zero()],
            ];
            for (end, parameter) in ends.iter_mut().zip(range) {
                let (uv, _) = match pcurve.sample(
                    NativeLocation {
                        chart: 0,
                        parameters: [parameter, 0.0],
                    },
                    dec,
                ) {
                    Ok(sampled) => sampled,
                    Err(error) if crate::tolerant::native::control(&error) => return Err(error),
                    Err(_) => return Ok(Err("its u/v line has no value at an end of the edge")),
                };
                *end = [exact(uv[0])?, exact(uv[1])?];
            }
            Ok(Ok(FileCopy::Straight(Box::new(ends))))
        }
        _ => Ok(Err("its u/v curve is neither a B-spline nor a line")),
    }
}

/// The precision a u/v value that is not a rational is enclosed at, as the
/// native samples are (`SAMPLE_BITS` of the native carriers).
const UV_BITS: u32 = 96;

/// A u/v value of a spline as binary64.
fn rounded(value: &Real, dec: &mut Decider) -> Result<f64, Error> {
    Ok(match value.rational_ref() {
        Some(exact) => exact.to_f64_nearest(),
        None => dec
            .enclose_budgeted_at(value, UV_BITS)?
            .mid()
            .to_f64_nearest(),
    })
}

/// Whether a u/v point lies inside the extent of the surface's charts. A
/// spline's hull holds it (Piegl and Tiller, "The NURBS Book", section 4.2,
/// for positive weights), so a spline whose controls are all inside stays on
/// the surface.
fn inside(extent: &[[Option<f64>; 2]; 2], uv: [f64; 2]) -> bool {
    (0..2).all(|axis| {
        extent[axis][0].is_none_or(|low| uv[axis] >= low)
            && extent[axis][1].is_none_or(|high| uv[axis] <= high)
    })
}

/// One edge use and the file u/v curve it carries.
pub(super) struct FileUse<'a> {
    /// The edge's stations on its 3D curve.
    pub edge: &'a EdgeSamples,
    /// The face's source surface and its native form.
    pub surface: (u64, &'a NativeSurface),
    /// The use's source PCURVE and its native form.
    pub pcurve: (u64, &'a NativeUvCurve),
    /// How far, in millimetres, the edge may run past the interval of its
    /// u/v curve: what the body allows a vertex off its edge's end.
    pub overshoot_mm: f64,
}

/// The gap rule for one edge use that carries a file u/v curve, or why the
/// checker does not read the use.
///
/// The use has one work allowance, the checker's cap for an edge. A use that
/// spends it refuses the body with `numerical.gap-work-cap`. The chart work
/// of every native evaluation is charged on top of the checker's own units.
pub(super) fn read_use(
    context: &mut Context<'_>,
    used: &FileUse<'_>,
    dec: &mut Decider,
) -> Result<Result<UseReading, &'static str>, Error> {
    let FileUse {
        edge,
        surface: (surface_entity, surface),
        pcurve: (pcurve_entity, pcurve),
        overshoot_mm,
    } = *used;
    let range = [
        edge.points[0].0.parameters[0],
        edge.points[STATIONS].0.parameters[0],
    ];
    let width = (range[1] - range[0]).abs();
    if !width.is_finite() || width <= 0.0 {
        return Ok(Err("its edge has no finite parameter range"));
    }
    let curve = match &edge.carrier {
        PhysicalCurve::Native(native) => {
            let reading = context.curve64(native, dec)?;
            let evaluator = Arc::clone(&reading);
            (
                gap::SampledCurve {
                    point: Arc::new(move |t: f64| evaluator.evaluate([t, 0.0])),
                    breaks: reading.borders(0),
                    degree: curve_degree(context.parsed, native.source(), dec)?,
                },
                Some(reading),
            )
        }
        PhysicalCurve::Line { origin, direction } => {
            let (origin, direction) = (*origin, *direction);
            (
                gap::SampledCurve {
                    point: Arc::new(move |t: f64| {
                        let point: [f64; 3] =
                            std::array::from_fn(|axis| origin[axis] + t * direction[axis]);
                        (point, 8.0 * f64::EPSILON * norm(point))
                    }),
                    breaks: Vec::new(),
                    degree: 1,
                },
                None,
            )
        }
        PhysicalCurve::Point(_) => return Ok(Err("its edge is a point")),
    };
    let surface64 = context.surface64(surface, dec)?;
    let extent = [surface64.extent(0), surface64.extent(1)];
    let evaluator = Arc::clone(&surface64);
    let lifted = gap::SampledSurface {
        point: Arc::new(move |uv: [f64; 2]| evaluator.evaluate(uv)),
        lines: [surface64.borders(0), surface64.borders(1)],
        periods: surface64.periods(),
        degree: surface_degree(context.parsed, surface_entity, dec)?,
    };
    let file = match file_copy(context.parsed, pcurve_entity, pcurve, range, dec)? {
        Ok(file) => file,
        Err(why) => return Ok(Err(why)),
    };
    if let FileCopy::Spline(spline) = &file {
        // The edge's range inside the spline's own interval. An end placed
        // on the 3D curve from its vertex may fall past it by what the body
        // allows a vertex off its edge's end, and no farther.
        let [low, high] = spline.domain().map(Rational::to_f64_nearest);
        let ends = [range[0].min(range[1]), range[0].max(range[1])];
        let past = (low - ends[0]).max(ends[1] - high).max(0.0);
        if past > 0.0 {
            let step = width / STATIONS as f64;
            let chord = [(0, 1), (STATIONS, STATIONS - 1)]
                .into_iter()
                .map(|(a, b)| {
                    norm(subtract(
                        edge.points[a].1.point_mm,
                        edge.points[b].1.point_mm,
                    ))
                })
                .fold(0.0, f64::max);
            if past * chord > overshoot_mm * step {
                return Ok(Err(
                    "the interval of its u/v curve does not hold the edge's parameter range",
                ));
            }
        }
        for control in spline.controls() {
            let uv = [rounded(&control.x, dec)?, rounded(&control.y, dec)?];
            if !inside(&extent, uv) {
                return Ok(Err("its u/v curve leaves the charts of its surface"));
            }
        }
        // The spline as read here against this reader's own evaluation of
        // the u/v curve, at the ends and the middle of the edge. They are the
        // same curve where both read the file's values in the same units.
        for share in [0.0, 0.5, 1.0] {
            let parameter = (range[0] + share * (range[1] - range[0])).clamp(low, high);
            let (native, uncertainty) = match pcurve.sample(
                NativeLocation {
                    chart: 0,
                    parameters: [parameter, 0.0],
                },
                dec,
            ) {
                Ok(sampled) => sampled,
                Err(error) if crate::tolerant::native::control(&error) => return Err(error),
                Err(_) => return Ok(Err("its u/v curve has no value inside the edge")),
            };
            let point = spline.point_at(&Real::from(exact(parameter)?), dec)?;
            let read = [rounded(&point.x, dec)?, rounded(&point.y, dec)?];
            let size = native.iter().fold(1.0f64, |size, x| size.max(x.abs()));
            if (0..2).any(|axis| {
                (read[axis] - native[axis]).abs() > uncertainty + 64.0 * f64::EPSILON * size
            }) {
                return Ok(Err(
                    "its u/v spline is not in the parameter units this reader evaluates",
                ));
            }
        }
    }
    if let FileCopy::Straight(ends) = &file
        && !ends
            .iter()
            .all(|uv| inside(&extent, uv.each_ref().map(Rational::to_f64_nearest)))
    {
        return Ok(Err("its u/v curve leaves the charts of its surface"));
    }
    let parameters = [exact(range[0])?, exact(range[1])?];
    let copies = gap::EdgeCopies {
        curve: gap::GapCurve::Sampled(&curve.0),
        parameters: parameters.each_ref(),
        surface: gap::GapSurface::Sampled(&lifted),
        face_copy: match &file {
            FileCopy::Spline(spline) => gap::FaceCopyRef::Native(spline),
            FileCopy::Straight(ends) => gap::FaceCopyRef::Linear(ends),
        },
    };
    let read = n::read_import_use(&copies, gap::Limits::default(), dec);
    // What the native charts charge for their evaluations, with or without
    // an answer. The rule's refusal stays the answer where there is one.
    let charts = curve.1.as_ref().map_or(0, |reading| reading.take_work()) + surface64.take_work();
    let charged = dec.charge_work(charts);
    let rule = match read {
        Ok(rule) => rule,
        Err(refusal) if resource_limit(&refusal) => return Err(refusal.into()),
        Err(_) => return Ok(Err("the checker has no reading of its two copies")),
    };
    charged?;
    if !rule.stored_mm().is_finite() || !rule.parametric.stored_mm().is_finite() {
        return Ok(Err("binary64 cannot place a sample of its two copies"));
    }
    let copy = match rule.new_copy() {
        None => None,
        Some(next) => {
            let mut controls = Vec::with_capacity(next.curve.control_count());
            for control in next.curve.controls() {
                let uv = [rounded(&control.x, dec)?, rounded(&control.y, dec)?];
                if !inside(&extent, uv) {
                    return Ok(Err(
                        "its copy in the edge's parameter leaves the charts of its surface",
                    ));
                }
                controls.push(uv);
            }
            Some(ReparametrisedTrace {
                degree: next.curve.degree(),
                knots: next
                    .curve
                    .knots()
                    .iter()
                    .map(Rational::to_f64_nearest)
                    .collect(),
                controls,
                spans: next.spans,
            })
        }
    };
    Ok(Ok(UseReading { rule, copy }))
}
