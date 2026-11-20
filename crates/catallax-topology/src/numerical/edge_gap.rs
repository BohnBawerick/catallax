//! Admission's use of the edge-gap checker: one measured, stored tolerance
//! per edge use.

use super::gap::{self, EdgeCopies, FaceCopyRef, GapSurface, Hints, Limits};
use super::*;
use catallax_numbers::Real;
use catallax_shapes::BSplineCurve;

/// Source: rule T4 of the tolerance rules: the stored tolerance is 1.1 times
/// the measured gap plus the measurement error. As a rational, for the one
/// case decided on an exact value.
fn margin() -> Rational {
    Rational::from_ratio(11, 10).expect("nonzero")
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| bad("trim", "an edge gap is not a finite number"))
}

/// A 128-bit upper bound of the gap at one parameter, from the exact
/// definitions. This is the evaluation admission made at every sample before
/// it measured; it now runs only where a measurement straddles the limit.
fn upper_at(
    face: &NumericalFace,
    edge: &NumericalEdge,
    copy: &FaceCopyRef<'_>,
    ends: [&Rational; 2],
    t: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let p = curve_at(&edge.curve, t, dec)?
        .ok_or_else(|| bad("unsupported", "edge has no parameter evaluator"))?;
    let at: [Real; 2] = match copy {
        FaceCopyRef::Native(curve) => {
            let point = curve.point_at(&Real::from(t.clone()), dec)?;
            [point.x, point.y]
        }
        FaceCopyRef::Linear(uv) => {
            let share = (t - ends[0]) / (ends[1] - ends[0]);
            std::array::from_fn(|axis| {
                Real::from(&uv[0][axis] + &(&share * &(&uv[1][axis] - &uv[0][axis])))
            })
        }
        FaceCopyRef::Projected => curved::parameters(&face.surface, &p, dec)?.map(Real::from),
    };
    let q = native::point(&face.surface, &at, dec)?;
    Ok(distance(&p, &q))
}

/// Measures one edge use and returns the record admission stores for it: the
/// parameter of the largest gap and the tolerance of rule T4, 1.1 times the
/// reported gap. `NumericalEvidence` then compares that tolerance with the
/// body's limit, as it compared every sample before.
///
/// Precision escalates in one case only (decision Q6, section 4.2 step 7):
/// the stored tolerance is over the limit, and it would be at or under it if
/// the binary64 rounding error of the measurement were taken off. The gap at
/// that parameter is then enclosed at 128 bits and the enclosure's upper end
/// replaces the binary64 value. What the samples left unresolved stays in the
/// tolerance; a higher precision does not remove it.
pub(super) fn stored_sample(
    face: &NumericalFace,
    trim: &NumericalTrim,
    edge: &NumericalEdge,
    native: Option<&BSplineCurve>,
    limit: &Rational,
    dec: &mut Decider,
) -> Result<NumericalSample, Refusal> {
    let ends = if trim.sense == Sense::Same {
        edge.parameters.each_ref()
    } else {
        [&edge.parameters[1], &edge.parameters[0]]
    };
    let projected = native.is_none() && analytic::projected_edge(face, edge);
    let copy = match native {
        Some(curve) => FaceCopyRef::Native(curve),
        None if projected => FaceCopyRef::Projected,
        None => FaceCopyRef::Linear(&trim.uv),
    };
    let copies = EdgeCopies {
        curve: &edge.curve,
        parameters: ends,
        surface: GapSurface::Kernel(&face.surface),
        face_copy: copy,
    };
    let mut found = gap::measure(&copies, &Hints::default(), Limits::default(), dec)?;
    let mut decided_on = copy;
    if projected {
        // A projected copy has no u/v of its own between its ends, but its
        // two declared end values are part of the trim and are measured.
        let declared = EdgeCopies {
            face_copy: FaceCopyRef::Linear(&trim.uv),
            ..copies
        };
        let at_ends = gap::measure_ends(&declared, dec)?;
        if at_ends.reported_mm() > found.reported_mm() {
            found = at_ends;
            decided_on = FaceCopyRef::Linear(&trim.uv);
        }
    }
    if let Some(exit) = found.domain_exit {
        return Err(bad(
            "trim",
            format!(
                "the face copy of edge {:?} leaves the surface of face {:?} at parameter {exit}",
                edge.key, face.key
            ),
        ));
    }
    let parameter = exact(found.parameter)?
        .max(ends[0].min(ends[1]).clone())
        .min(ends[0].max(ends[1]).clone());
    let stored = exact(found.stored_mm())?;
    let unresolved = exact((found.error_mm - found.noise_mm).max(0.0))?;
    let least = margin() * (exact(found.raw_mm - found.noise_mm)? + unresolved.clone());
    let distance_upper_mm = if stored <= *limit || least > *limit {
        stored
    } else {
        margin() * (upper_at(face, edge, &decided_on, ends, &parameter, dec)? + unresolved)
    };
    // A straight face copy has exact u/v at the reported parameter; they
    // follow it in the record, as they followed every sample before.
    let mut parameters = vec![parameter.clone()];
    if let FaceCopyRef::Linear(uv) = decided_on {
        let share = (&parameter - ends[0]) / (ends[1] - ends[0]);
        parameters
            .extend((0..2).map(|axis| &uv[0][axis] + &(&share * &(&uv[1][axis] - &uv[0][axis]))));
    }
    Ok(NumericalSample {
        parameters,
        distance_upper_mm,
    })
}
