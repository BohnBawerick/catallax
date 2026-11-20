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

/// The tolerance stored for one reading (decision Q6, section 4.2 steps 7
/// and 8): 1.1 times the reported gap.
///
/// Precision escalates in one case only: the stored tolerance is over the
/// limit, and it would be at or under it if the binary64 rounding error of
/// the measurement were taken off. `enclose` then gives the gap at the
/// reported parameter at 128 bits, and the upper end of that enclosure
/// replaces the binary64 value. What the samples left unresolved stays in the
/// tolerance; a higher precision does not remove it.
fn settled(
    found: &gap::GapMeasure,
    limit: &Rational,
    enclose: impl FnOnce() -> Result<Rational, Refusal>,
) -> Result<Rational, Refusal> {
    let stored = exact(found.stored_mm())?;
    let unresolved = exact((found.error_mm - found.noise_mm).max(0.0))?;
    let least = margin() * (exact(found.raw_mm - found.noise_mm)? + unresolved.clone());
    if stored <= *limit || least > *limit {
        return Ok(stored);
    }
    Ok(margin() * (enclose()? + unresolved))
}

/// Settles every reading of one edge use and returns the one with the largest
/// stored tolerance, with the parameter it was read at. A reading that looks
/// largest in binary64 can fall under the limit at higher precision while
/// another stays over it, so none is dropped before it is settled.
fn largest_settled<'a, T: Copy>(
    candidates: &'a [(T, gap::GapMeasure)],
    mut settle: impl FnMut(&T, &'a gap::GapMeasure) -> Result<(Rational, Rational), Refusal>,
) -> Result<(T, Rational, Rational), Refusal> {
    let mut largest: Option<(T, Rational, Rational)> = None;
    for (which, found) in candidates {
        let (parameter, stored) = settle(which, found)?;
        if largest.as_ref().is_none_or(|(.., best)| stored > *best) {
            largest = Some((*which, parameter, stored));
        }
    }
    largest.ok_or_else(|| bad("trim", "an edge use has no gap reading"))
}

/// Measures one edge use and returns the record admission stores for it: the
/// parameter of the largest gap and the tolerance of rule T4, 1.1 times the
/// reported gap. `NumericalEvidence` then compares that tolerance with the
/// body's limit, as it compared every sample before.
///
/// A plane-projected use has two readings, the projected curve and its two
/// declared ends. Both are settled by [`settled`] and the larger is stored.
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
    let found = gap::measure(&copies, &Hints::default(), Limits::default(), dec)?;
    if let Some(exit) = found.domain_exit {
        return Err(bad(
            "trim",
            format!(
                "the face copy of edge {:?} leaves the surface of face {:?} at parameter {exit}",
                edge.key, face.key
            ),
        ));
    }
    let mut candidates = vec![(copy, found)];
    if projected {
        // A projected copy has no u/v of its own between its ends, but its
        // two declared end values are part of the trim and are measured. It
        // is the same edge use, so this second reading draws on what the
        // first one left of the use's one allowance.
        let declared = EdgeCopies {
            face_copy: FaceCopyRef::Linear(&trim.uv),
            ..copies
        };
        let left = Limits {
            work: Limits::default().work.saturating_sub(candidates[0].1.work),
        };
        candidates.push((
            FaceCopyRef::Linear(&trim.uv),
            gap::measure_ends(&declared, left, dec)?,
        ));
    }
    // Every reading goes through the decision near the limit on its own, and
    // the largest result is stored with its own parameter.
    let (decided_on, parameter, distance_upper_mm) =
        largest_settled(&candidates, |copy, found| {
            let parameter = exact(found.parameter)?
                .max(ends[0].min(ends[1]).clone())
                .min(ends[0].max(ends[1]).clone());
            let upper = |dec: &mut Decider| upper_at(face, edge, copy, ends, &parameter, dec);
            let stored = settled(found, limit, || upper(dec))?;
            Ok((parameter, stored))
        })?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Rational;

    fn q(text: &str) -> Rational {
        text.parse().expect("a rational")
    }

    fn reading(raw_mm: f64, noise_mm: f64, parameter: f64) -> gap::GapMeasure {
        gap::GapMeasure {
            raw_mm,
            error_mm: noise_mm,
            noise_mm,
            parameter,
            domain_exit: None,
            method: gap::Method::EndsOnly,
            plan: Vec::new(),
            refinements: Vec::new(),
            evaluations: 0,
            work: 0,
        }
    }

    /// Two readings of one edge use, both over a 0.01 mm limit in binary64
    /// and both within their rounding error of it. The first has the larger
    /// rounding error, so it looks larger. At 128 bits it is 1e-16 mm under
    /// the gap that stores exactly the limit, and the second is 1e-16 mm
    /// over. The second must be the one stored, with its own parameter.
    #[test]
    fn a_reading_that_drops_under_the_limit_does_not_hide_one_that_stays_over() {
        let limit = q("0.01");
        let at_limit = q("1/110");
        let raw = 1.0 / 110.0;
        let candidates = [
            ("projected", reading(raw, 1e-13, 0.25)),
            ("declared ends", reading(raw, 1e-14, 1.0)),
        ];
        // In binary64 the first is the larger and both are over the limit.
        let stored: Vec<_> = candidates
            .iter()
            .map(|(_, found)| exact(found.stored_mm()).expect("finite"))
            .collect();
        assert!(stored[0] > stored[1] && stored[1] > limit);
        let mut enclosed = Vec::new();
        let (which, parameter, value) = largest_settled(&candidates, |which, found| {
            let value = settled(found, &limit, || {
                enclosed.push(*which);
                Ok(if *which == "projected" {
                    &at_limit - &q("1e-16")
                } else {
                    &at_limit + &q("1e-16")
                })
            })?;
            Ok((exact(found.parameter)?, value))
        })
        .expect("settled");
        // Both were decided at the higher precision, and the one still over
        // the limit is stored.
        assert_eq!(enclosed, ["projected", "declared ends"]);
        assert_eq!(which, "declared ends");
        assert_eq!(parameter, q("1"));
        assert!(value > limit);
        assert_eq!(value, margin() * (&at_limit + &q("1e-16")));
    }

    /// A reading that is over the limit by more than its rounding error is
    /// stored as measured, and one that is under it is too: neither is
    /// enclosed.
    #[test]
    fn a_reading_clear_of_the_limit_is_not_escalated() {
        let limit = q("0.01");
        for raw in [0.005, 0.02] {
            let found = reading(raw, 1e-13, 0.5);
            let value = settled(&found, &limit, || panic!("escalated")).expect("settled");
            assert_eq!(value, exact(found.stored_mm()).expect("finite"));
        }
    }
}
