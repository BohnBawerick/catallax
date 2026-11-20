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

/// The work allowance of one edge use across an import admission, with what
/// the checker has read of the use so far.
///
/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." An import
/// admission does more than one thing with a use. Admission reads its
/// same-parameter gap. The gap rule for imports (decision Q9) reads its
/// geometric gap and may give its face copy the edge's parameter. The body
/// is then admitted with the new copies. All of that takes its work from
/// this one account, and each step is given what the steps before it left.
///
/// Where a stored tolerance is over the cap, rule T3's refit and the
/// admissions around it draw on the same account.
///
/// A reading that is kept here is not taken again. A later step of the same
/// admission uses it as it stands, so no copy is read twice and no reading
/// starts a new cap.
#[derive(Clone, Debug)]
pub(super) struct UseAccount {
    /// The one allowance of the use.
    limit: u64,
    /// What every reading so far has counted against it.
    spent: u64,
    /// The checker's readings of the copy the use has now: the copy itself
    /// and, for a projected copy, its two declared ends. `None` before the
    /// use is read, and for a derived copy, which is not measured.
    readings: Option<Vec<gap::GapMeasure>>,
    /// The readings of the copy as imported, for an admission that goes
    /// back to it after a new copy was tried.
    imported: Option<Vec<gap::GapMeasure>>,
}

impl UseAccount {
    /// The account of a use that has not been read, with the allowance
    /// `allowance`.
    pub(super) fn new(allowance: Limits) -> Self {
        Self {
            limit: allowance.work,
            spent: 0,
            readings: None,
            imported: None,
        }
    }

    /// One account per edge use of `definition`, in definition order, each
    /// with the allowance `each_use`.
    pub(super) fn for_uses(definition: &NumericalDefinition, each_use: Limits) -> Vec<Self> {
        let uses = definition
            .faces
            .iter()
            .flat_map(|face| &face.loops)
            .map(Vec::len)
            .sum();
        vec![Self::new(each_use); uses]
    }

    /// What the use has left of its allowance.
    pub(super) fn left(&self) -> Limits {
        Limits {
            work: self.limit.saturating_sub(self.spent),
        }
    }

    /// The work counted against the allowance so far.
    pub(super) fn spent(&self) -> u64 {
        self.spent
    }

    /// Counts work that a reading took from [`Self::left`].
    pub(super) fn spend(&mut self, work: u64) {
        self.spent = self.spent.saturating_add(work);
    }

    /// The checker's reading of the copy the use has now, where one is
    /// kept.
    pub(super) fn reading(&self) -> Option<&gap::GapMeasure> {
        self.readings.as_ref().and_then(|kept| kept.first())
    }

    /// Counts and keeps readings of the copy as imported that were taken
    /// from [`Self::left`].
    pub(super) fn keep(&mut self, readings: &[gap::GapMeasure]) {
        for reading in readings {
            self.spend(reading.work);
        }
        self.readings = Some(readings.to_vec());
        self.imported = self.readings.clone();
    }

    /// Keeps `reading` as the reading of a new copy of the use, in place of
    /// the old copy's. Its work was counted when it was taken.
    pub(super) fn replace(&mut self, reading: gap::GapMeasure) {
        self.readings = Some(vec![reading]);
    }

    /// Goes back to the copy as imported: its readings are the kept ones
    /// again. Only a use whose copy was replaced changes. What was counted
    /// stays counted.
    pub(super) fn back_to_imported(&mut self) {
        self.readings = self.imported.clone();
    }
}

/// The edge's two parameters in the order this use runs through them.
pub(super) fn traversal<'a>(trim: &NumericalTrim, edge: &'a NumericalEdge) -> [&'a Rational; 2] {
    if trim.sense == Sense::Same {
        edge.parameters.each_ref()
    } else {
        [&edge.parameters[1], &edge.parameters[0]]
    }
}

/// The face copy of one edge use: the retained native curve where there is
/// one, the projection for a circle or spline on a plane, and the straight
/// u/v segment otherwise.
fn face_copy<'a>(
    face: &NumericalFace,
    trim: &'a NumericalTrim,
    edge: &NumericalEdge,
    native: Option<&'a BSplineCurve>,
) -> FaceCopyRef<'a> {
    match native {
        Some(curve) => FaceCopyRef::Native(curve),
        None if analytic::projected_edge(face, edge) => FaceCopyRef::Projected,
        None => FaceCopyRef::Linear(&trim.uv),
    }
}

/// The checker's reading of one edge use against its face copy, under the
/// checker's cap, for the audit of derived copies. Admission does not read
/// through here: it takes its readings in [`stored_sample`].
fn audit_reading<'a>(
    face: &'a NumericalFace,
    trim: &'a NumericalTrim,
    edge: &'a NumericalEdge,
    native: Option<&'a BSplineCurve>,
    dec: &mut Decider,
) -> Result<(FaceCopyRef<'a>, gap::GapMeasure), Refusal> {
    let copy = face_copy(face, trim, edge, native);
    let copies = EdgeCopies {
        curve: &edge.curve,
        parameters: traversal(trim, edge),
        surface: GapSurface::Kernel(&face.surface),
        face_copy: copy,
    };
    let found = gap::measure(&copies, &Hints::default(), Limits::default(), dec)?;
    Ok((copy, found))
}

/// The checker's readings of one measured edge use under `limits`, the use's
/// one allowance: the copy itself and, for a projected copy, its two
/// declared ends.
///
/// A projected copy has no u/v of its own between its ends, but its two
/// declared end values are part of the trim and are measured. It is the same
/// edge use, so that second reading draws on what the first one left. A copy
/// that leaves its surface's domain is refused by the caller, and its ends
/// are not read.
fn take_readings(
    copies: &EdgeCopies<'_>,
    trim: &NumericalTrim,
    projected: bool,
    limits: Limits,
    dec: &mut Decider,
) -> Result<Vec<gap::GapMeasure>, Refusal> {
    let found = gap::measure(copies, &Hints::default(), limits, dec)?;
    if !projected || found.domain_exit.is_some() {
        return Ok(vec![found]);
    }
    let declared = EdgeCopies {
        face_copy: FaceCopyRef::Linear(&trim.uv),
        ..*copies
    };
    let left = Limits {
        work: limits.work.saturating_sub(found.work),
    };
    let ends = gap::measure_ends(&declared, left, dec)?;
    Ok(vec![found, ends])
}

/// Refits the face copy of one edge use from its 3D curve (rules T3, T7 and
/// T8) until rule T4's stored tolerance is at most `target_mm`. `None` for a
/// curve projected on a plane, which has no u/v copy of its own.
///
/// `account` is the use's account in an import admission (rule T3). The
/// refit then draws on what the use has left of its one allowance, and its
/// work is counted there with or without a new copy. The refit of an own
/// build passes `None`.
pub(super) fn refit_use(
    face: &NumericalFace,
    trim: &NumericalTrim,
    edge: &NumericalEdge,
    native: Option<&BSplineCurve>,
    target_mm: f64,
    account: Option<&mut UseAccount>,
    dec: &mut Decider,
) -> Result<Option<gap::refit::Refit>, Refusal> {
    let copy = match native {
        Some(curve) => FaceCopyRef::Native(curve),
        None if analytic::projected_edge(face, edge) => return Ok(None),
        None => FaceCopyRef::Linear(&trim.uv),
    };
    let copies = EdgeCopies {
        curve: &edge.curve,
        parameters: traversal(trim, edge),
        surface: GapSurface::Kernel(&face.surface),
        face_copy: copy,
    };
    match account {
        None => gap::refit::refit(&copies, target_mm, dec),
        Some(account) => {
            let fit = gap::refit::refit_within(&copies, target_mm, account.left(), dec)?;
            account.spend(fit.work);
            Ok(fit.copy)
        }
    }
}

/// One derived copy as the full check reads it.
#[derive(Clone, Debug)]
pub struct DerivedCopyAudit {
    /// Key of the face of the use.
    pub face: String,
    /// Key of the edge.
    pub edge: String,
    /// The rule the copy's tag names.
    pub rule: CopyRule,
    /// The bound the tag states, in millimetres.
    pub bound_mm: Rational,
    /// The largest gap the checker found, in millimetres.
    pub measured_mm: Rational,
    /// The measurement error of that reading, in millimetres.
    pub error_mm: Rational,
}

impl DerivedCopyAudit {
    /// Whether the measured gap is inside the stated bound plus the
    /// measurement error. A copy that fails this carries a wrong tag.
    pub fn holds(&self) -> bool {
        self.measured_mm <= &self.bound_mm + &self.error_mm
    }
}

/// The audit of decision Q8 (section 4.4): runs the full check on every
/// derived copy of a body, the copies admission did not measure.
///
/// Returns one row per derived copy. A caller fails when a row does not
/// [`hold`](DerivedCopyAudit::holds). Plane-projected uses are left out:
/// admission measures them whatever their tag says.
pub fn audit_derived_copies(
    body: &NumericalBody,
    dec: &mut Decider,
) -> Result<Vec<DerivedCopyAudit>, Refusal> {
    let definition = body.definition();
    let index = native::Index::new(
        definition,
        body.native_trims(),
        body.sampled_carriers(),
        dec,
    )?;
    let mut rows = Vec::new();
    for (face_index, face) in definition.faces.iter().enumerate() {
        for (boundary, row) in face.loops.iter().enumerate() {
            for (use_index, trim) in row.iter().enumerate() {
                let edge = &definition.edges[trim.edge];
                let native = index.get(face_index, boundary, use_index);
                let tag = index
                    .tag(face_index, boundary, use_index)
                    .unwrap_or(&trim.tag);
                let CopyTag::Derived { rule, bound, .. } = tag else {
                    continue;
                };
                if native.is_none() && analytic::projected_edge(face, edge) {
                    continue;
                }
                let (_, found) = audit_reading(face, trim, edge, native, dec)?;
                rows.push(DerivedCopyAudit {
                    face: face.key.clone(),
                    edge: edge.key.clone(),
                    rule: *rule,
                    bound_mm: bound.clone(),
                    measured_mm: exact(found.raw_mm)?,
                    error_mm: exact(found.error_mm)?,
                });
            }
        }
    }
    Ok(rows)
}

/// The record admission stores for one edge use, and whether it is the bound
/// of a derived copy (`true`) or a measurement (`false`).
///
/// A measured record holds the parameter of the largest gap and the
/// tolerance of rule T4, 1.1 times the reported gap. `NumericalEvidence` then
/// compares that tolerance with the body's limit, as it compared every
/// sample before.
///
/// A plane-projected use has two readings, the projected curve and its two
/// declared ends. Both are settled by [`settled`] and the larger is stored.
///
/// `account` is the use's account in an import admission. The readings are
/// then taken under what the use has left of its one allowance, counted and
/// kept. Where the account already keeps the readings of this copy, they
/// are used as they stand and the checker does not read the use again.
/// Every other admission passes `None` and reads under the checker's cap.
#[allow(clippy::too_many_arguments)]
pub(super) fn stored_sample(
    face: &NumericalFace,
    trim: &NumericalTrim,
    edge: &NumericalEdge,
    native: Option<&BSplineCurve>,
    tag: &CopyTag,
    limit: &Rational,
    account: Option<&mut UseAccount>,
    dec: &mut Decider,
) -> Result<(bool, NumericalSample), Refusal> {
    let ends = traversal(trim, edge);
    let projected = native.is_none() && analytic::projected_edge(face, edge);
    if let (CopyTag::Derived { bound, .. }, false) = (tag, projected) {
        // A derived copy is not measured. Its rule states the bound, and
        // the bound is what is stored: the margin of rule T4 is for fresh
        // measurements only. A projected copy is always measured, whatever
        // its tag says, because that reading is the only check that the 3D
        // curve lies on the plane.
        if *bound < Rational::zero() {
            return Err(bad(
                "definition",
                format!(
                    "edge {:?} states a negative bound for a derived copy",
                    edge.key
                ),
            ));
        }
        // The record names the range of the edge the bound covers.
        return Ok((
            true,
            NumericalSample {
                parameters: ends.map(Clone::clone).to_vec(),
                distance_upper_mm: bound.clone(),
            },
        ));
    }
    let copy = face_copy(face, trim, edge, native);
    let copies = EdgeCopies {
        curve: &edge.curve,
        parameters: ends,
        surface: GapSurface::Kernel(&face.surface),
        face_copy: copy,
    };
    let readings = match account
        .as_deref()
        .and_then(|account| account.readings.clone())
    {
        Some(kept) => kept,
        None => {
            let limits = account
                .as_deref()
                .map_or_else(Limits::default, UseAccount::left);
            let taken = take_readings(&copies, trim, projected, limits, dec)?;
            if let Some(account) = account {
                account.keep(&taken);
            }
            taken
        }
    };
    if let Some(exit) = readings.first().and_then(|found| found.domain_exit) {
        return Err(bad(
            "trim",
            format!(
                "the face copy of edge {:?} leaves the surface of face {:?} at parameter {exit}",
                edge.key, face.key
            ),
        ));
    }
    // The copy itself, then the declared ends of a projected copy.
    let candidates: Vec<(FaceCopyRef<'_>, gap::GapMeasure)> = [copy, FaceCopyRef::Linear(&trim.uv)]
        .into_iter()
        .zip(readings)
        .collect();
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
    Ok((
        false,
        NumericalSample {
            parameters,
            distance_upper_mm,
        },
    ))
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
