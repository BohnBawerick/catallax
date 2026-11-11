// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the operations return new owned curves instead of
// editing a curve in place, statuses become `Fault`, and SMLib's notification
// calls and borrowed-curve bookkeeping are dropped. The zone tolerance is a
// parameter because SMLib reads it from the curve's edge or the system.
// Translated from SMLib SmBSplineCurve.cpp (CreateCanonical,
// EditParameterization), SmBSplineCurveNLib.cpp (Trim, the version compiled
// without THIS_IS_THE_OLD_ALGORITHM...), SmNurbsCrv.cpp (sm_FixupKnotVector)
// and SmTol.cpp (MapTo1d) at c9c5979f04cea8c38213004d4f04df0dcfb40d05;
// rust-brep has no port of these.
//! Building a B-spline curve from STEP data, changing its parameter range,
//! and trimming it.

use crate::Fault;
use crate::curve_cache::{is_clamped, trim as cut_curve};
use crate::extent::Extent1;
use crate::iges_curve::{RawCurve, iges_126_curve};
use crate::scalar::{EFF_ZERO, are_same, scaled_zero};
use crate::spline::{BSplineCurve, NO_Z, unique_knots};
use crate::vector::Vector3;

/// Why [`BSplineCurve::trim`] produced no curve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrimRefusal {
    /// SMLib's `SM_ERR` or `SM_ERR_INVALID_INPUT`.
    Failed(Fault),
    /// The cut leaves a span shorter than the zone tolerance. SMLib removes
    /// it with `N_CrvRemoveDegenSegs`, which a later port step brings.
    ShortSpan,
}

impl From<Fault> for TrimRefusal {
    fn from(fault: Fault) -> Self {
        Self::Failed(fault)
    }
}

/// A trimmed curve and the interval it was cut to.
#[derive(Clone, Debug, PartialEq)]
pub struct Trimmed {
    /// The curve over `interval`, or the whole curve when nothing was cut.
    pub curve: BSplineCurve,
    /// The requested interval after snapping to existing knots.
    pub interval: Extent1,
    /// Whether the curve was cut.
    pub cut: bool,
}

/// `sm_FixupKnotVector` at `SM_EFF_ZERO`: knots closer than the scaled zero
/// become equal. SMLib's second pass, which spreads distinct knots closer
/// than the tolerance, cannot act here: distinct knots differ by at least
/// the scaled zero, which is at least `SM_EFF_ZERO`.
fn fixup_knots(knots: &mut [f64]) {
    let scaled = scaled_zero(knots[0]).max(scaled_zero(knots[knots.len() - 1]));
    for i in 0..knots.len() - 1 {
        if knots[i + 1] - knots[i] < scaled {
            knots[i + 1] = knots[i];
        }
    }
}

/// `SmTol::MapTo1d` for a curve: the zone tolerance divided by the speed,
/// or a hundredth of the parameter range where the speed vanishes.
fn map_to_1d(zone: f64, derivative: Vector3, natural: Extent1) -> f64 {
    let speed = derivative.length();
    if speed < EFF_ZERO {
        natural.length() / 100.0
    } else {
        zone / speed
    }
}

impl BSplineCurve {
    /// `SmBSplineCurve::CreateCanonical`: the STEP form of a curve.
    ///
    /// `points` are Euclidean, with [`NO_Z`] in every `z` for a planar
    /// curve; `knots` are the distinct knots with their `multiplicities`;
    /// `weights` are one per point. A `trim` interval replaces the natural
    /// one. A three-dimensional request turns planar points into points on
    /// the `z = 0` plane. A dimension outside 1 to 3, a degree outside 1 to
    /// 32, or storage whose counts do not agree is [`Fault::InvalidInput`]
    /// (SMLib reads past its arrays when the counts disagree). Degree above
    /// 17, a weight outside `[0.001, 100]`, or a trim outside the knot range
    /// is [`Fault::Degenerate`], SMLib's `SM_ERR`. Non-finite coordinates
    /// are also refused here.
    pub fn create_canonical(
        dimension: usize,
        degree: usize,
        points: &[Vector3],
        multiplicities: &[usize],
        knots: &[f64],
        weights: Option<&[f64]>,
        trim: Option<Extent1>,
    ) -> Result<Self, Fault> {
        if !(1..=3).contains(&dimension) || !(1..=32).contains(&degree) {
            return Err(Fault::InvalidInput);
        }
        if weights.is_some_and(|w| w.len() != points.len()) || multiplicities.len() != knots.len() {
            return Err(Fault::InvalidInput);
        }
        // Sum the multiplicities before expanding any of them.
        let total = multiplicities
            .iter()
            .try_fold(0usize, |sum, &count| sum.checked_add(count));
        if points.is_empty() || total != Some(points.len() + degree + 1) {
            return Err(Fault::InvalidInput);
        }
        let expanded: Vec<f64> = knots
            .iter()
            .zip(multiplicities)
            .flat_map(|(&knot, &count)| std::iter::repeat_n(knot, count))
            .collect();
        let interval = trim.map_or(
            [expanded[degree], expanded[expanded.len() - degree - 1]],
            |t| [t.min, t.max],
        );
        let unit = vec![1.0; points.len()];
        let mut raw: RawCurve = iges_126_curve(
            degree,
            &expanded,
            weights.unwrap_or(&unit),
            points,
            interval,
        )?;
        if dimension == 3 && raw.points[0].coordinates[2] == NO_Z {
            for point in &mut raw.points {
                point.coordinates[2] = 0.0;
            }
        }
        fixup_knots(&mut raw.knots);
        Self::new(raw.degree, raw.knots, raw.points)
    }

    /// `SmBSplineCurve::EditParameterization`: the same curve over a new
    /// parameter range. A range within the scaled zero of the current one
    /// changes nothing. A range of no length, or a non-finite one, is
    /// [`Fault::InvalidInput`]; SMLib builds a curve with collapsed knots.
    /// Knots that are not clamped at both ends are [`Fault::InvalidInput`]
    /// too: SMLib's rescaling of the end knots reshapes such a curve.
    pub fn edit_parameterization(&self, range: Extent1) -> Result<Self, Fault> {
        let old = self.natural_interval();
        if are_same(old.min, range.min) && are_same(old.max, range.max) {
            return Ok(self.clone());
        }
        if !is_clamped(self) {
            return Err(Fault::InvalidInput);
        }
        // `N_BasisReparam`.
        let p = self.degree();
        let mut knots = self.knots().to_vec();
        let m = knots.len() - 1;
        if knots[0] != range.min || knots[m] != range.max {
            let first = knots[0];
            let factor = (range.max - range.min) / (knots[m] - knots[0]);
            for knot in &mut knots[..=p] {
                *knot = range.min;
            }
            for knot in &mut knots[p + 1..m - p] {
                *knot = range.min + factor * (*knot - first);
            }
            for knot in &mut knots[m - p..] {
                *knot = range.max;
            }
        }
        Self::new(p, knots, self.control_points().to_vec())
    }

    /// What `SmBSplineCurve::Trim` does with `requested` before it cuts:
    /// `None` where it returns success without changing anything, else the
    /// interval snapped to the curve's knots with the zone tolerance mapped
    /// to parameter space at each end.
    fn snap_trim(&self, requested: Extent1, zone: f64) -> Result<Option<Snapped>, Fault> {
        let natural = self.natural_interval();
        if natural.is_contained_by(requested, EFF_ZERO) {
            return Ok(None);
        }
        let (mut lower, mut upper) = (requested.min, requested.max);
        let tolerance_lower = 100.0 * scaled_zero(lower);
        let tolerance_upper = 100.0 * scaled_zero(upper);
        // A failed evaluation leaves SMLib's undefined derivative, which
        // maps to a hundredth of the range.
        let speed = |u: f64| {
            self.evaluate(u, 1, true, true)
                .map_or_else(|_| Vector3::default(), |values| values[1])
        };
        let zone_lower = map_to_1d(zone, speed(lower), natural);
        let zone_upper = map_to_1d(zone, speed(upper), natural);

        if !natural.contains_value(lower, 0.0) {
            if (natural.min - lower).abs() < tolerance_lower {
                lower = natural.clamp_value(lower);
            } else {
                return Err(Fault::Degenerate);
            }
        }
        if !natural.contains_value(upper, 0.0) {
            if (natural.max - upper).abs() < tolerance_upper {
                upper = natural.clamp_value(upper);
            } else {
                return Err(Fault::Degenerate);
            }
        }

        let unique = unique_knots(self.knots(), None);
        for &(knot, _) in &unique {
            if (knot - lower).abs() < tolerance_lower {
                lower = knot;
            }
            if (knot - upper).abs() < tolerance_upper {
                upper = knot;
                break;
            }
            if knot > upper + tolerance_upper {
                break;
            }
        }

        // A snap that collapses the interval falls back to the request.
        if (lower - upper).abs() < tolerance_lower + tolerance_upper {
            lower = requested.min;
            upper = requested.max;
            if upper - lower < EFF_ZERO {
                return Err(Fault::InvalidInput);
            }
        }

        // Both ends are already at the curve's ends.
        if (lower - unique[0].0).abs() < tolerance_lower
            && (upper - unique[unique.len() - 1].0).abs() < tolerance_upper
        {
            return Ok(None);
        }
        let mut interval = Extent1::default();
        interval.set_min_max(lower, upper)?;
        Ok(Some(Snapped {
            interval,
            zone_lower,
            zone_upper,
            distinct_knots: unique.len(),
        }))
    }

    /// `SmBSplineCurve::Trim` with its notification argument unsure: the
    /// interval the trim would cut to, without cutting. It is `requested`
    /// itself where SMLib returns success without editing the interval.
    pub fn snapped_trim_interval(&self, requested: Extent1, zone: f64) -> Result<Extent1, Fault> {
        Ok(self
            .snap_trim(requested, zone)?
            .map_or(requested, |snapped| snapped.interval))
    }

    /// `SmBSplineCurve::Trim`: the curve restricted to `requested`, which
    /// snaps to nearby knots, with `zone` the zone tolerance in model units
    /// (`1e-5` unless a model sets one).
    ///
    /// An interval that contains the natural one, or that snaps to it,
    /// leaves the curve whole. A cut that leaves a span shorter than the
    /// mapped tolerance is [`TrimRefusal::ShortSpan`].
    pub fn trim(&self, requested: Extent1, zone: f64) -> Result<Trimmed, TrimRefusal> {
        let Some(snapped) = self.snap_trim(requested, zone)? else {
            return Ok(Trimmed {
                curve: self.clone(),
                interval: requested,
                cut: false,
            });
        };
        // `sm_CrvTrim` assumes clamped knots.
        if !is_clamped(self) {
            return Err(Fault::InvalidInput.into());
        }
        let interval = snapped.interval;
        let curve = cut_curve(self, interval.min, interval.max).map_err(|_| Fault::Degenerate)?;
        let knots = curve.knots();
        let (p, last) = (curve.degree(), knots.len() - 1);
        if snapped.distinct_knots > 2
            && (knots[p + 1] - knots[0] < snapped.zone_lower
                || knots[last] - knots[last - p - 1] < snapped.zone_upper)
        {
            return Err(TrimRefusal::ShortSpan);
        }
        Ok(Trimmed {
            curve,
            interval,
            cut: true,
        })
    }
}

/// A trim interval snapped to knots, with what the cut needs afterwards.
struct Snapped {
    interval: Extent1,
    zone_lower: f64,
    zone_upper: f64,
    distinct_knots: usize,
}
