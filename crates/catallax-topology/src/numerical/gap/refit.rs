//! New face copies of an edge use, in the edge's own parameter.
//!
//! Two kinds are built here, by the same fit:
//!
//! - A refit from the 3D curve (rules T3, T7 and T8: "refit first", "refit
//!   the 2D copy", "refit edges above the receiver tolerance"). At each node
//!   the 3D point's foot on the surface is found from the old copy as a
//!   seed. The new copy lies under the 3D curve.
//! - A change of parameter of the copy the use has (decision Q9, and section
//!   4.2 step 2 of the Checkpoint 2 decisions: an import whose copies run at
//!   different speeds is reparametrised, and its tolerance is not raised). At
//!   each node the old copy is read where its lift is nearest the 3D point.
//!   The new copy runs along the old one and its two ends stay where they
//!   were. OpenCASCADE's `Approx_SameParameter` pairs the copies the same
//!   way, by projection, keeps the two ends, and asks the pairs to ascend.
//!
//! Either way the new copy is a u/v curve in the edge's parameter, so the
//! gap stays a same-parameter gap. The rate of change at a node comes from a
//! second reading a small step along the edge. Each span between nodes is
//! the cubic Hermite interpolant of the two values and their rates, written
//! as a Bezier piece (Farin, "Curves and Surfaces for CAGD", 5th edition,
//! section 7.5), so the copy is a cubic B-spline with triple interior knots.
//! The nodes are the breaks of the 3D curve and, for a change of parameter,
//! the places where the old copy's knots land on it; each piece is cut into
//! equal spans, and the spans double until the checker reads a stored
//! tolerance under the target or the span limit is reached. Every candidate
//! is measured by the A+ checker, so a copy is kept only on what the checker
//! reads.
//!
//! In an import admission the work of a fit is held to the one allowance of
//! its edge use (section 4.2 step 8: "A hard per-edge work cap refuses by
//! name"). The fit comes after other readings of the same use, so it is
//! handed what they left, and every foot, every reading of the old copy and
//! the checker's work on every candidate is drawn from that remainder
//! (`Drawn`). That holds for the change of parameter and for the refit of
//! rule T3 alike. The fit reports what it took, with or without a new copy,
//! and the checker's reading of the copy it returns, so that reading is not
//! taken again.

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{BSplineCurve, Curve, Point3};

use super::geometric::{self, Way};
use super::{EdgeCopies, FaceCopyRef, GapMeasure, Hints, Limits, lift, measure, work_cap};

/// Source: the span limit of the fitter of `numerical_embedding.rs`
/// (`MAX_SPANS`), which builds the same kind of copy.
const MAX_SPANS: usize = 256;

/// Source: judgement. The rate of a foot comes from a forward or backward
/// difference over a thousandth of its span: small against the span, so the
/// difference error (second order in the step) is far under a fit target,
/// and large against the parameter's rounding.
const RATE_STEP_SHARE: f64 = 1e-3;

/// What a refit found for one edge use.
#[derive(Clone, Debug)]
pub struct Refit {
    /// The new face copy, a u/v curve in the edge's parameter.
    pub curve: BSplineCurve,
    /// The tolerance rule T4 stores for it, as the checker read it.
    pub stored_mm: f64,
    /// Spans of the new copy.
    pub spans: usize,
    /// The checker's reading of the new copy, which `stored_mm` comes from.
    pub reading: GapMeasure,
}

/// What a fit under an allowance found for one edge use, and what it took.
#[derive(Clone, Debug)]
pub struct Fit {
    /// The new face copy. `None` for a copy that has no u/v of its own, and
    /// for a change of parameter where the old copy has none that is valid.
    pub copy: Option<Refit>,
    /// Work counted against the allowance the fit was given: one unit per
    /// evaluation while a foot is found, one per reading of the old copy in
    /// a change of parameter, and the checker's work on every candidate,
    /// kept or not.
    pub work: u64,
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| {
        Refusal::new(
            "numerical.invalid-trim",
            "a refitted face copy is not a finite number",
        )
    })
}

/// Where the nodes of a new face copy are read.
#[derive(Clone, Copy, Debug)]
enum Nodes {
    /// The foot of the 3D point on the surface. The limits are what the edge
    /// use has left of its one work allowance, in an import admission.
    SurfaceFoot(Option<Limits>),
    /// The old copy, where its lift is nearest the 3D point. The limits are
    /// what the edge use has left of its one work allowance.
    OldCopy(Limits),
}

/// The work one fit has counted, and what it may still count.
///
/// Section 4.2 step 8 gives an edge one hard work cap. In an import
/// admission a fit is a later reading of its edge use: a change of parameter
/// (decision Q9) comes after the use's same-parameter gap and its geometric
/// gap, and rule T3's refit after those and after any change of parameter.
/// So the fit starts from what they left. Every foot, every reading of the
/// old copy and the checker's work on every candidate is taken from that
/// one remainder, and where it runs out the use refuses by name. Nothing
/// restarts between candidates.
///
/// The refit of an own build (rules T7 and T8) is given no remainder: its
/// feet are charged to the operation, and each candidate is another copy
/// that the checker measures under its own cap. It only counts.
struct Drawn {
    /// Units counted so far.
    spent: u64,
    /// The most the fit may count: the remainder of the use's allowance,
    /// where it was given one.
    limit: Option<u64>,
}

impl Drawn {
    fn new(nodes_from: Nodes) -> Self {
        Self {
            spent: 0,
            limit: match nodes_from {
                Nodes::SurfaceFoot(limits) => limits.map(|limits| limits.work),
                Nodes::OldCopy(limits) => Some(limits.work),
            },
        }
    }

    /// What is left to count.
    fn left(&self) -> u64 {
        self.limit
            .map_or(u64::MAX, |limit| limit.saturating_sub(self.spent))
    }

    /// The limits the checker measures the next candidate under: what is
    /// left of the use's allowance, or the checker's own cap for the refit
    /// of an own build.
    fn candidate(&self) -> Limits {
        match self.limit {
            Some(_) => Limits { work: self.left() },
            None => Limits::default(),
        }
    }

    /// Runs `work`, which takes each unit from what is left before it spends
    /// it and gives `None` where that runs out. What it took is counted and
    /// charged to the operation whether or not it finished. Where the
    /// allowance ran out the named refusal is the answer, also when the
    /// operation's own allowance is gone as well.
    fn draw<T>(
        &mut self,
        dec: &mut Decider,
        work: impl FnOnce(&mut u64) -> Option<T>,
    ) -> Result<T, Refusal> {
        let before = self.left();
        let mut left = before;
        let found = work(&mut left);
        let taken = before - left;
        self.spent += taken;
        let charged = if taken == 0 {
            Ok(())
        } else {
            dec.charge_work(taken)
        };
        let found = found.ok_or_else(|| {
            work_cap(&format!(
                "an edge spent its work allowance before its new face copy was fitted; \
                 {} units were spent on the fit",
                self.spent
            ))
        })?;
        charged?;
        Ok(found)
    }
}

/// Refits the face copy of one edge use until rule T4's stored tolerance is
/// at most `target_mm`, and returns the best candidate the checker read.
/// `None` for a copy that has no u/v of its own (a curve projected on a
/// plane: its gap is the curve's distance from the plane, which no copy
/// changes).
pub fn refit(
    copies: &EdgeCopies<'_>,
    target_mm: f64,
    dec: &mut Decider,
) -> Result<Option<Refit>, Refusal> {
    fit(copies, target_mm, Nodes::SurfaceFoot(None), dec).map(|(copy, _)| copy)
}

/// [`refit`] in an import admission (rule T3: "above that, refit first").
///
/// `limits` is what the edge use has left of its one work allowance after
/// the readings that came first. One unit is counted per evaluation while a
/// foot is found, and the checker counts its own units on each candidate.
/// All of them are taken from `limits` and charged to the operation, as
/// [`refit`] charges them. Where `limits` runs out the use refuses with
/// `numerical.gap-work-cap`. The work taken is returned with or without a
/// new copy.
pub fn refit_within(
    copies: &EdgeCopies<'_>,
    target_mm: f64,
    limits: Limits,
    dec: &mut Decider,
) -> Result<Fit, Refusal> {
    fit(copies, target_mm, Nodes::SurfaceFoot(Some(limits)), dec)
        .map(|(copy, work)| Fit { copy, work })
}

/// Gives the face copy of one edge use the parameter of its 3D curve
/// (decision Q9): the new copy is the old one read where its lift is nearest
/// the 3D point of each parameter, with its two ends kept. The spans double
/// until rule T4's stored tolerance is at most `target_mm`; the best
/// candidate the checker read is returned.
///
/// `limits` is what the edge use has left of its one work allowance after
/// the readings that came first. One unit is counted per evaluation of
/// either copy while a foot is found and per reading of the old copy, and
/// the checker counts its own units on each candidate. All of them are taken
/// from `limits` and charged to the operation. Where `limits` runs out the
/// use refuses with `numerical.gap-work-cap`; no candidate read before that
/// is returned in its place.
///
/// The copy comes with the checker's reading of it, the last one the fit
/// took, so an admission of the copy need not read it again.
///
/// No copy is returned for one that has no u/v of its own, and where the old
/// copy has no valid change of parameter: the nearest points must ascend
/// along it as the 3D curve is run through, which is the order OpenCASCADE's
/// `Approx_SameParameter` asks of its projected pairs. The work taken is
/// returned either way.
pub fn reparametrise(
    copies: &EdgeCopies<'_>,
    target_mm: f64,
    limits: Limits,
    dec: &mut Decider,
) -> Result<Fit, Refusal> {
    fit(copies, target_mm, Nodes::OldCopy(limits), dec).map(|(copy, work)| Fit { copy, work })
}

/// Source: Dennis and Schnabel 1983, section 5.4: a forward difference takes
/// a step of the square root of the machine epsilon. Two breaks closer than
/// that share of the edge leave no room for the rate step of a span between
/// them, so they are one node.
fn same_node(a: f64, b: f64, width: f64) -> bool {
    (a - b).abs() <= f64::EPSILON.sqrt() * width
}

/// The new copy, where one was found, and the work the fit counted: the
/// evaluations of its feet and the checker's work on every candidate.
fn fit(
    copies: &EdgeCopies<'_>,
    target_mm: f64,
    nodes_from: Nodes,
    dec: &mut Decider,
) -> Result<(Option<Refit>, u64), Refusal> {
    if matches!(copies.face_copy, FaceCopyRef::Projected) {
        return Ok((None, 0));
    }
    let lifted = lift::lift(copies, dec)?;
    let mut scratch = lift::Scratch::default();
    let mut drawn = Drawn::new(nodes_from);
    let [lo, hi] = lifted.range();
    let mut breaks = vec![lo];
    if let Curve::BSpline(curve) = copies.curve {
        for knot in curve.knots() {
            let k = knot.to_f64_nearest();
            if k > lo && k < hi && breaks.last() != Some(&k) {
                breaks.push(k);
            }
        }
    }
    breaks.push(hi);
    if let (Nodes::OldCopy(_), FaceCopyRef::Native(old)) = (nodes_from, copies.face_copy) {
        // The old copy is smooth between its knots only, and so is the new
        // one between the places where those knots land on the 3D curve.
        let mut last = None;
        for knot in old.knots() {
            let k = knot.to_f64_nearest();
            if k <= lo || k >= hi || last == Some(k) {
                continue;
            }
            last = Some(k);
            let landed = drawn.draw(dec, |left| {
                geometric::foot(&lifted, &mut scratch, Way::Backward, k, left)
            })?;
            if !breaks.iter().any(|b| same_node(*b, landed.foot, hi - lo)) {
                breaks.push(landed.foot);
            }
        }
        breaks.sort_by(f64::total_cmp);
    }
    let mut best: Option<Refit> = None;
    let mut cuts = 1;
    loop {
        let spans = (breaks.len() - 1) * cuts;
        if spans > MAX_SPANS {
            break;
        }
        let mut nodes = Vec::with_capacity(spans + 1);
        for window in breaks.windows(2) {
            for k in 0..cuts {
                nodes.push(window[0] + (window[1] - window[0]) * k as f64 / cuts as f64);
            }
        }
        nodes.push(hi);
        let mut knots = Vec::with_capacity(3 * spans + 5);
        let mut controls: Vec<[f64; 2]> = Vec::with_capacity(3 * spans + 1);
        let mut ascending = true;
        for (index, window) in nodes.windows(2).enumerate() {
            let [a, b] = [window[0], window[1]];
            let step = (b - a) * RATE_STEP_SHARE;
            // The u/v of the new copy at `t`, and the old copy's parameter
            // it was read at.
            let mut read = |t: f64, dec: &mut Decider| -> Result<([f64; 2], f64), Refusal> {
                let old = |s: f64, scratch: &mut lift::Scratch| {
                    lifted.face_at(s, scratch).map_or([0.0, 0.0], |uv| uv.value)
                };
                match nodes_from {
                    Nodes::SurfaceFoot(_) => {
                        let seed = old(t, &mut scratch);
                        let uv = drawn
                            .draw(dec, |left| lifted.foot_within(t, seed, &mut scratch, left))?;
                        Ok((uv, t))
                    }
                    // The two ends of the old copy stay the ends of the new:
                    // one reading of the old copy.
                    Nodes::OldCopy(_) if t <= lo || t >= hi => {
                        drawn.draw(dec, lift::spend)?;
                        Ok((old(t, &mut scratch), t))
                    }
                    // The foot, then one reading of the old copy there.
                    Nodes::OldCopy(_) => {
                        let found = drawn.draw(dec, |left| {
                            let found =
                                geometric::foot(&lifted, &mut scratch, Way::Forward, t, left)?;
                            lift::spend(left)?;
                            Some(found)
                        })?;
                        Ok((old(found.foot, &mut scratch), found.foot))
                    }
                }
            };
            let ((pa, sa), (pa2, sa2)) = (read(a, dec)?, read(a + step, dec)?);
            let ((pb, sb), (pb2, sb2)) = (read(b, dec)?, read(b - step, dec)?);
            ascending &= sa < sb && sa <= sa2 && sa2 <= sb2 && sb2 <= sb;
            let third = (b - a) / 3.0;
            let rate_a = [(pa2[0] - pa[0]) / step, (pa2[1] - pa[1]) / step];
            let rate_b = [(pb[0] - pb2[0]) / step, (pb[1] - pb2[1]) / step];
            if index == 0 {
                knots.extend([a; 4]);
                controls.push(pa);
            } else {
                knots.extend([a; 3]);
            }
            controls.push([pa[0] + rate_a[0] * third, pa[1] + rate_a[1] * third]);
            controls.push([pb[0] - rate_b[0] * third, pb[1] - rate_b[1] * third]);
            controls.push(pb);
        }
        if !ascending {
            // The old copy is not run through once in order: this is no
            // change of parameter.
            break;
        }
        knots.extend([hi; 4]);
        let points = controls
            .iter()
            .map(|uv| {
                Ok(Point3::new(
                    Real::from(&exact(uv[0])? + &lifted.shift[0]),
                    Real::from(&exact(uv[1])? + &lifted.shift[1]),
                    Real::zero(),
                ))
            })
            .collect::<Result<Vec<_>, Refusal>>()?;
        // The end knots are the edge's exact parameters, so the new copy
        // covers the whole finite edge; the interior knots lie strictly
        // inside and keep their binary64 values.
        let (first, last) = {
            let [a, b] = copies.parameters;
            (a.clone().min(b.clone()), a.clone().max(b.clone()))
        };
        let count = knots.len();
        let knots = knots
            .into_iter()
            .enumerate()
            .map(|(i, knot)| match i {
                0..=3 => Ok(first.clone()),
                i if i + 4 >= count => Ok(last.clone()),
                _ => exact(knot),
            })
            .collect::<Result<Vec<_>, Refusal>>()?;
        let weights = vec![Real::one(); points.len()];
        let curve = BSplineCurve::new(3, knots, points, weights, dec)?;
        let found = measure(
            &EdgeCopies {
                face_copy: FaceCopyRef::Native(&curve),
                ..*copies
            },
            &Hints::default(),
            drawn.candidate(),
            dec,
        )?;
        drawn.spent += found.work;
        // A copy that leaves a spline surface's domain is refused by
        // admission, so it is no candidate.
        let stored_mm = if found.domain_exit.is_some() {
            f64::INFINITY
        } else {
            found.stored_mm()
        };
        let better = stored_mm.is_finite() && best.as_ref().is_none_or(|b| stored_mm < b.stored_mm);
        if better {
            best = Some(Refit {
                curve,
                stored_mm,
                spans,
                reading: found,
            });
        }
        if stored_mm <= target_mm || !better {
            break;
        }
        cuts *= 2;
    }
    Ok((best, drawn.spent))
}
