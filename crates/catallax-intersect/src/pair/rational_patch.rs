//! Complete finite analytic sections from the checked rational chart definition.

use super::{Contact, CurveForm, Plan, Spec};
use crate::math::sign;
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{
    BSplineCurve, PatchSectionBranch, RationalAtlas, RationalPatch, SectionSurface,
};

fn reversed(curve: &BSplineCurve, dec: &mut Decider) -> Result<BSplineCurve, Refusal> {
    let mut controls = curve.controls().cloned().collect::<Vec<_>>();
    let mut weights = curve.weights().cloned().collect::<Vec<_>>();
    controls.reverse();
    weights.reverse();
    let [lo, hi] = curve.domain();
    let sum = lo + hi;
    let knots = curve.knots().iter().rev().map(|knot| &sum - knot).collect();
    BSplineCurve::new(curve.degree(), knots, controls, weights, dec)
}

/// Sections by a plane, cylinder or sphere. With `swap` false the analytic
/// surface is the pair's first surface.
pub(super) fn section(
    surface: &SectionSurface,
    patch: &RationalPatch,
    swap: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let section = patch.sections(surface, dec)?;
    let carrier = surface.surface();
    let mut plan = Plan::empty();
    plan.points.extend_from_slice(section.points());
    push_branches(&mut plan, section.branches(), swap, dec)?;
    for (curve, crossing) in section.curves() {
        let [lo, hi] = curve.domain();
        let middle = Real::from_rational((lo + hi).mul_pow2(-1));
        let reverse = if crossing {
            let point = curve.point_at(&middle, dec)?;
            let normal = patch.normal_at(&point, dec)?;
            let analytic = carrier.normal_at(&point, dec)?;
            let along = if swap {
                normal.cross(&analytic)
            } else {
                analytic.cross(&normal)
            };
            let tangent = curve.first_at(&middle, dec)?;
            match sign(&tangent.dot(&along), dec)? {
                Sign::Positive => false,
                Sign::Negative => true,
                Sign::Zero => {
                    return Err(Refusal::new(
                        "intersect.rational-patch-section",
                        "a higher-order crossing has no certified pair-normal orientation",
                    ));
                }
            }
        } else {
            swap
        };
        let curve = if reverse {
            reversed(curve, dec)?
        } else {
            curve.clone()
        };
        plan.curves.push(Spec {
            form: CurveForm::BSpline(curve),
            witness: None,
            contact: if crossing {
                Contact::Crossing
            } else {
                Contact::Tangent
            },
            singular: Vec::new(),
        });
    }
    if plan.curves.is_empty() && plan.points.is_empty() {
        Ok(Plan::disjoint())
    } else {
        Ok(plan)
    }
}

/// Oriented exact branches: the first surface's normal cross the second's.
fn push_branches(
    plan: &mut Plan,
    branches: &[PatchSectionBranch],
    swap: bool,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for branch in branches {
        let branch = if branch.follows_surface_first(dec)? != swap {
            branch.clone()
        } else {
            branch.reversed()
        };
        plan.curves.push(Spec {
            form: CurveForm::PatchSection(branch),
            witness: None,
            contact: Contact::Crossing,
            singular: Vec::new(),
        });
    }
    Ok(())
}

/// Sections of a joined atlas face. Every component is a cell-cover branch
/// over the face's charts; there are no isolated contact points.
pub(super) fn atlas(
    surface: &SectionSurface,
    atlas: &RationalAtlas,
    swap: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let section = atlas.sections(surface, dec)?;
    if section.branches().is_empty() {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    push_branches(&mut plan, section.branches(), swap, dec)?;
    Ok(plan)
}
