//! Complete finite plane sections from the checked rational chart definition.

use super::{Contact, CurveForm, Plan, Spec};
use crate::math::sign;
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, Plane, RationalPatch};

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

pub(super) fn plane(
    plane: &Plane,
    patch: &RationalPatch,
    swap: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let section = patch.plane_sections(plane, dec)?;
    let mut plan = Plan::empty();
    plan.points.extend_from_slice(section.points());
    for (curve, crossing) in section.curves() {
        let [lo, hi] = curve.domain();
        let middle = Real::from_rational((lo + hi).mul_pow2(-1));
        let reverse = if crossing {
            let point = curve.point_at(&middle, dec)?;
            let normal = patch.normal_at(&point, dec)?;
            let along = if swap {
                normal.cross(&plane.normal)
            } else {
                plane.normal.cross(&normal)
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
