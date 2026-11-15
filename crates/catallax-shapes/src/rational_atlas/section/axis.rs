//! Reuse exact univariate roots when every component has physical border ends.

use super::*;
use crate::rational_patch::{axis_section, exact_sign};
use catallax_numbers::{Real, Sign};

pub(super) fn section(
    atlas: &RationalAtlas,
    surface: &SectionSurface,
    dec: &mut Decider,
) -> Result<Option<Vec<(BSplineCurve, bool)>>, Refusal> {
    let mut components = Vec::new();
    for chart in atlas.charts() {
        let Some(section) = axis_section(chart.patch(), surface, dec)? else {
            return Ok(None);
        };
        if !section.points().is_empty() || !section.branches().is_empty() {
            return Ok(None);
        }
        let curves = section.curves().collect::<Vec<_>>();
        if curves.is_empty() {
            continue;
        }
        // shortcut: components spanning several charts keep the complete cell
        // search until their parameter intervals can be joined explicitly.
        if !components.is_empty() {
            return Ok(None);
        }
        for (curve, crossing) in curves {
            for bound in curve.domain() {
                let point = curve.point_at(&Real::from_rational(bound.clone()), dec)?;
                let Some(at) = chart.patch().parameters_at(&point, dec)? else {
                    return Ok(None);
                };
                let mut physical = false;
                for (side, axis, end) in [
                    (SplineBoundary::UMin, 0, 0),
                    (SplineBoundary::UMax, 0, 1),
                    (SplineBoundary::VMin, 1, 0),
                    (SplineBoundary::VMax, 1, 1),
                ] {
                    if chart.neighbor(side).is_none()
                        && exact_sign(
                            &(&at[axis]
                                - &Real::from_rational(chart.patch().domain()[axis][end].clone())),
                            dec,
                        )? == Sign::Zero
                    {
                        physical = true;
                        break;
                    }
                }
                if !physical {
                    return Ok(None);
                }
            }
            components.push((curve.clone(), crossing));
        }
    }
    Ok(Some(components))
}
