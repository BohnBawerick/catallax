//! Fitted carriers keep their exact spline coefficients and recorded fit bound.

use super::*;
use catallax_shapes::RollingEnvelope;

fn proposals(
    cutter: &SectionSurface,
    source: &BSplineSurface,
    swap: bool,
    dec: &mut Decider,
) -> Result<Candidates, Refusal> {
    if let SectionSurface::Plane(plane) = cutter {
        return Candidates::plane_section(plane, source, swap, dec);
    }
    let (native, map) = match cutter {
        SectionSurface::Cylinder(cylinder) => {
            let direction = axis(&cylinder.axis, dec)?;
            let origin = point(&cylinder.origin, dec)?;
            let radius = rational(&cylinder.radius(), dec)?;
            let real_origin = origin.to_real();
            let real_direction = direction.to_real();
            let mut bound = 0.0_f64;
            for point in source.controls() {
                dec.check_deadline()?;
                let height = point.minus(&real_origin).dot(&real_direction).to_f64();
                if !height.is_finite() {
                    return Err(Refusal::new(
                        "intersect.trace-fitted-proposals",
                        "the unproved cylinder proposal window needs finite controls",
                    ));
                }
                bound = bound.max(height.abs());
            }
            // This approximate window belongs only to the proposer. The exact
            // section below cuts by the unrestricted original native equation.
            let bound = Rational::from_f64_exact(bound + 1.0).ok_or_else(|| {
                Refusal::new(
                    "intersect.trace-fitted-proposals",
                    "the unproved cylinder proposal window overflowed",
                )
            })?;
            (
                cylinder_source(&radius, &bound, dec)?,
                frame(&direction, &origin)?,
            )
        }
        SectionSurface::Torus(torus) => {
            let major = rational(&torus.major, dec)?;
            let minor = rational(&torus.minor, dec)?;
            if minor.sign() != Sign::Positive || major <= minor {
                return Err(Refusal::new(
                    "intersect.trace-fitted-proposals",
                    "the torus proposer needs a strict ring torus",
                ));
            }
            let direction = axis(&torus.axis, dec)?;
            let center = point(&torus.center, dec)?;
            (
                torus_source(&major, &minor, dec)?,
                frame(&direction, &center)?,
            )
        }
        _ => {
            return Err(Refusal::new(
                "intersect.trace-fitted-proposals",
                "this fitted section needs an analytic proposal-surface adapter",
            ));
        }
    };
    let identity = Isometry::identity();
    Candidates::generate(
        if swap {
            [(source, &identity), (&native, &map)]
        } else {
            [(&native, &map), (source, &identity)]
        },
        dec,
    )
}

pub(in crate::pair) fn section(
    cutter: &SectionSurface,
    envelope: &RollingEnvelope,
    swap: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let mut stage = "unproved analytic/spline branch proposals";
    let result = (|| -> Result<Plan, Refusal> {
        let source = envelope.carrier();
        let candidates = proposals(cutter, source, swap, dec)?;
        stage = "exact fitted-source embedding and complete parameter-chart cover";
        // Keep the original algebraic coefficients. Do not replace controls by
        // midpoints or attach the envelope's error bound to a different spline.
        let mut atlases = RationalAtlas::partition(source.clone(), [false, false], dec)?;
        if atlases.len() != 1 {
            return Err(Refusal::new(
                "intersect.trace-fitted-component-ownership",
                "the fitted source needs a single complete joined atlas before native section ownership",
            ));
        }
        let atlas = atlases.pop().expect("one complete atlas was checked");
        if !atlas.source().same_definition(source) {
            return Err(Refusal::new(
                "intersect.trace-fitted-source-binding",
                "the admitted atlas needs the exact retained rolling-envelope carrier",
            ));
        }
        stage = "complete exact analytic-section certificates and endpoint topology";
        let section = atlas.sections(cutter, dec)?;
        stage = "exact proposal selection of every certified fitted-section component";
        candidates.prove_sections(section.branches(), dec)?;
        if section.branches().is_empty() {
            return Ok(Plan::disjoint());
        }
        stage = "oriented retained fitted-section branches";
        let mut plan = Plan::empty();
        super::super::rational_patch::push_branches(&mut plan, section.branches(), swap, dec)?;
        Ok(plan)
    })();
    result.map_err(|error| {
        Refusal::new(
            error.code,
            format!("fitted rolling-envelope section {stage}: {}", error.detail),
        )
    })
}
