//! Finite control hulls exclude native spline/patch meetings.

use super::*;
use catallax_numbers::Ball;
use catallax_shapes::SplineSupport;

mod bilinear;

pub(super) fn meetings(
    geom: &mut GeomStore,
    curve: CurveId,
    surface: SurfaceId,
    dec: &mut Decider,
) -> Result<Option<CurveSurface>, Refusal> {
    bilinear::meetings(geom, curve, surface, dec)
}

fn hull<'a>(
    controls: impl Iterator<Item = &'a Point3>,
    bits: u32,
    dec: &mut Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    let mut result: Option<[Ball; 3]> = None;
    for point in controls {
        dec.charge_work(1)?;
        let next = [&point.x, &point.y, &point.z]
            .map(|value| dec.enclose_budgeted_at(value, bits))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        let next: [Ball; 3] = next.try_into().unwrap();
        result = Some(match result {
            None => next,
            Some(previous) => std::array::from_fn(|axis| previous[axis].hull(&next[axis])),
        });
    }
    Ok(result)
}

pub(super) fn disjoint(
    geom: &GeomStore,
    curve: CurveId,
    surface: SurfaceId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (Curve::BSpline(curve), Surface::BSpline(support)) =
        (geom.curve(curve), geom.surface(surface))
    else {
        return Ok(false);
    };
    let source = match support {
        SplineSupport::RationalPatch(patch) => patch.carrier(),
        SplineSupport::RationalAtlas(atlas) => atlas.source(),
        _ => return Ok(false),
    };
    let bits = dec.policy().max_precision_bits.min(128);
    if bits == 0 {
        return Ok(false);
    }
    // Admitted positive weights keep each complete finite source in its
    // control hull. Strictly separated enclosures therefore exclude meetings.
    let Some(first) = hull(curve.controls(), bits, dec)? else {
        return Ok(false);
    };
    let Some(second) = hull(source.controls().iter(), bits, dec)? else {
        return Ok(false);
    };
    for axis in 0..3 {
        dec.charge_work(1)?;
        if first[axis].hi() < second[axis].lo() || second[axis].hi() < first[axis].lo() {
            return Ok(true);
        }
    }
    Ok(false)
}
