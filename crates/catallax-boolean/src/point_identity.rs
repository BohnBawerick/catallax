//! Optional uniqueness certificates for existing cylinder-intersection witnesses.

use catallax_numbers::{BoxRoots, Decider, FactKind, Outcome, Rational, Refusal, RootBox};
use catallax_shapes::{
    Curve, GeomKey, GeomStore, MeetOf, Point3, PointDef, PointId, Position, Surface, SurfaceId,
};
use std::sync::Arc;

fn closed_position(geom: &GeomStore, point: PointId) -> Option<&Point3> {
    match geom.point(point) {
        PointDef::Explicit(position) => Some(position),
        PointDef::Meet(meeting) => match &meeting.position {
            Position::Closed(position) => Some(position),
            Position::Isolated(_) => None,
        },
    }
}

fn cylinder_sources(geom: &GeomStore, point: PointId) -> Option<[SurfaceId; 3]> {
    let PointDef::Meet(meeting) = geom.point(point) else {
        return None;
    };
    let MeetOf::CurveSurface(curve, surface) = &meeting.of else {
        return None;
    };
    let Curve::Intersection(intersection) = geom.curve(*curve) else {
        return None;
    };
    let surfaces = [intersection.surfaces[0], intersection.surfaces[1], *surface];
    if surfaces[0] == surfaces[1]
        || surfaces[0] == surfaces[2]
        || surfaces[1] == surfaces[2]
        || surfaces.iter().any(|&surface| {
            !matches!(geom.surface(surface), Surface::Cylinder(_))
                || !geom.surface(surface).has_rational_definition()
        })
    {
        None
    } else {
        Some(surfaces)
    }
}

/// Proves separation or equality; inconclusive bounds retain scalar equality.
/// Incidence establishes both roots, not which root, trim membership or orientation.
pub(crate) fn same_cylinder_point(
    geom: &GeomStore,
    existing: PointId,
    incoming: PointId,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    if dec.policy().max_precision_bits < 64 {
        return Ok(None);
    }
    let (Some(a), Some(b)) = (
        closed_position(geom, existing),
        closed_position(geom, incoming),
    ) else {
        return Ok(None);
    };
    let coordinates = [(&a.x, &b.x), (&a.y, &b.y), (&a.z, &b.z)];
    if coordinates
        .iter()
        .all(|(a, b)| a.rational_ref().is_some() && b.rational_ref().is_some())
        || coordinates.iter().any(|(a, b)| {
            let (a, b) = (a.approx(), b.approx());
            a.hi() < b.lo() || b.hi() < a.lo()
        })
    {
        return Ok(None);
    }
    let Some(surfaces) =
        cylinder_sources(geom, incoming).or_else(|| cylinder_sources(geom, existing))
    else {
        return Ok(None);
    };
    if surfaces.iter().any(|&surface| {
        [existing, incoming].into_iter().any(|point| {
            !geom.facts().holds(
                FactKind::Incident,
                GeomKey::Point(point),
                GeomKey::Surface(surface),
            )
        })
    }) {
        return Ok(None);
    }
    dec.check_deadline()?;
    let mut lo = std::array::from_fn::<_, 3, _>(|_| Rational::zero());
    let mut hi = std::array::from_fn::<_, 3, _>(|_| Rational::zero());
    let margin = Rational::one().mul_pow2(-32);
    for (axis, (a, b)) in coordinates.into_iter().enumerate() {
        let bounds = (|| {
            Ok::<_, Refusal>((
                dec.enclose_budgeted_at(a, 64)?,
                dec.enclose_budgeted_at(b, 64)?,
            ))
        })();
        let (a, b) = match bounds {
            Ok(bounds) => bounds,
            Err(error)
                if error.code == "numbers.deadline"
                    || error.code == "numbers.operation-budget"
                    || catallax_numbers::current_operation_work().exhausted =>
            {
                return Err(error);
            }
            Err(_) => return Ok(None),
        };
        if a.hi() < b.lo() || b.hi() < a.lo() {
            return Ok(Some(false));
        }
        lo[axis] = &a.lo().min(b.lo()) - &margin;
        hi[axis] = &a.hi().max(b.hi()) + &margin;
    }
    let Some(equations) = surfaces
        .iter()
        .map(|&surface| geom.implicit(surface)?.as_rational())
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let system = Arc::new(
        catallax_numbers::PolySystem::new(equations).expect("three cylinder equations in xyz"),
    );
    let candidate = RootBox::new(Vec::from(lo), Vec::from(hi)).expect("ordered coordinate hull");
    // Both retained witnesses are roots in this hull. Padding gives the
    // certificate an interior; it is not a coincidence tolerance. A hull with
    // two roots cannot prove uniqueness. One precision attempt leaves the
    // original scalar method available near singularities.
    match dec.certify_root_at(&system, &candidate, 64) {
        Outcome::Proven(BoxRoots::One(_)) => Ok(Some(true)),
        Outcome::Refused(error)
            if error.code == "numbers.root-not-certified"
                && !catallax_numbers::current_operation_work().exhausted =>
        {
            Ok(None)
        }
        Outcome::Refused(error) => Err(error),
        _ => Ok(None),
    }
}

#[cfg(test)]
mod tests;
