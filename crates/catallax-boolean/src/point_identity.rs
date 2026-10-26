//! Canonical point witnesses with request-local certified coordinate reuse.

use crate::geometry::{coincident, record_same_point};
use catallax_numbers::{
    BoxRoots, Decider, FactKind, Outcome, PolySystem, Rational, Real, Refusal, RootBox,
    current_operation_work,
};
use catallax_shapes::{
    Curve, GeomKey, GeomStore, MeetOf, Point3, PointDef, PointId, Position, Surface, SurfaceId,
};
use catallax_topology::Model;
use std::{
    collections::{HashMap, hash_map::Entry},
    sync::Arc,
};

type Coordinates = [(Rational, Rational); 3];

enum Boxes {
    Separated,
    Overlapping,
    Unavailable,
}

/// One canonicalization scope on a single append-only model. Never reuse this
/// pool with a replacement model: point handles are local to its geometry store.
#[derive(Default)]
pub(crate) struct CanonicalPoints {
    pub(crate) points: Vec<PointId>,
    bounds: HashMap<PointId, Option<Coordinates>>,
}

fn closed_position(geom: &GeomStore, point: PointId) -> Option<&Point3> {
    match geom.point(point) {
        PointDef::Explicit(position) => Some(position),
        PointDef::Meet(meeting) => match &meeting.position {
            Position::Closed(position) => Some(position),
            Position::Isolated(_) => None,
        },
    }
}

fn coordinate(value: &Real, dec: &mut Decider) -> Result<(Rational, Rational), Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok((value.clone(), value.clone()));
    }
    let bound = dec.enclose_at(value, 64)?;
    Ok((bound.lo(), bound.hi()))
}

impl CanonicalPoints {
    fn cache(
        &mut self,
        geom: &GeomStore,
        point: PointId,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if let Entry::Vacant(entry) = self.bounds.entry(point) {
            let position = geom.position(point);
            let answer: Result<Coordinates, Refusal> = (|| {
                Ok([
                    coordinate(&position.x, dec)?,
                    coordinate(&position.y, dec)?,
                    coordinate(&position.z, dec)?,
                ])
            })();
            let bounds = match answer {
                Ok(bounds) => Some(bounds),
                Err(refusal)
                    if refusal.code == "numbers.deadline"
                        || refusal.code == "numbers.operation-budget"
                        || current_operation_work().exhausted =>
                {
                    return Err(refusal);
                }
                // This optional filter has no answer; the original scalar
                // comparison below still owns any arithmetic refusal.
                Err(_) => None,
            };
            entry.insert(bounds);
        }
        Ok(())
    }

    fn boxes(
        &mut self,
        geom: &GeomStore,
        a: PointId,
        b: PointId,
        dec: &mut Decider,
    ) -> Result<Boxes, Refusal> {
        if dec.policy().max_precision_bits < 64 {
            return Ok(Boxes::Unavailable);
        }
        let ready = matches!(self.bounds.get(&a), Some(Some(_)))
            && matches!(self.bounds.get(&b), Some(Some(_)));
        if !ready {
            if let Some(point) = closed_position(geom, b)
                && [&point.x, &point.y, &point.z]
                    .iter()
                    .all(|x| x.rational_ref().is_some())
            {
                return Ok(Boxes::Unavailable);
            }
            if let (Some(a), Some(b)) = (closed_position(geom, a), closed_position(geom, b))
                && [&a.x, &a.y, &a.z]
                    .into_iter()
                    .zip([&b.x, &b.y, &b.z])
                    .any(|(a, b)| {
                        let (a, b) = (a.approx(), b.approx());
                        a.hi() < b.lo() || b.hi() < a.lo()
                    })
            {
                // Keep the existing cheap scalar filter; do not pay for new
                // enclosures when its cached double hint already separates.
                return Ok(Boxes::Unavailable);
            }
            self.cache(geom, b, dec)?;
            self.cache(geom, a, dec)?;
        }
        let (Some(Some(a)), Some(Some(b))) = (self.bounds.get(&a), self.bounds.get(&b)) else {
            return Ok(Boxes::Unavailable);
        };
        // Cached arithmetic does not make a quadratic scan unmetered.
        dec.check_deadline()?;
        Ok(
            if a.iter()
                .zip(b)
                .any(|((alo, ahi), (blo, bhi))| ahi < blo || bhi < alo)
            {
                Boxes::Separated
            } else {
                Boxes::Overlapping
            },
        )
    }

    pub(crate) fn insert(
        &mut self,
        model: &mut Model,
        point: PointId,
        dec: &mut Decider,
        code: &'static str,
    ) -> Result<PointId, Refusal> {
        let span = catallax_numbers::OperationSpan::begin(
            "canonical-point",
            &format!("{point:?}"),
            "intern",
            dec,
        );
        let result = (|| {
            let mut shared = SharedCylinderPoint::new(model.geom(), point);
            for index in 0..self.points.len() {
                let existing = self.points[index];
                let equal = if model.geom().facts().holds(
                    FactKind::Coincident,
                    GeomKey::Point(existing),
                    GeomKey::Point(point),
                ) {
                    dec.note_fact_decision();
                    true
                } else {
                    match self.boxes(model.geom(), existing, point, dec)? {
                        Boxes::Separated => false,
                        Boxes::Unavailable => {
                            coincident(model.geom().same_point(existing, point, dec), code)?
                        }
                        Boxes::Overlapping => {
                            let answer = match shared.as_mut() {
                                Some(shared) => shared.compare(
                                    existing,
                                    self.bounds[&existing].as_ref().unwrap(),
                                    self.bounds[&point].as_ref().unwrap(),
                                    dec,
                                )?,
                                None => None,
                            };
                            match answer {
                                Some(equal) => equal,
                                None => {
                                    coincident(model.geom().same_point(existing, point, dec), code)?
                                }
                            }
                        }
                    }
                };
                if equal {
                    drop(shared);
                    record_same_point(model, existing, point);
                    return Ok(existing);
                }
            }
            self.points.push(point);
            Ok(point)
        })();
        if result.is_err() {
            let source = match model.geom().point(point) {
                PointDef::Meet(meeting) => format!("{:?}", meeting.of),
                PointDef::Explicit(_) => "Explicit".to_owned(),
            };
            let failed = catallax_numbers::OperationSpan::begin(
                "canonical-source",
                &format!("{point:?}/{source}"),
                "refused",
                dec,
            );
            failed.finish(dec, None, None);
        }
        span.finish(dec, result.as_ref().err(), None);
        result
    }
}

struct SharedCylinderPoint<'a> {
    geom: &'a GeomStore,
    point: PointId,
    surfaces: Option<[SurfaceId; 3]>,
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

impl<'a> SharedCylinderPoint<'a> {
    fn new(geom: &'a GeomStore, point: PointId) -> Option<Self> {
        let PointDef::Meet(meeting) = geom.point(point) else {
            return None;
        };
        let MeetOf::CurveSurface(_, surface) = &meeting.of else {
            return None;
        };
        if !matches!(geom.surface(*surface), Surface::Cylinder(_))
            || !geom.surface(*surface).has_rational_definition()
        {
            return None;
        }
        Some(Self {
            geom,
            point,
            surfaces: cylinder_sources(geom, point),
        })
    }

    fn compare(
        &mut self,
        existing: PointId,
        a: &Coordinates,
        b: &Coordinates,
        dec: &mut Decider,
    ) -> Result<Option<bool>, Refusal> {
        let Some(surfaces) = self
            .surfaces
            .or_else(|| cylinder_sources(self.geom, existing))
        else {
            return Ok(None);
        };
        if surfaces.iter().any(|&surface| {
            [existing, self.point].into_iter().any(|point| {
                !self.geom.facts().holds(
                    FactKind::Incident,
                    GeomKey::Point(point),
                    GeomKey::Surface(surface),
                )
            })
        }) {
            return Ok(None);
        }
        let span = catallax_numbers::OperationSpan::begin(
            "point-identity",
            &format!("{existing:?}/{:?}", self.point),
            "compare",
            dec,
        );
        let result = (|| {
            let Some(equations) = surfaces
                .iter()
                .map(|&surface| self.geom.implicit(surface)?.as_rational())
                .collect::<Option<Vec<_>>>()
            else {
                return Ok(None);
            };
            let margin = Rational::one().mul_pow2(-32);
            let lo = a
                .iter()
                .zip(b)
                .map(|((alo, _), (blo, _))| alo.min(blo) - &margin)
                .collect();
            let hi = a
                .iter()
                .zip(b)
                .map(|((_, ahi), (_, bhi))| ahi.max(bhi) + &margin)
                .collect();
            let system =
                Arc::new(PolySystem::new(equations).expect("three cylinder equations in xyz"));
            let candidate = RootBox::new(lo, hi).expect("ordered coordinate hull");
            // Both points are proved roots of this system and lie in this hull.
            // Padding supplies an interior for the certificate, not a coincidence
            // tolerance: two distinct roots cannot pass the uniqueness proof.
            match dec.certify_root(&system, &candidate) {
                Outcome::Proven(BoxRoots::One(_)) => Ok(Some(true)),
                Outcome::Refused(refusal)
                    if refusal.code == "numbers.root-not-certified"
                        && !current_operation_work().exhausted =>
                {
                    Ok(None)
                }
                Outcome::Refused(refusal) => Err(refusal),
                _ => Ok(None),
            }
        })();
        span.finish(dec, result.as_ref().err(), None);
        let outcome = match &result {
            Ok(Some(true)) => "equal",
            Ok(Some(false)) => "apart",
            Ok(None) => "fallback",
            Err(_) => "refused",
        };
        let state =
            catallax_numbers::OperationSpan::begin("point-identity-result", &outcome, "state", dec);
        state.finish(dec, None, None);
        result
    }
}

#[cfg(test)]
mod tests;
