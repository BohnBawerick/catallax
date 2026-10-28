//! Exact, export-only cuts of periodic carrier charts.

use catallax_intersect::{Contact, CurveSurface};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{PointId, Surface, Vector3};
use catallax_topology::{EdgeId, FaceId, Model};
use std::collections::BTreeSet;

use crate::{
    number,
    trim::{self, Geometry, Location},
};

pub(crate) struct Crossing {
    pub edge: EdgeId,
    pub point: PointId,
    pub height: Real,
    /// Signed distance from the torus core circle along the meridian
    /// reference; zero on other carriers.
    pub offset: Real,
    /// The normalized, surface-CCW boundary enters the meridian's upper side.
    pub enters: bool,
    pub increasing: bool,
}

pub(crate) struct Meridian {
    pub reference: Vector3,
    pub crossings: Vec<Crossing>,
}

/// A generic meridian avoids source vertices, contained boundary curves, and
/// tangent contacts. Consequently every accepted event has two distinct chart
/// sides and a certified entering/exiting direction.
pub(crate) fn meridian(
    geometry: &mut Geometry,
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Meridian, Refusal> {
    let mut major = None;
    let (origin, axis, reference) = match geometry.store.surface(model.face(face).surface) {
        Surface::Cylinder(s) => (s.origin.clone(), s.axis.clone(), s.ref_dir.clone()),
        Surface::Sphere(s) => (s.center.clone(), s.axis.clone(), s.ref_dir.clone()),
        Surface::Cone(s) => (s.apex.clone(), s.axis.clone(), s.ref_dir.clone()),
        Surface::Torus(s) => {
            major = Some(s.major.clone());
            (s.center.clone(), s.axis.clone(), s.ref_dir.clone())
        }
        _ => unreachable!("meridian on a non-periodic carrier"),
    };
    let axis = number::unit(&axis);
    let reference = number::unit(&reference);
    let across = axis.cross(&reference);
    let boundaries = geometry.boundaries(model, face, dec)?;
    let attempts = boundaries.len().saturating_add(1).saturating_mul(128);
    let mut fallback = None;
    for attempt in 0..attempts.saturating_add(4) {
        if attempt == 4
            && let Some(meridian) = fallback.take()
        {
            return Ok(meridian);
        }
        dec.check_deadline()?;
        // Prefer an exterior chart cut before splitting a contractible trim.
        // The antipodal ray shares the first plane/intersection cache and
        // keeps ordinary partial walls and caps as one face.
        let reference = match attempt {
            0 => reference.clone(),
            1 => reference.scaled(&(-Real::one())),
            2 => across.clone(),
            3 => across.scaled(&(-Real::one())),
            _ => {
                let serial = attempt - 1;
                let magnitude = serial.div_ceil(2) as i64;
                let t = Real::from_i64(if serial % 2 == 0 {
                    -magnitude
                } else {
                    magnitude
                });
                reference
                    .scaled(&(Real::one() - &t * &t))
                    .plus(&across.scaled(&(&Real::from_i64(2) * &t)))
            }
        };
        let normal = axis.cross(&reference);
        let plane = geometry.plane(origin.clone(), normal, axis.clone(), dec)?;
        let mut crossings = Vec::new();
        let mut seen = BTreeSet::new();
        let mut regular = true;
        for boundary in boundaries.iter() {
            let intersections = geometry.intersections(boundary.curve, plane, dec)?;
            match intersections.as_ref() {
                CurveSurface::Disjoint => {}
                CurveSurface::Contained => {
                    regular = false;
                    break;
                }
                CurveSurface::Meets(hits) => {
                    for hit in hits {
                        let location = geometry.on_arc(boundary, hit.point, dec)?;
                        if location == Location::Outside {
                            continue;
                        }
                        let radial = geometry.store.position(hit.point).minus(&origin);
                        if trim::sign(&radial.dot(&reference), dec)? != Sign::Positive {
                            continue;
                        }
                        if location == Location::Boundary || hit.contact != Contact::Crossing {
                            regular = false;
                            break;
                        }
                        if !seen.insert((boundary.edge, hit.point)) {
                            regular = false;
                            break;
                        }
                        let tangent = catallax_intersect::tangent_at(
                            &geometry.store,
                            boundary.curve,
                            hit.point,
                            dec,
                        )?;
                        let direction = trim::sign(&axis.dot(&radial.cross(&tangent[1])), dec)?;
                        if direction == Sign::Zero {
                            regular = false;
                            break;
                        }
                        let offset = major.as_ref().map_or_else(Real::zero, |major| {
                            &radial.dot(&number::unit(&reference)) - major
                        });
                        crossings.push(Crossing {
                            edge: boundary.edge,
                            point: hit.point,
                            height: radial.dot(&axis),
                            offset,
                            enters: (direction == Sign::Positive) == boundary.forward,
                            increasing: direction == Sign::Positive,
                        });
                    }
                }
            }
            if !regular {
                break;
            }
        }
        if !regular {
            continue;
        }
        let torus = major.is_some();
        trim::try_sort(&mut crossings, |a, b| order(a, b, torus, dec))?;
        for pair in crossings.windows(2) {
            if order(&pair[0], &pair[1], torus, dec)? == std::cmp::Ordering::Equal {
                regular = false;
                break;
            }
        }
        if regular {
            let meridian = Meridian {
                reference,
                crossings,
            };
            if meridian.crossings.is_empty() || attempt >= 4 {
                return Ok(meridian);
            }
            if fallback.is_none() {
                fallback = Some(meridian);
            }
        }
    }
    Err(number::budget(
        "no transverse periodic meridian certificate",
    ))
}

/// The written toroidal chart puts its v seam on the outer equator. Accept a
/// torus face only when that whole circle is outside the face: no boundary arc
/// meets it, and the meridian's crossings alternate from entering to exiting
/// in tube-angle order, so the meridian arc through the equator is exterior.
/// A disjoint connected circle with one exterior point is exterior throughout.
pub(crate) fn torus_seam_outside(
    geometry: &mut Geometry,
    model: &Model,
    face: FaceId,
    meridian: &Meridian,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let Surface::Torus(torus) = geometry.store.surface(model.face(face).surface).clone() else {
        unreachable!("torus seam check on another carrier")
    };
    if meridian.crossings.is_empty() || !meridian.crossings.len().is_multiple_of(2) {
        return Ok(false);
    }
    for (index, crossing) in meridian.crossings.iter().enumerate() {
        if crossing.enters != index.is_multiple_of(2)
            || (trim::sign(&crossing.height, dec)? == Sign::Zero
                && trim::sign(&crossing.offset, dec)? != Sign::Negative)
        {
            return Ok(false);
        }
    }
    let axis = number::unit(&torus.axis);
    let major_sq = torus.major.square();
    let outer = |geometry: &Geometry, point: PointId, dec: &mut Decider| {
        let radial = geometry.store.position(point).minus(&torus.center);
        let height = radial.dot(&axis);
        let distance = &(&radial.dot(&radial) - &(&height * &height)) - &major_sq;
        Ok::<_, Refusal>(trim::sign(&distance, dec)? != Sign::Negative)
    };
    let reference = number::unit(&torus.ref_dir);
    let equator = geometry.plane(torus.center.clone(), axis.clone(), reference, dec)?;
    for boundary in geometry.boundaries(model, face, dec)?.iter() {
        match geometry
            .intersections(boundary.curve, equator, dec)?
            .as_ref()
        {
            CurveSurface::Disjoint => {}
            CurveSurface::Contained => {
                // A contained trim curve is one of the two equator circles.
                let point = match boundary.from {
                    Some(point) => point,
                    None => {
                        let origin = geometry
                            .store
                            .curve(boundary.curve)
                            .parameter_origin(&geometry.store);
                        geometry.point(origin)
                    }
                };
                if outer(geometry, point, dec)? {
                    return Ok(false);
                }
            }
            CurveSurface::Meets(hits) => {
                for hit in hits {
                    if geometry.on_arc(boundary, hit.point, dec)? != Location::Outside
                        && outer(geometry, hit.point, dec)?
                    {
                        return Ok(false);
                    }
                }
            }
        }
    }
    Ok(true)
}

/// Order crossings along the meridian. A torus meridian is a tube circle, so
/// its crossings are ordered by the tube angle in `[0, 2pi)`, starting on the
/// outer equator; other carriers order by height along the axis.
fn order(
    a: &Crossing,
    b: &Crossing,
    torus: bool,
    dec: &mut Decider,
) -> Result<std::cmp::Ordering, Refusal> {
    let ordering = |sign| match sign {
        Sign::Negative => std::cmp::Ordering::Less,
        Sign::Zero => std::cmp::Ordering::Equal,
        Sign::Positive => std::cmp::Ordering::Greater,
    };
    if !torus {
        return Ok(ordering(trim::sign(&(&a.height - &b.height), dec)?));
    }
    let mut half = |crossing: &Crossing| -> Result<u8, Refusal> {
        Ok(match trim::sign(&crossing.height, dec)? {
            Sign::Positive => 1,
            Sign::Negative => 3,
            Sign::Zero if trim::sign(&crossing.offset, dec)? == Sign::Negative => 2,
            Sign::Zero => 0,
        })
    };
    let (first, second) = (half(a)?, half(b)?);
    if first != second || first % 2 == 0 {
        return Ok(first.cmp(&second));
    }
    // Within one open half circle the earlier angle turns positively to the later.
    let turn = &(&a.offset * &b.height) - &(&a.height * &b.offset);
    Ok(ordering(trim::sign(&turn, dec)?).reverse())
}

pub(crate) struct Plan {
    pub faces: std::collections::BTreeMap<FaceId, crate::wire::ChartFace>,
    pub surfaces: BTreeSet<catallax_shapes::SurfaceId>,
    pub curves: BTreeSet<catallax_shapes::CurveId>,
    pub uses: std::collections::BTreeMap<EdgeId, Vec<(FaceId, catallax_topology::CoedgeId)>>,
    pub written_vertices:
        std::collections::BTreeMap<catallax_topology::VertexId, crate::number::Point>,
    pub added_edges: usize,
    pub added_vertices: usize,
}

impl Plan {
    pub fn new(
        model: &Model,
        faces: std::collections::BTreeMap<FaceId, crate::wire::ChartFace>,
        added_edges: usize,
        added_vertices: usize,
    ) -> Self {
        let mut surfaces = BTreeSet::new();
        let mut curves = BTreeSet::new();
        let mut uses = std::collections::BTreeMap::<EdgeId, Vec<_>>::new();
        for &face in faces.keys() {
            surfaces.insert(model.face(face).surface);
            for loop_id in &model.face(face).loops {
                let catallax_topology::Loop::Coedges(coedges) = model.loop_(*loop_id) else {
                    unreachable!()
                };
                for &coedge in coedges {
                    let edge = model.coedge(coedge).edge;
                    curves.insert(model.edge(edge).curve);
                    uses.entry(edge).or_default().push((face, coedge));
                }
            }
        }
        Self {
            faces,
            surfaces,
            curves,
            uses,
            added_edges,
            added_vertices,
            written_vertices: Default::default(),
        }
    }
}
