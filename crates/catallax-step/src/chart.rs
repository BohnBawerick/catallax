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
    let (origin, axis, reference) = match geometry.store.surface(model.face(face).surface) {
        Surface::Cylinder(s) => (s.origin.clone(), s.axis.clone(), s.ref_dir.clone()),
        Surface::Sphere(s) => (s.center.clone(), s.axis.clone(), s.ref_dir.clone()),
        Surface::Cone(s) => (s.apex.clone(), s.axis.clone(), s.ref_dir.clone()),
        _ => unreachable!("meridian on a non-periodic carrier"),
    };
    let axis = number::unit(&axis);
    let reference = number::unit(&reference);
    let across = axis.cross(&reference);
    let boundaries = geometry.boundaries(model, face, dec)?;
    let attempts = boundaries.len().saturating_add(1).saturating_mul(128);
    for attempt in 0..attempts {
        dec.check_deadline()?;
        let magnitude = attempt.div_ceil(2) as i64;
        let t = Real::from_i64(if attempt % 2 == 0 {
            -magnitude
        } else {
            magnitude
        });
        let reference = reference
            .scaled(&(Real::one() - &t * &t))
            .plus(&across.scaled(&(&Real::from_i64(2) * &t)));
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
                        crossings.push(Crossing {
                            edge: boundary.edge,
                            point: hit.point,
                            height: radial.dot(&axis),
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
        trim::try_sort(&mut crossings, |a, b| {
            Ok(match trim::sign(&(&a.height - &b.height), dec)? {
                Sign::Negative => std::cmp::Ordering::Less,
                Sign::Zero => std::cmp::Ordering::Equal,
                Sign::Positive => std::cmp::Ordering::Greater,
            })
        })?;
        for pair in crossings.windows(2) {
            if trim::sign(&(&pair[0].height - &pair[1].height), dec)? == Sign::Zero {
                regular = false;
                break;
            }
        }
        if regular {
            return Ok(Meridian {
                reference,
                crossings,
            });
        }
    }
    Err(number::budget(
        "no transverse periodic meridian certificate",
    ))
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
