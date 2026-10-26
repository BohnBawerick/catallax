//! Export-only cuts for complete cylindrical and conical bands.

use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, FactKind, Outcome, Real, Refusal, StrictOrder, StrictSign};
use catallax_shapes::{CurveId, GeomKey, Point3, Surface, SurfaceId, Vector3};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

use crate::number;

pub(crate) struct Association {
    pub surface: SurfaceId,
    pub v: Real,
    pub v_scale: Real,
    pub increasing: bool,
}

pub(crate) struct Ring {
    pub origin: Point3,
    pub radius: Real,
    pub associations: Vec<Association>,
}

pub(crate) struct Band {
    pub surface: SurfaceId,
    pub edges: [EdgeId; 2],
    pub forward: [bool; 2],
    pub low_v: Real,
    pub v_scale: Real,
    pub radius_bound: Real,
    pub seam_direction: Vector3,
}

#[derive(Default)]
pub(crate) struct Plan {
    pub faces: BTreeMap<FaceId, Rc<Band>>,
    pub rings: BTreeMap<EdgeId, Rc<Ring>>,
    pub curves: BTreeMap<CurveId, CurveForm>,
    pub references: BTreeMap<SurfaceId, Vector3>,
}

struct SourceRing {
    edge: EdgeId,
    curve: CurveId,
    center: Point3,
    normal: Vector3,
    radius_sq: Real,
    radius: Real,
    v: Real,
    increasing: bool,
}

struct SourceBand {
    face: FaceId,
    surface: SurfaceId,
    axis: Vector3,
    reference: Vector3,
    v_scale: Real,
    rings: [SourceRing; 2],
}

impl Plan {
    pub(crate) fn prepare(
        model: &Model,
        bodies: impl Iterator<Item = BodyId> + Clone,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut has_band = false;
        for body in bodies.clone() {
            for lump in &model.body(body).lumps {
                for shell in &model.lump(*lump).shells {
                    for &face in &model.shell(*shell).faces {
                        let definition = model.face(face);
                        match model.geom().surface(definition.surface) {
                            Surface::Plane(_) => {}
                            Surface::Sphere(_) | Surface::Torus(_)
                                if definition.loops.is_empty() => {}
                            Surface::Cylinder(_) | Surface::Cone(_)
                                if ring_edges(model, face).is_some() =>
                            {
                                has_band = true
                            }
                            _ => return Ok(Self::default()),
                        }
                    }
                }
            }
        }
        if !has_band {
            return Ok(Self::default());
        }
        let mut sources = Vec::new();
        let mut seen = BTreeSet::new();
        for body in bodies {
            for lump in &model.body(body).lumps {
                for shell in &model.lump(*lump).shells {
                    for &face in &model.shell(*shell).faces {
                        if !seen.insert(face) {
                            continue;
                        }
                        let definition = model.face(face);
                        match model.geom().surface(definition.surface) {
                            Surface::Plane(_) => {}
                            Surface::Sphere(_) | Surface::Torus(_)
                                if definition.loops.is_empty() => {}
                            Surface::Cylinder(_) | Surface::Cone(_) => {
                                let Some(band) = source_band(model, face, dec)? else {
                                    // Do not partly change a file that still needs the
                                    // arbitrary-trim or singular-apex export contract.
                                    return Ok(Self::default());
                                };
                                sources.push(band);
                            }
                            _ => return Ok(Self::default()),
                        }
                    }
                }
            }
        }
        // Connected bands must cut shared rings at the same physical meridian.
        // A reference rotation changes the parametrisation, not the carrier.
        let mut parents: Vec<_> = (0..sources.len()).collect();
        let mut edge_owners = BTreeMap::new();
        let mut curve_owners = BTreeMap::new();
        let mut surface_owners = BTreeMap::new();
        for (index, band) in sources.iter().enumerate() {
            if let Some(previous) = surface_owners.insert(band.surface, index) {
                join(&mut parents, index, previous);
            }
            for ring in &band.rings {
                if let Some(previous) = edge_owners.insert(ring.edge, index) {
                    join(&mut parents, index, previous);
                }
                if let Some(previous) = curve_owners.insert(ring.curve, index) {
                    join(&mut parents, index, previous);
                }
            }
        }
        let mut plan = Self::default();
        let mut rings: BTreeMap<EdgeId, Ring> = BTreeMap::new();
        for (index, source) in sources.iter().enumerate() {
            let reference = &sources[root(&mut parents, index)].reference;
            plan.references.insert(source.surface, reference.clone());
            // Use the certified carrier frame, not a normalized difference
            // of algebraic endpoint heights that must be re-proved unit.
            let seam_direction = match model.geom().surface(source.surface) {
                Surface::Cylinder(_) => source.axis.clone(),
                Surface::Cone(cone) => source.axis.plus(&reference.scaled(&cone.slope)),
                _ => unreachable!(),
            };
            plan.faces.insert(
                source.face,
                Rc::new(Band {
                    surface: source.surface,
                    edges: source.rings.each_ref().map(|r| r.edge),
                    forward: source
                        .rings
                        .each_ref()
                        .map(|r| r.increasing == (model.edge(r.edge).sense == Sense::Same)),
                    low_v: source.rings[0].v.clone(),
                    v_scale: source.v_scale.clone(),
                    radius_bound: &source.rings[0].radius + &source.rings[1].radius,
                    seam_direction,
                }),
            );
            for ring in &source.rings {
                let origin = number::point_plus(&ring.center, reference, &ring.radius);
                rings
                    .entry(ring.edge)
                    .or_insert_with(|| Ring {
                        origin,
                        radius: ring.radius.clone(),
                        associations: Vec::with_capacity(2),
                    })
                    .associations
                    .push(Association {
                        surface: source.surface,
                        v: ring.v.clone(),
                        v_scale: source.v_scale.clone(),
                        increasing: ring.increasing,
                    });
                plan.curves
                    .entry(ring.curve)
                    .or_insert_with(|| CurveForm::Circle {
                        center: ring.center.clone(),
                        normal: ring.normal.clone(),
                        radius_sq: ring.radius_sq.clone(),
                        ref_dir: reference.clone(),
                    });
            }
        }
        plan.rings = rings
            .into_iter()
            .map(|(edge, ring)| (edge, Rc::new(ring)))
            .collect();
        Ok(plan)
    }
}

fn ring_edges(model: &Model, face: FaceId) -> Option<[EdgeId; 2]> {
    let [a, b] = model.face(face).loops.as_slice() else {
        return None;
    };
    let edge = |loop_id| {
        let Loop::Coedges(coedges) = model.loop_(loop_id) else {
            return None;
        };
        let [coedge] = coedges.as_slice() else {
            return None;
        };
        let edge = model.coedge(*coedge).edge;
        let definition = model.edge(edge);
        (definition.start.is_none() && definition.end.is_none()).then_some(edge)
    };
    Some([edge(*a)?, edge(*b)?])
}

pub(crate) fn is_complete(model: &Model, face: FaceId, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(source_band(model, face, dec)?.is_some())
}

fn source_band(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<SourceBand>, Refusal> {
    let face_definition = model.face(face);
    let Some(edges) = ring_edges(model, face) else {
        return Ok(None);
    };
    let surface = model.geom().surface(face_definition.surface);
    let (origin, axis, reference, slant) = match surface {
        Surface::Cylinder(s) => (
            &s.origin,
            number::unit(&s.axis),
            number::unit(&s.ref_dir),
            Real::one(),
        ),
        Surface::Cone(s) => (
            &s.apex,
            number::unit(&s.axis),
            number::unit(&s.ref_dir),
            number::sqrt_real(&(Real::one() + &s.slope * &s.slope)),
        ),
        _ => return Ok(None),
    };
    let mut rings = [None, None];
    for (index, edge) in edges.into_iter().enumerate() {
        let definition = model.edge(edge);
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = curve_form(model.geom(), definition.curve, dec)?
        else {
            return Ok(None);
        };
        let offset = center.minus(origin);
        let height = offset.dot(&axis);
        let curve_key = GeomKey::Curve(definition.curve);
        let surface_key = GeomKey::Surface(face_definition.surface);
        if model
            .geom()
            .facts()
            .holds(FactKind::Coaxial, curve_key, surface_key)
        {
            dec.note_fact_decision();
        } else if !zero(&normal.cross(&axis).norm_squared(), dec)?
            || !zero(&offset.cross(&axis).norm_squared(), dec)?
        {
            return Ok(None);
        }
        if matches!(surface, Surface::Cone(_))
            && !matches!(dec.sign(&height), Outcome::Proven(StrictSign::Positive))
        {
            return Ok(None);
        }
        // Whole-carrier incidence already proves the circular section's
        // radius equation. Re-eliminating an algebraic drill angle here can
        // exhaust the request after construction has certified that identity.
        if model
            .geom()
            .facts()
            .holds(FactKind::Incident, curve_key, surface_key)
        {
            dec.note_fact_decision();
        } else {
            let expected_radius_sq = match surface {
                Surface::Cylinder(s) => s.radius_sq.clone(),
                Surface::Cone(s) => (&s.slope * &height).square(),
                _ => unreachable!(),
            };
            if !zero(&(&radius_sq - &expected_radius_sq), dec)? {
                return Ok(None);
            }
        }
        let increasing = match dec.sign(&normal.dot(&axis)) {
            Outcome::Proven(sign) => sign == StrictSign::Positive,
            Outcome::Refused(refusal) => return Err(refusal),
            _ => return Ok(None),
        };
        rings[index] = Some(SourceRing {
            edge,
            curve: definition.curve,
            center,
            normal,
            radius: number::sqrt_real(&radius_sq),
            radius_sq,
            v: height,
            increasing,
        });
    }
    let mut rings = rings.map(Option::unwrap);
    match dec.compare(&rings[0].v, &rings[1].v) {
        Outcome::Proven(StrictOrder::Greater) => rings.swap(0, 1),
        Outcome::Proven(StrictOrder::Less) => {}
        Outcome::Refused(refusal) => return Err(refusal),
        _ => return Ok(None),
    }
    Ok(Some(SourceBand {
        face,
        surface: face_definition.surface,
        axis,
        reference,
        v_scale: slant,
        rings,
    }))
}

fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::Refused(refusal) => Err(refusal),
        _ => Ok(false),
    }
}

fn root(parents: &mut [usize], mut index: usize) -> usize {
    while parents[index] != index {
        parents[index] = parents[parents[index]];
        index = parents[index];
    }
    index
}

fn join(parents: &mut [usize], a: usize, b: usize) {
    let a = root(parents, a);
    let b = root(parents, b);
    parents[b] = a;
}
