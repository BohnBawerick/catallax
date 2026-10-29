//! Export-only sphere frames. The frame only parameterizes the sphere, so a
//! trimmed face may be written with another exact frame of the same centre and
//! radius. The point set, the face sense and the loops are unchanged, so each
//! trim keeps its incidence.

use crate::{number, trim};
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Real, Refusal};
use catallax_shapes::{GeomKey, Sphere, Surface, SurfaceId, Vector3};
use catallax_topology::{BodyId, Face, FaceId, Loop, Model, Shell};
use std::collections::{BTreeMap, BTreeSet};

pub(crate) struct Reframed {
    pub model: Model,
    pub bodies: Vec<BodyId>,
    /// Source face to its face on the new frame.
    pub faces: BTreeMap<FaceId, FaceId>,
}

/// Every trimmed spherical face moved onto its sphere's perpendicular frame,
/// whose axis is the source reference direction. Old poles lie on the new
/// frame's equator, away from trims that pass close to them.
pub(crate) fn perpendicular(model: &Model, bodies: &[BodyId]) -> Reframed {
    let mut frames = BTreeMap::new();
    for (surface, _) in trimmed(model, bodies, |_| true) {
        let Surface::Sphere(sphere) = model.geom().surface(surface) else {
            unreachable!("trimmed faces are spherical")
        };
        frames.insert(
            surface,
            Sphere {
                axis: sphere.ref_dir.clone(),
                ref_dir: sphere.axis.clone(),
                ..sphere.clone()
            },
        );
    }
    rewrite(model, bodies, frames)
}

/// A stored frame whose pole lies on a trim has no admissible chart: the pole
/// is neither outside the face nor strictly inside it, where an export-only
/// cap could hold it. Every face on `spheres` is moved onto the first candidate
/// frame whose poles are certified outside each of those faces, or inside when
/// `split_poles` allows a cap. `None` keeps the legacy encoding.
pub(crate) fn reframe(
    model: &Model,
    bodies: &[BodyId],
    spheres: &BTreeSet<SurfaceId>,
    split_poles: bool,
    geometry: &mut trim::Geometry,
    dec: &mut Decider,
) -> Result<Option<Reframed>, Refusal> {
    let mut frames = BTreeMap::new();
    'spheres: for (surface, faces) in trimmed(model, bodies, |surface| spheres.contains(&surface)) {
        let Surface::Sphere(sphere) = geometry.store.surface(surface).clone() else {
            unreachable!("trimmed faces are spherical")
        };
        'candidates: for (axis, ref_dir) in candidates(&sphere) {
            let radius = number::unit(&axis).scaled(&sphere.radius);
            let poles = [
                geometry.point(sphere.center.plus(&radius.scaled(&(-Real::one())))),
                geometry.point(sphere.center.plus(&radius)),
            ];
            // Membership depends only on the boundary, not on the stored frame.
            for &face in &faces {
                for pole in poles {
                    match geometry.sphere_contains(model, face, pole, dec)? {
                        Some(trim::Location::Outside) => {}
                        Some(trim::Location::Inside) if split_poles => {}
                        _ => continue 'candidates,
                    }
                }
            }
            frames.insert(
                surface,
                Sphere {
                    axis,
                    ref_dir,
                    ..sphere
                },
            );
            continue 'spheres;
        }
        return Ok(None);
    }
    Ok(Some(rewrite(model, bodies, frames)))
}

/// Exact frames in a fixed order: the stored axes permuted, then their
/// diagonals. Every pair is perpendicular because the stored pair is.
fn candidates(sphere: &Sphere) -> Vec<(Vector3, Vector3)> {
    let a = &sphere.axis;
    let r = &sphere.ref_dir;
    let c = a.cross(r);
    vec![
        (r.clone(), a.clone()),
        (c.clone(), a.clone()),
        (a.plus(r), c.clone()),
        (a.minus(r), c.clone()),
        (a.plus(&c), r.clone()),
        (a.minus(&c), r.clone()),
        (r.plus(&c), a.clone()),
        (r.minus(&c), a.clone()),
    ]
}

/// Trimmed spherical faces of `bodies` on the selected surfaces, by surface.
fn trimmed(
    model: &Model,
    bodies: &[BodyId],
    select: impl Fn(SurfaceId) -> bool,
) -> BTreeMap<SurfaceId, Vec<FaceId>> {
    let mut faces = BTreeMap::<SurfaceId, Vec<FaceId>>::new();
    for &body in bodies {
        for face in model.faces_of_body(body) {
            let definition = model.face(face);
            if !definition.loops.is_empty()
                && matches!(model.geom().surface(definition.surface), Surface::Sphere(_))
                && select(definition.surface)
            {
                let faces = faces.entry(definition.surface).or_default();
                if !faces.contains(&face) {
                    faces.push(face);
                }
            }
        }
    }
    faces
}

/// Every trimmed face on a key of `frames` moved onto that frame.
fn rewrite(model: &Model, bodies: &[BodyId], frames: BTreeMap<SurfaceId, Sphere>) -> Reframed {
    let mut moved = model.clone();
    let mut faces = BTreeMap::new();
    let mut trims = trimmed(model, bodies, |surface| frames.contains_key(&surface));
    for (surface, sphere) in frames {
        let frame = moved.geom_mut().add_surface(Surface::Sphere(sphere));
        for face in trims.remove(&surface).into_iter().flatten() {
            let source = model.face(face).clone();
            // The implicit sphere is unchanged, so every source incidence holds.
            for &id in &source.loops {
                let Loop::Coedges(coedges) = model.loop_(id) else {
                    continue;
                };
                for &coedge in coedges {
                    let curve = GeomKey::Curve(model.edge(model.coedge(coedge).edge).curve);
                    if model.geom().facts().holds(
                        FactKind::Incident,
                        curve,
                        GeomKey::Surface(surface),
                    ) {
                        moved.geom_mut().record_fact(Fact {
                            kind: FactKind::Incident,
                            a: curve,
                            b: GeomKey::Surface(frame),
                            at: None,
                            origin: FactOrigin::Construction("step.sphere-frame"),
                        });
                    }
                }
            }
            faces.insert(
                face,
                moved.add_face(Face {
                    surface: frame,
                    ..source
                }),
            );
        }
    }
    let mut replacements = Vec::with_capacity(bodies.len());
    for &body in bodies {
        let mut body = model.body(body).clone();
        for id in &mut body.lumps {
            let mut lump = model.lump(*id).clone();
            for id in &mut lump.shells {
                let shell = model
                    .shell(*id)
                    .faces
                    .iter()
                    .map(|face| faces.get(face).copied().unwrap_or(*face))
                    .collect();
                *id = moved.add_shell(Shell { faces: shell });
            }
            *id = moved.add_lump(lump);
        }
        replacements.push(moved.add_body(body));
    }
    Reframed {
        model: moved,
        bodies: replacements,
        faces,
    }
}
