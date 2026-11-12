//! Export-only caps around spherical parameter poles. The raw translator needs
//! a degenerate edge for a pole inside a face, and STEP has no such edge. A
//! latitude ring strictly inside the face moves each interior pole onto its own
//! cap, written with a perpendicular frame whose poles lie outside that cap.
//! A pole that the model already marks with a vertex loop is that same pole:
//! the loop is dropped, and the cap carries the pole instead.

use crate::{number, trim};
use catallax_intersect::CurveSurface;
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Rational, Real, Refusal, Sign};
use catallax_shapes::{Circle, Curve, GeomKey, Point3, Sphere, Surface, SurfaceId};
use catallax_topology::{BodyId, Coedge, Edge, Face, FaceId, Loop, Model, Sense, Shell};
use std::collections::BTreeMap;

/// Ring heights tried per pole, from a half-radius cap towards the pole.
const CAP_HEIGHTS: u32 = 24;

pub(crate) struct Split {
    pub model: Model,
    pub bodies: Vec<BodyId>,
    /// Source face to its holed replacement.
    pub faces: BTreeMap<FaceId, FaceId>,
    /// Export-only ring edges added around poles.
    pub rings: usize,
}

/// `poles` lists each face with the poles certified inside it, low then high
/// along the sphere's axis. `None` keeps the legacy encoding.
pub(crate) fn split(
    model: &Model,
    bodies: &[BodyId],
    poles: &[(FaceId, [bool; 2])],
    geometry: &mut trim::Geometry,
    dec: &mut Decider,
) -> Result<Option<Split>, Refusal> {
    let mut rings = Vec::new();
    for &(face, inside) in poles {
        let Surface::Sphere(sphere) = geometry.store.surface(model.face(face).surface).clone()
        else {
            unreachable!("pole on a non-spherical face")
        };
        for (side, inside) in [(-1, inside[0]), (1, inside[1])] {
            if !inside {
                continue;
            }
            let Some(lambda) = ring(model, face, &sphere, side, geometry, dec)? else {
                return Ok(None);
            };
            rings.push((face, sphere.clone(), side, lambda));
        }
    }
    let mut split = model.clone();
    *split.geom_mut() = geometry.store.clone();
    let mut frames = BTreeMap::<SurfaceId, SurfaceId>::new();
    let mut holes = BTreeMap::<FaceId, Vec<_>>::new();
    let mut caps = BTreeMap::<FaceId, Vec<FaceId>>::new();
    let count = rings.len();
    for (face, sphere, side, lambda) in rings {
        let source = model.face(face).clone();
        let frame = *frames.entry(source.surface).or_insert_with(|| {
            split.geom_mut().add_surface(Surface::Sphere(Sphere {
                center: sphere.center.clone(),
                radius: sphere.radius.clone(),
                axis: sphere.ref_dir.clone(),
                ref_dir: sphere.axis.clone(),
            }))
        });
        let lambda = Real::from_rational(lambda);
        let height = &(&lambda * &lambda) * &sphere.axis.norm_squared();
        let circle = split.geom_mut().add_curve(Curve::Circle(Circle {
            center: sphere.center.plus(&sphere.axis.scaled(&lambda)),
            normal: sphere.axis.clone(),
            radius_sq: &(&sphere.radius * &sphere.radius) - &height,
            ref_dir: sphere.ref_dir.clone(),
        }));
        // The ring is the sphere's section by its own cutting plane.
        for surface in [source.surface, frame] {
            split.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(circle),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("step.pole-cap-ring"),
            });
        }
        let edge = split.add_edge(Edge {
            curve: circle,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        // The circle turns counterclockwise about the axis. Its cap lies to
        // the left of travel, seen against the face normal, when the pole side
        // and the face sense agree.
        let cap_sense = if (side > 0) == (source.sense == Sense::Same) {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let cap_loop = split.add_coedge(Coedge {
            edge,
            sense: cap_sense,
        });
        let cap_loop = split.add_loop(Loop::Coedges(vec![cap_loop]));
        caps.entry(face).or_default().push(split.add_face(Face {
            surface: frame,
            sense: source.sense,
            loops: vec![cap_loop],
        }));
        let hole = split.add_coedge(Coedge {
            edge,
            sense: match cap_sense {
                Sense::Same => Sense::Reversed,
                Sense::Reversed => Sense::Same,
            },
        });
        holes
            .entry(face)
            .or_default()
            .push(split.add_loop(Loop::Coedges(vec![hole])));
    }
    let mut faces = BTreeMap::new();
    for (face, mut loops) in holes {
        let mut definition = model.face(face).clone();
        definition.loops.append(&mut loops);
        faces.insert(face, split.add_face(definition));
    }
    let mut replacements = Vec::with_capacity(bodies.len());
    for &body in bodies {
        let mut body = model.body(body).clone();
        for id in &mut body.lumps {
            let mut lump = model.lump(*id).clone();
            for id in &mut lump.shells {
                let mut shell = Vec::new();
                for face in &model.shell(*id).faces {
                    shell.push(faces.get(face).copied().unwrap_or(*face));
                    shell.extend(caps.get(face).into_iter().flatten().copied());
                }
                *id = split.add_shell(Shell { faces: shell });
            }
            *id = split.add_lump(lump);
        }
        replacements.push(split.add_body(body));
    }
    Ok(Some(Split {
        model: split,
        bodies: replacements,
        faces,
        rings: count,
    }))
}

/// The model without the vertex loops of spherical faces, each of which
/// must sit on a pole of its sphere beside a boundary of edges. `None` when
/// no face has one, or when one is not at a pole.
pub(crate) fn strip(
    model: &Model,
    bodies: &[BodyId],
    dec: &mut Decider,
) -> Result<Option<Split>, Refusal> {
    let mut stripped = BTreeMap::new();
    let mut replaced = model.clone();
    let mut seen = std::collections::BTreeSet::new();
    for &body in bodies {
        for face in model.faces_of_body(body) {
            if !seen.insert(face) {
                continue;
            }
            let definition = model.face(face);
            let Surface::Sphere(sphere) = model.geom().surface(definition.surface) else {
                continue;
            };
            let mut loops = Vec::with_capacity(definition.loops.len());
            let mut vertices = Vec::new();
            for id in &definition.loops {
                match model.loop_(*id) {
                    Loop::Vertex(vertex) => vertices.push(*vertex),
                    Loop::Coedges(_) => loops.push(*id),
                }
            }
            if vertices.is_empty() {
                continue;
            }
            if loops.is_empty() {
                return Ok(None);
            }
            let radius = number::unit(&sphere.axis).scaled(&sphere.radius);
            for vertex in vertices {
                let at = model.geom().position(model.vertex(vertex).point);
                let mut on_pole = false;
                for sign in [-1, 1] {
                    let pole = sphere.center.plus(&radius.scaled(&Real::from_i64(sign)));
                    on_pole |= trim::sign(&at.minus(&pole).norm_squared(), dec)? == Sign::Zero;
                }
                if !on_pole {
                    return Ok(None);
                }
            }
            let mut definition = definition.clone();
            definition.loops = loops;
            stripped.insert(face, replaced.add_face(definition));
        }
    }
    if stripped.is_empty() {
        return Ok(None);
    }
    let mut replacements = Vec::with_capacity(bodies.len());
    for &body in bodies {
        let mut body = model.body(body).clone();
        for id in &mut body.lumps {
            let mut lump = model.lump(*id).clone();
            for id in &mut lump.shells {
                let faces = model
                    .shell(*id)
                    .faces
                    .iter()
                    .map(|face| stripped.get(face).copied().unwrap_or(*face))
                    .collect();
                *id = replaced.add_shell(Shell { faces });
            }
            *id = replaced.add_lump(lump);
        }
        replacements.push(replaced.add_body(body));
    }
    Ok(Some(Split {
        model: replaced,
        bodies: replacements,
        faces: stripped,
        rings: 0,
    }))
}

/// The largest tried ring `axis·(x - center) = lambda·|axis|²` whose closed cap
/// on `side` meets no boundary of `face`. The pole is certified inside `face`,
/// so that cap lies strictly inside the face.
fn ring(
    model: &Model,
    face: FaceId,
    sphere: &Sphere,
    side: i64,
    geometry: &mut trim::Geometry,
    dec: &mut Decider,
) -> Result<Option<Rational>, Refusal> {
    let norm_squared = sphere.axis.norm_squared();
    // A rational multiple of the axis keeps every cut plane free of the
    // axis length's square root. Any rational below radius/|axis| works.
    let ratio = &sphere.radius * &(Real::one() / number::sqrt_real(&norm_squared));
    let coarse = dec.enclose(&ratio, &Rational::one())?.hi().abs();
    let scale = dec
        .enclose(
            &ratio,
            &(coarse * Rational::from_ratio(1, 1 << 20).unwrap()),
        )?
        .lo();
    if scale.sign() != Sign::Positive {
        return Ok(None);
    }
    let boundaries = geometry.boundaries(model, face, dec)?;
    'heights: for k in 1..=CAP_HEIGHTS {
        dec.check_deadline()?;
        let fraction = Rational::one() - Rational::from_ratio(1, 1 << k).unwrap();
        let lambda = &(&scale * &fraction) * &Rational::from_i64(side);
        let offset = Real::from_rational(lambda.clone());
        let plane = geometry.plane(
            sphere.center.plus(&sphere.axis.scaled(&offset)),
            sphere.axis.clone(),
            sphere.ref_dir.clone(),
            dec,
        )?;
        let level = &offset * &norm_squared;
        for boundary in boundaries.iter() {
            match geometry.intersections(boundary.curve, plane, dec)?.as_ref() {
                CurveSurface::Contained => continue 'heights,
                CurveSurface::Disjoint => {}
                CurveSurface::Meets(hits) => {
                    for hit in hits {
                        if geometry.on_arc(boundary, hit.point, dec)? != trim::Location::Outside {
                            continue 'heights;
                        }
                    }
                }
            }
            // The arc does not reach the ring, so one point decides its side.
            let height = sphere
                .axis
                .dot(&sample(geometry, boundary).minus(&sphere.center));
            let below = trim::sign(&(&(&height - &level) * &Real::from_i64(side)), dec)?;
            if below != Sign::Negative {
                continue 'heights;
            }
        }
        return Ok(Some(lambda));
    }
    Ok(None)
}

fn sample(geometry: &trim::Geometry, boundary: &trim::Boundary) -> Point3 {
    if let Some(point) = boundary.from {
        return geometry.store.position(point);
    }
    match geometry.store.curve(boundary.curve) {
        Curve::Intersection(curve) => geometry.store.position(curve.witness),
        curve => curve.parameter_origin(&geometry.store),
    }
}
