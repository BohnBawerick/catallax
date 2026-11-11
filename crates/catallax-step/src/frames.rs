//! Export-only sphere frames. The frame only parameterizes the sphere, so a
//! trimmed face may be written with another exact frame of the same centre and
//! radius. The point set, the face sense and the loops are unchanged, so each
//! trim keeps its incidence.

use crate::{number, trim};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, Outcome, Rational, Real, Refusal, StrictSign,
};
use catallax_shapes::{
    BoundaryRow, Cylinder, GeomKey, SaddleRing, Sphere, SplineSupport, Surface, SurfaceId, Vector3,
};
use catallax_topology::{
    BodyId, EdgeId, Face, FaceId, Loop, Model, RollingMitreBand, Shell, pipe_rows, saddle_rows,
};
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
    let spherical = |surface| matches!(model.geom().surface(surface), Surface::Sphere(_));
    for (surface, _) in trimmed(model, bodies, spherical) {
        let Surface::Sphere(sphere) = model.geom().surface(surface) else {
            unreachable!("trimmed faces are spherical")
        };
        frames.insert(
            surface,
            Surface::Sphere(Sphere {
                axis: sphere.ref_dir.clone(),
                ref_dir: sphere.axis.clone(),
                ..sphere.clone()
            }),
        );
    }
    rewrite(model, bodies, frames)
}

/// Every wall under a band, moved onto the frame its chart needs. A rolling
/// pipe's wall takes the band's `t = 0` direction, so its first chart cut meets
/// the wall rail at the band's seam. A saddle's narrow rail winds round its
/// cylinder and takes the band's seam at `t = -pi/2`. A saddle's wide rail
/// stays strictly on one side of the plane through the wide axis along
/// `across`, so a cut along `across` misses every ring on that cylinder. A loop
/// winds round neither cylinder, so each wall's cut goes where its rail never
/// does. A band over a mitre has rail `i` on round `i` at `R (cos t w + sin t
/// y)` from its axis, for `t` in `[0, pi/2]`, with `w` the foot's unit radial
/// direction and `y` the apex's, square to each other. A cut plane through the
/// axis along `d` misses that closed quarter both ways when `(d . w) (d . y) <
/// 0`, so the wall takes `d = y - w` unless an earlier band's cut already
/// misses it. Its end section's wall takes the band's side direction `z`.
/// `None` when the bodies have no band.
pub(crate) fn band_walls(
    model: &Model,
    bodies: &[BodyId],
    dec: &mut Decider,
) -> Result<Option<Reframed>, Refusal> {
    let refused = || {
        Refusal::new(
            "step.unsupported-geometry",
            "a trimmed band has no STEP form",
        )
    };
    let mut frames = BTreeMap::new();
    let frame = |rail: EdgeId, direction: Vector3, frames: &mut BTreeMap<_, _>| {
        for coedge in model.coedges_of_edge(rail) {
            let Some((_, neighbour)) = model.owner_of_coedge(coedge) else {
                continue;
            };
            let surface = model.face(neighbour).surface;
            match model.geom().surface(surface) {
                Surface::Cylinder(wall) => {
                    frames.insert(
                        surface,
                        Surface::Cylinder(Cylinder {
                            ref_dir: direction.clone(),
                            ..wall.clone()
                        }),
                    );
                }
                // A cone's band puts its wall rail's seam vertex `across`
                // from the axis too.
                Surface::Cone(wall) => {
                    frames.insert(
                        surface,
                        Surface::Cone(catallax_shapes::Cone {
                            ref_dir: direction.clone(),
                            ..wall.clone()
                        }),
                    );
                }
                _ => {}
            }
        }
    };
    for &body in bodies {
        for face in model.faces_of_body(body) {
            match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RollingPipe(pipe)) => {
                    let [_, (rail, _)] = pipe_rows(model, face, dec)?.ok_or_else(refused)?;
                    frame(rail, pipe.law().across.clone(), &mut frames);
                }
                Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => {
                    let [(wide, _), (narrow, _)] = saddle_rows(model, face).ok_or_else(refused)?;
                    let law = saddle.definition().law();
                    if matches!(law.ring, SaddleRing::Crossing { .. }) {
                        // Both rails meet the common perpendicular at their
                        // poles, so neither crosses this cut in its interior.
                        frame(wide, law.across.clone(), &mut frames);
                        frame(narrow, law.across.clone(), &mut frames);
                        continue;
                    }
                    if law.ring == SaddleRing::Overlap {
                        // A loop's wide rail never reaches `-e across`, where
                        // `X = -P`, and its narrow rail never reaches `e across`
                        // from the foot, where `X - g = p`. So each wall's cut
                        // goes there and the rail misses it.
                        let facing = Real::from_i64(law.gap_sign);
                        frame(wide, law.across.scaled(&-&facing), &mut frames);
                        frame(narrow, law.across.scaled(&facing), &mut frames);
                        continue;
                    }
                    frame(wide, law.across.clone(), &mut frames);
                    frame(narrow, law.side_axis().scaled(&-Real::one()), &mut frames);
                }
                Surface::BSpline(SplineSupport::RollingMitre(band)) => {
                    let checked = RollingMitreBand::checked(model, face, dec)?
                        .map_err(|undecided| refused().keeping_any(undecided))?;
                    let [out, up, side] = band.frame();
                    let apex = up.scaled(band.round());
                    let (along, across) = (band.along().square(), band.across().square());
                    for (edge, row, _, _) in checked.rows() {
                        if row == BoundaryRow::Negative {
                            // The section's sector `a m1 + b m2`, `a, b >= 0`,
                            // has `x` component `(a + b) sa > 0`. A cut
                            // `d = p x + q z` misses it both ways when
                            // `p^2 sc^2 < q^2 sa^2`, as `z` does.
                            mitre_wall(
                                model,
                                edge,
                                side,
                                &mut frames,
                                |cut, dec| {
                                    let p = cut.dot(out).square();
                                    let q = cut.dot(side).square();
                                    dec.sign(&(&(&p * &across) - &(&q * &along)))
                                },
                                dec,
                            )?;
                            continue;
                        }
                        let i = usize::from(row == BoundaryRow::Mouth);
                        let axis = &band.axes()[i];
                        let foot = band.rail(i, dec)?.major_dir;
                        let foot = foot.minus(&axis.scaled(&foot.dot(axis)));
                        let cut = apex.minus(&foot);
                        mitre_wall(
                            model,
                            edge,
                            &cut,
                            &mut frames,
                            |cut, dec| dec.sign(&(&cut.dot(&foot) * &cut.dot(&apex))),
                            dec,
                        )?;
                    }
                }
                Surface::BSpline(SplineSupport::RollingSphere(band)) => {
                    let [(round, _), (tube, _)] = saddle_rows(model, face).ok_or_else(refused)?;
                    let def = band.definition();
                    // The cylinder rail's angle is the chart's own, which
                    // starts at -pi/2.
                    frame(tube, def.along().scaled(&-Real::one()), &mut frames);
                    // The sphere turns about the cylinder's axis, and its cut
                    // meets the sphere rail where that starts. Both rings of
                    // one sphere start at the same angle, or neither is cut.
                    let start = band.rails()[0].control(0).minus(def.centre());
                    let start = start.minus(&def.axis().scaled(&start.dot(def.axis())));
                    for coedge in model.coedges_of_edge(round) {
                        let Some((_, neighbour)) = model.owner_of_coedge(coedge) else {
                            continue;
                        };
                        let surface = model.face(neighbour).surface;
                        let Surface::Sphere(sphere) = model.geom().surface(surface) else {
                            continue;
                        };
                        if let Some(Surface::Sphere(earlier)) = frames.get(&surface) {
                            // The same start direction and axis, checked in
                            // order; the first that fails refuses and keeps
                            // its decision if it was undecided.
                            let checks = [
                                (earlier.ref_dir.cross(&start).norm_squared(), None),
                                (earlier.ref_dir.dot(&start), Some(StrictSign::Positive)),
                                (earlier.axis.cross(def.axis()).norm_squared(), None),
                            ];
                            for (value, wanted) in checks {
                                let outcome = dec.sign(&value);
                                let holds = match (&outcome, wanted) {
                                    (Outcome::ProvenCoincident(_), None) => true,
                                    (Outcome::Proven(sign), Some(wanted)) => *sign == wanted,
                                    _ => false,
                                };
                                if !holds {
                                    return Err(refused().keeping_undecided(outcome));
                                }
                            }
                            continue;
                        }
                        frames.insert(
                            surface,
                            Surface::Sphere(Sphere {
                                axis: def.axis().clone(),
                                ref_dir: start.clone(),
                                ..sphere.clone()
                            }),
                        );
                    }
                }
                _ => {}
            }
        }
    }
    if frames.is_empty() {
        return Ok(None);
    }
    Ok(Some(rewrite(model, bodies, frames)))
}

/// Moves the cylinder under a mitre band's boundary `edge` onto `cut`, or
/// keeps an earlier frame when `misses` proves its cut clear of that
/// boundary. Any other frame there refuses.
fn mitre_wall(
    model: &Model,
    edge: EdgeId,
    cut: &Vector3,
    frames: &mut BTreeMap<SurfaceId, Surface>,
    // The sign of a quantity that is negative exactly when the cut misses.
    misses: impl Fn(&Vector3, &mut Decider) -> Outcome<StrictSign>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for coedge in model.coedges_of_edge(edge) {
        let Some((_, neighbour)) = model.owner_of_coedge(coedge) else {
            continue;
        };
        let surface = model.face(neighbour).surface;
        let Surface::Cylinder(wall) = model.geom().surface(surface) else {
            continue;
        };
        let different = || {
            Refusal::new(
                "step.unsupported-geometry",
                "two bands need different cuts on one wall",
            )
        };
        match frames.get(&surface) {
            Some(Surface::Cylinder(framed)) => match misses(&framed.ref_dir, dec) {
                Outcome::Proven(StrictSign::Negative) => {}
                other => return Err(different().keeping_undecided(other)),
            },
            Some(_) => return Err(different()),
            None => {
                frames.insert(
                    surface,
                    Surface::Cylinder(Cylinder {
                        ref_dir: cut.clone(),
                        ..wall.clone()
                    }),
                );
            }
        }
    }
    Ok(())
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
        let triangular = faces.iter().all(|&face| {
            let loops = &model.face(face).loops;
            loops.len() == 1
                && matches!(model.loop_(loops[0]), Loop::Coedges(edges) if edges.len() == 3)
        });
        'candidates: for (axis, ref_dir) in candidates(&sphere, triangular) {
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
                Surface::Sphere(Sphere {
                    axis,
                    ref_dir,
                    ..sphere
                }),
            );
            continue 'spheres;
        }
        return Ok(None);
    }
    Ok(Some(rewrite(model, bodies, frames)))
}

/// Stored axes first, with rational turns before the diagonals on triangle trims.
/// Every pair is perpendicular because the stored pair is.
fn candidates(sphere: &Sphere, rational_first: bool) -> Vec<(Vector3, Vector3)> {
    let a = &sphere.axis;
    let r = &sphere.ref_dir;
    let c = a.cross(r);
    let mut frames = vec![
        (r.clone(), a.clone()),
        (c.clone(), a.clone()),
        (a.plus(r), c.clone()),
        (a.minus(r), c.clone()),
        (a.plus(&c), r.clone()),
        (a.minus(&c), r.clone()),
        (r.plus(&c), a.clone()),
        (r.minus(&c), a.clone()),
    ];
    if !rational_first {
        return frames;
    }
    let ua = number::unit(a);
    let ur = number::unit(r);
    let uc = ua.cross(&ur);
    // Rational turns avoid an extra square root in chart proofs.
    let three = Real::from_rational(Rational::from_ratio(3, 5).unwrap());
    let four = Real::from_rational(Rational::from_ratio(4, 5).unwrap());
    let a3 = ua.scaled(&three);
    let r3 = ur.scaled(&three);
    let r4 = ur.scaled(&four);
    let c4 = uc.scaled(&four);
    let rational = [
        (a3.plus(&r4), uc.clone()),
        (a3.minus(&r4), uc),
        (a3.plus(&c4), ur.clone()),
        (a3.minus(&c4), ur),
        (r3.plus(&c4), ua.clone()),
        (r3.minus(&c4), ua),
    ];
    // Three-arc corners benefit; other trims keep their earlier chart cuts.
    drop(frames.splice(2..2, rational));
    frames
}

/// Trimmed faces of `bodies` on the selected surfaces, by surface.
fn trimmed(
    model: &Model,
    bodies: &[BodyId],
    select: impl Fn(SurfaceId) -> bool,
) -> BTreeMap<SurfaceId, Vec<FaceId>> {
    let mut faces = BTreeMap::<SurfaceId, Vec<FaceId>>::new();
    for &body in bodies {
        for face in model.faces_of_body(body) {
            let definition = model.face(face);
            if !definition.loops.is_empty() && select(definition.surface) {
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
fn rewrite(model: &Model, bodies: &[BodyId], frames: BTreeMap<SurfaceId, Surface>) -> Reframed {
    let mut moved = model.clone();
    let mut faces = BTreeMap::new();
    let mut trims = trimmed(model, bodies, |surface| frames.contains_key(&surface));
    for (surface, framed) in frames {
        let frame = moved.geom_mut().add_surface(framed);
        for face in trims.remove(&surface).into_iter().flatten() {
            let source = model.face(face).clone();
            // The implicit surface is unchanged, so every source incidence holds.
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
                            origin: FactOrigin::Construction("step.surface-frame"),
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{RPoint3, RVector3};

    #[test]
    fn rational_turns_keep_unit_perpendicular_frames_for_scaled_source_axes() {
        for scale in [1, 3, -2] {
            let sphere = Sphere {
                center: RPoint3::new(0.into(), 0.into(), 0.into()).into(),
                radius: Real::one(),
                axis: RVector3::new(0.into(), 0.into(), (2 * scale).into()).into(),
                ref_dir: RVector3::new(scale.into(), 0.into(), 0.into()).into(),
            };
            let legacy_first = candidates(&sphere, false);
            assert_eq!(legacy_first.len(), 8);
            assert_eq!(
                legacy_first[2]
                    .0
                    .minus(&sphere.axis.plus(&sphere.ref_dir))
                    .norm_squared()
                    .as_rational(),
                Some(Rational::zero())
            );
            let frames = candidates(&sphere, true);
            assert_eq!(frames.len(), 14);
            for (axis, reference) in frames.into_iter().skip(2).take(6) {
                assert_eq!(axis.norm_squared().as_rational(), Some(Rational::one()));
                assert_eq!(
                    reference.norm_squared().as_rational(),
                    Some(Rational::one())
                );
                assert_eq!(axis.dot(&reference).as_rational(), Some(Rational::zero()));
                for vector in [&axis, &reference] {
                    assert!(
                        [&vector.x, &vector.y, &vector.z]
                            .iter()
                            .all(|c| c.as_rational().is_some())
                    );
                }
            }
        }
    }
}
