//! Finite normal images and projection sectors of selected convex rotational faces.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Circle, Surface, Torus};

use crate::{
    cone_skin::{ConeSkin, reduced},
    cone_sweep::ConeSweep,
    cylinder_skin::CylinderSkin,
    profile::sign,
};

/// A circular projection sector between two selected finite normal images.
#[derive(Clone, Debug)]
pub struct RoundJoin {
    /// Exact ring-torus carrier.
    pub surface: Surface,
    /// The selected input edge, not a boundary rim.
    pub origin: Origin,
    /// Removed rim supports and reversed material senses. The sector lies outside both images.
    pub cuts: [(Surface, bool); 2],
    /// Far circular contacts; true is a positive use in the removed rim.
    pub contacts: [(Circle, bool); 2],
}

/// Selected finite images, including only interior-edge projection sectors.
#[derive(Clone, Debug)]
pub struct RotationalThicken {
    /// Face images. Unnamed rim patches are internal construction boundaries.
    pub sweeps: Vec<ConeSweep>,
    /// Front sectors. Convex source skins have no back edge sectors.
    pub joins: Vec<RoundJoin>,
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.unsupported-envelope", detail)
}

/// Whether a generated rim lies between two selected source faces.
pub fn interior_rim(origin: &Origin, selected: &[&[Origin]]) -> bool {
    matches!(&origin.role, Role::Rim { edge } if edge.between().iter().all(|o| selected.iter().any(|face| face.contains(o))))
}

impl RotationalThicken {
    /// Normal images of complete conical bands and perpendicular planar caps.
    pub fn conical(
        skin: &ConeSkin,
        selected: &[bool],
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut sweeps = Vec::new();
        for (i, yes) in selected.iter().enumerate() {
            if *yes {
                if let Surface::Cylinder(cylinder) = &skin.faces[i].surface {
                    let mut rims = Vec::new();
                    for boundary in &skin.boundaries {
                        if !boundary.faces.contains(&i) {
                            continue;
                        }
                        let catallax_shapes::Curve::Circle(circle) = &boundary.curve else {
                            return Err(bad(
                                "a cylindrical selected neighbour requires complete circular trims",
                            ));
                        };
                        let neighbour = boundary
                            .faces
                            .iter()
                            .copied()
                            .find(|j| *j != i)
                            .ok_or_else(|| bad("a cylindrical rim has no neighbour"))?;
                        rims.push((circle.center.clone(), skin.faces[neighbour].origins.clone()));
                    }
                    let [a, b] = rims.as_mut_slice() else {
                        return Err(bad(
                            "a cylindrical selected neighbour requires two circular rims",
                        ));
                    };
                    if sign(&b.0.minus(&a.0).dot(&skin.axis), dec)? == Sign::Negative {
                        std::mem::swap(a, b);
                    }
                    let mut cylinder = cylinder.clone();
                    cylinder.axis = skin.axis.clone();
                    let band = CylinderSkin {
                        cylinder,
                        ends: [a.0.clone(), b.0.clone()],
                        origins: [a.1.clone(), b.1.clone(), skin.faces[i].origins.clone()],
                    };
                    let mut image = Self::cylindrical(
                        &band,
                        &[false, false, true],
                        thickness,
                        side,
                        step,
                        dec,
                    )?;
                    sweeps.append(&mut image.sweeps);
                } else {
                    sweeps.push(skin.face_sweep(i, thickness, side, step, dec)?);
                }
            }
        }
        let origins: Vec<_> = skin
            .faces
            .iter()
            .zip(selected)
            .filter_map(|(f, yes)| yes.then_some(f.origins.as_slice()))
            .collect();
        Self::new(sweeps, &origins, thickness, side, step, dec)
    }

    /// Normal images of complete cylindrical walls and perpendicular planar caps.
    pub fn cylindrical(
        skin: &CylinderSkin,
        selected: &[bool],
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let length = reduced(skin.cylinder.axis.norm_squared().sqrt(), dec)?;
        let axis = skin.cylinder.axis.scaled(&(&Real::one() / &length));
        let mut sweeps = Vec::new();
        for (i, yes) in selected.iter().enumerate() {
            if !yes {
                continue;
            }
            let band = skin.thicken(i, thickness, side, step, dec)?;
            let low = reduced(
                band.outer.ends[0].minus(&skin.cylinder.origin).dot(&axis),
                dec,
            )?;
            let high = reduced(
                band.outer.ends[1].minus(&skin.cylinder.origin).dot(&axis),
                dec,
            )?;
            let outer = reduced(band.outer.cylinder.radius(), dec)?;
            let (inner, inner_names) = match &band.inner {
                Some(inner) => (
                    reduced(inner.cylinder.radius(), dec)?,
                    inner.origins[2].clone(),
                ),
                None => (Real::zero(), Vec::new()),
            };
            sweeps.push(ConeSweep::new(
                skin.cylinder.origin.clone(),
                axis.clone(),
                skin.cylinder.ref_dir.clone(),
                vec![
                    [inner.clone(), low.clone()],
                    [outer.clone(), low],
                    [outer, high.clone()],
                    [inner, high],
                ],
                vec![
                    band.outer.origins[0].clone(),
                    band.outer.origins[2].clone(),
                    band.outer.origins[1].clone(),
                    inner_names,
                ],
                dec,
            )?);
        }
        let origins: Vec<_> = skin
            .origins
            .iter()
            .zip(selected)
            .filter_map(|(f, yes)| yes.then_some(f.as_slice()))
            .collect();
        Self::new(sweeps, &origins, thickness, side, step, dec)
    }

    fn new(
        mut sweeps: Vec<ConeSweep>,
        selected: &[&[Origin]],
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut rims: Vec<(EdgeName, Vec<(usize, usize)>)> = Vec::new();
        for (i, sweep) in sweeps.iter().enumerate() {
            for (j, face) in sweep.faces.iter().enumerate() {
                let Some(origin) = face.origins.iter().find(|o| interior_rim(o, selected)) else {
                    continue;
                };
                let Role::Rim { edge } = &origin.role else {
                    unreachable!()
                };
                if let Some((_, uses)) = rims.iter_mut().find(|(name, _)| name == edge) {
                    uses.push((i, j));
                } else {
                    rims.push((edge.clone(), vec![(i, j)]));
                }
            }
        }
        let mut joins = Vec::new();
        for (edge, uses) in rims {
            if matches!(side, ThickenSide::Front) {
                let [a, b] = uses.as_slice() else {
                    return Err(bad(
                        "a circular selected edge needs exactly two finite images",
                    ));
                };
                let cuts = [a, b].map(|(i, j)| {
                    let face = &sweeps[*i].faces[*j];
                    (face.surface.clone(), face.reversed)
                });
                let contact = |(s, f): (usize, usize)| -> Result<(Circle, bool, Circle), Refusal> {
                    let sweep = &sweeps[s];
                    let face = &sweep.faces[f];
                    for (end, vertex) in face.ends.iter().enumerate() {
                        if sweep.faces.iter().any(|other| {
                            other.ends.contains(vertex)
                                && other
                                    .origins
                                    .iter()
                                    .any(|o| matches!(o.role, Role::Offset { .. }))
                        }) {
                            return Ok((
                                sweep.circle(*vertex),
                                end == 0,
                                sweep.circle(face.ends[1 - end]),
                            ));
                        }
                    }
                    Err(bad(
                        "a circular selected edge has no far normal-image contact",
                    ))
                };
                let (ca, sa, base) = contact(*a)?;
                let (cb, sb, other) = contact(*b)?;
                if sign(&base.center.minus(&other.center).norm_squared(), dec)? != Sign::Zero
                    || sign(&(&base.radius_sq - &other.radius_sq), dec)? != Sign::Zero
                {
                    return Err(bad("selected circular rims do not share their source edge"));
                }
                let radius = reduced(base.radius_sq.sqrt(), dec)?;
                if sign(&(&radius - thickness), dec)? != Sign::Positive {
                    return Err(Refusal::new(
                        "offset.spindle-carrier",
                        "the circular projection sector requires a ring-torus carrier",
                    ));
                }
                for (surface, _) in &cuts {
                    if let Surface::Cone(c) = surface {
                        let axial = reduced(base.center.minus(&c.apex).dot(&c.axis).square(), dec)?;
                        if sign(&(&axial - &thickness.square()), dec)? != Sign::Positive {
                            return Err(bad(
                                "a circular sector crossing its rim-cone apex plane needs a directed meridian proof",
                            ));
                        }
                    }
                }
                joins.push(RoundJoin {
                    surface: Surface::Torus(Torus {
                        center: base.center,
                        axis: base.normal,
                        major: radius,
                        minor: thickness.clone(),
                        ref_dir: base.ref_dir,
                    }),
                    origin: Origin {
                        step: step.clone(),
                        role: Role::Round { edge },
                    },
                    cuts,
                    contacts: [(ca, sa), (cb, sb)],
                });
            }
            for (i, j) in uses {
                sweeps[i].faces[j]
                    .origins
                    .retain(|o| !interior_rim(o, selected));
            }
        }
        Ok(Self { sweeps, joins })
    }
}

impl RotationalThicken {
    pub(crate) fn build(
        &self,
        model: &mut catallax_topology::Model,
        side: ThickenSide,
        dec: &mut Decider,
    ) -> Result<catallax_topology::BodyId, Refusal> {
        use catallax_shapes::Curve;
        use catallax_topology::{Edge, EdgeId, Sense};
        if matches!(side, ThickenSide::Back) || self.joins.is_empty() {
            let mut result = catallax_topology::empty(model);
            for sweep in &self.sweeps {
                let piece = crate::cone_build::sweep(model, sweep, dec)?;
                result = catallax_boolean::boolean(
                    model,
                    result,
                    piece,
                    catallax_boolean::BoolOp::Unite,
                    dec,
                )
                .map_err(|e| e.refusal)?
                .body;
            }
            return Ok(result);
        }
        fn ring(
            model: &mut catallax_topology::Model,
            rings: &mut Vec<(Circle, EdgeId)>,
            circle: Circle,
            dec: &mut Decider,
        ) -> Result<EdgeId, Refusal> {
            for (known, edge) in rings.iter() {
                if sign(&known.center.minus(&circle.center).norm_squared(), dec)? == Sign::Zero
                    && sign(&(&known.radius_sq - &circle.radius_sq), dec)? == Sign::Zero
                    && sign(&known.normal.cross(&circle.normal).norm_squared(), dec)? == Sign::Zero
                {
                    if sign(&known.normal.dot(&circle.normal), dec)? != Sign::Positive {
                        return Err(bad(
                            "circular normal-image contacts have opposite parameter axes",
                        ));
                    }
                    return Ok(*edge);
                }
            }
            let curve = model.geom_mut().add_curve(Curve::Circle(circle.clone()));
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            });
            rings.push((circle, edge));
            Ok(edge)
        }
        let mut rings = Vec::new();
        let mut faces = Vec::new();
        for sweep in &self.sweeps {
            for patch in &sweep.faces {
                if patch.origins.is_empty() {
                    continue;
                }
                let mut paths = Vec::new();
                for (i, end) in patch.ends.iter().enumerate() {
                    if sign(&sweep.vertices[*end][0], dec)? == Sign::Zero {
                        if !matches!(patch.surface, Surface::Plane(_)) {
                            return Err(bad(
                                "a sharp rotational image reaching an axis needs its singular boundary proof",
                            ));
                        }
                        continue;
                    }
                    let edge = ring(model, &mut rings, sweep.circle(*end), dec)?;
                    paths.push(vec![(
                        edge,
                        if i == 0 { Sense::Same } else { Sense::Reversed },
                    )]);
                }
                faces.push(crate::extrude::face(
                    model,
                    patch.surface.clone(),
                    if patch.reversed {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    },
                    paths,
                    None,
                    patch.origins.clone(),
                ));
            }
        }
        for join in &self.joins {
            let mut paths = Vec::new();
            for (circle, first) in &join.contacts {
                let edge = ring(model, &mut rings, circle.clone(), dec)?;
                paths.push(vec![(
                    edge,
                    if *first { Sense::Same } else { Sense::Reversed },
                )]);
            }
            faces.push(crate::extrude::face(
                model,
                join.surface.clone(),
                Sense::Same,
                paths,
                None,
                vec![join.origin.clone()],
            ));
        }
        assemble(model, faces, &self.joins[0].origin.step)
    }
}

pub(crate) fn assemble(
    model: &mut catallax_topology::Model,
    faces: Vec<catallax_topology::FaceId>,
    step: &StepId,
) -> Result<catallax_topology::BodyId, Refusal> {
    use catallax_topology::{Body, BodyKind, Loop, Lump, Shell};
    use std::collections::BTreeMap;
    let mut owners: BTreeMap<_, Vec<usize>> = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        for lp in &model.face(*face).loops {
            if let Loop::Coedges(coedges) = model.loop_(*lp) {
                for coedge in coedges {
                    owners
                        .entry(model.coedge(*coedge).edge)
                        .or_default()
                        .push(i);
                }
            }
        }
    }
    let mut neighbours = vec![Vec::new(); faces.len()];
    for uses in owners.values() {
        for i in uses {
            neighbours[*i].extend(uses.iter().copied().filter(|j| j != i));
        }
    }
    let mut seen = vec![false; faces.len()];
    let mut shells = Vec::new();
    for i in 0..faces.len() {
        if seen[i] {
            continue;
        }
        let mut stack = vec![i];
        let mut component = Vec::new();
        seen[i] = true;
        while let Some(j) = stack.pop() {
            component.push(faces[j]);
            for k in &neighbours[j] {
                if !seen[*k] {
                    seen[*k] = true;
                    stack.push(*k);
                }
            }
        }
        let outer = component.iter().any(|face| {
            model.face_origins(*face).iter().any(|o| {
                o.step == *step
                    && matches!(
                        o.role,
                        Role::Offset { .. } | Role::Round { .. } | Role::RoundCorner { .. }
                    )
            })
        });
        shells.push((outer, model.add_shell(Shell { faces: component })));
    }
    shells.sort_by_key(|(outer, _)| !outer);
    let outer_count = shells.iter().filter(|(outer, _)| *outer).count();
    let lumps = if outer_count == 1 {
        vec![model.add_lump(Lump {
            shells: shells.into_iter().map(|(_, shell)| shell).collect(),
        })]
    } else {
        if outer_count != shells.len() {
            return Err(bad(
                "disconnected front normal fans need a proved cavity-to-component assignment",
            ));
        }
        shells
            .into_iter()
            .map(|(_, shell)| {
                model.add_lump(Lump {
                    shells: vec![shell],
                })
            })
            .collect()
    };
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
        marks: Default::default(),
    }))
}
