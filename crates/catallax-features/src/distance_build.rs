//! Distance feature construction from immutable input boundaries.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Origin, StepId, ThickenSide};
use catallax_shapes::Surface;
use catallax_topology::{Body, BodyId, BodyKind, Face, Loop, Lump, Model, Sense, Shell};

use crate::{
    polyhedron::{ConvexCell, PlaneFace},
    profile::sign,
    radial::RadialBand,
};

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.unsupported-envelope", detail)
}
fn positive(thickness: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    if sign(thickness, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "distance.non-positive-thickness",
            "thickness must be positive",
        ));
    }
    Ok(())
}
fn radial(model: &Model, body: BodyId) -> Option<(Surface, Vec<Origin>)> {
    let faces = model.faces_of_body(body);
    let [face] = faces.as_slice() else {
        return None;
    };
    let f = model.face(*face);
    let surface = model.geom().surface(f.surface);
    if f.sense != Sense::Same
        || !f.loops.is_empty()
        || !matches!(surface, Surface::Sphere(_) | Surface::Torus(_))
    {
        return None;
    }
    Some((surface.clone(), model.face_origins(*face).to_vec()))
}
fn closed(model: &mut Model, band: RadialBand) -> BodyId {
    let mut shells = Vec::new();
    for ((surface, origins), sense) in std::iter::once((band.outer, Sense::Same))
        .chain(band.inner.map(|inner| (inner, Sense::Reversed)))
    {
        let surface = model.geom_mut().add_surface(surface);
        let face = model.add_named_face(
            Face {
                surface,
                sense,
                loops: Vec::new(),
            },
            origins,
        );
        shells.push(model.add_shell(Shell { faces: vec![face] }));
    }
    let lump = model.add_lump(Lump { shells });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    })
}
fn cylindrical(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<Option<crate::cylinder_skin::CylinderSkin>, Refusal> {
    let mut surfaces = Vec::new();
    for face in model.faces_of_body(body) {
        let f = model.face(face);
        if f.sense != Sense::Same {
            return Ok(None);
        }
        for lp in &f.loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                return Ok(None);
            };
            let [coedge] = cs.as_slice() else {
                return Ok(None);
            };
            let edge = model.edge(model.coedge(*coedge).edge);
            if edge.start.is_some()
                || edge.end.is_some()
                || !matches!(
                    model.geom().curve(edge.curve),
                    catallax_shapes::Curve::Circle(_)
                )
            {
                return Ok(None);
            }
        }
        surfaces.push((
            model.geom().surface(f.surface).clone(),
            model.face_origins(face).to_vec(),
        ));
    }
    crate::cylinder_skin::CylinderSkin::new(&surfaces, dec)
}

fn cylinder(model: &mut Model, skin: &crate::cylinder_skin::CylinderSkin) -> BodyId {
    use catallax_shapes::{Circle, Curve};
    use catallax_topology::Edge;
    let rings = skin.ends.each_ref().map(|center| {
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: center.clone(),
            normal: skin.cylinder.axis.clone(),
            radius_sq: skin.cylinder.radius_sq.clone(),
            ref_dir: skin.cylinder.ref_dir.clone(),
        }));
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        })
    });
    let surfaces = skin.surfaces();
    let faces = vec![
        crate::extrude::face(
            model,
            surfaces[0].clone(),
            Sense::Same,
            vec![vec![(rings[0], Sense::Reversed)]],
            None,
            skin.origins[0].clone(),
        ),
        crate::extrude::face(
            model,
            surfaces[1].clone(),
            Sense::Same,
            vec![vec![(rings[1], Sense::Same)]],
            None,
            skin.origins[1].clone(),
        ),
        crate::extrude::face(
            model,
            surfaces[2].clone(),
            Sense::Same,
            vec![
                vec![(rings[0], Sense::Same)],
                vec![(rings[1], Sense::Reversed)],
            ],
            None,
            skin.origins[2].clone(),
        ),
    ];
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    })
}
fn convex(model: &Model, body: BodyId, dec: &mut Decider) -> Result<ConvexCell, Refusal> {
    let mut planes = Vec::new();
    let mut points = Vec::new();
    for face in model.faces_of_body(body) {
        let f = model.face(face);
        let Surface::Plane(plane) = model.geom().surface(f.surface) else {
            return Err(bad(
                "curved trimmed skins require their offset and edge-envelope arrangement",
            ));
        };
        let mut plane = plane.clone();
        if f.sense == Sense::Reversed {
            plane.normal = plane.normal.scaled(&Real::from_i64(-1));
        }
        planes.push(PlaneFace {
            plane,
            origins: model.face_origins(face).to_vec(),
        });
        for lp in &f.loops {
            let Loop::Coedges(coedges) = model.loop_(*lp) else {
                return Err(bad("a convex planar face needs a polygonal boundary"));
            };
            for coedge in coedges {
                let edge = model.edge(model.coedge(*coedge).edge);
                if !matches!(
                    model.geom().curve(edge.curve),
                    catallax_shapes::Curve::Line(_)
                ) {
                    return Err(bad("a convex planar skin requires straight polygon edges"));
                }
                let Some(vertex) = edge.start else {
                    return Err(bad("a polygon edge requires finite endpoints"));
                };
                points.push(model.geom().position(model.vertex(vertex).point));
            }
        }
    }
    let cell = ConvexCell::new(planes, dec)?;
    if cell.is_empty() || !cell.contains_vertices(&points, dec)? {
        return Err(bad(
            "nonconvex planar skins need the complete trimmed-face envelope arrangement",
        ));
    }
    Ok(cell)
}
fn closed_selection(origins: &[Origin], refs: &[FaceRef]) -> Result<bool, Refusal> {
    for reference in refs {
        if !reference.touching.is_empty() || !origins.iter().any(|o| reference.face.matches(o)) {
            return Err(Refusal::new(
                "reference.lost-face",
                format!("no closed input face resolves {reference:?}"),
            ));
        }
    }
    Ok(!refs.is_empty())
}

pub(crate) fn offset(
    model: &mut Model,
    body: BodyId,
    distance: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let direction = sign(distance, dec)?;
    if direction == Sign::Zero || model.body(body).kind == BodyKind::Empty {
        return Ok(body);
    }
    if let Some((surface, origins)) = radial(model, body) {
        return Ok(match surface.parallel_support(distance, dec)? {
            Some(surface) => closed(
                model,
                RadialBand {
                    outer: (surface, origins),
                    inner: None,
                },
            ),
            None => catallax_topology::empty(model),
        });
    }
    if direction == Sign::Positive
        && let Some(skin) = crate::circular_build::skin(model, body, dec)?
    {
        return crate::circular_build::outward(model, &skin, distance, step, dec);
    }
    if let Some(skin) = cylindrical(model, body, dec)? {
        return Ok(match skin.eroded(&(-distance), dec)? {
            Some(skin) => cylinder(model, &skin),
            None => catallax_topology::empty(model),
        });
    }
    if let Some(skin) = crate::cone_build::skin(model, body, dec)? {
        if direction == Sign::Positive {
            return crate::cone_rounded::build(model, &skin, distance, step, dec);
        }
        return Ok(match skin.eroded(&(-distance), dec)? {
            Some(skin) => crate::cone_build::build(model, &skin),
            None => catallax_topology::empty(model),
        });
    }
    let cell = convex(model, body, dec)?;
    if direction == Sign::Positive {
        let marks = model.body(body).marks.clone();
        return crate::rounded_polyhedron::build(model, &cell, distance, step, marks, dec);
    }
    Ok(cell.eroded(&(-distance), dec)?.build(model))
}

pub(crate) fn shell(
    model: &mut Model,
    body: BodyId,
    open: &[FaceRef],
    thickness: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    positive(thickness, dec)?;
    if model.body(body).kind == BodyKind::Empty {
        if open.is_empty() {
            return Ok(body);
        }
        return Err(Refusal::new(
            "reference.lost-face",
            "an empty body has no open face",
        ));
    }
    if let Some((surface, origins)) = radial(model, body) {
        if closed_selection(&origins, open)? {
            return Ok(catallax_topology::empty(model));
        }
        return Ok(closed(
            model,
            RadialBand::shell(&surface, &origins, thickness, step, dec)?,
        ));
    }
    if let Some(skin) = cylindrical(model, body, dec)? {
        let selected = skin.select(open)?;
        if selected.iter().all(|v| *v) {
            return Ok(catallax_topology::empty(model));
        }
        let Some(cavity) = skin.cavity(&selected, thickness, Some(step), dec)? else {
            return Ok(body);
        };
        let cavity = cylinder(model, &cavity);
        return Ok(catallax_boolean::boolean(
            model,
            body,
            cavity,
            catallax_boolean::BoolOp::Subtract,
            dec,
        )
        .map_err(|e| e.refusal)?
        .body);
    }
    if let Some(skin) = crate::cone_build::skin(model, body, dec)? {
        let selected = skin.select(open)?;
        if selected.iter().all(|v| *v) {
            return Ok(catallax_topology::empty(model));
        }
        let cavity = match skin.cavity(&selected, thickness, step, dec) {
            Ok(Some(cavity)) => cavity,
            Ok(None) => return Ok(body),
            Err(e) => return Err(e),
        };
        let cavity = crate::cone_build::build(model, &cavity);
        return Ok(catallax_boolean::boolean(
            model,
            body,
            cavity,
            catallax_boolean::BoolOp::Subtract,
            dec,
        )
        .map_err(|e| e.refusal)?
        .body);
    }
    if let Some(skin) = crate::circular_build::skin(model, body, dec)? {
        let selected = skin.select(open)?;
        if selected.iter().all(|v| *v) {
            return Ok(catallax_topology::empty(model));
        }
        let Some(cavity) = skin.cavity(&selected, thickness, step, dec)? else {
            return Ok(body);
        };
        let cavity = crate::cone_build::sections(model, &cavity.faces, &cavity.boundaries, None);
        return Ok(catallax_boolean::boolean(
            model,
            body,
            cavity,
            catallax_boolean::BoolOp::Subtract,
            dec,
        )
        .map_err(|e| e.refusal)?
        .body);
    }
    let cell = convex(model, body, dec)?;
    let selected = cell.select(open, dec)?;
    if selected.iter().all(|v| *v) {
        return Ok(catallax_topology::empty(model));
    }
    let cavity = match cell.cavity(&selected, thickness, step, dec) {
        Ok(cavity) => cavity,
        Err(e) if e.code == "distance.planar-domain" => {
            let mut walls = catallax_topology::empty(model);
            for (i, opened) in selected.iter().enumerate() {
                if *opened || cell.polygons[i].is_empty() {
                    continue;
                }
                let tube = crate::rounded_polyhedron::plate(model, &cell, i, thickness, step, dec)?;
                walls = catallax_boolean::boolean(
                    model,
                    walls,
                    tube,
                    catallax_boolean::BoolOp::Unite,
                    dec,
                )
                .map_err(|e| e.refusal)?
                .body;
            }
            return Ok(catallax_boolean::boolean(
                model,
                body,
                walls,
                catallax_boolean::BoolOp::Intersect,
                dec,
            )
            .map_err(|e| e.refusal)?
            .body);
        }
        Err(e) => return Err(e),
    };
    if cavity.is_empty() {
        return Ok(body);
    }
    let cavity = cavity.build(model);
    Ok(
        catallax_boolean::boolean(model, body, cavity, catallax_boolean::BoolOp::Subtract, dec)
            .map_err(|e| e.refusal)?
            .body,
    )
}

pub(crate) fn thicken(
    model: &mut Model,
    body: BodyId,
    refs: &[FaceRef],
    thickness: &Real,
    side: ThickenSide,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    positive(thickness, dec)?;
    if refs.is_empty() || model.body(body).kind == BodyKind::Empty {
        return Err(Refusal::new(
            "reference.lost-face",
            "thicken needs a surviving selected face",
        ));
    }
    if let Some((surface, origins)) = radial(model, body) {
        if !closed_selection(&origins, refs)? {
            return Err(Refusal::new(
                "reference.lost-face",
                "thicken needs a selected face",
            ));
        }
        return Ok(closed(
            model,
            RadialBand::thicken(&surface, &origins, thickness, side, step, dec)?,
        ));
    }
    if let Some(skin) = cylindrical(model, body, dec)? {
        let selected = skin.select(refs)?;
        let indices: Vec<_> = selected
            .iter()
            .enumerate()
            .filter_map(|(i, v)| v.then_some(i))
            .collect();
        let [index] = indices.as_slice() else {
            return Err(bad(
                "selected cylinder joins require their projection normal sectors",
            ));
        };
        let band = skin.thicken(*index, thickness, side, step, dec)?;
        let outer = cylinder(model, &band.outer);
        return if let Some(inner) = band.inner {
            let inner = cylinder(model, &inner);
            Ok(catallax_boolean::boolean(
                model,
                outer,
                inner,
                catallax_boolean::BoolOp::Subtract,
                dec,
            )
            .map_err(|e| e.refusal)?
            .body)
        } else {
            Ok(outer)
        };
    }
    if let Some(skin) = crate::cone_build::skin(model, body, dec)? {
        let selected = skin.select(refs)?;
        let indices: Vec<_> = selected
            .iter()
            .enumerate()
            .filter_map(|(i, v)| v.then_some(i))
            .collect();
        let [index] = indices.as_slice() else {
            return Err(bad(
                "sharp selected conical joins require their projection normal sectors",
            ));
        };
        let sweep = skin.face_sweep(*index, thickness, side, step, dec)?;
        return crate::cone_build::sweep(model, &sweep, dec);
    }
    if let Some(skin) = crate::circular_build::skin(model, body, dec)? {
        let selected = skin.select(refs)?;
        let indices: Vec<_> = selected
            .iter()
            .enumerate()
            .filter_map(|(i, v)| v.then_some(i))
            .collect();
        let [index] = indices.as_slice() else {
            return Err(bad(
                "sharp selected circular joins require their projection normal sectors",
            ));
        };
        let sweep = skin.face_sweep(*index, thickness, side, step, dec)?;
        return crate::cone_build::sweep(model, &sweep, dec);
    }
    let cell = convex(model, body, dec)?;
    let selected = cell.select(refs, dec)?;
    let indices: Vec<_> = selected
        .iter()
        .enumerate()
        .filter_map(|(i, v)| v.then_some(i))
        .collect();
    let [index] = indices.as_slice() else {
        return Err(bad(
            "sharp selected joins require their projection normal sectors",
        ));
    };
    Ok(cell
        .face_sweep(*index, thickness, side, step, dec)?
        .build(model))
}

pub(crate) fn inherit_marks(model: &mut Model, source: BodyId, result: BodyId) -> BodyId {
    if source == result
        || model.body(source).marks.is_empty()
        || model.body(source).marks == model.body(result).marks
    {
        return result;
    }
    let mut body = model.body(result).clone();
    body.marks.extend(&model.body(source).marks);
    model.add_body(body)
}
