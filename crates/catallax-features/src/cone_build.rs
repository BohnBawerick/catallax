//! Complete-section recognition and shared analytic conical topology.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Refusal};
use catallax_shapes::{Curve, PointDef, Surface};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Loop, Lump, Model, Sense, Shell, Vertex};

use crate::cone_skin::ConeSkin;

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.cone-domain", detail)
}

pub(crate) fn skin(
    model: &mut Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<Option<ConeSkin>, Refusal> {
    let faces = model.faces_of_body(body);
    if !faces.iter().any(|f| {
        matches!(
            model.geom().surface(model.face(*f).surface),
            Surface::Cone(_)
        )
    }) {
        return Ok(None);
    }
    if model.body(body).lumps.len() != 1 || model.lump(model.body(body).lumps[0]).shells.len() != 1
    {
        return Err(bad(
            "multiple conical cells need their complete distance-envelope arrangement",
        ));
    }
    let mut edges = BTreeMap::new();
    let mut surfaces = Vec::new();
    let mut apex_count = 0;
    for (i, face) in faces.iter().enumerate() {
        let f = model.face(*face);
        if f.sense != Sense::Same {
            return Err(bad(
                "concave conical faces need their edge-distance envelope arrangement",
            ));
        }
        for lp in &f.loops {
            match model.loop_(*lp) {
                Loop::Coedges(cs) => {
                    let [coedge] = cs.as_slice() else {
                        return Err(bad(
                            "partial conical sections need a finite trim arrangement",
                        ));
                    };
                    let edge = model.coedge(*coedge).edge;
                    let e = model.edge(edge);
                    if e.start.is_some()
                        || e.end.is_some()
                        || !matches!(
                            catallax_intersect::curve_form(model.geom(), e.curve, dec)?,
                            catallax_intersect::CurveForm::Circle { .. }
                                | catallax_intersect::CurveForm::Ellipse { .. }
                        )
                    {
                        return Err(bad(
                            "conical distance construction requires complete closed conic trims",
                        ));
                    }
                    edges.entry(edge).or_insert_with(Vec::new).push(i);
                }
                Loop::Vertex(v) if matches!(model.geom().surface(f.surface), Surface::Cone(_)) => {
                    let Surface::Cone(c) = model.geom().surface(f.surface) else {
                        unreachable!()
                    };
                    let point = model.geom().position(model.vertex(*v).point);
                    if crate::profile::sign(&point.minus(&c.apex).norm_squared(), dec)?
                        != catallax_numbers::Sign::Zero
                    {
                        return Err(bad("a singular trim is not the retained cone apex"));
                    }
                    apex_count += 1;
                }
                _ => return Err(bad("only the cone tip can be a singular conical trim")),
            }
        }
        surfaces.push((
            model.geom().surface(f.surface).clone(),
            model.face_origins(*face).to_vec(),
        ));
    }
    let Some(skin) = ConeSkin::new(&surfaces, dec)? else {
        return Err(bad(
            "the conical support inventory is not a bounded convex cell",
        ));
    };
    if skin.faces.len() != faces.len()
        || skin.boundaries.len() != edges.len()
        || apex_count != usize::from(skin.apex.is_some())
    {
        return Err(bad(
            "the complete conical sections do not match the input trim inventory",
        ));
    }
    for boundary in &skin.boundaries {
        let curve = model.geom_mut().add_curve(boundary.curve.clone());
        let named: Vec<_> = boundary
            .faces
            .iter()
            .map(|i| {
                surfaces
                    .iter()
                    .position(|(s, _)| s.same_definition(&skin.faces[*i].surface))
                    .expect("retained source face")
            })
            .collect();
        let mut matched = false;
        for (edge, incident) in &edges {
            if incident.len() == 2
                && named.iter().all(|i| incident.contains(i))
                && catallax_intersect::same_curve(model.geom(), curve, model.edge(*edge).curve, dec)
                    .is_coincident()
            {
                matched = true;
                break;
            }
        }
        if !matched {
            return Err(bad(
                "a conical source trim differs from its full analytic section",
            ));
        }
    }
    Ok(Some(skin))
}

pub(crate) fn build(model: &mut Model, skin: &ConeSkin) -> BodyId {
    let mut paths = vec![Vec::new(); skin.faces.len()];
    for boundary in &skin.boundaries {
        let curve = model.geom_mut().add_curve(boundary.curve.clone());
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        paths[boundary.faces[0]].push(vec![(edge, Sense::Reversed)]);
        paths[boundary.faces[1]].push(vec![(edge, Sense::Same)]);
    }
    let apex = skin.apex.as_ref().map(|(p, face)| {
        let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
        (model.add_vertex(Vertex { point }), *face)
    });
    let faces = skin
        .faces
        .iter()
        .zip(paths)
        .enumerate()
        .map(|(i, (f, paths))| {
            crate::extrude::face(
                model,
                f.surface.clone(),
                Sense::Same,
                paths,
                apex.filter(|a| a.1 == i).map(|a| a.0),
                f.origins.clone(),
            )
        })
        .collect();
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

pub(crate) fn sweep(
    model: &mut Model,
    sweep: &crate::cone_sweep::ConeSweep,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let mut rings = Vec::new();
    let mut vertices = Vec::new();
    for (i, point) in sweep.vertices.iter().enumerate() {
        if crate::profile::sign(&point[0], dec)? == catallax_numbers::Sign::Zero {
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit(sweep.circle(i).center));
            vertices.push(Some(model.add_vertex(Vertex { point })));
            rings.push(None);
        } else {
            vertices.push(None);
            let curve = model.geom_mut().add_curve(Curve::Circle(sweep.circle(i)));
            rings.push(Some(model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            })));
        }
    }
    let mut faces = Vec::new();
    for f in &sweep.faces {
        let mut paths = Vec::new();
        let mut apex = None;
        for (i, index) in f.ends.iter().enumerate() {
            if let Some(edge) = rings[*index] {
                paths.push(vec![(
                    edge,
                    if i == 0 { Sense::Same } else { Sense::Reversed },
                )]);
            } else if !matches!(f.surface, Surface::Plane(_)) {
                apex = vertices[*index];
            }
        }
        faces.push(crate::extrude::face(
            model,
            f.surface.clone(),
            if f.reversed {
                Sense::Reversed
            } else {
                Sense::Same
            },
            paths,
            apex,
            f.origins.clone(),
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    }))
}
