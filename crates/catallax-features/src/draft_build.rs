//! Immutable topology construction from drafted support arrangements.

use crate::draft::{
    DraftDefinition,
    cell::{PlanarCell, PlanePatch},
    round::RoundCell,
    sign,
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Origin, StepId};
use catallax_shapes::{Circle, Curve, Line, Point3, PointDef, Surface};
use catallax_topology::{
    Body, BodyId, BodyKind, Edge, FaceId, Loop, Lump, Model, Sense, Shell, Vertex,
};
use std::collections::{BTreeMap, BTreeSet};

/// Removed face history, separate from selectable surviving origins.
#[derive(Clone, Debug)]
pub struct ConsumedFace {
    /// Draft operation that consumed the source face.
    pub step: StepId,
    /// Complete original origin set of the removed face.
    pub source: Vec<Origin>,
    /// Origins of replacing drafted faces, never a new draft origin.
    pub by: Vec<Origin>,
}

enum Cell {
    Planar(PlanarCell),
    Round(RoundCell),
}

impl Cell {
    fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Planar(cell) => cell.contains_interior(point, dec),
            Self::Round(cell) => cell.contains_interior(point, dec),
        }
    }

    fn anchored(&self, definition: &DraftDefinition, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Planar(cell) => cell.anchored(&definition.neutral, dec),
            Self::Round(cell) => cell.anchored(&definition.neutral, dec),
        }
    }

    fn origins(&self) -> BTreeSet<Origin> {
        match self {
            Self::Planar(cell) => cell
                .supports
                .iter()
                .zip(&cell.cycles)
                .filter(|(_, cycle)| !cycle.is_empty())
                .flat_map(|(support, _)| support.origins.iter().cloned())
                .collect(),
            Self::Round(cell) => cell
                .origins
                .iter()
                .chain(cell.caps.iter().flatten().flat_map(|(_, origins)| origins))
                .cloned()
                .collect(),
        }
    }

    fn radius(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        let point_bound = |point: &Point3, dec: &mut Decider| -> Result<Rational, Refusal> {
            [&point.x, &point.y, &point.z]
                .into_iter()
                .try_fold(Rational::zero(), |sum, x| {
                    let interval = dec.enclose_at(x, 128)?;
                    Ok(sum + interval.lo().abs().max(interval.hi().abs()))
                })
        };
        match self {
            Self::Planar(cell) => cell
                .vertices
                .iter()
                .try_fold(Rational::zero(), |bound, point| {
                    Ok(bound.max(point_bound(point, dec)?))
                }),
            Self::Round(cell) => {
                let mut bound = Rational::zero();
                for i in 0..2 {
                    // Each coordinate's radial excursion is at most the radius.
                    let radius = dec.enclose_at(&cell.radii[i], 128)?.hi();
                    bound = bound
                        .max(point_bound(&cell.centre(i), dec)? + Rational::from_i64(3) * radius);
                }
                Ok(bound)
            }
        }
    }

    fn assemble(&self, model: &mut Model) -> BodyId {
        match self {
            Self::Planar(cell) => assemble(model, cell),
            Self::Round(cell) => assemble_round(model, cell),
        }
    }
}

fn round(
    model: &Model,
    faces: &[FaceId],
    chosen: &BTreeSet<FaceId>,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<RoundCell, Refusal> {
    let mut cylinder = None;
    let mut caps = Vec::new();
    for face in faces {
        let f = model.face(*face);
        match model.geom().surface(f.surface) {
            Surface::Cylinder(c) if cylinder.is_none() && f.sense == Sense::Same => {
                if !chosen.contains(face) {
                    return Err(Refusal::new(
                        "draft.unsupported-arrangement",
                        "unchanged curved neighbours require mixed support cells",
                    ));
                }
                cylinder = Some((c, model.face_origins(*face).to_vec()));
            }
            Surface::Plane(p) => {
                if chosen.contains(face) {
                    definition.support(&Surface::Plane(p.clone()), f.sense == Sense::Same, dec)?;
                    return Err(Refusal::new(
                        "draft.unsupported-round-cap",
                        "simultaneously drafted round caps require conic trimmed-cell construction",
                    ));
                }
                let mut p = p.clone();
                if f.sense == Sense::Reversed {
                    p.normal = p.normal.scaled(&Real::from_i64(-1));
                }
                caps.push((p, model.face_origins(*face).to_vec()));
            }
            _ => {
                return Err(Refusal::new(
                    "draft.unsupported-arrangement",
                    "round-cell construction needs one cylindrical side and two planar caps",
                ));
            }
        }
    }
    let (cylinder, origins) = cylinder.ok_or_else(|| {
        Refusal::new(
            "draft.unsupported-arrangement",
            "no axial cylindrical support",
        )
    })?;
    let caps = caps.try_into().map_err(|_| {
        Refusal::new(
            "draft.unsupported-arrangement",
            "round-cell construction needs two complete planar caps",
        )
    })?;
    RoundCell::new(cylinder, origins, caps, definition, dec)
}

fn assemble_round(model: &mut Model, cell: &RoundCell) -> BodyId {
    let mut paths = Vec::new();
    let mut apex = None;
    let mut faces = Vec::new();
    for i in 0..2 {
        let centre = cell.centre(i);
        if let Some((plane, origins)) = &cell.caps[i] {
            let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                center: centre,
                normal: cell.axis.clone(),
                radius_sq: cell.radii[i].square(),
                ref_dir: cell.reference.clone(),
            }));
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            });
            faces.push(crate::extrude::face(
                model,
                Surface::Plane(plane.clone()),
                Sense::Same,
                vec![vec![(
                    edge,
                    if i == 0 { Sense::Reversed } else { Sense::Same },
                )]],
                None,
                origins.clone(),
            ));
            paths.push(vec![(
                edge,
                if i == 0 { Sense::Same } else { Sense::Reversed },
            )]);
        } else {
            let point = model.geom_mut().add_point(PointDef::Explicit(centre));
            apex = Some(model.add_vertex(Vertex { point }));
        }
    }
    faces.push(crate::extrude::face(
        model,
        cell.surface.clone(),
        Sense::Same,
        paths,
        apex,
        cell.origins.clone(),
    ));
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

fn selected(model: &Model, body: BodyId, refs: &[FaceRef]) -> Result<BTreeSet<FaceId>, Refusal> {
    let faces = model.faces_of_body(body);
    let mut edges = BTreeMap::<_, Vec<FaceId>>::new();
    for face in &faces {
        for lp in &model.face(*face).loops {
            if let Loop::Coedges(cs) = model.loop_(*lp) {
                for c in cs {
                    edges.entry(model.coedge(*c).edge).or_default().push(*face);
                }
            }
        }
    }
    let mut result = BTreeSet::new();
    for reference in refs {
        let mut matches = Vec::new();
        for face in &faces {
            if !model
                .face_origins(*face)
                .iter()
                .any(|o| reference.face.matches(o))
            {
                continue;
            }
            if reference.touching.iter().all(|pattern| {
                edges.values().any(|uses| {
                    uses.contains(face)
                        && uses.iter().any(|f| {
                            *f != *face && model.face_origins(*f).iter().any(|o| pattern.matches(o))
                        })
                })
            }) {
                matches.push(*face);
            }
        }
        if matches.is_empty() {
            return Err(Refusal::new(
                "draft.lost-face",
                format!("no input face resolves {reference:?}"),
            ));
        }
        if reference.near.is_some() && matches.len() > 1 {
            return Err(Refusal::new(
                "draft.unsupported-near",
                "nearest-face draft selection needs exact trimmed-face distances",
            ));
        }
        result.extend(matches);
    }
    Ok(result)
}

fn planar(
    model: &Model,
    faces: &[FaceId],
    chosen: &BTreeSet<FaceId>,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<PlanarCell, Refusal> {
    let mut original = Vec::with_capacity(faces.len());
    let mut supports = Vec::with_capacity(faces.len());
    for face in faces {
        let f = model.face(*face);
        let surface = model.geom().surface(f.surface);
        let Surface::Plane(plane) = surface else {
            return Err(Refusal::new(
                "draft.unsupported-arrangement",
                "this support arrangement requires curved cell construction",
            ));
        };
        let mut plane = plane.clone();
        if f.sense == Sense::Reversed {
            plane.normal = plane.normal.scaled(&Real::from_i64(-1));
        }
        original.push(plane.clone());
        if chosen.contains(face) {
            let replacement = definition.support(surface, f.sense == Sense::Same, dec)?;
            let Surface::Plane(replacement) = replacement.surface else {
                unreachable!("plane draft remains planar")
            };
            plane = replacement;
        }
        supports.push(PlanePatch {
            plane,
            origins: model.face_origins(*face).to_vec(),
        });
    }
    // Convexity is proved against the immutable input skin, not guessed from
    // face normals or the new arrangement's successful closure.
    let mut vertices = BTreeSet::new();
    for face in faces {
        for lp in &model.face(*face).loops {
            if let Loop::Coedges(cs) = model.loop_(*lp) {
                for c in cs {
                    let edge = model.edge(model.coedge(*c).edge);
                    vertices.extend([edge.start, edge.end].into_iter().flatten());
                }
            }
        }
    }
    for vertex in vertices {
        let point = model.geom().position(model.vertex(vertex).point);
        for plane in &original {
            if sign(&point.minus(&plane.origin).dot(&plane.normal), dec)? == Sign::Positive {
                return Err(Refusal::new(
                    "draft.unsupported-arrangement",
                    "nonconvex planar input requires a multi-cell support arrangement",
                ));
            }
        }
    }
    PlanarCell::new(supports, dec)?.ok_or_else(|| {
        Refusal::new(
            "draft.reversed-solid",
            "replaced boundary has no full-dimensional material cell",
        )
    })
}

fn assemble(model: &mut Model, cell: &PlanarCell) -> BodyId {
    let vertices: Vec<_> = cell
        .vertices
        .iter()
        .map(|point| {
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit(point.clone()));
            model.add_vertex(Vertex { point })
        })
        .collect();
    let mut edges = BTreeMap::new();
    let mut faces = Vec::new();
    for (support, cycle) in cell.supports.iter().zip(&cell.cycles) {
        if cycle.is_empty() {
            continue;
        }
        let mut path = Vec::with_capacity(cycle.len());
        for i in 0..cycle.len() {
            let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
            let pair = (a.min(b), a.max(b));
            let edge = *edges.entry(pair).or_insert_with(|| {
                let curve = model.geom_mut().add_curve(Curve::Line(Line {
                    origin: cell.vertices[pair.0].clone(),
                    dir: cell.vertices[pair.1].minus(&cell.vertices[pair.0]),
                }));
                model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: Some(vertices[pair.0]),
                    end: Some(vertices[pair.1]),
                })
            });
            path.push((edge, if a < b { Sense::Same } else { Sense::Reversed }));
        }
        faces.push(crate::extrude::face(
            model,
            Surface::Plane(support.plane.clone()),
            Sense::Same,
            vec![path],
            None,
            support.origins.clone(),
        ));
    }
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

pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    refs: &[FaceRef],
    definition: &DraftDefinition,
    keep_regions: Option<&[catallax_recipe::Point3]>,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<ConsumedFace>, Rational), Refusal> {
    let chosen = selected(model, body, refs)?;
    let input = model.body(body).clone();
    let mut components = Vec::new();
    for lump in &input.lumps {
        let shells = &model.lump(*lump).shells;
        if shells.len() != 1 {
            return Err(Refusal::new(
                "draft.unsupported-arrangement",
                "multiple-shell draft needs explicit oriented void cells",
            ));
        }
        let faces = model.shell(shells[0]).faces.clone();
        let cell = if faces.iter().all(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Plane(_)
            )
        }) {
            Cell::Planar(planar(model, &faces, &chosen, definition, dec)?)
        } else {
            Cell::Round(round(model, &faces, &chosen, definition, dec)?)
        };
        components.push((faces, cell));
    }
    let mut keep = vec![false; components.len()];
    if let Some(points) = keep_regions {
        for point in points {
            let [x, y, z] = point.0.each_ref().map(|n| n.value().clone().into());
            let point = Point3::new(x, y, z);
            let mut found = false;
            for (i, (_, cell)) in components.iter().enumerate() {
                if cell.contains_interior(&point, dec)? {
                    keep[i] = true;
                    found = true;
                }
            }
            if !found {
                return Err(Refusal::new(
                    "draft.invalid-keep-region",
                    "saved region point is on a boundary or outside all candidate material",
                ));
            }
        }
    } else {
        for (i, (_, cell)) in components.iter().enumerate() {
            keep[i] = cell.anchored(definition, dec)?;
        }
    }
    let mut result = input;
    result.lumps.clear();
    let mut consumed = Vec::new();
    let by: Vec<_> = chosen
        .iter()
        .flat_map(|f| model.face_origins(*f).iter().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut radius = Rational::zero();
    for ((faces, cell), keep) in components.into_iter().zip(keep) {
        let surviving = if keep {
            cell.origins()
        } else {
            BTreeSet::new()
        };
        for face in faces {
            let source = model.face_origins(face);
            if source.iter().all(|o| !surviving.contains(o)) {
                consumed.push(ConsumedFace {
                    step: step.clone(),
                    source: source.to_vec(),
                    by: by.clone(),
                });
            }
        }
        if !keep {
            continue;
        }
        radius = radius.max(cell.radius(dec)?);
        let component = cell.assemble(model);
        result.lumps.extend_from_slice(&model.body(component).lumps);
    }
    result.kind = if result.lumps.is_empty() {
        BodyKind::Empty
    } else {
        BodyKind::Solid
    };
    let result = model.add_body(result);
    let faults = model.check(result);
    if !faults.is_empty() {
        return Err(Refusal::new(
            "draft.invalid-topology",
            format!("{faults:?}"),
        ));
    }
    consumed.sort_by(|a, b| a.source.cmp(&b.source));
    consumed.dedup_by(|a, b| a.source == b.source);
    Ok((result, consumed, radius))
}
