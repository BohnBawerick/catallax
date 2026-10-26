//! Immutable topology construction from drafted support arrangements.

use crate::draft::{
    DraftDefinition,
    annular::{AnnularCell, AnnularWall},
    arrangement::{Patch, Polyhedron, Skin},
    capped::CappedCylinder,
    cell::{PlanarCell, PlanePatch},
    near,
    round::RoundCell,
};
use catallax_numbers::{Decider, Rational, Real, Refusal};
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
    Annular(Box<AnnularCell>),
    Planar(PlanarCell),
    Round(Box<RoundCell>),
    Polyhedral(Box<Polyhedron>),
    Capped(Box<CappedCylinder>),
}

impl Cell {
    fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Annular(cell) => cell.contains_interior(point, dec),
            Self::Planar(cell) => cell.contains_interior(point, dec),
            Self::Round(cell) => cell.contains_interior(point, dec),
            Self::Polyhedral(cell) => cell.contains_interior(point, dec),
            Self::Capped(cell) => cell.contains_interior(point, dec),
        }
    }

    fn anchored(&self, definition: &DraftDefinition, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Annular(cell) => cell.anchored(&definition.neutral, dec),
            Self::Planar(cell) => cell.anchored(&definition.neutral, dec),
            Self::Round(cell) => cell.anchored(&definition.neutral, dec),
            Self::Polyhedral(cell) => cell.anchored(&definition.neutral, dec),
            Self::Capped(cell) => cell.anchored(&definition.neutral, dec),
        }
    }

    fn origins(&self) -> BTreeSet<Origin> {
        match self {
            Self::Annular(cell) => cell
                .outer
                .origins
                .iter()
                .chain(&cell.inner.origins)
                .chain(
                    cell.outer
                        .caps
                        .iter()
                        .flatten()
                        .flat_map(|(_, origins)| origins),
                )
                .cloned()
                .collect(),
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
            Self::Polyhedral(cell) => cell
                .faces
                .iter()
                .flat_map(|f| f.origins.iter().cloned())
                .collect(),
            Self::Capped(cell) => cell
                .origins
                .iter()
                .chain(cell.caps.iter().flat_map(|(_, names)| names))
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
            Self::Annular(cell) => Self::Round(Box::new(cell.outer.clone())).radius(dec),
            Self::Planar(cell) => cell
                .vertices
                .iter()
                .try_fold(Rational::zero(), |bound, point| {
                    Ok(bound.max(point_bound(point, dec)?))
                }),
            Self::Polyhedral(cell) => cell
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
            Self::Capped(cell) => cell.bound_radius(dec),
        }
    }

    fn assemble(&self, model: &mut Model, dec: &mut Decider) -> Result<BodyId, Refusal> {
        Ok(match self {
            Self::Annular(cell) => assemble_annular(model, cell),
            Self::Planar(cell) => assemble(model, cell),
            Self::Round(cell) => assemble_round(model, cell),
            Self::Polyhedral(cell) => assemble_polyhedral(model, cell),
            Self::Capped(cell) => assemble_capped(model, cell, dec)?,
        })
    }
}

fn round(
    model: &Model,
    faces: &[FaceId],
    chosen: &BTreeSet<FaceId>,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<Cell, Refusal> {
    if faces
        .iter()
        .filter(|f| {
            matches!(
                model.geom().surface(model.face(**f).surface),
                Surface::Cylinder(_)
            )
        })
        .count()
        == 2
    {
        return annular(model, faces, chosen, definition, dec).map(|c| Cell::Annular(Box::new(c)));
    }
    let mut cylinder = None;
    let mut caps = Vec::new();
    for face in faces {
        let f = model.face(*face);
        match model.geom().surface(f.surface) {
            Surface::Cylinder(c) if cylinder.is_none() && f.sense == Sense::Same => {
                cylinder = Some((c, model.face_origins(*face).to_vec(), chosen.contains(face)));
            }
            Surface::Plane(p) => {
                let mut p = if chosen.contains(face) {
                    let support = definition.support(
                        &Surface::Plane(p.clone()),
                        f.sense == Sense::Same,
                        dec,
                    )?;
                    let Surface::Plane(mut plane) = support.surface else {
                        unreachable!("a drafted plane is planar")
                    };
                    if !support.same_sense {
                        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
                    }
                    caps.push((plane, model.face_origins(*face).to_vec()));
                    continue;
                } else {
                    p.clone()
                };
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
    let (cylinder, origins, move_side) = cylinder.ok_or_else(|| {
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
    if move_side {
        RoundCell::new(cylinder, origins, caps, definition, dec)
            .map(|cell| Cell::Round(Box::new(cell)))
    } else {
        CappedCylinder::new(cylinder.clone(), origins, caps, dec)
            .map(|cell| Cell::Capped(Box::new(cell)))
    }
}

fn annular(
    model: &Model,
    faces: &[FaceId],
    chosen: &BTreeSet<FaceId>,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<AnnularCell, Refusal> {
    let unsupported = || {
        Refusal::new(
            "draft.unsupported-arrangement",
            "annular construction needs two complete cylindrical walls and two complete annular caps",
        )
    };
    let mut walls = Vec::new();
    let mut caps = Vec::new();
    if faces.len() != 4 {
        return Err(unsupported());
    }
    for id in faces {
        let face = model.face(*id);
        if face.loops.len() != 2 {
            return Err(unsupported());
        }
        for lp in &face.loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                return Err(unsupported());
            };
            let [c] = cs.as_slice() else {
                return Err(unsupported());
            };
            let edge = model.edge(model.coedge(*c).edge);
            if edge.start.is_some()
                || edge.end.is_some()
                || !matches!(model.geom().curve(edge.curve), Curve::Circle(_))
            {
                return Err(unsupported());
            }
        }
        let origins = model.face_origins(*id).to_vec();
        match model.geom().surface(face.surface) {
            Surface::Cylinder(cylinder) => walls.push(AnnularWall {
                cylinder: cylinder.clone(),
                origins,
                same_sense: face.sense == Sense::Same,
                selected: chosen.contains(id),
            }),
            Surface::Plane(plane) => {
                let mut plane = plane.clone();
                if chosen.contains(id) {
                    definition.support(
                        &Surface::Plane(plane.clone()),
                        face.sense == Sense::Same,
                        dec,
                    )?;
                    return Err(unsupported());
                }
                if face.sense == Sense::Reversed {
                    plane.normal = plane.normal.scaled(&Real::from_i64(-1));
                }
                caps.push((plane, origins));
            }
            _ => return Err(unsupported()),
        }
    }
    AnnularCell::new(
        walls.try_into().map_err(|_| unsupported())?,
        caps.try_into().map_err(|_| unsupported())?,
        definition,
        dec,
    )
}

fn assemble_annular(model: &mut Model, cell: &AnnularCell) -> BodyId {
    let mut paths = [Vec::new(), Vec::new()];
    let mut faces = Vec::new();
    for i in 0..2 {
        let mut cap_paths = Vec::new();
        for (j, wall) in [&cell.outer, &cell.inner].into_iter().enumerate() {
            let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                center: wall.centre(i),
                normal: wall.axis.clone(),
                radius_sq: wall.radii[i].square(),
                ref_dir: wall.reference.clone(),
            }));
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            });
            let cap_sense = if i == j { Sense::Reversed } else { Sense::Same };
            cap_paths.push(vec![(edge, cap_sense)]);
            paths[j].push(vec![(edge, crate::extrude::reverse(cap_sense))]);
        }
        let (plane, origins) = cell.outer.caps[i]
            .as_ref()
            .expect("certified nonvanishing annular caps");
        faces.push(crate::extrude::face(
            model,
            Surface::Plane(plane.clone()),
            Sense::Same,
            cap_paths,
            None,
            origins.clone(),
        ));
    }
    for (j, (wall, paths)) in [&cell.outer, &cell.inner]
        .into_iter()
        .zip(paths)
        .enumerate()
    {
        faces.push(crate::extrude::face(
            model,
            wall.surface.clone(),
            if j == 0 { Sense::Same } else { Sense::Reversed },
            paths,
            None,
            wall.origins.clone(),
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

fn assemble_capped(
    model: &mut Model,
    cell: &CappedCylinder,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    use crate::draft::capped::FootprintBoundary;
    let (vertices, chord) = match cell.footprint_boundary(dec)? {
        FootprintBoundary::Closed => (None, None),
        FootprintBoundary::Chord(q) => {
            let points = q.each_ref().map(|q| cell.point(q, &cell.heights_at(q)[0]));
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                origin: points[0].clone(),
                dir: points[1].minus(&points[0]),
            }));
            let vertices = points.map(|point| {
                let point = model.geom_mut().add_point(PointDef::Explicit(point));
                model.add_vertex(Vertex { point })
            });
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(vertices[0]),
                end: Some(vertices[1]),
            });
            (Some(vertices), Some(edge))
        }
        FootprintBoundary::Tangent(q) => {
            let point = cell.point(&q, &cell.heights_at(&q)[0]);
            let point = model.geom_mut().add_point(PointDef::Explicit(point));
            let vertex = model.add_vertex(Vertex { point });
            (Some([vertex, vertex]), None)
        }
    };
    let right_handed = crate::draft::sign(&cell.u.cross(&cell.v).dot(&cell.axis), dec)?
        == catallax_numbers::Sign::Positive;
    let mut faces = Vec::with_capacity(3);
    let mut side_arcs = Vec::with_capacity(2);
    for i in 0..2 {
        let curve = model.geom_mut().add_curve(cell.cap_curve(i, dec)?);
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: vertices.map(|v| v[0]),
            end: vertices.map(|v| v[1]),
        });
        let cap_sense = if (i == 1) == right_handed {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let side_sense = if cap_sense == Sense::Same {
            Sense::Reversed
        } else {
            Sense::Same
        };
        let mut boundary = vec![(edge, cap_sense)];
        if let Some(chord) = chord {
            boundary.push((chord, side_sense));
        }
        faces.push(crate::extrude::face(
            model,
            Surface::Plane(cell.caps[i].0.clone()),
            Sense::Same,
            vec![boundary],
            None,
            cell.caps[i].1.clone(),
        ));
        side_arcs.push((edge, side_sense));
    }
    // A cap contact identifies the annulus ends: both arcs then form one wire.
    let side_loops = if vertices.is_some() {
        vec![side_arcs]
    } else {
        side_arcs.into_iter().map(|arc| vec![arc]).collect()
    };
    faces.push(crate::extrude::face(
        model,
        Surface::Cylinder(cell.cylinder.clone()),
        Sense::Same,
        side_loops,
        None,
        cell.origins.clone(),
    ));
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

fn selected(
    model: &Model,
    body: BodyId,
    refs: &[FaceRef],
    dec: &mut Decider,
) -> Result<BTreeSet<FaceId>, Refusal> {
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
                "recipe.missing-reference",
                format!("no input face resolves {reference:?}"),
            ));
        }
        if let Some(near) = &reference.near
            && matches.len() > 1
        {
            matches = vec![nearest_axial(model, &faces, &matches, near, dec)?];
        }
        result.extend(matches);
    }
    Ok(result)
}

/// Picks the uniquely nearest closed trimmed face of a plane-capped axial
/// cylinder, or of a coaxial annulus with two walls and two annular caps.
fn nearest_axial(
    model: &Model,
    faces: &[FaceId],
    matches: &[FaceId],
    near: &catallax_recipe::Point3,
    dec: &mut Decider,
) -> Result<FaceId, Refusal> {
    let unsupported = || {
        Refusal::new(
            "draft.unsupported-near",
            "nearest-face draft selection needs one or two coaxial cylindrical sides and two planar caps",
        )
    };
    let mut cylinders = Vec::new();
    let mut caps = Vec::new();
    for face in faces {
        match model.geom().surface(model.face(*face).surface) {
            Surface::Cylinder(c) if cylinders.len() < 2 => cylinders.push(c),
            Surface::Plane(p) => caps.push(p),
            _ => return Err(unsupported()),
        }
    }
    let [first, second] = caps.as_slice() else {
        return Err(unsupported());
    };
    let walls = match cylinders.as_slice() {
        [_] => None,
        [a, b] => Some(near::annulus_walls(a, b, dec)?),
        _ => return Err(unsupported()),
    };
    let [x, y, z] = near.0.each_ref().map(|n| n.value().clone().into());
    let point = Point3::new(x, y, z);
    let mut distances = Vec::with_capacity(matches.len());
    for face in matches {
        distances.push(match model.geom().surface(model.face(*face).surface) {
            Surface::Plane(cap) => match walls {
                Some(walls) => near::annular_cap_distance_squared(walls, cap, &point, dec)?,
                None => near::cap_distance_squared(cylinders[0], cap, &point, dec)?,
            },
            Surface::Cylinder(wall) => {
                near::wall_distance_squared(wall, [first, second], &point, dec)?
            }
            _ => return Err(unsupported()),
        });
    }
    Ok(matches[near::unique_nearest(&distances, dec)?])
}

fn planar_skin(model: &Model, faces: &[FaceId]) -> Result<Skin, Refusal> {
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
    let indices: BTreeMap<_, _> = vertices.iter().enumerate().map(|(i, v)| (*v, i)).collect();
    let points = vertices
        .iter()
        .map(|v| model.geom().position(model.vertex(*v).point))
        .collect();
    let mut patches = Vec::with_capacity(faces.len());
    for face in faces {
        let f = model.face(*face);
        let Surface::Plane(plane) = model.geom().surface(f.surface) else {
            return Err(Refusal::new(
                "draft.unsupported-arrangement",
                "planar construction needs planar carriers",
            ));
        };
        let mut plane = plane.clone();
        if f.sense == Sense::Reversed {
            plane.normal = plane.normal.scaled(&Real::from_i64(-1));
        }
        let mut loops = Vec::new();
        for lp in &f.loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                return Err(Refusal::new(
                    "draft.unsupported-arrangement",
                    "planar skin requires finite straight wires",
                ));
            };
            let mut cycle = Vec::with_capacity(cs.len());
            for c in cs {
                let coedge = model.coedge(*c);
                let edge = model.edge(coedge.edge);
                if !matches!(model.geom().curve(edge.curve), Curve::Line(_)) {
                    return Err(Refusal::new(
                        "draft.unsupported-arrangement",
                        "planar skin requires line intersections",
                    ));
                }
                let vertex = if coedge.sense == Sense::Same {
                    edge.start
                } else {
                    edge.end
                };
                let Some(vertex) = vertex else {
                    return Err(Refusal::new(
                        "draft.nonclosing-boundary",
                        "planar skin wire has no endpoint",
                    ));
                };
                cycle.push(indices[&vertex]);
            }
            loops.push(cycle);
        }
        patches.push(Patch {
            plane,
            origins: model.face_origins(*face).to_vec(),
            loops,
        });
    }
    Ok(Skin {
        vertices: points,
        patches,
    })
}

fn planar(
    skin: &Skin,
    chosen: &BTreeSet<usize>,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<Vec<Cell>, Refusal> {
    if skin.is_convex(dec)? {
        let mut original = Vec::with_capacity(skin.patches.len());
        for (index, patch) in skin.patches.iter().enumerate() {
            let plane = if chosen.contains(&index) {
                let Surface::Plane(plane) = definition
                    .support(&Surface::Plane(patch.plane.clone()), true, dec)?
                    .surface
                else {
                    unreachable!("plane draft remains planar")
                };
                plane
            } else {
                patch.plane.clone()
            };
            original.push(PlanePatch {
                plane,
                origins: patch.origins.clone(),
            });
        }
        let cell = PlanarCell::new(original, dec)?.ok_or_else(|| {
            Refusal::new(
                "draft.reversed-solid",
                "replaced boundary has no full-dimensional material cell",
            )
        })?;
        return Ok(vec![Cell::Planar(cell)]);
    }
    skin.draft(chosen, definition, dec).map(|cells| {
        cells
            .into_iter()
            .map(|c| Cell::Polyhedral(Box::new(c)))
            .collect()
    })
}

fn assemble_polyhedral(model: &mut Model, cell: &Polyhedron) -> BodyId {
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
    for face in &cell.faces {
        let mut paths = Vec::new();
        for cycle in &face.loops {
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
            paths.push(path);
        }
        faces.push(crate::extrude::face(
            model,
            Surface::Plane(cell.plane(face)),
            Sense::Same,
            paths,
            None,
            face.origins.clone(),
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
    let source_faces = model.faces_of_body(body);
    let (chosen, mut components) = if !source_faces.is_empty()
        && source_faces.iter().all(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Plane(_)
            )
        }) {
        let skin = planar_skin(model, &source_faces)?;
        let selected = skin.resolve_faces(refs, dec)?;
        let chosen = selected.iter().map(|i| source_faces[*i]).collect();
        (chosen, planar(&skin, &selected, definition, dec)?)
    } else {
        (selected(model, body, refs, dec)?, Vec::new())
    };
    let input = model.body(body).clone();
    if components.is_empty() && !source_faces.is_empty() {
        for lump in &input.lumps {
            let shells = &model.lump(*lump).shells;
            if shells.len() != 1 {
                return Err(Refusal::new(
                    "draft.unsupported-arrangement",
                    "mixed curved shells require oriented curved cells",
                ));
            }
            let faces = &model.shell(shells[0]).faces;
            components.push(round(model, faces, &chosen, definition, dec)?);
        }
    }
    if components.len() > 1
        && components
            .iter()
            .any(|c| matches!(c, Cell::Capped(_) | Cell::Annular(_)))
    {
        return Err(Refusal::new(
            "draft.unsupported-arrangement",
            "several plane-capped cylinders need a certified curved component arrangement",
        ));
    }
    for (i, component) in components.iter().enumerate() {
        let Cell::Round(a) = component else {
            continue;
        };
        for other in &components[..i] {
            if let Cell::Round(b) = other
                && !a.certified_disjoint(b, dec)?
            {
                return Err(Refusal::new(
                    "draft.unsupported-arrangement",
                    "intersecting round cells need reconstructed curved boundary patches",
                ));
            }
        }
    }
    let needs_maximal = components.iter().any(|c| matches!(c, Cell::Polyhedral(_)));
    let mut keep = vec![false; components.len()];
    if let Some(points) = keep_regions {
        for point in points {
            let [x, y, z] = point.0.each_ref().map(|n| n.value().clone().into());
            let point = Point3::new(x, y, z);
            let mut found = false;
            for (i, cell) in components.iter().enumerate() {
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
        for (i, cell) in components.iter().enumerate() {
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
    let mut surviving = BTreeSet::new();
    for (cell, keep) in components.into_iter().zip(keep) {
        if keep {
            surviving.extend(cell.origins());
        }
        if !keep {
            continue;
        }
        radius = radius.max(cell.radius(dec)?);
        let component = cell.assemble(model, dec)?;
        result.lumps.extend_from_slice(&model.body(component).lumps);
    }
    for face in source_faces {
        let source = model.face_origins(face);
        if source.iter().all(|o| !surviving.contains(o)) {
            consumed.push(ConsumedFace {
                step: step.clone(),
                source: source.to_vec(),
                by: by.clone(),
            });
        }
    }
    result.kind = if result.lumps.is_empty() {
        BodyKind::Empty
    } else {
        BodyKind::Solid
    };
    let result = model.add_body(result);
    let result = if needs_maximal {
        catallax_boolean::make_maximal(model, result, dec)?.value
    } else {
        result
    };
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
