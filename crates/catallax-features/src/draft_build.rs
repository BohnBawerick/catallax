//! Immutable topology construction from drafted support arrangements.

use crate::draft::{
    DraftDefinition,
    annular::{AnnularCell, AnnularWall},
    arrangement::{Patch, Polyhedron, Skin},
    bored::{self, Bore, BoredCell},
    bossed::{self, Boss, BossedCell},
    capped::CappedCylinder,
    cell::{PlanarCell, PlanePatch},
    near,
    round::RoundCell,
};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_recipe::{FaceRef, Origin, StepId};
use catallax_shapes::{Circle, Curve, Line, Point3, PointDef, Surface};
use catallax_topology::{
    Body, BodyId, BodyKind, Edge, EdgeId, FaceId, Loop, LoopId, Lump, Model, Sense, Shell, Vertex,
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

/// One vertex-free circular wire shared by a cap hole and the bore wall.
#[derive(Clone)]
struct BoreRing {
    curve: Curve,
    edge_sense: Sense,
    cap_sense: Sense,
    wall_sense: Sense,
    cap: FaceId,
}

/// One inward bore wall and the two rings that trim it.
struct SourceBore {
    wall: FaceId,
    rings: [BoreRing; 2],
}

/// A planar envelope's source wires and its unselected bore walls.
struct BoredSource {
    skin: Skin,
    planes: Vec<FaceId>,
    bores: Vec<SourceBore>,
}

/// One preserved wall's sense, and its rings bound to their cap supports.
struct BoundBore {
    rings: [(BoreRing, usize); 2],
    sense: Sense,
}

/// A planar base's source wires around one outward wall. The wall and top are
/// unchanged, so only the base's planes enter the envelope.
struct BossedSource {
    skin: Skin,
    planes: Vec<FaceId>,
    wall: FaceId,
    /// The wall's rings, shoulder first, then top.
    rings: [BoreRing; 2],
    /// A separate bore through the base, if any.
    bores: Vec<SourceBore>,
}

enum Cell {
    Annular(Box<AnnularCell>),
    Bored(Box<BoredCell>, Vec<BoundBore>),
    Bossed(Box<BossedCell>, Box<[BoreRing; 2]>, Vec<BoundBore>),
    Planar(PlanarCell),
    Round(Box<RoundCell>),
    Polyhedral(Box<Polyhedron>),
    Capped(Box<CappedCylinder>),
}

impl Cell {
    fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Annular(cell) => cell.contains_interior(point, dec),
            Self::Bored(cell, ..) => cell.contains_interior(point, dec),
            Self::Bossed(cell, ..) => cell.contains_interior(point, dec),
            Self::Planar(cell) => cell.contains_interior(point, dec),
            Self::Round(cell) => cell.contains_interior(point, dec),
            Self::Polyhedral(cell) => cell.contains_interior(point, dec),
            Self::Capped(cell) => cell.contains_interior(point, dec),
        }
    }

    fn anchored(&self, definition: &DraftDefinition, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Annular(cell) => cell.anchored(&definition.neutral, dec),
            // Strict clearance leaves positive area around the removed disk.
            Self::Bored(cell, ..) => cell.outer.anchored(&definition.neutral, dec),
            Self::Bossed(cell, ..) => cell.anchored(&definition.neutral, dec),
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
            Self::Bored(cell, ..) => Self::Planar(cell.outer.clone())
                .origins()
                .into_iter()
                .chain(
                    cell.bores
                        .iter()
                        .flat_map(|bore| bore.origins.iter().cloned()),
                )
                .collect(),
            Self::Bossed(cell, ..) => Self::Planar(cell.outer.clone())
                .origins()
                .into_iter()
                .chain(cell.boss.origins.iter().cloned())
                .chain(cell.boss.top_origins.iter().cloned())
                .chain(
                    cell.bored
                        .iter()
                        .flat_map(|bored| &bored.bores)
                        .flat_map(|bore| bore.origins.iter().cloned()),
                )
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
            Self::Bored(cell, ..) => Self::Planar(cell.outer.clone()).radius(dec),
            Self::Bossed(cell, ..) => {
                // The shoulder ring lies in the base's hull. Only the top reaches out.
                let boss = &cell.boss;
                let centre = bored::axis_point(&boss.cylinder, &boss.top);
                let radius = dec.enclose_at(&boss.cylinder.radius_sq.sqrt(), 128)?.hi();
                Ok(Self::Planar(cell.outer.clone())
                    .radius(dec)?
                    .max(point_bound(&centre, dec)? + Rational::from_i64(3) * radius))
            }
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
            Self::Bored(cell, bores) => assemble_bored(model, cell, bores),
            Self::Bossed(cell, rings, bores) => assemble_bossed(model, cell, rings, bores),
            Self::Planar(cell) => assemble(model, cell, &[], Vec::new()),
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

fn planar_skin(model: &Model, faces: &[FaceId], holes: &BTreeSet<LoopId>) -> Result<Skin, Refusal> {
    let mut vertices = BTreeSet::new();
    for face in faces {
        for lp in &model.face(*face).loops {
            if holes.contains(lp) {
                continue;
            }
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
            if holes.contains(lp) {
                continue;
            }
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

fn convex_cell(
    skin: &Skin,
    chosen: &BTreeSet<usize>,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<PlanarCell, Refusal> {
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
    PlanarCell::new(original, dec)?.ok_or_else(|| {
        Refusal::new(
            "draft.reversed-solid",
            "replaced boundary has no full-dimensional material cell",
        )
    })
}

fn planar(
    skin: &Skin,
    chosen: &BTreeSet<usize>,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<Vec<Cell>, Refusal> {
    if skin.is_convex(dec)? {
        return Ok(vec![Cell::Planar(convex_cell(
            skin, chosen, definition, dec,
        )?)]);
    }
    skin.draft(chosen, definition, dec).map(|cells| {
        cells
            .into_iter()
            .map(|c| Cell::Polyhedral(Box::new(c)))
            .collect()
    })
}

/// A wall's two complete rings, each with the one planar loop that shares it.
/// `None` when the wall has any other boundary.
fn wall_rings(model: &Model, wall: FaceId, planes: &[FaceId]) -> Option<[(BoreRing, LoopId); 2]> {
    let mut rings = Vec::with_capacity(2);
    for lp in &model.face(wall).loops {
        let Loop::Coedges(cs) = model.loop_(*lp) else {
            return None;
        };
        let [coedge] = cs.as_slice() else {
            return None;
        };
        let coedge = model.coedge(*coedge);
        let edge = model.edge(coedge.edge);
        let curve = model.geom().curve(edge.curve);
        if edge.start.is_some() || edge.end.is_some() || !matches!(curve, Curve::Circle(_)) {
            return None;
        }
        let mut uses = planes.iter().flat_map(|face| {
            model
                .face(*face)
                .loops
                .iter()
                .filter_map(|lp| match model.loop_(*lp) {
                    Loop::Coedges(cs) => match cs.as_slice() {
                        [c] if model.coedge(*c).edge == coedge.edge => {
                            Some((*face, *lp, model.coedge(*c).sense))
                        }
                        _ => None,
                    },
                    Loop::Vertex(_) => None,
                })
        });
        let (Some((cap, hole, cap_sense)), None) = (uses.next(), uses.next()) else {
            return None;
        };
        rings.push((
            BoreRing {
                curve: curve.clone(),
                edge_sense: edge.sense,
                cap_sense,
                wall_sense: coedge.sense,
                cap,
            },
            hole,
        ));
    }
    rings.try_into().ok()
}

/// Recognises planar faces around inward walls, each trimmed by two complete
/// rings. Any other shape returns `None`, keeping its established route.
fn bored_source(model: &Model, faces: &[FaceId]) -> Result<Option<BoredSource>, Refusal> {
    let unsupported = || {
        Refusal::new(
            "draft.unsupported-arrangement",
            "a preserved bore needs one wall between two complete cap rings",
        )
    };
    let mut walls = Vec::new();
    let mut planes = Vec::new();
    for face in faces {
        match model.geom().surface(model.face(*face).surface) {
            Surface::Cylinder(_) if model.face(*face).sense == Sense::Reversed => walls.push(*face),
            Surface::Plane(_) => planes.push(*face),
            _ => return Ok(None),
        }
    }
    if walls.is_empty() {
        return Ok(None);
    }
    let mut holes = BTreeSet::new();
    let mut caps = BTreeSet::new();
    let mut bores = Vec::with_capacity(walls.len());
    for wall in walls {
        let [(a, a_hole), (b, b_hole)] =
            wall_rings(model, wall, &planes).ok_or_else(unsupported)?;
        if a.cap == b.cap {
            return Err(unsupported());
        }
        holes.extend([a_hole, b_hole]);
        caps.extend([a.cap, b.cap]);
        bores.push(SourceBore {
            wall,
            rings: [a, b],
        });
    }
    // A blind floor or a shoulder has only rings. Every bore needs the same two
    // separate caps, each holding one outer wire and exactly these rings.
    if caps.len() != 2
        || caps
            .iter()
            .any(|cap| model.face(*cap).loops.len() != bores.len() + 1)
    {
        return Err(unsupported());
    }
    let skin = planar_skin(model, &planes, &holes)?;
    Ok(Some(BoredSource {
        skin,
        planes,
        bores,
    }))
}

/// Recognises planar faces around one outward wall that rises from a hole in one
/// face to the sole wire of another, and at most one separate inward wall between
/// that face and a third. Any other shape returns `None`, keeping its established
/// route and refusal.
fn bossed_source(model: &Model, faces: &[FaceId]) -> Result<Option<BossedSource>, Refusal> {
    let mut wall = None;
    let mut bore = None;
    let mut planes = Vec::new();
    for face in faces {
        let sense = model.face(*face).sense;
        match model.geom().surface(model.face(*face).surface) {
            Surface::Cylinder(_) if sense == Sense::Same && wall.is_none() => wall = Some(*face),
            Surface::Cylinder(_) if sense == Sense::Reversed && bore.is_none() => {
                bore = Some(*face);
            }
            Surface::Plane(_) => planes.push(*face),
            _ => return Ok(None),
        }
    }
    let Some(wall) = wall else {
        return Ok(None);
    };
    let Some([a, b]) = wall_rings(model, wall, &planes) else {
        return Ok(None);
    };
    let wires = 2 + usize::from(bore.is_some());
    let loops = |ring: &(BoreRing, LoopId)| model.face(ring.0.cap).loops.len();
    let ((shoulder, hole), (top, _)) = match (loops(&a), loops(&b)) {
        (n, 1) if n == wires => (a, b),
        (1, n) if n == wires => (b, a),
        _ => return Ok(None),
    };
    let mut holes = BTreeSet::from([hole]);
    let mut bores = Vec::new();
    if let Some(bore) = bore {
        let Some([(c, c_hole), (d, d_hole)]) = wall_rings(model, bore, &planes) else {
            return Ok(None);
        };
        // One mouth opens on the shoulder and the other on a separate cap.
        let floor = match (c.cap == shoulder.cap, d.cap == shoulder.cap) {
            (true, false) => d.cap,
            (false, true) => c.cap,
            _ => return Ok(None),
        };
        if model.face(floor).loops.len() != 2 {
            return Ok(None);
        }
        holes.extend([c_hole, d_hole]);
        bores.push(SourceBore {
            wall: bore,
            rings: [c, d],
        });
    }
    planes.retain(|face| *face != top.cap);
    let skin = planar_skin(model, &planes, &holes)?;
    Ok(Some(BossedSource {
        skin,
        planes,
        wall,
        rings: [shoulder, top],
        bores,
    }))
}

/// The unchanged definitions of unselected source bores.
fn preserved_bores(
    model: &Model,
    source: &[SourceBore],
    refs: &[FaceRef],
) -> Result<Vec<Bore>, Refusal> {
    let mut bores = Vec::with_capacity(source.len());
    for bore in source {
        let origins = model.face_origins(bore.wall).to_vec();
        if refs
            .iter()
            .any(|r| origins.iter().any(|origin| r.face.matches(origin)))
        {
            return Err(Refusal::new(
                "draft.unsupported-arrangement",
                "a selected bore needs a curved arrangement",
            ));
        }
        let Surface::Cylinder(cylinder) = model.geom().surface(model.face(bore.wall).surface)
        else {
            unreachable!("bored sources keep a cylindrical wall")
        };
        bores.push(Bore {
            cylinder: cylinder.clone(),
            origins,
        });
    }
    Ok(bores)
}

/// Binds each source bore's rings to the certified cell's cap supports.
fn bind_bores(
    model: &Model,
    source: &[SourceBore],
    cell: &BoredCell,
) -> Result<Vec<BoundBore>, Refusal> {
    let lost = || {
        Refusal::new(
            "draft.unsupported-arrangement",
            "a bore ring lost its cap support",
        )
    };
    let mut bound = Vec::with_capacity(source.len());
    for bore in source {
        let mut rings = bore.rings.clone().map(|ring| {
            let support = cell.outer.supports.iter().position(|support| {
                model
                    .face_origins(ring.cap)
                    .iter()
                    .all(|origin| support.origins.contains(origin))
            });
            (ring, support)
        });
        rings.sort_by_key(|(_, support)| *support);
        let [(a, Some(lower)), (b, Some(upper))] = rings else {
            return Err(lost());
        };
        if lower == upper || !cell.caps.contains(&lower) || !cell.caps.contains(&upper) {
            return Err(lost());
        }
        bound.push(BoundBore {
            rings: [(a, lower), (b, upper)],
            sense: model.face(bore.wall).sense,
        });
    }
    Ok(bound)
}

/// The drafted envelope minus its unchanged bores, proved before any topology.
fn bored_cell(
    model: &Model,
    source: &BoredSource,
    refs: &[FaceRef],
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<(BTreeSet<FaceId>, Cell), Refusal> {
    let unsupported = |detail| Refusal::new("draft.unsupported-arrangement", detail);
    let bores = preserved_bores(model, &source.bores, refs)?;
    let selected = source.skin.resolve_faces(refs, dec)?;
    let chosen = selected.iter().map(|i| source.planes[*i]).collect();
    if !source.skin.is_convex(dec)? {
        return Err(unsupported(
            "a preserved bore needs a convex planar envelope",
        ));
    }
    let cylinders: Vec<_> = bores.iter().map(|bore| &bore.cylinder).collect();
    let caps = bored::certify_bores(
        source.skin.patches.iter().map(|patch| &patch.plane),
        &cylinders,
        dec,
    )?;
    if caps.iter().any(|cap| selected.contains(cap)) {
        return Err(unsupported("a selected bore cap moves the cap-hole trim"));
    }
    let outer = convex_cell(&source.skin, &selected, definition, dec)?;
    let cell = BoredCell::new(outer, bores, definition, dec)?;
    let bound = bind_bores(model, &source.bores, &cell)?;
    Ok((chosen, Cell::Bored(Box::new(cell), bound)))
}

/// The planar base under an unchanged boss and around any unchanged bore, proved
/// before any topology.
fn bossed_cell(
    model: &Model,
    source: &BossedSource,
    refs: &[FaceRef],
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<(BTreeSet<FaceId>, Cell), Refusal> {
    let unsupported = |detail| Refusal::new("draft.unsupported-arrangement", detail);
    let [shoulder, end] = &source.rings;
    let origins = model.face_origins(source.wall).to_vec();
    let top_origins = model.face_origins(end.cap).to_vec();
    if refs.iter().any(|r| {
        origins
            .iter()
            .chain(&top_origins)
            .any(|origin| r.face.matches(origin))
    }) {
        return Err(unsupported("a selected boss needs a curved arrangement"));
    }
    let bores = preserved_bores(model, &source.bores, refs)?;
    let selected = source.skin.resolve_faces(refs, dec)?;
    let chosen = selected.iter().map(|i| source.planes[*i]).collect();
    if !source.skin.is_convex(dec)? {
        return Err(unsupported(
            "a preserved boss needs a convex planar envelope",
        ));
    }
    let top = model.face(end.cap);
    let (Surface::Cylinder(cylinder), Surface::Plane(plane)) = (
        model.geom().surface(model.face(source.wall).surface),
        model.geom().surface(top.surface),
    ) else {
        unreachable!("bossed sources keep a cylindrical wall and a planar top")
    };
    let mut plane = plane.clone();
    if top.sense == Sense::Reversed {
        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
    }
    let boss = Boss {
        cylinder: cylinder.clone(),
        origins,
        top: plane,
        top_origins,
    };
    let cap = bossed::certify_footprint(
        source.skin.patches.iter().map(|patch| &patch.plane),
        &boss,
        dec,
    )?;
    if source.planes[cap] != shoulder.cap {
        return Err(unsupported("a boss ring lost its cap support"));
    }
    if selected.contains(&cap) {
        return Err(unsupported("a selected shoulder moves the boss's ring"));
    }
    if !bores.is_empty() {
        let cylinders: Vec<_> = bores.iter().map(|bore| &bore.cylinder).collect();
        let caps = bored::certify_bores(
            source.skin.patches.iter().map(|patch| &patch.plane),
            &cylinders,
            dec,
        )?;
        if caps.iter().any(|cap| selected.contains(cap)) {
            return Err(unsupported("a selected bore cap moves the cap-hole trim"));
        }
    }
    let outer = convex_cell(&source.skin, &selected, definition, dec)?;
    let cell = BossedCell::new(outer, boss, bores, definition, dec)?;
    let bound = match &cell.bored {
        Some(bored) => bind_bores(model, &source.bores, bored)?,
        None => Vec::new(),
    };
    Ok((
        chosen,
        Cell::Bossed(Box::new(cell), Box::new(source.rings.clone()), bound),
    ))
}

fn ring_edge(model: &mut Model, ring: &BoreRing) -> EdgeId {
    let curve = model.geom_mut().add_curve(ring.curve.clone());
    model.add_edge(Edge {
        curve,
        sense: ring.edge_sense,
        start: None,
        end: None,
    })
}

/// Each bore's wall, and its rings as holes in the two caps.
fn bore_faces(
    model: &mut Model,
    bores: &[Bore],
    bound: &[BoundBore],
) -> (Vec<HolePath>, Vec<WallFace>) {
    let mut holes = Vec::with_capacity(2 * bound.len());
    let mut walls = Vec::with_capacity(bound.len());
    for (bore, bound) in bores.iter().zip(bound) {
        let mut paths = Vec::with_capacity(2);
        for (ring, support) in &bound.rings {
            let edge = ring_edge(model, ring);
            holes.push((*support, vec![(edge, ring.cap_sense)]));
            paths.push(vec![(edge, ring.wall_sense)]);
        }
        walls.push((
            Surface::Cylinder(bore.cylinder.clone()),
            bound.sense,
            paths,
            bore.origins.clone(),
        ));
    }
    (holes, walls)
}

fn assemble_bossed(
    model: &mut Model,
    cell: &BossedCell,
    rings: &[BoreRing; 2],
    bound: &[BoundBore],
) -> BodyId {
    let [lower, upper] = rings.each_ref().map(|ring| ring_edge(model, ring));
    let [shoulder, end] = rings;
    let boss = &cell.boss;
    let bores = cell.bored.as_ref().map_or(&[][..], |bored| &bored.bores);
    let (mut holes, bore_walls) = bore_faces(model, bores, bound);
    let mut walls = vec![
        (
            Surface::Cylinder(boss.cylinder.clone()),
            Sense::Same,
            vec![
                vec![(lower, shoulder.wall_sense)],
                vec![(upper, end.wall_sense)],
            ],
            boss.origins.clone(),
        ),
        (
            Surface::Plane(boss.top.clone()),
            Sense::Same,
            vec![vec![(upper, end.cap_sense)]],
            boss.top_origins.clone(),
        ),
    ];
    walls.extend(bore_walls);
    holes.push((cell.shoulder, vec![(lower, shoulder.cap_sense)]));
    assemble(model, &cell.outer, &holes, walls)
}

fn assemble_bored(model: &mut Model, cell: &BoredCell, bound: &[BoundBore]) -> BodyId {
    let (holes, walls) = bore_faces(model, &cell.bores, bound);
    assemble(model, &cell.outer, &holes, walls)
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

type HolePath = (usize, Vec<(EdgeId, Sense)>);
type WallFace = (Surface, Sense, Vec<Vec<(EdgeId, Sense)>>, Vec<Origin>);

fn assemble(
    model: &mut Model,
    cell: &PlanarCell,
    holes: &[HolePath],
    walls: Vec<WallFace>,
) -> BodyId {
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
    for (index, (support, cycle)) in cell.supports.iter().zip(&cell.cycles).enumerate() {
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
        let mut paths = vec![path];
        paths.extend(
            holes
                .iter()
                .filter(|(cap, _)| *cap == index)
                .map(|(_, hole)| hole.clone()),
        );
        faces.push(crate::extrude::face(
            model,
            Surface::Plane(support.plane.clone()),
            Sense::Same,
            paths,
            None,
            support.origins.clone(),
        ));
    }
    for (surface, sense, paths, origins) in walls {
        faces.push(crate::extrude::face(
            model, surface, sense, paths, None, origins,
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
        let skin = planar_skin(model, &source_faces, &BTreeSet::new())?;
        let selected = skin.resolve_faces(refs, dec)?;
        let chosen = selected.iter().map(|i| source_faces[*i]).collect();
        (chosen, planar(&skin, &selected, definition, dec)?)
    } else if let Some(source) = bored_source(model, &source_faces)? {
        let (chosen, cell) = bored_cell(model, &source, refs, definition, dec)?;
        (chosen, vec![cell])
    } else if let Some(source) = bossed_source(model, &source_faces)? {
        let (chosen, cell) = bossed_cell(model, &source, refs, definition, dec)?;
        (chosen, vec![cell])
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
