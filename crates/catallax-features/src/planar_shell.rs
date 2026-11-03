//! Finite planar shell carriers classified together inside the immutable source.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::{Containment, PatchSelection};
use catallax_numbers::{Decider, MarkLog, Outcome, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, FaceRef, Origin, Role, StepId, VertexName};
use catallax_shapes::distance::PlanarBoundary;
use catallax_shapes::{
    Circle, Curve, Cylinder, Line, Plane, Point3, PointDef, Sphere, Surface, Vector3,
};
use catallax_topology::{BodyId, Edge, Face, FaceId, Model, Sense, Vertex};

use crate::{draft::arrangement::Skin, profile::sign};

#[path = "planar_shell/nearest.rs"]
mod nearest;

struct EdgeIncident {
    face: usize,
    lateral: Vector3,
}
struct VertexIncident {
    face: usize,
    toward: [Vector3; 2],
    convex: bool,
}
fn winding(skin: &Skin, face: usize, dec: &mut Decider) -> Result<Real, Refusal> {
    let patch = &skin.patches[face];
    let cycle = &patch.loops[0];
    if cycle.len() < 3 {
        return Err(bad("a planar nearest region needs a polygonal outer loop"));
    }
    let origin = &skin.vertices[cycle[0]];
    let mut area = Real::zero();
    for i in 1..cycle.len() - 1 {
        area = area
            + skin.vertices[cycle[i]]
                .minus(origin)
                .cross(&skin.vertices[cycle[i + 1]].minus(origin))
                .dot(&patch.plane.normal);
    }
    match sign(&area, dec)? {
        Sign::Positive => Ok(Real::one()),
        Sign::Negative => Ok(Real::from_i64(-1)),
        Sign::Zero => Err(bad("a physical planar face has zero outer-loop area")),
    }
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("distance.unsupported-envelope", detail)
}

#[derive(Default)]
struct Sides {
    minus: bool,
    plus: bool,
    tangent: bool,
}
impl Sides {
    fn admit(
        &mut self,
        delta: &Vector3,
        direction: &Vector3,
        squared: &Real,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        match sign(&(&delta.norm_squared() - squared), dec)? {
            Sign::Negative => {
                self.minus = true;
                self.plus = true;
            }
            Sign::Positive => {}
            Sign::Zero => match sign(&delta.dot(direction), dec)? {
                Sign::Negative => self.plus = true,
                Sign::Positive => self.minus = true,
                Sign::Zero => self.tangent = true,
            },
        }
        Ok(())
    }
}

struct FiniteDistance {
    faces: Vec<PlanarBoundary>,
    edges: Vec<[Point3; 2]>,
    squared: Real,
}
impl FiniteDistance {
    fn sides(
        &self,
        point: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<(bool, bool)>, Refusal> {
        let mut sides = Sides::default();
        for face in &self.faces {
            let height = point.minus(&face.plane.origin).dot(&face.plane.normal);
            let projection = point.plus(
                &face
                    .plane
                    .normal
                    .scaled(&(-height / face.plane.normal.norm_squared())),
            );
            if face.contains(&projection, dec)? {
                sides.admit(&point.minus(&projection), direction, &self.squared, dec)?;
                if sides.minus && sides.plus {
                    return Ok(Some((true, true)));
                }
            }
        }
        for [a, b] in &self.edges {
            let edge = b.minus(a);
            let along = point.minus(a).dot(&edge);
            let length_squared = edge.norm_squared();
            let foot = if sign(&along, dec)? != Sign::Positive {
                a.clone()
            } else if sign(&(&along - &length_squared), dec)? != Sign::Negative {
                b.clone()
            } else {
                a.plus(&edge.scaled(&(along / length_squared)))
            };
            sides.admit(&point.minus(&foot), direction, &self.squared, dec)?;
            if sides.minus && sides.plus {
                return Ok(Some((true, true)));
            }
        }
        // A zero first derivative needs another strict cell witness. Never
        // choose a side by a floating perturbation or an enlarged tolerance.
        Ok((!sides.tangent).then_some((sides.minus, sides.plus)))
    }
}

fn normal(model: &Model, face: FaceId, point: &Point3) -> Vector3 {
    let face = model.face(face);
    let normal = match model.geom().surface(face.surface) {
        Surface::Plane(plane) => plane.normal.clone(),
        Surface::Cylinder(cylinder) => {
            let delta = point.minus(&cylinder.origin);
            delta.minus(
                &cylinder
                    .axis
                    .scaled(&(delta.dot(&cylinder.axis) / cylinder.axis.norm_squared())),
            )
        }
        Surface::Sphere(sphere) => point.minus(&sphere.center),
        _ => unreachable!("finite planar carriers are planes, cylinders or spheres"),
    };
    if face.sense == Sense::Same {
        normal
    } else {
        normal.scaled(&Real::from_i64(-1))
    }
}

fn parallel_face(
    model: &mut Model,
    skin: &Skin,
    index: usize,
    thickness: &Real,
    direction: i64,
    step: &StepId,
) -> FaceId {
    let patch = &skin.patches[index];
    let scale =
        &(thickness / &patch.plane.normal.norm_squared().sqrt()) * &Real::from_i64(direction);
    let shift = patch.plane.normal.scaled(&scale);
    let mut plane = patch.plane.clone();
    plane.origin = plane.origin.plus(&shift);
    plane.normal = plane.normal.scaled(&Real::from_i64(direction));
    let mut vertices = BTreeMap::new();
    let mut paths = Vec::new();
    for cycle in &patch.loops {
        let order: Vec<_> = if direction > 0 {
            cycle.clone()
        } else {
            cycle.iter().rev().copied().collect()
        };
        let mut path = Vec::with_capacity(order.len());
        for i in 0..order.len() {
            let a = order[i];
            let b = order[(i + 1) % order.len()];
            for vertex in [a, b] {
                vertices.entry(vertex).or_insert_with(|| {
                    let point = model
                        .geom_mut()
                        .add_point(PointDef::Explicit(skin.vertices[vertex].plus(&shift)));
                    model.add_vertex(Vertex { point })
                });
            }
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                origin: skin.vertices[a].plus(&shift),
                dir: skin.vertices[b].minus(&skin.vertices[a]),
            }));
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(vertices[&a]),
                end: Some(vertices[&b]),
            });
            path.push((edge, Sense::Same));
        }
        paths.push(path);
    }
    crate::extrude::face(
        model,
        Surface::Plane(plane),
        Sense::Same,
        paths,
        None,
        patch
            .origins
            .iter()
            .map(|of| Origin {
                step: step.clone(),
                role: Role::Inner {
                    of: Box::new(of.clone()),
                },
            })
            .collect(),
    )
}

fn cylinder_face(
    model: &mut Model,
    cylinder: Cylinder,
    ends: [&Point3; 2],
    origins: Vec<Origin>,
) -> FaceId {
    let rings = ends.map(|center| {
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: center.clone(),
            normal: cylinder.axis.clone(),
            radius_sq: cylinder.radius_sq.clone(),
            ref_dir: cylinder.ref_dir.clone(),
        }));
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        })
    });
    crate::extrude::face(
        model,
        Surface::Cylinder(cylinder),
        Sense::Same,
        vec![
            vec![(rings[0], Sense::Same)],
            vec![(rings[1], Sense::Reversed)],
        ],
        None,
        origins,
    )
}

pub(super) fn build(
    model: &mut Model,
    source: BodyId,
    open: &[FaceRef],
    thickness: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let source_faces = model.faces_of_body(source);
    let skin = crate::draft_build::planar_skin(model, &source_faces, &BTreeSet::new())
        .map_err(|error| bad(error.detail))?;
    let selected = skin.resolve_faces(open, dec)?;
    if selected.len() == skin.patches.len() {
        return Ok(catallax_topology::empty(model));
    }
    let kept = |i: usize| !selected.contains(&i);
    let mut source_cells = nearest::Source::new(&skin, dec)?;
    let mut trim_marks = MarkLog::new();
    let mut candidates = Vec::new();
    let mut edges: BTreeMap<(usize, usize), Vec<EdgeIncident>> = BTreeMap::new();
    let mut vertices: BTreeMap<usize, Vec<VertexIncident>> = BTreeMap::new();
    let mut distance_faces = Vec::new();
    for (i, patch) in skin.patches.iter().enumerate() {
        if kept(i) {
            distance_faces.push(PlanarBoundary::polygon(
                patch.plane.clone(),
                &skin.vertices,
                &patch.loops,
            ));
            for direction in [-1, 1] {
                candidates.push(parallel_face(model, &skin, i, thickness, direction, step));
            }
        }
        let sense = winding(&skin, i, dec)?;
        for cycle in &patch.loops {
            if cycle.len() < 3 {
                return Err(bad("a planar nearest region needs closed polygonal loops"));
            }
            for j in 0..cycle.len() {
                let a = cycle[j];
                let b = cycle[(j + 1) % cycle.len()];
                let previous = cycle[(j + cycle.len() - 1) % cycle.len()];
                let forward = skin.vertices[b].minus(&skin.vertices[a]);
                let backward = skin.vertices[previous].minus(&skin.vertices[a]);
                let lateral = forward.cross(&patch.plane.normal).scaled(&sense);
                edges
                    .entry((a.min(b), a.max(b)))
                    .or_default()
                    .push(EdgeIncident { face: i, lateral });
                let turn = backward
                    .scaled(&Real::from_i64(-1))
                    .cross(&forward)
                    .dot(&patch.plane.normal)
                    * sense.clone();
                let convex = sign(&turn, dec)? == Sign::Positive;
                vertices.entry(a).or_default().push(VertexIncident {
                    face: i,
                    toward: [backward, forward],
                    convex,
                });
            }
        }
    }
    let mut distance_edges = Vec::new();
    for ((a, b), incident) in &edges {
        let [first, second] = incident.as_slice() else {
            return Err(bad("a physical source edge needs two incident faces"));
        };
        if !kept(first.face) && !kept(second.face) {
            continue;
        }
        let ends = [&skin.vertices[*a], &skin.vertices[*b]];
        let axis = ends[1].minus(ends[0]);
        let reference = &skin.patches[first.face].plane.normal;
        let mut constraints = Vec::with_capacity(4);
        constraints.push(Plane {
            origin: ends[0].clone(),
            normal: axis.scaled(&Real::from_i64(-1)),
            ref_dir: reference.clone(),
        });
        constraints.push(Plane {
            origin: ends[1].clone(),
            normal: axis.clone(),
            ref_dir: reference.clone(),
        });
        for face in incident.iter().filter(|face| kept(face.face)) {
            constraints.push(Plane {
                origin: ends[0].clone(),
                normal: face.lateral.scaled(&Real::from_i64(-1)),
                ref_dir: axis.clone(),
            });
        }
        let clipped = source_cells.clip(
            model,
            &constraints,
            |model, dec| {
                let mut origins = Vec::new();
                for a in &skin.patches[first.face].origins {
                    for b in &skin.patches[second.face].origins {
                        let edge = EdgeName::new(a.clone(), b.clone())
                            .map_err(|error| bad(error.to_string()))?;
                        origins.push(Origin {
                            step: step.clone(),
                            role: Role::Round { edge },
                        });
                    }
                }
                if origins.is_empty() {
                    return Err(bad("a physical edge envelope needs named source faces"));
                }
                let cylinder = Cylinder::new(
                    ends[0].clone(),
                    axis,
                    thickness.clone(),
                    reference.clone(),
                    dec,
                )?;
                Ok(cylinder_face(model, cylinder, ends, origins))
            },
            dec,
        )?;
        candidates.extend(clipped.value);
        trim_marks.extend(&clipped.marks);
        // The exact minimum still contains every kept physical edge, even
        // when its curved boundary stratum is absent from the final source.
        distance_edges.push([ends[0].clone(), ends[1].clone()]);
    }
    for (vertex, incident) in vertices {
        let mut active = 0;
        let mut excluded = false;
        for face in &incident {
            if kept(face.face) {
                active += 1;
                excluded |= !face.convex;
            }
        }
        if active == 0 || excluded {
            continue;
        }
        let mut constraints = Vec::with_capacity(active * 2);
        for face in incident.iter().filter(|face| kept(face.face)) {
            for toward in &face.toward {
                constraints.push(Plane {
                    origin: skin.vertices[vertex].clone(),
                    normal: toward.clone(),
                    ref_dir: skin.patches[face.face].plane.normal.clone(),
                });
            }
        }
        let clipped = source_cells.clip(
            model,
            &constraints,
            |model, dec| {
                let mut origins: Vec<_> = incident
                    .iter()
                    .flat_map(|face| skin.patches[face.face].origins.iter().cloned())
                    .collect();
                origins.sort();
                origins.dedup();
                let vertex_name =
                    VertexName::new(origins).map_err(|error| bad(error.to_string()))?;
                let frame = &skin.patches[incident[0].face].plane;
                let sphere = Sphere::new(
                    skin.vertices[vertex].clone(),
                    thickness.clone(),
                    frame.normal.clone(),
                    frame.ref_dir.clone(),
                    dec,
                )?;
                let surface = model.geom_mut().add_surface(Surface::Sphere(sphere));
                Ok(model.add_named_face(
                    Face {
                        surface,
                        sense: Sense::Same,
                        loops: Vec::new(),
                    },
                    vec![Origin {
                        step: step.clone(),
                        role: Role::RoundCorner {
                            vertex: vertex_name,
                        },
                    }],
                ))
            },
            dec,
        )?;
        candidates.extend(clipped.value);
        trim_marks.extend(&clipped.marks);
    }
    let source = if trim_marks.is_empty() {
        source
    } else {
        let mut marked = model.body(source).clone();
        marked.marks.extend(&trim_marks);
        model.add_body(marked)
    };
    let distance = FiniteDistance {
        faces: distance_faces,
        edges: distance_edges,
        squared: thickness.square(),
    };
    let source_ids: BTreeSet<_> = source_faces.iter().copied().collect();
    catallax_boolean::arrange(
        model,
        source,
        &candidates,
        |model, owner, point, dec| {
            let source_patch = source_ids.contains(&owner);
            if !source_patch {
                match catallax_boolean::classify_point(model, source, point, dec) {
                    Outcome::Proven(Containment::Outside) => {
                        return Ok(Some(PatchSelection::Discard));
                    }
                    Outcome::Proven(Containment::Inside) => {}
                    Outcome::ProvenCoincident(_) => {
                        for &face in &source_faces {
                            match model.geom().same_surface(
                                model.face(owner).surface,
                                model.face(face).surface,
                                dec,
                            ) {
                                Outcome::ProvenCoincident(_) => {
                                    if catallax_boolean::face_contains(model, face, point, dec)? {
                                        return Ok(Some(PatchSelection::Discard));
                                    }
                                }
                                Outcome::Refused(error) => return Err(error),
                                _ => {}
                            }
                        }
                        return Ok(None);
                    }
                    Outcome::NotExact(error) => return Err(bad(error.detail)),
                    Outcome::Refused(error) => return Err(error),
                }
            }
            let direction = normal(model, owner, point);
            let Some((minus, plus)) = distance.sides(point, &direction, dec)? else {
                return Ok(None);
            };
            Ok(Some(if source_patch {
                if minus {
                    PatchSelection::Keep(Sense::Same)
                } else {
                    PatchSelection::Discard
                }
            } else {
                match (minus, plus) {
                    (true, false) => PatchSelection::Keep(Sense::Same),
                    (false, true) => PatchSelection::Keep(Sense::Reversed),
                    _ => PatchSelection::Discard,
                }
            }))
        },
        dec,
    )
    .map(|answer| answer.body)
    .map_err(|error| error.refusal)
}
