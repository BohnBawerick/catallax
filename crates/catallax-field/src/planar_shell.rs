//! Finite-face distance tubes from independently reconstructed planar definitions.

use std::collections::BTreeMap;

use catallax_features::{
    draft::arrangement::{Patch, Skin},
    polyhedron::{ConvexCell, PlaneFace},
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, EdgeName, FaceRef, Origin, Role, StepId, VertexName};
use catallax_shapes::{Cylinder, Plane, Point3, Sphere, Surface, Vector3};

use super::bad;
use crate::{Plan, geometry};

fn cache_rational(value: &mut Real, dec: &mut Decider) -> Result<(), Refusal> {
    if value.rational_ref().is_none()
        && let Some(rational) = dec.rational_value(value)?
    {
        *value = Real::from(rational);
    }
    Ok(())
}
fn cache_point_rationals(point: &mut Point3, dec: &mut Decider) -> Result<(), Refusal> {
    for value in [&mut point.x, &mut point.y, &mut point.z] {
        cache_rational(value, dec)?;
    }
    Ok(())
}
fn cache_vector_rationals(vector: &mut Vector3, dec: &mut Decider) -> Result<(), Refusal> {
    for value in [&mut vector.x, &mut vector.y, &mut vector.z] {
        cache_rational(value, dec)?;
    }
    Ok(())
}
fn cache_plane_rationals(plane: &mut Plane, dec: &mut Decider) -> Result<(), Refusal> {
    cache_point_rationals(&mut plane.origin, dec)?;
    cache_vector_rationals(&mut plane.normal, dec)?;
    cache_vector_rationals(&mut plane.ref_dir, dec)
}
fn cache_skin_rationals(skin: &mut Skin, dec: &mut Decider) -> Result<(), Refusal> {
    for point in &mut skin.vertices {
        cache_point_rationals(point, dec)?;
    }
    for patch in &mut skin.patches {
        cache_plane_rationals(&mut patch.plane, dec)?;
    }
    Ok(())
}

fn face_tube(skin: &Skin, patch: &Patch, thickness: &Real, step: &StepId) -> Skin {
    let shift = patch
        .plane
        .normal
        .scaled(&(thickness / &patch.plane.normal.norm_squared().sqrt()));
    let backwards = shift.scaled(&Real::from_i64(-1));
    let count: usize = patch.loops.iter().map(Vec::len).sum();
    let mut vertices = Vec::with_capacity(2 * count);
    let mut indices = BTreeMap::new();
    for &vertex in patch.loops.iter().flatten() {
        indices.entry(vertex).or_insert_with(|| {
            let index = vertices.len();
            vertices.push(skin.vertices[vertex].plus(&backwards));
            vertices.push(skin.vertices[vertex].plus(&shift));
            index
        });
    }
    let mut origins: Vec<_> = patch
        .origins
        .iter()
        .map(|of| Origin {
            step: step.clone(),
            role: Role::Inner { of: Box::new(of.clone()) },
        })
        .collect();
    let mut patches = Vec::with_capacity(count + 2);
    for side in [0, 1] {
        let direction = Real::from_i64(if side == 0 { -1 } else { 1 });
        let mut plane = patch.plane.clone();
        plane.origin = plane.origin.plus(if side == 0 { &backwards } else { &shift });
        plane.normal = plane.normal.scaled(&direction);
        let loops = patch
            .loops
            .iter()
            .map(|cycle| {
                if side == 0 {
                    cycle.iter().rev().map(|v| indices[v]).collect()
                } else {
                    cycle.iter().map(|v| indices[v] + 1).collect()
                }
            })
            .collect();
        let cap_origins = if side == 0 { origins.clone() } else { std::mem::take(&mut origins) };
        patches.push(Patch { plane, loops, origins: cap_origins });
    }
    for cycle in &patch.loops {
        for i in 0..cycle.len() {
            let a = cycle[i];
            let b = cycle[(i + 1) % cycle.len()];
            patches.push(Patch {
                plane: Plane {
                    origin: vertices[indices[&a]].clone(),
                    normal: skin.vertices[b].minus(&skin.vertices[a]).cross(&patch.plane.normal),
                    ref_dir: patch.plane.normal.clone(),
                },
                loops: vec![vec![indices[&a], indices[&b], indices[&b] + 1, indices[&a] + 1]],
                origins: Vec::new(),
            });
        }
    }
    Skin { vertices, patches }
}

fn tube_leaf(plan: &mut Plan, mut tube: Skin, bound: &Rational, dec: &mut Decider) -> Result<usize, Refusal> {
    cache_skin_rationals(&mut tube, dec)?;
    if tube.is_convex(dec)? {
        let faces = tube
            .patches
            .into_iter()
            .map(|patch| PlaneFace { plane: patch.plane, origins: patch.origins })
            .collect();
        return plan.distance_cell(ConvexCell::new(faces, dec)?, bound.clone(), dec);
    }
    let mut result = 0;
    for component in tube.shell_field_components(dec)? {
        for term in component.terms {
            let faces = term
                .sides
                .into_iter()
                .map(|(support, same)| {
                    let mut plane = component.planes[support].clone();
                    if !same {
                        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
                    }
                    cache_plane_rationals(&mut plane, dec)?;
                    let origins = component
                        .faces
                        .iter()
                        .filter(|face| face.support == support && face.same_sense == same)
                        .flat_map(|face| face.origins.iter().cloned())
                        .collect();
                    Ok(PlaneFace { plane, origins })
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            let part = plan.distance_cell(ConvexCell::new(faces, dec)?, bound.clone(), dec)?;
            result = plan.boolean(BooleanKind::Union, result, part, dec)?;
        }
    }
    Ok(result)
}

struct EdgeIncident {
    face: usize,
    outward: Vector3,
}

struct VertexIncident {
    face: usize,
    toward: Option<[Vector3; 2]>,
}

fn region_leaf(
    plan: &mut Plan,
    surface: Surface,
    origins: Vec<Origin>,
    planes: Vec<Plane>,
    bound: &Rational,
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    let faces = origins.into_iter()
        .map(|origin| geometry::face(surface.clone(), origin, dec))
        .collect::<Result<Vec<_>, _>>()?;
    let mut surfaces = Vec::with_capacity(planes.len() + 1);
    surfaces.push(surface);
    for mut plane in planes {
        cache_plane_rationals(&mut plane, dec)?;
        surfaces.push(Surface::Plane(plane));
    }
    plan.leaf(
        geometry::Leaf::Cell(geometry::Cell::new(surfaces, Some(bound.clone()))?),
        faces, dec,
    )
}

struct Vertex {
    faces: Vec<VertexIncident>,
    kept: bool,
    convex: bool,
}

impl Plan {
    pub(crate) fn planar_distance_shell(
        &mut self,
        (node, source_bound): (usize, Rational),
        mut skin: Skin,
        open: &[FaceRef],
        thickness: &Rational,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        cache_skin_rationals(&mut skin, dec)?;
        let selected = skin.resolve_faces(open, dec)?;
        if selected.len() == skin.patches.len() {
            return Ok(0);
        }
        let bound = &source_bound + thickness;
        let thickness = Real::from(thickness.clone());
        let mut walls = 0;
        let mut edges: BTreeMap<(usize, usize), Vec<EdgeIncident>> = BTreeMap::new();
        let mut vertices: BTreeMap<usize, Vertex> = BTreeMap::new();
        for (face, patch) in skin.patches.iter().enumerate() {
            let kept = !selected.contains(&face);
            if kept {
                let prism = tube_leaf(self, face_tube(&skin, patch, &thickness, step), &bound, dec)?;
                walls = self.boolean(BooleanKind::Union, walls, prism, dec)?;
            }
            for cycle in &patch.loops {
                for i in 0..cycle.len() {
                    let a = cycle[i];
                    let b = cycle[(i + 1) % cycle.len()];
                    let outgoing = skin.vertices[b].minus(&skin.vertices[a]);
                    edges.entry((a.min(b), a.max(b))).or_default().push(EdgeIncident {
                        face,
                        outward: outgoing.cross(&patch.plane.normal),
                    });
                    let vertex = vertices.entry(a).or_insert_with(|| Vertex {
                        faces: Vec::new(), kept: false, convex: true,
                    });
                    let toward = if kept {
                        vertex.kept = true;
                        let previous = cycle[(i + cycle.len() - 1) % cycle.len()];
                        let incoming = skin.vertices[a].minus(&skin.vertices[previous]);
                        if vertex.convex {
                            vertex.convex = geometry::sign(
                                &incoming.cross(&outgoing).dot(&patch.plane.normal), dec,
                            )? == Sign::Positive;
                        }
                        Some([incoming.scaled(&Real::from_i64(-1)), outgoing])
                    } else {
                        None
                    };
                    vertex.faces.push(VertexIncident { face, toward });
                }
            }
        }
        for ((a, b), incident) in edges {
            let [first, second] = incident.as_slice() else {
                return Err(bad("a physical planar edge needs two incident source faces"));
            };
            if selected.contains(&first.face) && selected.contains(&second.face) {
                continue;
            }
            let mut origins = Vec::with_capacity(
                skin.patches[first.face].origins.len() * skin.patches[second.face].origins.len(),
            );
            for a in &skin.patches[first.face].origins {
                for b in &skin.patches[second.face].origins {
                    let edge = EdgeName::new(a.clone(), b.clone()).map_err(|error| bad(&error.to_string()))?;
                    origins.push(Origin { step: step.clone(), role: Role::Round { edge } });
                }
            }
            let ends = [skin.vertices[a].clone(), skin.vertices[b].clone()];
            let axis = ends[1].minus(&ends[0]);
            let reference = &skin.patches[first.face].plane.normal;
            let mut planes = Vec::with_capacity(incident.len() + 2);
            planes.push(Plane {
                origin: ends[0].clone(),
                normal: axis.scaled(&Real::from_i64(-1)),
                ref_dir: reference.clone(),
            });
            planes.push(Plane {
                origin: ends[1].clone(), normal: axis.clone(), ref_dir: reference.clone(),
            });
            for face in incident.iter().filter(|face| !selected.contains(&face.face)) {
                // These planes bound the angular nearest-edge sector. At reflex
                // edges the finite-face conormals select the opposite quadrant.
                planes.push(Plane {
                    origin: ends[0].clone(),
                    normal: face.outward.scaled(&Real::from_i64(-1)),
                    ref_dir: axis.clone(),
                });
            }
            let mut cylinder = Cylinder::new(
                ends[0].clone(), axis, thickness.clone(), reference.clone(), dec,
            )?;
            cache_point_rationals(&mut cylinder.origin, dec)?;
            cache_vector_rationals(&mut cylinder.axis, dec)?;
            cache_vector_rationals(&mut cylinder.ref_dir, dec)?;
            cache_rational(&mut cylinder.radius_sq, dec)?;
            let part = region_leaf(self, Surface::Cylinder(cylinder), origins, planes, &bound, dec)?;
            walls = self.boolean(BooleanKind::Union, walls, part, dec)?;
        }
        for (index, vertex) in vertices {
            // A flat or reflex kept-face corner has no open vertex-nearest region;
            // its face prism and adjacent edge tubes already cover the entire ball.
            if !vertex.kept || !vertex.convex {
                continue;
            }
            let mut origins: Vec<_> = vertex.faces.iter()
                .flat_map(|face| skin.patches[face.face].origins.iter().cloned()).collect();
            origins.sort();
            origins.dedup();
            let name = VertexName::new(origins).map_err(|error| bad(&error.to_string()))?;
            let frame = &skin.patches[vertex.faces[0].face].plane;
            let mut sphere = Sphere::new(
                skin.vertices[index].clone(), thickness.clone(),
                frame.normal.clone(), frame.ref_dir.clone(), dec,
            )?;
            cache_point_rationals(&mut sphere.center, dec)?;
            cache_vector_rationals(&mut sphere.axis, dec)?;
            cache_vector_rationals(&mut sphere.ref_dir, dec)?;
            cache_rational(&mut sphere.radius, dec)?;
            let mut planes = Vec::with_capacity(vertex.faces.len() * 2);
            for face in vertex.faces {
                if let Some(toward) = face.toward {
                    // The polar cone of the kept finite-face tangent directions
                    // is its vertex-nearest normal cone, including both sides.
                    for normal in toward {
                        planes.push(Plane {
                            origin: skin.vertices[index].clone(), normal,
                            ref_dir: skin.patches[face.face].plane.normal.clone(),
                        });
                    }
                }
            }
            let part = region_leaf(self, Surface::Sphere(sphere), vec![Origin {
                step: step.clone(), role: Role::RoundCorner { vertex: name },
            }], planes, &bound, dec)?;
            walls = self.boolean(BooleanKind::Union, walls, part, dec)?;
        }
        self.boolean(BooleanKind::Intersect, node, walls, dec)
    }
}
