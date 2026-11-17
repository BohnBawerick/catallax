use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::Refusal;
use catallax_shapes::{Curve, Point3, Surface};
use catallax_topology::{EdgeId, FaceId, Loop, Model, VertexId};

use crate::{EdgeMesh, FaceMesh, Mesh, MeshOptions, Triangle, edges::Boundaries, failure};

type Chords = BTreeSet<(usize, usize)>;

pub(crate) struct Previous<'a> {
    model: &'a Model,
    mesh: &'a Mesh,
    edges: BTreeMap<EdgeId, &'a EdgeMesh>,
    chords: BTreeMap<(usize, usize), usize>,
}

fn same_point(a: &Point3, b: &Point3) -> bool {
    a.x.same_definition(&b.x) && a.y.same_definition(&b.y) && a.z.same_definition(&b.z)
}

fn same_vertex(a: &Model, av: Option<VertexId>, b: &Model, bv: Option<VertexId>) -> bool {
    match (av, bv) {
        (Some(av), Some(bv)) => same_point(
            &a.geom().position(a.vertex(av).point),
            &b.geom().position(b.vertex(bv).point),
        ),
        (None, None) => true,
        _ => false,
    }
}

fn same_surface(a: &Surface, b: &Surface) -> bool {
    match (a, b) {
        // A spline carrier alone does not identify its finite parameter cover.
        (Surface::BSpline(a), Surface::BSpline(b)) => a.shares_definition_nodes(b),
        _ => a.same_definition(b),
    }
}

fn same_curve(a: &Model, ae: EdgeId, b: &Model, be: EdgeId) -> bool {
    let ac = a.geom().curve(a.edge(ae).curve);
    let bc = b.geom().curve(b.edge(be).curve);
    if !ac.same_definition(bc) {
        return false;
    }
    match (ac, bc) {
        // Intersection construction keys contain model-local surface handles.
        (Curve::Intersection(ac), Curve::Intersection(bc)) => ac
            .surfaces
            .iter()
            .zip(bc.surfaces.iter())
            .all(|(&a_id, &b_id)| same_surface(a.geom().surface(a_id), b.geom().surface(b_id))),
        _ => true,
    }
}

impl<'a> Previous<'a> {
    pub(crate) fn new(model: &'a Model, mesh: &'a Mesh) -> Self {
        let mut boundary = vec![false; mesh.positions.len()];
        for edge in &mesh.edges {
            for &vertex in &edge.vertices {
                boundary[vertex] = true;
            }
        }
        let mut owner = vec![usize::MAX; mesh.positions.len()];
        for (index, face) in mesh.faces.iter().enumerate() {
            for triangle in &mesh.triangles[face.triangles.clone()] {
                for vertex in triangle.vertices {
                    if owner[vertex] != usize::MAX && owner[vertex] != index {
                        boundary[vertex] = true;
                    }
                    owner[vertex] = index;
                }
            }
        }
        let mut chords = BTreeMap::new();
        for (index, face) in mesh.faces.iter().enumerate() {
            let mut once = Chords::new();
            for triangle in &mesh.triangles[face.triangles.clone()] {
                for i in 0..3 {
                    let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                    if boundary[a] && boundary[b] {
                        let pair = (a.min(b), a.max(b));
                        if !once.insert(pair) {
                            chords.entry(pair).or_insert(index);
                        }
                    }
                }
            }
        }
        Self {
            model,
            chords,
            mesh,
            edges: mesh.edges.iter().map(|edge| (edge.edge, edge)).collect(),
        }
    }

    fn boundary_map(
        &self,
        old: &FaceMesh,
        model: &Model,
        face: FaceId,
        boundaries: &Boundaries,
    ) -> Option<BTreeMap<usize, usize>> {
        let before = self.model.face(old.face);
        let after = model.face(face);
        if before.sense != after.sense
            || before.loops.len() != after.loops.len()
            || !same_surface(
                self.model.geom().surface(before.surface),
                model.geom().surface(after.surface),
            )
        {
            return None;
        }
        let mut remap = BTreeMap::new();
        let mut reverse = BTreeMap::new();
        let mut bind = |old: usize, new: usize| -> Option<()> {
            if self.mesh.positions.get(old)? != boundaries.mesh.positions.get(new)?
                || remap.get(&old).is_some_and(|&value| value != new)
                || reverse.get(&new).is_some_and(|&value| value != old)
            {
                return None;
            }
            remap.insert(old, new);
            reverse.insert(new, old);
            Some(())
        };
        for (&old_loop, &new_loop) in before.loops.iter().zip(&after.loops) {
            match (self.model.loop_(old_loop), model.loop_(new_loop)) {
                (Loop::Vertex(a), Loop::Vertex(b)) => {
                    if !same_vertex(self.model, Some(*a), model, Some(*b)) {
                        return None;
                    }
                    let new = boundaries.vertices[b];
                    let position = boundaries.mesh.positions[new];
                    // Singleton trims have no EdgeMesh entry. Refuse ambiguous
                    // coordinate aliases instead of welding distinct vertices.
                    let mut candidates = self
                        .mesh
                        .positions
                        .iter()
                        .enumerate()
                        .filter(|(_, point)| **point == position);
                    let old = candidates.next()?.0;
                    if candidates.next().is_some() {
                        return None;
                    }
                    bind(old, new)?;
                }
                (Loop::Coedges(a), Loop::Coedges(b)) if a.len() == b.len() => {
                    for (&a, &b) in a.iter().zip(b) {
                        let a = self.model.coedge(a);
                        let b = model.coedge(b);
                        let ae = self.model.edge(a.edge);
                        let be = model.edge(b.edge);
                        if a.sense != b.sense
                            || ae.sense != be.sense
                            || !same_curve(self.model, a.edge, model, b.edge)
                            || !same_vertex(self.model, ae.start, model, be.start)
                            || !same_vertex(self.model, ae.end, model, be.end)
                        {
                            return None;
                        }
                        let old = self.edges.get(&a.edge)?;
                        let new = &boundaries.mesh.edges[boundaries.edges[&b.edge]];
                        if old.vertices.len() != new.vertices.len()
                            || old.max_chord_error != new.max_chord_error
                        {
                            return None;
                        }
                        for (&old, &new) in old.vertices.iter().zip(&new.vertices) {
                            bind(old, new)?;
                        }
                    }
                }
                _ => return None,
            }
        }
        Some(remap)
    }

    fn prior_matches(
        &self,
        face_index: usize,
        remap: &BTreeMap<usize, usize>,
        prior: &Chords,
    ) -> bool {
        let old_prior: Chords = self
            .chords
            .iter()
            .filter(|&(_, &first_face)| first_face < face_index)
            .filter_map(|(&(a, b), _)| {
                let (&a, &b) = (remap.get(&a)?, remap.get(&b)?);
                Some((a.min(b), a.max(b)))
            })
            .collect();
        let vertices: BTreeSet<_> = remap.values().copied().collect();
        old_prior.iter().copied().eq(prior
            .iter()
            .copied()
            .filter(|(a, b)| vertices.contains(a) && vertices.contains(b)))
    }

    pub(crate) fn append(
        &self,
        model: &Model,
        face: FaceId,
        boundaries: &mut Boundaries,
        prior: &Chords,
        options: &MeshOptions,
    ) -> Result<bool, Refusal> {
        for (index, old) in self.mesh.faces.iter().enumerate() {
            let Some(mut remap) = self.boundary_map(old, model, face, boundaries) else {
                continue;
            };
            if !self.prior_matches(index, &remap, prior) {
                continue;
            }
            let triangles = &self.mesh.triangles[old.triangles.clone()];
            let interior: BTreeSet<_> = triangles
                .iter()
                .flat_map(|triangle| triangle.vertices)
                .filter(|vertex| !remap.contains_key(vertex))
                .collect();
            if interior.len()
                > options
                    .max_vertices
                    .saturating_sub(boundaries.mesh.positions.len())
            {
                return Err(failure("reused face triangulation reached max_vertices"));
            }
            // Preserve the original position order, including periodic chart
            // copies that survived compaction. Only boundary indices are shared.
            for vertex in interior {
                let next = boundaries.mesh.positions.len();
                boundaries.mesh.positions.push(self.mesh.positions[vertex]);
                remap.insert(vertex, next);
            }
            let start = boundaries.mesh.triangles.len();
            boundaries
                .mesh
                .triangles
                .extend(triangles.iter().map(|triangle| Triangle {
                    vertices: triangle.vertices.map(|vertex| remap[&vertex]),
                    face,
                }));
            let mut report = old.clone();
            report.face = face;
            report.triangles = start..boundaries.mesh.triangles.len();
            boundaries.mesh.faces.push(report);
            boundaries.mesh.unavailable_controls.extend(
                self.mesh
                    .unavailable_controls
                    .iter()
                    .filter(|control| control.face == old.face)
                    .map(|control| {
                        let mut control = control.clone();
                        control.face = face;
                        control
                    }),
            );
            return Ok(true);
        }
        Ok(false)
    }
}
