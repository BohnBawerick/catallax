use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Refusal};
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
    pub(crate) fn new(
        model: &'a Model,
        mesh: &'a Mesh,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        // Reserve both vertex-sized initializations before allocating. Each
        // remaining record visit is charged before its scan or copy.
        dec.charge_work((mesh.positions.len() as u64).saturating_mul(2))?;
        let mut boundary = vec![false; mesh.positions.len()];
        for edge in &mesh.edges {
            dec.charge_work(1)?;
            for &vertex in &edge.vertices {
                dec.charge_work(1)?;
                boundary[vertex] = true;
            }
        }
        dec.check_deadline()?;
        let mut owner = vec![usize::MAX; mesh.positions.len()];
        for (index, face) in mesh.faces.iter().enumerate() {
            dec.charge_work(1)?;
            for triangle in &mesh.triangles[face.triangles.clone()] {
                dec.charge_work(3)?;
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
            dec.charge_work(1)?;
            let mut once = Chords::new();
            for triangle in &mesh.triangles[face.triangles.clone()] {
                dec.charge_work(3)?;
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
        let mut edges = BTreeMap::new();
        for edge in &mesh.edges {
            dec.charge_work(1)?;
            edges.insert(edge.edge, edge);
        }
        Ok(Self {
            model,
            chords,
            mesh,
            edges,
        })
    }

    fn boundary_map(
        &self,
        old: &FaceMesh,
        model: &Model,
        face: FaceId,
        boundaries: &Boundaries,
        dec: &mut Decider,
    ) -> Result<Option<BTreeMap<usize, usize>>, Refusal> {
        let before = self.model.face(old.face);
        let after = model.face(face);
        if before.sense != after.sense
            || before.loops.len() != after.loops.len()
            || !same_surface(
                self.model.geom().surface(before.surface),
                model.geom().surface(after.surface),
            )
        {
            return Ok(None);
        }
        let mut remap = BTreeMap::new();
        let mut reverse = BTreeMap::new();
        let mut bind = |old: usize, new: usize, dec: &mut Decider| -> Result<bool, Refusal> {
            dec.charge_work(1)?;
            let (Some(before), Some(after)) = (
                self.mesh.positions.get(old),
                boundaries.mesh.positions.get(new),
            ) else {
                return Ok(false);
            };
            if before
                .iter()
                .zip(after)
                .any(|(a, b)| a.to_bits() != b.to_bits())
                || remap.get(&old).is_some_and(|&value| value != new)
                || reverse.get(&new).is_some_and(|&value| value != old)
            {
                return Ok(false);
            }
            remap.insert(old, new);
            reverse.insert(new, old);
            Ok(true)
        };
        for (&old_loop, &new_loop) in before.loops.iter().zip(&after.loops) {
            dec.charge_work(1)?;
            match (self.model.loop_(old_loop), model.loop_(new_loop)) {
                (Loop::Vertex(a), Loop::Vertex(b)) => {
                    if !same_vertex(self.model, Some(*a), model, Some(*b)) {
                        return Ok(None);
                    }
                    // Singleton trims bound the chart. Their poles are emitted
                    // as face-local samples, not the unused topological point
                    // reserved by edges::sample. Keep that generation order.
                }
                (Loop::Coedges(a), Loop::Coedges(b)) if a.len() == b.len() => {
                    for (&a, &b) in a.iter().zip(b) {
                        dec.charge_work(1)?;
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
                            return Ok(None);
                        }
                        let Some(old) = self.edges.get(&a.edge) else {
                            return Ok(None);
                        };
                        let new = &boundaries.mesh.edges[boundaries.edges[&b.edge]];
                        if old.vertices.len() != new.vertices.len()
                            || old.max_chord_error != new.max_chord_error
                        {
                            return Ok(None);
                        }
                        for (&old, &new) in old.vertices.iter().zip(&new.vertices) {
                            if !bind(old, new, dec)? {
                                return Ok(None);
                            }
                        }
                    }
                }
                _ => return Ok(None),
            }
        }
        Ok(Some(remap))
    }

    fn prior_matches(
        &self,
        face_index: usize,
        remap: &BTreeMap<usize, usize>,
        prior: &Chords,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut old_prior = Chords::new();
        for (&(a, b), &first_face) in &self.chords {
            dec.charge_work(1)?;
            if first_face < face_index
                && let (Some(&a), Some(&b)) = (remap.get(&a), remap.get(&b))
            {
                old_prior.insert((a.min(b), a.max(b)));
            }
        }
        let mut vertices = BTreeSet::new();
        for &vertex in remap.values() {
            dec.charge_work(1)?;
            vertices.insert(vertex);
        }
        for &(a, b) in prior {
            dec.charge_work(1)?;
            if vertices.contains(&a) && vertices.contains(&b) && !old_prior.remove(&(a, b)) {
                return Ok(false);
            }
        }
        Ok(old_prior.is_empty())
    }

    pub(crate) fn append(
        &self,
        model: &Model,
        face: FaceId,
        boundaries: &mut Boundaries,
        prior: &Chords,
        options: &MeshOptions,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for (index, old) in self.mesh.faces.iter().enumerate() {
            dec.charge_work(1)?;
            let Some(mut remap) = self.boundary_map(old, model, face, boundaries, dec)? else {
                continue;
            };
            if !self.prior_matches(index, &remap, prior, dec)? {
                continue;
            }
            let triangles = &self.mesh.triangles[old.triangles.clone()];
            let mut interior = BTreeSet::new();
            for triangle in triangles {
                dec.charge_work(3)?;
                for vertex in triangle.vertices {
                    if !remap.contains_key(&vertex) {
                        interior.insert(vertex);
                    }
                }
            }
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
                dec.charge_work(1)?;
                let next = boundaries.mesh.positions.len();
                boundaries.mesh.positions.push(self.mesh.positions[vertex]);
                remap.insert(vertex, next);
            }
            let start = boundaries.mesh.triangles.len();
            dec.charge_work((triangles.len() as u64).saturating_mul(3))?;
            boundaries.mesh.triangles.reserve(triangles.len());
            for triangle in triangles {
                dec.check_deadline()?;
                boundaries.mesh.triangles.push(Triangle {
                    vertices: triangle.vertices.map(|vertex| remap[&vertex]),
                    face,
                });
            }
            dec.charge_work(1)?;
            let mut report = old.clone();
            report.face = face;
            report.triangles = start..boundaries.mesh.triangles.len();
            boundaries.mesh.faces.push(report);
            for control in &self.mesh.unavailable_controls {
                dec.charge_work(1)?;
                if control.face == old.face {
                    let mut control = control.clone();
                    control.face = face;
                    boundaries.mesh.unavailable_controls.push(control);
                }
            }
            return Ok(true);
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expired_deadline_refuses_prior_mesh_indexing() {
        let model = Model::new();
        let mesh = Mesh::default();
        let mut dec = Decider::default();
        dec.set_deadline(Some(
            std::time::Instant::now() - std::time::Duration::from_secs(1),
        ));
        let refusal = Previous::new(&model, &mesh, &mut dec).err().unwrap();
        assert_eq!(refusal.code, "numbers.deadline");
    }
}
