use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use catallax_numbers::{Decider, Refusal};
use catallax_shapes::{Curve, Point3, Surface};
use catallax_topology::{EdgeId, FaceId, Loop, Model, VertexId};

use crate::{
    EdgeMesh, FaceMesh, Mesh, MeshOptions, Triangle,
    assembly::{Assembly, Chords, FaceResult},
    edges::Boundaries,
};

/// Certified face provenance from one [`crate::mesh_body_cached`] result.
///
/// Positions live in the public mesh. Expanded bodies retain the physical cell
/// reports, triangles and edges that were merged in the public picking view.
#[derive(Clone, Debug)]
pub struct FaceCache {
    options: MeshOptions,
    data: Arc<Capture>,
}

impl FaceCache {
    pub(crate) fn new(options: MeshOptions, data: Capture) -> Self {
        Self {
            options,
            data: Arc::new(data),
        }
    }

    /// Number of physical faces served from the previous cache.
    pub fn reused_faces(&self) -> usize {
        self.data.reused
    }

    /// Number of physical faces triangulated for this request.
    pub fn meshed_faces(&self) -> usize {
        self.data.records.len() - self.data.reused
    }

    pub(crate) fn input<'a>(
        &'a self,
        model: &'a Model,
        mesh: &'a Mesh,
        options: &MeshOptions,
    ) -> Option<Input<'a>> {
        (&self.options == options).then_some(Input {
            model,
            geometry: mesh,
            source: self.data.expanded.as_ref().unwrap_or(mesh),
            records: Some(&self.data.records),
        })
    }
}

#[derive(Debug, Default)]
pub(crate) struct Capture {
    records: BTreeMap<FaceId, Record>,
    pub(crate) expanded: Option<Mesh>,
    reused: usize,
}

impl Capture {
    pub(crate) fn face(&mut self, result: &FaceResult, options: &MeshOptions, reused: bool) {
        self.records.insert(
            result.face,
            Record {
                chord: options.chord_tolerance.clone(),
                peak: result.peak,
                uncompacted: result.uncompacted,
            },
        );
        self.reused += usize::from(reused);
    }
}

#[derive(Debug)]
pub(crate) struct Record {
    chord: catallax_numbers::Rational,
    peak: usize,
    uncompacted: usize,
}

#[derive(Clone, Copy)]
pub(crate) struct Input<'a> {
    pub(crate) model: &'a Model,
    pub(crate) geometry: &'a Mesh,
    pub(crate) source: &'a Mesh,
    pub(crate) records: Option<&'a BTreeMap<FaceId, Record>>,
}

pub(crate) struct Edit<'a, 'out> {
    pub(crate) previous: Option<Input<'a>>,
    pub(crate) capture: Option<&'out mut Capture>,
}

pub(crate) struct Previous<'a> {
    model: &'a Model,
    mesh: &'a Mesh,
    positions: &'a [[f64; 3]],
    records: Option<&'a BTreeMap<FaceId, Record>>,
    edges: BTreeMap<EdgeId, &'a EdgeMesh>,
    chords: BTreeMap<(usize, usize), usize>,
}

pub(crate) struct ReusePlan<'a> {
    face: FaceId,
    base: usize,
    peak: usize,
    uncompacted: usize,
    report: &'a FaceMesh,
    positions: &'a [[f64; 3]],
    triangles: &'a [Triangle],
    controls: &'a [crate::UnavailableControl],
    remap: BTreeMap<usize, usize>,
    interior: BTreeSet<usize>,
    shared: BTreeSet<usize>,
    relies_on: Chords,
}

impl ReusePlan<'_> {
    /// Copies a matching face into the new body without evaluating its surface.
    pub(crate) fn emit(self, dec: &mut Decider) -> Result<FaceResult, Refusal> {
        let Self {
            face,
            base,
            peak,
            uncompacted,
            report: old,
            positions: source_positions,
            triangles,
            controls,
            mut remap,
            interior,
            shared,
            relies_on,
        } = self;
        // Preserve the original position order, including periodic chart
        // copies that survived compaction. Only boundary indices are shared.
        let mut positions = Vec::with_capacity(interior.len());
        for vertex in interior {
            dec.charge_work(1)?;
            remap.insert(vertex, base + positions.len());
            positions.push(source_positions[vertex]);
        }
        dec.charge_work((triangles.len() as u64).saturating_mul(3))?;
        let mut copied = Vec::with_capacity(triangles.len());
        for triangle in triangles {
            dec.check_deadline()?;
            copied.push(Triangle {
                vertices: triangle.vertices.map(|vertex| remap[&vertex]),
                face,
            });
        }
        dec.charge_work(1)?;
        let mut report = old.clone();
        report.face = face;
        report.triangles = 0..copied.len();
        let mut unavailable_controls = Vec::new();
        for control in controls {
            dec.charge_work(1)?;
            if control.face == old.face {
                let mut control = control.clone();
                control.face = face;
                unavailable_controls.push(control);
            }
        }
        Ok(FaceResult {
            face,
            base,
            peak,
            uncompacted,
            positions,
            triangles: copied,
            report,
            unavailable_controls,
            shared,
            relies_on,
        })
    }
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
    pub(crate) fn new(input: Input<'a>, dec: &mut Decider) -> Result<Self, Refusal> {
        let Input {
            model,
            geometry,
            source: mesh,
            records,
        } = input;
        let positions = &geometry.positions;
        // Reserve both vertex-sized initializations before allocating. Each
        // remaining record visit is charged before its scan or copy.
        dec.charge_work((positions.len() as u64).saturating_mul(2))?;
        let mut boundary = vec![false; positions.len()];
        for edge in &mesh.edges {
            dec.charge_work(1)?;
            for &vertex in &edge.vertices {
                dec.charge_work(1)?;
                boundary[vertex] = true;
            }
        }
        dec.check_deadline()?;
        let mut owner = vec![usize::MAX; positions.len()];
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
            positions,
            records,
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
            let (Some(before), Some(after)) =
                (self.positions.get(old), boundaries.mesh.positions.get(new))
            else {
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
        // Sorted triangle and canonical cyclic order survive only when the
        // relative order of the shared indices is unchanged.
        if !remap.values().is_sorted() {
            return Ok(None);
        }
        Ok(Some(remap))
    }

    /// Chords that earlier old faces emitted between this old face's samples,
    /// in new indices: what the cached triangulation was built against.
    fn relied_on(
        &self,
        face_index: usize,
        remap: &BTreeMap<usize, usize>,
        dec: &mut Decider,
    ) -> Result<Chords, Refusal> {
        let mut relied_on = Chords::new();
        for (&old, &new) in remap {
            for (&(_, b), &first_face) in self.chords.range((old, 0)..=(old, usize::MAX)) {
                dec.charge_work(1)?;
                if first_face < face_index
                    && let Some(&b) = remap.get(&b)
                {
                    relied_on.insert((new, b));
                }
            }
        }
        Ok(relied_on)
    }

    /// Finds a face with matching definition, shared samples and earlier chords
    /// without copying its point or triangle buffers.
    pub(crate) fn plan_face(
        &self,
        model: &Model,
        face: FaceId,
        boundaries: &Boundaries,
        assembly: &Assembly,
        options: &MeshOptions,
        dec: &mut Decider,
    ) -> Result<Option<ReusePlan<'_>>, Refusal> {
        for (index, old) in self.mesh.faces.iter().enumerate() {
            dec.charge_work(1)?;
            let record = if let Some(records) = self.records {
                let Some(record) = records.get(&old.face) else {
                    continue;
                };
                if record.chord != options.chord_tolerance {
                    continue;
                }
                Some(record)
            } else {
                None
            };
            let Some(remap) = self.boundary_map(old, model, face, boundaries, dec)? else {
                continue;
            };
            let relies_on = self.relied_on(index, &remap, dec)?;
            let mut shared = BTreeSet::new();
            for &vertex in remap.values() {
                dec.charge_work(1)?;
                shared.insert(vertex);
            }
            if !prior_matches(&relies_on, &shared, assembly.prior(), dec)? {
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
            let base = assembly.base();
            let peak = record.map_or(interior.len(), |record| record.peak);
            let uncompacted = record.map_or(interior.len(), |record| record.uncompacted);
            if uncompacted < interior.len() || peak > options.max_vertices.saturating_sub(base) {
                continue;
            }
            return Ok(Some(ReusePlan {
                face,
                base,
                peak,
                uncompacted,
                report: old,
                positions: self.positions,
                triangles,
                controls: &self.mesh.unavailable_controls,
                remap,
                interior,
                shared,
                relies_on,
            }));
        }
        Ok(None)
    }
}

/// Whether the chords emitted so far between `shared` samples are exactly
/// `relies_on`. Charged per chord, as the old face loop was.
fn prior_matches(
    relies_on: &Chords,
    shared: &BTreeSet<usize>,
    prior: &Chords,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut count = 0;
    for &a in shared {
        for &pair in prior.range((a, 0)..=(a, usize::MAX)) {
            dec.charge_work(1)?;
            if shared.contains(&pair.1) {
                if !relies_on.contains(&pair) {
                    return Ok(false);
                }
                count += 1;
            }
        }
    }
    Ok(count == relies_on.len())
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
        let input = Input {
            model: &model,
            geometry: &mesh,
            source: &mesh,
            records: None,
        };
        let refusal = Previous::new(input, &mut dec).err().unwrap();
        assert_eq!(refusal.code, "numbers.deadline");
    }
}
