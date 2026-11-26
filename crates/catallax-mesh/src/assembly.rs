//! One face's mesh as its own result, and the ordered body assembly.
//!
//! Edge sampling fixes every shared boundary vertex before any face is built.
//! A face builder then writes into a [`FaceDraft`]. The draft reads the shared
//! samples and appends only the face's own vertices, triangles and control
//! reports. It finishes as a [`FaceResult`]. Fresh and reused faces both reach
//! the body mesh through [`Assembly::push`], which checks the result's vertex
//! references and the chords it relied on before placing it in face order.

use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Index,
};

use catallax_numbers::{Ball, MarkLog, Refusal};
use catallax_topology::{EdgeId, FaceId, Loop, Model, VertexId};

use crate::{
    EdgeMesh, FaceMesh, Mesh, Triangle, UnavailableControl,
    edges::{Boundaries, Position, Rails},
    failure,
};

/// Pairs of shared boundary vertices, smaller index first.
pub(crate) type Chords = BTreeSet<(usize, usize)>;

/// The positions one face builder sees: the shared boundary samples, then the
/// face's own vertices numbered from `base`. Indices in between belong to
/// other faces and are never readable here.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Positions<'a> {
    shared: &'a [[f64; 3]],
    base: usize,
    local: Vec<[f64; 3]>,
}

impl<'a> Positions<'a> {
    fn new(shared: &'a [[f64; 3]], base: usize) -> Self {
        debug_assert!(base >= shared.len());
        Self {
            shared,
            base,
            local: Vec::new(),
        }
    }

    /// Positions with no shared prefix, numbered from zero.
    pub(crate) fn detached(local: Vec<[f64; 3]>) -> Positions<'static> {
        Positions {
            shared: &[],
            base: 0,
            local,
        }
    }

    /// The face's own vertices, in generation order.
    pub(crate) fn into_local(self) -> Vec<[f64; 3]> {
        self.local
    }

    /// One past the largest readable index.
    pub(crate) fn len(&self) -> usize {
        self.base + self.local.len()
    }

    pub(crate) fn push(&mut self, position: [f64; 3]) {
        self.local.push(position);
    }

    /// Drops the face's own vertices from index `len` on.
    pub(crate) fn truncate(&mut self, len: usize) {
        self.local.truncate(len.saturating_sub(self.base));
    }
}

impl Index<usize> for Positions<'_> {
    type Output = [f64; 3];

    fn index(&self, index: usize) -> &[f64; 3] {
        if index < self.shared.len() {
            return &self.shared[index];
        }
        let local = index
            .checked_sub(self.base)
            .expect("a face reads only shared samples and its own vertices");
        &self.local[local]
    }
}

/// What one face builder appends, with the shared samples it reads.
#[derive(Debug, PartialEq)]
pub(crate) struct DraftMesh<'a> {
    pub(crate) positions: Positions<'a>,
    /// Shared B-rep boundary polylines, read only.
    pub(crate) edges: &'a [EdgeMesh],
    pub(crate) triangles: Vec<Triangle>,
    pub(crate) faces: Vec<FaceMesh>,
    pub(crate) unavailable_controls: Vec<UnavailableControl>,
}

/// A face builder's workspace. Everything shared is read only, except a cache
/// of exact boundary-vertex enclosures that planar faces fill on demand.
pub(crate) struct FaceDraft<'a> {
    pub(crate) mesh: DraftMesh<'a>,
    pub(crate) edges: &'a BTreeMap<EdgeId, usize>,
    pub(crate) vertices: &'a BTreeMap<VertexId, usize>,
    pub(crate) exact_positions: &'a BTreeMap<usize, Position>,
    pub(crate) enclosures: &'a mut BTreeMap<usize, [Ball; 3]>,
    pub(crate) rails: &'a Rails,
}

impl<'a> FaceDraft<'a> {
    /// A draft whose own vertices are numbered from `base`.
    pub(crate) fn new(
        shared: &'a Boundaries,
        base: usize,
        enclosures: &'a mut BTreeMap<usize, [Ball; 3]>,
    ) -> Self {
        Self {
            mesh: DraftMesh {
                positions: Positions::new(&shared.mesh.positions, base),
                edges: &shared.mesh.edges,
                triangles: Vec::new(),
                faces: Vec::new(),
                unavailable_controls: Vec::new(),
            },
            edges: &shared.edges,
            vertices: &shared.vertices,
            exact_positions: &shared.exact_positions,
            enclosures,
            rails: &shared.rails,
        }
    }

    /// The finished face, recording the earlier chords it was built against.
    pub(crate) fn finish(
        self,
        model: &Model,
        face: FaceId,
        prior: &Chords,
    ) -> Result<FaceResult, Refusal> {
        let DraftMesh {
            positions,
            triangles,
            faces,
            unavailable_controls,
            ..
        } = self.mesh;
        let Ok([report]) = <[FaceMesh; 1]>::try_from(faces) else {
            return Err(failure("a face builder must report its face exactly once"));
        };
        let mut shared = BTreeSet::new();
        for &boundary in &model.face(face).loops {
            if let Loop::Coedges(coedges) = model.loop_(boundary) {
                for &coedge in coedges {
                    if let Some(&edge) = self.edges.get(&model.coedge(coedge).edge) {
                        shared.extend(self.mesh.edges[edge].vertices.iter().copied());
                    }
                }
            }
        }
        let relies_on = chords_among(prior, &shared);
        Ok(FaceResult {
            face,
            base: positions.base,
            positions: positions.local,
            triangles,
            report,
            unavailable_controls,
            shared,
            relies_on,
        })
    }
}

/// One face's complete mesh, independent of the faces built before it.
///
/// Triangle indices below the boundary sample count refer to shared edge
/// samples. Indices from `base` refer to `positions`. The report's triangle
/// range is relative to `triangles` until assembly.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct FaceResult {
    pub(crate) face: FaceId,
    pub(crate) base: usize,
    pub(crate) positions: Vec<[f64; 3]>,
    pub(crate) triangles: Vec<Triangle>,
    /// Chord, aspect and angular evidence for the face.
    pub(crate) report: FaceMesh,
    pub(crate) unavailable_controls: Vec<UnavailableControl>,
    /// Samples of the shared edges that bound the face.
    pub(crate) shared: BTreeSet<usize>,
    /// Chords between `shared` samples that earlier faces had emitted when
    /// this face was built. The face's triangulation avoided each of them.
    pub(crate) relies_on: Chords,
}

impl FaceResult {
    /// Chords between shared samples that this face uses on two triangles.
    fn chords(&self, boundary: usize) -> Chords {
        let mut once = Chords::new();
        let mut twice = Chords::new();
        for triangle in &self.triangles {
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                if a < boundary && b < boundary {
                    let pair = (a.min(b), a.max(b));
                    if !once.insert(pair) {
                        twice.insert(pair);
                    }
                }
            }
        }
        twice
    }
}

fn chords_among(chords: &Chords, vertices: &BTreeSet<usize>) -> Chords {
    let mut among = Chords::new();
    for &a in vertices {
        for &(a, b) in chords.range((a, 0)..=(a, usize::MAX)) {
            if vertices.contains(&b) {
                among.insert((a, b));
            }
        }
    }
    among
}

/// The body mesh after its shared edge samples, built face by face in order.
pub(crate) struct Assembly {
    boundary: usize,
    positions: Vec<[f64; 3]>,
    triangles: Vec<Triangle>,
    faces: Vec<FaceMesh>,
    unavailable_controls: Vec<UnavailableControl>,
    prior: Chords,
}

impl Assembly {
    pub(crate) fn new(boundaries: &Boundaries) -> Self {
        Self {
            boundary: boundaries.mesh.positions.len(),
            positions: Vec::new(),
            triangles: Vec::new(),
            faces: Vec::new(),
            unavailable_controls: Vec::new(),
            prior: Chords::new(),
        }
    }

    /// The index the next face's first own vertex receives.
    pub(crate) fn base(&self) -> usize {
        self.boundary + self.positions.len()
    }

    /// Chords between shared samples emitted by the faces placed so far.
    pub(crate) fn prior(&self) -> &Chords {
        &self.prior
    }

    /// Places a face after the earlier ones. Returns false, placing nothing,
    /// when the chords it relied on differ from those emitted so far; the
    /// caller then builds the face again in order. A triangle index outside
    /// the shared samples and the face's own vertices is refused.
    pub(crate) fn push(&mut self, result: FaceResult) -> Result<bool, Refusal> {
        if chords_among(&self.prior, &result.shared) != result.relies_on {
            return Ok(false);
        }
        let base = self.base();
        let local = result.base..result.base + result.positions.len();
        let start = self.triangles.len();
        self.triangles.reserve(result.triangles.len());
        for triangle in &result.triangles {
            let mut vertices = triangle.vertices;
            for vertex in &mut vertices {
                if *vertex >= self.boundary {
                    if !local.contains(vertex) {
                        self.triangles.truncate(start);
                        return Err(failure(
                            "a face mesh references a vertex outside its shared edges and its own samples",
                        ));
                    }
                    *vertex = *vertex - result.base + base;
                }
            }
            self.triangles.push(Triangle {
                vertices,
                face: triangle.face,
            });
        }
        self.prior.extend(result.chords(self.boundary));
        self.positions.extend(result.positions);
        let mut report = result.report;
        report.triangles = start + report.triangles.start..start + report.triangles.end;
        self.faces.push(report);
        self.unavailable_controls
            .extend(result.unavailable_controls);
        Ok(true)
    }

    /// The assembled mesh, after the shared samples and edge polylines.
    pub(crate) fn finish(self, boundaries: Boundaries) -> (Mesh, MarkLog) {
        let Boundaries {
            mut mesh, marks, ..
        } = boundaries;
        debug_assert_eq!(mesh.positions.len(), self.boundary);
        mesh.positions.extend(self.positions);
        mesh.triangles.extend(self.triangles);
        mesh.faces.extend(self.faces);
        mesh.unavailable_controls.extend(self.unavailable_controls);
        (mesh, marks)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Rational;
    use catallax_shapes::Pose;

    fn face(face: FaceId, base: usize, own: &[[f64; 3]], triangles: &[[usize; 3]]) -> FaceResult {
        FaceResult {
            face,
            base,
            positions: own.to_vec(),
            triangles: triangles
                .iter()
                .map(|&vertices| Triangle { vertices, face })
                .collect(),
            report: FaceMesh {
                face,
                triangles: 0..triangles.len(),
                max_chord_error: Rational::zero(),
                max_aspect_ratio: None,
                max_angular_deviation: None,
            },
            unavailable_controls: Vec::new(),
            shared: [0, 1, 2, 3].into(),
            relies_on: Chords::new(),
        }
    }

    #[test]
    fn assembly_checks_relied_chords_and_relocates_a_face_built_earlier() {
        let mut model = Model::new();
        let body = catallax_topology::block(
            &mut model,
            [1.into(), 1.into(), 1.into()],
            &Pose::identity(),
        )
        .unwrap();
        let faces = model.faces_of_body(body);
        let shared = Boundaries {
            mesh: Mesh {
                positions: vec![[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                ..Mesh::default()
            },
            edges: BTreeMap::new(),
            vertices: BTreeMap::new(),
            exact_positions: BTreeMap::new(),
            rails: Rails::default(),
            marks: MarkLog::new(),
        };
        let mut assembly = Assembly::new(&shared);
        // Two faces drafted against the same empty prefix, as workers would.
        let first = face(faces[0], 4, &[[5.0; 3]], &[[0, 1, 2], [0, 2, 3]]);
        let second = face(faces[1], 4, &[[6.0; 3]], &[[1, 3, 4]]);
        assert!(assembly.push(first).unwrap());
        assert_eq!(assembly.prior(), &Chords::from([(0, 2)]));
        // The second face was drafted before the first emitted chord 0-2.
        assert!(!assembly.push(second.clone()).unwrap());
        let mut rebuilt = second;
        rebuilt.relies_on = Chords::from([(0, 2)]);
        assert!(assembly.push(rebuilt).unwrap());
        // A vertex of the first face is neither shared nor the third face's own.
        let mut stray = face(faces[2], 6, &[], &[[0, 1, 4]]);
        stray.relies_on = Chords::from([(0, 2)]);
        assert_eq!(
            assembly.push(stray).unwrap_err().code,
            "mesh.precision-or-budget"
        );
        let (mesh, _) = assembly.finish(shared);
        assert_eq!(mesh.positions.len(), 6);
        assert_eq!(mesh.positions[4..], [[5.0; 3], [6.0; 3]]);
        assert_eq!(mesh.triangles[2].vertices, [1, 3, 5]);
        assert_eq!(
            mesh.faces
                .iter()
                .map(|f| f.triangles.clone())
                .collect::<Vec<_>>(),
            [0..2, 2..3]
        );
    }
}
