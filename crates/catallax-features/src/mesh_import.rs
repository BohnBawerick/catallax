//! The B-rep of an imported closed mesh, RECIPE.md section 4.19.
//!
//! [`ClosedMesh::read`] has already proved the surface closed, consistently
//! oriented and free of crossings, exactly and once for the whole mesh. The
//! construction here therefore proves nothing face by face. It merges each
//! edge-connected set of coplanar triangles into one planar face, records the
//! exact incidences it creates, and costs time linear in the triangle count.

use std::collections::{BTreeMap, HashMap};

use catallax_import::mesh::ClosedMesh;
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{Origin, Role, StepId};
use catallax_shapes::{Curve, Isometry, Line, Plane, Point3, PointDef, RPoint3, Surface};
use catallax_topology::{
    Body, BodyId, BodyKind, EdgeId, Lump, Model, Sense, Shell, Vertex, VertexId,
};

struct Builder<'a> {
    mesh: &'a ClosedMesh,
    points: Vec<RPoint3>,
    /// The same points as shared exact reals: every use clones the handle.
    reals: Vec<Point3>,
    vertices: Vec<VertexId>,
    edges: HashMap<(u32, u32), EdgeId>,
    origins: Vec<Origin>,
}

impl Builder<'_> {
    fn edge(&mut self, model: &mut Model, a: u32, b: u32) -> (EdgeId, Sense) {
        let key = (a.min(b), a.max(b));
        let (points, reals, vertices) = (&self.points, &self.reals, &self.vertices);
        let edge = *self.edges.entry(key).or_insert_with(|| {
            let (start, end) = (&points[key.0 as usize], &points[key.1 as usize]);
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                origin: reals[key.0 as usize].clone(),
                dir: end.minus(start).to_real(),
            }));
            model.add_edge(catallax_topology::Edge {
                curve,
                sense: Sense::Same,
                start: Some(vertices[key.0 as usize]),
                end: Some(vertices[key.1 as usize]),
            })
        });
        (edge, if a < b { Sense::Same } else { Sense::Reversed })
    }

    /// One planar face through triangle `t`'s plane, bounded by vertex loops.
    fn face(&mut self, model: &mut Model, t: u32, loops: &[Vec<u32>]) -> catallax_topology::FaceId {
        let corners = self.mesh.triangles()[t as usize];
        let [a, b, c] = corners.map(|v| &self.points[v as usize]);
        let (ab, ac) = (b.minus(a), c.minus(a));
        let plane = Plane {
            origin: self.reals[corners[0] as usize].clone(),
            normal: ab.cross(&ac).to_real(),
            ref_dir: ab.to_real(),
        };
        let paths = loops
            .iter()
            .map(|lp| {
                (0..lp.len())
                    .map(|i| self.edge(model, lp[i], lp[(i + 1) % lp.len()]))
                    .collect()
            })
            .collect();
        crate::extrude::face(
            model,
            Surface::Plane(plane),
            Sense::Same,
            paths,
            None,
            self.origins.clone(),
        )
    }
}

fn root(parents: &mut [u32], mut t: u32) -> u32 {
    while parents[t as usize] != t {
        let up = parents[parents[t as usize] as usize];
        parents[t as usize] = up;
        t = up;
    }
    t
}

/// Boundary loops of one coplanar region, outer loop first, or `None` when a
/// boundary vertex is passed twice and the region is built triangle by triangle.
fn region_loops(
    mesh: &ClosedMesh,
    region: &[u32],
    of: &[u32],
    dec: &mut Decider,
) -> Result<Option<Vec<Vec<u32>>>, Refusal> {
    let mut next: BTreeMap<u32, u32> = BTreeMap::new();
    let mut starts = Vec::new();
    for &t in region {
        let tri = mesh.triangles()[t as usize];
        for slot in 0..3 {
            if of[mesh.twin(t, slot) as usize] != of[t as usize] {
                let (a, b) = (tri[slot], tri[(slot + 1) % 3]);
                if next.insert(a, b).is_some() {
                    return Ok(None);
                }
                starts.push(a);
            }
        }
    }
    let mut loops = Vec::new();
    for start in starts {
        if !next.contains_key(&start) {
            continue;
        }
        let mut lp = vec![start];
        let mut at = next.remove(&start).expect("present");
        while at != start {
            lp.push(at);
            let Some(following) = next.remove(&at) else {
                return Ok(None);
            };
            at = following;
        }
        loops.push(lp);
    }
    let (k, sign) = mesh.projection(region[0]);
    let mut outer = Vec::new();
    for (i, lp) in loops.iter().enumerate() {
        if mesh.loop_area_sign(k, lp, dec)? == sign {
            outer.push(i);
        }
    }
    if outer.len() != 1 {
        return Ok(None);
    }
    loops.swap(0, outer[0]);
    Ok(Some(loops))
}

/// Builds the solid. Returns the body and the sum over axes of the largest
/// absolute coordinate, which bounds the body for later operations.
pub(crate) fn build(
    model: &mut Model,
    mesh: &ClosedMesh,
    placement: &Isometry,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Rational), Refusal> {
    let triangles = mesh.triangles();
    if triangles.is_empty() {
        return Ok((catallax_topology::empty(model), Rational::zero()));
    }
    let points = (0..mesh.vertex_count() as u32)
        .map(|v| {
            // Scaling and placing one vertex is about thirty exact operations
            // on numbers the size of its coordinates.
            let size = mesh
                .file_vertex(v)
                .iter()
                .map(|x| (x.numerator().bit_length() + x.denominator().bit_length()) / 64 + 1)
                .max()
                .unwrap_or(1);
            dec.charge_work(30u64.saturating_mul(size.saturating_mul(size)))?;
            let [x, y, z] = mesh.vertex(v);
            placement
                .apply_rational_point(&RPoint3::new(x, y, z))
                .ok_or_else(|| {
                    Refusal::new(
                        "import.mesh-pose",
                        "a mesh pose must be an exact rational motion",
                    )
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let radius = (0..3).fold(Rational::zero(), |sum, axis| {
        let largest = points
            .iter()
            .map(|p| [&p.x, &p.y, &p.z][axis].abs())
            .max()
            .expect("a nonempty mesh has vertices");
        sum + largest
    });
    let reals: Vec<Point3> = points.iter().map(RPoint3::to_real).collect();
    let vertices = reals
        .iter()
        .map(|p| {
            let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
            model.add_vertex(Vertex { point })
        })
        .collect();
    let mut builder = Builder {
        mesh,
        points,
        reals,
        vertices,
        edges: HashMap::new(),
        origins: vec![Origin {
            step: step.clone(),
            role: Role::Mesh,
        }],
    };

    let mut parents: Vec<u32> = (0..triangles.len() as u32).collect();
    for t in 0..triangles.len() as u32 {
        for slot in 0..3 {
            let u = mesh.twin(t, slot);
            if t < u && mesh.coplanar_across(t, slot, dec)? {
                let (a, b) = (root(&mut parents, t), root(&mut parents, u));
                parents[a.max(b) as usize] = a.min(b);
            }
        }
    }
    let of: Vec<u32> = (0..triangles.len() as u32)
        .map(|t| root(&mut parents, t))
        .collect();
    let mut regions: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (t, &r) in of.iter().enumerate() {
        regions.entry(r).or_default().push(t as u32);
    }

    let mut shell_of = vec![0usize; triangles.len()];
    for (i, shell) in mesh.shells().iter().enumerate() {
        for &t in &shell.triangles {
            shell_of[t as usize] = i;
        }
    }
    let mut shell_faces = vec![Vec::new(); mesh.shells().len()];
    for region in regions.values() {
        // Face construction is linear in the region's triangles.
        dec.charge_work(region.len() as u64)?;
        let faces = &mut shell_faces[shell_of[region[0] as usize]];
        match region_loops(mesh, region, &of, dec)? {
            Some(loops) => faces.push(builder.face(model, region[0], &loops)),
            None => {
                for &t in region {
                    let loops = [triangles[t as usize].to_vec()];
                    faces.push(builder.face(model, t, &loops));
                }
            }
        }
    }
    let shells: Vec<_> = shell_faces
        .into_iter()
        .map(|faces| model.add_shell(Shell { faces }))
        .collect();
    let lumps = mesh
        .lumps()
        .into_iter()
        .map(|(outer, cavities)| {
            let shells = std::iter::once(outer)
                .chain(cavities)
                .map(|s| shells[s])
                .collect();
            model.add_lump(Lump { shells })
        })
        .collect();
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
        marks: Default::default(),
    });
    Ok((body, radius))
}
