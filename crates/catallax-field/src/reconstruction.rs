//! Opt-in, bounded M23 recognition of an axis box with one coaxial through-bore.
//! The reconstructed definition is exact; fidelity to the input is NOT certified.
mod deviation;
use crate::Plan;
use catallax_import::mesh::ClosedMesh;
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{BooleanKind, Number, Operation, Pose, Recipe, Step, StepId};
pub use deviation::{SampledDeviation, Sampling, filled_envelope_control, sampled_deviation};
use std::collections::{BTreeMap, BTreeSet};

type P = [Rational; 3];
/// Maximum source size for this new recognizer; does not change M22's limit.
pub const MAX_TRIANGLES: usize = 4096;
/// Reconstructed solid and its separate imported reference identity.
pub struct Reconstruction {
    /// Original file bytes' SHA-256; never the identity of the fitted solid.
    source_sha256: String,
    fidelity: SampledDeviation,
    /// Millimetres per source unit.
    source_scale: Rational,
    /// Inferred ordinary definition, retained independently of the source mesh.
    recipe: Recipe,
    /// Evaluated analytic B-rep, not the original faceted import.
    brep: catallax_features::Built,
    /// Usable field of the inferred mathematical solid, not the source scan.
    field: Plan,
    /// Seven connected regions: six planes and one fitted bore.
    region_triangles: [usize; 7],
    /// Inferred circle centres, common to both z levels, in mm.
    centre: [Rational; 2],
    /// Bottom/top radii in mm; equal means cylinder, unequal means cone.
    radii: [Rational; 2],
    /// Inferred outer box in mm.
    bounds: [[Rational; 2]; 3],
}
impl Reconstruction {
    /// Mandatory unproved bidirectional fidelity report at fixed 4/8/128 sampling.
    pub fn fidelity(&self) -> &SampledDeviation {
        &self.fidelity
    }
    /// Hash of the original imported bytes, distinct from the inferred recipe.
    pub fn source_sha256(&self) -> &str {
        &self.source_sha256
    }
    /// Exact millimetres per source unit used for this fit.
    pub fn source_scale(&self) -> &Rational {
        &self.source_scale
    }
    /// The immutable inferred exact definition.
    pub fn recipe(&self) -> &Recipe {
        &self.recipe
    }
    /// The validated reconstructed B-rep.
    pub fn brep(&self) -> &catallax_features::Built {
        &self.brep
    }
    /// Field of the inferred solid, not a certificate about the scan.
    pub fn field(&self) -> &Plan {
        &self.field
    }
    /// Complete per-region triangle counts (six planes, then bore).
    pub fn region_triangles(&self) -> &[usize; 7] {
        &self.region_triangles
    }
    /// Inferred common circle centre in mm.
    pub fn centre(&self) -> &[Rational; 2] {
        &self.centre
    }
    /// Inferred lower/upper radii in mm.
    pub fn radii(&self) -> &[Rational; 2] {
        &self.radii
    }
    /// Inferred box extents in mm.
    pub fn bounds(&self) -> &[[Rational; 2]; 3] {
        &self.bounds
    }
}

fn fail(message: &str) -> Refusal {
    Refusal::new("reconstruction.unsupported", message)
}
fn id(s: &str) -> StepId {
    StepId::new(s).expect("static id")
}
fn cross(a: &[Rational; 2], b: &[Rational; 2]) -> Rational {
    &a[0] * &b[1] - &a[1] * &b[0]
}
fn circle(points: &[P]) -> Result<([Rational; 2], Rational), Refusal> {
    if points.len() < 3 {
        return Err(fail("a bore ring needs three noncollinear points"));
    }
    let a = &points[0];
    let b = &points[1];
    let c = &points[2];
    let u = [&b[0] - &a[0], &b[1] - &a[1]];
    let v = [&c[0] - &a[0], &c[1] - &a[1]];
    let square = |p: &P| &p[0] * &p[0] + &p[1] * &p[1];
    let ub = (square(b) - square(a)).mul_pow2(-1);
    let vc = (square(c) - square(a)).mul_pow2(-1);
    let det = cross(&u, &v);
    let centre = [
        (&ub * &v[1] - &u[1] * &vc).checked_div(&det),
        (&u[0] * &vc - &ub * &v[0]).checked_div(&det),
    ];
    let [Some(x), Some(y)] = centre else {
        return Err(fail("ambiguous circle fit"));
    };
    let centre = [x, y];
    let radius_sq = |p: &P| {
        (0..2).fold(Rational::zero(), |s, i| {
            let d = &p[i] - &centre[i];
            s + &d * &d
        })
    };
    let r2 = radius_sq(a);
    if points.iter().any(|p| radius_sq(p) != r2) {
        return Err(fail("ring vertices do not fit one exact circle"));
    }
    let radius = crate::profile::sqrt(&r2)
        .filter(|r| r > &Rational::zero())
        .ok_or_else(|| fail("only positive rational radii are supported"))?;
    Ok((centre, radius))
}

/// Imports and validates the entire closed mesh through M22, then reconstructs.
/// No ground-truth recipe or dimensions are input. All fits are exact on vertices;
/// the faceted surface differs from its fitted cylinder/cone between vertices.
/// Only one connected, outward shell, six axis-aligned planar regions and a
/// two-ring cylindrical/conical through-bore are supported. No triangles are dropped.
pub fn reconstruct(
    bytes: &[u8],
    scale: &Rational,
    dec: &mut Decider,
) -> Result<Reconstruction, Refusal> {
    let mesh = ClosedMesh::read_with_limit(bytes, scale, MAX_TRIANGLES, dec)?;
    if mesh.shells().len() != 1 || !mesh.shells()[0].outward {
        return Err(fail("one outward connected shell is required"));
    }
    let bounds = mesh.bounds().ok_or_else(|| fail("empty source"))?;
    if bounds.iter().any(|b| b[0] >= b[1]) {
        return Err(fail("positive box extents required"));
    }
    let vertices: Vec<_> = (0..mesh.vertex_count())
        .map(|i| mesh.vertex(i as u32))
        .collect();
    let triangles = mesh.triangles();
    let labels: Vec<_> = triangles
        .iter()
        .map(|t| {
            (0..6)
                .find(|&f| {
                    t.iter()
                        .all(|&v| vertices[v as usize][f / 2] == bounds[f / 2][f % 2])
                })
                .unwrap_or(6)
        })
        .collect();
    // Deterministic region growing on actual M22 triangle adjacency, not normal
    // clustering that could join disconnected patches across a missing feature.
    let mut seen = vec![false; triangles.len()];
    let mut regions = [0; 7];
    let mut counts = [0; 7];
    for seed in 0..triangles.len() {
        if seen[seed] {
            continue;
        }
        let label = labels[seed];
        regions[label] += 1;
        let mut queue = vec![seed];
        seen[seed] = true;
        while let Some(t) = queue.pop() {
            dec.charge_work(1)?;
            counts[label] += 1;
            for slot in 0..3 {
                let next = mesh.twin(t as u32, slot) as usize;
                if !seen[next] && labels[next] == label {
                    seen[next] = true;
                    queue.push(next);
                }
            }
        }
    }
    if regions != [1; 7] {
        return Err(fail(
            "expected six connected planes and one connected bore region",
        ));
    }
    let mut rings: [BTreeSet<u32>; 2] = Default::default();
    let mut edges: [BTreeSet<(u32, u32)>; 2] = Default::default();
    for (t, &label) in triangles.iter().zip(&labels) {
        if label == 6 {
            for &v in t {
                let z = &vertices[v as usize][2];
                let level = if z == &bounds[2][0] {
                    0
                } else if z == &bounds[2][1] {
                    1
                } else {
                    return Err(fail("bore must have exactly two end rings"));
                };
                rings[level].insert(v);
            }
            for j in 0..3 {
                let (a, b) = (t[j], t[(j + 1) % 3]);
                if vertices[a as usize][2] == vertices[b as usize][2] {
                    let level = usize::from(vertices[a as usize][2] == bounds[2][1]);
                    edges[level].insert((a.min(b), a.max(b)));
                }
            }
        }
    }
    let fits: Vec<_> = rings
        .iter()
        .map(|ring| {
            circle(
                &ring
                    .iter()
                    .map(|&v| vertices[v as usize].clone())
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Result<_, _>>()?;
    if fits[0].0 != fits[1].0 {
        return Err(fail("bore rings must have the same fitted centre"));
    }
    let centre = fits[0].0.clone();
    let radii = [fits[0].1.clone(), fits[1].1.clone()];
    for r in &radii {
        for a in 0..2 {
            if &centre[a] - r <= bounds[a][0] || &centre[a] + r >= bounds[a][1] {
                return Err(fail("bore must lie strictly within all side planes"));
            }
        }
    }
    let mut rays: [BTreeMap<[Rational; 2], u32>; 2] = Default::default();
    for level in 0..2 {
        for &v in &rings[level] {
            let ray = std::array::from_fn(|a| {
                (&vertices[v as usize][a] - &centre[a])
                    .checked_div(&radii[level])
                    .unwrap()
            });
            if rays[level].insert(ray, v).is_some() {
                return Err(fail("duplicate bore ray"));
            }
        }
    }
    if rays[0].keys().ne(rays[1].keys()) {
        return Err(fail("end rings need matching angular samples"));
    }
    let mut adjacency: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for &(a, b) in &edges[0] {
        adjacency.entry(a).or_default().push(b);
        adjacency.entry(b).or_default().push(a);
    }
    if adjacency.len() != rings[0].len() || adjacency.values().any(|v| v.len() != 2) {
        return Err(fail("bore boundary is not a simple cycle"));
    }
    let start = *rings[0].first().unwrap();
    let mut cycle = vec![start];
    let mut previous = start;
    let mut current = adjacency[&start][0];
    while current != start {
        if cycle.contains(&current) {
            return Err(fail("multiple bore boundary cycles"));
        }
        cycle.push(current);
        let next = adjacency[&current]
            .iter()
            .find(|&&v| v != previous)
            .copied()
            .unwrap();
        previous = current;
        current = next;
    }
    if cycle.len() != rings[0].len() {
        return Err(fail("unmatched bore boundary component"));
    }
    let ray_for = |v: u32| {
        std::array::from_fn(|a| {
            (&vertices[v as usize][a] - &centre[a])
                .checked_div(&radii[0])
                .unwrap()
        })
    };
    let first_sign = cross(&ray_for(cycle[0]), &ray_for(cycle[1])).sign();
    if first_sign == catallax_numbers::Sign::Zero {
        return Err(fail("ambiguous angular span"));
    }
    let mut accepted = BTreeSet::new();
    let mut top_edges = BTreeSet::new();
    for i in 0..cycle.len() {
        dec.charge_work(triangles.len() as u64)?;
        let a = cycle[i];
        let b = cycle[(i + 1) % cycle.len()];
        if cross(&ray_for(a), &ray_for(b)).sign() != first_sign {
            return Err(fail("bore cycle must wind once convexly about its centre"));
        }
        let c = rays[1][&ray_for(a)];
        let d = rays[1][&ray_for(b)];
        top_edges.insert((c.min(d), c.max(d)));
        let corners = [a, b, c, d];
        let mut strip = Vec::new();
        let mut used = BTreeSet::new();
        for (j, t) in triangles.iter().enumerate() {
            if labels[j] == 6 && t.iter().all(|v| corners.contains(v)) {
                strip.push(j);
                used.extend(t.iter().copied());
            }
        }
        if strip.len() != 2 || used.len() != 4 {
            return Err(fail("bore strip must be a complete two-triangle quad"));
        }
        for j in strip {
            if !accepted.insert(j) {
                return Err(fail("ambiguous bore strip assignment"));
            }
        }
    }
    if accepted.len() != counts[6] || top_edges != edges[1] {
        return Err(fail("unmatched bore triangles or end edges"));
    }
    // Every inter-region cap edge must be on the rectangular perimeter or on the
    // recognized bore ring. Combined with M22 closedness/nonintersection and one
    // connected region per plane, this excludes missing caps, islands and extra holes.
    for (i, t) in triangles.iter().enumerate() {
        if labels[i] == 4 || labels[i] == 5 {
            let level = labels[i] - 4;
            for slot in 0..3 {
                if labels[mesh.twin(i as u32, slot) as usize] != labels[i] {
                    let (a, b) = (t[slot], t[(slot + 1) % 3]);
                    let pair = (a.min(b), a.max(b));
                    let perimeter = (0..4).any(|f| {
                        vertices[a as usize][f / 2] == bounds[f / 2][f % 2]
                            && vertices[b as usize][f / 2] == bounds[f / 2][f % 2]
                    });
                    if !perimeter && !edges[level].contains(&pair) {
                        return Err(fail("unmatched cap boundary"));
                    }
                }
            }
        }
    }
    let mut outer_pose = Pose::identity();
    outer_pose.translation.0 = bounds.clone().map(|b| Number::new(b[0].clone()));
    let mut bore_pose = Pose::identity();
    bore_pose.translation.0 =
        [centre[0].clone(), centre[1].clone(), bounds[2][0].clone()].map(Number::new);
    let height = &bounds[2][1] - &bounds[2][0];
    let bore = if radii[0] == radii[1] {
        Operation::Cylinder {
            radius: radii[0].clone().into(),
            height: height.clone().into(),
            pose: bore_pose,
        }
    } else {
        Operation::Cone {
            base_radius: radii[0].clone().into(),
            top_radius: radii[1].clone().into(),
            height: height.into(),
            pose: bore_pose,
        }
    };
    let recipe = Recipe::new(
        "mesh reconstruction",
        vec![
            Step {
                id: id("envelope"),
                op: Operation::Block {
                    size: bounds.clone().map(|b| (&b[1] - &b[0]).into()),
                    pose: outer_pose,
                },
            },
            Step {
                id: id("bore"),
                op: bore,
            },
            Step {
                id: id("part"),
                op: Operation::Boolean {
                    kind: BooleanKind::Subtract,
                    target: id("envelope"),
                    tool: id("bore"),
                },
            },
        ],
        id("part"),
    );
    let brep = catallax_features::evaluate(&recipe, dec).map_err(|e| match e {
        catallax_features::BuildError::Geometry { refusal, .. } => refusal,
        _ => fail("invalid inferred recipe"),
    })?;
    if !brep.model.check(brep.result).is_empty() {
        return Err(fail("reconstructed topology failed structural validation"));
    }
    let field = Plan::compile_checked(&recipe, dec)?;
    let mut result = Reconstruction {
        source_sha256: catallax_import::sha256_hex(bytes),
        fidelity: SampledDeviation {
            scan_to_part_mm: 0.,
            part_to_scan_mm: 0.,
            scan_samples: 0,
            part_samples: 0,
        },
        source_scale: scale.clone(),
        recipe,
        brep,
        field,
        region_triangles: counts,
        centre,
        radii,
        bounds,
    };
    result.fidelity = sampled_deviation(
        bytes,
        &result,
        Sampling {
            triangle_steps: 4,
            surface_steps: 8,
            angles: 128,
        },
        dec,
    )?;
    Ok(result)
}
