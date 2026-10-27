//! Watertight field meshes on a conforming Kuhn tetrahedral grid.
//!
//! Interval exclusion visits every possible boundary cell. A mixed cell has a
//! true boundary witness between differently classified corners. Every deferred
//! cell must neighbour a mixed cell, so both boundary-to-mesh and mesh-to-boundary
//! distances are at most four cell edges. Binary32 coordinates are admitted only
//! when exact, so STL rounding cannot invalidate the bound.

use crate::{Instruction, Location, Plan, Workspace, geometry, lattice::BoxState, metric::Kind};
use catallax_numbers::{Decider, Interval, Rational, Refusal};
use catallax_recipe::BooleanKind;
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

pub mod brep;
mod membership;

type Index = [u32; 3];
type Terms = Vec<([u32; 3], Interval)>;

/// Print-output controls. These do not change the model or its tolerances.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    /// Two-sided distance limit in millimetres.
    pub tolerance: Rational,
    /// Maximum octree cells visited across all refinement attempts.
    pub max_cells: usize,
    /// Maximum mesh vertices.
    pub max_vertices: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            tolerance: Rational::from_ratio(1, 100).unwrap(),
            max_cells: 2_000_000,
            max_vertices: 1_000_000,
        }
    }
}

/// A closed, consistently oriented triangular field approximation.
#[derive(Clone, Debug)]
pub struct Mesh {
    /// Exactly representable binary32 millimetre coordinates, shared by facets.
    pub positions: Vec<[f32; 3]>,
    /// Outward-oriented facet vertex indices.
    pub triangles: Vec<[usize; 3]>,
    /// Proven two-sided distance bound, including output representation.
    pub max_deviation: Rational,
    /// Total visited octree cells, including unsuccessful refinements.
    pub cells: usize,
}
/// Bytes gathered before one write. Unbuffered sinks such as files otherwise
/// pay a system call per field.
const CHUNK: usize = 1 << 16;
fn drain(out: &mut impl Write, chunk: &mut Vec<u8>, limit: usize) -> std::io::Result<()> {
    if chunk.len() >= limit {
        out.write_all(chunk)?;
        chunk.clear();
    }
    Ok(())
}
fn stl_normal(points: [[f64; 3]; 3]) -> std::io::Result<[f32; 3]> {
    let invalid = || {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "STL facet has non-finite coordinates or no finite nonzero normal",
        )
    };
    if points.iter().flatten().any(|value| !value.is_finite()) {
        return Err(invalid());
    }
    let normalise = |normal: [f64; 3]| {
        let length = crate::metric::Scalar::sqrt(&normal.iter().map(|n| n * n).sum::<f64>());
        normal.map(|component| (component / length) as f32)
    };
    let a: [f64; 3] = std::array::from_fn(|i| points[1][i] - points[0][i]);
    let b: [f64; 3] = std::array::from_fn(|i| points[2][i] - points[0][i]);
    let normal = normalise(std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        a[j] * b[k] - a[k] * b[j]
    }));
    if normal.iter().all(|component| component.is_finite()) {
        // Preserve existing finite field normals, including their exact bytes.
        return Ok(normal);
    }
    // Subtract the encoded coordinates exactly before forming the cross product.
    // Binary32 differences/products have exponents within binary64's range, so a
    // nonzero exact component survives conversion and its squared norm is finite.
    let points = points.map(|point| {
        point.map(|value| Rational::from_f64_exact(value).expect("finite encoded coordinate"))
    });
    let a: [Rational; 3] = std::array::from_fn(|i| &points[1][i] - &points[0][i]);
    let b: [Rational; 3] = std::array::from_fn(|i| &points[2][i] - &points[0][i]);
    let normal = normalise(std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        (&(&a[j] * &b[k]) - &(&a[k] * &b[j])).to_f64_nearest()
    }));
    if normal.iter().all(|component| component.is_finite()) {
        Ok(normal)
    } else {
        Err(invalid())
    }
}

impl Mesh {
    /// The number of edge-connected components. A closed mesh with no triangles has none.
    pub fn component_count(&self) -> usize {
        let mut parent: Vec<usize> = (0..self.positions.len()).collect();
        fn root(parent: &mut [usize], mut v: usize) -> usize {
            while parent[v] != v {
                parent[v] = parent[parent[v]];
                v = parent[v];
            }
            v
        }
        for &[a, b, c] in &self.triangles {
            for (x, y) in [(a, b), (b, c)] {
                let (x, y) = (root(&mut parent, x), root(&mut parent, y));
                parent[x] = y;
            }
        }
        let used: BTreeSet<usize> = self.triangles.iter().flatten().copied().collect();
        used.into_iter()
            .map(|v| root(&mut parent, v))
            .collect::<BTreeSet<_>>()
            .len()
    }
    /// Writes a binary STL with finite advisory normals. All topology checks
    /// finish before publication. Invalid coordinates or degenerate facets
    /// return [`std::io::ErrorKind::InvalidData`], never NaN normal bytes.
    pub fn write_stl(&self, out: &mut impl Write) -> std::io::Result<()> {
        let count = u32::try_from(self.triangles.len())
            .map_err(|_| std::io::Error::other("too many STL triangles"))?;
        let mut chunk = Vec::with_capacity(CHUNK + 50);
        let mut header = [0u8; 80];
        let label = b"catallax field mesh; millimetres";
        header[..label.len()].copy_from_slice(label);
        chunk.extend_from_slice(&header);
        chunk.extend_from_slice(&count.to_le_bytes());
        for triangle in &self.triangles {
            let p = triangle.map(|i| self.positions[i].map(f64::from));
            for component in stl_normal(p)? {
                chunk.extend_from_slice(&component.to_le_bytes());
            }
            for &vertex in triangle {
                for coordinate in self.positions[vertex] {
                    chunk.extend_from_slice(&coordinate.to_le_bytes());
                }
            }
            chunk.extend_from_slice(&[0u8; 2]);
            drain(out, &mut chunk, CHUNK)?;
        }
        drain(out, &mut chunk, 1)
    }
    /// Writes the same rounded vertices and facets for display.
    pub fn write_obj(&self, out: &mut impl Write) -> std::io::Result<()> {
        let mut chunk = Vec::with_capacity(CHUNK + 256);
        writeln!(chunk, "# catallax field print mesh; millimetres")?;
        writeln!(chunk, "# proven_two_sided_deviation {}", self.max_deviation)?;
        for p in &self.positions {
            writeln!(chunk, "v {} {} {}", p[0] as f64, p[1] as f64, p[2] as f64)?;
            drain(out, &mut chunk, CHUNK)?;
        }
        for t in &self.triangles {
            writeln!(chunk, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1)?;
            drain(out, &mut chunk, CHUNK)?;
        }
        drain(out, &mut chunk, 1)
    }
}
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn refusal(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}
fn charge(dec: &mut Decider) -> Result<(), Refusal> {
    dec.check_deadline()
}
fn range(terms: &Terms, xyz: &[Interval; 3]) -> Interval {
    terms.iter().fold(
        Interval::point(0.0).unwrap(),
        |sum, (powers, coefficient)| {
            let mut term = *coefficient;
            for axis in 0..3 {
                for _ in 0..powers[axis] {
                    term = term * xyz[axis];
                }
            }
            sum + term
        },
    )
}
struct Primitive {
    kind: Kind,
    matrix: [[Interval; 3]; 3],
    translation: [Interval; 3],
    parameters: [Interval; 3],
}
fn state(bounds: &[Interval]) -> BoxState {
    if bounds.iter().any(|b| b.lo() > 0.0) {
        BoxState::Outside
    } else if bounds.iter().all(|b| b.hi() < 0.0) {
        BoxState::Inside
    } else {
        BoxState::Unknown
    }
}
impl Primitive {
    fn classify(&self, world: &[Interval; 3]) -> BoxState {
        let local: [Interval; 3] = std::array::from_fn(|i| {
            (0..3).fold(self.translation[i], |sum, j| {
                sum + self.matrix[i][j] * world[j]
            })
        });
        let zero = Interval::point(0.0).unwrap();
        match self.kind {
            Kind::Block => state(&[
                zero - local[0],
                local[0] - self.parameters[0],
                zero - local[1],
                local[1] - self.parameters[1],
                zero - local[2],
                local[2] - self.parameters[2],
            ]),
            Kind::Round => {
                let radius = self.parameters[0]
                    + (self.parameters[1] - self.parameters[0]) * local[2] / self.parameters[2];
                state(&[
                    zero - local[2],
                    local[2] - self.parameters[2],
                    local[0].square() + local[1].square() - radius.square(),
                ])
            }
            Kind::Sphere => {
                let value = local[0].square() + local[1].square() + local[2].square()
                    - self.parameters[0].square();
                state(&[value])
            }
            Kind::Torus => {
                let radial = (local[0].square() + local[1].square())
                    .sqrt()
                    .expect("sum of squares");
                let value = (radial - self.parameters[0]).square() + local[2].square()
                    - self.parameters[1].square();
                state(&[value])
            }
        }
    }
}
enum Definition {
    Primitive(Box<Primitive>),
    Cell(Vec<Terms>),
    Lattice,
}
pub(crate) struct Classifier<'a> {
    plan: &'a Plan,
    definitions: Vec<Option<Definition>>,
    states: Vec<Option<BoxState>>,
    pending: Vec<usize>,
}
impl<'a> Classifier<'a> {
    pub(crate) fn new(plan: &'a Plan, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut definitions = Vec::with_capacity(plan.leaves.len());
        for (i, leaf) in plan.leaves.iter().enumerate() {
            // Compaction retains source leaves for face history. Keep their
            // indices, but prepare only definitions referenced by the live tape.
            if plan.leaf_registers[i].is_none() {
                definitions.push(None);
                continue;
            }
            match leaf {
                geometry::Leaf::Primitive(primitive) => {
                    let inverse = primitive.inverse_placement();
                    let mut matrix = [[Interval::point(0.0).unwrap(); 3]; 3];
                    for (i, row) in matrix.iter_mut().enumerate() {
                        for (j, value) in row.iter_mut().enumerate() {
                            *value =
                                geometry::enclose(&inverse.matrix()[i][j], 192, dec)?.to_interval();
                        }
                    }
                    let translation = [
                        &inverse.translation.x,
                        &inverse.translation.y,
                        &inverse.translation.z,
                    ];
                    let mut offsets = [Interval::point(0.0).unwrap(); 3];
                    for i in 0..3 {
                        offsets[i] = geometry::enclose(translation[i], 192, dec)?.to_interval();
                    }
                    definitions.push(Some(Definition::Primitive(Box::new(Primitive {
                        kind: primitive.kind,
                        matrix,
                        translation: offsets,
                        parameters: primitive.parameters.each_ref().map(Interval::from_rational),
                    }))));
                }
                geometry::Leaf::Cell(_) => {
                    let mut supports = Vec::new();
                    for constraint in leaf.constraints() {
                        let mut terms = Vec::new();
                        for (monomial, coefficient) in constraint.polynomial.terms() {
                            terms.push((
                                monomial.exponents().try_into().expect("three coordinates"),
                                geometry::enclose(&coefficient, 192, dec)?.to_interval(),
                            ));
                        }
                        supports.push(terms);
                    }
                    definitions.push(Some(Definition::Cell(supports)));
                }
                geometry::Leaf::Lattice(_) => definitions.push(Some(Definition::Lattice)),
                _ => {
                    return Err(refusal(
                        "field.print-coverage",
                        "print interval coverage requires analytic solid cells or periodic beams",
                    ));
                }
            }
        }
        Ok(Self {
            plan,
            definitions,
            states: vec![None; plan.instructions.len()],
            pending: Vec::with_capacity(plan.instructions.len()),
        })
    }
    fn classify(&mut self, xyz: &[Interval; 3], dec: &mut Decider) -> Result<BoxState, Refusal> {
        self.classify_with(xyz, dec, |_| None)
    }
    pub(crate) fn classify_with(
        &mut self,
        xyz: &[Interval; 3],
        dec: &mut Decider,
        mut exact_leaf: impl FnMut(usize) -> Option<BoxState>,
    ) -> Result<BoxState, Refusal> {
        use BoxState::*;
        charge(dec)?;
        self.states.fill(None);
        self.pending.clear();
        self.pending.push(self.plan.root);
        while let Some(&node) = self.pending.last() {
            if self.states[node].is_some() {
                self.pending.pop();
                continue;
            }
            let state = match self.plan.instructions[node] {
                Instruction::Empty => Outside,
                Instruction::Leaf(i) => {
                    if let Some(state) = exact_leaf(i) {
                        state
                    } else {
                        match self.definitions[i].as_ref().expect("live leaf definition") {
                            Definition::Primitive(primitive) => primitive.classify(xyz),
                            Definition::Lattice => {
                                let geometry::Leaf::Lattice(lattice) = &self.plan.leaves[i] else {
                                    unreachable!("compiled lattice")
                                };
                                lattice.box_state(xyz, dec)?
                            }
                            Definition::Cell(supports) => {
                                let mut inside = true;
                                let mut outside = false;
                                for support in supports {
                                    let bound = range(support, xyz);
                                    outside |= bound.lo() > 0.0;
                                    inside &= bound.hi() < 0.0;
                                }
                                if outside {
                                    Outside
                                } else if inside {
                                    Inside
                                } else {
                                    Unknown
                                }
                            }
                        }
                    }
                }
                Instruction::Boolean { kind, a, b } => {
                    let Some(left) = self.states[a] else {
                        self.pending.push(a);
                        continue;
                    };
                    match (kind, left) {
                        (BooleanKind::Union, Inside) => Inside,
                        (BooleanKind::Intersect, Outside) | (BooleanKind::Subtract, Outside) => {
                            Outside
                        }
                        _ => {
                            let Some(right) = self.states[b] else {
                                self.pending.push(b);
                                continue;
                            };
                            match (kind, left, right) {
                                (BooleanKind::Union, _, Inside) => Inside,
                                (BooleanKind::Union, Outside, Outside) => Outside,
                                (BooleanKind::Intersect, _, Outside) => Outside,
                                (BooleanKind::Intersect, Inside, Inside) => Inside,
                                (BooleanKind::Subtract, _, Inside) => Outside,
                                (BooleanKind::Subtract, Inside, Outside) => Inside,
                                _ => Unknown,
                            }
                        }
                    }
                }
            };
            self.states[node] = Some(state);
            self.pending.pop();
        }
        Ok(self.states[self.plan.root].expect("classified root"))
    }
    fn inside(
        &mut self,
        point: &[Rational; 3],
        workspace: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        match self.classify(&point.each_ref().map(Interval::from_rational), dec)? {
            BoxState::Inside => Ok(true),
            BoxState::Outside => Ok(false),
            BoxState::Unknown => match self.corner_location(point, workspace, dec)? {
                Location::Inside => Ok(true),
                Location::Outside | Location::Boundary => Ok(false),
                Location::Indeterminate => Err(refusal(
                    "field.print-classification",
                    "an exact grid-corner classification is unresolved",
                )),
            },
        }
    }
}
pub(crate) fn bounds(plan: &Plan, dec: &mut Decider) -> Result<[[Rational; 2]; 3], Refusal> {
    let mut tape: Vec<Option<[[Rational; 2]; 3]>> = Vec::with_capacity(plan.instructions.len());
    for instruction in &plan.instructions {
        let bound = match *instruction {
            Instruction::Empty => Some(std::array::from_fn(|_| [q(0), q(0)])),
            Instruction::Leaf(i) => plan.leaves[i]
                .display_box(dec)?
                .map(|b| b.map(|x| [x.lo(), x.hi()])),
            Instruction::Boolean { kind, a, b } => match kind {
                BooleanKind::Subtract => tape[a].clone(),
                BooleanKind::Union => match (&tape[a], &tape[b]) {
                    (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                        [
                            a[i][0].clone().min(b[i][0].clone()),
                            a[i][1].clone().max(b[i][1].clone()),
                        ]
                    })),
                    _ => None,
                },
                BooleanKind::Intersect => match (&tape[a], &tape[b]) {
                    (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                        [
                            a[i][0].clone().max(b[i][0].clone()),
                            a[i][1].clone().min(b[i][1].clone()),
                        ]
                    })),
                    (a, b) => a.clone().or_else(|| b.clone()),
                },
            },
        };
        tape.push(bound);
    }
    tape[plan.root].clone().ok_or_else(|| {
        refusal(
            "field.print-unbounded",
            "print output needs a bounded region",
        )
    })
}
struct Builder<'a> {
    classifier: Classifier<'a>,
    workspace: Workspace,
    origin: [Rational; 3],
    edge: Rational,
    corner_signs: BTreeMap<Index, bool>,
    edge_vertices: BTreeMap<(Index, Index), usize>,
    grid_vertices: Vec<Index>,
    mesh: Mesh,
    options: &'a Options,
}
impl Builder<'_> {
    fn point(&self, index: Index) -> [Rational; 3] {
        std::array::from_fn(|i| &self.origin[i] + &(&self.edge * &q(index[i] as i64)))
    }
    fn corners(
        &mut self,
        index: Index,
        dec: &mut Decider,
    ) -> Result<([Index; 8], [bool; 8]), Refusal> {
        let corners: [Index; 8] = std::array::from_fn(|bits| {
            std::array::from_fn(|i| index[i] + ((bits >> i) & 1) as u32)
        });
        let mut signs = [false; 8];
        for (i, corner) in corners.iter().enumerate() {
            if let Some(&sign) = self.corner_signs.get(corner) {
                signs[i] = sign;
            } else {
                signs[i] =
                    self.classifier
                        .inside(&self.point(*corner), &mut self.workspace, dec)?;
                self.corner_signs.insert(*corner, signs[i]);
            }
        }
        Ok((corners, signs))
    }
    fn vertex(&mut self, a: Index, b: Index) -> Result<usize, Refusal> {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(&id) = self.edge_vertices.get(&key) {
            return Ok(id);
        }
        if self.mesh.positions.len() >= self.options.max_vertices {
            return Err(refusal("field.print-budget", "print vertex limit reached"));
        }
        let position: [Rational; 3] = std::array::from_fn(|i| {
            &self.origin[i] + &(&self.edge * &q(a[i] as i64 + b[i] as i64)).mul_pow2(-1)
        });
        let rounded = position.each_ref().map(|r| r.to_f64_nearest() as f32);
        for i in 0..3 {
            if !rounded[i].is_finite()
                || Rational::from_f64_exact(rounded[i] as f64).as_ref() != Some(&position[i])
            {
                return Err(refusal(
                    "field.print-coordinate",
                    "the chosen binary32 print grid is not exactly representable at this placement",
                ));
            }
        }
        let id = self.mesh.positions.len();
        self.mesh.positions.push(rounded);
        self.grid_vertices
            .push(std::array::from_fn(|i| a[i] + b[i]));
        self.edge_vertices.insert(key, id);
        Ok(id)
    }
    fn triangle(&mut self, mut vertices: [usize; 3], direction: [i64; 3]) -> Result<(), Refusal> {
        let p = vertices.map(|i| self.grid_vertices[i]);
        let a: [i64; 3] = std::array::from_fn(|i| p[1][i] as i64 - p[0][i] as i64);
        let b: [i64; 3] = std::array::from_fn(|i| p[2][i] as i64 - p[0][i] as i64);
        let normal: [i64; 3] = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            a[j] * b[k] - a[k] * b[j]
        });
        let dot: i64 = (0..3).map(|i| normal[i] * direction[i]).sum();
        if dot < 0 {
            vertices.swap(1, 2);
        } else if dot == 0 {
            return Err(refusal(
                "field.print-topology",
                "a tetrahedral facet has no outward orientation",
            ));
        }
        self.mesh.triangles.push(vertices);
        Ok(())
    }
    fn polygonise(&mut self, index: Index, dec: &mut Decider) -> Result<bool, Refusal> {
        let (corners, signs) = self.corners(index, dec)?;
        if signs.iter().all(|s| *s) || signs.iter().all(|s| !*s) {
            return Ok(false);
        }
        for [a, b] in [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]] {
            let tetra = [0, 1 << a, (1 << a) | (1 << b), 7];
            let mut inside = [0usize; 4];
            let mut outside = [0usize; 4];
            let (mut ni, mut no) = (0, 0);
            for i in tetra {
                if signs[i] {
                    inside[ni] = i;
                    ni += 1;
                } else {
                    outside[no] = i;
                    no += 1;
                }
            }
            let inside = &inside[..ni];
            let outside = &outside[..no];
            if inside.is_empty() || outside.is_empty() {
                continue;
            }
            let direction: [i64; 3] = std::array::from_fn(|axis| {
                let out: i64 = outside.iter().map(|i| corners[*i][axis] as i64).sum();
                let inn: i64 = inside.iter().map(|i| corners[*i][axis] as i64).sum();
                out * inside.len() as i64 - inn * outside.len() as i64
            });
            if inside.len() == 1 || outside.len() == 1 {
                let (one, others) = if inside.len() == 1 {
                    (inside[0], &outside)
                } else {
                    (outside[0], &inside)
                };
                let vertices = [
                    self.vertex(corners[one], corners[others[0]])?,
                    self.vertex(corners[one], corners[others[1]])?,
                    self.vertex(corners[one], corners[others[2]])?,
                ];
                self.triangle(vertices, direction)?;
            } else {
                let [i, j] = [inside[0], inside[1]];
                let [k, l] = [outside[0], outside[1]];
                let v = [
                    self.vertex(corners[i], corners[k])?,
                    self.vertex(corners[i], corners[l])?,
                    self.vertex(corners[j], corners[k])?,
                    self.vertex(corners[j], corners[l])?,
                ];
                self.triangle([v[0], v[1], v[3]], direction)?;
                self.triangle([v[0], v[3], v[2]], direction)?;
            }
        }
        Ok(true)
    }
}
fn check_topology(mesh: &Mesh) -> Result<(), Refusal> {
    let mut edges: BTreeMap<(usize, usize), (usize, i64)> = BTreeMap::new();
    let mut links: Vec<Vec<(usize, usize)>> = vec![Vec::new(); mesh.positions.len()];
    for &[a, b, c] in &mesh.triangles {
        for (v, x, y) in [(a, b, c), (b, c, a), (c, a, b)] {
            links[v].push((x, y));
            let (key, sense) = if v < x { ((v, x), 1) } else { ((x, v), -1) };
            let entry = edges.entry(key).or_default();
            entry.0 += 1;
            entry.1 += sense;
        }
    }
    if edges
        .values()
        .any(|&(count, sense)| count != 2 || sense != 0)
    {
        return Err(refusal(
            "field.print-topology",
            "the field mesh has an open or inconsistently oriented edge",
        ));
    }
    for link in links {
        let Some(&(start, _)) = link.first() else {
            return Err(refusal(
                "field.print-topology",
                "an unused print vertex remains",
            ));
        };
        let mut adjacency: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (a, b) in link {
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
        }
        if adjacency.values().any(|n| n.len() != 2) {
            return Err(refusal(
                "field.print-topology",
                "a print vertex link is not a cycle",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut pending = vec![start];
        while let Some(v) = pending.pop() {
            if seen.insert(v) {
                pending.extend(&adjacency[&v]);
            }
        }
        if seen.len() != adjacency.len() {
            return Err(refusal(
                "field.print-topology",
                "a print vertex joins disconnected sheets",
            ));
        }
    }
    Ok(())
}

/// Builds the print output of a plan whose result must be one connected solid, such as
/// a lattice. A closed mesh with several components is refused as
/// `field.print-disconnected` rather than published.
pub fn mesh_plan_connected(
    plan: &Plan,
    options: &Options,
    dec: &mut Decider,
) -> Result<Mesh, Refusal> {
    let mesh = mesh_plan(plan, options, dec)?;
    match mesh.component_count() {
        1 => Ok(mesh),
        n => Err(refusal(
            "field.print-disconnected",
            format!("the print output has {n} connected components, expected one"),
        )),
    }
}

/// Builds a watertight output from the field alone. No B-rep builder is called.
pub fn mesh_plan(plan: &Plan, options: &Options, dec: &mut Decider) -> Result<Mesh, Refusal> {
    if options.tolerance <= q(0) || options.max_cells == 0 || options.max_vertices < 4 {
        return Err(refusal(
            "field.print-options",
            "positive tolerance and cell limit, and at least four vertices required",
        ));
    }
    if plan.is_proven_empty() {
        return Ok(Mesh {
            positions: vec![],
            triangles: vec![],
            max_deviation: q(0),
            cells: 0,
        });
    }
    let bounds = bounds(plan, dec)?;
    let mut side = q(1);
    let extent = (0..3).map(|i| &bounds[i][1] - &bounds[i][0]).max().unwrap();
    while side <= &extent + &q(2) {
        side = side.mul_pow2(1);
    }
    let mut edge = side.clone();
    let mut depth = 0u32;
    while &edge * &q(4) > options.tolerance {
        edge = edge.mul_pow2(-1);
        depth += 1;
        if depth > 29 {
            return Err(refusal(
                "field.print-grid",
                "print grid index range exhausted",
            ));
        }
    }
    let origin: [Rational; 3] = std::array::from_fn(|i| {
        &(&edge * &Rational::from_integer((&bounds[i][0] / &edge).floor())) - &edge
    });
    while (0..3).any(|i| &origin[i] + &side <= bounds[i][1]) {
        side = side.mul_pow2(1);
        depth += 1;
        if depth > 29 {
            return Err(refusal(
                "field.print-grid",
                "print grid index range exhausted",
            ));
        }
    }
    let classifier = Classifier::new(plan, dec)?;
    let mut builder = Builder {
        classifier,
        workspace: Workspace::default(),
        origin,
        edge,
        corner_signs: BTreeMap::new(),
        edge_vertices: BTreeMap::new(),
        grid_vertices: Vec::new(),
        mesh: Mesh {
            positions: vec![],
            triangles: vec![],
            max_deviation: q(0),
            cells: 0,
        },
        options,
    };
    loop {
        let count = 1u32 << depth;
        let mut pending = vec![([0u32; 3], count)];
        let mut unknown = Vec::new();
        while let Some((index, width)) = pending.pop() {
            builder.mesh.cells += 1;
            if builder.mesh.cells > options.max_cells {
                return Err(refusal(
                    "field.print-budget",
                    "print interval refinement cell limit reached",
                ));
            }
            let lo = builder.point(index);
            let hi: [Rational; 3] =
                std::array::from_fn(|i| &lo[i] + &(&builder.edge * &q(width as i64)));
            let cover = std::array::from_fn(|i| {
                Interval::new(
                    Interval::from_rational(&lo[i]).lo(),
                    Interval::from_rational(&hi[i]).hi(),
                )
                .unwrap()
            });
            if builder.classifier.classify(&cover, dec)? != BoxState::Unknown {
                continue;
            }
            if width == 1 {
                unknown.push(index);
            } else {
                let half = width / 2;
                for bits in 0..8 {
                    pending.push((
                        std::array::from_fn(|i| index[i] + half * ((bits >> i) & 1)),
                        half,
                    ));
                }
            }
        }
        let mut mixed = BTreeSet::new();
        for &cell in &unknown {
            if builder.polygonise(cell, dec)? {
                mixed.insert(cell);
            }
        }
        let covered = unknown.iter().all(|cell| {
            (-1i64..=1).any(|x| {
                (-1i64..=1).any(|y| {
                    (-1i64..=1).any(|z| {
                        let candidate =
                            [cell[0] as i64 + x, cell[1] as i64 + y, cell[2] as i64 + z];
                        candidate.iter().all(|c| *c >= 0 && *c < count as i64)
                            && mixed.contains(&candidate.map(|n| n as u32))
                    })
                })
            })
        });
        if covered && !builder.mesh.triangles.is_empty() {
            check_topology(&builder.mesh)?;
            builder.mesh.max_deviation = &builder.edge * &q(4);
            return Ok(builder.mesh);
        }
        depth += 1;
        if depth > 29 {
            return Err(refusal(
                "field.print-unresolved",
                "boundary coverage did not settle on the print grid",
            ));
        }
        builder.edge = builder.edge.mul_pow2(-1);
        builder.mesh.positions.clear();
        builder.mesh.triangles.clear();
        builder.edge_vertices.clear();
        builder.corner_signs.clear();
        builder.grid_vertices.clear();
    }
}
