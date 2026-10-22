//! Adaptive display directly from field definitions and exact trim predicates.
//!
//! Analytic parameter cells cover every candidate support within a certified
//! result bound. Interval constraint signs prove constant trim membership.
//! Cells crossing trims shrink until a surviving boundary witness certifies
//! their diameter, or a neighbouring mesh vertex covers the entire cell. An
//! unresolved cell refuses; point samples never establish absence or coverage.

use crate::{FaceId, FaceState, Instruction, Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Interval, Rational, Refusal};
use catallax_recipe::BooleanKind;
use catallax_shapes::RPoint3;
use std::{collections::BTreeMap, io::Write};

mod chart;
use chart::Chart;

/// Deterministic display limits, in millimetres and counts.
#[derive(Clone, Debug)]
pub struct Options {
    /// Two-sided distance limit between the mesh and the defined boundary.
    pub tolerance: Rational,
    /// Maximum distinct binary64 display vertices.
    pub max_vertices: usize,
    /// Maximum visited parameter cells, including rejected and deferred cells.
    pub max_cells: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            tolerance: Rational::from_ratio(1, 100).unwrap(),
            max_vertices: 1_000_000,
            max_cells: 2_000_000,
        }
    }
}
/// Triangle tagged with a field registry handle, not a B-rep face.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Triangle {
    /// Indices into [`Mesh::positions`].
    pub vertices: [usize; 3],
    /// A surviving named field support at this patch's certified witness.
    pub face: FaceId,
}
/// A display approximation, never an exact solid or STEP export.
#[derive(Clone, Debug, PartialEq)]
pub struct Mesh {
    /// Binary64 millimetre coordinates.
    pub positions: Vec<[f64; 3]>,
    /// Definition-native triangles. Adaptive trim cells need not sew watertight.
    pub triangles: Vec<Triangle>,
    /// Proven two-sided distance upper bound to the defined boundary.
    pub max_deviation: Rational,
    /// Deterministic number of visited parameter cells.
    pub cells: usize,
}
impl Mesh {
    /// Writes the actual display artifact, with field registry picking groups.
    pub fn write_obj(&self, out: &mut impl Write) -> std::io::Result<()> {
        writeln!(out, "# catallax field display mesh; millimetres")?;
        writeln!(out, "# proven_two_sided_deviation {}", self.max_deviation)?;
        for [x, y, z] in &self.positions {
            writeln!(out, "v {x} {y} {z}")?;
        }
        let mut previous = None;
        for t in &self.triangles {
            if previous != Some(t.face) {
                writeln!(out, "g field_face_{}", t.face.0)?;
                previous = Some(t.face);
            }
            let [a, b, c] = t.vertices.map(|v| v + 1);
            writeln!(out, "f {a} {b} {c}")?;
        }
        Ok(())
    }
}
fn failure(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}

#[derive(Clone)]
struct Cell {
    range: [[Rational; 2]; 2],
    depth: u32,
}
impl Cell {
    fn centre(&self) -> [Rational; 2] {
        self.range.each_ref().map(|r| (&r[0] + &r[1]).mul_pow2(-1))
    }
    fn span(&self) -> [Rational; 2] {
        self.range.each_ref().map(|r| &r[1] - &r[0])
    }
    fn corners(&self) -> [[Rational; 2]; 4] {
        let [u, v] = &self.range;
        [
            [u[0].clone(), v[0].clone()],
            [u[1].clone(), v[0].clone()],
            [u[1].clone(), v[1].clone()],
            [u[0].clone(), v[1].clone()],
        ]
    }
    fn split(&self, axis: usize) -> [Self; 2] {
        let mid = (&self.range[axis][0] + &self.range[axis][1]).mul_pow2(-1);
        let mut a = self.clone();
        let mut b = self.clone();
        a.range[axis][1] = mid.clone();
        b.range[axis][0] = mid;
        a.depth += 1;
        b.depth += 1;
        [a, b]
    }
}
struct Constraint {
    terms: Vec<([u32; 3], Interval)>,
    support: catallax_numbers::RealPoly,
}
impl Constraint {
    fn range(&self, p: &[Interval; 3]) -> Interval {
        self.terms
            .iter()
            .fold(Interval::point(0.0).unwrap(), |sum, (powers, c)| {
                let mut t = *c;
                for i in 0..3 {
                    let mut exponent = powers[i];
                    while exponent >= 2 {
                        t = t * p[i].square();
                        exponent -= 2;
                    }
                    if exponent == 1 {
                        t = t * p[i];
                    }
                }
                sum + t
            })
    }
}
struct Patch {
    chart: Chart,
    faces: Vec<usize>,
    coincident: Vec<Vec<bool>>,
}
struct Builder<'a> {
    plan: &'a Plan,
    options: &'a Options,
    constraints: Vec<Vec<Constraint>>,
    workspace: Workspace,
    locations: Vec<(Location, bool)>,
    vertices: BTreeMap<[u64; 3], usize>,
    mesh: Mesh,
}

/// Constructs the field's display without calling a B-rep builder or mesher.
/// Every unsupported chart, unresolved trim or exhausted limit returns a stable
/// refusal rather than publishing incomplete geometry or falling back to B-rep.
pub fn mesh_plan(plan: &Plan, options: &Options, dec: &mut Decider) -> Result<Mesh, Refusal> {
    if options.tolerance <= q(0) || options.max_vertices < 4 || options.max_cells == 0 {
        return Err(failure(
            "field.display-options",
            "positive tolerance and cell limit, and at least four vertices required",
        ));
    }
    let mut builder = Builder {
        plan,
        options,
        constraints: Vec::new(),
        workspace: Workspace::default(),
        locations: Vec::new(),
        vertices: BTreeMap::new(),
        mesh: Mesh {
            positions: Vec::new(),
            triangles: Vec::new(),
            max_deviation: q(0),
            cells: 0,
        },
    };
    if plan.is_proven_empty() {
        return Ok(builder.mesh);
    }
    let mut bounds: Vec<Option<Rational>> = Vec::with_capacity(plan.instructions.len());
    for instruction in &plan.instructions {
        let bound = match *instruction {
            Instruction::Empty => Some(q(0)),
            Instruction::Leaf(i) => plan.leaves[i].display_radius(dec)?,
            Instruction::Boolean { kind, a, b } => match kind {
                BooleanKind::Subtract => bounds[a].clone(),
                BooleanKind::Union => match (&bounds[a], &bounds[b]) {
                    (Some(a), Some(b)) => Some(a.max(b).clone()),
                    _ => None,
                },
                BooleanKind::Intersect => match (&bounds[a], &bounds[b]) {
                    (Some(a), Some(b)) => Some(a.min(b).clone()),
                    (a, b) => a.clone().or_else(|| b.clone()),
                },
            },
        };
        bounds.push(bound);
    }
    let bound = bounds[plan.root].clone().ok_or_else(|| {
        failure(
            "field.display-unbounded",
            "result has no certified compact display domain",
        )
    })? + q(1);
    for leaf in &plan.leaves {
        let mut constraints = Vec::new();
        for c in leaf.constraints() {
            let (support, _) = crate::geometry::support_geometry(&c.polynomial, dec)?;
            let mut terms = Vec::new();
            for (m, c) in c.polynomial.terms() {
                terms.push((
                    m.exponents().try_into().expect("three coordinates"),
                    crate::geometry::enclose(&c, 128, dec)?.to_interval(),
                ));
            }
            constraints.push(Constraint { terms, support });
        }
        builder.constraints.push(constraints);
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (i, face) in plan.faces.iter().enumerate() {
        if plan.leaf_registers[face.leaf].is_none() {
            continue;
        }
        let mut found = None;
        for (j, group) in groups.iter().enumerate() {
            if crate::same_support(face, &plan.faces[group[0]], dec)? {
                found = Some(j);
                break;
            }
        }
        if let Some(j) = found {
            groups[j].push(i);
        } else {
            groups.push(vec![i]);
        }
    }
    let mut patches = Vec::new();
    for faces in groups {
        let support = &plan.faces[faces[0]].support;
        let mut coincident = Vec::new();
        for constraints in &builder.constraints {
            let mut same = Vec::new();
            for c in constraints {
                same.push(support.equals(&c.support, dec)?);
            }
            coincident.push(same);
        }
        for chart in Chart::all(&plan.faces[faces[0]].surface, &bound)? {
            patches.push(Patch {
                chart,
                faces: faces.clone(),
                coincident: coincident.clone(),
            });
        }
    }
    let mut stack: Vec<_> = patches
        .iter()
        .enumerate()
        .rev()
        .map(|(i, p)| (i, p.chart.domain.clone()))
        .collect();
    let mut deferred = Vec::new();
    loop {
        while let Some((id, cell)) = stack.pop() {
            dec.check_deadline()?;
            builder.mesh.cells += 1;
            if builder.mesh.cells > options.max_cells || cell.depth > 64 {
                return Err(failure(
                    "field.display-budget",
                    format!(
                        "refinement limit: {} cells, depth {}, support {:?}",
                        builder.mesh.cells, cell.depth, plan.faces[patches[id].faces[0]].origin
                    ),
                ));
            }
            let patch = &patches[id];
            let cover = patch.chart.cover(&cell);
            let point = patch.chart.point(&cell.centre());
            let (location, stable) = builder.classify(patch, &cover, &point, dec)?;
            if matches!(location, Location::Inside | Location::Outside) {
                continue;
            }
            let diameter = cover.iter().fold(q(0), |sum, b| sum + b.hi() - b.lo());
            let small = diameter <= options.tolerance.mul_pow2(-1);
            if stable || small {
                let mut present = None;
                for &face in &patch.faces {
                    let state = plan
                        .face_query_checked(FaceId(face), &point, &mut builder.workspace, dec)?
                        .state;
                    if state == FaceState::OnFace {
                        present = Some(face);
                        break;
                    }
                    if stable && state == FaceState::Indeterminate {
                        break;
                    }
                }
                if let Some(face) = present {
                    let chord = patch.chart.chord(&cell);
                    if (stable && chord <= options.tolerance.mul_pow2(-1)) || small {
                        builder.emit(
                            patch,
                            &cell,
                            face,
                            if stable { chord } else { diameter },
                            dec,
                        )?;
                        continue;
                    }
                } else if stable {
                    // A sign-stable patch has constant trim membership. An exact
                    // absence is valid only if every candidate settles absent.
                    let mut absent = true;
                    for &face in &patch.faces {
                        absent &= plan
                            .face_query_checked(FaceId(face), &point, &mut builder.workspace, dec)?
                            .state
                            == FaceState::NotOnFace;
                    }
                    if absent {
                        continue;
                    }
                }
                if small {
                    deferred.push((id, cell, cover));
                    continue;
                }
            }
            let span = cell.span();
            let axis = usize::from(
                &patch.chart.lipschitz[1] * &span[1] > &patch.chart.lipschitz[0] * &span[0],
            );
            let [a, b] = cell.split(axis);
            stack.push((id, b));
            stack.push((id, a));
        }
        if deferred.is_empty() {
            break;
        }
        let bins = builder.bins();
        for (id, cell, cover) in deferred.drain(..) {
            if builder.covered(&cover, &bins) {
                continue;
            }
            let patch = &patches[id];
            let span = cell.span();
            let axis = usize::from(
                &patch.chart.lipschitz[1] * &span[1] > &patch.chart.lipschitz[0] * &span[0],
            );
            let [a, b] = cell.split(axis);
            stack.push((id, b));
            stack.push((id, a));
        }
        if stack.is_empty() {
            break;
        }
    }
    builder.mesh.max_deviation = options.tolerance.clone();
    Ok(builder.mesh)
}
impl Builder<'_> {
    fn classify(
        &mut self,
        patch: &Patch,
        cover: &[Ball; 3],
        point: &RPoint3,
        dec: &mut Decider,
    ) -> Result<(Location, bool), Refusal> {
        self.locations.clear();
        let intervals = cover.each_ref().map(Ball::to_interval);
        for instruction in &self.plan.instructions {
            let value = match *instruction {
                Instruction::Empty => (Location::Outside, true),
                Instruction::Leaf(i) => {
                    let mut stable = true;
                    let mut outside = false;
                    for (j, c) in self.constraints[i].iter().enumerate() {
                        if patch.coincident[i][j] {
                            continue;
                        }
                        let range = c.range(&intervals);
                        stable &= !range.contains_zero();
                        if matches!(self.plan.leaves[i], crate::geometry::Leaf::Primitive(_))
                            && range.lo() > 0.0
                        {
                            outside = true;
                            break;
                        }
                    }
                    if outside {
                        (Location::Outside, true)
                    } else if stable {
                        (
                            self.plan.leaves[i].location(
                                &[point.x.clone(), point.y.clone(), point.z.clone()],
                                dec,
                            )?,
                            true,
                        )
                    } else {
                        (Location::Indeterminate, false)
                    }
                }
                Instruction::Boolean { kind, a, b } => {
                    let (a, sa) = self.locations[a];
                    let (b, sb) = self.locations[b];
                    let location = crate::combine(kind, a, b);
                    (
                        location,
                        (sa && sb) || matches!(location, Location::Inside | Location::Outside),
                    )
                }
            };
            self.locations.push(value);
        }
        Ok(self.locations[self.plan.root])
    }
    fn emit(
        &mut self,
        patch: &Patch,
        cell: &Cell,
        face: usize,
        error: Rational,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let mut indices = [0; 4];
        let mut rounding = q(0);
        for (i, uv) in cell.corners().iter().enumerate() {
            let p = patch.chart.point(uv);
            let exact = [p.x, p.y, p.z];
            let position = exact.each_ref().map(|x| {
                let f = x.to_f64_nearest();
                if f == 0.0 { 0.0 } else { f }
            });
            if position.iter().any(|x| !x.is_finite()) {
                return Err(failure(
                    "field.display-precision",
                    "non-finite display coordinate",
                ));
            }
            let residual = (0..3).fold(q(0), |sum, i| {
                sum + (&Rational::from_f64_exact(position[i]).unwrap() - &exact[i]).abs()
            });
            rounding = rounding.max(residual);
            let key = position.map(f64::to_bits);
            indices[i] = if let Some(&id) = self.vertices.get(&key) {
                id
            } else {
                if self.mesh.positions.len() >= self.options.max_vertices {
                    return Err(failure(
                        "field.display-budget",
                        "display vertex limit reached",
                    ));
                }
                let id = self.mesh.positions.len();
                self.mesh.positions.push(position);
                self.vertices.insert(key, id);
                id
            };
        }
        if error + rounding > self.options.tolerance {
            return Err(failure(
                "field.display-precision",
                "coordinate rounding exceeds the certified display tolerance",
            ));
        }
        let point = patch.chart.point(&cell.centre());
        let direction = self
            .plan
            .face_query_checked(FaceId(face), &point, &mut self.workspace, dec)?
            .direction;
        let [a, b, c] = [indices[0], indices[1], indices[2]].map(|i| self.mesh.positions[i]);
        let u = std::array::from_fn::<_, 3, _>(|i| b[i] - a[i]);
        let v = std::array::from_fn::<_, 3, _>(|i| c[i] - a[i]);
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let reverse = if let crate::Direction::Unique(d) = direction {
            (0..3).map(|i| n[i] * (d[i].lo() + d[i].hi())).sum::<f64>() < 0.0
        } else {
            false
        };
        for mut vertices in [
            [indices[0], indices[1], indices[2]],
            [indices[0], indices[2], indices[3]],
        ] {
            if reverse {
                vertices.swap(1, 2);
            }
            self.mesh.triangles.push(Triangle {
                vertices,
                face: FaceId(face),
            });
        }
        Ok(())
    }
    fn bin(&self, p: [f64; 3]) -> [i64; 3] {
        let scale = self.options.tolerance.to_f64_nearest();
        p.map(|x| (x / scale).floor() as i64)
    }
    fn bins(&self) -> BTreeMap<[i64; 3], Vec<usize>> {
        let mut bins = BTreeMap::new();
        for (i, &p) in self.mesh.positions.iter().enumerate() {
            bins.entry(self.bin(p)).or_insert_with(Vec::new).push(i);
        }
        bins
    }
    fn covered(&self, cover: &[Ball; 3], bins: &BTreeMap<[i64; 3], Vec<usize>>) -> bool {
        let key = self.bin(cover.each_ref().map(|b| b.mid().to_f64_nearest()));
        for x in -2..=2 {
            for y in -2..=2 {
                for z in -2..=2 {
                    let Some(ids) = bins.get(&[
                        key[0].saturating_add(x),
                        key[1].saturating_add(y),
                        key[2].saturating_add(z),
                    ]) else {
                        continue;
                    };
                    for &i in ids {
                        let p = self.mesh.positions[i];
                        let distance = (0..3).fold(q(0), |sum, j| {
                            let x = Rational::from_f64_exact(p[j]).unwrap();
                            sum + (&x - &cover[j].lo()).abs().max((&x - &cover[j].hi()).abs())
                        });
                        if distance <= self.options.tolerance {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }
}
