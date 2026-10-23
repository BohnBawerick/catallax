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
mod diagnostics;
use chart::Chart;
use diagnostics::Snapshot;
pub use diagnostics::{Cost, Diagnostics, FaceWork, SupportWork};

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
    linear: Option<([Rational; 3], Rational)>,
    anchor: [Interval; 3],
}
impl Constraint {
    fn range(&self, p: &[Interval; 3]) -> Interval {
        let p: [Interval; 3] = std::array::from_fn(|i| p[i] - self.anchor[i]);
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

    fn interior_sign(
        &self,
        chart: &Chart,
        cell: &Cell,
        corners: &mut Option<[RPoint3; 4]>,
    ) -> Option<catallax_numbers::Sign> {
        let (normal, constant) = self.linear.as_ref()?;
        if !chart.corner_signs_suffice(normal) {
            return None;
        }
        let corners = corners.get_or_insert_with(|| cell.corners().map(|uv| chart.point(&uv)));
        let mut positive = false;
        let mut negative = false;
        for p in corners {
            let value =
                constant.clone() + &normal[0] * &p.x + &normal[1] * &p.y + &normal[2] * &p.z;
            positive |= value > q(0);
            negative |= value < q(0);
        }
        match (positive, negative) {
            (true, false) => Some(catallax_numbers::Sign::Positive),
            (false, true) => Some(catallax_numbers::Sign::Negative),
            _ => None,
        }
    }
}
struct Patch {
    chart: Chart,
    faces: Vec<usize>,
    coincident: Vec<Vec<bool>>,
    group: usize,
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
    mesh_plan_inner(plan, options, dec, None)
}

/// Builds the same mesh and retains non-charging per-face/support observations,
/// including a refusal's executed prefix. No partial mesh is returned on error.
pub fn mesh_plan_reported(
    plan: &Plan,
    options: &Options,
    dec: &mut Decider,
    report: &mut Diagnostics,
) -> Result<Mesh, Refusal> {
    *report = Diagnostics::default();
    let before = Snapshot::new(dec);
    let result = mesh_plan_inner(plan, options, dec, Some(report));
    report.total.record(before, dec, 0, false);
    report.total.cells = report.supports.iter().map(|s| s.total.cells).sum();
    report.total.max_depth = report
        .supports
        .iter()
        .map(|s| s.total.max_depth)
        .max()
        .unwrap_or(0);
    report.common.work = report
        .total
        .work
        .saturating_sub(report.supports.iter().map(|s| s.total.work).sum());
    report.common.decisions = report
        .total
        .decisions
        .saturating_sub(report.supports.iter().map(|s| s.total.decisions).sum());
    result
}

fn mesh_plan_inner(
    plan: &Plan,
    options: &Options,
    dec: &mut Decider,
    mut report: Option<&mut Diagnostics>,
) -> Result<Mesh, Refusal> {
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
            // Translate exactly once, before interval evaluation. Expanding
            // (x-c)^2 in world coordinates needlessly correlates x^2 and 2cx.
            let anchor = c
                .reference_point()
                .and_then(|p| Some([p.x.as_rational()?, p.y.as_rational()?, p.z.as_rational()?]))
                .unwrap_or_else(|| [q(0), q(0), q(0)]);
            let mut local = None;
            for (axis, shift) in anchor.iter().enumerate() {
                if !shift.is_zero() {
                    let variable = catallax_numbers::RealPoly::var(3, axis);
                    let constant = catallax_numbers::RealPoly::constant(
                        3,
                        catallax_numbers::Real::from(shift.clone()),
                    );
                    local = Some(
                        local
                            .as_ref()
                            .unwrap_or(&c.polynomial)
                            .substitute(axis, &(&variable + &constant)),
                    );
                }
            }
            let mut terms = Vec::new();
            for (m, c) in local.as_ref().unwrap_or(&c.polynomial).terms() {
                terms.push((
                    m.exponents().try_into().expect("three coordinates"),
                    crate::geometry::enclose(&c, 128, dec)?.to_interval(),
                ));
            }
            let linear = c.polynomial.rational().and_then(|p| {
                let mut normal = [q(0), q(0), q(0)];
                let mut constant = q(0);
                for (m, coefficient) in p.terms() {
                    let degree: u32 = m.exponents().iter().sum();
                    if degree > 1 {
                        return None;
                    }
                    if degree == 0 {
                        constant = constant + coefficient;
                    } else {
                        let axis = m.exponents().iter().position(|&e| e == 1).unwrap();
                        normal[axis] = &normal[axis] + &coefficient;
                    }
                }
                Some((normal, constant))
            });
            constraints.push(Constraint {
                terms,
                support,
                linear,
                anchor: anchor.each_ref().map(Interval::from_rational),
            });
        }
        builder.constraints.push(constraints);
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (i, face) in plan.faces.iter().enumerate() {
        if !face.in_result || plan.leaf_registers[face.leaf].is_none() {
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
    if let Some(report) = report.as_deref_mut() {
        report.faces = plan
            .faces
            .iter()
            .enumerate()
            .map(|(i, _)| FaceWork {
                face: FaceId(i),
                candidate_cells: 0,
                max_depth: 0,
                queries: Cost::default(),
            })
            .collect();
        report.supports = groups
            .iter()
            .map(|faces| SupportWork {
                faces: faces.iter().copied().map(FaceId).collect(),
                total: Cost::default(),
                classification: Cost::default(),
                queries: Cost::default(),
                emission: Cost::default(),
            })
            .collect();
    }
    let source_boxes: Vec<_> = plan
        .leaves
        .iter()
        .map(|leaf| leaf.display_box(dec))
        .collect::<Result<_, _>>()?;
    let mut patches = Vec::new();
    let mut face_work: Vec<(usize, usize, u32)> = groups.iter().map(|g| (g[0], 0, 0)).collect();
    for (group, faces) in groups.into_iter().enumerate() {
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
                group,
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
            let before = report.as_ref().map(|_| Snapshot::new(dec));
            let previous_cells = builder.mesh.cells;
            let mut deferred_cover = None;
            let result = (|| -> Result<Visit, Refusal> {
                dec.check_deadline()?;
                builder.mesh.cells += 1;
                let work = &mut face_work[patches[id].group];
                work.1 += 1;
                work.2 = work.2.max(cell.depth);
                if let Some(report) = report.as_deref_mut() {
                    for &face in &patches[id].faces {
                        report.faces[face].candidate_cells += 1;
                        report.faces[face].max_depth = report.faces[face].max_depth.max(cell.depth);
                    }
                }
                if builder.mesh.cells > options.max_cells || cell.depth > 64 {
                    return Err(failure(
                        "field.display-budget",
                        format!(
                            "refinement limit: {} cells; per support (origin, cells, maximum depth): {:?}",
                            builder.mesh.cells,
                            face_work
                                .iter()
                                .map(|(face, cells, depth)| (
                                    &plan.faces[*face].origin,
                                    *cells,
                                    *depth
                                ))
                                .collect::<Vec<_>>()
                        ),
                    ));
                }
                let patch = &patches[id];
                let point = patch.chart.point(&cell.centre());
                let cover = patch.chart.cover(&cell, &point);
                if patch.faces.iter().all(|&face| {
                    source_boxes[plan.faces[face].leaf]
                        .as_ref()
                        .is_some_and(|bounds| {
                            (0..3).any(|i| {
                                cover[i].hi() < bounds[i].lo() || cover[i].lo() > bounds[i].hi()
                            })
                        })
                }) {
                    return Ok(Visit::Skip);
                }
                let stage = report.as_ref().map(|_| Snapshot::new(dec));
                let classified = builder.classify(patch, &cell, &cover, &point, dec);
                if let (Some(report), Some(stage)) = (report.as_deref_mut(), stage) {
                    report.supports[patch.group]
                        .classification
                        .record(stage, dec, cell.depth, true);
                }
                let (location, stable) = classified?;
                if matches!(location, Location::Inside | Location::Outside) {
                    return Ok(Visit::Skip);
                }
                let diameter = cover.iter().fold(q(0), |sum, b| sum + b.hi() - b.lo());
                let small = diameter <= options.tolerance.mul_pow2(-1);
                if stable || small {
                    let mut present = None;
                    let mut absent = true;
                    for &face in &patch.faces {
                        let stage = report.as_ref().map(|_| Snapshot::new(dec));
                        let queried = builder.face_at(face, &point, stable, dec);
                        if let (Some(report), Some(stage)) = (report.as_deref_mut(), stage) {
                            report.faces[face]
                                .queries
                                .record(stage, dec, cell.depth, true);
                            report.supports[patch.group]
                                .queries
                                .record(stage, dec, cell.depth, true);
                        }
                        let query = queried?;
                        if query.state == FaceState::OnFace {
                            present = Some((face, query.direction));
                            break;
                        }
                        absent &= query.state == FaceState::NotOnFace;
                    }
                    if let Some((face, direction)) = present {
                        let chord = patch.chart.chord(&cell);
                        if (stable && chord <= options.tolerance.mul_pow2(-1)) || small {
                            let stage = report.as_ref().map(|_| Snapshot::new(dec));
                            let emitted = builder.emit(
                                patch,
                                &cell,
                                face,
                                if stable { chord } else { diameter },
                                direction,
                            );
                            if let (Some(report), Some(stage)) = (report.as_deref_mut(), stage) {
                                report.supports[patch.group]
                                    .emission
                                    .record(stage, dec, cell.depth, true);
                            }
                            emitted?;
                            return Ok(Visit::Skip);
                        }
                    } else if stable && absent {
                        // Only the open cell is classified. Closed surviving
                        // neighbours cover its aligned trim edges without area.
                        return Ok(Visit::Skip);
                    }
                    if small {
                        deferred_cover = Some(cover);
                        return Ok(Visit::Defer);
                    }
                }
                Ok(Visit::Split)
            })();
            if let (Some(report), Some(before)) = (report.as_deref_mut(), before) {
                report.supports[patches[id].group].total.record(
                    before,
                    dec,
                    cell.depth,
                    builder.mesh.cells != previous_cells,
                );
            }
            match result? {
                Visit::Skip => continue,
                Visit::Defer => {
                    deferred.push((id, cell, deferred_cover.expect("deferred cell cover")));
                    continue;
                }
                Visit::Split => {
                    let patch = &patches[id];
                    let span = cell.span();
                    let axis = usize::from(
                        &patch.chart.lipschitz[1] * &span[1] > &patch.chart.lipschitz[0] * &span[0],
                    );
                    let [a, b] = cell.split(axis);
                    stack.push((id, b));
                    stack.push((id, a));
                }
            }
        }
        if deferred.is_empty() {
            break;
        }
        let bins = builder.bins();
        for (id, cell, cover) in deferred.drain(..) {
            let before = report.as_ref().map(|_| Snapshot::new(dec));
            let covered = builder.covered(&cover, &bins);
            if let (Some(report), Some(before)) = (report.as_deref_mut(), before) {
                report.supports[patches[id].group]
                    .total
                    .record(before, dec, cell.depth, false);
            }
            if covered {
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
enum Visit {
    Skip,
    Defer,
    Split,
}

impl Builder<'_> {
    fn ray_at(
        &mut self,
        point: &[Rational; 3],
        direction: &catallax_shapes::Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        self.workspace.rays.clear();
        for (slot, instruction) in self.plan.instructions.iter().enumerate() {
            let location = self.locations[slot].0;
            let location = if matches!(location, Location::Inside | Location::Outside) {
                // Strict regional classification is locally constant at the
                // cell's interior witness, including under an infinitesimal ray.
                location
            } else {
                match *instruction {
                    Instruction::Empty => Location::Outside,
                    Instruction::Leaf(i) => self.plan.leaves[i].ray(point, direction, dec)?,
                    Instruction::Boolean { kind, a, b } => {
                        crate::combine(kind, self.workspace.rays[a], self.workspace.rays[b])
                    }
                }
            };
            self.workspace.rays.push(location);
        }
        Ok(self.workspace.rays[self.plan.root])
    }

    fn face_at(
        &mut self,
        id: usize,
        point: &RPoint3,
        stable: bool,
        dec: &mut Decider,
    ) -> Result<crate::FaceQuery, Refusal> {
        use crate::{Direction, FaceQuery};
        let absent = FaceQuery {
            state: FaceState::NotOnFace,
            direction: Direction::Undefined,
        };
        let unknown = FaceQuery {
            state: FaceState::Indeterminate,
            direction: Direction::Ambiguous,
        };
        let plan = self.plan;
        let face = &plan.faces[id];
        let p = [point.x.clone(), point.y.clone(), point.z.clone()];
        let slot = plan.leaf_registers[face.leaf].expect("live display support");
        let location = match self.locations[slot].0 {
            Location::Indeterminate => plan.leaves[face.leaf].location(&p, dec)?,
            location => location,
        };
        match location {
            Location::Inside | Location::Outside => return Ok(absent),
            Location::Indeterminate => return Ok(unknown),
            Location::Boundary => (),
        }
        match plan.leaves[face.leaf].on_face(plan.definition_role(id), &p, dec)? {
            Some(false) => return Ok(absent),
            None => return Ok(unknown),
            Some(true) => (),
        }
        let n = crate::normal(face, &p);
        if crate::geometry::sign(&n.norm_squared(), dec)? == catallax_numbers::Sign::Zero {
            return Ok(unknown);
        }
        let minus_n = n.scaled(&catallax_numbers::Real::from_i64(-1));
        let plus = self.ray_at(&p, &n, dec)?;
        let minus = self.ray_at(&p, &minus_n, dec)?;
        let n = match (plus, minus) {
            (Location::Outside, Location::Inside) => n,
            (Location::Inside, Location::Outside) => minus_n,
            (Location::Inside, Location::Inside) | (Location::Outside, Location::Outside)
                if stable =>
            {
                return Ok(absent);
            }
            _ => return Ok(unknown),
        };
        Ok(FaceQuery {
            state: FaceState::OnFace,
            direction: crate::unit_direction(n, dec)?,
        })
    }

    fn classify(
        &mut self,
        patch: &Patch,
        cell: &Cell,
        cover: &[Ball; 3],
        point: &RPoint3,
        dec: &mut Decider,
    ) -> Result<(Location, bool), Refusal> {
        self.locations.clear();
        let intervals = cover.each_ref().map(Ball::to_interval);
        let mut corners = None;
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
                        let sign = range
                            .sign()
                            .map(|s| match s {
                                catallax_numbers::StrictSign::Positive => {
                                    catallax_numbers::Sign::Positive
                                }
                                catallax_numbers::StrictSign::Negative => {
                                    catallax_numbers::Sign::Negative
                                }
                            })
                            .or_else(|| c.interior_sign(&patch.chart, cell, &mut corners));
                        stable &= sign.is_some();
                        if matches!(self.plan.leaves[i], crate::geometry::Leaf::Primitive(_))
                            && sign == Some(catallax_numbers::Sign::Positive)
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
        direction: crate::Direction,
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
