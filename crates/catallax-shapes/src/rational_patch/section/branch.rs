//! Complete regular zero sets of a cut equation in patch parameters.

use super::*;

mod tube;
pub use tube::{SectionBranchApproximation, SectionCandidate};

/// Counts for the candidate-independent complete cell search.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SectionSearch {
    /// Cells certified to contain no zero.
    pub empty_cells: usize,
    /// Cells certified to contain exactly one regular arc.
    pub arc_cells: usize,
    /// Deepest binary split, counted from each knot span.
    pub depth: u32,
    /// Deterministic decision work, including exact boundary roots and joins.
    pub work: u64,
}

#[derive(Clone, Debug)]
struct ArcCell {
    span: usize,
    domain: Domain,
    solved: usize,
    ends: [[Real; 2]; 2],
}

/// One exact connected regular component, selected by its ordered cell cover.
/// The fitted candidate is never used as the geometric definition.
#[derive(Clone, Debug)]
pub struct PatchSectionBranch {
    patch: RationalPatch,
    surface: SectionSurface,
    cells: Arc<[ArcCell]>,
    closed: bool,
}

impl PatchSectionBranch {
    pub(crate) fn same_definition(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.cells, &other.cells)
    }
    /// Exact carrier retained beside the branch selector.
    pub fn patch(&self) -> &RationalPatch {
        &self.patch
    }
    /// Exact cutting surface.
    pub fn surface(&self) -> &SectionSurface {
        &self.surface
    }
    /// Whether the exact endpoint graph is a cycle.
    pub fn is_closed(&self) -> bool {
        self.closed
    }
    /// Endpoints of an open component, or its common seam point when closed.
    pub fn endpoints(&self) -> [Point3; 2] {
        [
            self.patch.data.spans[self.cells[0].span].point(&self.cells[0].ends[0]),
            self.patch.data.spans[self.cells.last().unwrap().span]
                .point(&self.cells.last().unwrap().ends[1]),
        ]
    }
    /// Parameter endpoints of the retained component.
    pub fn parameter_endpoints(&self) -> [[Real; 2]; 2] {
        [
            self.cells[0].ends[0].clone(),
            self.cells.last().unwrap().ends[1].clone(),
        ]
    }
    /// Exact cell-join witnesses in traversal order, including the first end.
    pub fn witnesses(&self) -> impl Iterator<Item = Point3> + '_ {
        self.cells
            .iter()
            .map(|cell| self.patch.data.spans[cell.span].point(&cell.ends[0]))
    }
    /// Reverse the natural traversal without changing the point set.
    pub fn reversed(&self) -> Self {
        let cells = self
            .cells
            .iter()
            .rev()
            .cloned()
            .map(|mut c| {
                c.ends.swap(0, 1);
                c
            })
            .collect::<Vec<_>>();
        Self {
            cells: cells.into(),
            ..self.clone()
        }
    }
    /// Whether the retained traversal follows surface-normal cross
    /// patch-normal. The cutting surface's normal is its equation gradient.
    pub fn follows_surface_first(&self, dec: &mut Decider) -> Result<bool, Refusal> {
        let cell = &self.cells[0];
        let equation = self.surface.equation(&self.patch.data.spans[cell.span]);
        let derivative = equation.partial(cell.solved).eval_real(&cell.ends[0]);
        let delta = &cell.ends[1][1 - cell.solved] - &cell.ends[0][1 - cell.solved];
        let sense = exact_sign(&(&derivative * &delta), dec)?;
        if sense == Sign::Zero {
            return Err(failure(
                "a regular branch needs a nonzero pair-normal orientation",
            ));
        }
        Ok((sense == Sign::Positive) == (cell.solved == 1))
    }
    fn locate(&self, point: &Point3, dec: &mut Decider) -> Result<Option<(usize, Real)>, Refusal> {
        for (i, cell) in self.cells.iter().enumerate() {
            let span = &self.patch.data.spans[cell.span];
            for end in 0..2 {
                if exact_sign(
                    &span.point(&cell.ends[end]).minus(point).norm_squared(),
                    dec,
                )? == Sign::Zero
                {
                    return Ok(Some((i, Real::from_i64(end as i64))));
                }
            }
        }
        if exact_sign(&self.surface.value_at(point), dec)? != Sign::Zero {
            return Ok(None);
        }
        let Some(at) = self.patch.parameters_at(point, dec)? else {
            return Ok(None);
        };
        for (i, cell) in self.cells.iter().enumerate() {
            if contains(&cell.domain, &at, dec)? {
                let axis = 1 - cell.solved;
                let t =
                    (&at[axis] - &cell.ends[0][axis]) / (&cell.ends[1][axis] - &cell.ends[0][axis]);
                return Ok(Some((i, t)));
            }
        }
        Ok(None)
    }
    /// Exact membership through the finite patch parameters and component cells.
    pub fn contains(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        self.locate(point, dec).map(|v| v.is_some())
    }
    /// Traversal order on an open component. Equal points compare equal.
    pub fn compare(&self, a: &Point3, b: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        if self.closed {
            return Err(failure("linear order requires an open patch branch"));
        }
        let (i, a) = self
            .locate(a, dec)?
            .ok_or_else(|| failure("order point is outside the patch branch"))?;
        let (j, b) = self
            .locate(b, dec)?
            .ok_or_else(|| failure("order point is outside the patch branch"))?;
        if i < j {
            Ok(Sign::Positive)
        } else if i > j {
            Ok(Sign::Negative)
        } else {
            exact_sign(&(b - a), dec)
        }
    }
    /// Sign selecting cyclic order for three points on a closed component.
    /// Positive means a,b,c follow the retained traversal. Zero means coincidence.
    pub fn cyclic_sign(&self, points: &[Point3; 3], dec: &mut Decider) -> Result<Sign, Refusal> {
        if !self.closed {
            return Err(failure("cyclic order requires a closed patch branch"));
        }
        let mut positions = Vec::new();
        for point in points {
            let (i, t) = self
                .locate(point, dec)?
                .ok_or_else(|| failure("cyclic-order point is outside the patch branch"))?;
            positions.push(Real::from_i64(i as i64) + t);
        }
        exact_sign(
            &(&(&positions[1] - &positions[0])
                * &(&(&positions[2] - &positions[0]) * &(&positions[2] - &positions[1]))),
            dec,
        )
    }
    /// Exact arriving and leaving tangents, retaining both sides of knot joins.
    pub fn tangents(&self, point: &Point3, dec: &mut Decider) -> Result<[Vector3; 2], Refusal> {
        let mut arriving = None;
        let mut leaving = None;
        for cell in self.cells.iter() {
            let span = &self.patch.data.spans[cell.span];
            for end in 0..2 {
                if exact_sign(
                    &span.point(&cell.ends[end]).minus(point).norm_squared(),
                    dec,
                )? == Sign::Zero
                {
                    let tangent = self.cell_tangent(cell, &cell.ends[end]);
                    if end == 0 {
                        leaving = Some(tangent);
                    } else {
                        arriving = Some(tangent);
                    }
                }
            }
        }
        match (arriving, leaving) {
            (Some(a), Some(b)) => Ok([a, b]),
            (Some(t), None) | (None, Some(t)) => Ok([t.clone(), t]),
            (None, None) => {
                let (i, _) = self
                    .locate(point, dec)?
                    .ok_or_else(|| failure("tangent point is outside the patch branch"))?;
                let at = self
                    .patch
                    .parameters_at(point, dec)?
                    .expect("membership checked");
                let tangent = self.cell_tangent(&self.cells[i], &at);
                Ok([tangent.clone(), tangent])
            }
        }
    }
    fn cell_tangent(&self, cell: &ArcCell, at: &[Real; 2]) -> Vector3 {
        let span = &self.patch.data.spans[cell.span];
        let equation = self.surface.equation(span);
        let fu = equation.partial(0).eval_real(at);
        let fv = equation.partial(1).eval_real(at);
        let free = 1 - cell.solved;
        let along = &cell.ends[1][free] - &cell.ends[0][free];
        let velocity = if cell.solved == 1 {
            [Real::one(), -(&fu / &fv)]
        } else {
            [-(&fv / &fu), Real::one()]
        };
        let [du, dv] = span.derivatives(at);
        du.scaled(&(&velocity[0] * &along))
            .plus(&dv.scaled(&(&velocity[1] * &along)))
    }
}

#[derive(Clone)]
struct Net {
    count: [usize; 2],
    values: Vec<Real>,
    domain: Domain,
    depth: u32,
}
impl Net {
    fn edges(&self) -> [Vec<Real>; 4] {
        let [nu, nv] = self.count;
        [
            self.values[..nv].to_vec(),
            self.values[(nu - 1) * nv..].to_vec(),
            (0..nu).map(|i| self.values[i * nv].clone()).collect(),
            (0..nu)
                .map(|i| self.values[i * nv + nv - 1].clone())
                .collect(),
        ]
    }
    fn derivative(&self, axis: usize, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let [nu, nv] = self.count;
        let mut result = Vec::new();
        for i in 0..nu - usize::from(axis == 0) {
            for j in 0..nv - usize::from(axis == 1) {
                let k = i * nv + j;
                let next = &self.values[k + if axis == 0 { nv } else { 1 }];
                reserve_values(next, &self.values[k], false, dec)?;
                result.push(next - &self.values[k]);
            }
        }
        Ok(result)
    }
    fn halves(&self, dec: &mut Decider) -> Result<[Self; 2], Refusal> {
        let axis = usize::from(
            (&self.domain[0][1] - &self.domain[0][0]) < (&self.domain[1][1] - &self.domain[1][0]),
        );
        let mut left = self.clone();
        let mut right = self.clone();
        let degree = self.count[axis] - 1;
        for fixed in 0..self.count[1 - axis] {
            let at = |k| {
                if axis == 0 {
                    k * self.count[1] + fixed
                } else {
                    fixed * self.count[1] + k
                }
            };
            let mut work = (0..=degree)
                .map(|k| self.values[at(k)].clone())
                .collect::<Vec<_>>();
            for level in 1..=degree {
                for k in 0..=degree - level {
                    reserve_values(&work[k], &work[k + 1], false, dec)?;
                    let sum = &work[k] + &work[k + 1];
                    reserve_values(&sum, &Real::from_i64(2), true, dec)?;
                    work[k] = sum / Real::from_i64(2);
                }
                left.values[at(level)] = work[0].clone();
                right.values[at(degree - level)] = work[degree - level].clone();
            }
        }
        let mid = (&self.domain[axis][0] + &self.domain[axis][1]).mul_pow2(-1);
        left.domain[axis][1] = mid.clone();
        right.domain[axis][0] = mid;
        left.depth += 1;
        right.depth += 1;
        Ok([left, right])
    }
}
fn reserve_values(a: &Real, b: &Real, multiply: bool, dec: &mut Decider) -> Result<(), Refusal> {
    if let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) {
        dec.reserve_exact_arithmetic(a, b, multiply)?;
    }
    Ok(())
}
fn uniform(values: &[Real], dec: &mut Decider) -> Result<Option<Sign>, Refusal> {
    let mut sign = None;
    for value in values {
        let next = exact_sign(value, dec)?;
        if next == Sign::Zero || sign.is_some_and(|s| s != next) {
            return Ok(None);
        }
        sign = Some(next);
    }
    Ok(sign)
}
fn variation(values: &[Real], dec: &mut Decider) -> Result<Option<usize>, Refusal> {
    let mut previous = None;
    let mut changes = 0;
    for (i, value) in values.iter().enumerate() {
        let sign = exact_sign(value, dec)?;
        if sign == Sign::Zero {
            if i == 0 || i + 1 == values.len() {
                return Ok(None);
            }
            continue;
        }
        if previous.is_some_and(|s| s != sign) {
            changes += 1;
        }
        previous = Some(sign);
    }
    Ok(Some(changes))
}
fn same_uv(a: &[Real; 2], b: &[Real; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    for i in 0..2 {
        if !a[i].same_definition(&b[i]) && exact_sign(&(&a[i] - &b[i]), dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

fn edge_root(
    equation: &RealPoly,
    domain: &Domain,
    edge: usize,
    dec: &mut Decider,
) -> Result<[Real; 2], Refusal> {
    let fixed = edge / 2;
    let bound = &domain[fixed][edge % 2];
    let restricted = equation.substitute(
        fixed,
        &RealPoly::constant(2, Real::from_rational(bound.clone())),
    );
    let polynomial = axis_polynomial(&restricted, 1 - fixed, dec)?;
    let mut found = None;
    for root in polynomial.real_roots(dec)? {
        if in_range(&root, &domain[1 - fixed], dec)? {
            if found.is_some() {
                return Err(failure("a certified cell edge has multiple isolated roots"));
            }
            let mut at = [
                Real::from_rational(bound.clone()),
                Real::from_rational(bound.clone()),
            ];
            at[1 - fixed] = root;
            found = Some(at);
        }
    }
    found.ok_or_else(|| failure("a certified edge crossing has no isolated root"))
}

pub(super) fn search(
    patch: &RationalPatch,
    surface: &SectionSurface,
    dec: &mut Decider,
) -> Result<RationalPatchSection, Refusal> {
    let start = dec.stats().work;
    let mut record = SectionSearch::default();
    let mut arcs = Vec::new();
    let [du, dv] = patch.carrier().degrees();
    for (si, span) in patch.data.spans.iter().enumerate() {
        let h = patch.carrier().span_bezier_controls(span.indices, dec)?;
        let (values, count) = surface.net(&h, [du, dv], dec)?;
        let equation = surface.equation(span);
        let mut pending = vec![Net {
            count,
            values,
            domain: span.domain.clone(),
            depth: 0,
        }];
        while let Some(net) = pending.pop() {
            dec.check_deadline()?;
            record.depth = record.depth.max(net.depth);
            if uniform(&net.values, dec)?.is_some() {
                record.empty_cells += 1;
                continue;
            }
            let mut classified = false;
            for solved in [1, 0] {
                if uniform(&net.derivative(solved, dec)?, dec)?.is_none() {
                    continue;
                }
                let mut crossings = Vec::new();
                let mut unresolved = false;
                for (edge, values) in net.edges().iter().enumerate() {
                    match variation(values, dec)? {
                        Some(0) => (),
                        Some(1) => crossings.push(edge),
                        _ => {
                            unresolved = true;
                            break;
                        }
                    }
                }
                if unresolved {
                    break;
                }
                if crossings.is_empty() {
                    record.empty_cells += 1;
                    classified = true;
                    break;
                }
                if crossings.len() == 2 {
                    let a = edge_root(&equation, &net.domain, crossings[0], dec)?;
                    let b = edge_root(&equation, &net.domain, crossings[1], dec)?;
                    arcs.push(ArcCell {
                        span: si,
                        domain: net.domain.clone(),
                        solved,
                        ends: [a, b],
                    });
                    record.arc_cells += 1;
                    classified = true;
                    break;
                }
            }
            if classified {
                continue;
            }
            let at = midpoint(&net.domain).map(Real::from_rational);
            if exact_sign(&equation.eval_real(&at), dec)? == Sign::Zero
                && exact_sign(&equation.partial(0).eval_real(&at), dec)? == Sign::Zero
                && exact_sign(&equation.partial(1).eval_real(&at), dec)? == Sign::Zero
            {
                let point = span.point(&at);
                let display = |v: &Real| {
                    v.rational_ref()
                        .map(ToString::to_string)
                        .unwrap_or_else(|| format!("{v:?}"))
                };
                return Err(Refusal::new(
                    surface.codes().tangent,
                    format!(
                        "nonregular contact at ({}, {}, {}), parameters ({}, {})",
                        display(&point.x),
                        display(&point.y),
                        display(&point.z),
                        display(&at[0]),
                        display(&at[1])
                    ),
                ));
            }
            if net.depth >= MAX_PROOF_DEPTH
                || record.empty_cells + record.arc_cells + pending.len() + 2 > MAX_PROOF_CELLS
            {
                return Err(failure(
                    "the complete section cell search reached its unchanged proof allowance",
                ));
            }
            let [left, right] = net.halves(dec)?;
            pending.push(right);
            pending.push(left);
        }
    }
    let mut vertices: Vec<[Real; 2]> = Vec::new();
    let mut links: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut endpoints = Vec::new();
    for (i, arc) in arcs.iter().enumerate() {
        let mut ids = [0, 0];
        for (end, at) in arc.ends.iter().enumerate() {
            let mut existing = None;
            for (j, old) in vertices.iter().enumerate() {
                if same_uv(old, at, dec)? {
                    existing = Some(j);
                    break;
                }
            }
            let j = existing.unwrap_or_else(|| {
                vertices.push(at.clone());
                links.push(Vec::new());
                vertices.len() - 1
            });
            links[j].push((i, end));
            ids[end] = j;
        }
        endpoints.push(ids);
    }
    for (at, adjacent) in vertices.iter().zip(&links) {
        if adjacent.len() == 2 {
            continue;
        }
        let border = (0..2).any(|axis| {
            patch.data.domain[axis]
                .iter()
                .any(|b| at[axis].rational_ref() == Some(b))
        });
        if adjacent.len() != 1 || !border {
            return Err(failure(
                "section endpoints do not form a complete manifold cover with exact patch-border ends",
            ));
        }
    }
    let mut used = vec![false; arcs.len()];
    let mut branches = Vec::new();
    while let Some(seed) = used.iter().position(|v| !v) {
        let mut start = endpoints[seed][0];
        // Find an open end before traversing, otherwise retain the first seam.
        let mut stack = vec![start];
        let mut visited = vec![false; vertices.len()];
        while let Some(v) = stack.pop() {
            if visited[v] {
                continue;
            }
            visited[v] = true;
            if links[v].len() == 1 {
                start = v;
                break;
            }
            for &(i, end) in &links[v] {
                stack.push(endpoints[i][1 - end]);
            }
        }
        let mut here = start;
        let mut ordered = Vec::new();
        while let Some(&(i, end)) = links[here].iter().find(|(i, _)| !used[*i]) {
            used[i] = true;
            let mut arc = arcs[i].clone();
            if end == 1 {
                arc.ends.swap(0, 1);
            }
            ordered.push(arc);
            here = endpoints[i][1 - end];
        }
        branches.push(PatchSectionBranch {
            patch: patch.clone(),
            surface: surface.clone(),
            cells: ordered.into(),
            closed: here == start,
        });
    }
    record.work = dec.stats().work - start;
    Ok(RationalPatchSection {
        patch: patch.clone(),
        surface: surface.clone(),
        curves: Vec::new(),
        points: Vec::new(),
        branches,
        search: Some(record),
    })
}
