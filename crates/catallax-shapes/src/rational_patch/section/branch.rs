//! Complete regular zero sets of a cut equation in patch parameters.

use super::*;

mod hits;
pub use hits::SectionHits;
mod germ;
mod tube;
pub use tube::{SectionBranchApproximation, SectionCandidate};
mod trace;
pub use trace::SectionTracePiece;
mod contact;

/// Counts for the candidate-independent complete cell search.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SectionSearch {
    /// Cells certified to contain no zero.
    pub empty_cells: usize,
    /// Cells certified to contain exactly one regular arc.
    pub arc_cells: usize,
    /// Cells whose only zeros are retained isolated corner contacts.
    pub contact_cells: usize,
    /// Deepest binary split, counted from each knot span.
    pub depth: u32,
    /// Deterministic decision work, including exact boundary roots and joins.
    pub work: u64,
}

#[derive(Clone, Debug)]
struct ArcCell {
    chart: usize,
    span: usize,
    domain: Domain,
    solved: usize,
    ends: [[Real; 2]; 2],
}

/// One exact connected regular component, selected by its ordered cell cover.
/// The fitted candidate is never used as the geometric definition.
#[derive(Clone, Debug)]
pub struct PatchSectionBranch {
    charts: Arc<[RationalPatch]>,
    surface: SectionSurface,
    cells: Arc<[ArcCell]>,
    closed: bool,
}

impl PatchSectionBranch {
    pub(crate) fn same_definition(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.cells, &other.cells)
    }
    /// Exact carrier charts retained beside the branch selector. A single
    /// patch is one chart; a joined atlas keeps every chart of its face, in
    /// the source's own parameters.
    pub fn charts(&self) -> &[RationalPatch] {
        &self.charts
    }
    fn span(&self, cell: &ArcCell) -> &Span {
        &self.charts[cell.chart].data.spans[cell.span]
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
        let last = self.cells.last().unwrap();
        [
            self.span(&self.cells[0]).point(&self.cells[0].ends[0]),
            self.span(last).point(&last.ends[1]),
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
            .map(|cell| self.span(cell).point(&cell.ends[0]))
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
        let equation = self.surface.equation(self.span(cell));
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
    /// The cell holding a point, its traversal fraction there and, away from
    /// exact cell ends, its chart parameters.
    #[allow(clippy::type_complexity)]
    fn locate(
        &self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<(usize, Real, Option<[Real; 2]>)>, Refusal> {
        for (i, cell) in self.cells.iter().enumerate() {
            let span = self.span(cell);
            for (end, parameters) in cell.ends.iter().enumerate() {
                if exact_sign(&span.point(parameters).minus(point).norm_squared(), dec)?
                    == Sign::Zero
                {
                    return Ok(Some((i, Real::from_i64(end as i64), None)));
                }
            }
        }
        if exact_sign(&self.surface.value_at(point), dec)? != Sign::Zero {
            return Ok(None);
        }
        for (chart, patch) in self.charts.iter().enumerate() {
            if self.cells.iter().all(|cell| cell.chart != chart) {
                continue;
            }
            let Some(at) = patch.parameters_at(point, dec)? else {
                continue;
            };
            for (i, cell) in self.cells.iter().enumerate() {
                if cell.chart == chart && contains(&cell.domain, &at, dec)? {
                    let axis = 1 - cell.solved;
                    let t = (&at[axis] - &cell.ends[0][axis])
                        / (&cell.ends[1][axis] - &cell.ends[0][axis]);
                    return Ok(Some((i, t, Some(at))));
                }
            }
        }
        Ok(None)
    }
    /// Whether two branches have the same point set. Components of one
    /// complete section are maximal and pairwise disjoint, and the search is
    /// exact and deterministic, so with equal carrier definitions one shared
    /// cell end decides. On the same charts with
    /// different cutting surfaces, a cell end where the other cut polynomial
    /// is nonzero proves the sets differ. `None` when neither proof applies.
    pub fn same_component(&self, other: &Self, dec: &mut Decider) -> Result<Option<bool>, Refusal> {
        let same_charts = self.charts.len() == other.charts.len()
            && self
                .charts
                .iter()
                .zip(other.charts.iter())
                .all(|(a, b)| a.carrier().same_definition(b.carrier()));
        if !same_charts {
            return Ok(None);
        }
        if self
            .surface
            .surface()
            .same_definition(&other.surface.surface())
        {
            // The search is exact and deterministic, so equal definitions give
            // equal cells: one component has the same cell ends in both.
            let start = &self.cells[0].ends[0];
            for cell in other.cells.iter() {
                for end in &cell.ends {
                    if same_uv(start, end, dec)? {
                        return Ok(Some(true));
                    }
                }
            }
            return Ok(Some(false));
        }
        for (a, b) in [(self, other), (other, self)] {
            for cell in a.cells.iter() {
                let equation = b.surface.equation(a.span(cell));
                for end in &cell.ends {
                    if exact_sign(&equation.eval_real(end), dec)? != Sign::Zero {
                        return Ok(Some(false));
                    }
                }
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
        let (i, a, _) = self
            .locate(a, dec)?
            .ok_or_else(|| failure("order point is outside the patch branch"))?;
        let (j, b, _) = self
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
            let (i, t, _) = self
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
            let span = self.span(cell);
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
                let (i, _, at) = self
                    .locate(point, dec)?
                    .ok_or_else(|| failure("tangent point is outside the patch branch"))?;
                let at = at.expect("cell ends were matched above");
                let tangent = self.cell_tangent(&self.cells[i], &at);
                Ok([tangent.clone(), tangent])
            }
        }
    }
    fn cell_tangent(&self, cell: &ArcCell, at: &[Real; 2]) -> Vector3 {
        let span = self.span(cell);
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
    for value in values {
        let sign = exact_sign(value, dec)?;
        if sign == Sign::Zero {
            continue;
        }
        if previous.is_some_and(|s| s != sign) {
            changes += 1;
        }
        previous = Some(sign);
    }
    Ok(previous.map(|_| changes))
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

/// Charts searched together, in one shared parameter space.
pub(crate) struct Cover {
    pub(crate) charts: Arc<[RationalPatch]>,
    /// Whether each chart side, ordered `u` low, `u` high, `v` low, `v` high,
    /// is a physical border rather than a join to another chart.
    pub(crate) physical: Vec<[bool; 4]>,
    /// Exactly identified end rows of each periodic parameter.
    pub(crate) wrap: [Option<[Rational; 2]>; 2],
}

pub(super) fn search(
    patch: &RationalPatch,
    surface: &SectionSurface,
    dec: &mut Decider,
) -> Result<RationalPatchSection, Refusal> {
    let cover = Cover {
        charts: Arc::from([patch.clone()]),
        physical: vec![[true; 4]],
        wrap: [None, None],
    };
    let section = search_cover(&cover, surface, false, dec)?;
    Ok(RationalPatchSection {
        patch: patch.clone(),
        surface: surface.clone(),
        curves: Vec::new(),
        points: section.points,
        branches: section.branches,
        search: Some(section.search),
    })
}

pub(crate) struct CoverSection {
    pub branches: Vec<PatchSectionBranch>,
    pub points: Vec<Point3>,
    pub search: SectionSearch,
}

/// The complete regular section over every chart of a cover. Arc ends join
/// by exact shared parameters; an open end must lie on a physical border.
pub(crate) fn search_cover(
    cover: &Cover,
    surface: &SectionSurface,
    trace_ready: bool,
    dec: &mut Decider,
) -> Result<CoverSection, Refusal> {
    enum Crossing {
        Edge(usize),
        Corner([Real; 2]),
    }
    let start = dec.stats().work;
    let mut record = SectionSearch::default();
    let mut arcs = Vec::new();
    let mut points: Vec<Point3> = Vec::new();
    let mut corner_contacts = Vec::new();
    for (ci, patch) in cover.charts.iter().enumerate() {
        let [du, dv] = patch.carrier().degrees();
        for (si, span) in patch.data.spans.iter().enumerate() {
            let h = patch.carrier().span_bezier_controls(span.indices, dec)?;
            if !surface.retains_span(&h, dec)? {
                record.empty_cells += 1;
                continue;
            }
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
                if matches!(surface, SectionSurface::Plane(_))
                    && let Some(corners) = contact::corners(&net, dec)?
                {
                    let isolated = !corners.is_empty();
                    for at in corners {
                        let point = span.point(&at);
                        let mut duplicate = false;
                        for previous in &points {
                            if exact_sign(&point.minus(previous).norm_squared(), dec)? == Sign::Zero
                            {
                                duplicate = true;
                                break;
                            }
                        }
                        if !duplicate {
                            points.push(point);
                        }
                    }
                    if isolated {
                        record.contact_cells += 1;
                    } else {
                        record.empty_cells += 1;
                    }
                    continue;
                }
                let mut classified = false;
                for solved in [1, 0] {
                    if uniform(&net.derivative(solved, dec)?, dec)?.is_none() {
                        continue;
                    }
                    if trace_ready {
                        // Retained trace jets also need a strict power-form bound.
                        let derivative = equation.partial(solved).eval_ball(
                            &box_balls(&net.domain, dec.policy().max_precision_bits.min(128)),
                            dec,
                        )?;
                        if derivative.sign().is_none() {
                            continue;
                        }
                    }
                    let mut crossings = Vec::new();
                    let mut unresolved = false;
                    for (edge, values) in net.edges().iter().enumerate() {
                        let first = exact_sign(&values[0], dec)? == Sign::Zero;
                        let last = exact_sign(values.last().unwrap(), dec)? == Sign::Zero;
                        match variation(values, dec)? {
                            Some(0) => {
                                for (end, zero) in [first, last].into_iter().enumerate() {
                                    if zero {
                                        let mut at = net
                                            .domain
                                            .each_ref()
                                            .map(|bounds| Real::from(bounds[edge % 2].clone()));
                                        at[1 - edge / 2] =
                                            Real::from(net.domain[1 - edge / 2][end].clone());
                                        let mut duplicate = false;
                                        for old in &crossings {
                                            if let Crossing::Corner(old) = old {
                                                duplicate |= same_uv(old, &at, dec)?;
                                            }
                                        }
                                        if !duplicate {
                                            crossings.push(Crossing::Corner(at));
                                        }
                                    }
                                }
                            }
                            Some(1) if !first && !last => crossings.push(Crossing::Edge(edge)),
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
                        let point = |crossing, dec: &mut Decider| match crossing {
                            Crossing::Edge(edge) => edge_root(&equation, &net.domain, edge, dec),
                            Crossing::Corner(point) => Ok(point),
                        };
                        let a = point(crossings.remove(0), dec)?;
                        let b = point(crossings.remove(0), dec)?;
                        if exact_sign(&(&a[1 - solved] - &b[1 - solved]), dec)? == Sign::Zero {
                            continue;
                        }
                        arcs.push(ArcCell {
                            chart: ci,
                            span: si,
                            domain: net.domain.clone(),
                            solved,
                            ends: [a, b],
                        });
                        record.arc_cells += 1;
                        classified = true;
                        break;
                    }
                    if crossings.len() == 1 {
                        let mut side = None;
                        let mut one_side = true;
                        for value in &net.values {
                            let sign = exact_sign(value, dec)?;
                            if sign != Sign::Zero {
                                if side.is_some_and(|old| old != sign) {
                                    one_side = false;
                                    break;
                                }
                                side = Some(sign);
                            }
                        }
                        if one_side && let Crossing::Corner(contact) = crossings.remove(0) {
                            corner_contacts.push(contact);
                            record.empty_cells += 1;
                            classified = true;
                            break;
                        }
                    }
                }
                if classified {
                    continue;
                }
                let at = midpoint(&net.domain).map(Real::from_rational);
                if exact_sign(&equation.eval_real(&at), dec)? == Sign::Zero
                    && exact_sign(&equation.partial(0).eval_real(&at), dec)? == Sign::Zero
                    && exact_sign(&equation.partial(1).eval_real(&at), dec)? == Sign::Zero
                    && !(matches!(surface, SectionSurface::Plane(_))
                        && contact::definite_hessian(&equation, &at, dec)?)
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
                    || record.empty_cells
                        + record.arc_cells
                        + record.contact_cells
                        + pending.len()
                        + 2
                        > MAX_PROOF_CELLS
                {
                    return Err(surface.unresolved(
                        "the complete section cell search reached its unchanged proof allowance",
                    ));
                }
                let [left, right] = net.halves(dec)?;
                pending.push(right);
                pending.push(left);
            }
        }
    }
    let mut vertices: Vec<[Real; 2]> = Vec::new();
    let mut links: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut endpoints = Vec::new();
    for (i, arc) in arcs.iter().enumerate() {
        let mut ids = [0, 0];
        for (end, at) in arc.ends.iter().enumerate() {
            let mut at = at.clone();
            for (axis, wrap) in cover.wrap.iter().enumerate() {
                if let Some([lo, hi]) = wrap
                    && at[axis].rational_ref() == Some(hi)
                {
                    at[axis] = Real::from_rational(lo.clone());
                }
            }
            let at = &at;
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
    for mut contact in corner_contacts {
        for (axis, wrap) in cover.wrap.iter().enumerate() {
            if let Some([lo, hi]) = wrap
                && contact[axis].rational_ref() == Some(hi)
            {
                contact[axis] = Real::from(lo.clone());
            }
        }
        let mut joined = false;
        for vertex in &vertices {
            joined |= same_uv(&contact, vertex, dec)?;
        }
        if !joined {
            return Err(failure(
                "an isolated corner contact needs a separate contact certificate",
            ));
        }
    }
    for adjacent in &links {
        if adjacent.len() == 2 {
            continue;
        }
        let border = adjacent.len() == 1 && {
            let (arc, end) = adjacent[0];
            let chart = arcs[arc].chart;
            let domain = cover.charts[chart].domain();
            let at = &arcs[arc].ends[end];
            (0..4).any(|side| {
                cover.physical[chart][side]
                    && at[side / 2].rational_ref() == Some(&domain[side / 2][side % 2])
            })
        };
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
            charts: cover.charts.clone(),
            surface: surface.clone(),
            cells: ordered.into(),
            closed: here == start,
        });
    }
    let mut isolated = Vec::new();
    for point in points {
        let mut covered = false;
        for branch in &branches {
            if branch.contains(&point, dec)? {
                covered = true;
                break;
            }
        }
        if !covered {
            isolated.push(point);
        }
    }
    record.work = dec.stats().work - start;
    Ok(CoverSection {
        branches,
        points: isolated,
        search: record,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BSplineSurface, Vector3};

    #[test]
    fn a_regular_section_joins_exact_cell_corners() {
        let mut dec = Decider::default();
        let knots = vec![
            Rational::zero(),
            Rational::zero(),
            Rational::one(),
            Rational::one(),
        ];
        let patch = RationalPatch::new(
            BSplineSurface::new(
                [1, 1],
                [knots.clone(), knots],
                [2, 2],
                [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
                    .map(|[x, y, z]| Point3::new(x.into(), y.into(), z.into()))
                    .to_vec(),
                vec![Real::one(); 4],
                &mut dec,
            )
            .unwrap(),
            &mut dec,
        )
        .unwrap();
        let plane = |x: i64| {
            SectionSurface::Plane(Plane {
                origin: Point3::new(x.into(), Real::zero(), Real::zero()),
                normal: Vector3::new(Real::one(), Real::one(), Real::zero()),
                ref_dir: Vector3::new(Real::zero(), Real::zero(), Real::one()),
            })
        };
        let section = search(&patch, &plane(1), &mut dec).unwrap();
        let [branch] = section.branches() else {
            panic!("one open diagonal");
        };
        assert!(!branch.is_closed());
        assert!(
            branch
                .contains(
                    &Point3::new(
                        Real::from(Rational::from_ratio(1, 2).unwrap()),
                        Real::from(Rational::from_ratio(1, 2).unwrap()),
                        Real::zero(),
                    ),
                    &mut dec
                )
                .unwrap()
        );
        assert!(search(&patch, &plane(0), &mut dec).is_err());
    }
}
