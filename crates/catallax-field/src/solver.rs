//! F7 no-mesh solver adapter. Certified geometry/integration enclosures are not
//! a certified elasticity solve. See `docs/SOLVER-QUERIES.md` for coverage.
use crate::{
    Instruction, Location, Plan, Query, Workspace, distance, geometry,
    printing::{self, BoxClassifier, BoxState, BoxStorage},
};
use catallax_measure::Enclosure;
use catallax_numbers::{
    Decider, OperationAllowance, Rational, Real, Refusal, operation_budget_active,
};
use catallax_recipe::{Origin, Recipe};
use catallax_shapes::RPoint3;

/// Exact axis intervals in millimetres.
pub type Bounds = [[Rational; 2]; 3];
fn q(n: i64) -> Rational {
    Rational::from(n)
}
fn refusal(code: &'static str, text: &str) -> Refusal {
    Refusal::new(code, text)
}
fn volume(b: &Bounds) -> Rational {
    (0..3).fold(q(1), |v, i| v * (&b[i][1] - &b[i][0]))
}
fn valid(b: &Bounds) -> bool {
    b.iter().all(|[a, b]| a < b)
}
fn intersect(a: &Bounds, b: &Bounds) -> Option<Bounds> {
    let c = std::array::from_fn(|i| {
        [
            a[i][0].clone().max(b[i][0].clone()),
            a[i][1].clone().min(b[i][1].clone()),
        ]
    });
    valid(&c).then_some(c)
}

/// Explicit homogeneous material or a complete ordered slab partition of the
/// compact bounding box. Material IDs carry no implicit constitutive law.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    /// Stable region identity.
    pub id: String,
    /// Explicit material identity.
    pub material: String,
    /// Lower slab coordinate in millimetres.
    pub lo: Rational,
    /// Upper slab coordinate in millimetres.
    pub hi: Rational,
}
/// Explicit material assignment; gaps, overlaps and duplicate region names refuse.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Materials {
    /// One material throughout the solid.
    Homogeneous {
        /// Stable region name.
        region: String,
        /// Explicit material identity, not a constitutive model.
        material: String,
    },
    /// Partition by planes normal to one world axis.
    Slabs {
        /// World axis: 0=x, 1=y, 2=z.
        axis: usize,
        /// Ordered contiguous intervals covering the compact bounds.
        regions: Vec<Region>,
    },
}
/// Certified material membership, preserving interfaces and uncertainty.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MaterialAt {
    /// Certified outside the part.
    Void,
    /// Interior to one region.
    Region {
        /// Stable region name.
        region: String,
        /// Explicit material identity.
        material: String,
    },
    /// Internal interface, with both adjacent (region, material) identities.
    Interface(Vec<(String, String)>),
    /// Part boundary; never silently assigned an interior material.
    Boundary,
    /// Point membership or assignment did not settle.
    Unresolved,
}

/// Owns the immutable compiled recipe; its revision cannot be supplied separately
/// from the plan. This first cut accepts native recipes, not external file maps.
pub struct Definition {
    plan: Plan,
    revision: String,
    materials: Materials,
}
impl Definition {
    /// Compile native geometry and validate its explicit material partition.
    pub fn compile(
        recipe: &Recipe,
        materials: Materials,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let plan = Plan::compile_checked(recipe, dec)?;
        let b = printing::bounds(&plan, dec)?;
        match &materials {
            Materials::Homogeneous { region, material }
                if !region.is_empty() && !material.is_empty() => {}
            Materials::Slabs { axis, regions } if *axis < 3 && !regions.is_empty() => {
                let mut end = b[*axis][0].clone();
                let mut ids = std::collections::BTreeSet::new();
                for r in regions {
                    dec.charge_work(1)?;
                    if r.id.is_empty()
                        || r.material.is_empty()
                        || !ids.insert(&r.id)
                        || r.lo != end
                        || r.hi <= r.lo
                    {
                        return Err(refusal(
                            "solver.material-partition",
                            "regions must be ordered, disjoint, named and gap-free",
                        ));
                    }
                    end = r.hi.clone();
                }
                if end != b[*axis][1] {
                    return Err(refusal(
                        "solver.material-partition",
                        "slabs must cover the complete compact bounds",
                    ));
                }
            }
            _ => {
                return Err(refusal(
                    "solver.material-partition",
                    "invalid or unsupported material assignment",
                ));
            }
        }
        Ok(Self {
            plan,
            revision: identity(recipe, &materials),
            materials,
        })
    }
    /// Identity of this geometry and material assignment.
    pub fn revision(&self) -> &str {
        &self.revision
    }
    /// Prepare F6 coverage and the native-box shortcut without constructing a B-rep.
    pub fn prepare(&self, dec: &mut Decider) -> Result<Queries<'_>, Refusal> {
        let bounds = printing::bounds(&self.plan, dec)?;
        let exact_box = match self.plan.instructions[self.plan.root] {
            Instruction::Leaf(id) => match &self.plan.leaves[id] {
                geometry::Leaf::Primitive(p) => crate::mass::axis_block(p),
                _ => None,
            },
            _ => None,
        };
        Ok(Queries {
            definition: self,
            classifier: BoxClassifier::new(&self.plan, dec)?,
            bounds,
            exact_box,
        })
    }
}

/// Prepared shared geometry. Workspaces and Deciders remain local to each caller.
pub struct Queries<'a> {
    definition: &'a Definition,
    classifier: BoxClassifier<'a>,
    bounds: Bounds,
    exact_box: Option<Bounds>,
}
/// Integration pieces are closed boxes whose overlaps have zero volume. Each
/// piece is either proven material or unresolved, never an estimated occupancy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Piece {
    /// Exact closed axis-aligned bounds in millimetres.
    pub bounds: Bounds,
    /// True only for a certified material piece; false retains unresolved measure.
    pub certain: bool,
    /// Stable region identity.
    pub region: String,
    /// Explicit material identity.
    pub material: String,
}
/// Attributable finite-cell integration data; no stiffness or solver certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct Cell {
    /// SHA-256 of canonical recipe and material assignment.
    pub revision: String,
    /// Exact closed axis-aligned bounds in millimetres.
    pub bounds: Bounds,
    /// Certified closed-cell state; Unknown is not a proved cut.
    pub state: BoxState,
    /// Enclosed occupied fraction in [0,1]; width is unresolved measure divided by cell volume.
    pub fraction: Enclosure,
    /// Disjoint-interior integration boxes, partitioned by material region.
    pub pieces: Vec<Piece>,
    /// True only for a nonempty single native axis-box intersection.
    pub connected_box: bool,
}
impl Cell {
    /// All 27 moments integral xi_x^a xi_y^b xi_z^c dV, index 9a+3b+c.
    /// The polynomial is nonnegative on the cell. Unresolved pieces contribute
    /// [0, exact box integral]. Region interfaces have zero measure.
    pub fn moments(&self, region: &str, dec: &mut Decider) -> Result<[Enclosure; 27], Refusal> {
        dec.charge_work(self.pieces.len() as u64)?;
        let mut out = std::array::from_fn(|_| Enclosure { lo: q(0), hi: q(0) });
        for piece in self.pieces.iter().filter(|p| p.region == region) {
            dec.charge_work(27)?;
            let factors: [[Rational; 3]; 3] = std::array::from_fn(|i| {
                let h = &self.bounds[i][1] - &self.bounds[i][0];
                let l = (&piece.bounds[i][0] - &self.bounds[i][0])
                    .checked_div(&h)
                    .unwrap();
                let u = (&piece.bounds[i][1] - &self.bounds[i][0])
                    .checked_div(&h)
                    .unwrap();
                let mut lp = q(1);
                let mut up = q(1);
                std::array::from_fn(|a| {
                    lp = &lp * &l;
                    up = &up * &u;
                    (&h * &(&up - &lp)).checked_div(&q(a as i64 + 1)).unwrap()
                })
            });
            for a in 0..3 {
                for b in 0..3 {
                    for c in 0..3 {
                        let v = &(&factors[0][a] * &factors[1][b]) * &factors[2][c];
                        let m = &mut out[9 * a + 3 * b + c];
                        m.hi = &m.hi + &v;
                        if piece.certain {
                            m.lo = &m.lo + &v;
                        }
                    }
                }
            }
        }
        if !self.pieces.iter().any(|p| p.region == region) {
            return Err(refusal(
                "solver.region-absent",
                "no integration pieces for this region in this cell",
            ));
        }
        Ok(out)
    }
}
impl Queries<'_> {
    /// Certified compact whole-part bounds.
    pub fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    /// Q1: exact point membership and labelled field value.
    pub fn membership(
        &self,
        p: &RPoint3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Query, Refusal> {
        self.definition.plan.proven_checked(p, storage, dec)
    }
    /// Q2: classify every point of a closed rational cell.
    pub fn whole_cell(
        &self,
        b: &Bounds,
        storage: &mut BoxStorage,
        dec: &mut Decider,
    ) -> Result<BoxState, Refusal> {
        if !valid(b) {
            return Err(refusal(
                "solver.cell",
                "cell must have strictly positive rational extents",
            ));
        }
        self.classifier.classify(
            &std::array::from_fn(|i| b[i][0].clone()),
            &std::array::from_fn(|i| b[i][1].clone()),
            storage,
            dec,
        )
    }
    /// Q6: material membership from geometry and the explicit region partition.
    pub fn material(
        &self,
        p: &RPoint3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<MaterialAt, Refusal> {
        match self.membership(p, storage, dec)?.location {
            Location::Outside => Ok(MaterialAt::Void),
            Location::Boundary => Ok(MaterialAt::Boundary),
            Location::Inside => match &self.definition.materials {
                Materials::Homogeneous { region, material } => Ok(MaterialAt::Region {
                    region: region.clone(),
                    material: material.clone(),
                }),
                Materials::Slabs { axis, regions } => {
                    dec.charge_work(regions.len() as u64)?;
                    let x = [&p.x, &p.y, &p.z][*axis];
                    let ids: Vec<_> = regions
                        .iter()
                        .filter(|r| &r.lo <= x && x <= &r.hi)
                        .map(|r| (r.id.clone(), r.material.clone()))
                        .collect();
                    Ok(match ids.as_slice() {
                        [(r, m)] => MaterialAt::Region {
                            region: r.clone(),
                            material: m.clone(),
                        },
                        [] => MaterialAt::Unresolved,
                        _ => MaterialAt::Interface(ids),
                    })
                }
            },
            _ => Ok(MaterialAt::Unresolved),
        }
    }
    fn piece(
        &self,
        b: Bounds,
        certain: bool,
        out: &mut Vec<Piece>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let scans = match &self.definition.materials {
            Materials::Homogeneous { .. } => 1,
            Materials::Slabs { regions, .. } => regions.len(),
        };
        dec.charge_work(scans as u64)?;
        match &self.definition.materials {
            Materials::Homogeneous { region, material } => out.push(Piece {
                bounds: b,
                certain,
                region: region.clone(),
                material: material.clone(),
            }),
            Materials::Slabs { axis, regions } => {
                for r in regions {
                    let mut s = b.clone();
                    s[*axis] = [r.lo.clone(), r.hi.clone()];
                    if let Some(bounds) = intersect(&b, &s) {
                        out.push(Piece {
                            bounds,
                            certain,
                            region: r.id.clone(),
                            material: r.material.clone(),
                        });
                    }
                }
            }
        }
        Ok(())
    }
    /// Bounded one-level octant refinement of unknown cells. No requested width
    /// is implied: fraction.width() and each moment.width() are achieved widths.
    pub fn integration(
        &self,
        b: &Bounds,
        storage: &mut BoxStorage,
        dec: &mut Decider,
    ) -> Result<Cell, Refusal> {
        dec.charge_work(1)?;
        let state = self.whole_cell(b, storage, dec)?;
        let mut pieces = Vec::new();
        if let Some(solid) = &self.exact_box {
            if let Some(clipped) = intersect(b, solid) {
                self.piece(clipped, true, &mut pieces, dec)?;
            }
        } else if state == BoxState::Inside {
            self.piece(b.clone(), true, &mut pieces, dec)?;
        } else if state == BoxState::Unknown {
            let mid: [Rational; 3] =
                std::array::from_fn(|i| (&b[i][0] + &b[i][1]).checked_div(&q(2)).unwrap());
            for k in 0..8 {
                dec.charge_work(1)?;
                let sub = std::array::from_fn(|i| {
                    if k >> i & 1 == 0 {
                        [b[i][0].clone(), mid[i].clone()]
                    } else {
                        [mid[i].clone(), b[i][1].clone()]
                    }
                });
                let s = self.whole_cell(&sub, storage, dec)?;
                if s != BoxState::Outside {
                    self.piece(sub, s == BoxState::Inside, &mut pieces, dec)?;
                }
            }
        }
        let mut v = Enclosure { lo: q(0), hi: q(0) };
        for p in &pieces {
            let w = volume(&p.bounds);
            v.hi = &v.hi + &w;
            if p.certain {
                v.lo = &v.lo + &w;
            }
        }
        let full = volume(b);
        Ok(Cell {
            revision: self.definition.revision.clone(),
            bounds: b.clone(),
            state,
            fraction: Enclosure {
                lo: v.lo.checked_div(&full).unwrap(),
                hi: v.hi.checked_div(&full).unwrap(),
            },
            connected_box: self.exact_box.is_some() && !pieces.is_empty(),
            pieces,
        })
    }
    /// Prepared trimmed-boundary distance and owning-name services. This separate
    /// preparation may refuse unsupported trims; cell integration remains usable.
    /// The search covers named faces only, so it refuses `solver.distance-coverage`
    /// when a live leaf (gyroid, sheet, lattice or mesh) has no named face: its
    /// boundary would be missed and a named-face bound is not a whole-boundary one.
    pub fn distances(&self, dec: &mut Decider) -> Result<distance::Prepared<'_>, Refusal> {
        // One pass over the compiled plan, already charged at compile; leaving it
        // uncharged keeps recorded distance-probe work unchanged.
        let plan = &self.definition.plan;
        let mut named = vec![false; plan.leaves.len()];
        for face in &plan.faces {
            named[face.leaf] = true;
        }
        if (0..plan.leaves.len()).any(|i| plan.leaf_registers[i].is_some() && !named[i]) {
            return Err(refusal(
                "solver.distance-coverage",
                "a live field leaf has no named boundary faces; whole-boundary distance is unproved",
            ));
        }
        let p = plan.prepare_distances(dec)?;
        p.prepare_nearest(dec)?;
        Ok(p)
    }
    /// Constant-size exact closest-point construction for the native axis box.
    /// The squared distance and displacement are exact. A unique nearest point
    /// need not have a unique surface normal (edge/corner); zero displacement has
    /// no direction. Tied interior faces retain every closest witness.
    pub fn box_distance(&self, p: &RPoint3, dec: &mut Decider) -> Result<BoxDistance, Refusal> {
        dec.charge_work(6)?;
        let b=self.exact_box.as_ref().ok_or_else(||refusal("solver.distance-coverage","fast closest construction supports one native axis box only; use trimmed distances for other geometry"))?;
        let x = [p.x.clone(), p.y.clone(), p.z.clone()];
        let clamped: [Rational; 3] =
            std::array::from_fn(|i| x[i].clone().max(b[i][0].clone()).min(b[i][1].clone()));
        let mut witnesses = Vec::new();
        if clamped != x {
            witnesses.push(clamped);
        } else {
            let min = (0..6)
                .map(|k| (&x[k / 2] - &b[k / 2][k % 2]).abs())
                .min()
                .unwrap();
            for i in 0..3 {
                for side in &b[i] {
                    if (&x[i] - side).abs() == min {
                        let mut w = x.clone();
                        w[i] = side.clone();
                        if !witnesses.contains(&w) {
                            witnesses.push(w);
                        }
                    }
                }
            }
        }
        let deltas: Vec<[Rational; 3]> = witnesses
            .iter()
            .map(|w| std::array::from_fn(|i| &w[i] - &x[i]))
            .collect();
        let squared = deltas[0].iter().fold(q(0), |s, d| s + d * d);
        Ok(BoxDistance {
            squared,
            displacements: deltas,
        })
    }
    /// A named boundary point is admitted only by the result's trimmed-face query.
    pub fn named_point(
        &self,
        name: &Origin,
        p: &RPoint3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<crate::FaceQuery, Refusal> {
        let id = self
            .definition
            .plan
            .select(name)
            .ok_or_else(|| refusal("solver.face-name", "name is absent from this revision"))?;
        self.definition.plan.face_query_checked(id, p, storage, dec)
    }
    /// Full-face load/support patch, confined to a single native axis box.
    /// Caller supplies a face-interior point; exact trim and unique normal are
    /// checked. No general curved patch integration is claimed.
    pub fn box_patch(
        &self,
        name: &Origin,
        p: &RPoint3,
        condition: Condition,
        dec: &mut Decider,
    ) -> Result<Patch, Refusal> {
        let b = self.exact_box.as_ref().ok_or_else(|| {
            refusal(
                "solver.patch-coverage",
                "full-face integration supports a native axis box only",
            )
        })?;
        let answer = self.named_point(name, p, &mut Workspace::default(), dec)?;
        if answer.state != crate::FaceState::OnFace {
            return Err(refusal(
                "solver.patch-membership",
                "point does not prove membership on selected face",
            ));
        }
        let x = [&p.x, &p.y, &p.z];
        let faces: Vec<_> = (0..3)
            .filter(|&i| x[i] == &b[i][0] || x[i] == &b[i][1])
            .collect();
        if faces.len() != 1
            || (0..3).any(|i| i != faces[0] && (x[i] <= &b[i][0] || x[i] >= &b[i][1]))
        {
            return Err(refusal(
                "solver.patch-ambiguous",
                "supply a strict face-interior point",
            ));
        }
        let axis = faces[0];
        let area = (0..3)
            .filter(|&i| i != axis)
            .fold(q(1), |s, i| s * (&b[i][1] - &b[i][0]));
        let centroid = std::array::from_fn(|i| {
            if i == axis {
                x[i].clone()
            } else {
                (&b[i][0] + &b[i][1]).checked_div(&q(2)).unwrap()
            }
        });
        Ok(Patch {
            revision: self.definition.revision.clone(),
            origin: name.clone(),
            area,
            centroid,
            condition,
        })
    }
}
/// Exact nearest boundary witness construction for one native axis box.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoxDistance {
    /// Exact squared unsigned boundary distance in square millimetres.
    pub squared: Rational,
    /// All nearest witness displacements in millimetres, ordered by axis/side on ties. A zero vector has no direction.
    pub displacements: Vec<[Rational; 3]>,
}
impl BoxDistance {
    /// Enclose the square root; no floating scalar is promoted to exact distance.
    pub fn distance(&self, width: &Rational, dec: &mut Decider) -> Result<Enclosure, Refusal> {
        if width <= &q(0) || self.squared < q(0) {
            return Err(refusal(
                "solver.distance-options",
                "distance needs positive width and nonnegative squared distance",
            ));
        }
        let b = dec.enclose(&Real::from(self.squared.clone()).sqrt(), width)?;
        Ok(Enclosure {
            lo: b.lo().max(q(0)),
            hi: b.hi(),
        })
    }
}
/// Boundary data in mm, N, N/mm²; no solver result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Condition {
    /// Zero displacement over the full patch.
    Fixed,
    /// Uniform world-space traction vector in N/mm².
    Traction([Rational; 3]),
}
/// An exactly integrated full named box face.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Patch {
    /// SHA-256 of canonical recipe and material assignment.
    pub revision: String,
    /// Persistent recipe face name.
    pub origin: Origin,
    /// Exact full-face area in square millimetres.
    pub area: Rational,
    /// Exact face area centroid in millimetres.
    pub centroid: [Rational; 3],
    /// Explicit boundary condition.
    pub condition: Condition,
}

/// Fixed canonical jobs of 128 cells, independent of worker count. Reserved work
/// is charged to the existing operation before spawning; unused units are never
/// refunded. Exhaustion discards the entire grid, not just its uncertain cells.
#[derive(Clone, Debug)]
pub struct GridOptions {
    /// Cells per axis, in 1..=64.
    pub divisions: usize,
    /// None uses every available core, bounded by the number of nonempty jobs.
    pub workers: Option<usize>,
    /// Units reserved per cell; allocation is fixed across worker counts.
    pub work_per_cell: u64,
}
impl Default for GridOptions {
    fn default() -> Self {
        Self {
            divisions: 64,
            workers: None,
            work_per_cell: 256,
        }
    }
}
/// Compact occupancy result; retrieve integration pieces by querying its bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct GridCell {
    /// Certified closed-cell state; Unknown is not a proved cut.
    pub state: BoxState,
    /// Enclosed occupied fraction in [0,1]; width is unresolved measure divided by cell volume.
    pub fraction: Enclosure,
}
/// Canonical ordered grid with a conservative whole-part volume enclosure.
#[derive(Clone, Debug, PartialEq)]
pub struct Grid {
    /// SHA-256 of canonical recipe and material assignment.
    pub revision: String,
    /// Exact closed axis-aligned bounds in millimetres.
    pub bounds: Bounds,
    /// Cells per axis, in 1..=64.
    pub divisions: usize,
    /// x varies fastest, then y, then z. Re-query integration for attributable
    /// material pieces and moments using cell_bounds(index).
    pub cells: Vec<GridCell>,
    /// Sum of enclosed fractions times exact cell volume, in cubic millimetres.
    pub volume: Enclosure,
    /// Actual numeric/classifier work performed inside reserved jobs, including refused jobs.
    pub actual_work: u64,
    /// Non-refundable work charged to the caller operation for all jobs.
    pub reserved_work: u64,
}
impl Grid {
    /// Exact bounds for one canonical cell, or None for an invalid index.
    pub fn cell_bounds(&self, index: usize) -> Option<Bounds> {
        (index < self.cells.len()).then(|| grid_bounds(&self.bounds, self.divisions, index))
    }
}
fn grid_bounds(bounds: &Bounds, n: usize, index: usize) -> Bounds {
    let indices = [index % n, index / n % n, index / (n * n)];
    std::array::from_fn(|i| {
        let h = (&bounds[i][1] - &bounds[i][0])
            .checked_div(&q(n as i64))
            .unwrap();
        [
            &bounds[i][0] + &(&h * &q(indices[i] as i64)),
            &bounds[i][0] + &(&h * &q(indices[i] as i64 + 1)),
        ]
    })
}
/// Execution metadata is deliberately separate from the semantic answer.
pub struct GridRun {
    /// Complete semantic answer or typed refusal; no partial grid escapes.
    pub answer: Result<Grid, Refusal>,
    /// Actual number of successfully spawned, nonempty workers.
    pub spawned_workers: usize,
    /// Actual numeric/classifier work performed inside reserved jobs, including refused jobs.
    pub actual_work: u64,
}
impl Queries<'_> {
    /// Build the complete grid using deterministic fixed jobs and an existing allowance.
    pub fn grid(&self, options: &GridOptions, dec: &mut Decider) -> GridRun {
        let mut spawned = 0;
        let mut actual_work = 0;
        let answer = (|| {
            if !operation_budget_active() {
                return Err(refusal(
                    "solver.operation-scope",
                    "grid requires the caller's existing operation allowance",
                ));
            }
            if options.divisions == 0
                || options.divisions > 64
                || options.workers == Some(0)
                || options.work_per_cell == 0
                || !valid(&self.bounds)
            {
                return Err(refusal(
                    "solver.grid-options",
                    "require 1..64 divisions, positive workers/work, and nonempty compact bounds",
                ));
            }
            dec.check_deadline()?;
            let count = options.divisions.pow(3);
            let reserved = options
                .work_per_cell
                .checked_mul(count as u64)
                .ok_or_else(|| refusal("solver.grid-options", "work allocation overflow"))?;
            let mut allowance = OperationAllowance::reserve(reserved)?;
            let jobs = count.div_ceil(128);
            let workers = options
                .workers
                .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, usize::from))
                .min(jobs);
            let mut assignments: Vec<Vec<_>> = (0..workers).map(|_| Vec::new()).collect();
            for job in 0..jobs {
                let start = job * 128;
                let end = (start + 128).min(count);
                let token = allowance
                    .split_off(options.work_per_cell * (end - start) as u64)
                    .unwrap();
                assignments[job % workers].push((job, start, end, token));
            }
            let policy = *dec.policy();
            let mut results = std::thread::scope(|scope| {
                let mut handles = Vec::new();
                let mut spawn_error = None;
                for assignment in assignments {
                    let task = move || {
                        assignment.into_iter().map(|(job, start, end, token)| {
                            let (result, work) = token.run(|| {
                                let mut dec = Decider::new(policy);
                                let mut storage = BoxStorage::default();
                                (start..end).map(|index| {
                                    let b = grid_bounds(&self.bounds, options.divisions, index);
                                    let c = self.integration(&b, &mut storage, &mut dec)?;
                                    Ok(GridCell { state: c.state, fraction: c.fraction })
                                }).collect::<Result<Vec<_>, Refusal>>()
                            });
                            let result = match (result, work.exhausted) {
                                // Preserve the original numeric refusal and its
                                // attached undecided evidence, including on exhaustion.
                                (Err(error), _) => Err(error),
                                (Ok(_), true) => Err(refusal("solver.job-budget",
                                    "canonical grid job exhausted its reserved work; no partial grid accepted")),
                                (Ok(value), false) => Ok(value),
                            };
                            (job, result, work.used)
                        }).collect::<Vec<_>>()
                    };
                    match std::thread::Builder::new().spawn_scoped(scope, task) {
                        Ok(h) => {
                            spawned += 1;
                            handles.push(h);
                        }
                        Err(_) => {
                            spawn_error = Some(refusal(
                                "solver.worker-spawn",
                                "unable to start a grid worker",
                            ));
                            break;
                        }
                    }
                }
                let mut result = Vec::new();
                for h in handles {
                    match h.join() {
                        Ok(v) => result.extend(v),
                        Err(_) => {
                            spawn_error = Some(refusal(
                                "solver.worker-panic",
                                "grid worker panicked; no answer accepted",
                            ))
                        }
                    }
                }
                spawn_error.map_or(Ok(result), Err)
            })?;
            results.sort_by_key(|r| r.0);
            actual_work = results.iter().map(|r| r.2).sum();
            let mut cells = Vec::with_capacity(count);
            let mut actual = 0;
            for (_, result, work) in results {
                actual += work;
                cells.extend(result?);
            }
            dec.check_deadline()?;
            let cell_volume = volume(&self.bounds).checked_div(&q(count as i64)).unwrap();
            let mut sum = Enclosure { lo: q(0), hi: q(0) };
            dec.charge_work(count as u64)?;
            for c in &cells {
                sum.lo = &sum.lo + &c.fraction.lo;
                sum.hi = &sum.hi + &c.fraction.hi;
            }
            Ok(Grid {
                revision: self.definition.revision.clone(),
                bounds: self.bounds.clone(),
                divisions: options.divisions,
                cells,
                volume: Enclosure {
                    lo: &sum.lo * &cell_volume,
                    hi: &sum.hi * &cell_volume,
                },
                actual_work: actual,
                reserved_work: reserved,
            })
        })();
        GridRun {
            answer,
            spawned_workers: spawned,
            actual_work,
        }
    }
}

fn identity(recipe: &Recipe, materials: &Materials) -> String {
    let mut text = recipe.to_text();
    let mut add = |s: &str| {
        text.push_str(&format!("{}:", s.len()));
        text.push_str(s);
    };
    match materials {
        Materials::Homogeneous { region, material } => {
            add("homogeneous");
            add(region);
            add(material);
        }
        Materials::Slabs { axis, regions } => {
            add("slabs");
            add(&axis.to_string());
            for r in regions {
                add(&r.id);
                add(&r.material);
                add(&r.lo.to_string());
                add(&r.hi.to_string());
            }
        }
    }
    catallax_import::sha256_hex(text.as_bytes())
}
