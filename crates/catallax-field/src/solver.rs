//! F7 no-mesh solver adapter. Certified geometry/integration enclosures are not
//! a certified elasticity solve. See `docs/SOLVER-QUERIES.md` for coverage.
use crate::{
    Direction, FaceId, Instruction, Location, Plan, Query, Workspace, distance, geometry,
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

/// Explicit homogeneous material, an ordered slab partition of the compact
/// bounding box, or body-step material regions declared in the recipe.
/// Material IDs carry no implicit constitutive law.
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
    /// Use the body-step `materials` regions declared directly on the [`Recipe`].
    Recipe,
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

struct CompiledRegion {
    id: String,
    material: String,
    plan: Plan,
    exact_box: Option<Bounds>,
}

fn node_bounds(plan: &Plan, node: usize, dec: &mut Decider) -> Result<Option<Bounds>, Refusal> {
    let mut tape: Vec<Option<Bounds>> = Vec::with_capacity(node + 1);
    for instruction in &plan.instructions[..=node] {
        let bound = match *instruction {
            Instruction::Empty => Some(std::array::from_fn(|_| [q(0), q(0)])),
            Instruction::Leaf(i) => plan.leaves[i]
                .display_box(dec)?
                .map(|b| b.map(|x| [x.lo(), x.hi()])),
            Instruction::Boolean { kind, a, b } => match kind {
                catallax_recipe::BooleanKind::Subtract => tape[a].clone(),
                catallax_recipe::BooleanKind::Union => match (&tape[a], &tape[b]) {
                    (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                        [
                            a[i][0].clone().min(b[i][0].clone()),
                            a[i][1].clone().max(b[i][1].clone()),
                        ]
                    })),
                    _ => None,
                },
                catallax_recipe::BooleanKind::Intersect => match (&tape[a], &tape[b]) {
                    (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                        [
                            a[i][0].clone().max(b[i][0].clone()),
                            a[i][1].clone().min(b[i][1].clone()),
                        ]
                    })),
                    (Some(a), None) | (None, Some(a)) => Some(a.clone()),
                    _ => None,
                },
            },
        };
        tape.push(bound);
    }
    Ok(tape[node].clone())
}

fn plan_exact_box(plan: &Plan, dec: &mut Decider) -> Result<Option<Bounds>, Refusal> {
    fn node_box(plan: &Plan, node: usize, dec: &mut Decider) -> Result<Option<Bounds>, Refusal> {
        match plan.instructions[node] {
            Instruction::Leaf(id) => Ok(match &plan.leaves[id] {
                geometry::Leaf::Primitive(p) => crate::mass::axis_block(p),
                _ => None,
            }),
            Instruction::Boolean {
                kind: catallax_recipe::BooleanKind::Union,
                a,
                b,
            } => {
                for (sub_node, tool_node) in [(a, b), (b, a)] {
                    if let Instruction::Boolean {
                        kind: catallax_recipe::BooleanKind::Subtract,
                        a: whole,
                        b: subtracted,
                    } = plan.instructions[sub_node]
                        && subtracted == tool_node
                        && let Some(whole_box) = node_box(plan, whole, dec)?
                        && let Some(tool_b) = node_bounds(plan, tool_node, dec)?
                        && (0..3).all(|i| {
                            whole_box[i][0] <= tool_b[i][0] && tool_b[i][1] <= whole_box[i][1]
                        })
                    {
                        return Ok(Some(whole_box));
                    }
                }
                Ok(None)
            }
            _ => Ok(None),
        }
    }
    node_box(plan, plan.root, dec)
}

/// Owns the immutable compiled recipe; its revision cannot be supplied separately
/// from the plan. This first cut accepts native recipes, not external file maps.
pub struct Definition {
    plan: Plan,
    revision: String,
    materials: Materials,
    recipe_regions: Vec<CompiledRegion>,
    all_recipe_boxes: bool,
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
        let mut recipe_regions = Vec::new();
        let mut all_recipe_boxes = false;
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
            Materials::Recipe if !recipe.materials.is_empty() => {
                let mut ids = std::collections::BTreeSet::new();
                for entry in &recipe.materials {
                    dec.charge_work(1)?;
                    if entry.id.as_str().is_empty()
                        || entry.material.trim().is_empty()
                        || !ids.insert(entry.id.as_str())
                    {
                        return Err(refusal(
                            "solver.material-partition",
                            "recipe material regions must be uniquely named with non-empty materials",
                        ));
                    }
                    let mut sub = recipe.clone();
                    sub.result = entry.region.clone();
                    sub.materials.clear();
                    let r_plan = Plan::compile_checked(&sub, dec)?;
                    if r_plan.is_proven_empty() {
                        return Err(refusal(
                            "solver.material-partition",
                            "recipe material region is empty",
                        ));
                    }
                    let r_bounds = printing::bounds(&r_plan, dec)?;
                    if intersect(&b, &r_bounds).is_none() {
                        return Err(refusal(
                            "solver.material-partition",
                            "recipe material region lies outside the solid bounds",
                        ));
                    }
                    let exact_box = plan_exact_box(&r_plan, dec)?;
                    recipe_regions.push(CompiledRegion {
                        id: entry.id.as_str().to_owned(),
                        material: entry.material.clone(),
                        plan: r_plan,
                        exact_box,
                    });
                }
                all_recipe_boxes = recipe_regions.iter().all(|r| r.exact_box.is_some());
                if let Some(part_box) = plan_exact_box(&plan, dec)?
                    && all_recipe_boxes
                {
                    let mut clipped_boxes = Vec::with_capacity(recipe_regions.len());
                    let mut total_vol = q(0);
                    for r in &recipe_regions {
                        let Some(c_box) = intersect(&part_box, r.exact_box.as_ref().unwrap())
                        else {
                            return Err(refusal(
                                "solver.material-partition",
                                "recipe material region does not intersect the solid",
                            ));
                        };
                        for prev in &clipped_boxes {
                            if intersect(&c_box, prev).is_some() {
                                return Err(refusal(
                                    "solver.material-partition",
                                    "recipe material regions overlap inside the solid",
                                ));
                            }
                        }
                        total_vol = &total_vol + &volume(&c_box);
                        clipped_boxes.push(c_box);
                    }
                    if total_vol != volume(&part_box) {
                        return Err(refusal(
                            "solver.material-partition",
                            "recipe material regions leave an unassigned gap inside the solid",
                        ));
                    }
                } else {
                    let part_cls = BoxClassifier::new(&plan, dec)?;
                    let reg_cls: Vec<_> = recipe_regions
                        .iter()
                        .map(|r| BoxClassifier::new(&r.plan, dec))
                        .collect::<Result<_, _>>()?;
                    let mut storage = BoxStorage::default();
                    let mut ws = Workspace::default();
                    let n = 4usize;
                    for idx in 0..n.pow(3) {
                        dec.charge_work(1)?;
                        let cell_b = grid_bounds(&b, n, idx);
                        let lo = std::array::from_fn(|i| cell_b[i][0].clone());
                        let hi = std::array::from_fn(|i| cell_b[i][1].clone());
                        let p_state = part_cls.classify(&lo, &hi, &mut storage, dec)?;
                        if p_state == BoxState::Outside {
                            continue;
                        }
                        let mut inside_count = 0usize;
                        let mut all_outside = true;
                        for cls in &reg_cls {
                            let rs = cls.classify(&lo, &hi, &mut storage, dec)?;
                            if rs != BoxState::Outside {
                                all_outside = false;
                            }
                            if rs == BoxState::Inside {
                                inside_count += 1;
                            }
                        }
                        if inside_count > 1 {
                            return Err(refusal(
                                "solver.material-partition",
                                "recipe material regions overlap inside the solid",
                            ));
                        }
                        if p_state == BoxState::Inside && all_outside {
                            return Err(refusal(
                                "solver.material-partition",
                                "recipe material regions leave an unassigned gap inside the solid",
                            ));
                        }
                        let center = RPoint3::new(
                            (&lo[0] + &hi[0]).mul_pow2(-1),
                            (&lo[1] + &hi[1]).mul_pow2(-1),
                            (&lo[2] + &hi[2]).mul_pow2(-1),
                        );
                        if plan.proven_checked(&center, &mut ws, dec)?.location == Location::Inside
                        {
                            let mut pt_inside = 0usize;
                            let mut pt_all_outside = true;
                            for r in &recipe_regions {
                                let rl = r.plan.proven_checked(&center, &mut ws, dec)?.location;
                                if rl != Location::Outside {
                                    pt_all_outside = false;
                                }
                                if rl == Location::Inside {
                                    pt_inside += 1;
                                }
                            }
                            if pt_inside > 1 {
                                return Err(refusal(
                                    "solver.material-partition",
                                    "recipe material regions overlap at an interior point",
                                ));
                            }
                            if pt_all_outside {
                                return Err(refusal(
                                    "solver.material-partition",
                                    "recipe material regions leave an unassigned gap at an interior point",
                                ));
                            }
                        }
                    }
                    for (i, r) in recipe_regions.iter().enumerate() {
                        let rb = printing::bounds(&r.plan, dec)?;
                        let center = RPoint3::new(
                            (&rb[0][0] + &rb[0][1]).mul_pow2(-1),
                            (&rb[1][0] + &rb[1][1]).mul_pow2(-1),
                            (&rb[2][0] + &rb[2][1]).mul_pow2(-1),
                        );
                        if r.plan.proven_checked(&center, &mut ws, dec)?.location
                            == Location::Inside
                            && plan.proven_checked(&center, &mut ws, dec)?.location
                                == Location::Inside
                        {
                            for (j, other) in recipe_regions.iter().enumerate() {
                                if i != j
                                    && other.plan.proven_checked(&center, &mut ws, dec)?.location
                                        == Location::Inside
                                {
                                    return Err(refusal(
                                        "solver.material-partition",
                                        "recipe material regions overlap inside the solid",
                                    ));
                                }
                            }
                        }
                    }
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
            recipe_regions,
            all_recipe_boxes,
        })
    }
    /// Compile native geometry using the `materials` regions declared on the [`Recipe`].
    pub fn from_recipe(recipe: &Recipe, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::compile(recipe, Materials::Recipe, dec)
    }
    /// Identity of this geometry and material assignment.
    pub fn revision(&self) -> &str {
        &self.revision
    }
    /// Prepare F6 coverage and the native-box shortcut without constructing a B-rep.
    pub fn prepare(&self, dec: &mut Decider) -> Result<Queries<'_>, Refusal> {
        let bounds = printing::bounds(&self.plan, dec)?;
        let exact_box = plan_exact_box(&self.plan, dec)?;
        let region_classifiers = self
            .recipe_regions
            .iter()
            .map(|r| BoxClassifier::new(&r.plan, dec))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Queries {
            definition: self,
            classifier: BoxClassifier::new(&self.plan, dec)?,
            region_classifiers,
            bounds,
            exact_box,
        })
    }
}

/// Prepared shared geometry. Workspaces and Deciders remain local to each caller.
pub struct Queries<'a> {
    definition: &'a Definition,
    classifier: BoxClassifier<'a>,
    region_classifiers: Vec<BoxClassifier<'a>>,
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
        let query = self.membership(p, storage, dec)?;
        match query.location {
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
                Materials::Recipe => self.recipe_material_at(p, true, storage, dec),
            },
            Location::Indeterminate => match &self.definition.materials {
                Materials::Recipe => self.recipe_material_at(p, false, storage, dec),
                _ => Ok(MaterialAt::Unresolved),
            },
        }
    }
    fn recipe_material_at(
        &self,
        p: &RPoint3,
        part_proven_inside: bool,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<MaterialAt, Refusal> {
        dec.charge_work(self.definition.recipe_regions.len() as u64)?;
        let mut inside = Vec::new();
        let mut boundary = Vec::new();
        let mut boundary_dirs = Vec::new();
        let mut indeterminate = false;
        for r in &self.definition.recipe_regions {
            let rq = r.plan.proven_checked(p, storage, dec)?;
            match rq.location {
                Location::Inside => inside.push((r.id.clone(), r.material.clone())),
                Location::Boundary => {
                    boundary.push((r.id.clone(), r.material.clone()));
                    boundary_dirs.push(rq.direction);
                }
                Location::Outside => (),
                Location::Indeterminate => indeterminate = true,
            }
        }
        if inside.len() > 1 {
            return Err(refusal(
                "solver.material-partition",
                "point lies inside multiple material regions",
            ));
        }
        if let [(region, material)] = inside.as_slice()
            && boundary.is_empty()
            && !indeterminate
        {
            return Ok(MaterialAt::Region {
                region: region.clone(),
                material: material.clone(),
            });
        }
        if inside.is_empty() && boundary.len() >= 2 && !indeterminate {
            if part_proven_inside {
                return Ok(MaterialAt::Interface(boundary));
            }
            // For a union of regions at a smooth curved interface, opposite
            // outward unit normals prove that the two regions cover both sides
            // of the surface, so the point is on an internal material interface.
            let mut opposite = false;
            for i in 0..boundary_dirs.len() {
                for j in (i + 1)..boundary_dirs.len() {
                    if let (Direction::Unique(a), Direction::Unique(b)) =
                        (&boundary_dirs[i], &boundary_dirs[j])
                    {
                        let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                        if dot.hi() < 0.0 {
                            opposite = true;
                            break;
                        }
                    }
                }
            }
            if opposite {
                return Ok(MaterialAt::Interface(boundary));
            }
        }
        if part_proven_inside && inside.is_empty() && boundary.is_empty() && !indeterminate {
            return Err(refusal(
                "solver.material-partition",
                "interior point is outside every declared material region",
            ));
        }
        Ok(MaterialAt::Unresolved)
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
            Materials::Recipe => self.definition.recipe_regions.len(),
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
            Materials::Recipe => {
                for r in &self.definition.recipe_regions {
                    if let Some(r_box) = &r.exact_box
                        && let Some(bounds) = intersect(&b, r_box)
                    {
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
    fn recipe_sub_pieces(
        &self,
        b: &Bounds,
        part_inside: bool,
        storage: &mut BoxStorage,
        out: &mut Vec<Piece>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        dec.charge_work(self.definition.recipe_regions.len() as u64)?;
        let lo = std::array::from_fn(|i| b[i][0].clone());
        let hi = std::array::from_fn(|i| b[i][1].clone());
        let mut states = Vec::with_capacity(self.region_classifiers.len());
        let mut inside_count = 0usize;
        let mut active_count = 0usize;
        for cls in &self.region_classifiers {
            let s = cls.classify(&lo, &hi, storage, dec)?;
            if s == BoxState::Inside {
                inside_count += 1;
            }
            if s != BoxState::Outside {
                active_count += 1;
            }
            states.push(s);
        }
        if inside_count > 1 {
            return Err(refusal(
                "solver.material-partition",
                "recipe material regions overlap on a cell",
            ));
        }
        if part_inside && active_count == 0 {
            return Err(refusal(
                "solver.material-partition",
                "recipe material regions leave an unassigned gap on a cell",
            ));
        }
        if part_inside && inside_count == 1 && active_count == 1 {
            for (r, s) in self.definition.recipe_regions.iter().zip(states) {
                if s == BoxState::Inside {
                    out.push(Piece {
                        bounds: b.clone(),
                        certain: true,
                        region: r.id.clone(),
                        material: r.material.clone(),
                    });
                }
            }
            return Ok(true);
        }
        for (r, s) in self.definition.recipe_regions.iter().zip(states) {
            if s != BoxState::Outside {
                out.push(Piece {
                    bounds: b.clone(),
                    certain: part_inside && s == BoxState::Inside,
                    region: r.id.clone(),
                    material: r.material.clone(),
                });
            }
        }
        Ok(false)
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
        let curved_recipe_regions = matches!(self.definition.materials, Materials::Recipe)
            && !self.definition.all_recipe_boxes;
        let mut v = Enclosure { lo: q(0), hi: q(0) };
        if !curved_recipe_regions {
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
            for p in &pieces {
                let w = volume(&p.bounds);
                v.hi = &v.hi + &w;
                if p.certain {
                    v.lo = &v.lo + &w;
                }
            }
        } else if let Some(solid) = &self.exact_box {
            if let Some(clipped) = intersect(b, solid) {
                let full_v = volume(&clipped);
                v = Enclosure {
                    lo: full_v.clone(),
                    hi: full_v,
                };
                let mut trial = Vec::new();
                if self.recipe_sub_pieces(&clipped, true, storage, &mut trial, dec)? {
                    pieces = trial;
                } else {
                    let mid: [Rational; 3] = std::array::from_fn(|i| {
                        (&clipped[i][0] + &clipped[i][1])
                            .checked_div(&q(2))
                            .unwrap()
                    });
                    for k in 0..8 {
                        dec.charge_work(1)?;
                        let sub = std::array::from_fn(|i| {
                            if k >> i & 1 == 0 {
                                [clipped[i][0].clone(), mid[i].clone()]
                            } else {
                                [mid[i].clone(), clipped[i][1].clone()]
                            }
                        });
                        self.recipe_sub_pieces(&sub, true, storage, &mut pieces, dec)?;
                    }
                }
            }
        } else if state == BoxState::Inside {
            let full_v = volume(b);
            v = Enclosure {
                lo: full_v.clone(),
                hi: full_v,
            };
            let mut trial = Vec::new();
            if self.recipe_sub_pieces(b, true, storage, &mut trial, dec)? {
                pieces = trial;
            } else {
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
                    self.recipe_sub_pieces(&sub, true, storage, &mut pieces, dec)?;
                }
            }
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
                    let w = volume(&sub);
                    v.hi = &v.hi + &w;
                    if s == BoxState::Inside {
                        v.lo = &v.lo + &w;
                    }
                    self.recipe_sub_pieces(&sub, s == BoxState::Inside, storage, &mut pieces, dec)?;
                }
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
    /// Q4: general unsigned true-boundary distance and direction semantics,
    /// reusing a prepared distance service.
    pub fn boundary_distance_with(
        &self,
        prepared: &distance::Prepared<'_>,
        p: &RPoint3,
        options: &distance::NearestOptions,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<BoundaryDistance, Refusal> {
        if options.width <= q(0) || options.max_cells < 2 {
            return Err(refusal(
                "solver.distance-options",
                "boundary distance requires positive width and at least two cells",
            ));
        }
        let query = self.membership(p, storage, dec)?;
        if query.location == Location::Boundary {
            return Ok(BoundaryDistance {
                revision: self.definition.revision.clone(),
                location: Some(Location::Boundary),
                lower: q(0),
                upper: Some(q(0)),
                witness: Some(p.clone().to_real()),
                converged: true,
                direction: BoundaryDirection::BoundaryNormal(query.direction),
                cells: 0,
            });
        }
        let point = p.clone().to_real();
        let ownership = prepared.nearest_face(&point, options, dec)?;
        match ownership {
            distance::Ownership::Proved(owner) => {
                let mut lower = owner.lower;
                let mut upper = owner.upper;
                let mut witness = owner.witness;
                let mut cells = owner.cells;
                if &upper - &lower > options.width {
                    let refined = prepared.point_to_face(
                        owner.face,
                        &point,
                        &catallax_shapes::distance::Options {
                            width: options.width.clone(),
                            max_cells: options.max_cells,
                        },
                        dec,
                    )?;
                    cells += refined.cells;
                    lower = lower.max(refined.lower);
                    if let Some(u) = refined.upper
                        && u < upper
                    {
                        upper = u;
                        if let Some([_, w]) = refined.witnesses {
                            witness = w;
                        }
                    }
                }
                let converged = &upper - &lower <= options.width;
                Ok(BoundaryDistance {
                    revision: self.definition.revision.clone(),
                    location: Some(query.location),
                    lower,
                    upper: Some(upper),
                    witness: Some(witness),
                    converged,
                    direction: BoundaryDirection::UniqueNearest {
                        owner: owner.origin,
                        toward_witness: owner.direction,
                        margin: owner.margin,
                    },
                    cells,
                })
            }
            distance::Ownership::Refused(ambiguity) => {
                let bounds = prepared.point_to_boundary(
                    &point,
                    &catallax_shapes::distance::Options {
                        width: options.width.clone(),
                        max_cells: options.max_cells,
                    },
                    dec,
                )?;
                let witness = bounds.witnesses.map(|[_, w]| w);
                let direction = if let Some(w) = &witness {
                    let offset = w.minus(&point);
                    let toward = if geometry::sign(&offset.norm_squared(), dec)?
                        == catallax_numbers::Sign::Zero
                    {
                        Direction::Undefined
                    } else {
                        crate::unit_direction(offset, dec)?
                    };
                    BoundaryDirection::WitnessOnly {
                        toward_witness: toward,
                        contenders: ambiguity.contenders.into_iter().map(|c| c.origin).collect(),
                    }
                } else {
                    BoundaryDirection::Undefined
                };
                Ok(BoundaryDistance {
                    revision: self.definition.revision.clone(),
                    location: Some(query.location),
                    lower: bounds.lower,
                    upper: bounds.upper,
                    witness,
                    converged: bounds.converged,
                    direction,
                    cells: ambiguity.cells + bounds.cells,
                })
            }
        }
    }
    /// Q4: general unsigned true-boundary distance and direction semantics.
    pub fn boundary_distance(
        &self,
        p: &RPoint3,
        options: &distance::NearestOptions,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<BoundaryDistance, Refusal> {
        let prepared = self.distances(dec)?;
        self.boundary_distance_with(&prepared, p, options, storage, dec)
    }
    /// Q3/Q5: integrate a load or support condition on any surviving named face,
    /// reusing a prepared distance service.
    pub fn face_patch_with(
        &self,
        prepared: &distance::Prepared<'_>,
        name: &Origin,
        condition: Condition,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<FacePatch, Refusal> {
        let integral = prepared.integrate_origin(name, width, dec)?;
        let zero = || Enclosure { lo: q(0), hi: q(0) };
        let mul_enc = |a: &Enclosure, b: &Enclosure| -> Enclosure {
            let products = [&a.lo * &b.lo, &a.lo * &b.hi, &a.hi * &b.lo, &a.hi * &b.hi];
            Enclosure {
                lo: products.iter().min().unwrap().clone(),
                hi: products.iter().max().unwrap().clone(),
            }
        };
        let sub_enc = |a: &Enclosure, b: &Enclosure| -> Enclosure {
            Enclosure {
                lo: &a.lo - &b.hi,
                hi: &a.hi - &b.lo,
            }
        };
        let cross_enc = |a: &[Enclosure; 3], b: &[Enclosure; 3]| -> [Enclosure; 3] {
            [
                sub_enc(&mul_enc(&a[1], &b[2]), &mul_enc(&a[2], &b[1])),
                sub_enc(&mul_enc(&a[2], &b[0]), &mul_enc(&a[0], &b[2])),
                sub_enc(&mul_enc(&a[0], &b[1]), &mul_enc(&a[1], &b[0])),
            ]
        };
        let exact_enc = |r: &Rational| Enclosure {
            lo: r.clone(),
            hi: r.clone(),
        };
        let (resultant_force, resultant_moment) = match &condition {
            Condition::Fixed => ([zero(), zero(), zero()], [zero(), zero(), zero()]),
            Condition::Traction(t) => {
                let t_enc = t.each_ref().map(exact_enc);
                let force = std::array::from_fn(|k| mul_enc(&integral.area, &t_enc[k]));
                let moment = cross_enc(&integral.first_moment, &t_enc);
                (force, moment)
            }
            Condition::TotalForce(f) => {
                let f_enc = f.each_ref().map(exact_enc);
                let inv_area = Enclosure {
                    lo: q(1).checked_div(&integral.area.hi).unwrap(),
                    hi: q(1).checked_div(&integral.area.lo).unwrap(),
                };
                let traction: [Enclosure; 3] =
                    std::array::from_fn(|k| mul_enc(&f_enc[k], &inv_area));
                let force = std::array::from_fn(|k| mul_enc(&integral.area, &traction[k]));
                let moment = cross_enc(&integral.centroid, &f_enc);
                (force, moment)
            }
            Condition::Pressure(p) => {
                let neg_p = exact_enc(&(-p));
                let force = std::array::from_fn(|k| mul_enc(&integral.normal_integral[k], &neg_p));
                let moment = std::array::from_fn(|k| mul_enc(&integral.normal_moment[k], &neg_p));
                (force, moment)
            }
        };
        Ok(FacePatch {
            revision: self.definition.revision.clone(),
            origin: name.clone(),
            faces: integral.faces,
            pieces: integral.pieces,
            area: integral.area,
            centroid: integral.centroid,
            normal_integral: integral.normal_integral,
            unit_normal: integral.unit_normal,
            condition,
            resultant_force,
            resultant_moment,
        })
    }
    /// Q3/Q5: integrate a load or support condition on any surviving named face
    /// (planar or curved, including split faces).
    pub fn face_patch(
        &self,
        name: &Origin,
        condition: Condition,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<FacePatch, Refusal> {
        let prepared = self.distances(dec)?;
        self.face_patch_with(&prepared, name, condition, width, dec)
    }
}

/// Direction semantics for a general boundary-distance query (Q4).
#[derive(Clone, Debug, PartialEq)]
pub enum BoundaryDirection {
    /// Certified boundary point (zero distance): carries the outward surface normal
    /// (`Direction::Unique` on a smooth face interior, `Direction::Ambiguous` on an edge/corner).
    BoundaryNormal(Direction),
    /// Off-surface point with one proved nearest named face: unit vector from the query
    /// point toward the admitted trimmed witness on that face, plus the certified
    /// lower bound on distance to the nearest ownership border when other faces exist.
    UniqueNearest {
        /// Proved nearest face name.
        owner: Origin,
        /// Unit vector from the query point toward the witness.
        toward_witness: Direction,
        /// Certified lower bound on distance to the nearest ownership border.
        margin: Option<Rational>,
    },
    /// Off-surface point where an admitted witness is available, but multiple named faces
    /// remain unseparated within the requested width (e.g. outside a convex edge or on a
    /// medial tie locus). The ray toward the witness is distinct from a unique nearest direction.
    WitnessOnly {
        /// Unit vector from the query point toward the best admitted witness.
        toward_witness: Direction,
        /// Unseparated contending face names.
        contenders: Vec<Origin>,
    },
    /// No certified direction (e.g. no admitted witness or unresolved boundary).
    Undefined,
}

/// General unsigned true-boundary distance and direction record (Q4).
#[derive(Clone, Debug)]
pub struct BoundaryDistance {
    /// SHA-256 of canonical recipe and material assignment.
    pub revision: String,
    /// Certified point classification when the query point is rational.
    pub location: Option<Location>,
    /// Proven unsigned Euclidean distance lower bound in millimetres.
    pub lower: Rational,
    /// Proven unsigned Euclidean distance upper bound in millimetres, when a witness is admitted.
    pub upper: Option<Rational>,
    /// Exact admitted boundary witness attaining `upper`, when available.
    pub witness: Option<catallax_shapes::Point3>,
    /// Whether `upper - lower <= options.width`.
    pub converged: bool,
    /// Explicit direction semantics: boundary outward normal, unique-owner witness ray,
    /// tied/ambiguous witness ray, or undefined.
    pub direction: BoundaryDirection,
    /// Total chart cells visited across boundary distance and ownership separation.
    pub cells: usize,
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
    /// Total world-space force vector in N, distributed as uniform traction over the face area.
    TotalForce([Rational; 3]),
    /// Uniform inward normal pressure in N/mm² (positive acts along -n_outward).
    Pressure(Rational),
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
/// Integrated boundary condition on any surviving named face (planar or curved,
/// including split faces carrying the same recipe origin).
#[derive(Clone, Debug, PartialEq)]
pub struct FacePatch {
    /// SHA-256 of canonical recipe and material assignment.
    pub revision: String,
    /// Persistent recipe face name.
    pub origin: Origin,
    /// Surviving registry face handles contributing to this named selection.
    pub faces: Vec<FaceId>,
    /// Number of connected boundary pieces integrated for this origin.
    pub pieces: usize,
    /// Certified total surface area enclosure in square millimetres.
    pub area: Enclosure,
    /// Certified area centroid enclosure in millimetres.
    pub centroid: [Enclosure; 3],
    /// Certified integral of the outward unit normal over the face, in square millimetres.
    pub normal_integral: [Enclosure; 3],
    /// Unit outward normal when constant across the patch (planar faces).
    pub unit_normal: Option<[Enclosure; 3]>,
    /// Explicit boundary condition.
    pub condition: Condition,
    /// Assembled resultant force vector in newtons (zero for `Fixed`).
    pub resultant_force: [Enclosure; 3],
    /// Assembled resultant moment about the world origin in N·mm (zero for `Fixed`).
    pub resultant_moment: [Enclosure; 3],
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
    /// Number of nonempty assignments admitted to the shared field pool.
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
                .unwrap_or_else(crate::workers::worker_count)
                .min(crate::workers::worker_count())
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
            let assignments = assignments
                .into_iter()
                .map(|assignment| std::sync::Mutex::new(Some(assignment)))
                .collect::<Vec<_>>();
            spawned = assignments.len();
            let mut results = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                crate::workers::with_workers(workers, || {
                    crate::workers::map(&assignments, |assignment| {
                        assignment.lock().unwrap().take().expect("single-use grid jobs")
                            .into_iter().map(|(job, start, end, token)| {
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
                                    (Err(error), _) => Err(error),
                                    (Ok(_), true) => Err(refusal("solver.job-budget",
                                        "canonical grid job exhausted its reserved work; no partial grid accepted")),
                                    (Ok(value), false) => Ok(value),
                                };
                                (job, result, work.used)
                            }).collect::<Vec<_>>()
                    })
                })
            })).map_err(|_| refusal("solver.worker-panic", "grid worker panicked; no answer accepted"))??
                .into_iter().flatten().collect::<Vec<_>>();
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
        Materials::Recipe => {
            add("recipe");
            for r in &recipe.materials {
                add(r.id.as_str());
                add(r.region.as_str());
                add(&r.material);
            }
        }
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
