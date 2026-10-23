//! Recipe-defined analytic fields in safe Rust.
//!
//! [`Plan`] owns a fixed-order evaluation tape and a separate named [`Face`]
//! registry. [`Plan::fast`] is tentative display work. [`Plan::proven`] uses
//! exact support and profile predicates and outward-rounded binary balls, with
//! no platform math library. Sweep and Boolean residuals are not true distance.
//! See `docs/FIELD.md` for the distance proof and unresolved-boundary contract.

#![forbid(unsafe_code)]

mod certified_prism;
/// Definition-native, certified display approximations. No B-rep is constructed.
pub mod display;
/// Certified trimmed-face and boundary minimum-distance queries.
pub mod distance;
mod draft;
mod features;
mod fillet;
mod geometry;
pub use certified_prism::CertifiedPrismField;
/// Proven two-sided deviation of B-rep faces from their own definitions.
pub mod measurement;
mod metric;
mod offset_spline;
mod planar;
mod profile;
mod samples;
mod support;

use std::collections::BTreeMap;

use catallax_numbers::{Ball, Decider, Interval, Rational, Real, RealPoly, Refusal, Sign};
use catallax_recipe::{BooleanKind, Operation, Origin, Recipe};
use catallax_shapes::{Isometry, Pose, RPoint3, RVector3, Surface, Vector3};

use metric::{PRECISION, Scalar};

/// Certified point classification. No tolerance changes a boundary into material.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Location {
    /// In the interior of the regularised solid.
    Inside,
    /// Outside the closed regularised solid.
    Outside,
    /// On the solid's boundary, not merely a residual's zero set.
    Boundary,
    /// The exact trim or deterministic local boundary proof did not settle.
    Indeterminate,
}

/// The numerical meaning of an evaluated field.
#[derive(Clone, Copy, PartialEq, Debug)]
#[non_exhaustive]
pub enum Value {
    /// Enclosure of true signed Euclidean distance, negative inside, in millimetres.
    TrueDistance(Interval),
    /// Enclosure of a sweep or Boolean residual in millimetres. Not true distance.
    Residual(Interval),
}

impl Value {
    /// The enclosure, without changing its meaning.
    pub fn enclosure(self) -> Interval {
        match self {
            Self::TrueDistance(v) | Self::Residual(v) => v,
        }
    }
}

/// Enclosure of unsigned distance to the result boundary, in millimetres.
/// An infinite upper end is explicitly unbounded, not a 1 nm certificate.
/// The empty solid has distance `[+infinity, +infinity]`.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct DistanceBound {
    /// Non-negative lower and upper bounds.
    pub millimetres: Interval,
}

/// Outward unit surface direction. Defined only at boundary points.
#[derive(Clone, Copy, PartialEq, Debug)]
#[non_exhaustive]
pub enum Direction {
    /// A unique outward unit normal, enclosed component by component.
    Unique([Interval; 3]),
    /// Edge, corner, apex, or a local configuration with no certified unique normal.
    Ambiguous,
    /// Not a certified boundary point.
    Undefined,
}

/// One authoritative point query.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Query {
    /// Exact classification or an explicit unresolved local boundary.
    pub location: Location,
    /// A true distance or a labelled residual.
    pub value: Value,
    /// A geometric unsigned distance enclosure, distinct from the residual.
    pub distance: DistanceBound,
    /// Surface direction where certified.
    pub direction: Direction,
}

/// Tentative display result. It carries no proof or geometric error certificate.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Approximate {
    /// Negative inside, positive outside. Zero need not be a Boolean boundary.
    pub value: f64,
    /// Whether the real-arithmetic formula is true signed distance.
    /// Floating rounding still makes `value` only an approximation.
    pub true_distance_formula: bool,
}

/// A selectable recipe face candidate. Surviving trimmed domains are queried
/// against the result, never inferred from a winning tape instruction alone.
#[derive(Clone, Debug)]
pub struct Face {
    /// Stable name from the creating recipe operation.
    pub origin: Origin,
    /// Exact world-space supporting surface, independent of the evaluation tape.
    pub surface: Surface,
    polynomial: RealPoly,
    support: RealPoly,
    gradient: [RealPoly; 3],
    leaf: usize,
    in_result: bool,
}

/// Handle into one plan's immutable face registry. Not a persistent recipe name.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct FaceId(usize);

impl FaceId {
    /// Index in the owning plan's face registry, never a persistent recipe name.
    pub fn index(self) -> usize {
        self.0
    }
}

/// Whether a supporting-surface point survives Boolean trimming.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum FaceState {
    /// Certified point on the surviving face.
    OnFace,
    /// Not on the supporting surface or outside the surviving trimmed domain.
    NotOnFace,
    /// Degenerate trim or edge incidence needs further geometric proof.
    Indeterminate,
}

/// Selection query on a named face, including subtraction's orientation reversal.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct FaceQuery {
    /// Trim membership, never merely a supporting-surface hit.
    pub state: FaceState,
    /// Outward normal of this face in the result, when unique and proven.
    pub direction: Direction,
}

#[derive(Clone, Copy, Debug)]
enum Instruction {
    Empty,
    Leaf(usize),
    Boolean {
        kind: BooleanKind,
        a: usize,
        b: usize,
    },
}

/// Reusable per-worker storage. One immutable plan can serve any number of workers.
/// Fast point evaluation allocates nothing after this storage has warmed up.
#[derive(Default, Debug)]
pub struct Workspace {
    fast: Vec<f64>,
    values: Vec<Ball>,
    locations: Vec<Location>,
    rays: Vec<Location>,
    active: Vec<usize>,
    planes: Vec<Vector3>,
}

/// Reusable batch storage. Results always have the input points' order.
#[derive(Default, Debug)]
pub struct Batch {
    workspace: Workspace,
    proven: Vec<Query>,
    fast: Vec<Approximate>,
}

/// An owned, fixed-order field plan compiled from a checked recipe.
/// Face identities remain outside the tape, including folded-away candidates.
#[derive(Clone, Debug)]
pub struct Plan {
    instructions: Vec<Instruction>,
    leaves: Vec<geometry::Leaf>,
    faces: Vec<Face>,
    root: usize,
    leaf_registers: Vec<Option<usize>>,
    // Only explicit canonical expansions rename roles; domain predicates retain
    // the generated feature role, independently of the public stock name.
    definition_roles: BTreeMap<usize, catallax_recipe::Role>,
    consumed: Vec<catallax_features::ConsumedFace>,
}

impl Plan {
    /// Removed-face provenance from draft operations, in recipe order.
    /// These records never make a vanished face selectable.
    pub fn consumed_faces(&self) -> &[catallax_features::ConsumedFace] {
        &self.consumed
    }

    /// Compiles stock solids, sketch sweeps and regularised booleans.
    /// Unsupported operations and impossible dimensions report a stable refusal.
    /// Sketch validation is shared with features; no B-rep builder is called.
    pub fn compile(recipe: &Recipe) -> Result<Self, Refusal> {
        Self::compile_checked(recipe, &mut Decider::default())
    }

    /// Compiles under the shared deterministic proof policy. Defining domains
    /// and support-normalisation divisors are proved before publishing a plan.
    pub fn compile_checked(recipe: &Recipe, dec: &mut Decider) -> Result<Self, Refusal> {
        recipe
            .check()
            .map_err(|e| Refusal::new("field.recipe", e.to_string()))?;
        let mut plan = Self {
            instructions: vec![Instruction::Empty],
            leaves: vec![],
            faces: vec![],
            root: 0,
            leaf_registers: vec![],
            definition_roles: BTreeMap::new(),
            consumed: Vec::new(),
        };
        let mut steps = BTreeMap::new();
        let mut profiles = BTreeMap::new();
        let mut step_faces: BTreeMap<catallax_recipe::StepId, Vec<usize>> = BTreeMap::new();
        for step in &recipe.steps {
            let span =
                catallax_numbers::OperationSpan::begin("field", &step.id, step.op.name(), dec);
            let first_face = plan.faces.len();
            let result = (|| -> Result<Option<usize>, Refusal> {
                let node = match &step.op {
                    Operation::Sketch { plane, loops } => {
                        profiles.insert(step.id.clone(), profile::Profile::new(plane, loops, dec)?);
                        return Ok(None);
                    }
                    Operation::Boolean { kind, target, tool } => {
                        plan.boolean(*kind, steps[target], steps[tool], dec)?
                    }
                    Operation::Extrude {
                        profile,
                        distance,
                        direction,
                    } => {
                        let built = features::Feature::extrude(
                            &profiles[profile],
                            &step.id,
                            Some(distance.value()),
                            direction.as_ref(),
                            false,
                            dec,
                        )?;
                        plan.feature(built, dec)?
                    }
                    Operation::ExtrudeThrough {
                        target,
                        profile,
                        direction,
                        side,
                    } => {
                        let built = features::Feature::extrude(
                            &profiles[profile],
                            &step.id,
                            None,
                            direction.as_ref(),
                            steps[target] == 0,
                            dec,
                        )?;
                        let mut tool = plan.feature(built, dec)?;
                        if let Some(
                            side @ (catallax_recipe::ThroughSide::Forward
                            | catallax_recipe::ThroughSide::Backward),
                        ) = side
                        {
                            let profile = &profiles[profile];
                            let normal = profile::rational_vector(
                                &profile.pose,
                                &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
                            );
                            let direction = direction
                                .as_ref()
                                .map(|v| {
                                    let [x, y, z] = v.0.each_ref().map(|n| n.value().clone());
                                    RVector3::new(x, y, z)
                                })
                                .unwrap_or_else(|| normal.clone());
                            let positive = direction.dot(&normal).sign() == Sign::Positive;
                            let forward = *side == catallax_recipe::ThroughSide::Forward;
                            let plane = Surface::Plane(catallax_shapes::Plane {
                                origin: profile
                                    .pose
                                    .apply_point(&catallax_shapes::Point3::origin()),
                                normal: normal
                                    .scaled(&Rational::from_i64(if positive == forward {
                                        -1
                                    } else {
                                        1
                                    }))
                                    .into(),
                                ref_dir: profile.pose.apply_vector(
                                    &RVector3::new(
                                        Rational::one(),
                                        Rational::zero(),
                                        Rational::zero(),
                                    )
                                    .to_real(),
                                ),
                            });
                            let face = geometry::face(
                                plane.clone(),
                                Origin {
                                    step: step.id.clone(),
                                    role: catallax_recipe::Role::Start,
                                },
                                dec,
                            )?;
                            let cap = plan.leaf(
                                geometry::Leaf::Cell(geometry::Cell::new(vec![plane], None)?),
                                vec![face],
                                dec,
                            )?;
                            tool = plan.boolean(BooleanKind::Intersect, tool, cap, dec)?;
                        }
                        plan.boolean(BooleanKind::Subtract, steps[target], tool, dec)?
                    }
                    Operation::Revolve {
                        profile,
                        axis,
                        angle,
                    } => {
                        let built = features::Feature::revolve(
                            &profiles[profile],
                            &step.id,
                            axis,
                            angle.value(),
                            dec,
                        )?;
                        plan.feature(built, dec)?
                    }
                    Operation::Pattern { body, layout } => {
                        catallax_features::pattern_placement(layout, 0, dec)?;
                        let source = steps[body];
                        let mut result = source;
                        for instance in 1..layout.count() {
                            let map = catallax_features::pattern_placement(layout, instance, dec)?;
                            let copy = plan.copy(
                                source,
                                &step_faces[body],
                                &map,
                                &step.id,
                                instance,
                                dec,
                            )?;
                            result = plan.boolean(BooleanKind::Union, result, copy, dec)?;
                        }
                        result
                    }
                    Operation::Mirror {
                        body,
                        plane,
                        keep_original,
                    } => {
                        let map = catallax_features::mirror_placement(plane)?;
                        let source = steps[body];
                        let copy = plan.copy(source, &step_faces[body], &map, &step.id, 1, dec)?;
                        if *keep_original {
                            plan.boolean(BooleanKind::Union, source, copy, dec)?
                        } else {
                            copy
                        }
                    }
                    Operation::Draft {
                        body,
                        faces,
                        pull,
                        neutral,
                        angle,
                        keep_regions,
                    } => {
                        let [x, y, z] = pull.0.each_ref().map(|n| n.value().clone());
                        let definition = catallax_features::draft::DraftDefinition::new(
                            RVector3::new(x, y, z).to_real(),
                            catallax_features::draft::neutral_plane(neutral),
                            angle.value(),
                            dec,
                        )?;
                        plan.draft(
                            (steps[body], &step_faces[body]),
                            faces,
                            &definition,
                            keep_regions.as_deref(),
                            &step.id,
                            dec,
                        )?
                    }
                    Operation::Fillet {
                        body,
                        edges,
                        radius,
                    } => plan.fillet(
                        steps[body],
                        &step_faces[body],
                        edges,
                        &radius.value().clone().into(),
                        &step.id,
                        dec,
                    )?,
                    Operation::ExtrudeTo {
                        target,
                        profile,
                        direction,
                        ..
                    }
                    | Operation::ExtrudeNext {
                        target,
                        profile,
                        direction,
                    } => {
                        let limit = if let Operation::ExtrudeTo { limit, .. } = &step.op {
                            Some(limit)
                        } else {
                            None
                        };
                        let cap = catallax_features::end_condition::Capped::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            limit,
                            dec,
                        )?;
                        let built = features::Feature::extrude(
                            &profiles[profile],
                            &step.id,
                            None,
                            direction.as_ref(),
                            false,
                            dec,
                        )?;
                        let tool = plan.feature(built, dec)?;
                        let mut envelope = 0;
                        for cell in cap.cells.iter().filter(|c| c.active) {
                            let planes = cap.cell_planes(cell, dec)?;
                            let mut faces = Vec::new();
                            for (plane, role) in &planes {
                                if let Some(role) = role {
                                    faces.push(geometry::face(
                                        Surface::Plane(plane.clone()),
                                        Origin {
                                            step: step.id.clone(),
                                            role: role.clone(),
                                        },
                                        dec,
                                    )?);
                                }
                            }
                            let surfaces =
                                planes.into_iter().map(|(p, _)| Surface::Plane(p)).collect();
                            let node = plan.leaf(
                                geometry::Leaf::Cell(geometry::Cell::new(
                                    surfaces,
                                    Some(cap.column_radius.clone()),
                                )?),
                                faces,
                                dec,
                            )?;
                            envelope = plan.boolean(BooleanKind::Union, envelope, node, dec)?;
                        }
                        plan.boolean(BooleanKind::Intersect, tool, envelope, dec)?
                    }
                    Operation::Hole {
                        target,
                        at,
                        diameter,
                        depth,
                        tip,
                        head,
                    } => {
                        catallax_features::hole::validate(
                            diameter.value(),
                            depth,
                            tip,
                            head.as_ref(),
                        )?;
                        let target = steps[target];
                        let bound = plan.node_radius(target, dec)?;
                        let mut tool_node = 0;
                        for placement in at {
                            let tool = catallax_features::hole::HoleTool::new(
                                placement,
                                diameter.value(),
                                depth,
                                tip,
                                head.as_ref(),
                                &bound,
                                dec,
                            )?;
                            for band in &tool.bands {
                                let side = tool.surface(band);
                                let lower = tool.cap(&band.depth[0], false);
                                let upper = tool.cap(&band.depth[1], true);
                                let mut faces = vec![geometry::face(
                                    side.clone(),
                                    Origin {
                                        step: step.id.clone(),
                                        role: catallax_recipe::Role::Hole {
                                            at: placement.id.clone(),
                                            part: band.part,
                                        },
                                    },
                                    dec,
                                )?];
                                for part in &band.end {
                                    faces.push(geometry::face(
                                        upper.clone(),
                                        Origin {
                                            step: step.id.clone(),
                                            role: catallax_recipe::Role::Hole {
                                                at: placement.id.clone(),
                                                part: *part,
                                            },
                                        },
                                        dec,
                                    )?);
                                }
                                let mut radius = placement
                                    .point
                                    .0
                                    .iter()
                                    .fold(Rational::zero(), |s, n| s + n.value().abs());
                                let mut axial = Rational::zero();
                                let mut radial = Rational::zero();
                                for x in &band.depth {
                                    let b = geometry::enclose(x, 128, dec)?;
                                    axial = axial.max(b.lo().abs().max(b.hi().abs()));
                                }
                                for x in &band.radius {
                                    radial = radial.max(geometry::enclose(x, 128, dec)?.hi());
                                }
                                radius = radius + axial + radial;
                                let node = plan.leaf(
                                    geometry::Leaf::Cell(geometry::Cell::new(
                                        vec![lower, upper, side],
                                        Some(radius),
                                    )?),
                                    faces,
                                    dec,
                                )?;
                                tool_node =
                                    plan.boolean(BooleanKind::Union, tool_node, node, dec)?;
                            }
                        }
                        plan.boolean(BooleanKind::Subtract, target, tool_node, dec)?
                    }
                    _ => {
                        let (leaf, faces) = geometry::compile(&step.op, &step.id, dec)?;
                        plan.leaf(leaf, faces, dec)?
                    }
                };
                Ok(Some(node))
            })();
            span.finish(dec, result.as_ref().err(), None);
            if let Some(node) = result? {
                let mut faces = match &step.op {
                    Operation::Boolean { target, tool, .. } => step_faces[target]
                        .iter()
                        .chain(&step_faces[tool])
                        .copied()
                        .collect(),
                    Operation::ExtrudeThrough { target, .. } | Operation::Hole { target, .. } => {
                        step_faces[target].clone()
                    }
                    Operation::Pattern { body, .. }
                    | Operation::Fillet { body, .. }
                    | Operation::Mirror {
                        body,
                        keep_original: true,
                        ..
                    } => step_faces[body].clone(),
                    _ => Vec::new(),
                };
                faces.extend(first_face..plan.faces.len());
                faces.sort_unstable();
                faces.dedup();
                step_faces.insert(step.id.clone(), faces);
                steps.insert(step.id.clone(), node);
            }
        }
        plan.root = steps[&recipe.result];
        for id in &step_faces[&recipe.result] {
            plan.faces[*id].in_result = true;
        }
        plan.compact();
        Ok(plan)
    }

    /// Compiles a stock solid through its canonical sketch and sweep, with the
    /// same public origins as the direct analytic compiler.
    pub fn compile_canonical_stock(
        operation: &Operation,
        step: &catallax_recipe::StepId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let Some(expansion) = catallax_recipe::CanonicalStock::new(operation, step)? else {
            return Ok(None);
        };
        let mut plan = Self::compile_checked(expansion.recipe(), dec)?;
        for (index, face) in plan.faces.iter_mut().enumerate() {
            let origin = expansion.origin(&face.origin);
            if origin != &face.origin {
                plan.definition_roles
                    .insert(index, face.origin.role.clone());
                face.origin = origin.clone();
            }
        }
        Ok(Some(plan))
    }

    fn definition_role(&self, face: usize) -> &catallax_recipe::Role {
        self.definition_roles
            .get(&face)
            .unwrap_or(&self.faces[face].origin.role)
    }

    fn node_radius(&self, node: usize, dec: &mut Decider) -> Result<Rational, Refusal> {
        let mut bounds: Vec<Option<Rational>> = Vec::with_capacity(node + 1);
        for instruction in &self.instructions[..=node] {
            let bound = match *instruction {
                Instruction::Empty => Some(Rational::zero()),
                Instruction::Leaf(i) => self.leaves[i].display_radius(dec)?,
                Instruction::Boolean { kind, a, b } => match kind {
                    BooleanKind::Subtract => bounds[a].clone(),
                    BooleanKind::Union => bounds[a]
                        .as_ref()
                        .zip(bounds[b].as_ref())
                        .map(|(a, b)| a.clone().max(b.clone())),
                    BooleanKind::Intersect => match (&bounds[a], &bounds[b]) {
                        (Some(a), Some(b)) => Some(a.clone().min(b.clone())),
                        (Some(a), None) | (None, Some(a)) => Some(a.clone()),
                        _ => None,
                    },
                },
            };
            bounds.push(bound);
        }
        bounds[node].clone().ok_or_else(|| {
            Refusal::new(
                "field.unbounded-target",
                "hole target has no certified finite bound",
            )
        })
    }

    fn feature(
        &mut self,
        built: Option<(features::Feature, Vec<Face>)>,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        match built {
            Some((feature, faces)) => self.leaf(geometry::Leaf::Feature(feature), faces, dec),
            None => Ok(0),
        }
    }

    fn leaf(
        &mut self,
        leaf: geometry::Leaf,
        mut faces: Vec<Face>,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut found = None;
        for (id, previous) in self.leaves.iter().enumerate() {
            if previous.same_geometry(&leaf, dec)? {
                found = Some(id);
                break;
            }
        }
        let leaf_id = found.unwrap_or_else(|| {
            let id = self.leaves.len();
            self.leaves.push(leaf);
            id
        });
        for face in &mut faces {
            face.leaf = leaf_id;
        }
        self.faces.extend(faces);
        if let Some(id) = self
            .instructions
            .iter()
            .position(|n| matches!(n, Instruction::Leaf(i) if *i == leaf_id))
        {
            Ok(id)
        } else {
            self.instructions.push(Instruction::Leaf(leaf_id));
            Ok(self.instructions.len() - 1)
        }
    }

    /// Applies an exact rigid motion without changing the recipe or metric.
    pub fn transformed(&self, pose: &Pose, dec: &mut Decider) -> Result<Self, Refusal> {
        let pose = Pose::new(pose.rotation.clone(), pose.translation.clone(), dec)?;
        let map = Isometry::from(&pose);
        let mut result = self.clone();
        for leaf in &mut result.leaves {
            leaf.transform(&map)?;
        }
        for face in &mut result.faces {
            geometry::face_geometry(face, face.surface.transformed(&pose), dec)?;
        }
        Ok(result)
    }

    fn copy(
        &mut self,
        source: usize,
        source_faces: &[usize],
        map: &Isometry,
        step: &catallax_recipe::StepId,
        instance: u32,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let end = source + 1;
        let mut live = vec![false; end];
        live[source] = true;
        for i in (0..end).rev() {
            if live[i]
                && let Instruction::Boolean { a, b, .. } = self.instructions[i]
            {
                live[a] = true;
                live[b] = true;
            }
        }
        let mut remap = vec![0; end];
        for i in 1..end {
            if !live[i] {
                continue;
            }
            remap[i] = match self.instructions[i] {
                Instruction::Empty => 0,
                Instruction::Boolean { kind, a, b } => {
                    self.boolean(kind, remap[a], remap[b], dec)?
                }
                Instruction::Leaf(id) => {
                    let mut leaf = self.leaves[id].clone();
                    leaf.transform(map)?;
                    let mut faces = source_faces
                        .iter()
                        .map(|i| &self.faces[*i])
                        .filter(|f| f.leaf == id)
                        .cloned()
                        .collect::<Vec<_>>();
                    for face in &mut faces {
                        geometry::face_geometry(face, map.surface(&face.surface), dec)?;
                        face.origin = catallax_features::copy_origin(step, instance, &face.origin);
                    }
                    self.leaf(leaf, faces, dec)?
                }
            };
        }
        Ok(remap[source])
    }

    fn boolean(
        &mut self,
        kind: BooleanKind,
        a: usize,
        b: usize,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        use BooleanKind::*;
        if kind == Intersect
            && a != b
            && let (Instruction::Leaf(i), Instruction::Leaf(j)) =
                (self.instructions[a], self.instructions[b])
            && self.leaves[i].separated_by_plane(&self.leaves[j], dec)?
        {
            return Ok(0);
        }
        Ok(match kind {
            Union if a == b || b == 0 => a,
            Union if a == 0 => b,
            Intersect if a == 0 || b == 0 => 0,
            Intersect if a == b => a,
            Subtract if a == 0 || a == b => 0,
            Subtract if b == 0 => a,
            _ => {
                self.instructions.push(Instruction::Boolean { kind, a, b });
                self.instructions.len() - 1
            }
        })
    }

    fn compact(&mut self) {
        let mut live = vec![false; self.instructions.len()];
        live[0] = true;
        live[self.root] = true;
        for i in (0..self.instructions.len()).rev() {
            if live[i]
                && let Instruction::Boolean { a, b, .. } = self.instructions[i]
            {
                live[a] = true;
                live[b] = true;
            }
        }
        let mut remap = vec![0; live.len()];
        let mut tape = Vec::new();
        self.leaf_registers = vec![None; self.leaves.len()];
        for (i, instruction) in self.instructions.iter().copied().enumerate() {
            if !live[i] {
                continue;
            }
            remap[i] = tape.len();
            let instruction = match instruction {
                Instruction::Boolean { kind, a, b } => Instruction::Boolean {
                    kind,
                    a: remap[a],
                    b: remap[b],
                },
                Instruction::Leaf(leaf) => {
                    self.leaf_registers[leaf] = Some(tape.len());
                    instruction
                }
                Instruction::Empty => instruction,
            };
            tape.push(instruction);
        }
        self.root = remap[self.root];
        self.instructions = tape;
    }

    /// Whether expression identities or separating primitive halfspaces proved emptiness.
    /// `false` makes no non-emptiness claim for general Boolean expressions.
    pub fn is_proven_empty(&self) -> bool {
        self.root == 0
    }

    /// All source faces, including candidates later trimmed away. Call
    /// [`Self::face_query`] to establish membership in the selected result.
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    /// Resolves a stable recipe name without depending on tape indices.
    /// When a support is replaced without changing its origin, the result's
    /// candidate takes precedence over the retained construction-history entry.
    pub fn select(&self, origin: &Origin) -> Option<FaceId> {
        self.faces
            .iter()
            .position(|f| f.in_result && &f.origin == origin)
            .or_else(|| self.faces.iter().position(|f| &f.origin == origin))
            .map(FaceId)
    }

    /// Tentative scalar evaluation. Rejects NaN, infinities and finite inputs that
    /// overflow the approximate arithmetic. Signed zero is canonicalised to +0.
    pub fn fast(&self, point: [f64; 3], storage: &mut Workspace) -> Result<Approximate, Refusal> {
        if point.iter().any(|p| !p.is_finite()) {
            return Err(Refusal::new(
                "field.non-finite-point",
                "fast coordinates must be finite",
            ));
        }
        storage.fast.clear();
        for instruction in &self.instructions {
            let value = match *instruction {
                Instruction::Empty => f64::INFINITY,
                Instruction::Leaf(i) => {
                    let value = self.leaves[i].fast_distance(point);
                    if !value.is_finite() {
                        return Err(Refusal::new(
                            "field.fast-range",
                            "use the proven evaluator for this numeric range",
                        ));
                    }
                    value
                }
                Instruction::Boolean { kind, a, b } => match kind {
                    BooleanKind::Union => storage.fast[a].min(storage.fast[b]),
                    BooleanKind::Intersect => storage.fast[a].max(storage.fast[b]),
                    BooleanKind::Subtract => storage.fast[a].max(-storage.fast[b]),
                },
            };
            storage.fast.push(if value == 0.0 { 0.0 } else { value });
        }
        Ok(Approximate {
            value: storage.fast[self.root],
            true_distance_formula: matches!(self.instructions[self.root], Instruction::Leaf(i) if self.leaves[i].true_distance()),
        })
    }

    /// Authoritative evaluation at an exact rational point. Precision and local
    /// witness order are fixed, with no deadline or schedule-dependent refinement.
    /// Difficult Boolean zero sets return `Indeterminate`, never a false proof.
    pub fn proven(&self, point: &RPoint3, storage: &mut Workspace) -> Query {
        self.proven_checked(point, storage, &mut Decider::default())
            .unwrap_or_else(|_| unresolved_query())
    }

    /// Authoritative query under the shared policy, preserving numeric failures.
    pub fn proven_checked(
        &self,
        point: &RPoint3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Query, Refusal> {
        if self.root == 0 {
            let infinity = Interval::point(f64::INFINITY).unwrap();
            return Ok(Query {
                location: Location::Outside,
                value: Value::Residual(infinity),
                distance: DistanceBound {
                    millimetres: infinity,
                },
                direction: Direction::Undefined,
            });
        }
        let coordinates = [point.x.clone(), point.y.clone(), point.z.clone()];
        storage.values.clear();
        storage.locations.clear();
        for instruction in &self.instructions {
            let (value, location) = match *instruction {
                // Empty is folded out of every non-empty expression, so its
                // numeric register is unused and needs no infinite binary ball.
                Instruction::Empty => (Ball::integer(0), Location::Outside),
                Instruction::Leaf(i) => {
                    let location = self.leaves[i].location(&coordinates, dec)?;
                    (
                        self.leaves[i].proven_distance(point, location, dec)?,
                        location,
                    )
                }
                Instruction::Boolean { kind, a, b } => {
                    let value = match kind {
                        BooleanKind::Union => storage.values[a].min(&storage.values[b]),
                        BooleanKind::Intersect => storage.values[a].max(&storage.values[b]),
                        BooleanKind::Subtract => storage.values[a].max(&(-&storage.values[b])),
                    };
                    (
                        value,
                        combine(kind, storage.locations[a], storage.locations[b]),
                    )
                }
            };
            storage.values.push(value);
            storage.locations.push(location);
        }
        let primitive = matches!(self.instructions[self.root], Instruction::Leaf(i) if self.leaves[i].true_distance());
        let mut location = storage.locations[self.root];
        storage.active.clear();
        for (id, face) in self.faces.iter().enumerate() {
            if face.in_result
                && self.leaf_registers[face.leaf]
                    .is_some_and(|slot| storage.locations[slot] == Location::Boundary)
                && polynomial_sign(&face.polynomial, &coordinates, dec)? == Sign::Zero
                && self.leaves[face.leaf].on_face(self.definition_role(id), &coordinates, dec)?
                    != Some(false)
            {
                storage.active.push(id);
            }
        }
        if location == Location::Indeterminate {
            if let Some(proven) = self.planar_location(&coordinates, storage, dec)? {
                location = proven;
            } else if self.boundary_witness(&coordinates, storage, dec)? {
                location = Location::Boundary;
            }
        }
        let direction = if location == Location::Boundary {
            self.direction(&coordinates, storage, dec)?
        } else {
            Direction::Undefined
        };
        let enclosure = storage.values[self.root].to_interval();
        let absolute = unsigned(enclosure);
        let distance = if location == Location::Boundary {
            Interval::point(0.0).unwrap()
        } else if primitive {
            absolute
        } else {
            Interval::new(absolute.lo(), f64::INFINITY).unwrap()
        };
        Ok(Query {
            location,
            value: if primitive {
                Value::TrueDistance(enclosure)
            } else {
                Value::Residual(enclosure)
            },
            distance: DistanceBound {
                millimetres: distance,
            },
            direction,
        })
    }

    fn ray(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        storage: &mut Vec<Location>,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        storage.clear();
        for instruction in &self.instructions {
            let location = match *instruction {
                Instruction::Empty => Location::Outside,
                Instruction::Leaf(i) => self.leaves[i].ray(p, d, dec)?,
                Instruction::Boolean { kind, a, b } => combine(kind, storage[a], storage[b]),
            };
            storage.push(location);
        }
        Ok(storage[self.root])
    }

    fn planar_location(
        &self,
        p: &[Rational; 3],
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Option<Location>, Refusal> {
        storage.planes.clear();
        for (slot, instruction) in self.instructions.iter().enumerate() {
            if let Instruction::Leaf(i) = *instruction {
                match storage.locations[slot] {
                    Location::Boundary => {
                        if !self.leaves[i].planes_at(p, &mut storage.planes, dec)? {
                            return Ok(None);
                        }
                    }
                    Location::Indeterminate => return Ok(None),
                    _ => (),
                }
            }
        }
        // Most genuine corners have cheap transverse witnesses. Exhaustive
        // plane cells are needed only when that search cannot prove boundary.
        if storage.planes.len() > 1 && self.boundary_witness(p, storage, dec)? {
            return Ok(Some(Location::Boundary));
        }
        planar::classify(&storage.planes, dec, |n, dec| {
            self.ray(p, n, &mut storage.rays, dec)
        })
    }

    fn boundary_witness(
        &self,
        p: &[Rational; 3],
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut inside = false;
        let mut outside = false;
        for &id in &storage.active {
            let n = normal(&self.faces[id], p);
            if geometry::sign(&n.norm_squared(), dec)? == Sign::Zero {
                continue;
            }
            for n in [n.clone(), n.scaled(&Real::from_i64(-1))] {
                let loc = self.ray(p, &n, &mut storage.rays, dec)?;
                inside |= loc == Location::Inside;
                outside |= loc == Location::Outside;
            }
            if inside && outside {
                return Ok(true);
            }
        }
        // A fixed rational corpus also reaches box corners and transverse edges.
        for x in -1..=1 {
            for y in -1..=1 {
                for z in -1..=1 {
                    if x == 0 && y == 0 && z == 0 {
                        continue;
                    }
                    let n = RVector3::new(x.into(), y.into(), z.into()).to_real();
                    let loc = self.ray(p, &n, &mut storage.rays, dec)?;
                    inside |= loc == Location::Inside;
                    outside |= loc == Location::Outside;
                    if inside && outside {
                        return Ok(true);
                    }
                }
            }
        }
        Ok(false)
    }

    fn oriented_normal(
        &self,
        face: &Face,
        p: &[Rational; 3],
        storage: &mut Vec<Location>,
        dec: &mut Decider,
    ) -> Result<Option<Vector3>, Refusal> {
        let n = normal(face, p);
        if geometry::sign(&n.norm_squared(), dec)? == Sign::Zero {
            return Ok(None);
        }
        let plus = self.ray(p, &n, storage, dec)?;
        let minus = self.ray(p, &n.scaled(&Real::from_i64(-1)), storage, dec)?;
        Ok(match (plus, minus) {
            (Location::Outside, Location::Inside) => Some(n),
            (Location::Inside, Location::Outside) => Some(n.scaled(&Real::from_i64(-1))),
            _ => None,
        })
    }

    fn direction(
        &self,
        p: &[Rational; 3],
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Direction, Refusal> {
        let mut result: Option<Vector3> = None;
        for &id in &storage.active {
            let Some(n) = self.oriented_normal(&self.faces[id], p, &mut storage.rays, dec)? else {
                return Ok(Direction::Ambiguous);
            };
            if let Some(previous) = &result {
                if geometry::sign(&n.cross(previous).norm_squared(), dec)? != Sign::Zero
                    || geometry::sign(&n.dot(previous), dec)? != Sign::Positive
                {
                    return Ok(Direction::Ambiguous);
                }
            } else {
                result = Some(n);
            }
        }
        result.map_or(Ok(Direction::Ambiguous), |n| unit_direction(n, dec))
    }

    /// Tests a selected face's actual domain. Invalid handles return a refusal;
    /// buried faces do not become selectable boundary patches by name alone.
    pub fn face_query(
        &self,
        id: FaceId,
        point: &RPoint3,
        storage: &mut Workspace,
    ) -> Result<FaceQuery, Refusal> {
        self.face_query_checked(id, point, storage, &mut Decider::default())
    }

    /// Named-face query with an explicit shared proof policy.
    pub fn face_query_checked(
        &self,
        id: FaceId,
        point: &RPoint3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<FaceQuery, Refusal> {
        let origin = &self
            .faces
            .get(id.0)
            .ok_or_else(|| Refusal::new("field.face-id", "face handle is not in this plan"))?
            .origin;
        let mut chosen: Option<(usize, FaceQuery)> = None;
        let mut indeterminate = false;
        for (index, face) in self
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| &f.origin == origin)
        {
            let query = self.face_candidate_query(FaceId(index), point, storage, dec)?;
            match query.state {
                FaceState::OnFace => {
                    if let Some((old, previous)) = &mut chosen {
                        if !same_support(&self.faces[*old], face, dec)? {
                            previous.direction = Direction::Ambiguous;
                        }
                    } else {
                        chosen = Some((index, query));
                    }
                }
                FaceState::Indeterminate => indeterminate = true,
                FaceState::NotOnFace => (),
            }
        }
        if let Some((_, mut query)) = chosen {
            if indeterminate {
                query.direction = Direction::Ambiguous;
            }
            return Ok(query);
        }
        Ok(FaceQuery {
            state: if indeterminate {
                FaceState::Indeterminate
            } else {
                FaceState::NotOnFace
            },
            direction: if indeterminate {
                Direction::Ambiguous
            } else {
                Direction::Undefined
            },
        })
    }

    fn face_candidate_query(
        &self,
        id: FaceId,
        point: &RPoint3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<FaceQuery, Refusal> {
        let face = self
            .faces
            .get(id.0)
            .ok_or_else(|| Refusal::new("field.face-id", "face handle is not in this plan"))?;
        let p = [point.x.clone(), point.y.clone(), point.z.clone()];
        let absent = FaceQuery {
            state: FaceState::NotOnFace,
            direction: Direction::Undefined,
        };
        if !face.in_result || self.leaf_registers[face.leaf].is_none() {
            return Ok(absent);
        }
        if polynomial_sign(&face.polynomial, &p, dec)? != Sign::Zero {
            return Ok(absent);
        }
        match self.leaves[face.leaf].location(&p, dec)? {
            Location::Inside | Location::Outside => return Ok(absent),
            Location::Indeterminate => {
                return Ok(FaceQuery {
                    state: FaceState::Indeterminate,
                    direction: Direction::Ambiguous,
                });
            }
            Location::Boundary => (),
        }
        match self.leaves[face.leaf].on_face(self.definition_role(id.0), &p, dec)? {
            Some(false) => return Ok(absent),
            None => {
                return Ok(FaceQuery {
                    state: FaceState::Indeterminate,
                    direction: Direction::Ambiguous,
                });
            }
            Some(true) => (),
        }
        let query = self.proven_checked(point, storage, dec)?;
        if matches!(query.location, Location::Inside | Location::Outside) {
            return Ok(absent);
        }
        if !matches!(self.instructions[self.root], Instruction::Leaf(_)) {
            for &other in &storage.active {
                if !same_support(face, &self.faces[other], dec)? {
                    return Ok(FaceQuery {
                        state: FaceState::Indeterminate,
                        direction: Direction::Ambiguous,
                    });
                }
            }
        }
        if let Some(n) = self.oriented_normal(face, &p, &mut storage.rays, dec)? {
            return Ok(FaceQuery {
                state: FaceState::OnFace,
                direction: unit_direction(n, dec)?,
            });
        }
        if matches!(self.instructions[self.root], Instruction::Leaf(i) if i == face.leaf) {
            return Ok(FaceQuery {
                state: FaceState::OnFace,
                direction: Direction::Ambiguous,
            });
        }
        Ok(FaceQuery {
            state: FaceState::Indeterminate,
            direction: Direction::Ambiguous,
        })
    }

    /// Evaluates a list with reusable registers and output capacity. The caller
    /// may partition points between workers; there is no cross-point reduction.
    pub fn proven_batch<'a>(&self, points: &[RPoint3], storage: &'a mut Batch) -> &'a [Query] {
        storage.proven.clear();
        for point in points {
            storage
                .proven
                .push(self.proven(point, &mut storage.workspace));
        }
        &storage.proven
    }

    /// Tentative batch evaluation. On failure no partial results are retained.
    pub fn fast_batch<'a>(
        &self,
        points: &[[f64; 3]],
        storage: &'a mut Batch,
    ) -> Result<&'a [Approximate], Refusal> {
        storage.fast.clear();
        for &point in points {
            match self.fast(point, &mut storage.workspace) {
                Ok(value) => storage.fast.push(value),
                Err(error) => {
                    storage.fast.clear();
                    return Err(error);
                }
            }
        }
        Ok(&storage.fast)
    }
}

fn combine(kind: BooleanKind, a: Location, b: Location) -> Location {
    use Location::*;
    match kind {
        BooleanKind::Union => {
            if a == Inside || b == Inside {
                Inside
            } else if a == Outside && b == Outside {
                Outside
            } else {
                Indeterminate
            }
        }
        BooleanKind::Intersect => {
            if a == Outside || b == Outside {
                Outside
            } else if a == Inside && b == Inside {
                Inside
            } else {
                Indeterminate
            }
        }
        BooleanKind::Subtract => {
            if a == Outside || b == Inside {
                Outside
            } else if a == Inside && b == Outside {
                Inside
            } else {
                Indeterminate
            }
        }
    }
}
fn unsigned(value: Interval) -> Interval {
    Interval::new(
        if value.contains_zero() {
            0.0
        } else {
            value.lo().abs().min(value.hi().abs())
        },
        value.lo().abs().max(value.hi().abs()),
    )
    .unwrap()
}
fn polynomial_sign(poly: &RealPoly, p: &[Rational; 3], dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(poly) = poly.rational() {
        return Ok(poly.eval(p).sign());
    }
    geometry::sign(&poly.eval_real(&p.clone().map(Real::from)), dec)
}
fn normal(face: &Face, p: &[Rational; 3]) -> Vector3 {
    let values = p.clone().map(Real::from);
    Vector3::new(
        face.gradient[0].eval_real(&values),
        face.gradient[1].eval_real(&values),
        face.gradient[2].eval_real(&values),
    )
}
fn unit_direction(n: Vector3, dec: &mut Decider) -> Result<Direction, Refusal> {
    if let Some(n) = n.as_rational() {
        let norm = Ball::from_rational(&n.norm_squared(), PRECISION)
            .sqrt()
            .expect("nonzero rational norm");
        let components = [n.x, n.y, n.z].map(|r| {
            Ball::from_rational(&r, PRECISION)
                .checked_div(&norm)
                .expect("positive norm")
                .to_interval()
        });
        return Ok(Direction::Unique(components));
    }
    let length = n.norm_squared().sqrt();
    Ok(Direction::Unique([
        geometry::enclose(&(&n.x / &length), PRECISION, dec)?.to_interval(),
        geometry::enclose(&(&n.y / &length), PRECISION, dec)?.to_interval(),
        geometry::enclose(&(&n.z / &length), PRECISION, dec)?.to_interval(),
    ]))
}
fn same_support(a: &Face, b: &Face, dec: &mut Decider) -> Result<bool, Refusal> {
    a.support.equals(&b.support, dec)
}
fn unresolved_query() -> Query {
    Query {
        location: Location::Indeterminate,
        value: Value::Residual(Interval::entire()),
        distance: DistanceBound {
            millimetres: Interval::new(0.0, f64::INFINITY).unwrap(),
        },
        direction: Direction::Ambiguous,
    }
}
