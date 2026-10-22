//! Recipe-defined analytic fields in safe Rust.
//!
//! [`Plan`] owns a fixed-order evaluation tape and a separate named [`Face`]
//! registry. [`Plan::fast`] is tentative display work. [`Plan::proven`] uses
//! exact support and profile predicates and outward-rounded binary balls, with
//! no platform math library. Sweep and Boolean residuals are not true distance.
//! See `docs/FIELD.md` for the distance proof and unresolved-boundary contract.

#![forbid(unsafe_code)]

mod features;
mod geometry;
mod metric;
mod profile;
mod samples;
mod support;

use std::collections::BTreeMap;

use catallax_numbers::{Ball, Interval, MPoly, Rational, Refusal, Sign};
use catallax_recipe::{BooleanKind, Operation, Origin, Recipe};
use catallax_shapes::{RPoint3, RVector3, Surface};

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
    polynomial: MPoly,
    support: MPoly,
    gradient: [MPoly; 3],
    leaf: usize,
}

/// Handle into one plan's immutable face registry. Not a persistent recipe name.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct FaceId(usize);

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
}

impl Plan {
    /// Compiles stock solids, sketch sweeps and regularised booleans.
    /// Unsupported operations and impossible dimensions report a stable refusal.
    /// Sketch validation is shared with features; no B-rep builder is called.
    pub fn compile(recipe: &Recipe) -> Result<Self, Refusal> {
        recipe
            .check()
            .map_err(|e| Refusal::new("field.recipe", e.to_string()))?;
        let mut plan = Self {
            instructions: vec![Instruction::Empty],
            leaves: vec![],
            faces: vec![],
            root: 0,
            leaf_registers: vec![],
        };
        let mut steps = BTreeMap::new();
        let mut profiles = BTreeMap::new();
        for step in &recipe.steps {
            let node = match &step.op {
                Operation::Sketch { plane, loops } => {
                    profiles.insert(step.id.clone(), profile::Profile::new(plane, loops)?);
                    continue;
                }
                Operation::Boolean { kind, target, tool } => {
                    plan.boolean(*kind, steps[target], steps[tool])
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
                    )?;
                    plan.feature(built)
                }
                Operation::ExtrudeThrough {
                    target,
                    profile,
                    direction,
                } => {
                    let built = features::Feature::extrude(
                        &profiles[profile],
                        &step.id,
                        None,
                        direction.as_ref(),
                        steps[target] == 0,
                    )?;
                    let tool = plan.feature(built);
                    plan.boolean(BooleanKind::Subtract, steps[target], tool)
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
                    )?;
                    plan.feature(built)
                }
                _ => {
                    let (leaf, faces) = geometry::compile(&step.op, &step.id)?;
                    plan.leaf(leaf, faces)
                }
            };
            steps.insert(step.id.clone(), node);
        }
        plan.root = steps[&recipe.result];
        plan.compact();
        Ok(plan)
    }

    fn feature(&mut self, built: Option<(features::Feature, Vec<Face>)>) -> usize {
        built.map_or(0, |(feature, faces)| {
            self.leaf(geometry::Leaf::Feature(feature), faces)
        })
    }

    fn leaf(&mut self, leaf: geometry::Leaf, mut faces: Vec<Face>) -> usize {
        let leaf_id = self
            .leaves
            .iter()
            .position(|l| l.same_geometry(&leaf))
            .unwrap_or_else(|| {
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
            id
        } else {
            self.instructions.push(Instruction::Leaf(leaf_id));
            self.instructions.len() - 1
        }
    }

    fn boolean(&mut self, kind: BooleanKind, a: usize, b: usize) -> usize {
        use BooleanKind::*;
        match kind {
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
        }
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

    /// Whether exact expression identities proved that the result is empty.
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
    pub fn select(&self, origin: &Origin) -> Option<FaceId> {
        self.faces
            .iter()
            .position(|f| &f.origin == origin)
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
        if self.root == 0 {
            let infinity = Interval::point(f64::INFINITY).unwrap();
            return Query {
                location: Location::Outside,
                value: Value::Residual(infinity),
                distance: DistanceBound {
                    millimetres: infinity,
                },
                direction: Direction::Undefined,
            };
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
                    let location = self.leaves[i].location(&coordinates);
                    (self.leaves[i].proven_distance(point, location), location)
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
            if self.leaf_registers[face.leaf]
                .is_some_and(|slot| storage.locations[slot] == Location::Boundary)
                && face.polynomial.eval(&coordinates).is_zero()
                && self.leaves[face.leaf].on_face(&face.origin.role, &coordinates) != Some(false)
            {
                storage.active.push(id);
            }
        }
        if location == Location::Indeterminate && self.boundary_witness(&coordinates, storage) {
            location = Location::Boundary;
        }
        let direction = if location == Location::Boundary {
            self.direction(&coordinates, storage)
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
        Query {
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
        }
    }

    fn ray(&self, p: &[Rational; 3], d: &RVector3, storage: &mut Vec<Location>) -> Location {
        storage.clear();
        let direction = [d.x.clone(), d.y.clone(), d.z.clone()];
        for instruction in &self.instructions {
            let location = match *instruction {
                Instruction::Empty => Location::Outside,
                Instruction::Leaf(i) => self.leaves[i].ray(p, &direction),
                Instruction::Boolean { kind, a, b } => combine(kind, storage[a], storage[b]),
            };
            storage.push(location);
        }
        storage[self.root]
    }

    fn boundary_witness(&self, p: &[Rational; 3], storage: &mut Workspace) -> bool {
        let mut inside = false;
        let mut outside = false;
        for &id in &storage.active {
            let n = normal(&self.faces[id], p);
            if n.is_zero() {
                continue;
            }
            for n in [n.clone(), n.scaled(&Rational::from_i64(-1))] {
                let loc = self.ray(p, &n, &mut storage.rays);
                inside |= loc == Location::Inside;
                outside |= loc == Location::Outside;
            }
            if inside && outside {
                return true;
            }
        }
        // A fixed rational corpus also reaches box corners and transverse edges.
        for x in -1..=1 {
            for y in -1..=1 {
                for z in -1..=1 {
                    if x == 0 && y == 0 && z == 0 {
                        continue;
                    }
                    let n = RVector3::new(x.into(), y.into(), z.into());
                    let loc = self.ray(p, &n, &mut storage.rays);
                    inside |= loc == Location::Inside;
                    outside |= loc == Location::Outside;
                    if inside && outside {
                        return true;
                    }
                }
            }
        }
        false
    }

    fn oriented_normal(
        &self,
        face: &Face,
        p: &[Rational; 3],
        storage: &mut Vec<Location>,
    ) -> Option<RVector3> {
        let n = normal(face, p);
        if n.is_zero() {
            return None;
        }
        let plus = self.ray(p, &n, storage);
        let minus = self.ray(p, &n.scaled(&Rational::from_i64(-1)), storage);
        match (plus, minus) {
            (Location::Outside, Location::Inside) => Some(n),
            (Location::Inside, Location::Outside) => Some(n.scaled(&Rational::from_i64(-1))),
            _ => None,
        }
    }

    fn direction(&self, p: &[Rational; 3], storage: &mut Workspace) -> Direction {
        let mut result: Option<RVector3> = None;
        for &id in &storage.active {
            let Some(n) = self.oriented_normal(&self.faces[id], p, &mut storage.rays) else {
                return Direction::Ambiguous;
            };
            if let Some(previous) = &result {
                if !n.cross(previous).is_zero() || n.dot(previous).sign() != Sign::Positive {
                    return Direction::Ambiguous;
                }
            } else {
                result = Some(n);
            }
        }
        result.map_or(Direction::Ambiguous, unit_direction)
    }

    /// Tests a selected face's actual domain. Invalid handles return a refusal;
    /// buried faces do not become selectable boundary patches by name alone.
    pub fn face_query(
        &self,
        id: FaceId,
        point: &RPoint3,
        storage: &mut Workspace,
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
        if self.leaf_registers[face.leaf].is_none() {
            return Ok(absent);
        }
        if !face.polynomial.eval(&p).is_zero() {
            return Ok(absent);
        }
        match self.leaves[face.leaf].location(&p) {
            Location::Inside | Location::Outside => return Ok(absent),
            Location::Indeterminate => {
                return Ok(FaceQuery {
                    state: FaceState::Indeterminate,
                    direction: Direction::Ambiguous,
                });
            }
            Location::Boundary => (),
        }
        match self.leaves[face.leaf].on_face(&face.origin.role, &p) {
            Some(false) => return Ok(absent),
            None => {
                return Ok(FaceQuery {
                    state: FaceState::Indeterminate,
                    direction: Direction::Ambiguous,
                });
            }
            Some(true) => (),
        }
        let query = self.proven(point, storage);
        if matches!(query.location, Location::Inside | Location::Outside) {
            return Ok(absent);
        }
        if !matches!(self.instructions[self.root], Instruction::Leaf(_))
            && storage
                .active
                .iter()
                .any(|&other| !same_support(face, &self.faces[other]))
        {
            return Ok(FaceQuery {
                state: FaceState::Indeterminate,
                direction: Direction::Ambiguous,
            });
        }
        if let Some(n) = self.oriented_normal(face, &p, &mut storage.rays) {
            return Ok(FaceQuery {
                state: FaceState::OnFace,
                direction: unit_direction(n),
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
fn normal(face: &Face, p: &[Rational; 3]) -> RVector3 {
    RVector3::new(
        face.gradient[0].eval(p),
        face.gradient[1].eval(p),
        face.gradient[2].eval(p),
    )
}
fn unit_direction(n: RVector3) -> Direction {
    let norm = Ball::from_rational(&n.norm_squared(), PRECISION)
        .sqrt()
        .expect("nonzero rational norm");
    let components = [n.x, n.y, n.z].map(|r| {
        Ball::from_rational(&r, PRECISION)
            .checked_div(&norm)
            .expect("positive norm")
            .to_interval()
    });
    Direction::Unique(components)
}

fn same_support(a: &Face, b: &Face) -> bool {
    a.support == b.support
}
