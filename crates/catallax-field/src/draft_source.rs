//! Named planar source expressions for draft, independent of B-rep topology.
//!
//! The distance tape shares equal geometry. This temporary expression does not
//! share differently named creators: Boolean boundary provenance must survive
//! that optimisation, including copies of an entire Boolean expression.

use crate::{Instruction, Location, Plan, draft::Cell, geometry};
use catallax_features::draft::{
    arrangement::{Patch, Skin},
    cell::{PlanarCell, PlanePatch},
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, Operation, Origin, Recipe, StepId};
use catallax_shapes::{Isometry, Point3, Surface, Vector3};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone)]
struct Atom {
    evaluation: std::sync::Arc<[usize]>,
    root: usize,
    inverse: Isometry,
    skin: Skin,
}
#[derive(Clone)]
enum Node {
    Empty,
    Atom(Atom),
    Boolean {
        kind: BooleanKind,
        a: usize,
        b: usize,
    },
}
struct Source {
    nodes: Vec<Node>,
    root: usize,
    active: Vec<bool>,
    offsets: Vec<usize>,
}
#[derive(Default)]
struct State {
    minus: bool,
    plus: bool,
    patches: Vec<usize>,
}
#[derive(Default)]
struct Workspace {
    locations: Vec<Location>,
    states: Vec<State>,
}
fn failure(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-source-boundary", detail)
}
fn operation(kind: BooleanKind, a: bool, b: bool) -> bool {
    match kind {
        BooleanKind::Union => a || b,
        BooleanKind::Intersect => a && b,
        BooleanKind::Subtract => a && !b,
    }
}
fn live(nodes: &[Node], root: usize) -> Vec<bool> {
    let mut active = vec![false; root + 1];
    active[root] = true;
    for i in (0..=root).rev() {
        if active[i]
            && let Node::Boolean { a, b, .. } = nodes[i]
        {
            active[a] = true;
            active[b] = true;
        }
    }
    active
}
fn cell_skin(cell: &PlanarCell) -> Skin {
    Skin {
        vertices: cell.vertices.clone(),
        patches: cell
            .supports
            .iter()
            .zip(&cell.cycles)
            .filter(|(_, cycle)| !cycle.is_empty())
            .map(|(support, cycle)| Patch {
                plane: support.plane.clone(),
                loops: vec![cycle.clone()],
                origins: support.origins.clone(),
            })
            .collect(),
    }
}
fn leaf_location(
    leaf: &geometry::Leaf,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Location, Refusal> {
    if let geometry::Leaf::Feature(feature) = leaf {
        return feature.classify(point, dec);
    }
    if let geometry::Leaf::Draft(region) = leaf
        && let Cell::Polyhedral(cell) = &region.cell
    {
        return Ok(if !cell.contains(point, None, dec)? {
            Location::Outside
        } else if cell.on_boundary(point, None, dec)? {
            Location::Boundary
        } else {
            Location::Inside
        });
    }
    let values = [point.x.clone(), point.y.clone(), point.z.clone()];
    let mut boundary = false;
    for constraint in leaf.constraints() {
        match geometry::sign(&constraint.polynomial.eval_real(&values), dec)? {
            Sign::Positive => return Ok(Location::Outside),
            Sign::Zero => boundary = true,
            Sign::Negative => (),
        }
    }
    Ok(if boundary {
        Location::Boundary
    } else {
        Location::Inside
    })
}
impl Atom {
    fn location(
        &self,
        plan: &Plan,
        point: &Point3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let point = self.inverse.apply_point(point);
        storage
            .locations
            .resize(plan.instructions.len(), Location::Outside);
        for i in self.evaluation.iter() {
            storage.locations[*i] = match plan.instructions[*i] {
                Instruction::Empty => Location::Outside,
                Instruction::Leaf(leaf) => leaf_location(&plan.leaves[leaf], &point, dec)?,
                Instruction::Boolean { kind, a, b } => {
                    crate::combine(kind, storage.locations[a], storage.locations[b])
                }
            };
        }
        Ok(storage.locations[self.root])
    }
}
impl Source {
    fn boolean(&mut self, kind: BooleanKind, a: usize, b: usize) -> usize {
        let index = self.nodes.len();
        self.nodes.push(Node::Boolean { kind, a, b });
        index
    }
    fn copy(&mut self, root: usize, map: &Isometry, step: &StepId, instance: u32) -> usize {
        let active = live(&self.nodes, root);
        let mut remap = vec![0; root + 1];
        for i in 0..=root {
            if !active[i] {
                continue;
            }
            let node = match &self.nodes[i] {
                Node::Empty => continue,
                Node::Boolean { kind, a, b } => Node::Boolean {
                    kind: *kind,
                    a: remap[*a],
                    b: remap[*b],
                },
                Node::Atom(atom) => {
                    let mut atom = atom.clone();
                    atom.inverse = map.inverse().then(&atom.inverse);
                    for point in &mut atom.skin.vertices {
                        *point = map.apply_point(point);
                    }
                    for patch in &mut atom.skin.patches {
                        let Surface::Plane(plane) =
                            map.surface(&Surface::Plane(patch.plane.clone()))
                        else {
                            unreachable!("an isometry preserves planes")
                        };
                        patch.plane = plane;
                        if map.reversing() {
                            for cycle in &mut patch.loops {
                                cycle.reverse();
                            }
                        }
                        for origin in &mut patch.origins {
                            *origin = catallax_features::copy_origin(step, instance, origin);
                        }
                    }
                    Node::Atom(atom)
                }
            };
            remap[i] = self.nodes.len();
            self.nodes.push(node);
        }
        remap[root]
    }
    fn atom(
        &mut self,
        plan: &Plan,
        root: usize,
        faces: &[usize],
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if root == 0 {
            return Ok(0);
        }
        let mut active = vec![false; root + 1];
        active[root] = true;
        for i in (0..=root).rev() {
            if active[i]
                && let Instruction::Boolean { a, b, .. } = plan.instructions[i]
            {
                active[a] = true;
                active[b] = true;
            }
        }
        let evaluation: Vec<_> = (0..=root).filter(|i| active[*i]).collect();
        let mut skin = Skin {
            vertices: Vec::new(),
            patches: Vec::new(),
        };
        for register in &evaluation {
            let Instruction::Leaf(leaf) = plan.instructions[*register] else {
                continue;
            };
            let origins: Vec<_> = faces
                .iter()
                .filter_map(|f| (plan.faces[*f].leaf == leaf).then_some(&plan.faces[*f].origin))
                .collect();
            if origins.is_empty() {
                return Err(failure("a source leaf has no scoped face definitions"));
            }
            let part = match &plan.leaves[leaf] {
                geometry::Leaf::Lattice(_) => {
                    return Err(failure(
                        "field-only beams are not planar draft source faces",
                    ));
                }
                geometry::Leaf::Feature(feature) => {
                    feature.planar_skin(&origins, Some(bound), dec)?
                }
                geometry::Leaf::Draft(region) => {
                    let mut skin = match &region.cell {
                        Cell::Planar(cell) => cell_skin(cell),
                        Cell::Polyhedral(cell) => cell.skin(),
                        Cell::Round(_) => {
                            return Err(failure("a curved source needs curved arrangement cells"));
                        }
                    };
                    crate::draft::remap_skin_origins(&mut skin, &origins)?;
                    skin
                }
                geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_) => {
                    let witness = plan.leaves[leaf].draft_witness(dec)?;
                    let mut supports = Vec::new();
                    for face in faces
                        .iter()
                        .map(|i| &plan.faces[*i])
                        .filter(|f| f.leaf == leaf)
                    {
                        let Surface::Plane(mut plane) = face.surface.clone() else {
                            return Err(failure("a curved source needs curved arrangement cells"));
                        };
                        match geometry::sign(&witness.minus(&plane.origin).dot(&plane.normal), dec)?
                        {
                            Sign::Positive => {
                                plane.normal = plane.normal.scaled(&Real::from_i64(-1))
                            }
                            Sign::Negative => (),
                            Sign::Zero => {
                                return Err(failure(
                                    "a source carrier passes through its strict witness",
                                ));
                            }
                        }
                        supports.push(PlanePatch {
                            plane,
                            origins: vec![face.origin.clone()],
                        });
                    }
                    let cell = PlanarCell::new(supports, dec)?.ok_or_else(|| {
                        failure("source halfspaces have no full-dimensional cell")
                    })?;
                    cell_skin(&cell)
                }
            };
            append(&mut skin, &part);
        }
        let index = self.nodes.len();
        self.nodes.push(Node::Atom(Atom {
            evaluation: evaluation.into(),
            root,
            inverse: Isometry::identity(),
            skin,
        }));
        Ok(index)
    }
    fn joined(&mut self) -> Skin {
        self.active = live(&self.nodes, self.root);
        self.offsets.resize(self.nodes.len(), 0);
        let mut skin = Skin {
            vertices: Vec::new(),
            patches: Vec::new(),
        };
        for (i, node) in self.nodes.iter().enumerate().take(self.root + 1) {
            if self.active[i]
                && let Node::Atom(atom) = node
            {
                self.offsets[i] = skin.patches.len();
                append(&mut skin, &atom.skin);
            }
        }
        skin
    }
    fn classify(
        &self,
        plan: &Plan,
        point: &Point3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        storage.states.resize_with(self.nodes.len(), State::default);
        for i in 0..=self.root {
            if !self.active[i] {
                continue;
            }
            let inside = match &self.nodes[i] {
                Node::Empty => false,
                Node::Atom(atom) => match atom.location(plan, point, storage, dec)? {
                    Location::Inside => true,
                    Location::Outside => false,
                    _ => {
                        return Err(failure(
                            "a full-dimensional source cell lies on an unresolved boundary",
                        ));
                    }
                },
                Node::Boolean { kind, a, b } => {
                    operation(*kind, storage.states[*a].minus, storage.states[*b].minus)
                }
            };
            storage.states[i].minus = inside;
        }
        Ok(storage.states[self.root].minus)
    }
    fn boundary(
        &self,
        plan: &Plan,
        skin: &Skin,
        point: &Point3,
        normal: &Vector3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Vec<Origin>, Refusal> {
        storage.states.resize_with(self.nodes.len(), State::default);
        for i in 0..=self.root {
            if !self.active[i] {
                continue;
            }
            storage.states[i].patches.clear();
            match &self.nodes[i] {
                Node::Empty => {
                    storage.states[i].minus = false;
                    storage.states[i].plus = false;
                }
                Node::Atom(atom) => {
                    let location = atom.location(plan, point, storage, dec)?;
                    let state = &mut storage.states[i];
                    match location {
                        Location::Inside | Location::Outside => {
                            state.minus = location == Location::Inside;
                            state.plus = state.minus;
                        }
                        Location::Boundary => {
                            let mut sense = None;
                            for (patch_index, patch) in atom.skin.patches.iter().enumerate() {
                                if !atom.skin.contains_patch(patch, point, dec)? {
                                    continue;
                                }
                                let aligned =
                                    match geometry::sign(&patch.plane.normal.dot(normal), dec)? {
                                        Sign::Positive => true,
                                        Sign::Negative => false,
                                        Sign::Zero => {
                                            return Err(failure(
                                                "source facet has no transverse outward direction",
                                            ));
                                        }
                                    };
                                if sense.is_some_and(|previous| previous != aligned) {
                                    return Err(failure(
                                        "source facet has contradictory outward carriers",
                                    ));
                                }
                                sense = Some(aligned);
                                state.patches.push(self.offsets[i] + patch_index);
                            }
                            state.minus = sense.ok_or_else(|| {
                                failure("source boundary has no physical defining patch")
                            })?;
                            state.plus = !state.minus;
                        }
                        Location::Indeterminate => {
                            return Err(failure("source boundary classification is unresolved"));
                        }
                    }
                }
                Node::Boolean { kind, a, b } => {
                    let (previous, current) = storage.states.split_at_mut(i);
                    let state = &mut current[0];
                    state.minus = operation(*kind, previous[*a].minus, previous[*b].minus);
                    state.plus = operation(*kind, previous[*a].plus, previous[*b].plus);
                    if state.minus != state.plus {
                        merge(
                            &previous[*a].patches,
                            &previous[*b].patches,
                            &mut state.patches,
                        );
                    }
                }
            }
        }
        let result = &storage.states[self.root];
        if !result.minus || result.plus {
            return Err(failure(
                "classified source boundary and its outward germ disagree",
            ));
        }
        Ok(result
            .patches
            .iter()
            .flat_map(|i| skin.patches[*i].origins.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect())
    }
}
fn append(target: &mut Skin, source: &Skin) {
    let offset = target.vertices.len();
    target.vertices.extend(source.vertices.iter().cloned());
    target.patches.extend(source.patches.iter().map(|patch| {
        Patch {
            plane: patch.plane.clone(),
            loops: patch
                .loops
                .iter()
                .map(|cycle| cycle.iter().map(|i| offset + i).collect())
                .collect(),
            origins: patch.origins.clone(),
        }
    }));
}
fn merge(a: &[usize], b: &[usize], output: &mut Vec<usize>) {
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        let value = match (a.get(i), b.get(j)) {
            (Some(a), Some(b)) => (*a).min(*b),
            (Some(a), None) => *a,
            (None, Some(b)) => *b,
            (None, None) => unreachable!("one sequence remains"),
        };
        output.push(value);
        if a.get(i) == Some(&value) {
            i += 1;
        }
        if b.get(j) == Some(&value) {
            j += 1;
        }
    }
}
impl Plan {
    pub(crate) fn draft_source_skin(
        &self,
        recipe: &Recipe,
        body: &StepId,
        registers: &BTreeMap<StepId, usize>,
        faces: &BTreeMap<StepId, Vec<usize>>,
        dec: &mut Decider,
    ) -> Result<Skin, Refusal> {
        if faces[body]
            .iter()
            .any(|f| !matches!(self.faces[*f].surface, Surface::Plane(_)))
        {
            return Err(failure(
                "mixed curved input needs oriented curved arrangement cells",
            ));
        }
        // Historical finite leaves are included deliberately. A copied source
        // expression is evaluated in its original frame; its finite clipping
        // domain must also contain the original, possibly remote, material.
        let mut bound = Rational::one();
        for leaf in &self.leaves {
            if let Some(radius) = leaf.display_radius(dec)? {
                bound = bound.max(radius);
            }
        }
        let mut needed = BTreeSet::from([body.clone()]);
        for step in recipe.steps.iter().rev() {
            if !needed.contains(&step.id) {
                continue;
            }
            match &step.op {
                Operation::Boolean { target, tool, .. } => {
                    needed.insert(target.clone());
                    needed.insert(tool.clone());
                }
                Operation::ExtrudeThrough { target, .. } => {
                    needed.insert(target.clone());
                }
                Operation::Pattern { body, .. } | Operation::Mirror { body, .. } => {
                    needed.insert(body.clone());
                }
                _ => (),
            }
        }
        let mut source = Source {
            nodes: vec![Node::Empty],
            root: 0,
            active: Vec::new(),
            offsets: Vec::new(),
        };
        let mut roots = BTreeMap::new();
        for step in &recipe.steps {
            if !needed.contains(&step.id) {
                continue;
            }
            let root = match &step.op {
                Operation::Boolean { kind, target, tool } => {
                    source.boolean(*kind, roots[target], roots[tool])
                }
                Operation::ExtrudeThrough { target, .. } => {
                    let tool_faces: Vec<_> = faces[&step.id]
                        .iter()
                        .copied()
                        .filter(|f| self.faces[*f].origin.step == step.id)
                        .collect();
                    let tool = if let Some(face) = tool_faces.first() {
                        let leaf = self.faces[*face].leaf;
                        let root = self
                            .instructions
                            .iter()
                            .position(|v| matches!(v, Instruction::Leaf(i) if *i == leaf))
                            .expect("registered faces have a leaf instruction");
                        source.atom(self, root, &tool_faces, &bound, dec)?
                    } else {
                        0
                    };
                    source.boolean(BooleanKind::Subtract, roots[target], tool)
                }
                Operation::Pattern { body, layout } => {
                    let seed = roots[body];
                    let mut root = seed;
                    for instance in 1..layout.count() {
                        let map = catallax_features::pattern_placement(layout, instance, dec)?;
                        let copy = source.copy(seed, &map, &step.id, instance);
                        root = source.boolean(BooleanKind::Union, root, copy);
                    }
                    root
                }
                Operation::Mirror {
                    body,
                    plane,
                    keep_original,
                } => {
                    let map = catallax_features::mirror_placement(plane)?;
                    let copy = source.copy(roots[body], &map, &step.id, 1);
                    if *keep_original {
                        source.boolean(BooleanKind::Union, roots[body], copy)
                    } else {
                        copy
                    }
                }
                _ => source.atom(self, registers[&step.id], &faces[&step.id], &bound, dec)?,
            };
            roots.insert(step.id.clone(), root);
        }
        source.root = roots[body];
        let joined = source.joined();
        let mut classifier = Workspace::default();
        let mut resolver = Workspace::default();
        let components = joined.partition(
            |point, dec| source.classify(self, point, &mut classifier, dec),
            |point, normal, dec| source.boundary(self, &joined, point, normal, &mut resolver, dec),
            dec,
        )?;
        let mut skin = Skin {
            vertices: Vec::new(),
            patches: Vec::new(),
        };
        for component in components {
            append(&mut skin, &component.skin());
        }
        Ok(skin)
    }
}
