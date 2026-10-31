//! Named planar source expressions for draft, independent of B-rep topology.
//!
//! The distance tape shares equal geometry. This temporary expression does not
//! share differently named creators: Boolean boundary provenance must survive
//! that optimisation, including copies of an entire Boolean expression.

use crate::{
    Instruction, Location, Plan,
    draft::{Cell, PreservedBore},
    geometry,
};
use catallax_features::draft::{
    DraftDefinition,
    arrangement::{Patch, Skin},
    bored,
    cell::{PlanarCell, PlanePatch},
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, FaceRef, Operation, Origin, Recipe, StepId};
use catallax_shapes::{Cylinder, Isometry, Plane, Point3, Surface, Vector3};
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
pub(crate) fn failure(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-source-boundary", detail)
}
/// Keeps this evaluator's code while a numeric refusal keeps its own.
pub(crate) fn source_refusal(refusal: Refusal) -> Refusal {
    if refusal.code == "draft.unsupported-arrangement" {
        failure(&refusal.detail)
    } else {
        refusal
    }
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
                geometry::Leaf::Lattice(_) | geometry::Leaf::Gyroid(_) => {
                    return Err(failure(
                        "field-only surfaces are not planar draft source faces",
                    ));
                }
                geometry::Leaf::Feature(feature) => {
                    feature.planar_skin(&origins, Some(bound), dec)?
                }
                geometry::Leaf::Draft(region) => {
                    let mut skin = match &region.cell {
                        Cell::Planar(cell) => cell_skin(cell),
                        Cell::Polyhedral(cell) => cell.skin(),
                        Cell::Round(_) | Cell::Annular(_) => {
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
    /// Republishes the unchanged wall as faces of the drafted result.
    pub(crate) fn preserve_faces(&mut self, faces: &[usize]) {
        for face in faces {
            self.faces.push(self.faces[*face].clone());
        }
    }

    pub(crate) fn register_location(
        &self,
        register: usize,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        Ok(match self.instructions[register] {
            Instruction::Empty => Location::Outside,
            Instruction::Leaf(leaf) => leaf_location(&self.leaves[leaf], point, dec)?,
            Instruction::Boolean { kind, a, b } => crate::combine(
                kind,
                self.register_location(a, point, dec)?,
                self.register_location(b, point, dec)?,
            ),
        })
    }

    /// The tool's stacked bands, each bounded only by end planes and one wall.
    fn tool_bands(
        &self,
        register: usize,
        cylinder: &Cylinder,
        planes: &mut Vec<Plane>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let leaf = match self.instructions[register] {
            Instruction::Leaf(leaf) => leaf,
            Instruction::Boolean {
                kind: BooleanKind::Union,
                a,
                b,
            } => {
                self.tool_bands(a, cylinder, planes, dec)?;
                return self.tool_bands(b, cylinder, planes, dec);
            }
            _ => return Err(failure("a preserved bore needs a union of tool bands")),
        };
        for constraint in self.leaves[leaf].constraints() {
            let axial = match constraint.surface() {
                Surface::Plane(plane) => {
                    planes.push(plane.clone());
                    plane.normal.cross(&cylinder.axis)
                }
                Surface::Cylinder(band) => {
                    if geometry::sign(&(&band.radius_sq - &cylinder.radius_sq), dec)? != Sign::Zero
                        || geometry::sign(
                            &band
                                .origin
                                .minus(&cylinder.origin)
                                .cross(&cylinder.axis)
                                .norm_squared(),
                            dec,
                        )? != Sign::Zero
                    {
                        return Err(failure("preserved bore bands need one cylinder"));
                    }
                    band.axis.cross(&cylinder.axis)
                }
                _ => return Err(failure("a preserved bore tool has an unsupported face")),
            };
            // A parallel flat or an oblique end would leave a noncircular mouth.
            if geometry::sign(&axial.norm_squared(), dec)? != Sign::Zero {
                return Err(failure("a preserved bore tool has an unsupported face"));
            }
        }
        Ok(())
    }

    /// Proves the tool holds the axis over every interval between its end planes
    /// inside the caps. Bands share one wall, so the whole section is removed.
    fn check_through(
        &self,
        tool: usize,
        cylinder: &Cylinder,
        caps: [&Plane; 2],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let mut planes = Vec::new();
        self.tool_bands(tool, cylinder, &mut planes, dec)?;
        let lower = bored::axis_point(cylinder, caps[0]);
        let upper = bored::axis_point(cylinder, caps[1]);
        let height = |p: &Point3| p.minus(&lower).dot(&cylinder.axis);
        let mut stops = vec![lower.clone(), upper];
        for plane in &planes {
            let stop = bored::axis_point(cylinder, plane);
            let mut index = 0;
            while index < stops.len()
                && geometry::sign(&(height(&stops[index]) - height(&stop)), dec)? == Sign::Negative
            {
                index += 1;
            }
            if index == 0
                || index == stops.len()
                || geometry::sign(&(height(&stops[index]) - height(&stop)), dec)? == Sign::Zero
            {
                continue;
            }
            stops.insert(index, stop);
        }
        for pair in stops.windows(2) {
            let half = pair[1]
                .minus(&pair[0])
                .scaled(&(Real::one() / Real::from_i64(2)));
            if self.register_location(tool, &pair[0].plus(&half), dec)? != Location::Inside {
                return Err(failure("a preserved bore must pass through both caps"));
            }
        }
        Ok(())
    }

    /// Recognises a recipe-defined subtraction of one complete cylinder from a
    /// planar target. The proof reads recipe faces and exact leaves, never a
    /// B-rep. Any other curved source returns `None` and keeps its refusal.
    pub(crate) fn preserved_bore(
        &self,
        recipe: &Recipe,
        body: &StepId,
        refs: &[FaceRef],
        definition: &DraftDefinition,
        (registers, faces): (&BTreeMap<StepId, usize>, &BTreeMap<StepId, Vec<usize>>),
        dec: &mut Decider,
    ) -> Result<Option<(StepId, Skin, PreservedBore)>, Refusal> {
        let target = match recipe.step(body).map(|s| &s.op) {
            Some(Operation::Boolean {
                kind: BooleanKind::Subtract,
                target,
                ..
            })
            | Some(Operation::Hole { target, .. })
            | Some(Operation::ExtrudeThrough { target, .. }) => target,
            _ => return Ok(None),
        };
        let plane = |f: &usize| matches!(self.faces[*f].surface, Surface::Plane(_));
        if faces[body].iter().all(plane) || !faces[target].iter().all(plane) {
            return Ok(None);
        }
        let Instruction::Boolean {
            kind: BooleanKind::Subtract,
            a,
            b,
        } = self.instructions[registers[body]]
        else {
            return Ok(None);
        };
        if a != registers[target] {
            return Ok(None);
        }
        let added: Vec<_> = faces[body]
            .iter()
            .copied()
            .filter(|f| faces[target].binary_search(f).is_err())
            .collect();
        let mut walls = Vec::new();
        let mut cylinder = None;
        for f in &added {
            if let Surface::Cylinder(wall) = &self.faces[*f].surface {
                walls.push(*f);
                cylinder.get_or_insert_with(|| wall.clone());
            }
        }
        let Some(cylinder) = cylinder else {
            return Err(failure("a preserved bore needs one cylindrical wall"));
        };
        if refs
            .iter()
            .any(|r| added.iter().any(|f| r.face.matches(&self.faces[*f].origin)))
        {
            return Err(failure("a selected bore needs a curved arrangement"));
        }
        let skin = self.draft_source_skin(recipe, target, registers, faces, dec)?;
        if !skin.is_convex(dec)? {
            return Err(failure("a preserved bore needs a convex planar envelope"));
        }
        bored::check_axial(&cylinder, definition, dec).map_err(source_refusal)?;
        let caps = bored::certify_clear(
            skin.patches.iter().map(|patch| &patch.plane),
            &cylinder,
            dec,
        )
        .map_err(source_refusal)?;
        let selected = skin.resolve_faces(refs, dec)?;
        if caps.iter().any(|cap| selected.contains(cap)) {
            return Err(failure("a selected bore cap moves the cap-hole trim"));
        }
        self.check_through(b, &cylinder, caps.map(|cap| &skin.patches[cap].plane), dec)?;
        Ok(Some((
            target.clone(),
            skin,
            PreservedBore {
                register: b,
                walls,
                cylinder,
            },
        )))
    }

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
