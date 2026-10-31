//! Named planar source expressions for draft, independent of B-rep topology.
//!
//! The distance tape shares equal geometry. This temporary expression does not
//! share differently named creators: Boolean boundary provenance must survive
//! that optimisation, including copies of an entire Boolean expression.

use crate::{
    Instruction, Location, Plan,
    draft::{Cell, Preserved, PreservedBore, PreservedBores, PreservedBoss, PreservedHead},
    geometry,
};
use catallax_features::draft::{
    DraftDefinition,
    arrangement::{Patch, Skin},
    bored,
    bossed::{self, Boss},
    cell::{PlanarCell, PlanePatch},
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, FaceRef, Operation, Origin, Recipe, StepId};
use catallax_shapes::{Cone, Cylinder, Isometry, Plane, Point3, Surface, Vector3};
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
/// The tape registers and faces of a recipe's compiled steps.
type Steps<'a> = (
    &'a BTreeMap<StepId, usize>,
    &'a BTreeMap<StepId, Vec<usize>>,
);
/// One subtraction's tool register, the faces it added, and whether a recipe hole
/// defined it.
type ToolChain = (usize, Vec<usize>, bool);
/// One bore inside a subtraction's tool, and the end planes of its bands.
struct ToolBore {
    cylinder: Cylinder,
    planes: Vec<Plane>,
    /// The leaves carrying this cylinder.
    leaves: Vec<usize>,
    /// The cone opening a countersunk bore.
    head: Option<ToolHead>,
    /// Wider bands of the same tool, each of which must stay beyond the head's cap.
    entries: Vec<ToolBore>,
}
/// One conical band of a subtraction's tool.
struct ToolHead {
    cone: Cone,
    planes: Vec<Plane>,
    leaf: usize,
    /// Whether a recipe hole defines this band and its wall from one station. Their
    /// shared ring then needs no exact comparison with that station's plane.
    stationed: bool,
}
/// One bore with its subtraction's tool register and added faces.
type ChainBore<'a> = (usize, &'a [usize], ToolBore);
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
/// Whether two bands lie on one cylinder, whichever way their axes point.
fn same_wall(a: &Cylinder, b: &Cylinder, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        geometry::sign(&(&a.radius_sq - &b.radius_sq), dec)? == Sign::Zero
            && geometry::sign(&a.axis.cross(&b.axis).norm_squared(), dec)? == Sign::Zero
            && geometry::sign(
                &b.origin.minus(&a.origin).cross(&a.axis).norm_squared(),
                dec,
            )? == Sign::Zero,
    )
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
    /// Republishes unchanged faces as faces of the drafted result.
    fn preserve_faces(&mut self, faces: &[usize]) {
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

    /// Groups a tool's walls by bore. Each band is bounded only by end planes
    /// and its cylinders or its one cone, and stacked bands of one bore share a
    /// cylinder.
    fn tool_bores(
        &self,
        register: usize,
        bores: &mut Vec<ToolBore>,
        heads: &mut Vec<ToolHead>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let leaf = match self.instructions[register] {
            Instruction::Leaf(leaf) => leaf,
            Instruction::Boolean {
                kind: BooleanKind::Union,
                a,
                b,
            } => {
                self.tool_bores(a, bores, heads, dec)?;
                return self.tool_bores(b, bores, heads, dec);
            }
            _ => return Err(failure("a preserved bore needs a union of tool bands")),
        };
        let constraints = self.leaves[leaf].constraints();
        let mut cones = constraints.iter().filter_map(|c| match c.surface() {
            Surface::Cone(cone) => Some(cone),
            _ => None,
        });
        if let Some(cone) = cones.next() {
            let mut planes = Vec::new();
            for constraint in constraints {
                match constraint.surface() {
                    Surface::Plane(plane)
                        if geometry::sign(&plane.normal.cross(&cone.axis).norm_squared(), dec)?
                            == Sign::Zero =>
                    {
                        planes.push(plane.clone());
                    }
                    Surface::Cone(_) if cones.next().is_none() => (),
                    _ => return Err(failure("a preserved bore tool has an unsupported face")),
                }
            }
            heads.push(ToolHead {
                cone: cone.clone(),
                planes,
                leaf,
                stationed: false,
            });
            return Ok(());
        }
        let mut owners = Vec::new();
        for constraint in constraints {
            let Surface::Cylinder(band) = constraint.surface() else {
                continue;
            };
            let mut owner = None;
            for (index, bore) in bores.iter().enumerate() {
                if same_wall(&bore.cylinder, band, dec)? {
                    owner = Some(index);
                    break;
                }
            }
            owners.push(owner.unwrap_or_else(|| {
                bores.push(ToolBore {
                    cylinder: band.clone(),
                    planes: Vec::new(),
                    leaves: Vec::new(),
                    head: None,
                    entries: Vec::new(),
                });
                bores.len() - 1
            }));
            bores[*owners.last().expect("just pushed")]
                .leaves
                .push(leaf);
        }
        if owners.is_empty() {
            return Err(failure("a preserved bore tool has an unsupported face"));
        }
        for constraint in constraints {
            match constraint.surface() {
                Surface::Plane(plane) => {
                    for owner in &owners {
                        let bore = &mut bores[*owner];
                        // A parallel flat or an oblique end would leave a noncircular mouth.
                        let axial = plane.normal.cross(&bore.cylinder.axis);
                        if geometry::sign(&axial.norm_squared(), dec)? != Sign::Zero {
                            return Err(failure("a preserved bore tool has an unsupported face"));
                        }
                        bore.planes.push(plane.clone());
                    }
                }
                Surface::Cylinder(_) => (),
                _ => return Err(failure("a preserved bore tool has an unsupported face")),
            }
        }
        Ok(())
    }

    /// Reads one subtraction's tool. A tool with one cone is one countersunk bore: its
    /// narrowest cylinder is the wall, and anything wider must be its entry.
    fn tool_walls(
        &self,
        tool: usize,
        stationed: bool,
        dec: &mut Decider,
    ) -> Result<Vec<ToolBore>, Refusal> {
        let (mut bores, mut heads) = (Vec::new(), Vec::new());
        self.tool_bores(tool, &mut bores, &mut heads, dec)?;
        if heads.is_empty() {
            return Ok(bores);
        }
        if heads.len() > 1 || bores.is_empty() {
            return Err(failure("a preserved bore tool has an unsupported face"));
        }
        let mut narrowest = 0;
        for (index, bore) in bores.iter().enumerate() {
            let excess = &bore.cylinder.radius_sq - &bores[narrowest].cylinder.radius_sq;
            if geometry::sign(&excess, dec)? == Sign::Negative {
                narrowest = index;
            }
        }
        let mut wall = bores.swap_remove(narrowest);
        wall.head = heads.pop().map(|head| ToolHead { stationed, ..head });
        wall.entries = bores;
        Ok(vec![wall])
    }

    /// Proves a countersunk tool removes only its two bands from the source. Wider
    /// bands lie wholly beyond the cap the cone opens on, so they bound no material.
    fn check_head(
        &self,
        bore: &ToolBore,
        head: &ToolHead,
        skin: &Skin,
        caps: [usize; 2],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let planes = skin.patches.iter().map(|patch| &patch.plane);
        let mouth = bored::certify_head(planes, &bore.cylinder, &head.cone, caps, dec)
            .map_err(source_refusal)?;
        let opening = &skin.patches[mouth].plane;
        for entry in &bore.entries {
            for plane in &entry.planes {
                let height = plane.origin.minus(&opening.origin).dot(&opening.normal);
                if geometry::sign(&height, dec)? == Sign::Negative {
                    return Err(failure("a countersink's entry reaches into its target"));
                }
            }
        }
        Ok(())
    }

    /// Proves the tool holds the axis over every interval between its end planes
    /// inside the caps. Bands share one wall, so the whole section is removed.
    /// Another bore's band cannot be the witness: the bores are proved apart. Under
    /// a head, the cone's own band holds every interval past their shared ring and a
    /// band of the wall every interval short of it, so each height keeps its radius.
    /// That ring is a stop, already present for a recipe hole's station.
    /// Returns one axis point per interval.
    fn check_through(
        &self,
        tool: usize,
        bore: &ToolBore,
        caps: [&Plane; 2],
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        let cylinder = &bore.cylinder;
        let [mut lower, mut upper] = caps.map(|cap| bored::axis_point(cylinder, cap));
        if geometry::sign(&upper.minus(&lower).dot(&cylinder.axis), dec)? != Sign::Positive {
            std::mem::swap(&mut lower, &mut upper);
        }
        let height = |p: &Point3| p.minus(&lower).dot(&cylinder.axis);
        let mut stops = vec![lower.clone(), upper];
        let head_planes = bore.head.iter().flat_map(|head| &head.planes);
        let mut bounds: Vec<_> = bore
            .planes
            .iter()
            .chain(head_planes)
            .map(|plane| bored::axis_point(cylinder, plane))
            .collect();
        let seam = bore
            .head
            .as_ref()
            .map(|head| bored::seam(cylinder, &head.cone));
        if bore.head.as_ref().is_some_and(|head| !head.stationed) {
            bounds.extend(seam.clone());
        }
        for stop in bounds {
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
        let mut middles = Vec::with_capacity(stops.len() - 1);
        for pair in stops.windows(2) {
            let half = pair[1]
                .minus(&pair[0])
                .scaled(&(Real::one() / Real::from_i64(2)));
            let middle = pair[0].plus(&half);
            let inside = |leaf: usize, dec: &mut Decider| {
                Ok::<_, Refusal>(
                    leaf_location(&self.leaves[leaf], &middle, dec)? == Location::Inside,
                )
            };
            let held = match (&bore.head, &seam) {
                (Some(head), Some(seam)) => {
                    let past = middle.minus(seam).dot(&head.cone.axis);
                    if geometry::sign(&past, dec)? == Sign::Positive {
                        inside(head.leaf, dec)?
                    } else {
                        let mut held = false;
                        for leaf in &bore.leaves {
                            held |= inside(*leaf, dec)?;
                        }
                        held
                    }
                }
                _ => self.register_location(tool, &middle, dec)? == Location::Inside,
            };
            if !held {
                return Err(failure("a preserved bore must pass through both caps"));
            }
            middles.push(middle);
        }
        Ok(middles)
    }

    /// Walks a recipe's subtractions of curved tools back from `body`, newest
    /// first, to the first step that is not one. The walk reads recipe faces and
    /// exact leaves, never a B-rep. `None` keeps a mismatched step's refusal.
    fn bore_chain<'a>(
        &self,
        recipe: &'a Recipe,
        body: &'a StepId,
        (registers, faces): Steps,
    ) -> Option<(&'a StepId, Vec<ToolChain>)> {
        let plane = |f: &usize| matches!(self.faces[*f].surface, Surface::Plane(_));
        let mut chain = Vec::new();
        let mut step = body;
        loop {
            let target = match recipe.step(step).map(|s| &s.op) {
                Some(Operation::Boolean {
                    kind: BooleanKind::Subtract,
                    target,
                    ..
                })
                | Some(Operation::Hole { target, .. })
                | Some(Operation::ExtrudeThrough { target, .. }) => target,
                _ => return Some((step, chain)),
            };
            if faces[step].iter().all(plane) {
                return Some((step, chain));
            }
            let Instruction::Boolean {
                kind: BooleanKind::Subtract,
                a,
                b,
            } = self.instructions[registers[step]]
            else {
                return None;
            };
            if a != registers[target] {
                return None;
            }
            let added: Vec<_> = faces[step]
                .iter()
                .copied()
                .filter(|f| faces[target].binary_search(f).is_err())
                .collect();
            let hole = matches!(
                recipe.step(step).map(|s| &s.op),
                Some(Operation::Hole { .. })
            );
            chain.push((b, added, hole));
            step = target;
        }
    }

    /// Splits a recipe's union into its body and the one cylinder leaf that
    /// could stand on it. A body under bores is curved too, but never one leaf.
    fn boss_union<'a>(
        &self,
        recipe: &'a Recipe,
        step: &'a StepId,
        (registers, faces): Steps,
    ) -> Option<(&'a StepId, &'a StepId)> {
        let Some(Operation::Boolean {
            kind: BooleanKind::Union,
            target,
            tool,
        }) = recipe.step(step).map(|s| &s.op)
        else {
            return None;
        };
        let planar = |step: &StepId| {
            faces[step]
                .iter()
                .all(|f| matches!(self.faces[*f].surface, Surface::Plane(_)))
        };
        let leaf =
            |step: &StepId| matches!(self.instructions[registers[step]], Instruction::Leaf(_));
        let (body, tool) = match (planar(target), planar(tool)) {
            (true, false) => (target, tool),
            (false, true) => (tool, target),
            (false, false) if leaf(tool) && !leaf(target) => (target, tool),
            (false, false) if leaf(target) && !leaf(tool) => (tool, target),
            _ => return None,
        };
        let (
            Instruction::Boolean {
                kind: BooleanKind::Union,
                ..
            },
            Instruction::Leaf(_),
        ) = (
            self.instructions[registers[step]],
            self.instructions[registers[tool]],
        )
        else {
            return None;
        };
        Some((body, tool))
    }

    /// One tool leaf's wall and two end planes. Anything else is no boss.
    fn boss_carriers(&self, register: usize) -> Option<(&Cylinder, [&Plane; 2])> {
        let Instruction::Leaf(leaf) = self.instructions[register] else {
            return None;
        };
        let mut wall = None;
        let mut ends = Vec::new();
        for constraint in self.leaves[leaf].constraints() {
            match constraint.surface() {
                Surface::Cylinder(cylinder) if wall.is_none() => wall = Some(cylinder),
                Surface::Plane(plane) => ends.push(plane),
                _ => return None,
            }
        }
        Some((wall?, ends.try_into().ok()?))
    }

    /// Reads each subtraction's bores and refuses a selected one.
    fn tool_chain_bores<'a>(
        &self,
        chain: &'a [ToolChain],
        refs: &[FaceRef],
        dec: &mut Decider,
    ) -> Result<Vec<ChainBore<'a>>, Refusal> {
        let mut bores = Vec::new();
        for (tool, added, hole) in chain {
            if refs
                .iter()
                .any(|r| added.iter().any(|f| r.face.matches(&self.faces[*f].origin)))
            {
                return Err(failure("a selected bore needs a curved arrangement"));
            }
            bores.extend(
                self.tool_walls(*tool, *hole, dec)?
                    .into_iter()
                    .map(|bore| (*tool, added.as_slice(), bore)),
            );
        }
        Ok(bores)
    }

    /// Proves complete cylinders of a tool chain unchanged by a draft of `skin`. A
    /// bore through `boss` runs past the shoulder to the boss's end.
    fn prove_bores(
        &self,
        chain: &[ToolChain],
        bores: Vec<ChainBore>,
        (skin, boss): (&Skin, Option<&PreservedBoss>),
        refs: &[FaceRef],
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<PreservedBores, Refusal> {
        for (.., bore) in &bores {
            bored::check_axial(&bore.cylinder, definition, dec).map_err(source_refusal)?;
        }
        let cylinders: Vec<_> = bores.iter().map(|(.., bore)| &bore.cylinder).collect();
        let caps = bored::certify_bores(
            skin.patches.iter().map(|patch| &patch.plane),
            &cylinders,
            dec,
        )
        .map_err(source_refusal)?;
        for (.., bore) in &bores {
            if let Some(head) = &bore.head {
                self.check_head(bore, head, skin, caps, dec)?;
            }
        }
        let selected = skin.resolve_faces(refs, dec)?;
        if caps.iter().any(|cap| selected.contains(cap)) {
            return Err(failure("a selected bore cap moves the cap-hole trim"));
        }
        let mut preserved = Vec::with_capacity(bores.len());
        for (tool, added, bore) in bores {
            let ends = caps.map(|cap| match boss {
                Some(boss) if boss.through && cap == boss.shoulder => &boss.boss.top,
                _ => &skin.patches[cap].plane,
            });
            let middles = self.check_through(tool, &bore, ends, dec)?;
            // A band beyond the caps, such as a hole's entry, bounds no source material.
            let mut walls = Vec::new();
            for wall in added {
                let Surface::Cylinder(carrier) = &self.faces[*wall].surface else {
                    continue;
                };
                if !same_wall(carrier, &bore.cylinder, dec)? {
                    continue;
                }
                let leaf = &self.leaves[self.faces[*wall].leaf];
                let mut within = false;
                for middle in &middles {
                    within |= leaf_location(leaf, middle, dec)? == Location::Inside;
                }
                if within {
                    walls.push(*wall);
                }
            }
            let head = bore.head.map(|head| PreservedHead {
                walls: added
                    .iter()
                    .copied()
                    .filter(|f| matches!(self.faces[*f].surface, Surface::Cone(_)))
                    .collect(),
                cone: head.cone,
            });
            preserved.push(PreservedBore {
                walls,
                cylinder: bore.cylinder,
                head,
            });
        }
        Ok(PreservedBores {
            tools: chain.iter().map(|(tool, ..)| *tool).collect(),
            bores: preserved,
        })
    }

    /// Proves one complete cylinder tool stands on a cap of `skin` and stays
    /// unchanged by its draft. The tool starts on a cap and rises beyond it. The
    /// B-rep sees only that cap's ring, so a tool sunk into the target stays outside
    /// this proof.
    fn prove_boss(
        &self,
        (register, tool_faces): (usize, &[usize]),
        skin: &Skin,
        refs: &[FaceRef],
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<PreservedBoss, Refusal> {
        let (cylinder, [first, second]) = self
            .boss_carriers(register)
            .expect("recognised boss tools keep their carriers");
        let value = |plane: &Plane, point: &Point3| point.minus(&plane.origin).dot(&plane.normal);
        let mut found = None;
        for (low, high) in [(first, second), (second, first)] {
            for (index, patch) in skin.patches.iter().enumerate() {
                let carrier = patch.plane.normal.cross(&low.normal).norm_squared();
                if geometry::sign(&carrier, dec)? == Sign::Zero
                    && geometry::sign(&value(&patch.plane, &low.origin), dec)? == Sign::Zero
                    && geometry::sign(&value(&patch.plane, &high.origin), dec)? == Sign::Positive
                {
                    found = Some((index, high));
                }
            }
        }
        let Some((shoulder, high)) = found else {
            return Err(failure(
                "a preserved boss must stand on a cap of its target",
            ));
        };
        let mut top = high.clone();
        let outward = top.normal.dot(&skin.patches[shoulder].plane.normal);
        if geometry::sign(&outward, dec)? == Sign::Negative {
            top.normal = top.normal.scaled(&Real::from_i64(-1));
        }
        let mut walls = Vec::new();
        let mut caps = Vec::new();
        for face in tool_faces {
            match &self.faces[*face].surface {
                Surface::Cylinder(_) => walls.push(*face),
                Surface::Plane(plane)
                    if geometry::sign(&value(&top, &plane.origin), dec)? == Sign::Zero =>
                {
                    caps.push(*face);
                }
                _ => (),
            }
        }
        let names = |faces: &[usize]| -> Vec<Origin> {
            faces
                .iter()
                .map(|f| self.faces[*f].origin.clone())
                .collect()
        };
        let (origins, top_origins) = (names(&walls), names(&caps));
        if refs.iter().any(|r| {
            origins
                .iter()
                .chain(&top_origins)
                .any(|origin| r.face.matches(origin))
        }) {
            return Err(failure("a selected boss needs a curved arrangement"));
        }
        bored::check_axial(cylinder, definition, dec).map_err(source_refusal)?;
        let boss = Boss {
            cylinder: cylinder.clone(),
            origins,
            top,
            top_origins,
        };
        let cap =
            bossed::certify_footprint(skin.patches.iter().map(|patch| &patch.plane), &boss, dec)
                .map_err(source_refusal)?;
        if cap != shoulder {
            return Err(failure(
                "a preserved boss must stand on a cap of its target",
            ));
        }
        if skin.resolve_faces(refs, dec)?.contains(&cap) {
            return Err(failure("a selected shoulder moves the boss's ring"));
        }
        Ok(PreservedBoss {
            tool: register,
            boss,
            wall: walls,
            top: caps,
            shoulder: cap,
            through: false,
        })
    }

    /// The planar target and the curved material to keep around its draft: a
    /// chain of subtractions of complete cylinders, one complete cylinder standing
    /// on a cap, or one of each in either history order. A bore newer than the boss
    /// may run through it. Any other curved source returns `None` and keeps its
    /// refusal.
    pub(crate) fn preserved_source(
        &self,
        recipe: &Recipe,
        body: &StepId,
        refs: &[FaceRef],
        definition: &DraftDefinition,
        steps: Steps,
        dec: &mut Decider,
    ) -> Result<Option<(StepId, Skin, Preserved)>, Refusal> {
        let (registers, faces) = steps;
        let planar = |step: &StepId| {
            faces[step]
                .iter()
                .all(|f| matches!(self.faces[*f].surface, Surface::Plane(_)))
        };
        let Some((stop, mut chain)) = self.bore_chain(recipe, body, steps) else {
            return Ok(None);
        };
        let mut boss_tool = None;
        // Only a bore made after the boss can have cut it.
        let mut cuts_boss = false;
        let target = if planar(stop) {
            stop
        } else {
            let Some((under, tool)) = self.boss_union(recipe, stop, steps) else {
                return Ok(None);
            };
            let Some((target, older)) = self.bore_chain(recipe, under, steps) else {
                return Ok(None);
            };
            if !planar(target) || self.boss_carriers(registers[tool]).is_none() {
                return Ok(None);
            }
            cuts_boss = !chain.is_empty();
            chain.extend(older);
            boss_tool = Some(tool);
            target
        };
        if chain.is_empty() && boss_tool.is_none() {
            return Ok(None);
        }
        chain.reverse();
        let bores = self.tool_chain_bores(&chain, refs, dec)?;
        // A boss beside several bores is a wider family than this proof claims.
        if boss_tool.is_some() && bores.len() > 1 {
            return Ok(None);
        }
        // A countersunk bore stands alone.
        if (boss_tool.is_some() || bores.len() > 1)
            && bores.iter().any(|(.., bore)| bore.head.is_some())
        {
            return Err(failure("a preserved bore tool has an unsupported face"));
        }
        let skin = self.draft_source_skin(recipe, target, registers, faces, dec)?;
        if !skin.is_convex(dec)? {
            return Err(failure(if boss_tool.is_some() {
                "a preserved boss needs a convex planar envelope"
            } else {
                "a preserved bore needs a convex planar envelope"
            }));
        }
        let boss = match boss_tool {
            Some(tool) => {
                let source = (registers[tool], faces[tool].as_slice());
                let mut boss = self.prove_boss(source, &skin, refs, definition, dec)?;
                if let (true, [(.., bore)]) = (cuts_boss, bores.as_slice()) {
                    boss.through = bored::axis_inside(&boss.boss.cylinder, &bore.cylinder, dec)?;
                }
                Some(Box::new(boss))
            }
            None => None,
        };
        let bores = if chain.is_empty() {
            None
        } else {
            let envelope = (&skin, boss.as_deref());
            Some(self.prove_bores(&chain, bores, envelope, refs, definition, dec)?)
        };
        Ok(Some((target.clone(), skin, Preserved { bores, boss })))
    }

    /// Applies the recipe's own tools to a drafted target, republishing the
    /// unchanged walls. Nothing returns once the target itself is gone. Proved
    /// apart, a boss and a bore commute, so either history restores this way. A bore
    /// through the boss is newer than it, so its whole tool cuts the restored boss too.
    pub(crate) fn restore_preserved(
        &mut self,
        drafted: usize,
        preserved: &Preserved,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if drafted == 0 {
            return Ok(0);
        }
        let mut result = drafted;
        if let Some(preserved) = &preserved.boss {
            self.preserve_faces(&preserved.wall);
            self.preserve_faces(&preserved.top);
            result = self.boolean(BooleanKind::Union, result, preserved.tool, dec)?;
        }
        if let Some(preserved) = &preserved.bores {
            for bore in &preserved.bores {
                self.preserve_faces(&bore.walls);
                if let Some(head) = &bore.head {
                    self.preserve_faces(&head.walls);
                }
            }
            for tool in &preserved.tools {
                result = self.boolean(BooleanKind::Subtract, result, *tool, dec)?;
            }
        }
        Ok(result)
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
