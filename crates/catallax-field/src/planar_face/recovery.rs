//! Recover a finite planar Boolean boundary from the active field definitions.
//! No B-rep construction, field display, rounded samples or winning-distance
//! instruction supplies a trim or an origin.

use crate::{Instruction, Location, Plan, geometry};
use catallax_features::draft::{
    arrangement::{Patch, Skin},
    cell::{PlanarCell, PlanePatch},
};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, Origin, Role};
use catallax_shapes::{Point3, Surface, Vector3};
use std::collections::{BTreeMap, BTreeSet};

pub(super) fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.planar-domain", detail)
}

fn operation(kind: BooleanKind, a: bool, b: bool) -> bool {
    match kind {
        BooleanKind::Union => a || b,
        BooleanKind::Intersect => a && b,
        BooleanKind::Subtract => a && !b,
    }
}

fn location(leaf: &geometry::Leaf, p: &Point3, dec: &mut Decider) -> Result<Location, Refusal> {
    if let geometry::Leaf::Feature(feature) = leaf {
        return feature.classify(p, dec);
    }
    let xyz = [p.x.clone(), p.y.clone(), p.z.clone()];
    let mut boundary = false;
    for constraint in leaf.constraints() {
        match geometry::sign(&constraint.polynomial.eval_real(&xyz), dec)? {
            Sign::Positive => return Ok(Location::Outside),
            Sign::Negative => (),
            Sign::Zero => boundary = true,
        }
    }
    Ok(if boundary {
        Location::Boundary
    } else {
        Location::Inside
    })
}

fn append(target: &mut Skin, part: &Skin) {
    let first = target.vertices.len();
    target.vertices.extend(part.vertices.iter().cloned());
    target.patches.extend(part.patches.iter().map(|p| {
        Patch {
            plane: p.plane.clone(),
            loops: p
                .loops
                .iter()
                .map(|c| c.iter().map(|i| first + i).collect())
                .collect(),
            origins: p.origins.clone(),
        }
    }));
}

#[derive(Default)]
struct State {
    minus: bool,
    plus: bool,
    origins: BTreeSet<Origin>,
}

pub(super) fn skin(plan: &Plan, root: usize, dec: &mut Decider) -> Result<Skin, Refusal> {
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
    let bound = plan.node_radius(root, dec)?;
    let mut sources = BTreeMap::new();
    let mut joined = Skin {
        vertices: Vec::new(),
        patches: Vec::new(),
    };
    for i in (0..=root).filter(|i| active[*i]) {
        let Instruction::Leaf(leaf) = plan.instructions[i] else {
            continue;
        };
        if sources.contains_key(&leaf) {
            continue;
        }
        let members: Vec<_> = plan.faces.iter().filter(|f| f.leaf == leaf).collect();
        // The unchanged fallback receives a tape node, not the recipe's named
        // snapshots. Geometry-folded creators cannot establish scoped ownership
        // here. Refuse that case rather than revive an unrelated alias.
        let creators: BTreeSet<_> = members.iter().map(|f| &f.origin.step).collect();
        if creators.len() != 1
            || members.iter().any(|f| {
                !matches!(f.origin.role, Role::Start | Role::End | Role::Side { .. })
                    || !matches!(f.surface, Surface::Plane(_))
            })
        {
            return Err(bad(
                "planar Boolean recovery needs unaliased elementary planar source definitions",
            ));
        }
        let origins: Vec<_> = members.iter().map(|f| &f.origin).collect();
        let part = match &plan.leaves[leaf] {
            geometry::Leaf::Feature(feature) => feature.planar_skin(&origins, Some(&bound), dec)?,
            geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_) => {
                let witness = plan.leaves[leaf].draft_witness(dec)?;
                let mut supports = Vec::new();
                for f in members {
                    let Surface::Plane(mut plane) = f.surface.clone() else {
                        unreachable!()
                    };
                    match geometry::sign(&witness.minus(&plane.origin).dot(&plane.normal), dec)? {
                        Sign::Positive => plane.normal = plane.normal.scaled(&Real::from_i64(-1)),
                        Sign::Negative => (),
                        Sign::Zero => {
                            return Err(bad("source plane passes through its strict witness"));
                        }
                    }
                    supports.push(PlanePatch {
                        plane,
                        origins: vec![f.origin.clone()],
                    });
                }
                let cell = PlanarCell::new(supports, dec)?
                    .ok_or_else(|| bad("source planes have no full-dimensional bounded cell"))?;
                Skin {
                    vertices: cell.vertices.clone(),
                    patches: cell
                        .supports
                        .iter()
                        .zip(&cell.cycles)
                        .filter(|(_, c)| !c.is_empty())
                        .map(|(s, c)| Patch {
                            plane: s.plane.clone(),
                            loops: vec![c.clone()],
                            origins: s.origins.clone(),
                        })
                        .collect(),
                }
            }
            _ => {
                return Err(bad(
                    "derived or curved source has no polygonal Boolean recovery proof",
                ));
            }
        };
        append(&mut joined, &part);
        sources.insert(leaf, part);
    }

    let components = joined.partition(
        |p, dec| {
            let mut inside = vec![false; root + 1];
            for i in (0..=root).filter(|i| active[*i]) {
                inside[i] = match plan.instructions[i] {
                    Instruction::Empty => false,
                    Instruction::Leaf(leaf) => match location(&plan.leaves[leaf], p, dec)? {
                        Location::Inside => true,
                        Location::Outside => false,
                        _ => {
                            return Err(bad(
                                "arrangement cell has no strict source membership proof",
                            ));
                        }
                    },
                    Instruction::Boolean { kind, a, b } => operation(kind, inside[a], inside[b]),
                };
            }
            Ok(inside[root])
        },
        |p: &Point3, normal: &Vector3, dec: &mut Decider| {
            let mut states: Vec<State> = (0..=root).map(|_| State::default()).collect();
            for i in (0..=root).filter(|i| active[*i]) {
                match plan.instructions[i] {
                    Instruction::Empty => (),
                    Instruction::Leaf(leaf) => {
                        let state = &mut states[i];
                        match location(&plan.leaves[leaf], p, dec)? {
                            loc @ (Location::Inside | Location::Outside) => {
                                state.minus = loc == Location::Inside;
                                state.plus = state.minus;
                            }
                            Location::Boundary => {
                                let source = &sources[&leaf];
                                let mut sense = None;
                                for patch in &source.patches {
                                    if !source.contains_patch(patch, p, dec)? {
                                        continue;
                                    }
                                    let aligned =
                                        match geometry::sign(&patch.plane.normal.dot(normal), dec)?
                                        {
                                            Sign::Positive => true,
                                            Sign::Negative => false,
                                            Sign::Zero => {
                                                return Err(bad(
                                                    "facet has no transverse source germ",
                                                ));
                                            }
                                        };
                                    if sense.is_some_and(|s| s != aligned) {
                                        return Err(bad(
                                            "source facet has contradictory material senses",
                                        ));
                                    }
                                    sense = Some(aligned);
                                    state.origins.extend(patch.origins.iter().cloned());
                                }
                                state.minus = sense.ok_or_else(|| {
                                    bad("source boundary has no finite defining patch")
                                })?;
                                state.plus = !state.minus;
                            }
                            Location::Indeterminate => {
                                return Err(bad("source facet classification is unresolved"));
                            }
                        }
                    }
                    Instruction::Boolean { kind, a, b } => {
                        let (previous, current) = states.split_at_mut(i);
                        let state = &mut current[0];
                        state.minus = operation(kind, previous[a].minus, previous[b].minus);
                        state.plus = operation(kind, previous[a].plus, previous[b].plus);
                        if state.minus != state.plus {
                            state.origins.extend(
                                previous[a]
                                    .origins
                                    .iter()
                                    .chain(&previous[b].origins)
                                    .cloned(),
                            );
                        }
                    }
                }
            }
            let state = &states[root];
            if !state.minus || state.plus || state.origins.is_empty() {
                return Err(bad(
                    "Boolean outward germ has no surviving named physical facet",
                ));
            }
            Ok(state.origins.iter().cloned().collect())
        },
        dec,
    )?;
    let mut result = Skin {
        vertices: Vec::new(),
        patches: Vec::new(),
    };
    for component in components {
        append(&mut result, &component.skin());
    }
    Ok(result)
}
