//! Recover a finite planar Boolean boundary from the active field definitions.
//! No B-rep construction, field display, rounded samples or winning-distance
//! instruction supplies a trim or an origin.

use crate::{FaceState, Instruction, Location, Plan, geometry};
use catallax_features::draft::{
    arrangement::{Patch, Skin},
    cell::{PlanarCell, PlanePatch},
};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, Origin, Role};
use catallax_shapes::{Plane, Point3, Surface, Vector3};
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
    origins: BTreeSet<(Origin, usize)>,
}

fn same_plane(a: &Plane, b: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        geometry::sign(&a.normal.cross(&b.normal).norm_squared(), dec)? == Sign::Zero
            && geometry::sign(&a.origin.minus(&b.origin).dot(&b.normal), dec)? == Sign::Zero,
    )
}

pub(crate) fn skin(
    plan: &Plan,
    root: usize,
    named: Option<(&[usize], &[usize])>,
    dec: &mut Decider,
) -> Result<Skin, Refusal> {
    let mut joined = Skin {
        vertices: Vec::new(),
        patches: Vec::new(),
    };
    if root == 0 {
        return Ok(joined);
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
    let bound = plan.node_radius(root, dec)?;
    let mut sources = BTreeMap::new();
    for i in (0..=root).filter(|i| active[*i]) {
        let Instruction::Leaf(leaf) = plan.instructions[i] else {
            continue;
        };
        if sources.contains_key(&leaf) {
            continue;
        }
        let members: Vec<_> = match named {
            Some((faces, domains)) => faces
                .iter()
                .zip(domains)
                .map(|(&face, &domain)| (&plan.faces[face], domain))
                .filter(|(f, _)| f.leaf == leaf)
                .collect(),
            None => plan
                .faces
                .iter()
                .filter(|f| f.leaf == leaf)
                .map(|f| (f, 0))
                .collect(),
        };
        // A tape-only caller cannot distinguish creators folded onto the
        // same leaf. A named snapshot has the ownership history to do so.
        if named.is_none()
            && (members
                .iter()
                .map(|(f, _)| &f.origin.step)
                .collect::<BTreeSet<_>>()
                .len()
                != 1
                || members.iter().any(|(f, _)| {
                    !matches!(f.origin.role, Role::Start | Role::End | Role::Side { .. })
                        || !matches!(f.surface, Surface::Plane(_))
                }))
        {
            return Err(bad(
                "planar Boolean recovery needs unaliased elementary planar source definitions",
            ));
        }
        let origins: Vec<_> = members.iter().map(|(f, _)| &f.origin).collect();
        let part = match &plan.leaves[leaf] {
            geometry::Leaf::Feature(feature) => feature.planar_skin(&origins, Some(&bound), dec)?,
            geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_) => {
                let mut supports = Vec::new();
                if named.is_some() {
                    // Include unnamed cuts so that a graded band's finite
                    // trim and its material side come from its actual cell.
                    for constraint in plan.leaves[leaf].constraints() {
                        let Surface::Plane(plane) = constraint.surface() else {
                            return Err(bad("a curved source has no planar Boolean recovery"));
                        };
                        let mut origins = Vec::new();
                        for (face, _) in &members {
                            if let Surface::Plane(candidate) = &face.surface
                                && same_plane(plane, candidate, dec)?
                            {
                                origins.push(face.origin.clone());
                            }
                        }
                        supports.push(PlanePatch {
                            plane: plane.clone(),
                            origins,
                        });
                    }
                } else {
                    let witness = plan.leaves[leaf].draft_witness(dec)?;
                    for (f, _) in &members {
                        let Surface::Plane(mut plane) = f.surface.clone() else {
                            unreachable!()
                        };
                        match geometry::sign(&witness.minus(&plane.origin).dot(&plane.normal), dec)?
                        {
                            Sign::Positive => {
                                plane.normal = plane.normal.scaled(&Real::from_i64(-1))
                            }
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
        sources.insert(leaf, (part, members));
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
                                let (source, members) = &sources[&leaf];
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
                                    if named.is_none() {
                                        state
                                            .origins
                                            .extend(patch.origins.iter().cloned().map(|o| (o, 0)));
                                        continue;
                                    }
                                    for (face, domain) in members {
                                        if patch.origins.contains(&face.origin)
                                            && let Surface::Plane(candidate) = &face.surface
                                            && same_plane(&patch.plane, candidate, dec)?
                                        {
                                            state.origins.insert((face.origin.clone(), *domain));
                                        }
                                    }
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
            if !state.minus || state.plus {
                return Err(bad("classified boundary and its outward germ disagree"));
            }
            if named.is_none() && !state.origins.is_empty() {
                return Ok(state
                    .origins
                    .iter()
                    .map(|(origin, _)| origin.clone())
                    .collect());
            }
            let location = |inside| {
                if inside {
                    Location::Inside
                } else {
                    Location::Outside
                }
            };
            let mut origins = BTreeSet::new();
            let mut domains = Vec::new();
            let plus: Vec<_> = states.iter().map(|state| location(state.plus)).collect();
            let minus: Vec<_> = states.iter().map(|state| location(state.minus)).collect();
            for (origin, domain) in &state.origins {
                if *domain == 0
                    || plan.domain_state(*domain, &plus, &minus, &mut domains) == FaceState::OnFace
                {
                    origins.insert(origin.clone());
                }
            }
            if origins.is_empty() {
                return Err(bad(
                    "Boolean outward germ has no surviving named physical facet",
                ));
            }
            Ok(origins.into_iter().collect())
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
