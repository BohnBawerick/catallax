//! Named-patch survival is separate from shared scalar geometry.
//!
//! A later coincident carrier cannot restore a name removed by an earlier
//! Boolean. Each result snapshot retains the boundary predicates on every
//! path by which its source faces reached that result.

use crate::{Direction, Face, FaceQuery, FaceState, Instruction, Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{Operation, StepId};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug)]
pub(crate) enum Domain {
    Source,
    Absent,
    Boundary { node: usize, parent: usize },
    Either { a: usize, b: usize },
}

pub(crate) type Snapshots = BTreeMap<StepId, Vec<usize>>;

impl Plan {
    pub(crate) fn compact_domains(&mut self, live: &[bool], registers: &[usize]) {
        let mut active = vec![false; self.domains.len()];
        active[0] = true;
        active[1] = true;
        for face in self.faces.iter().filter(|face| face.in_result) {
            active[face.domain] = true;
        }
        for id in (2..self.domains.len()).rev() {
            if !active[id] {
                continue;
            }
            match self.domains[id] {
                Domain::Boundary { node, parent } => {
                    // Folding can discard a whole branch only through an empty
                    // result, which snapshot_domains has already made Absent.
                    debug_assert!(live[node]);
                    active[parent] = true;
                }
                Domain::Either { a, b } => {
                    active[a] = true;
                    active[b] = true;
                }
                Domain::Source | Domain::Absent => (),
            }
        }
        let mut remap = vec![0; self.domains.len()];
        let mut domains = Vec::new();
        for (id, domain) in self.domains.iter().enumerate() {
            if !active[id] {
                continue;
            }
            remap[id] = domains.len();
            domains.push(match *domain {
                Domain::Boundary { node, parent } => Domain::Boundary {
                    node: registers[node],
                    parent: remap[parent],
                },
                Domain::Either { a, b } => Domain::Either {
                    a: remap[a],
                    b: remap[b],
                },
                constant => constant,
            });
        }
        for face in &mut self.faces {
            face.domain = remap[face.domain];
        }
        self.domains = domains;
    }

    fn either_domain(&mut self, a: usize, b: usize) -> usize {
        if a == b || b == 1 {
            return a;
        }
        if a == 1 {
            return b;
        }
        if a == 0 || b == 0 {
            return 0;
        }
        let id = self.domains.len();
        self.domains.push(Domain::Either { a, b });
        id
    }

    pub(crate) fn snapshot_domains(
        &mut self,
        op: &Operation,
        node: usize,
        faces: &[usize],
        step_faces: &BTreeMap<StepId, Vec<usize>>,
        snapshots: &Snapshots,
    ) -> Vec<usize> {
        let (a, b) = match op {
            Operation::Boolean { target, tool, .. } => (Some(target), Some(tool)),
            Operation::ExtrudeThrough { target, .. } | Operation::Hole { target, .. } => {
                (Some(target), None)
            }
            Operation::Lattice { region, .. } => (Some(region), None),
            Operation::Pattern { body, .. }
            | Operation::Fillet { body, .. }
            | Operation::Chamfer { body, .. }
            | Operation::Shell { body, .. }
            | Operation::Mirror {
                body,
                keep_original: true,
                ..
            } => (Some(body), None),
            _ => (None, None),
        };
        let inherited = |step: Option<&StepId>, face: usize| {
            step.and_then(|step| {
                step_faces[step]
                    .binary_search(&face)
                    .ok()
                    .map(|index| snapshots[step][index])
            })
        };
        let mut gates = BTreeMap::new();
        faces.iter().map(|&face| {
            let parent = match (inherited(a, face), inherited(b, face)) {
                (Some(a), Some(b)) => self.either_domain(a, b),
                (Some(a), None) | (None, Some(a)) => a,
                (None, None) => self.faces[face].domain,
            };
            if node == 0 || parent == 1 {
                return 1;
            }
            // The candidate's own membership and trim checks already prove a
            // leaf boundary. Keep inherited predicates even after tape folding.
            if matches!(self.instructions[node], Instruction::Leaf(leaf) if leaf == self.faces[face].leaf) {
                return parent;
            }
            *gates.entry(parent).or_insert_with(|| {
                let id = self.domains.len();
                self.domains.push(Domain::Boundary { node, parent });
                id
            })
        }).collect()
    }

    pub(crate) fn copy_domain(
        &mut self,
        domain: usize,
        remap: &[usize],
        copied: &mut BTreeMap<usize, usize>,
    ) -> usize {
        if domain < 2 {
            return domain;
        }
        if let Some(&id) = copied.get(&domain) {
            return id;
        }
        let predicate = match self.domains[domain] {
            Domain::Source | Domain::Absent => unreachable!("constant domain handled above"),
            Domain::Boundary { node, parent } => Domain::Boundary {
                node: remap[node],
                parent: self.copy_domain(parent, remap, copied),
            },
            Domain::Either { a, b } => Domain::Either {
                a: self.copy_domain(a, remap, copied),
                b: self.copy_domain(b, remap, copied),
            },
        };
        let id = self.domains.len();
        self.domains.push(predicate);
        copied.insert(domain, id);
        id
    }

    /// Strict point membership can exclude a historical boundary, never admit it.
    pub(crate) fn domain_excludes_point(
        &self,
        domain: usize,
        locations: &[Location],
        states: &mut Vec<FaceState>,
    ) -> bool {
        use FaceState::{Indeterminate, NotOnFace};
        states.clear();
        for predicate in &self.domains[..=domain] {
            states.push(match *predicate {
                Domain::Source => Indeterminate,
                Domain::Absent => NotOnFace,
                Domain::Boundary { node, parent } => {
                    if matches!(locations[node], Location::Inside | Location::Outside) {
                        NotOnFace
                    } else {
                        states[parent]
                    }
                }
                Domain::Either { a, b } => {
                    if states[a] == NotOnFace && states[b] == NotOnFace {
                        NotOnFace
                    } else {
                        Indeterminate
                    }
                }
            });
        }
        states[domain] == NotOnFace
    }

    pub(crate) fn domain_normal(
        &self,
        face: &Face,
        point: &[Rational; 3],
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Option<FaceQuery>, Refusal> {
        if face.domain == 0 {
            return self
                .oriented_normal(face, point, &mut storage.rays, dec)?
                .map(|n| {
                    Ok(FaceQuery {
                        state: FaceState::OnFace,
                        direction: crate::unit_direction(n, dec)?,
                    })
                })
                .transpose();
        }
        let normal = crate::normal(face, point);
        if crate::geometry::sign(&normal.norm_squared(), dec)? == Sign::Zero {
            return Ok(None);
        }
        let plus = self.ray(point, &normal, &mut storage.rays, dec)?;
        let minus = self.ray(
            point,
            &normal.scaled(&catallax_numbers::Real::from_i64(-1)),
            &mut storage.opposite_rays,
            dec,
        )?;
        let outward = match (plus, minus) {
            (Location::Outside, Location::Inside) => normal,
            (Location::Inside, Location::Outside) => {
                normal.scaled(&catallax_numbers::Real::from_i64(-1))
            }
            _ => return Ok(None),
        };
        let state = self.domain_state(
            face.domain,
            &storage.rays,
            &storage.opposite_rays,
            &mut storage.domains,
        );
        Ok(Some(FaceQuery {
            state,
            direction: match state {
                FaceState::OnFace => crate::unit_direction(outward, dec)?,
                FaceState::NotOnFace => Direction::Undefined,
                FaceState::Indeterminate => Direction::Ambiguous,
            },
        }))
    }

    pub(crate) fn domain_state(
        &self,
        domain: usize,
        plus: &[Location],
        minus: &[Location],
        states: &mut Vec<FaceState>,
    ) -> FaceState {
        use FaceState::*;
        states.clear();
        for predicate in &self.domains[..=domain] {
            states.push(match *predicate {
                Domain::Source => OnFace,
                Domain::Absent => NotOnFace,
                Domain::Boundary { node, parent } => {
                    let boundary = match (plus[node], minus[node]) {
                        (Location::Inside, Location::Outside)
                        | (Location::Outside, Location::Inside) => OnFace,
                        (Location::Inside, Location::Inside)
                        | (Location::Outside, Location::Outside) => NotOnFace,
                        _ => Indeterminate,
                    };
                    match (states[parent], boundary) {
                        (NotOnFace, _) | (_, NotOnFace) => NotOnFace,
                        (OnFace, OnFace) => OnFace,
                        _ => Indeterminate,
                    }
                }
                Domain::Either { a, b } => match (states[a], states[b]) {
                    (OnFace, _) | (_, OnFace) => OnFace,
                    (NotOnFace, NotOnFace) => NotOnFace,
                    _ => Indeterminate,
                },
            });
        }
        states[domain]
    }
}
