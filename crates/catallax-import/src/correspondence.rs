//! Conservative correspondence of complete accepted entity definitions.
//!
//! Entity numbers are addresses within one hash-fixed source, never matching
//! evidence. Equality includes the oriented surface, every finite boundary,
//! carrier, endpoint and sense. Different parameterizations remain unresolved.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{
    Decider, Rational, Refusal, operation_budget_active, with_operation_budget,
};
use catallax_topology::{BodyId, Model};
use serde::{Deserialize, Serialize};

use crate::{Face, Imported, OrientedEdge, Placement, Solid, StepFile, Surface};

/// Kind of physical source entity; oriented edge uses are deduplicated.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EntityKind {
    /// One finite, oriented face.
    Face,
    /// One finite edge, independent of its direction of use in a face loop.
    Edge,
}

/// An entity address within one immutable source. Solid membership is explicit.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct SourceEntity {
    /// Physical entity kind.
    pub kind: EntityKind,
    /// Source solid key.
    pub solid: String,
    /// Source face or EDGE_CURVE key.
    pub key: String,
}

/// One verified one-to-one equality of complete accepted definitions.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EntityCorrespondence {
    /// Address in the old source.
    pub old: SourceEntity,
    /// Address in the replacement source.
    pub new: SourceEntity,
}

/// Why an old entity has no verified replacement.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UnresolvedReason {
    /// No equal verified definition exists in the new import. This can also
    /// mean a different parameterization or an unverified replacement geometry.
    ChangedOrRemoved,
    /// Equal definitions occur more than once in either source.
    AmbiguousGeometry,
    /// An unavailable definition of the same kind could be another equal
    /// definition, so uniqueness across both imports is not established.
    UnverifiedMultiplicity,
    /// The solid was left open, or its builder repaired geometry internally
    /// without returning the resulting finite definitions to this comparator.
    UnverifiedGeometry,
}

/// An unresolved old entity. Candidates are information, never a chosen target.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnresolvedEntity {
    /// Address in the old source.
    pub old: SourceEntity,
    /// Explicit unresolved cause.
    pub reason: UnresolvedReason,
    /// Equal verified new definitions, including all ambiguity candidates.
    pub candidates: Vec<SourceEntity>,
}

/// Hash-bound correspondence report. This is not a document edit or authority
/// to skip import/rebuild checks; serialized reports can be edited by callers.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReplacementCorrespondence {
    /// Format of this report.
    pub format: String,
    /// Actual immutable old source digest retained by the import producer.
    pub old_source_sha256: String,
    /// Actual immutable replacement source digest.
    pub new_source_sha256: String,
    /// Unique complete-definition matches, in old-address order.
    pub verified: Vec<EntityCorrespondence>,
    /// Every old entity without a verified singleton match.
    pub unresolved: Vec<UnresolvedEntity>,
    /// New entities without a verified old partner. These are not necessarily
    /// additions: changed, ambiguous and unverified geometry appears here too.
    pub unmatched_new: Vec<SourceEntity>,
}

#[derive(Clone, Debug)]
pub(crate) struct SourceSnapshot {
    sha256: String,
    file: StepFile,
    unverified: BTreeSet<String>,
}

impl SourceSnapshot {
    pub(crate) fn complete_definition(&self) -> Option<&StepFile> {
        self.unverified.is_empty().then_some(&self.file)
    }

    pub(crate) fn new(
        sha256: String,
        file: StepFile,
        bodies: &[(String, BodyId)],
        left_open: &[String],
        model: &Model,
    ) -> Self {
        let unverified = left_open
            .iter()
            .cloned()
            // Spline incidence repair currently returns marks but keeps
            // its accepted derivative private inside the strict builder.
            // Raw source definitions must not stand in for that result.
            .chain(
                bodies
                    .iter()
                    .filter(|(_, body)| !model.body(*body).marks.is_empty())
                    .map(|(solid, _)| solid.clone()),
            )
            .collect();
        Self {
            sha256,
            file,
            unverified,
        }
    }
}

fn clear_edge_names(edge: &mut OrientedEdge) {
    edge.edge.clear();
    edge.start.name.clear();
    edge.end.name.clear();
}

fn signature(mut face: Face, dec: &mut Decider) -> Result<Vec<u8>, Refusal> {
    face.name.clear();
    for boundary in &mut face.loops {
        dec.charge_work(1)?;
        if let Some(vertex) = &mut boundary.vertex {
            vertex.name.clear();
        }
        for edge in &mut boundary.edges {
            clear_edge_names(edge);
        }
    }
    let bytes = crate::canonical_bytes(&StepFile {
        mm_per_unit: Rational::one(), // recognized coordinates are already mm
        solids: vec![Solid {
            name: String::new(),
            faces: vec![face],
        }],
    });
    // Charge the complete representation, including large rational operands
    // and spline nets. Equality uses the actual bytes, not a hash collision.
    dec.charge_work(u64::try_from(bytes.len()).unwrap_or(u64::MAX))?;
    Ok(bytes)
}

fn edge_signature(edge: &OrientedEdge, dec: &mut Decider) -> Result<Vec<u8>, Refusal> {
    let mut edge = edge.clone();
    edge.forward = true; // a physical edge's coedge direction is not its identity
    signature(
        Face {
            name: String::new(),
            surface: Surface::Plane(Placement {
                origin: std::array::from_fn(|_| Rational::zero()),
                axis: [Rational::zero(), Rational::zero(), Rational::one()],
                ref_dir: [Rational::one(), Rational::zero(), Rational::zero()],
            }),
            same_sense: true,
            loops: vec![crate::FaceLoop {
                outer: true,
                orientation: true,
                edges: vec![edge],
                vertex: None,
            }],
        },
        dec,
    )
}

type Inventory = BTreeMap<SourceEntity, Option<Vec<u8>>>;

fn inventory(source: &SourceSnapshot, dec: &mut Decider) -> Result<Inventory, Refusal> {
    let mut result = BTreeMap::new();
    for solid in &source.file.solids {
        let verified = !source.unverified.contains(&solid.name);
        for face in &solid.faces {
            dec.charge_work(1)?;
            let address = SourceEntity {
                kind: EntityKind::Face,
                solid: solid.name.clone(),
                key: face.name.clone(),
            };
            let geometry = verified.then(|| signature(face.clone(), dec)).transpose()?;
            if result.insert(address, geometry).is_some() {
                return Err(Refusal::new(
                    "import.correspondence-duplicate-face",
                    "the accepted source repeats a face address within one solid",
                ));
            }
            for edge in face.loops.iter().flat_map(|boundary| &boundary.edges) {
                dec.charge_work(1)?;
                let address = SourceEntity {
                    kind: EntityKind::Edge,
                    solid: solid.name.clone(),
                    key: edge.edge.clone(),
                };
                let geometry = verified.then(|| edge_signature(edge, dec)).transpose()?;
                if let Some(previous) = result.insert(address, geometry.clone())
                    && previous != geometry
                {
                    return Err(Refusal::new(
                        "import.correspondence-inconsistent-edge",
                        "one accepted edge address has different complete definitions",
                    ));
                }
            }
        }
    }
    Ok(result)
}

/// Compare two completed imports without mutating either model or source.
///
/// Matches ignore entity numbers and require unique exact equality of complete
/// accepted finite definitions in both files. They do not use nearest geometry,
/// area/centroid scores, source ordering or partial control-point samples.
/// Different loop starts, equivalent parameterizations and changed geometry
/// may remain unresolved. Open solids and internally repaired spline geometry
/// whose resulting definitions are unavailable always remain unresolved.
///
/// Planar/seam/cylinder repairs performed by the recorded importer compare the
/// actual retained derivative, not the inconsistent original source. Public
/// record edits cannot change the producer-owned definitions or source digests.
/// The caller's operation allowance applies, with the unchanged 20M default
/// when no request is active. Failure returns no partial report.
pub fn compare(
    old: &Imported,
    new: &Imported,
    dec: &mut Decider,
) -> Result<ReplacementCorrespondence, Refusal> {
    if operation_budget_active() {
        return compare_inner(&old.comparison, &new.comparison, dec);
    }
    let (result, work) = with_operation_budget(20_000_000, || {
        compare_inner(&old.comparison, &new.comparison, dec)
    });
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "source correspondence exhausted its operation allowance",
        ));
    }
    result
}

fn compare_inner(
    old: &SourceSnapshot,
    new: &SourceSnapshot,
    dec: &mut Decider,
) -> Result<ReplacementCorrespondence, Refusal> {
    let old_entities = inventory(old, dec)?;
    let new_entities = inventory(new, dec)?;
    let mut unavailable = BTreeSet::new();
    for (entity, geometry) in old_entities.iter().chain(&new_entities) {
        dec.charge_work(1)?;
        if geometry.is_none() {
            unavailable.insert(entity.kind);
        }
    }
    let group = |entities: &Inventory| {
        let mut groups = BTreeMap::<_, Vec<_>>::new();
        for (entity, geometry) in entities {
            if let Some(geometry) = geometry {
                groups
                    .entry((entity.kind, geometry.clone()))
                    .or_default()
                    .push(entity.clone());
            }
        }
        groups
    };
    let old_groups = group(&old_entities);
    let new_groups = group(&new_entities);
    let mut verified = Vec::new();
    let mut unresolved = Vec::new();
    let mut matched = BTreeSet::new();
    for (old, geometry) in old_entities {
        dec.charge_work(1)?;
        let Some(geometry) = geometry else {
            unresolved.push(UnresolvedEntity {
                old,
                reason: UnresolvedReason::UnverifiedGeometry,
                candidates: Vec::new(),
            });
            continue;
        };
        let key = (old.kind, geometry);
        let candidates = if let Some(group) = new_groups.get(&key) {
            dec.charge_work(u64::try_from(group.len()).unwrap_or(u64::MAX))?;
            group.clone()
        } else {
            Vec::new()
        };
        if candidates.len() == 1 && old_groups[&key].len() == 1 && !unavailable.contains(&old.kind)
        {
            let new = candidates.into_iter().next().expect("one candidate");
            matched.insert(new.clone());
            verified.push(EntityCorrespondence { old, new });
        } else {
            unresolved.push(UnresolvedEntity {
                reason: if unavailable.contains(&old.kind) {
                    UnresolvedReason::UnverifiedMultiplicity
                } else if candidates.is_empty() {
                    UnresolvedReason::ChangedOrRemoved
                } else {
                    UnresolvedReason::AmbiguousGeometry
                },
                old,
                candidates,
            });
        }
    }
    Ok(ReplacementCorrespondence {
        format: "catallax-import-correspondence/1".into(),
        old_source_sha256: old.sha256.clone(),
        new_source_sha256: new.sha256.clone(),
        verified,
        unresolved,
        unmatched_new: new_entities
            .into_keys()
            .filter(|entity| !matched.contains(entity))
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::incidence_features::common;
    use common::{BoxFixture, box_step};

    // Existing spline_certification.rs repair control, with a 0.001 mm middle pole.
    fn bowed_spline() -> String {
        let text = box_step(BoxFixture::new(["2", "3", "4"]));
        let source = crate::read_step(text.as_bytes()).unwrap();
        let edge = &source.solids[0].faces[0].loops[0].edges[0];
        let edge_line = text
            .lines()
            .find(|line| line.starts_with(&format!("{}=", edge.edge)))
            .unwrap();
        let curve = edge_line.split(',').nth(3).unwrap();
        let curve_line = text
            .lines()
            .find(|line| line.starts_with(&format!("{curve}=")))
            .unwrap();
        let a = &edge.start.point;
        let b = &edge.end.point;
        let mut mid = std::array::from_fn::<_, 3, _>(|i| (&a[i] + &b[i]).mul_pow2(-1));
        mid[0] = &mid[0] + &common::r("0.001");
        let mut points = String::new();
        for (id, point) in [a, &mid, b].into_iter().enumerate() {
            points.push_str(&format!(
                "#{}=CARTESIAN_POINT('',({},{},{}));\n",
                900001 + id,
                common::dec(&point[0]),
                common::dec(&point[1]),
                common::dec(&point[2])
            ));
        }
        text.replace(curve_line, &format!("{curve}=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900001,#900002,#900003),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);"))
            .replace("DATA;", &format!("DATA;\n{points}"))
    }

    #[test]
    fn internal_spline_derivative_is_not_compared_as_raw_source() {
        let (report, work) = with_operation_budget(20_000_000, || {
            let source = bowed_spline();
            let file = crate::read_step(source.as_bytes()).unwrap();
            let mut model = Model::new();
            let bodies = crate::build_file(&mut model, &file, None).unwrap();
            assert!(
                bodies
                    .iter()
                    .any(|body| !model.body(*body).marks.is_empty())
            );
            let bodies = file
                .solids
                .iter()
                .zip(bodies)
                .map(|(solid, body)| (solid.name.clone(), body))
                .collect::<Vec<_>>();
            let snapshot = SourceSnapshot::new(
                crate::sha256_hex(source.as_bytes()),
                file,
                &bodies,
                &[],
                &model,
            );
            compare_inner(&snapshot, &snapshot, &mut Decider::default())
        });
        assert!(!work.exhausted);
        let report = report.unwrap();
        assert!(report.verified.is_empty());
        assert_eq!(report.unresolved.len(), 18);
        assert!(
            report
                .unresolved
                .iter()
                .all(|entry| entry.reason == UnresolvedReason::UnverifiedGeometry)
        );
        assert_eq!(report.unmatched_new.len(), 18);
    }
}
