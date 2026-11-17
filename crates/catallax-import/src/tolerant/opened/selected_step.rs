//! Compose selected body artifacts with their retained source occurrence graph.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Refusal};
use serde::Serialize;

use super::{BodyView, OpenedStep, Tier};
use crate::geometry_source::{ParsedStep, StepGraphKind, step_reachable};
use crate::part21::{Entity, Graph, Param};

/// Binding between one source prototype and its actual selected STEP geometry.
#[derive(Clone, Debug, Serialize)]
pub struct SelectedStepBody {
    /// Index in the default opening's source-body inventory.
    pub source_index: usize,
    /// Original source prototype entity, retained as its mapping item.
    pub source_entity: u64,
    /// Actual emitted solid or surface-model roots for this prototype.
    pub roots: Vec<u64>,
    /// Selected opening tier, independent of exact versus numerical proof.
    pub tier: Tier,
    /// Geometry authority of the selected body writer.
    pub proof: &'static str,
    /// Whether the opening records an approximate source interpretation.
    pub approximate: bool,
    /// Emitted face, edge and vertex counts, in that order.
    pub counts: [usize; 3],
    /// Exact selected local tolerance in millimetres, not a global body bound.
    pub tolerance_mm: String,
    /// Length, plane-angle and solid-angle factors of the emitted prototype.
    pub unit_factors: [f64; 3],
}

/// Whole-source STEP with selected geometry and retained occurrence bindings.
#[derive(Debug)]
pub struct SelectedStepFile {
    bytes: Vec<u8>,
    bodies: Vec<SelectedStepBody>,
}

impl SelectedStepFile {
    /// Actual selected-body document bytes.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
    /// Source prototype bindings and separately labelled geometry authority.
    pub fn bodies(&self) -> &[SelectedStepBody] {
        &self.bodies
    }
    /// Consume the document without copying its bytes.
    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("import.selected-document", detail.into())
}

fn source_error(error: impl std::fmt::Debug) -> Refusal {
    bad(format!("{error:?}"))
}

fn add(
    graph: &mut Graph,
    next: &mut u64,
    kind: &str,
    parameters: Vec<Param>,
) -> Result<u64, Refusal> {
    let id = *next;
    *next = next
        .checked_add(1)
        .ok_or_else(|| bad("STEP entity index overflow"))?;
    graph.insert(
        id,
        Entity {
            parts: vec![(kind.into(), parameters)],
        },
    );
    Ok(id)
}

fn remap(value: &mut Param, ids: &BTreeMap<u64, u64>, dec: &mut Decider) -> Result<(), Refusal> {
    dec.charge_work(1)?;
    match value {
        Param::Ref(id) => {
            *id = *ids
                .get(id)
                .ok_or_else(|| bad("selected geometry has an unbound reference"))?
        }
        Param::List(values) | Param::Typed(_, values) => {
            for value in values {
                remap(value, ids, dec)?;
            }
        }
        _ => {}
    }
    Ok(())
}

fn insert_prototype(
    source: &ParsedStep,
    bytes: &[u8],
    source_entity: u64,
    changes: &mut Graph,
    next: &mut u64,
    dec: &mut Decider,
) -> Result<(Vec<u64>, [usize; 3]), Refusal> {
    let written = ParsedStep::read(bytes, dec).map_err(source_error)?;
    let roots: Vec<_> = written
        .model_declarations(dec)
        .map_err(source_error)?
        .into_iter()
        .map(|model| model.entity)
        .collect();
    if roots.is_empty() {
        return Err(bad("selected writer emitted no prototype geometry"));
    }
    let mut representation = None;
    for id in written.entities() {
        dec.charge_work(1)?;
        if let Ok(frame) = written.representation(id)
            && roots.iter().all(|root| frame.items.contains(root))
            && written
                .representation_context(frame.context)
                .is_ok_and(|context| context.dimension == Some(3))
        {
            representation = Some(id);
            break;
        }
    }
    let representation = representation
        .ok_or_else(|| bad("selected prototype has no complete 3D representation"))?;
    let reached = step_reachable(
        written.source_graph(),
        &Graph::new(),
        &[representation],
        dec,
    )?;
    let mut ids = BTreeMap::new();
    for old in &reached {
        let id = *next;
        *next = next
            .checked_add(1)
            .ok_or_else(|| bad("STEP entity index overflow"))?;
        ids.insert(*old, id);
    }
    let mut counts = [0; 3];
    for old in reached {
        let mut entity = written.source_graph()[&old].clone();
        for (index, names) in [
            &["ADVANCED_FACE", "FACE_SURFACE", "FACE"][..],
            &["EDGE_CURVE"][..],
            &["VERTEX_POINT"][..],
        ]
        .iter()
        .enumerate()
        {
            counts[index] += usize::from(
                entity
                    .parts
                    .iter()
                    .any(|(kind, _)| names.contains(&kind.as_str())),
            );
        }
        for (_, parameters) in &mut entity.parts {
            for value in parameters {
                remap(value, &ids, dec)?;
            }
        }
        changes.insert(ids[&old], entity);
    }
    let origin = add(
        changes,
        next,
        "CARTESIAN_POINT",
        vec![
            Param::Str(String::new()),
            Param::List(vec![Param::Number("0.".into()); 3]),
        ],
    )?;
    let frame = add(
        changes,
        next,
        "AXIS2_PLACEMENT_3D",
        vec![
            Param::Str(String::new()),
            Param::Ref(origin),
            Param::Unset,
            Param::Unset,
        ],
    )?;
    let map = add(
        changes,
        next,
        "REPRESENTATION_MAP",
        vec![Param::Ref(frame), Param::Ref(ids[&representation])],
    )?;
    let transform = add(
        changes,
        next,
        "CARTESIAN_TRANSFORMATION_OPERATOR_3D",
        vec![
            Param::Str(String::new()),
            Param::Unset,
            Param::Unset,
            Param::Ref(origin),
            Param::Number("1.".into()),
            Param::Unset,
        ],
    )?;
    let name = source
        .source_graph()
        .get(&source_entity)
        .and_then(|entity| entity.parts.first())
        .and_then(|(_, parameters)| parameters.first())
        .filter(|value| matches!(value, Param::Str(_)))
        .cloned()
        .unwrap_or_else(|| Param::Str(String::new()));
    changes.insert(
        source_entity,
        Entity {
            parts: vec![(
                "MAPPED_ITEM".into(),
                vec![name, Param::Ref(map), Param::Ref(transform)],
            )],
        },
    );
    Ok((roots.into_iter().map(|root| ids[&root]).collect(), counts))
}

impl OpenedStep {
    /// Export the default opening's selected bodies without another admission.
    /// Original source occurrence links and transforms are retained. Each body
    /// keeps its writer's own dimensional context through a representation map.
    /// Every source prototype must have a retained representation binding.
    /// A missing source graph or any retained surface body refuses atomically.
    /// The caller's current allowance covers every body writer and composition.
    pub fn export_step_document(&self, dec: &mut Decider) -> Result<SelectedStepFile, Refusal> {
        let source = self
            .parsed_source()
            .ok_or_else(|| bad("opening has no retained source occurrence graph"))?;
        if self.bodies().is_empty()
            || (0..self.bodies().len())
                .any(|index| matches!(self.body(index), Some(BodyView::Surfaces { .. }) | None))
        {
            return Err(bad(
                "whole-document export requires every source body to have selected geometry",
            ));
        }
        let mut changes = Graph::new();
        let source_entities: Vec<u64> = self
            .source()
            .bodies()
            .iter()
            .map(|body| {
                body.name()
                    .strip_prefix('#')
                    .and_then(|name| name.parse().ok())
                    .ok_or_else(|| bad("source prototype has no entity binding"))
            })
            .collect::<Result<_, _>>()?;
        let source_roots: BTreeSet<_> = source_entities.iter().copied().collect();
        let mut represented = BTreeSet::new();
        for id in source.entities() {
            dec.charge_work(1)?;
            if let Ok(representation) = source.representation(id) {
                represented.extend(representation.items);
            }
        }
        if !source_roots.is_subset(&represented) {
            return Err(bad(
                "source prototype has no representation binding; export body-local geometry separately",
            ));
        }
        let mut next = source
            .entities()
            .last()
            .unwrap_or(0)
            .checked_add(1)
            .ok_or_else(|| bad("STEP entity index overflow"))?;
        let mut bodies = Vec::with_capacity(self.bodies().len());
        for index in 0..self.bodies().len() {
            let source_entity = source_entities[index];
            let (bytes, proof, tolerance_mm, unit_factors) = match self.body(index).unwrap() {
                BodyView::Certified { model, body } => {
                    let mut bytes = Vec::new();
                    let options = catallax_step::StepOptions::default();
                    catallax_step::write_step(model, &[body], &options, dec, &mut bytes)?;
                    (
                        bytes,
                        "certified",
                        options.curve_tolerance.to_string(),
                        [1.; 3],
                    )
                }
                BodyView::Numerical(body) => {
                    let file = catallax_step::numerical::build(body, &Default::default(), dec)?;
                    file.check_written(file.bytes(), dec)?;
                    (
                        file.bytes().to_vec(),
                        "numerical",
                        body.options().tolerance_mm.to_string(),
                        [1.; 3],
                    )
                }
                BodyView::Tolerant { import, body } => {
                    let file = import.export_body(body, dec).map_err(source_error)?;
                    let selected = &import.bodies()[body];
                    (
                        file.into_bytes(),
                        "native_sampled",
                        selected.definition().maximum_tolerance_mm().to_string(),
                        import.source_unit_factors(),
                    )
                }
                BodyView::Surfaces { .. } => {
                    unreachable!("all source bodies were checked before writing")
                }
            };
            let (roots, counts) =
                insert_prototype(source, &bytes, source_entity, &mut changes, &mut next, dec)?;
            bodies.push(SelectedStepBody {
                source_index: index,
                source_entity,
                roots,
                tier: self.bodies()[index].tier(),
                proof,
                approximate: self.bodies()[index].approximate(),
                counts,
                tolerance_mm,
                unit_factors,
            });
        }
        for (&id, entity) in source.source_graph() {
            for (kind, parameters) in &entity.parts {
                if kind.ends_with("SHAPE_REPRESENTATION")
                    && kind != "SHAPE_REPRESENTATION"
                    && let Some(Param::List(items)) = parameters.get(1)
                    && items
                        .iter()
                        .any(|item| matches!(item, Param::Ref(root) if source_roots.contains(root)))
                {
                    let mut replacement = entity.clone();
                    for (name, _) in &mut replacement.parts {
                        if name == kind {
                            *name = "SHAPE_REPRESENTATION".into();
                        }
                    }
                    changes.insert(id, replacement);
                }
            }
        }
        let evidence = serde_json::to_string(&serde_json::json!({
            "export_scope":"selected_source_document","source_sha256":source.source_sha256(),
            "bodies":bodies,"occurrence_links":"retained_source_entities"}))
        .map_err(|error| bad(error.to_string()))?;
        let bytes = source.write_step_graph(
            &changes,
            source.schemas().map_err(source_error)?,
            &evidence,
            None,
            StepGraphKind::Selected,
            dec,
        )?;
        Ok(SelectedStepFile { bytes, bodies })
    }
}
