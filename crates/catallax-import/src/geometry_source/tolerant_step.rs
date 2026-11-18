//! Serialize labelled body graphs with their retained source context.
//! Geometry authority remains separate for native and selected-body artifacts.

use std::collections::BTreeSet;
use std::fmt::Write;

use catallax_numbers::{Decider, Refusal};

use super::ParsedStep;
use crate::part21::{Entity, Graph, Param};

pub(crate) enum StepGraphKind {
    Native,
    Selected,
}

impl ParsedStep {
    /// Immutable declarations for the native tolerant topology mapper.
    pub(crate) fn source_graph(&self) -> &Graph {
        &self.graph
    }

    /// Serialize actual parsed geometry and repaired topology, not input bytes.
    /// The evidence travels in ignored report comments; no exact authority is inferred.
    pub(crate) fn write_step_graph(
        &self,
        changes: &Graph,
        schemas: &[String],
        evidence_json: &str,
        roots: Option<&[u64]>,
        kind: StepGraphKind,
        dec: &mut Decider,
    ) -> Result<Vec<u8>, Refusal> {
        let ids = match roots {
            Some(roots) => step_reachable(&self.graph, changes, roots, dec)?,
            None => self.graph.keys().chain(changes.keys()).copied().collect(),
        };
        let mut output = String::with_capacity(evidence_json.len() + ids.len() * 80);
        output.push_str("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((");
        let (description, file_name) = match kind {
            StepGraphKind::Native => (
                "CATALLAX TOLERANT: sampled native B-rep, not certified; per-entity millimetre tolerances and repair/check record: CATALLAX_REPORT_JSON comments",
                "catallax-tolerant.step",
            ),
            StepGraphKind::Selected => (
                "CATALLAX SELECTED BODY DOCUMENT: retained source occurrences and separately labelled selected body exports; CATALLAX_REPORT_JSON comments",
                "catallax-selected.step",
            ),
        };
        quoted(&mut output, description);
        output.push_str("),'2;1');\nFILE_NAME(");
        quoted(&mut output, file_name);
        output.push_str(",'',(''),(''),'Catallax','Catallax','');\nFILE_SCHEMA((");
        if schemas.is_empty() {
            quoted(&mut output, "AUTOMOTIVE_DESIGN");
        } else {
            for (index, schema) in schemas.iter().enumerate() {
                if index > 0 {
                    output.push(',');
                }
                quoted(&mut output, schema);
            }
        }
        output.push_str("));\nENDSEC;\n");
        // Concatenating comment payloads reconstructs one JSON report. Escaped
        // JSON slashes cannot start or terminate STEP comments, even in names.
        let mut remaining = evidence_json;
        while !remaining.is_empty() {
            dec.check_deadline()?;
            let mut end = remaining.len().min(2048);
            while !remaining.is_char_boundary(end) {
                end -= 1;
            }
            output.push_str("/* CATALLAX_REPORT_JSON ");
            for character in remaining[..end].chars() {
                if character == '/' {
                    output.push('\\');
                }
                output.push(character);
            }
            output.push_str(" */\n");
            remaining = &remaining[end..];
        }
        output.push_str("DATA;\n");
        for id in ids {
            dec.check_deadline()?;
            dec.charge_work(1)?;
            let entity = changes
                .get(&id)
                .or_else(|| self.graph.get(&id))
                .expect("union key");
            validate(entity, &self.graph, changes, dec)?;
            write!(&mut output, "#{id}=").expect("string formatting");
            if entity.parts.len() == 1 {
                let (kind, parameters) = &entity.parts[0];
                output.push_str(kind);
                parameters_text(&mut output, parameters);
            } else {
                output.push('(');
                for (kind, parameters) in &entity.parts {
                    output.push_str(kind);
                    parameters_text(&mut output, parameters);
                }
                output.push(')');
            }
            output.push_str(";\n");
        }
        output.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        dec.charge_work(output.len() as u64)?;
        Ok(output.into_bytes())
    }
}

pub(crate) fn step_reachable(
    source: &Graph,
    changes: &Graph,
    roots: &[u64],
    dec: &mut Decider,
) -> Result<BTreeSet<u64>, Refusal> {
    fn references(value: &Param, pending: &mut Vec<u64>, dec: &mut Decider) -> Result<(), Refusal> {
        dec.charge_work(1)?;
        match value {
            Param::Ref(id) => pending.push(*id),
            Param::List(values) | Param::Typed(_, values) => {
                for value in values {
                    references(value, pending, dec)?;
                }
            }
            _ => {}
        }
        Ok(())
    }
    let mut pending = roots.to_vec();
    let mut seen = BTreeSet::new();
    while let Some(id) = pending.pop() {
        if !seen.insert(id) {
            continue;
        }
        let entity = changes
            .get(&id)
            .or_else(|| source.get(&id))
            .ok_or_else(|| {
                Refusal::new(
                    "import.tolerant-broken-reference",
                    format!("absent entity #{id}"),
                )
            })?;
        for (_, parameters) in &entity.parts {
            for value in parameters {
                references(value, &mut pending, dec)?;
            }
        }
    }
    Ok(seen)
}

pub(crate) fn product_representations(
    graph: &Graph,
    dec: &mut Decider,
) -> Result<Vec<u64>, Refusal> {
    fn product_bound(graph: &Graph, root: u64, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut pending = vec![root];
        let mut seen = BTreeSet::new();
        while let Some(id) = pending.pop() {
            dec.charge_work(1)?;
            if !seen.insert(id) {
                continue;
            }
            let Some(entity) = graph.get(&id) else {
                continue;
            };
            if entity.parts.iter().any(|(kind, _)| {
                kind == "PRODUCT_DEFINITION"
                    || kind == "PRODUCT_DEFINITION_RELATIONSHIP"
                    || kind.ends_with("USAGE_OCCURRENCE")
            }) {
                return Ok(true);
            }
            for (kind, parameters) in &entity.parts {
                if matches!(
                    kind.as_str(),
                    "PRODUCT_DEFINITION_SHAPE" | "PROPERTY_DEFINITION" | "SHAPE_ASPECT"
                ) && let Some(Param::Ref(target)) = parameters.get(2)
                {
                    pending.push(*target);
                }
            }
        }
        Ok(false)
    }
    let mut roots = Vec::new();
    for entity in graph.values() {
        dec.charge_work(1)?;
        if entity.part("SHAPE_DEFINITION_REPRESENTATION").is_none() {
            continue;
        }
        let parameters = entity
            .part("SHAPE_DEFINITION_REPRESENTATION")
            .filter(|parameters| parameters.len() >= 2)
            .or_else(|| entity.part("PROPERTY_DEFINITION_REPRESENTATION"));
        if let Some([Param::Ref(definition), Param::Ref(representation), ..]) = parameters
            && product_bound(graph, *definition, dec)?
        {
            roots.push(*representation);
        }
    }
    Ok(roots)
}

fn validate(
    entity: &Entity,
    source: &Graph,
    changes: &Graph,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    fn visit(
        value: &Param,
        source: &Graph,
        changes: &Graph,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        dec.charge_work(1)?;
        match value {
            Param::Ref(id) if !source.contains_key(id) && !changes.contains_key(id) => {
                Err(Refusal::new(
                    "import.tolerant-broken-reference",
                    format!("repaired STEP references absent entity #{id}"),
                ))
            }
            Param::List(values) | Param::Typed(_, values) => {
                for value in values {
                    visit(value, source, changes, dec)?;
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }
    for (_, parameters) in &entity.parts {
        for parameter in parameters {
            visit(parameter, source, changes, dec)?;
        }
    }
    Ok(())
}

fn quoted(output: &mut String, value: &str) {
    output.push('\'');
    for character in value.chars() {
        output.push(character);
        if character == '\'' {
            output.push('\'');
        }
    }
    output.push('\'');
}

fn parameters_text(output: &mut String, values: &[Param]) {
    output.push('(');
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            output.push(',');
        }
        match value {
            Param::Ref(id) => write!(output, "#{id}").expect("string formatting"),
            Param::Number(text) => output.push_str(text),
            Param::Str(text) => quoted(output, text),
            Param::Enum(text) => {
                output.push('.');
                output.push_str(text);
                output.push('.');
            }
            Param::List(values) => parameters_text(output, values),
            Param::Typed(name, values) => {
                output.push_str(name);
                parameters_text(output, values);
            }
            Param::Unset => output.push('$'),
            Param::Derived => output.push('*'),
        }
    }
    output.push(')');
}
