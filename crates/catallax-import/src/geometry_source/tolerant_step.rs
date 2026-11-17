//! Rebuild a tolerant STEP document from its parsed entity model and repairs.
//! This is serialization of a lower-tier B-rep, not checked-Model export.

use std::collections::BTreeSet;
use std::fmt::Write;

use catallax_numbers::{Decider, Refusal};

use super::ParsedStep;
use crate::part21::{Entity, Graph, Param};

impl ParsedStep {
    /// Immutable declarations for the native tolerant topology mapper.
    pub(crate) fn source_graph(&self) -> &Graph {
        &self.graph
    }

    /// Serialize actual parsed geometry and repaired topology, not input bytes.
    /// The evidence travels in FILE_DESCRIPTION; no exact authority is inferred.
    pub(crate) fn write_tolerant_step(
        &self,
        changes: &Graph,
        schemas: &[String],
        evidence_json: &str,
        dec: &mut Decider,
    ) -> Result<Vec<u8>, Refusal> {
        let mut output = String::with_capacity(evidence_json.len() + self.graph.len() * 80);
        output.push_str("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((");
        quoted(&mut output, "CATALLAX TOLERANT: sampled fits, measured local tolerances; not certified");
        output.push(',');
        quoted(&mut output, evidence_json);
        output.push_str("),'2;1');\nFILE_NAME('catallax-tolerant.step','',(''),(''),'Catallax','Catallax','');\nFILE_SCHEMA((");
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
        output.push_str("));\nENDSEC;\nDATA;\n");
        let ids: BTreeSet<_> = self.graph.keys().chain(changes.keys()).copied().collect();
        for id in ids {
            dec.check_deadline()?;
            dec.charge_work(1)?;
            let entity = changes.get(&id).or_else(|| self.graph.get(&id)).expect("union key");
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

fn validate(entity: &Entity, source: &Graph, changes: &Graph, dec: &mut Decider) -> Result<(), Refusal> {
    fn visit(value: &Param, source: &Graph, changes: &Graph, dec: &mut Decider) -> Result<(), Refusal> {
        dec.charge_work(1)?;
        match value {
            Param::Ref(id) if !source.contains_key(id) && !changes.contains_key(id) => Err(Refusal::new(
                "import.tolerant-broken-reference",
                format!("repaired STEP references absent entity #{id}"),
            )),
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
