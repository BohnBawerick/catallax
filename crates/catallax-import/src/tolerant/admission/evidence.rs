//! Borrowed provenance views: a body-local derivative must not disclose siblings.
use serde::{Serialize, Serializer, ser::SerializeStruct};

use super::{TolerantEvidence, healing::HealingChange};

pub(super) struct Scope<'a> {
    evidence: &'a TolerantEvidence,
    body: Option<&'a str>,
}

impl<'a> Scope<'a> {
    pub(super) fn new(evidence: &'a TolerantEvidence, body: Option<&'a str>) -> Self {
        Self { evidence, body }
    }
}

struct Selected<'a, T> {
    records: &'a [T],
    body: &'a str,
    belongs: fn(&T, &str) -> bool,
}

impl<T: Serialize> Serialize for Selected<'_, T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.collect_seq(
            self.records
                .iter()
                .filter(|record| (self.belongs)(record, self.body)),
        )
    }
}

impl Serialize for Scope<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let Some(body) = self.body else {
            return self.evidence.serialize(serializer);
        };
        let evidence = self.evidence;
        let mut result = serializer.serialize_struct("TolerantEvidence", 11)?;
        result.serialize_field("method", &evidence.method)?;
        result.serialize_field("source_sha256", &evidence.source_sha256)?;
        result.serialize_field("authority", &evidence.authority)?;
        result.serialize_field(
            "tolerances",
            &Selected {
                records: &evidence.tolerances,
                body,
                belongs: |record, body| record.body == body,
            },
        )?;
        result.serialize_field(
            "sampled_checks",
            &Selected {
                records: &evidence.sampled_checks,
                body,
                belongs: |record, body| record.body == body,
            },
        )?;
        result.serialize_field(
            "repairs",
            &Selected {
                records: &evidence.repairs,
                body,
                belongs: |record: &HealingChange, selected| record.body() == selected,
            },
        )?;
        result.serialize_field(
            "topology_checks",
            &Selected {
                records: &evidence.topology_checks,
                body,
                belongs: |record, body| record.body == body,
            },
        )?;
        result.serialize_field(
            "topology_measurements",
            &Selected {
                records: &evidence.topology_measurements,
                body,
                belongs: |record, body| record.body == body,
            },
        )?;
        result.serialize_field(
            "warnings",
            &Selected {
                records: &evidence.warnings,
                body,
                belongs: |record, body| record.body.as_deref().is_none_or(|owner| owner == body),
            },
        )?;
        result.serialize_field(
            "native_trims",
            &Selected {
                records: &evidence.native_trims,
                body,
                belongs: |record, body| record.body == body,
            },
        )?;
        result.serialize_field(
            "native_domains",
            &Selected {
                records: &evidence.native_domains,
                body,
                belongs: |record, body| record.body == body,
            },
        )?;
        result.end()
    }
}
